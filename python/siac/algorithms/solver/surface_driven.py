"""Surface-driven aerosol solver.

An alternative to the Bayesian multi-grid inversion: instead of an L-BFGS-B
descent, sweep a fixed AOT axis at the prior TCWV, atmospherically correct the
TOA to BOA at each node, and score that against the surface prior. The per-node
cost cube is built here (vectorised numpy), then the heavy spatial median-pool
+ Gaussian-backstop + per-pixel arg-min is delegated to the Rust kernel
``surface_driven_pool_argmin``. This is the cheap, robust scheme the
experimental harness validated; it is selected by
``algorithms.solver.method = "surface_driven"``.

TCWV is held at the prior for the sweep and carried through to the solved state
(the surface-driven retrieval constrains AOT, not column water vapour).
"""

from __future__ import annotations

import logging
import os
from dataclasses import replace
from pathlib import Path
from typing import TYPE_CHECKING, Any, cast

import numpy as np
import xarray as xr

from siac._rust_compat import surface_driven_pool_argmin
from siac.algorithms.solver._grid_search import (
    build_candidate_coeff_provider,
    prepare_grid_search_observations,
)
from siac.algorithms.solver.multigrid import SolverResult, _log_aot_axis
from siac.domain.protocols import rt_optional_capability
from siac.runtime import GeometryAngles

if TYPE_CHECKING:
    from siac.config.algorithms import SolverAlgorithmConfig
    from siac.domain import SensorBand
    from siac.domain.protocols import RTModelBackend
    from siac.runtime import (
        AtmosphericState,
        SolverInputBundle,
        SurfacePrior,
    )

#: BOA-uncertainty floor for the per-band chi-square (matches the harness
#: ``UNC_FLOOR``; keeps a single clean band from dominating the cost).
_MIN_BOA_UNC = 0.006


logger = logging.getLogger(__name__)

_CANONICAL_6S_AEROSOL_PROFILES = (
    "continental",
    "maritime",
    "urban",
    "desert",
    "biomass_burning",
)


def _coord_mean(template: xr.DataArray, names: tuple[str, ...]) -> float | None:
    for name in names:
        if name not in template.coords:
            continue
        values = np.asarray(template.coords[name].values, dtype=np.float64)
        finite = values[np.isfinite(values)]
        if finite.size:
            return float(np.mean(finite))
    return None


def _scene_center_lon_lat(template: xr.DataArray) -> tuple[float, float]:
    lon = _coord_mean(template, ("lon", "longitude"))
    lat = _coord_mean(template, ("lat", "latitude"))
    if lon is not None and lat is not None:
        if lon > 180.0:
            lon = ((lon + 180.0) % 360.0) - 180.0
        return float(lon), float(lat)

    x = _coord_mean(template, ("x",))
    y = _coord_mean(template, ("y",))
    if x is None or y is None:
        raise ValueError(
            "surface_driven_aerosol_species requires solver-grid lon/lat coordinates or x/y "
            "coordinates with CRS metadata."
        )

    try:
        import rioxarray  # noqa: F401

        crs = template.rio.crs
    except Exception:
        crs = None
    if crs is not None:
        from pyproj import Transformer

        lon_t, lat_t = Transformer.from_crs(crs, "EPSG:4326", always_xy=True).transform(x, y)
        return float(lon_t), float(lat_t)

    if -180.0 <= x <= 360.0 and -90.0 <= y <= 90.0:
        lon = ((x + 180.0) % 360.0) - 180.0 if x > 180.0 else x
        return float(lon), float(y)

    raise ValueError(
        "surface_driven_aerosol_species requires CRS metadata to derive scene lon/lat from "
        "projected solver-grid coordinates."
    )


def _rt_month(rt_model: object) -> int:
    observation_time = getattr(rt_model, "observation_time", None)
    if observation_time is not None and hasattr(observation_time, "month"):
        return int(observation_time.month)
    rt_config = getattr(rt_model, "_config", None)
    return int(getattr(rt_config, "month", 1))


def _species_candidate_rt_models(
    *,
    rt_model: object,
    config: object,
    template: xr.DataArray,
) -> list[RTModelBackend]:
    mode = str(getattr(config, "surface_driven_aerosol_species", "none"))
    if mode == "none":
        return []
    if mode not in {"cci_climatology", "cci_climatology_exact", "canonical_6s"}:
        raise ValueError(f"Unsupported surface-driven aerosol species mode {mode!r}.")
    if str(getattr(rt_model, "backend_name", "")) != "sixs":
        raise ValueError(
            "surface_driven_aerosol_species requires algorithms.rt.backend=\"sixs\"."
        )

    clone_rt = rt_optional_capability(rt_model, "with_rt_setup")
    base_setup = getattr(rt_model, "rt_setup", None)
    if clone_rt is None or base_setup is None or not hasattr(base_setup, "model_copy"):
        raise ValueError(
            "Native 6S species mode requires an RT backend exposing rt_setup and with_rt_setup()."
        )

    if mode == "canonical_6s":
        from siac.config.algorithms import RTAerosolSetupConfig

        aerosol_setups = tuple(
            RTAerosolSetupConfig(profile=profile)
            for profile in _CANONICAL_6S_AEROSOL_PROFILES
        )
    elif mode == "cci_climatology_exact":
        lon, lat = _scene_center_lon_lat(template)
        month = _rt_month(rt_model)
        from siac.algorithms.rt.aerosol_species import climatology_cci_aerosol_setup

        aerosol_setups = (climatology_cci_aerosol_setup(lon, lat, month),)
    else:
        lon, lat = _scene_center_lon_lat(template)
        month = _rt_month(rt_model)
        n_candidates = int(getattr(config, "surface_driven_aerosol_species_candidates", 3))
        from siac.algorithms.rt.aerosol_species import candidate_cci_aerosol_setups

        aerosol_setups = candidate_cci_aerosol_setups(lon, lat, month, n=n_candidates)
    return [
        clone_rt(base_setup.model_copy(update={"aerosol": aerosol_setup}, deep=True))
        for aerosol_setup in aerosol_setups
    ]


def _acixthree_aot_axis() -> np.ndarray:
    """The validated fine non-uniform AOT axis (denser at low AOD).

    Reproduces the surface-driven harness axis: 0.01 steps to 0.2, 0.025 to 0.5,
    0.05 to 1.5, then 0.1 to 2.5, extended 0.25 steps to 4.0 (~68 nodes; the
    packaged LUT covers AOT 0.001-5.0). Without the extension, truths above the
    old 2.5 cap (e.g. wildfire-smoke matchups at 3.2-3.5) were unreachable by
    construction. The fine low-AOD spacing avoids the coarse-log quantisation
    that snaps clean retrievals onto sparse nodes.
    """
    return cast(
        "np.ndarray",
        np.concatenate(
            [
                np.arange(0.01, 0.2, 0.01),
                np.arange(0.2, 0.5, 0.025),
                np.arange(0.5, 1.5, 0.05),
                np.arange(1.5, 2.6, 0.1),
                np.arange(2.75, 4.01, 0.25),
            ]
        ).astype(np.float32),
    )


def _shape_cost_node(
    boa: np.ndarray,
    ref: np.ndarray,
    sig: np.ndarray,
    band_indices: list[int],
    anchor_index: int,
) -> np.ndarray:
    """Spectral-shape surface cost for one AOT node (NaN-propagating).

    Absolute brightness anchor on the reddest band + ratio (bluer/anchor) terms,
    with ratio sigma by error propagation. Cancels multiplicative surface-
    brightness errors and B01's absolute dark bias while keeping the deep-blue
    contrast that carries the AOD signal.
    """
    a = anchor_index
    sig_a = np.maximum(sig[a], 1e-6)
    ref_a = np.maximum(ref[a], 1e-4)
    boa_a = np.maximum(boa[a], 1e-4)
    with np.errstate(invalid="ignore", divide="ignore"):
        d = (boa[a] - ref[a]) ** 2 / sig_a**2
        for bi in band_indices:
            if bi == a:
                continue
            rm = boa[bi] / boa_a
            rr = ref[bi] / ref_a
            sr = np.abs(rr) * np.sqrt(
                (sig[bi] / np.maximum(ref[bi], 1e-4)) ** 2 + (sig[a] / ref_a) ** 2
            )
            d = d + (rm - rr) ** 2 / np.maximum(sr, 1e-3) ** 2
    return cast("np.ndarray", d)


def _abs_cost_node(
    boa: np.ndarray,
    ref: np.ndarray,
    inv_var: np.ndarray,
    band_indices: list[int],
) -> np.ndarray:
    """Absolute per-band chi-square for one AOT node (NaN-propagating)."""
    d = np.zeros(boa.shape[1:], dtype=np.float64)
    with np.errstate(invalid="ignore"):
        for bi in band_indices:
            diff = boa[bi] - ref[bi]
            d = d + diff * diff * inv_var[bi]
    return cast("np.ndarray", d)


def _trimmed_chi2_cost_node(
    boa: np.ndarray,
    ref: np.ndarray,
    inv_var: np.ndarray,
    band_indices: list[int],
) -> np.ndarray:
    """Sum chi-square terms after dropping the worst visible-band mismatch."""
    if len(band_indices) < 3:
        raise ValueError("trimmed chi-square cost requires at least three bands")

    terms = np.stack(
        [
            np.square(boa[band_index] - ref[band_index]) * inv_var[band_index]
            for band_index in band_indices
        ],
        axis=0,
    )
    keep = len(band_indices) - 1
    finite_count = np.count_nonzero(np.isfinite(terms), axis=0)
    ordered = np.partition(np.where(np.isfinite(terms), terms, np.inf), keep - 1, axis=0)
    cost = np.sum(ordered[:keep], axis=0)
    return cast("np.ndarray", np.where(finite_count >= keep, cost, np.nan))


def _additive_offset_cost_node(
    boa: np.ndarray,
    ref: np.ndarray,
    sig: np.ndarray,
    band_indices: list[int],
) -> np.ndarray:
    """Profile a spectrally flat surface-reflectance offset at one AOT node.

    A common additive mismatch can arise from residual surface normalization or
    local illumination. Profiling that one nuisance parameter leaves the
    wavelength-dependent residual, rather than absolute prior brightness, to
    constrain AOD. The surface source and visible spectral shape are unchanged.
    """
    if len(band_indices) < 2:
        raise ValueError("additive-offset cost requires at least two bands")

    numerator = np.zeros(boa.shape[1:], dtype=np.float64)
    denominator = np.zeros(boa.shape[1:], dtype=np.float64)
    with np.errstate(invalid="ignore", divide="ignore"):
        for band_index in band_indices:
            variance = np.square(np.maximum(sig[band_index], 1e-6))
            weight = 1.0 / variance
            numerator += (boa[band_index] - ref[band_index]) * weight
            denominator += weight
        offset = numerator / denominator
        cost = np.zeros(boa.shape[1:], dtype=np.float64)
        for band_index in band_indices:
            variance = np.square(np.maximum(sig[band_index], 1e-6))
            residual = boa[band_index] - ref[band_index] - offset
            cost += np.square(residual) / variance
    return cast("np.ndarray", cost)


def _profile_scale_cost_node(
    boa: np.ndarray,
    ref: np.ndarray,
    sig: np.ndarray,
    band_indices: list[int],
    *,
    scale_sigma: float,
) -> np.ndarray:
    """Profile a constrained common surface-brightness scale at one AOT node.

    The surface source remains unchanged.  A single multiplicative nuisance
    parameter absorbs small scene-wide brightness mismatch while its Gaussian
    prior prevents the surface term from explaining arbitrary atmospheric
    signal.
    """
    if scale_sigma <= 0.0:
        raise ValueError("scale_sigma must be positive")
    if not band_indices:
        raise ValueError("profile-scale cost requires at least one band")

    precision = 1.0 / float(scale_sigma) ** 2
    numerator = np.full(boa.shape[1:], precision, dtype=np.float64)
    denominator = np.full(boa.shape[1:], precision, dtype=np.float64)
    with np.errstate(invalid="ignore", divide="ignore"):
        for band_index in band_indices:
            variance = np.square(np.maximum(sig[band_index], 1e-6))
            weight = 1.0 / variance
            numerator += ref[band_index] * boa[band_index] * weight
            denominator += np.square(ref[band_index]) * weight
        scale = np.clip(numerator / denominator, 0.5, 1.5)
        cost = precision * np.square(scale - 1.0)
        for band_index in band_indices:
            variance = np.square(np.maximum(sig[band_index], 1e-6))
            residual = boa[band_index] - scale * ref[band_index]
            cost += np.square(residual) / variance
    return cast("np.ndarray", cost)


def _loo_scale_cost_node(
    boa: np.ndarray,
    ref: np.ndarray,
    sig: np.ndarray,
    band_indices: list[int],
    *,
    scale_sigma: float,
) -> np.ndarray:
    """Score each band after profiling brightness on the other visible bands.

    This is a leave-one-band-out predictive likelihood.  Candidate AOD is
    rewarded only when one common surface scale inferred without a band also
    predicts that band's corrected BOA, reducing the incentive to overfit a
    biased surface prior in the solve bands.
    """
    if scale_sigma <= 0.0:
        raise ValueError("scale_sigma must be positive")
    if len(band_indices) < 3:
        raise ValueError("leave-one-band-out scale cost requires at least three bands")

    precision = 1.0 / float(scale_sigma) ** 2
    cost = np.zeros(boa.shape[1:], dtype=np.float64)
    with np.errstate(invalid="ignore", divide="ignore"):
        for held_out in band_indices:
            numerator = np.full(boa.shape[1:], precision, dtype=np.float64)
            denominator = np.full(boa.shape[1:], precision, dtype=np.float64)
            for band_index in band_indices:
                if band_index == held_out:
                    continue
                variance = np.square(np.maximum(sig[band_index], 1e-6))
                weight = 1.0 / variance
                numerator += ref[band_index] * boa[band_index] * weight
                denominator += np.square(ref[band_index]) * weight
            scale = np.clip(numerator / denominator, 0.5, 1.5)
            held_variance = np.square(np.maximum(sig[held_out], 1e-6))
            held_residual = boa[held_out] - scale * ref[held_out]
            cost += np.square(held_residual) / held_variance
            cost += precision * np.square(scale - 1.0) / len(band_indices)
    return cast("np.ndarray", cost)


def _median_per_node(cube: np.ndarray) -> np.ndarray:
    flat = cube.reshape(cube.shape[0], -1)
    out = np.full(cube.shape[0], np.nan, dtype=np.float64)
    for index in range(flat.shape[0]):
        finite = flat[index][np.isfinite(flat[index])]
        if finite.size:
            out[index] = float(np.median(finite))
    return out


def _cost_curve_diagnostics(cube: np.ndarray, aot_axis: np.ndarray) -> dict[str, Any]:
    node_cost = _median_per_node(cube)
    finite_indices = np.flatnonzero(np.isfinite(node_cost))
    diagnostics: dict[str, Any] = {
        "surface_cost_curve_valid_nodes": int(finite_indices.size),
        "surface_cost_curve_curvature": None,
        "surface_cost_curve_min_at_edge": None,
    }
    if finite_indices.size == 0:
        return diagnostics

    ordered = finite_indices[np.argsort(node_cost[finite_indices])]
    best = int(ordered[0])
    min_at_edge = best == 0 or best == node_cost.size - 1
    diagnostics.update(
        {
            "surface_cost_curve_min_aot": float(aot_axis[best]),
            "surface_cost_curve_min_cost": float(node_cost[best]),
            "surface_cost_curve_min_at_edge": bool(min_at_edge),
        }
    )
    if ordered.size > 1:
        second = int(ordered[1])
        delta = float(node_cost[second] - node_cost[best])
        diagnostics.update(
            {
                "surface_cost_curve_second_aot": float(aot_axis[second]),
                "surface_cost_curve_second_delta": delta,
                "surface_cost_curve_relative_second_delta": delta
                / max(abs(float(node_cost[best])), 1e-6),
            }
        )
    if 0 < best < node_cost.size - 1:
        left = node_cost[best - 1]
        right = node_cost[best + 1]
        if np.isfinite(left) and np.isfinite(right):
            diagnostics["surface_cost_curve_curvature"] = float(
                left + right - 2.0 * node_cost[best]
            )
    return diagnostics


def _backstop_conflict_decision(
    *,
    surface_min_aot: object,
    aot_prior: np.ndarray,
    calibrated_unc: np.ndarray,
    z_threshold: float | None,
) -> tuple[bool, dict[str, Any]]:
    """Decide whether a tile-wide positive surface/prior conflict is material.

    The statistic uses only the surface likelihood and atmospheric prior that
    already exist in the retrieval.  It does not inspect a validation target.
    """
    configured = z_threshold is not None
    threshold = float(z_threshold) if configured else None
    if threshold is not None and threshold <= 0.0:
        raise ValueError("surface_driven_backstop_conflict_z must be positive.")

    try:
        surface_min = float(surface_min_aot)
    except (TypeError, ValueError):
        surface_min = float("nan")
    prior_values = np.asarray(aot_prior, dtype=np.float64)
    prior_values = prior_values[np.isfinite(prior_values)]
    sigma_values = np.asarray(calibrated_unc, dtype=np.float64)
    sigma_values = sigma_values[np.isfinite(sigma_values) & (sigma_values > 0.0)]
    prior_median = float(np.median(prior_values)) if prior_values.size else float("nan")
    sigma_median = float(np.median(sigma_values)) if sigma_values.size else float("nan")
    standardized = (
        (surface_min - prior_median) / sigma_median
        if np.isfinite(surface_min)
        and np.isfinite(prior_median)
        and np.isfinite(sigma_median)
        and sigma_median > 0.0
        else float("nan")
    )
    use_flat = bool(
        threshold is not None and np.isfinite(standardized) and standardized > threshold
    )
    diagnostics: dict[str, Any] = {
        "surface_backstop_conflict_gate_configured": configured,
        "surface_backstop_conflict_surface_min_aot": (
            surface_min if np.isfinite(surface_min) else None
        ),
        "surface_backstop_conflict_prior_aot_median": (
            prior_median if np.isfinite(prior_median) else None
        ),
        "surface_backstop_conflict_calibrated_sigma_median": (
            sigma_median if np.isfinite(sigma_median) else None
        ),
        "surface_backstop_conflict_standardized_positive": (
            standardized if np.isfinite(standardized) else None
        ),
        "surface_backstop_conflict_z_threshold": threshold,
        "surface_backstop_conflict_decision": "flat50" if use_flat else "calibrated",
    }
    return use_flat, diagnostics


def _select_auto2_solution(
    *,
    aot_main: np.ndarray,
    unc_main: np.ndarray,
    cost_main: np.ndarray,
    aot_abs: np.ndarray,
    unc_abs: np.ndarray,
    cost_abs: np.ndarray,
    clean_threshold: float,
    high_threshold: float,
    cost_gain: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, dict[str, Any]]:
    """Combine auto2 shape and absolute solves with a cost-aware materiality check.

    The legacy rule is to use shape only in ``clean_threshold .. high_threshold`` and
    absolute elsewhere. That can force a low-cost 0.01 result whenever the absolute
    solve is unstable at the regime boundary. This variant keeps the same regime split
    but requires absolute to be materially better (or shape to be unavailable).
    """
    aot_main_arr = np.asarray(aot_main, dtype=np.float32)
    unc_main_arr = np.asarray(unc_main, dtype=np.float32)
    cost_main_arr = np.asarray(cost_main, dtype=np.float32)
    aot_abs_arr = np.asarray(aot_abs, dtype=np.float32)
    unc_abs_arr = np.asarray(unc_abs, dtype=np.float32)
    cost_abs_arr = np.asarray(cost_abs, dtype=np.float32)

    if aot_main_arr.shape != aot_abs_arr.shape:
        raise ValueError(
            f"shape and abs auto2 solve arrays must match: main={aot_main_arr.shape!r}, "
            f"abs={aot_abs_arr.shape!r}"
        )
    if aot_main_arr.shape != cost_main_arr.shape or aot_main_arr.shape != cost_abs_arr.shape:
        raise ValueError("shape and abs auto2 solve arrays must all have matching shapes.")

    main_valid = np.isfinite(aot_main_arr) & np.isfinite(cost_main_arr)
    abs_valid = np.isfinite(aot_abs_arr) & np.isfinite(cost_abs_arr)

    # Legacy regime split (shape in clean→moderate, abs in tails).
    abs_tails = np.isfinite(aot_abs_arr) & (
        (aot_abs_arr < float(clean_threshold)) | (aot_abs_arr > float(high_threshold))
    )

    # Start from the legacy selection.
    use_abs = abs_tails & abs_valid

    # If both solutions are finite in the same pixel, only keep absolute when it is
    # materially better. This prevents hard-rail 0.01 carry-through in high-AOD
    # cases where the absolute branch is not clearly preferred.
    both_valid = abs_tails & main_valid & abs_valid
    if both_valid.any() and cost_gain > 0.0:
        both_main = cost_main_arr[both_valid]
        both_abs = cost_abs_arr[both_valid]
        abs_improves = both_abs <= both_main * (1.0 - float(cost_gain))
        # Allow the no-cost equal/near-zero tie to switch through only when close to
        # exact equal (keeps low-cost clean-tail behavior stable).
        zero_tie = (both_main <= 0.0) & (both_abs <= 1.0e-8)
        use_abs[both_valid] = abs_improves | zero_tie

    # Backward-compatible fallback: if absolute is NaN in a regime, keep shape even
    # in tails (except where shape is also NaN).
    both_nan = abs_tails & ~abs_valid
    use_abs[both_nan] = False

    selected_aot = np.where(use_abs, aot_abs_arr, aot_main_arr)
    selected_unc = np.where(use_abs, unc_abs_arr, unc_main_arr)
    selected_cost = np.where(use_abs, cost_abs_arr, cost_main_arr)

    diagnostics: dict[str, Any] = {
        "surface_auto2_abs_tail_pixels": int(np.count_nonzero(abs_tails)),
        "surface_auto2_abs_selected_pixels": int(np.count_nonzero(use_abs)),
        "surface_auto2_shape_selected_pixels": int(np.count_nonzero(~use_abs)),
        "surface_auto2_abs_cost_gain": float(cost_gain),
    }
    # Keep the most useful summary stats as stable scalar telemetry.
    selected_ratio = float(np.mean(use_abs)) if use_abs.size else 0.0
    diagnostics["surface_auto2_abs_selected_fraction"] = selected_ratio
    return selected_aot, selected_unc, selected_cost, diagnostics


def _band_curve_diagnostics(
    *,
    band_names: list[str],
    aot_axis: np.ndarray,
    band_node_cost: np.ndarray,
    band_node_residual: np.ndarray,
    final_aot: float,
) -> dict[str, Any]:
    diagnostics: dict[str, Any] = {}
    final_node = int(np.nanargmin(np.abs(aot_axis - final_aot))) if np.isfinite(final_aot) else -1
    if final_node >= 0:
        diagnostics["surface_final_aot_node"] = float(aot_axis[final_node])

    argmins: list[float] = []
    for band_index, band_name in enumerate(band_names):
        node_cost = band_node_cost[band_index]
        finite = np.flatnonzero(np.isfinite(node_cost))
        if finite.size:
            best = int(finite[np.argmin(node_cost[finite])])
            argmin_aot = float(aot_axis[best])
            argmins.append(argmin_aot)
            diagnostics[f"surface_band_{band_name}_argmin_aot"] = argmin_aot
            diagnostics[f"surface_band_{band_name}_argmin_cost"] = float(node_cost[best])
        if final_node >= 0:
            final_cost = float(band_node_cost[band_index, final_node])
            final_residual = float(band_node_residual[band_index, final_node])
            diagnostics[f"surface_band_{band_name}_cost_final_node"] = final_cost
            diagnostics[f"surface_band_{band_name}_residual_final_node"] = final_residual
    if len(argmins) >= 2:
        diagnostics["surface_band_argmin_spread"] = float(max(argmins) - min(argmins))
    return diagnostics


def _prefix_surface_diagnostics(diagnostics: dict[str, Any], prefix: str) -> dict[str, Any]:
    out: dict[str, Any] = {}
    for key, value in diagnostics.items():
        new_key = str(key)
        if new_key.startswith("surface_"):
            new_key = f"{prefix}_{new_key[len('surface_') :]}"
        else:
            new_key = f"{prefix}_{new_key}"
        out[new_key] = value
    return out


def _cost_field_window_indices(
    cube: np.ndarray,
    *,
    x: np.ndarray,
    y: np.ndarray,
    center_x: float,
    center_y: float,
    radius_m: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    cube_arr = np.asarray(cube, dtype=np.float64)
    if cube_arr.ndim != 3:
        raise ValueError(f"cost cube must have shape (aot, y, x); got {cube_arr.shape!r}")
    _n_aot, ny, nx = cube_arr.shape
    x_arr = np.asarray(x, dtype=np.float64).reshape(-1)
    y_arr = np.asarray(y, dtype=np.float64).reshape(-1)
    if x_arr.size != nx or y_arr.size != ny:
        raise ValueError(
            "cost-field x/y coordinates do not match cube shape: "
            f"x={x_arr.size}, y={y_arr.size}, cube_yx={(ny, nx)!r}"
        )
    if radius_m <= 0.0:
        raise ValueError("cost-field window radius must be positive.")

    x_idx = np.flatnonzero(np.abs(x_arr - float(center_x)) <= float(radius_m))
    y_idx = np.flatnonzero(np.abs(y_arr - float(center_y)) <= float(radius_m))
    if x_idx.size == 0 or y_idx.size == 0:
        raise ValueError(
            "cost-field window does not intersect the solver grid "
            f"(center=({center_x:.3f}, {center_y:.3f}), radius_m={radius_m:.1f})."
        )
    return cube_arr, x_arr, y_arr, x_idx, y_idx


def _local_cost_diagnostics(
    cube: np.ndarray,
    *,
    aot_axis: np.ndarray,
    x: np.ndarray,
    y: np.ndarray,
    center_x: float,
    center_y: float,
    radius_m: float,
) -> dict[str, Any]:
    cube_arr, _x_arr, _y_arr, x_idx, y_idx = _cost_field_window_indices(
        cube,
        x=x,
        y=y,
        center_x=center_x,
        center_y=center_y,
        radius_m=radius_m,
    )
    cube_win = np.ascontiguousarray(cube_arr[:, y_idx][:, :, x_idx], dtype=np.float64)
    diagnostics = _prefix_surface_diagnostics(
        _cost_curve_diagnostics(cube_win, aot_axis),
        "local",
    )
    diagnostics.update(
        {
            "local_window_radius_m": float(radius_m),
            "local_window_shape_y": int(y_idx.size),
            "local_window_shape_x": int(x_idx.size),
        }
    )
    return diagnostics


def _local_band_diagnostics(
    *,
    band_names: list[str],
    aot_axis: np.ndarray,
    band_cost_cube: np.ndarray | None,
    band_residual_cube: np.ndarray | None,
    x: np.ndarray,
    y: np.ndarray,
    center_x: float,
    center_y: float,
    radius_m: float,
    final_aot: float,
) -> dict[str, Any]:
    if band_cost_cube is None:
        return {}
    band_cost = np.asarray(band_cost_cube, dtype=np.float64)
    if band_cost.ndim != 4 or band_cost.shape[0] != len(band_names):
        return {}
    first_band = band_cost[0]
    _cube_arr, _x_arr, _y_arr, x_idx, y_idx = _cost_field_window_indices(
        first_band,
        x=x,
        y=y,
        center_x=center_x,
        center_y=center_y,
        radius_m=radius_m,
    )
    band_residual = (
        np.asarray(band_residual_cube, dtype=np.float64)
        if band_residual_cube is not None
        and np.asarray(band_residual_cube).shape == band_cost.shape
        else None
    )
    final_node = int(np.nanargmin(np.abs(aot_axis - final_aot))) if np.isfinite(final_aot) else -1
    diagnostics: dict[str, Any] = {}
    if final_node >= 0:
        diagnostics["local_final_aot_node"] = float(aot_axis[final_node])

    argmins: list[float] = []
    for band_index, band_name in enumerate(band_names):
        cost_win = np.ascontiguousarray(
            band_cost[band_index][:, y_idx][:, :, x_idx],
            dtype=np.float64,
        )
        node_cost = _median_per_node(cost_win)
        finite = np.flatnonzero(np.isfinite(node_cost))
        if finite.size:
            best = int(finite[np.argmin(node_cost[finite])])
            argmin_aot = float(aot_axis[best])
            argmins.append(argmin_aot)
            diagnostics[f"local_band_{band_name}_argmin_aot"] = argmin_aot
            diagnostics[f"local_band_{band_name}_argmin_cost"] = float(node_cost[best])
        if final_node >= 0:
            diagnostics[f"local_band_{band_name}_cost_final_node"] = float(node_cost[final_node])
            if band_residual is not None:
                residual_win = np.ascontiguousarray(
                    band_residual[band_index][:, y_idx][:, :, x_idx],
                    dtype=np.float64,
                )
                node_residual = _median_per_node(residual_win)
                diagnostics[f"local_band_{band_name}_residual_final_node"] = float(
                    node_residual[final_node]
                )
    if len(argmins) >= 2:
        diagnostics["local_band_argmin_spread"] = float(max(argmins) - min(argmins))
    return diagnostics


def _scalar_int(value: np.ndarray, default: int) -> int:
    arr = np.asarray(value).reshape(-1)
    if arr.size == 0:
        return default
    try:
        return int(arr[0])
    except (TypeError, ValueError):
        return default


def _resolve_cost_field_window(
    cube: np.ndarray,
    *,
    aot_axis: np.ndarray,
    aot_prior: np.ndarray,
    aot_prior_unc: np.ndarray,
    solve_valid: np.ndarray,
    x: np.ndarray,
    y: np.ndarray,
    center_x: float,
    center_y: float,
    radius_m: float,
    pool_window: int,
    min_count: int,
) -> dict[str, object]:
    cube_arr = np.asarray(cube, dtype=np.float64)
    if cube_arr.ndim != 3:
        raise ValueError(f"cost cube must have shape (aot, y, x); got {cube_arr.shape!r}")
    _n_aot, ny, nx = cube_arr.shape
    x_arr = np.asarray(x, dtype=np.float64).reshape(-1)
    y_arr = np.asarray(y, dtype=np.float64).reshape(-1)
    if x_arr.size != nx or y_arr.size != ny:
        raise ValueError(
            "cost-field x/y coordinates do not match cube shape: "
            f"x={x_arr.size}, y={y_arr.size}, cube_yx={(ny, nx)!r}"
        )
    if radius_m <= 0.0:
        raise ValueError("cost-field window radius must be positive.")

    x_idx = np.flatnonzero(np.abs(x_arr - float(center_x)) <= float(radius_m))
    y_idx = np.flatnonzero(np.abs(y_arr - float(center_y)) <= float(radius_m))
    if x_idx.size == 0 or y_idx.size == 0:
        raise ValueError(
            "cost-field window does not intersect the solver grid "
            f"(center=({center_x:.3f}, {center_y:.3f}), radius_m={radius_m:.1f})."
        )

    cube_win = np.ascontiguousarray(cube_arr[:, y_idx][:, :, x_idx], dtype=np.float64)
    prior_win = np.ascontiguousarray(np.asarray(aot_prior, dtype=np.float64)[y_idx][:, x_idx])
    unc_win = np.ascontiguousarray(np.asarray(aot_prior_unc, dtype=np.float64)[y_idx][:, x_idx])
    valid_win = np.ascontiguousarray(np.asarray(solve_valid, dtype=bool)[y_idx][:, x_idx])

    aod, unc, cost = surface_driven_pool_argmin(
        cube_win,
        np.asarray(aot_axis, dtype=np.float64),
        prior_win,
        unc_win,
        valid_win,
        int(pool_window),
        int(min_count),
    )
    aod_arr = np.asarray(aod, dtype=np.float64)
    unc_arr = np.asarray(unc, dtype=np.float64)
    cost_arr = np.asarray(cost, dtype=np.float64)
    finite = np.isfinite(aod_arr)
    if not np.any(finite):
        raise ValueError("integrated cost-field solve produced no finite pixels in the window.")

    cost_finite = np.isfinite(cost_arr) & finite
    unc_finite = np.isfinite(unc_arr) & finite
    return {
        "aod": float(np.nanmedian(aod_arr[finite])),
        "aod_unc": float(np.nanmedian(unc_arr[unc_finite])) if np.any(unc_finite) else None,
        "cost": float(np.nanmedian(cost_arr[cost_finite])) if np.any(cost_finite) else None,
        "n_finite": int(np.count_nonzero(finite)),
        "n_window": int(finite.size),
        "n_valid_input": int(np.count_nonzero(valid_win)),
        "window_shape": [int(y_idx.size), int(x_idx.size)],
        "x_range": [float(np.nanmin(x_arr[x_idx])), float(np.nanmax(x_arr[x_idx]))],
        "y_range": [float(np.nanmin(y_arr[y_idx])), float(np.nanmax(y_arr[y_idx]))],
    }


def integrated_cost_field_aod(
    *,
    cube: np.ndarray,
    aot_axis: np.ndarray,
    aot_prior: np.ndarray,
    aot_prior_unc: np.ndarray,
    solve_valid: np.ndarray,
    x: np.ndarray,
    y: np.ndarray,
    center_x: float,
    center_y: float,
    radius_m: float = 1500.0,
    pool_window: int = 1,
    min_count: int = 1,
    cube_abs: np.ndarray | None = None,
    band_cost_cube: np.ndarray | None = None,
    band_residual_cube: np.ndarray | None = None,
    band_names: list[str] | None = None,
    clean_threshold: float = 0.15,
    high_threshold: float = 0.6,
) -> dict[str, object]:
    """Resolve a point retrieval by solving a cropped spatial cost field.

    The production raster solve returns a per-pixel AOT map and validation often
    samples the nearest pixel. The research harness instead crops the raw cost
    cube to the target window, runs the spatial median-pool/argmin there, and
    returns the median AOT over finite solved pixels. This helper implements
    that integrated site-window retrieval with the same Rust kernel used by the
    raster solver.
    """
    main = _resolve_cost_field_window(
        cube,
        aot_axis=aot_axis,
        aot_prior=aot_prior,
        aot_prior_unc=aot_prior_unc,
        solve_valid=solve_valid,
        x=x,
        y=y,
        center_x=center_x,
        center_y=center_y,
        radius_m=radius_m,
        pool_window=pool_window,
        min_count=min_count,
    )
    abs_result: dict[str, object] | None = None
    selected = main
    selected_pass = "main"
    mode = "single"

    if cube_abs is not None:
        cube_abs_arr = np.asarray(cube_abs)
        if cube_abs_arr.shape == np.asarray(cube).shape:
            abs_result = _resolve_cost_field_window(
                cube_abs_arr,
                aot_axis=aot_axis,
                aot_prior=aot_prior,
                aot_prior_unc=aot_prior_unc,
                solve_valid=solve_valid,
                x=x,
                y=y,
                center_x=center_x,
                center_y=center_y,
                radius_m=radius_m,
                pool_window=pool_window,
                min_count=min_count,
            )
            mode = "auto2"
            selected_pass = "shape"
            abs_aod = float(abs_result["aod"])
            if not (float(clean_threshold) <= abs_aod <= float(high_threshold)):
                selected = abs_result
                selected_pass = "abs"

    selected_aod = float(selected["aod"])
    local_diagnostics = _local_cost_diagnostics(
        cube,
        aot_axis=aot_axis,
        x=x,
        y=y,
        center_x=center_x,
        center_y=center_y,
        radius_m=radius_m,
    )
    if band_names is not None:
        local_diagnostics.update(
            _local_band_diagnostics(
                band_names=[str(name) for name in band_names],
                aot_axis=aot_axis,
                band_cost_cube=band_cost_cube,
                band_residual_cube=band_residual_cube,
                x=x,
                y=y,
                center_x=center_x,
                center_y=center_y,
                radius_m=radius_m,
                final_aot=selected_aod,
            )
        )

    return {
        "mode": mode,
        "selected_pass": selected_pass,
        "aod": selected_aod,
        "aod_unc": selected.get("aod_unc"),
        "cost": selected.get("cost"),
        "radius_m": float(radius_m),
        "pool_window": int(pool_window),
        "min_count": int(min_count),
        "clean_threshold": float(clean_threshold),
        "high_threshold": float(high_threshold),
        "main": main,
        "abs": abs_result,
        "diagnostics": local_diagnostics,
    }


def integrated_cost_field_aod_from_npz(
    path: str,
    *,
    center_x: float,
    center_y: float,
    radius_m: float = 1500.0,
    clean_threshold: float = 0.15,
    high_threshold: float = 0.6,
) -> dict[str, object]:
    """Load a ``SIAC_DUMP_COST_CUBE`` archive and resolve its target window."""
    with np.load(path, allow_pickle=False) as z:
        cube_abs = z["cube_abs"] if "cube_abs" in z.files and z["cube_abs"].size else None
        band_cost_cube = (
            z["band_cost_cube"]
            if "band_cost_cube" in z.files and z["band_cost_cube"].size
            else None
        )
        band_residual_cube = (
            z["band_residual_cube"]
            if "band_residual_cube" in z.files and z["band_residual_cube"].size
            else None
        )
        band_names = (
            [str(name) for name in z["band_names"].astype(str).tolist()]
            if "band_names" in z.files
            else None
        )
        return integrated_cost_field_aod(
            cube=z["cube"],
            cube_abs=cube_abs,
            band_cost_cube=band_cost_cube,
            band_residual_cube=band_residual_cube,
            band_names=band_names,
            aot_axis=z["aot_axis"],
            aot_prior=z["aot_prior"],
            aot_prior_unc=z["aot_prior_unc"],
            solve_valid=z["solve_valid"],
            x=z["x"],
            y=z["y"],
            center_x=center_x,
            center_y=center_y,
            radius_m=radius_m,
            pool_window=_scalar_int(
                z["pool_window"] if "pool_window" in z.files else np.array([]), 1
            ),
            min_count=_scalar_int(z["min_count"] if "min_count" in z.files else np.array([]), 1),
            clean_threshold=clean_threshold,
            high_threshold=high_threshold,
        )


class SurfaceDrivenSolver:
    """Surface-prior-driven AOT solver over a fixed AOT axis."""

    def __init__(self, config: SolverAlgorithmConfig):
        self.config = config

    # -- backstop ----------------------------------------------------------
    def _backstop_unc(self, aot_prior: np.ndarray, aot_prior_unc_raw: np.ndarray) -> np.ndarray:
        """Per-pixel AOT-prior (CAMS) backstop sigma.

        Calibrated mode (default) widens the prior when clean and at the
        high-AOD tail and tightens it through the moderate band where the
        surface signal is shallow; otherwise a flat 50 % relative prior.
        """
        if not self.config.surface_driven_backstop_calibrated:
            unc = np.maximum(0.5 * aot_prior, 0.02)
        else:
            m = np.asarray(aot_prior, dtype=np.float64)
            loose = np.maximum(0.5 * m, 0.02)
            mid = np.maximum(0.07, 0.5 * m / (1.0 + np.exp(-(m - 0.5) / 0.15)))
            unc = np.where(m < 0.15, loose, mid)
        unc = np.where(np.isfinite(unc), unc, np.maximum(aot_prior_unc_raw, 0.02))
        unc = unc * float(getattr(self.config, "surface_driven_backstop_uncertainty_scale", 1.0))
        return cast("np.ndarray", np.asarray(unc, dtype=np.float32))

    # -- solve -------------------------------------------------------------
    def solve(
        self,
        toa: xr.DataArray,
        surface_prior: SurfacePrior,
        geometry: GeometryAngles,
        atmo_prior: AtmosphericState,
        rt_model: RTModelBackend,
        cloud_mask: xr.DataArray,
        bands: list[SensorBand],
        sharp_transition_mask: xr.DataArray | None = None,
        water_mask: xr.DataArray | None = None,
        _force_tau_prior: bool | None = None,
    ) -> SolverResult:
        # sharp_transition_mask is accepted for solver-interface parity but the
        # surface-driven sweep does not use it (no cross-edge smoothness term).
        del sharp_transition_mask
        ny = int(cloud_mask.sizes["y"])
        nx = int(cloud_mask.sizes["x"])
        shape = (ny, nx)
        n_bands = len(bands)

        # Cloud and water masks are enforced by default. Thick aerosol can be
        # misclassified as cloud, so an explicit cloud-only rescue switch is
        # available; it still masks water and requires complete finite support.
        # The older combined switch remains for backward-compatible experiments.
        ignore_cloud_water = bool(getattr(self.config, "surface_driven_ignore_cloud_water", False))
        allow_cloud_retrieval = ignore_cloud_water or bool(
            getattr(self.config, "surface_driven_allow_cloud_retrieval", False)
        )
        valid_mask: np.ndarray = (
            np.ones(shape, dtype=bool) if allow_cloud_retrieval else ~cloud_mask.values.astype(bool)
        )
        if water_mask is not None and not ignore_cloud_water:
            valid_mask = valid_mask & ~water_mask.values.astype(bool)

        aot_prior = atmo_prior.aot.values.astype(np.float32)
        tcwv_prior = atmo_prior.tcwv.values.astype(np.float32)
        aot_prior_unc_raw = atmo_prior.aot_unc.values.astype(np.float32)
        tcwv_prior_unc = np.maximum(atmo_prior.tcwv_unc.values.astype(np.float32), 1e-3)
        finite_tcwv_prior = tcwv_prior[np.isfinite(tcwv_prior)]
        tcwv_prior_is_spatial = bool(
            finite_tcwv_prior.size > 1
            and np.max(finite_tcwv_prior) - np.min(finite_tcwv_prior)
            > 1e-6 * max(1.0, float(np.max(np.abs(finite_tcwv_prior))))
        )

        # AOT axis: log-spaced over the configured bounds (default) or the
        # validated fine non-uniform "acixthree" axis (denser at low AOD).
        # TCWV remains fixed while this axis is swept: normally at its
        # per-pixel prior field, or at a configured scalar reference.
        if str(getattr(self.config, "surface_driven_aot_axis", "log")) == "acixthree":
            aot_axis = _acixthree_aot_axis()
        else:
            n_aot = max(3, int(self.config.grid_search_aot_points))
            aot_axis = _log_aot_axis(self.config.aot_bounds[0], self.config.aot_bounds[1], n_aot)
        n_aot = int(aot_axis.size)
        # A scalar reference (the validated harness ran the cost cube at
        # TCWV=2.0, sea-level) overrides the full prior field. With ``None``,
        # direct RT evaluation keeps the per-pixel field fixed and the median
        # below is only the scalar joint-LUT axis placeholder.
        reference_tcwv = getattr(self.config, "surface_driven_reference_tcwv", None)
        if reference_tcwv is not None:
            tcwv_fixed = float(reference_tcwv)
        else:
            tcwv_fixed = (
                float(np.median(finite_tcwv_prior))
                if finite_tcwv_prior.size
                else float(np.mean(self.config.tcwv_bounds))
            )

        (
            toa_values,
            _no_observation_mask,
            boa_prior,
            boa_unc,
            support_mask,
        ) = prepare_grid_search_observations(
            toa=toa,
            surface_prior=surface_prior,
            n_bands=n_bands,
            shape=shape,
            min_boa_unc=_MIN_BOA_UNC,
            aot_prior=aot_prior,
            tcwv_prior=tcwv_prior,
            aot_prior_unc=np.maximum(aot_prior_unc_raw, 1e-3),
            tcwv_prior_unc=tcwv_prior_unc,
        )
        solve_valid = valid_mask & support_mask

        # Coefficient provider with TCWV fixed while AOT varies. When a
        # reference is configured, replace the full field so every RT path uses
        # that scalar. Otherwise the direct path preserves the per-pixel prior;
        # the joint-LUT shortcut detects and declines non-uniform fixed fields.
        coeff_atmo_prior = atmo_prior
        if reference_tcwv is not None:
            coeff_atmo_prior = replace(
                atmo_prior,
                tcwv=xr.full_like(atmo_prior.tcwv, np.float32(tcwv_fixed)),
            )
        # Optionally evaluate the cost-cube RT at the single scene-mean geometry
        # (nanmean of each radian angle field), matching the harness. Per-pixel
        # view azimuth varies across the swath and shifts the shape-cube AOT
        # minimum; collapsing to the AOI mean removes that spread.
        coeff_geometry = geometry
        if getattr(self.config, "surface_driven_scene_mean_geometry", False):

            def _mean_field(da: xr.DataArray) -> xr.DataArray:
                vals = np.asarray(da.values, dtype=np.float64)
                finite = vals[np.isfinite(vals)]
                fill = float(np.mean(finite)) if finite.size else 0.0
                return xr.full_like(da, np.float32(fill))

            coeff_geometry = GeometryAngles(
                sza=_mean_field(geometry.sza),
                saa=_mean_field(geometry.saa),
                vza=_mean_field(geometry.vza),
                vaa=_mean_field(geometry.vaa),
            )

        # Surface-cost cube(s). f64 and the 1.0 (no-0.5) chi-square convention so
        # the kernel's backstop ``d^2/unc^2`` is on the same scale as the obs
        # cost — a bit-exact match to the reference numpy resolve.
        #
        # cost_mode:
        #   "chi2"  -> absolute per-band chi-square (lenient: partial-band sum;
        #              unchanged default behaviour).
        #   "shape" -> spectral-shape ratio cost (brightness-bias robust).
        #   "auto2" -> validated regime scheme: shape (all bands) + absolute
        #              (non-B01) cubes; classify by the absolute-solved AOD and
        #              use the shape solve only in the moderate band.
        #   "profile_scale" -> one constrained common surface-brightness scale.
        #   "loo_scale" -> leave-one-band-out prediction after fitting that scale
        #                  on the remaining visible bands.
        #   "additive_offset" -> profile one spectrally flat BOA offset and
        #                        retain only wavelength-dependent mismatch.
        #   "trimmed_chi2" -> absolute chi-square after dropping one inconsistent
        #                     visible band at each pixel/AOT node.
        cost_mode = str(getattr(self.config, "surface_driven_cost_mode", "chi2"))
        band_names = [str(getattr(b, "name", b)) for b in bands]
        anchor_index = band_names.index("B04") if "B04" in band_names else n_bands - 1
        abs_indices = [i for i, nm in enumerate(band_names) if nm != "B01"] or list(range(n_bands))
        all_indices = list(range(n_bands))
        two_pass = cost_mode == "auto2"
        use_shape_main = cost_mode in ("shape", "auto2")
        profile_scale_sigma = float(getattr(self.config, "surface_driven_profile_scale_sigma", 0.1))
        if profile_scale_sigma <= 0.0:
            raise ValueError("surface_driven_profile_scale_sigma must be positive")
        if cost_mode == "loo_scale" and len(all_indices) < 3:
            raise ValueError('surface_driven_cost_mode="loo_scale" requires at least 3 bands')
        if cost_mode == "trimmed_chi2" and len(all_indices) < 3:
            raise ValueError('surface_driven_cost_mode="trimmed_chi2" requires at least 3 bands')

        inv_var = 1.0 / np.square(np.maximum(boa_unc, 1e-6))  # (band, y, x)

        # -- tau-dependent surface prior (joint surface/AOD self-consistency) --
        # When enabled, the visible prior is re-predicted per AOT node from the
        # scene NIR/SWIR anchors corrected at THAT node's AOD, so the surface
        # term carries real information at high AOD where a fixed-anchor
        # prediction goes off the dictionary manifold.
        tau_payload = getattr(surface_prior, "tau_predictor", None)
        tau_config_on = (
            _force_tau_prior
            if _force_tau_prior is not None
            else bool(getattr(self.config, "surface_driven_tau_dependent_prior", False))
        )
        tau_enabled = (
            tau_config_on
            and isinstance(tau_payload, dict)
            and all(
                key in tau_payload
                for key in (
                    "trees",
                    "localizer_grid",
                    "anchor_toa_grid",
                    "anchor_sensor_bands",
                    "target_bands",
                )
            )
        )
        tau_target_idx: list[int] = []
        if tau_enabled:
            tau_target_idx = [
                band_names.index(name) for name in tau_payload["target_bands"] if name in band_names
            ]
            tau_enabled = bool(tau_target_idx)
        if tau_enabled:
            tau_trees = tau_payload["trees"]
            tau_loc = np.asarray(tau_payload["localizer_grid"], dtype=np.float64)
            tau_loc_flat = tau_loc.reshape(tau_loc.shape[0], -1).T  # (px, n_loc)
            tau_anchor_toa = np.asarray(tau_payload["anchor_toa_grid"], dtype=np.float64)
            tau_loc_valid = np.all(np.isfinite(tau_loc_flat), axis=1)
            tau_debias = dict(tau_payload.get("debias") or {})
            tau_debias_scale = float(tau_payload.get("debias_scale", 1.0))
            logger.info(
                "surface_driven: tau-dependent prior ENABLED (%d trees, targets %s)",
                len(tau_trees),
                ",".join(tau_payload["target_bands"]),
            )

        def _tau_prior_node(anchor_boa: np.ndarray, tau_value: float) -> np.ndarray:
            """Predict the visible prior for one AOT node from tau-corrected anchors."""
            flat = anchor_boa.reshape(anchor_boa.shape[0], -1).T  # (px, n_anchor)
            valid = tau_loc_valid & np.all(np.isfinite(flat) & (flat > 0.0) & (flat < 1.2), axis=1)
            prior_k = boa_prior.copy()
            if int(np.count_nonzero(valid)) < 50:
                return prior_k
            feats = np.column_stack([flat[valid], tau_loc_flat[valid]])
            preds = np.median(np.stack([tree.predict(feats) for tree in tau_trees]), axis=0)
            if preds.ndim == 1:
                preds = preds[:, np.newaxis]
            for out_index, name in enumerate(tau_payload["target_bands"]):
                if name not in band_names:
                    continue
                band_index = band_names.index(name)
                intercept, slope = tau_debias.get(name, (0.0, 0.0))
                values = preds[:, out_index] + tau_debias_scale * (
                    float(intercept) + float(slope) * float(tau_value)
                )
                plane = prior_k[band_index].reshape(-1)
                ok = np.isfinite(values) & (values > 0.001) & (values < 0.6)
                idx = np.flatnonzero(valid)[ok]
                plane[idx] = values[ok]
                prior_k[band_index] = plane.reshape(prior_k.shape[1:])
            return prior_k

        # Backstop centre fed to the pool-argmin is the raw aerosol prior. It
        # regularizes valid retrievals but is never used to fill invalid pixels.
        aot_prior_backstop = aot_prior
        aot_prior_unc = self._backstop_unc(aot_prior_backstop, aot_prior_unc_raw)
        # Full rolling-window edge length (pandas/xarray center convention) from
        # the ~1.2 km radius; min_periods matches the reference resolve.
        pool_window = max(
            1,
            int(
                round(
                    2.0
                    * self.config.surface_driven_pool_radius_m
                    / max(self.config.aerosol_resolution, 1e-3)
                )
            ),
        )
        min_count = (
            max(int(self.config.surface_driven_pool_min_count), 20, pool_window * pool_window // 5)
            if pool_window > 1
            else 1
        )

        def _resolve(
            cube_in: np.ndarray,
            *,
            backstop_unc: np.ndarray | None = None,
        ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
            prior_unc = aot_prior_unc if backstop_unc is None else backstop_unc
            a, u, c = surface_driven_pool_argmin(
                np.ascontiguousarray(cube_in, dtype=np.float64),
                aot_axis.astype(np.float64, copy=False),
                np.ascontiguousarray(aot_prior_backstop, dtype=np.float64),
                np.ascontiguousarray(prior_unc, dtype=np.float64),
                np.ascontiguousarray(solve_valid, dtype=bool),
                int(pool_window),
                int(min_count),
            )
            return (
                np.asarray(a, dtype=np.float32),
                np.asarray(u, dtype=np.float32),
                np.asarray(c, dtype=np.float32),
            )

        def _build_cost_cubes(
            candidate_rt_model: RTModelBackend,
        ) -> tuple[
            np.ndarray,
            np.ndarray | None,
            dict[str, Any],
            np.ndarray,
            np.ndarray,
            np.ndarray | None,
            np.ndarray | None,
            np.ndarray | None,
        ]:
            coeff_provider = build_candidate_coeff_provider(
                rt_model=candidate_rt_model,
                bands=bands,
                geometry=coeff_geometry,
                atmo_prior=coeff_atmo_prior,
                aot_axis=aot_axis,
                tcwv_axis=np.asarray([tcwv_fixed], dtype=np.float32),
                solve_aot=True,
                solve_tcwv=False,
                aot_bounds=self.config.aot_bounds,
                tcwv_bounds=self.config.tcwv_bounds,
                rt_sample_step=1,
                shape=shape,
            )
            anchor_coeff_provider = None
            if tau_enabled:
                anchor_coeff_provider = build_candidate_coeff_provider(
                    rt_model=candidate_rt_model,
                    bands=list(tau_payload["anchor_sensor_bands"]),
                    geometry=coeff_geometry,
                    atmo_prior=coeff_atmo_prior,
                    aot_axis=aot_axis,
                    tcwv_axis=np.asarray([tcwv_fixed], dtype=np.float32),
                    solve_aot=True,
                    solve_tcwv=False,
                    aot_bounds=self.config.aot_bounds,
                    tcwv_bounds=self.config.tcwv_bounds,
                    rt_sample_step=1,
                    shape=shape,
                )
            cube: np.ndarray = np.full((n_aot, ny, nx), np.nan, dtype=np.float64)
            cube_abs: np.ndarray | None = (
                np.full((n_aot, ny, nx), np.nan, dtype=np.float64) if two_pass else None
            )
            band_node_cost = np.full((n_bands, n_aot), np.nan, dtype=np.float64)
            band_node_residual = np.full((n_bands, n_aot), np.nan, dtype=np.float64)
            dump_band_cubes = bool(os.environ.get("SIAC_DUMP_COST_CUBE"))
            band_cost_cube: np.ndarray | None = (
                np.full((n_bands, n_aot, ny, nx), np.nan, dtype=np.float32)
                if dump_band_cubes
                else None
            )
            band_residual_cube: np.ndarray | None = (
                np.full((n_bands, n_aot, ny, nx), np.nan, dtype=np.float32)
                if dump_band_cubes
                else None
            )
            band_signed_residual_cube: np.ndarray | None = (
                np.full((n_bands, n_aot, ny, nx), np.nan, dtype=np.float32)
                if dump_band_cubes
                else None
            )
            for k in range(n_aot):
                xap, xbp, xcp = coeff_provider(float(aot_axis[k]), tcwv_fixed)
                y = np.asarray(xap, np.float64) * toa_values - np.asarray(xbp, np.float64)
                denom = 1.0 + np.asarray(xcp, np.float64) * y
                with np.errstate(invalid="ignore", divide="ignore"):
                    boa = np.where(np.abs(denom) > 1e-12, y / denom, np.nan)
                prior_node = boa_prior
                if anchor_coeff_provider is not None:
                    xa, xb, xc = anchor_coeff_provider(float(aot_axis[k]), tcwv_fixed)
                    ya = np.asarray(xa, np.float64) * tau_anchor_toa - np.asarray(xb, np.float64)
                    da = 1.0 + np.asarray(xc, np.float64) * ya
                    with np.errstate(invalid="ignore", divide="ignore"):
                        anchor_boa = np.where(np.abs(da) > 1e-12, ya / da, np.nan)
                    prior_node = _tau_prior_node(anchor_boa, float(aot_axis[k]))
                diff = boa - prior_node
                abs_term = inv_var * diff * diff  # (band, y, x)
                diag_valid = solve_valid
                for band_index in range(n_bands):
                    term_band = abs_term[band_index]
                    signed_resid_band = diff[band_index]
                    resid_band = np.abs(signed_resid_band)
                    if band_cost_cube is not None:
                        band_cost_cube[band_index, k] = np.where(
                            solve_valid, term_band, np.nan
                        ).astype(np.float32)
                    if band_residual_cube is not None:
                        band_residual_cube[band_index, k] = np.where(
                            solve_valid, resid_band, np.nan
                        ).astype(np.float32)
                    if band_signed_residual_cube is not None:
                        band_signed_residual_cube[band_index, k] = np.where(
                            solve_valid, signed_resid_band, np.nan
                        ).astype(np.float32)
                    valid_band = diag_valid & np.isfinite(term_band)
                    if np.any(valid_band):
                        band_node_cost[band_index, k] = float(np.median(term_band[valid_band]))
                    valid_resid = diag_valid & np.isfinite(resid_band)
                    if np.any(valid_resid):
                        band_node_residual[band_index, k] = float(
                            np.median(resid_band[valid_resid])
                        )
                if cost_mode == "profile_scale":
                    cube[k] = np.where(
                        solve_valid,
                        _profile_scale_cost_node(
                            boa,
                            prior_node,
                            boa_unc,
                            all_indices,
                            scale_sigma=profile_scale_sigma,
                        ),
                        np.nan,
                    )
                elif cost_mode == "additive_offset":
                    cube[k] = np.where(
                        solve_valid,
                        _additive_offset_cost_node(
                            boa,
                            prior_node,
                            boa_unc,
                            all_indices,
                        ),
                        np.nan,
                    )
                elif cost_mode == "loo_scale":
                    cube[k] = np.where(
                        solve_valid,
                        _loo_scale_cost_node(
                            boa,
                            prior_node,
                            boa_unc,
                            all_indices,
                            scale_sigma=profile_scale_sigma,
                        ),
                        np.nan,
                    )
                elif cost_mode == "trimmed_chi2":
                    cube[k] = np.where(
                        solve_valid,
                        _trimmed_chi2_cost_node(
                            boa,
                            prior_node,
                            inv_var,
                            all_indices,
                        ),
                        np.nan,
                    )
                elif use_shape_main:
                    cube[k] = np.where(
                        solve_valid,
                        _shape_cost_node(boa, prior_node, boa_unc, all_indices, anchor_index),
                        np.nan,
                    )
                else:
                    finite = np.isfinite(abs_term)
                    used = np.count_nonzero(finite, axis=0)
                    cost = np.where(finite, abs_term, 0.0).sum(axis=0)
                    cube[k] = np.where(solve_valid & (used > 0), cost, np.nan)
                if cube_abs is not None:
                    cube_abs[k] = np.where(
                        solve_valid,
                        _abs_cost_node(boa, prior_node, inv_var, abs_indices),
                        np.nan,
                    )
            diagnostics = {
                **_cost_curve_diagnostics(cube, aot_axis),
                "surface_cost_mode": cost_mode,
                "surface_tau_prior_enabled": bool(tau_enabled),
                "surface_tcwv_cost_mode": (
                    "scalar_reference"
                    if reference_tcwv is not None
                    else (
                        "fixed_spatial_prior"
                        if tcwv_prior_is_spatial
                        else "fixed_scalar_prior"
                    )
                ),
                "surface_tcwv_reference_cm": (
                    float(reference_tcwv) if reference_tcwv is not None else None
                ),
                "surface_aot_axis_min": float(np.nanmin(aot_axis)),
                "surface_aot_axis_max": float(np.nanmax(aot_axis)),
                "surface_aot_axis_nodes": int(aot_axis.size),
                "surface_solve_band_count": int(n_bands),
            }
            return (
                cube,
                cube_abs,
                diagnostics,
                band_node_cost,
                band_node_residual,
                band_cost_cube,
                band_residual_cube,
                band_signed_residual_cube,
            )

        def _dump_cost_cubes(
            cube: np.ndarray,
            cube_abs: np.ndarray | None,
            band_cost_cube: np.ndarray | None,
            band_residual_cube: np.ndarray | None,
            band_signed_residual_cube: np.ndarray | None,
            *,
            species_index: int | None,
        ) -> None:
            # Optional diagnostic dump of the raw per-node cost cube(s) + grid,
            # for offline analysis (env-gated; no-op by default). When species
            # candidates are active, each candidate gets its own suffixed file.
            dump_path_raw = os.environ.get("SIAC_DUMP_COST_CUBE")
            if not dump_path_raw:
                return
            dump_path = Path(dump_path_raw)
            if species_index is not None:
                suffix = dump_path.suffix or ".npz"
                dump_path = dump_path.with_name(
                    f"{dump_path.stem}.species{species_index}{suffix}"
                )
            template = atmo_prior.aot
            x_coord = template.coords.get("x")
            y_coord = template.coords.get("y")
            temporary_path = dump_path.with_name(f".{dump_path.name}.tmp.{os.getpid()}.npz")
            try:
                np.savez(
                    temporary_path,
                    cube=cube.astype(np.float32),
                    cube_abs=(
                        cube_abs.astype(np.float32) if cube_abs is not None else np.zeros(0)
                    ),
                    band_cost_cube=(
                        band_cost_cube.astype(np.float32)
                        if band_cost_cube is not None
                        else np.zeros(0, dtype=np.float32)
                    ),
                    band_residual_cube=(
                        band_residual_cube.astype(np.float32)
                        if band_residual_cube is not None
                        else np.zeros(0, dtype=np.float32)
                    ),
                    band_signed_residual_cube=(
                        band_signed_residual_cube.astype(np.float32)
                        if band_signed_residual_cube is not None
                        else np.zeros(0, dtype=np.float32)
                    ),
                    aot_axis=aot_axis.astype(np.float32),
                    aot_prior=aot_prior_backstop.astype(np.float32),
                    aot_prior_unc=aot_prior_unc.astype(np.float32),
                    solve_valid=solve_valid.astype(bool),
                    boa_prior=boa_prior.astype(np.float32),
                    boa_unc=boa_unc.astype(np.float32),
                    toa=toa_values.astype(np.float32),
                    x=(
                        np.asarray(x_coord.values, dtype=np.float64)
                        if x_coord is not None
                        else np.zeros(0)
                    ),
                    y=(
                        np.asarray(y_coord.values, dtype=np.float64)
                        if y_coord is not None
                        else np.zeros(0)
                    ),
                    band_names=np.asarray(band_names),
                    pool_window=np.asarray([pool_window]),
                    min_count=np.asarray([min_count]),
                    species_index=np.asarray([-1 if species_index is None else species_index]),
                )
                temporary_path.replace(dump_path)
            finally:
                temporary_path.unlink(missing_ok=True)

        def _resolve_candidate(
            candidate_rt_model: RTModelBackend,
            *,
            species_index: int | None,
        ) -> tuple[np.ndarray, np.ndarray, np.ndarray, dict[str, Any]]:
            escape_enabled = bool(
                getattr(self.config, "surface_driven_backstop_escape_enabled", False)
            )
            escape_low_aot = float(
                getattr(self.config, "surface_driven_backstop_escape_low_aot", 0.08)
            )
            escape_cost_threshold = float(
                getattr(self.config, "surface_driven_backstop_escape_cost_threshold", 20.0)
            )
            escape_delta_threshold = float(
                getattr(self.config, "surface_driven_backstop_escape_delta", 0.05)
            )
            escape_band_spread_threshold = float(
                getattr(
                    self.config,
                    "surface_driven_backstop_escape_band_spread",
                    0.25,
                )
            )
            escape_min_jump = float(
                getattr(self.config, "surface_driven_backstop_escape_min_jump", 0.35)
            )
            escape_cost_ratio = float(
                getattr(self.config, "surface_driven_backstop_escape_cost_ratio", 1.8)
            )
            auto2_abs_cost_gain = float(
                getattr(self.config, "surface_driven_auto2_abs_cost_gain", 0.2)
            )
            if not (0.0 <= auto2_abs_cost_gain <= 1.0):
                raise ValueError("surface_driven_auto2_abs_cost_gain must be in [0.0, 1.0].")

            (
                cube,
                cube_abs,
                diagnostics,
                band_node_cost,
                band_node_residual,
                band_cost_cube,
                band_residual_cube,
                band_signed_residual_cube,
            ) = _build_cost_cubes(candidate_rt_model)
            _dump_cost_cubes(
                cube,
                cube_abs,
                band_cost_cube,
                band_residual_cube,
                band_signed_residual_cube,
                species_index=species_index,
            )
            def _resolve_cost_mode(
                backstop_unc: np.ndarray | None = None,
            ) -> tuple[np.ndarray, np.ndarray, np.ndarray, dict[str, Any]]:
                aot_main, unc_main, cost_main = _resolve(cube, backstop_unc=backstop_unc)
                if not two_pass or cube_abs is None:
                    return aot_main, unc_main, cost_main, {}
                aot_abs, unc_abs, cost_abs = _resolve(cube_abs, backstop_unc=backstop_unc)
                # Regime by the absolute-solved AOD: use the shape solve only in the
                # moderate band; absolute (non-B01) at clean and high tails, but only
                # when the tail branch is materially better.
                clean_thr = float(getattr(self.config, "surface_driven_aod_clean", 0.15))
                high_thr = float(getattr(self.config, "surface_driven_aod_high", 0.6))
                return _select_auto2_solution(
                    aot_main=aot_main,
                    unc_main=unc_main,
                    cost_main=cost_main,
                    aot_abs=aot_abs,
                    unc_abs=unc_abs,
                    cost_abs=cost_abs,
                    clean_threshold=clean_thr,
                    high_threshold=high_thr,
                    cost_gain=auto2_abs_cost_gain,
                )

            aot_arr, aot_unc_arr, cost_arr, resolve_diagnostics = _resolve_cost_mode()
            diagnostics.update(resolve_diagnostics)
            conflict_z = getattr(self.config, "surface_driven_backstop_conflict_z", None)
            use_flat, conflict_diagnostics = _backstop_conflict_decision(
                surface_min_aot=diagnostics.get("surface_cost_curve_min_aot"),
                aot_prior=aot_prior_backstop,
                calibrated_unc=aot_prior_unc,
                z_threshold=conflict_z,
            )
            if use_flat and bool(self.config.surface_driven_backstop_calibrated):
                flat_unc = np.maximum(0.5 * aot_prior_backstop, 0.02)
                flat_unc = flat_unc * float(
                    getattr(self.config, "surface_driven_backstop_uncertainty_scale", 1.0)
                )
                aot_arr, aot_unc_arr, cost_arr, resolve_diagnostics = _resolve_cost_mode(
                    np.asarray(flat_unc, dtype=np.float32)
                )
                diagnostics.update(resolve_diagnostics)
                conflict_diagnostics["surface_backstop_conflict_relaxed"] = True
            else:
                conflict_diagnostics["surface_backstop_conflict_relaxed"] = False
            diagnostics.update(conflict_diagnostics)
            finite = np.isfinite(cost_arr)
            base_cost = float(np.nanmedian(cost_arr[finite])) if np.any(finite) else float("nan")
            finite_aot = aot_arr[np.isfinite(aot_arr)]
            final_aot = float(np.median(finite_aot)) if finite_aot.size else float("nan")

            # ------------------------------------------------------------------ #
            # Optional low-aod backstop escape:
            # ------------------------------------------------------------------ #
            if escape_enabled:
                try:
                    curve_delta = diagnostics.get("surface_cost_curve_relative_second_delta")
                    curve_delta_finite = (
                        float(curve_delta)
                        if np.isfinite(np.asarray(curve_delta, dtype=np.float64))
                        else float("nan")
                    )
                except TypeError:
                    curve_delta_finite = float("nan")
                band_spread = diagnostics.get("surface_band_argmin_spread")
                band_spread_finite = (
                    float(band_spread)
                    if np.isfinite(np.asarray(band_spread, dtype=np.float64))
                    else float("nan")
                )

                curve_min_at_edge = bool(diagnostics.get("surface_cost_curve_min_at_edge", False))
                suspicious_curve = bool(
                    curve_min_at_edge
                    or (
                        np.isfinite(curve_delta_finite)
                        and curve_delta_finite <= escape_delta_threshold
                    )
                    or (
                        np.isfinite(band_spread_finite)
                        and band_spread_finite <= escape_band_spread_threshold
                    )
                )
                if suspicious_curve and np.isfinite(final_aot) and final_aot <= escape_low_aot:
                    no_backstop_unc = np.full_like(aot_prior_unc, np.inf, dtype=np.float64)
                    aot_escape, aot_escape_unc, cost_escape = _resolve(
                        cube, backstop_unc=no_backstop_unc
                    )
                    finite_escape = np.isfinite(aot_escape)
                    escape_aot = (
                        float(np.nanmedian(aot_escape[finite_escape]))
                        if np.any(finite_escape)
                        else float("nan")
                    )
                    finite_escape_cost = np.isfinite(cost_escape)
                    escape_cost = (
                        float(np.nanmedian(cost_escape[finite_escape_cost]))
                        if np.any(finite_escape_cost)
                        else float("nan")
                    )
                    diagnostics.update(
                        {
                            "surface_backstop_escape_tested": True,
                            "surface_backstop_escape_candidate_aot_median": escape_aot,
                            "surface_backstop_escape_candidate_cost": escape_cost,
                            "surface_backstop_escape_final_delta": (
                                escape_aot - final_aot if np.isfinite(escape_aot) else float("nan")
                            ),
                        }
                    )
                    if (
                        np.isfinite(escape_cost)
                        and np.isfinite(escape_aot)
                        and np.isfinite(base_cost)
                        and escape_cost <= base_cost * escape_cost_ratio
                        and escape_aot >= final_aot + escape_min_jump
                        and base_cost >= escape_cost_threshold
                    ):
                        replace = (
                            finite_escape
                            & np.isfinite(aot_arr)
                            & (aot_arr <= final_aot + escape_min_jump)
                            & (aot_escape > aot_arr)
                            & ((aot_escape - aot_arr) >= escape_min_jump)
                        )
                        replaced = int(np.count_nonzero(replace))
                        if replaced > 0:
                            aot_arr = np.where(replace, aot_escape, aot_arr)
                            aot_unc_arr = np.where(replace, aot_escape_unc, aot_unc_arr)
                            cost_arr = np.where(replace, cost_escape, cost_arr)
                            finite_after = np.isfinite(aot_arr)
                            final_aot = (
                                float(np.median(aot_arr[finite_after]))
                                if np.any(finite_after)
                                else float("nan")
                            )
                            finite_cost_after = np.isfinite(cost_arr)
                            base_cost = (
                                float(np.nanmedian(cost_arr[finite_cost_after]))
                                if np.any(finite_cost_after)
                                else float("nan")
                            )
                            diagnostics["surface_backstop_escape_applied"] = True
                            diagnostics["surface_backstop_escape_replaced_pixels"] = int(replaced)
                        else:
                            diagnostics["surface_backstop_escape_applied"] = False
                            diagnostics["surface_backstop_escape_replaced_pixels"] = int(replaced)
                    else:
                        diagnostics["surface_backstop_escape_applied"] = False
                        diagnostics["surface_backstop_escape_replaced_pixels"] = 0
                else:
                    diagnostics["surface_backstop_escape_tested"] = False
                    diagnostics["surface_backstop_escape_applied"] = False
                    diagnostics["surface_backstop_escape_replaced_pixels"] = 0
                diagnostics["surface_backstop_escape_base_aot"] = final_aot
                diagnostics["surface_backstop_escape_base_cost"] = base_cost
            else:
                diagnostics["surface_backstop_escape_tested"] = False
                diagnostics["surface_backstop_escape_applied"] = False
                diagnostics["surface_backstop_escape_replaced_pixels"] = 0

            diagnostics["surface_final_aot_median"] = final_aot
            diagnostics["surface_rt_branch"] = (
                f"species_{int(species_index)}" if species_index is not None else "default"
            )
            diagnostics.update(
                _band_curve_diagnostics(
                    band_names=band_names,
                    aot_axis=aot_axis,
                    band_node_cost=band_node_cost,
                    band_node_residual=band_node_residual,
                    final_aot=final_aot,
                )
            )
            return aot_arr, aot_unc_arr, cost_arr, diagnostics

        species_rt_models = _species_candidate_rt_models(
            rt_model=rt_model,
            config=self.config,
            template=atmo_prior.aot,
        )
        species_mode = str(getattr(self.config, "surface_driven_aerosol_species", "none"))
        species_labels = (
            list(_CANONICAL_6S_AEROSOL_PROFILES)
            if species_mode == "canonical_6s"
            else [f"candidate_{index}" for index in range(len(species_rt_models))]
        )
        candidate_rt_models = species_rt_models or [rt_model]
        if len(candidate_rt_models) == 1:
            aot_arr, aot_unc_arr, cost_arr, diagnostics = _resolve_candidate(
                candidate_rt_models[0],
                species_index=0 if species_rt_models else None,
            )
            if species_rt_models:
                diagnostics["surface_species_candidate_labels"] = species_labels
                diagnostics["surface_species_selected_candidate"] = 0
                diagnostics["surface_species_selected_profile"] = species_labels[0]
        else:
            aot_arr = np.full(shape, np.nan, dtype=np.float32)
            aot_unc_arr = np.full(shape, np.nan, dtype=np.float32)
            cost_arr = np.full(shape, np.inf, dtype=np.float32)
            diagnostics: dict[str, Any] = {}
            best_candidate_cost = float("inf")
            best_candidate_index = -1
            candidate_median_costs: list[float | None] = []
            for species_index, candidate_rt_model in enumerate(candidate_rt_models):
                cand_aot, cand_unc, cand_cost, cand_diagnostics = _resolve_candidate(
                    candidate_rt_model,
                    species_index=species_index,
                )
                finite_cand_cost = cand_cost[np.isfinite(cand_cost)]
                median_cand_cost = (
                    float(np.median(finite_cand_cost)) if finite_cand_cost.size else float("inf")
                )
                candidate_median_costs.append(
                    median_cand_cost if np.isfinite(median_cand_cost) else None
                )
                if median_cand_cost < best_candidate_cost:
                    best_candidate_cost = median_cand_cost
                    best_candidate_index = species_index
                    aot_arr = np.asarray(cand_aot, dtype=np.float32)
                    aot_unc_arr = np.asarray(cand_unc, dtype=np.float32)
                    cost_arr = np.asarray(cand_cost, dtype=np.float32)
                    diagnostics = cand_diagnostics
            diagnostics["surface_rt_branch_diagnostic_candidate"] = diagnostics.get(
                "surface_rt_branch"
            )
            diagnostics["surface_species_candidate_median_costs"] = candidate_median_costs
            diagnostics["surface_species_candidate_labels"] = species_labels
            diagnostics["surface_species_selected_candidate"] = best_candidate_index
            diagnostics["surface_species_selected_profile"] = (
                species_labels[best_candidate_index] if best_candidate_index >= 0 else None
            )
            diagnostics["surface_species_selected_median_cost"] = (
                best_candidate_cost if np.isfinite(best_candidate_cost) else None
            )
            diagnostics["surface_rt_branch"] = "tile_species_min_median_cost"

        # Keep pixels without usable observations unresolved. Filling them from
        # the atmospheric prior would make patch means depend on cloud/water/no-data
        # coverage and would misrepresent a fallback value as a surface retrieval.
        solved = np.isfinite(aot_arr)

        template = atmo_prior.aot
        aot_da = xr.DataArray(
            np.asarray(aot_arr, dtype=np.float32),
            dims=template.dims,
            coords=template.coords,
            name="aot",
        )
        aot_unc_da = xr.DataArray(
            np.asarray(aot_unc_arr, dtype=np.float32),
            dims=template.dims,
            coords=template.coords,
            name="aot_unc",
        )

        finite_cost = cost_arr[np.isfinite(cost_arr)]
        final_cost = float(np.median(finite_cost)) if finite_cost.size else float("inf")
        n_solved = int(np.count_nonzero(solved))
        diagnostics.update(
            {
                "surface_valid_observation_count": int(np.count_nonzero(solve_valid)),
                "surface_valid_observation_fraction": float(np.mean(solve_valid)),
                "surface_solved_pixel_count": n_solved,
                "surface_solved_pixel_fraction": float(n_solved / solve_valid.size),
                "surface_cloud_mask_bypassed": bool(allow_cloud_retrieval),
                "surface_water_mask_bypassed": bool(ignore_cloud_water),
            }
        )

        solved_atmo = atmo_prior.with_updated_aot_tcwv(
            aot=aot_da,
            tcwv=atmo_prior.tcwv,
            aot_unc=aot_unc_da,
            tcwv_unc=atmo_prior.tcwv_unc,
        )
        return SolverResult(
            aot=solved_atmo.aot,
            tcwv=solved_atmo.tcwv,
            aot_unc=solved_atmo.aot_unc,
            tcwv_unc=solved_atmo.tcwv_unc,
            n_iterations=1,
            final_cost=final_cost,
            success=n_solved > 0,
            message=f"surface_driven: {n_solved}/{ny * nx} pixels solved",
            atmo_state=solved_atmo,
            diagnostics=diagnostics,
        )

    def solve_bundle(self, bundle: SolverInputBundle) -> SolverResult:
        # Cost-gated tau-dependent prior: run the static prior first; only where
        # its per-band misfit is catastrophic (the fixed-anchor prediction has
        # gone off the dictionary manifold, e.g. thick smoke) re-solve with the
        # tau-dependent self-consistency prior. The static solve is a reliable
        # self-diagnosis (cost/band ~1 when it fits, ~20+ when it fails), so the
        # gate needs no external label. gate=None (default) preserves behaviour.
        gate = getattr(self.config, "surface_driven_tau_gate_cost", None)
        payload = getattr(bundle.surface_prior, "tau_predictor", None)
        tau_available = isinstance(payload, dict) and "trees" in payload

        def _annotate_tau_gate(
            result: SolverResult,
            *,
            fired: bool,
            static_cost_per_band: float | None = None,
            tau_cost_per_band: float | None = None,
        ) -> SolverResult:
            diagnostics = dict(getattr(result, "diagnostics", {}) or {})
            diagnostics.update(
                {
                    "surface_tau_gate_configured": gate is not None,
                    "surface_tau_available": bool(tau_available),
                    "surface_tau_gate_fired": bool(fired),
                    "surface_tau_gate_threshold": None if gate is None else float(gate),
                }
            )
            if static_cost_per_band is not None:
                diagnostics["surface_static_cost_per_band"] = float(static_cost_per_band)
            if tau_cost_per_band is not None:
                diagnostics["surface_tau_cost_per_band"] = float(tau_cost_per_band)
            result.diagnostics = diagnostics
            return result

        if gate is None or not tau_available:
            result = self.solve(
                bundle.toa,
                bundle.surface_prior,
                bundle.geometry,
                bundle.atmo_prior,
                bundle.rt_model,
                bundle.cloud_mask,
                bundle.bands,
                sharp_transition_mask=bundle.sharp_transition_mask,
                water_mask=bundle.water_mask,
                _force_tau_prior=None,
            )
            return _annotate_tau_gate(result, fired=False)

        static_result = self.solve(
            bundle.toa,
            bundle.surface_prior,
            bundle.geometry,
            bundle.atmo_prior,
            bundle.rt_model,
            bundle.cloud_mask,
            bundle.bands,
            sharp_transition_mask=bundle.sharp_transition_mask,
            water_mask=bundle.water_mask,
            _force_tau_prior=False,
        )
        n_solve_bands = len(bundle.bands)
        cost_per_band = static_result.final_cost / max(1, n_solve_bands)
        if not np.isfinite(cost_per_band) or cost_per_band <= float(gate):
            return _annotate_tau_gate(
                static_result,
                fired=False,
                static_cost_per_band=cost_per_band,
            )
        logger.info(
            "surface_driven: static cost/band %.2f > gate %.2f — re-solving with "
            "tau-dependent prior.",
            cost_per_band,
            float(gate),
        )
        tau_result = self.solve(
            bundle.toa,
            bundle.surface_prior,
            bundle.geometry,
            bundle.atmo_prior,
            bundle.rt_model,
            bundle.cloud_mask,
            bundle.bands,
            sharp_transition_mask=bundle.sharp_transition_mask,
            water_mask=bundle.water_mask,
            _force_tau_prior=True,
        )
        tau_cost_per_band = tau_result.final_cost / max(1, n_solve_bands)
        return _annotate_tau_gate(
            tau_result,
            fired=True,
            static_cost_per_band=cost_per_band,
            tau_cost_per_band=tau_cost_per_band,
        )
