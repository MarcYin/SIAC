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
    if mode != "cci_climatology_exact":
        raise ValueError(f"Unsupported surface-driven aerosol species mode {mode!r}.")
    if str(getattr(rt_model, "backend_name", "")) != "sixs":
        raise ValueError('surface_driven_aerosol_species requires algorithms.rt.backend="sixs".')

    clone_rt = rt_optional_capability(rt_model, "with_rt_setup")
    base_setup = getattr(rt_model, "rt_setup", None)
    if clone_rt is None or base_setup is None or not hasattr(base_setup, "model_copy"):
        raise ValueError(
            "Native 6S species mode requires an RT backend exposing rt_setup and with_rt_setup()."
        )

    lon, lat = _scene_center_lon_lat(template)
    month = _rt_month(rt_model)
    from siac.algorithms.rt.aerosol_species import climatology_cci_aerosol_setup

    aerosol_setups = (climatology_cci_aerosol_setup(lon, lat, month),)
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


def _backstop_conflict_diagnostics(
    *,
    surface_min_aot: object,
    aot_prior: np.ndarray,
    calibrated_unc: np.ndarray,
) -> dict[str, Any]:
    """Report how far the surface solution sits above the aerosol prior.

    Reported in units of the calibrated backstop sigma, so a large positive
    value means the surface likelihood and the prior disagree about a thick
    scene. Diagnostic only — it uses just the surface likelihood and the prior
    that the retrieval already computed, and never inspects a validation target.
    """
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
    diagnostics: dict[str, Any] = {
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
    }
    return diagnostics


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
        # Per-band absolute chi-square, summed over whichever bands are finite
        # at this pixel (a partial-band sum rather than dropping the pixel).
        band_names = [str(getattr(b, "name", b)) for b in bands]
        inv_var = 1.0 / np.square(np.maximum(boa_unc, 1e-6))  # (band, y, x)

        # -- tau-dependent surface prior (joint surface/AOD self-consistency) --
        # When enabled, the visible prior is re-predicted per AOT node from the
        # scene NIR/SWIR anchors corrected at THAT node's AOD, so the surface
        # term carries real information at high AOD where a fixed-anchor
        # prediction goes off the dictionary manifold.
        tau_payload = getattr(surface_prior, "tau_predictor", None)
        tau_config_on = bool(getattr(self.config, "surface_driven_tau_dependent_prior", False))
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
                finite = np.isfinite(abs_term)
                used = np.count_nonzero(finite, axis=0)
                cost = np.where(finite, abs_term, 0.0).sum(axis=0)
                cube[k] = np.where(solve_valid & (used > 0), cost, np.nan)
            diagnostics = {
                **_cost_curve_diagnostics(cube, aot_axis),
                "surface_cost_mode": "chi2",
                "surface_tau_prior_enabled": bool(tau_enabled),
                "surface_tcwv_cost_mode": (
                    "scalar_reference"
                    if reference_tcwv is not None
                    else ("fixed_spatial_prior" if tcwv_prior_is_spatial else "fixed_scalar_prior")
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
                diagnostics,
                band_node_cost,
                band_node_residual,
                band_cost_cube,
                band_residual_cube,
                band_signed_residual_cube,
            )

        def _dump_cost_cubes(
            cube: np.ndarray,
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
                dump_path = dump_path.with_name(f"{dump_path.stem}.species{species_index}{suffix}")
            template = atmo_prior.aot
            x_coord = template.coords.get("x")
            y_coord = template.coords.get("y")
            temporary_path = dump_path.with_name(f".{dump_path.name}.tmp.{os.getpid()}.npz")
            try:
                np.savez(
                    temporary_path,
                    cube=cube.astype(np.float32),
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
            (
                cube,
                diagnostics,
                band_node_cost,
                band_node_residual,
                band_cost_cube,
                band_residual_cube,
                band_signed_residual_cube,
            ) = _build_cost_cubes(candidate_rt_model)
            _dump_cost_cubes(
                cube,
                band_cost_cube,
                band_residual_cube,
                band_signed_residual_cube,
                species_index=species_index,
            )
            aot_arr, aot_unc_arr, cost_arr = _resolve(cube)
            diagnostics.update(
                _backstop_conflict_diagnostics(
                    surface_min_aot=diagnostics.get("surface_cost_curve_min_aot"),
                    aot_prior=aot_prior_backstop,
                    calibrated_unc=aot_prior_unc,
                )
            )
            finite_aot = aot_arr[np.isfinite(aot_arr)]
            final_aot = float(np.median(finite_aot)) if finite_aot.size else float("nan")

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
        # The validated species mode resolves exactly one aerosol setup per
        # scene, so there is a single candidate to solve.
        candidate_rt_model = species_rt_models[0] if species_rt_models else rt_model
        aot_arr, aot_unc_arr, cost_arr, diagnostics = _resolve_candidate(
            candidate_rt_model,
            species_index=0 if species_rt_models else None,
        )

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
        return self.solve(
            bundle.toa,
            bundle.surface_prior,
            bundle.geometry,
            bundle.atmo_prior,
            bundle.rt_model,
            bundle.cloud_mask,
            bundle.bands,
            sharp_transition_mask=bundle.sharp_transition_mask,
            water_mask=bundle.water_mask,
        )
