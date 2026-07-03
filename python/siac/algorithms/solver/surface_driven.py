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

from dataclasses import replace
from typing import TYPE_CHECKING, cast

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
    if mode != "cci_climatology":
        raise ValueError(f"Unsupported surface-driven aerosol species mode {mode!r}.")
    if str(getattr(rt_model, "backend_name", "")) != "sixs":
        raise ValueError(
            'surface_driven_aerosol_species="cci_climatology" requires '
            'algorithms.rt.backend="sixs".'
        )

    clone_rt = rt_optional_capability(rt_model, "with_rt_setup")
    base_setup = getattr(rt_model, "rt_setup", None)
    if clone_rt is None or base_setup is None or not hasattr(base_setup, "model_copy"):
        raise ValueError(
            "Native 6S species mode requires an RT backend exposing rt_setup and with_rt_setup()."
        )

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
    0.05 to 1.5, then 0.1 to 2.5 (~62 nodes). The fine low-AOD spacing avoids the
    coarse-log quantisation that snaps clean retrievals onto sparse nodes.
    """
    return cast(
        "np.ndarray",
        np.concatenate(
            [
                np.arange(0.01, 0.2, 0.01),
                np.arange(0.2, 0.5, 0.025),
                np.arange(0.5, 1.5, 0.05),
                np.arange(1.5, 2.6, 0.1),
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

    return {
        "mode": mode,
        "selected_pass": selected_pass,
        "aod": float(selected["aod"]),
        "aod_unc": selected.get("aod_unc"),
        "cost": selected.get("cost"),
        "radius_m": float(radius_m),
        "pool_window": int(pool_window),
        "min_count": int(min_count),
        "clean_threshold": float(clean_threshold),
        "high_threshold": float(high_threshold),
        "main": main,
        "abs": abs_result,
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
        return integrated_cost_field_aod(
            cube=z["cube"],
            cube_abs=cube_abs,
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

        # Valid = not cloudy/invalid, optionally also not water. The validated
        # surface-driven recipe (harness GATE=0) does NOT exclude cloud/water in
        # the solve — it pools the per-pixel surface-cost cube spatially and
        # relies on the robust spatial median + support gate, so high-AOD biomass
        # / near-river stations (whose TOA is finite and whose cost has a real
        # high-AOT minimum) are still solved instead of falling back to the
        # under-estimating MAIAC prior. ``surface_driven_ignore_cloud_water``
        # opts into that behaviour (default keeps the conservative masking).
        if getattr(self.config, "surface_driven_ignore_cloud_water", False):
            valid_mask: np.ndarray = np.ones(shape, dtype=bool)
        else:
            valid_mask = ~cloud_mask.values.astype(bool)
            if water_mask is not None:
                valid_mask = valid_mask & ~water_mask.values.astype(bool)

        aot_prior = atmo_prior.aot.values.astype(np.float32)
        tcwv_prior = atmo_prior.tcwv.values.astype(np.float32)
        aot_prior_unc_raw = atmo_prior.aot_unc.values.astype(np.float32)
        tcwv_prior_unc = np.maximum(atmo_prior.tcwv_unc.values.astype(np.float32), 1e-3)

        # AOT axis: log-spaced over the configured bounds (default) or the
        # validated fine non-uniform "acixthree" axis (denser at low AOD). TCWV
        # is fixed at the scene-median prior for the sweep.
        if str(getattr(self.config, "surface_driven_aot_axis", "log")) == "acixthree":
            aot_axis = _acixthree_aot_axis()
        else:
            n_aot = max(3, int(self.config.grid_search_aot_points))
            aot_axis = _log_aot_axis(self.config.aot_bounds[0], self.config.aot_bounds[1], n_aot)
        n_aot = int(aot_axis.size)
        # TCWV held at a single value for the AOT sweep. A configured reference
        # (the validated harness ran the cost-cube RT at TCWV=2.0, sea-level)
        # overrides the real scene-median so the per-node TOA→BOA correction
        # matches the reference atmosphere that cancels the surface prior's blue
        # dark bias; ``None`` (default) keeps the real scene-median behaviour.
        reference_tcwv = getattr(self.config, "surface_driven_reference_tcwv", None)
        if reference_tcwv is not None:
            tcwv_fixed = float(reference_tcwv)
        else:
            finite_tcwv = tcwv_prior[np.isfinite(tcwv_prior)]
            tcwv_fixed = (
                float(np.median(finite_tcwv))
                if finite_tcwv.size
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

        # Coefficient provider with TCWV fixed at the prior; one shared
        # (band, y, x) stack per AOT node. When a reference TCWV is configured,
        # replace the prior's TCWV field so BOTH the joint-LUT path (tcwv_axis)
        # and the per-candidate compute path (which holds the unsolved TCWV at the
        # template field) run the cost-cube RT at the reference value.
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
        cost_mode = str(getattr(self.config, "surface_driven_cost_mode", "chi2"))
        band_names = [str(getattr(b, "name", b)) for b in bands]
        anchor_index = band_names.index("B04") if "B04" in band_names else n_bands - 1
        abs_indices = [i for i, nm in enumerate(band_names) if nm != "B01"] or list(range(n_bands))
        all_indices = list(range(n_bands))
        two_pass = cost_mode == "auto2"
        use_shape_main = cost_mode in ("shape", "auto2")

        inv_var = 1.0 / np.square(np.maximum(boa_unc, 1e-6))  # (band, y, x)

        # Backstop centre fed to the pool-argmin is the raw aerosol prior; the
        # no-observation fallback fill below also keeps the raw prior.
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

        def _resolve(cube_in: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
            a, u, c = surface_driven_pool_argmin(
                np.ascontiguousarray(cube_in, dtype=np.float64),
                aot_axis.astype(np.float64, copy=False),
                np.ascontiguousarray(aot_prior_backstop, dtype=np.float64),
                np.ascontiguousarray(aot_prior_unc, dtype=np.float64),
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
        ) -> tuple[np.ndarray, np.ndarray | None]:
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
            cube: np.ndarray = np.full((n_aot, ny, nx), np.nan, dtype=np.float64)
            cube_abs: np.ndarray | None = (
                np.full((n_aot, ny, nx), np.nan, dtype=np.float64) if two_pass else None
            )
            for k in range(n_aot):
                xap, xbp, xcp = coeff_provider(float(aot_axis[k]), tcwv_fixed)
                y = np.asarray(xap, np.float64) * toa_values - np.asarray(xbp, np.float64)
                denom = 1.0 + np.asarray(xcp, np.float64) * y
                with np.errstate(invalid="ignore", divide="ignore"):
                    boa = np.where(np.abs(denom) > 1e-12, y / denom, np.nan)
                if use_shape_main:
                    cube[k] = _shape_cost_node(boa, boa_prior, boa_unc, all_indices, anchor_index)
                else:
                    diff = boa - boa_prior
                    term = inv_var * diff * diff  # (band, y, x)
                    finite = np.isfinite(term)
                    used = np.count_nonzero(finite, axis=0)
                    cost = np.where(finite, term, 0.0).sum(axis=0)
                    cube[k] = np.where(used > 0, cost, np.nan)
                if cube_abs is not None:
                    cube_abs[k] = _abs_cost_node(boa, boa_prior, inv_var, abs_indices)
            return cube, cube_abs

        def _dump_cost_cubes(
            cube: np.ndarray,
            cube_abs: np.ndarray | None,
            *,
            species_index: int | None,
        ) -> None:
            # Optional diagnostic dump of the raw per-node cost cube(s) + grid,
            # for offline analysis (env-gated; no-op by default). When species
            # candidates are active, each candidate gets its own suffixed file.
            dump_path = __import__("os").environ.get("SIAC_DUMP_COST_CUBE")
            if not dump_path:
                return
            if species_index is not None:
                root, ext = __import__("os").path.splitext(dump_path)
                dump_path = f"{root}.species{species_index}{ext or '.npz'}"
            template = atmo_prior.aot
            x_coord = template.coords.get("x")
            y_coord = template.coords.get("y")
            __import__("numpy").savez(
                dump_path,
                cube=cube.astype(np.float32),
                cube_abs=(cube_abs.astype(np.float32) if cube_abs is not None else np.zeros(0)),
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

        def _resolve_candidate(
            candidate_rt_model: RTModelBackend,
            *,
            species_index: int | None,
        ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
            cube, cube_abs = _build_cost_cubes(candidate_rt_model)
            _dump_cost_cubes(cube, cube_abs, species_index=species_index)
            aot_arr, aot_unc_arr, cost_arr = _resolve(cube)
            if two_pass and cube_abs is not None:
                aot_abs, aot_unc_abs, cost_abs = _resolve(cube_abs)
                # Regime by the absolute-solved AOD: use the shape solve only in
                # the moderate band; absolute (non-B01) at clean and high tails.
                clean_thr = float(getattr(self.config, "surface_driven_aod_clean", 0.15))
                high_thr = float(getattr(self.config, "surface_driven_aod_high", 0.6))
                use_main = np.isfinite(aot_abs) & (aot_abs >= clean_thr) & (aot_abs <= high_thr)
                aot_arr = np.where(use_main, aot_arr, aot_abs).astype(np.float32)
                aot_unc_arr = np.where(use_main, aot_unc_arr, aot_unc_abs).astype(np.float32)
                cost_arr = np.where(use_main, cost_arr, cost_abs).astype(np.float32)
            return aot_arr, aot_unc_arr, cost_arr

        species_rt_models = _species_candidate_rt_models(
            rt_model=rt_model,
            config=self.config,
            template=atmo_prior.aot,
        )
        candidate_rt_models = species_rt_models or [rt_model]
        if len(candidate_rt_models) == 1:
            aot_arr, aot_unc_arr, cost_arr = _resolve_candidate(
                candidate_rt_models[0],
                species_index=0 if species_rt_models else None,
            )
        else:
            aot_arr = np.full(shape, np.nan, dtype=np.float32)
            aot_unc_arr = np.full(shape, np.nan, dtype=np.float32)
            cost_arr = np.full(shape, np.inf, dtype=np.float32)
            for species_index, candidate_rt_model in enumerate(candidate_rt_models):
                cand_aot, cand_unc, cand_cost = _resolve_candidate(
                    candidate_rt_model,
                    species_index=species_index,
                )
                better = np.isfinite(cand_cost) & (cand_cost < cost_arr)
                aot_arr = np.where(better, cand_aot, aot_arr).astype(np.float32)
                aot_unc_arr = np.where(better, cand_unc, aot_unc_arr).astype(np.float32)
                cost_arr = np.where(better, cand_cost, cost_arr).astype(np.float32)

        # Pixels with no usable observation fall back to the prior (backstop).
        solved = np.isfinite(aot_arr)
        aot_filled = np.where(solved, aot_arr, aot_prior).astype(np.float32)
        aot_unc_filled = np.where(solved, aot_unc_arr, aot_prior_unc).astype(np.float32)

        template = atmo_prior.aot
        aot_da = xr.DataArray(aot_filled, dims=template.dims, coords=template.coords, name="aot")
        aot_unc_da = xr.DataArray(
            aot_unc_filled, dims=template.dims, coords=template.coords, name="aot_unc"
        )

        finite_cost = cost_arr[np.isfinite(cost_arr)]
        final_cost = float(np.median(finite_cost)) if finite_cost.size else float("inf")
        n_solved = int(np.count_nonzero(solved))

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
