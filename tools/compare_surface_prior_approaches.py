#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path
from typing import Any

import numpy as np
import xarray as xr
from rasterio.enums import Resampling

from siac.algorithms.solver.multigrid import build_solver_valid_mask
from siac.algorithms.surface.swir_refine import resample_geometry_for_surface_prior
from siac.algorithms.surface.kernel_model import KernelModelDeriver
from siac.adapters.earthdata import earthaccess_cache_dir
from siac.adapters.earthdata_common import build_target_template
from siac.app._assembly_surface import (
    _select_visible_surface_prior_bands,
    _surface_prior_brdf_request,
    _surface_prior_mapping_state,
)
from siac.app.assembly import (
    _prepare_monthly_surface_prior_runtime,
    _query_monthly_surface_prior,
    resolve_brdf_provider,
    resolve_monthly_composite_provider,
)
from siac.app.planning import build_execution_plan
from siac.app.requests import SceneProcessRequest
from siac.config import SIACConfig
from siac.runtime import AtmosphericState, GeometryAngles, SurfacePrior
from siac.workflows.pipeline import _call_grid_assembler, _resample_field_to_template, _select_band_slice

LOGGER = logging.getLogger("compare_surface_prior_approaches")

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except ModuleNotFoundError:
    matplotlib = None
    plt = None


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run a detailed surface-prior experiment that compares "
            "monthly_database/SWIR_refine against direct same-day MCD43/kernel-model priors."
        )
    )
    parser.add_argument("--config", type=Path, required=True, help="SIAC TOML config path.")
    parser.add_argument(
        "--input",
        type=Path,
        default=None,
        help="Scene input path. If omitted, the script tries to infer a single SAFE under the S2 cache.",
    )
    parser.add_argument(
        "--reference-output-dir",
        type=Path,
        default=None,
        help=(
            "Existing SIAC output directory for the same scene. "
            "If provided, the script compares against stored surface_prior.nc and "
            "reconstructs the run's solved atmosphere from auxiliary.nc."
        ),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="Directory where experiment outputs will be written.",
    )
    parser.add_argument(
        "--sensor",
        default="s2",
        help="Sensor override for config resolution. Defaults to s2.",
    )
    parser.add_argument(
        "--monthly-source",
        default="auto",
        choices=("auto", "reference", "recompute"),
        help=(
            "Where to source the monthly_database prior. "
            "'auto' prefers reference-output-dir/surface_prior.nc when available."
        ),
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=("DEBUG", "INFO", "WARNING", "ERROR"),
        help="Logging verbosity.",
    )
    return parser.parse_args()


def _setup_logging(level: str) -> None:
    logging.basicConfig(
        level=getattr(logging, level.upper()),
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )


def _infer_input_path(config: Any, explicit_input: Path | None) -> Path:
    if explicit_input is not None:
        return explicit_input
    if getattr(config.run, "input_path", None) is not None:
        return Path(config.run.input_path)
    s2_cache = getattr(config.paths.caches, "s2", None) or getattr(config.s2_data, "cache_dir", None)
    if s2_cache is None:
        raise ValueError("Could not infer input path: no explicit --input and no S2 cache configured.")
    safes = sorted(Path(s2_cache).glob("*.SAFE"))
    if len(safes) != 1:
        raise ValueError(
            f"Could not infer input path from {s2_cache}: expected 1 SAFE, found {len(safes)}."
        )
    return safes[0]


def _to_dataset(data: xr.DataArray, *, prefix: str) -> xr.Dataset:
    if "band" not in data.dims:
        return xr.Dataset({prefix: data.astype(np.float32)})
    band_values = [str(value) for value in np.asarray(data.coords["band"].values)]
    return xr.Dataset(
        {
            f"{prefix}_{band_name}": data.sel(band=band_name, drop=True).astype(np.float32)
            for band_name in band_values
        }
    )


def _mask_for_band(mask: xr.DataArray, band_name: str) -> xr.DataArray:
    if "band" not in mask.dims:
        return mask.astype(bool)
    band_values = [str(value) for value in np.asarray(mask.coords["band"].values)]
    if band_name in band_values:
        return mask.sel(band=band_name, drop=True).astype(bool)
    return mask.isel(band=0, drop=True).astype(bool)


def _corrcoef(a: np.ndarray, b: np.ndarray) -> float:
    if a.size < 2 or b.size < 2:
        return float("nan")
    if np.allclose(a, a[0]) or np.allclose(b, b[0]):
        return float("nan")
    return float(np.corrcoef(a, b)[0, 1])


def _band_metrics(
    left: xr.DataArray,
    right: xr.DataArray,
    *,
    left_mask: xr.DataArray,
    right_mask: xr.DataArray,
) -> dict[str, Any]:
    valid = (
        left_mask.values.astype(bool)
        & right_mask.values.astype(bool)
        & np.isfinite(left.values)
        & np.isfinite(right.values)
    )
    count = int(np.count_nonzero(valid))
    if count == 0:
        return {"valid_count": 0}

    left_values = np.asarray(left.values[valid], dtype=np.float64)
    right_values = np.asarray(right.values[valid], dtype=np.float64)
    delta = right_values - left_values

    return {
        "valid_count": count,
        "left_mean": float(np.mean(left_values)),
        "right_mean": float(np.mean(right_values)),
        "bias_right_minus_left": float(np.mean(delta)),
        "median_right_minus_left": float(np.median(delta)),
        "mae_right_minus_left": float(np.mean(np.abs(delta))),
        "rmse_right_minus_left": float(np.sqrt(np.mean(delta**2))),
        "p05_right_minus_left": float(np.quantile(delta, 0.05)),
        "p50_right_minus_left": float(np.quantile(delta, 0.50)),
        "p95_right_minus_left": float(np.quantile(delta, 0.95)),
        "corr": _corrcoef(left_values, right_values),
    }


def _surface_prior_metrics(
    monthly_prior: SurfacePrior,
    direct_prior: SurfacePrior,
) -> dict[str, Any]:
    band_values = [str(value) for value in np.asarray(monthly_prior.boa.coords["band"].values)]
    metrics: dict[str, Any] = {}
    for band_name in band_values:
        monthly_band = monthly_prior.boa.sel(band=band_name, drop=True)
        direct_band = direct_prior.boa.sel(band=band_name, drop=True)
        metrics[band_name] = _band_metrics(
            monthly_band,
            direct_band,
            left_mask=_mask_for_band(monthly_prior.mask, band_name),
            right_mask=_mask_for_band(direct_prior.mask, band_name),
        )
    return metrics


def _same_spatial_grid(left: SurfacePrior, right: SurfacePrior) -> bool:
    left_boa = left.boa
    right_boa = right.boa
    if left_boa.sizes.get("y") != right_boa.sizes.get("y") or left_boa.sizes.get("x") != right_boa.sizes.get("x"):
        return False
    for axis in ("y", "x"):
        left_coord = np.asarray(left_boa.coords[axis].values)
        right_coord = np.asarray(right_boa.coords[axis].values)
        if left_coord.shape != right_coord.shape:
            return False
        if not np.allclose(left_coord, right_coord):
            return False
    return True


def _difference_dataset(
    monthly_prior: SurfacePrior,
    direct_prior: SurfacePrior,
) -> xr.Dataset:
    monthly = monthly_prior.boa.astype(np.float32)
    direct = direct_prior.boa.astype(np.float32)
    delta = (direct - monthly).astype(np.float32)
    ds = xr.Dataset()
    ds.update(_to_dataset(monthly, prefix="monthly"))
    ds.update(_to_dataset(direct, prefix="direct"))
    ds.update(_to_dataset(delta, prefix="direct_minus_monthly"))
    ds["monthly_mask"] = monthly_prior.mask.astype(np.int16)
    ds["direct_mask"] = direct_prior.mask.astype(np.int16)
    return ds


def _surface_prior_from_reference(path: Path) -> SurfacePrior:
    ds = xr.load_dataset(path)
    boa = ds.to_array(dim="band").astype(np.float32)
    mask = xr.DataArray(
        np.all(np.isfinite(boa.values), axis=0),
        dims=("y", "x"),
        coords={"y": boa.coords["y"], "x": boa.coords["x"]},
    )
    boa_unc = xr.zeros_like(boa, dtype=np.float32)
    return SurfacePrior(boa=boa, boa_unc=boa_unc, kernels=None, mask=mask)


def _align_single_band(
    data: xr.DataArray,
    template: xr.DataArray,
    *,
    resampling: Resampling,
) -> xr.DataArray:
    if tuple(data.shape) == tuple(template.shape) and all(
        np.array_equal(np.asarray(data.coords[dim].values), np.asarray(template.coords[dim].values))
        for dim in ("y", "x")
        if dim in data.coords and dim in template.coords
    ):
        return data

    try:
        import rioxarray  # noqa: F401

        src = data
        dst = template
        if "x" in src.dims and "y" in src.dims:
            src = src.rio.set_spatial_dims(x_dim="x", y_dim="y", inplace=False)
        if "x" in dst.dims and "y" in dst.dims:
            dst = dst.rio.set_spatial_dims(x_dim="x", y_dim="y", inplace=False)
        if src.rio.crs is not None and dst.rio.crs is not None:
            return src.rio.reproject_match(dst, resampling=resampling)
    except Exception:
        pass

    method = "nearest" if resampling == Resampling.nearest else "linear"
    return data.interp(x=template.coords["x"], y=template.coords["y"], method=method)


def _align_banded_to_template(
    data: xr.DataArray,
    template: xr.DataArray,
    *,
    resampling: Resampling,
) -> xr.DataArray:
    if "band" not in data.dims:
        return _align_single_band(data, template, resampling=resampling)
    band_values = [str(value) for value in np.asarray(data.coords["band"].values)]
    aligned = [
        _align_single_band(data.sel(band=band_name, drop=True), template, resampling=resampling).assign_coords(
            band=band_name
        )
        for band_name in band_values
    ]
    return xr.concat(aligned, dim="band").astype(np.float32)


def _align_surface_prior_to_reference(
    prior: SurfacePrior,
    reference_prior: SurfacePrior,
) -> SurfacePrior:
    template = reference_prior.boa.isel(band=0, drop=True) if "band" in reference_prior.boa.dims else reference_prior.boa
    aligned_boa = _align_banded_to_template(prior.boa, template, resampling=Resampling.bilinear)
    aligned_unc = _align_banded_to_template(prior.boa_unc, template, resampling=Resampling.bilinear)
    aligned_mask = _align_banded_to_template(prior.mask.astype(np.float32), template, resampling=Resampling.nearest) > 0.5
    return SurfacePrior(
        boa=aligned_boa.astype(np.float32),
        boa_unc=aligned_unc.astype(np.float32),
        kernels=prior.kernels,
        mask=aligned_mask.astype(bool),
        monthly_composites=prior.monthly_composites,
    )


def _align_geometry_to_template(
    geometry: GeometryAngles,
    template: xr.DataArray,
) -> GeometryAngles:
    return GeometryAngles(
        sza=_align_single_band(geometry.sza, template, resampling=Resampling.bilinear),
        saa=_align_single_band(geometry.saa, template, resampling=Resampling.bilinear),
        vza=_align_single_band(geometry.vza, template, resampling=Resampling.bilinear),
        vaa=_align_single_band(geometry.vaa, template, resampling=Resampling.bilinear),
    )


def _reconstruct_run_atmosphere(reference_output_dir: Path, atmo_prior: AtmosphericState) -> AtmosphericState:
    aux = xr.load_dataset(reference_output_dir / "auxiliary.nc")
    template = aux["aot"].astype(np.float32)
    prior_aot = _resample_field_to_template(atmo_prior.aot, template).astype(np.float32)
    prior_tcwv = _resample_field_to_template(atmo_prior.tcwv, template).astype(np.float32)
    return AtmosphericState(
        aot=xr.where(np.isfinite(template), template, prior_aot).astype(np.float32),
        tcwv=xr.where(np.isfinite(aux["tcwv"]), aux["tcwv"], prior_tcwv).astype(np.float32),
        tco3=_resample_field_to_template(atmo_prior.tco3, template).astype(np.float32),
        aot_unc=_resample_field_to_template(atmo_prior.aot_unc, template).astype(np.float32),
        tcwv_unc=_resample_field_to_template(atmo_prior.tcwv_unc, template).astype(np.float32),
        tco3_unc=_resample_field_to_template(atmo_prior.tco3_unc, template).astype(np.float32),
        elevation=_resample_field_to_template(atmo_prior.elevation, template).astype(np.float32),
    )


def _prepare_direct_runtime(
    plan: Any,
    obs: Any,
) -> dict[str, Any]:
    visible_bands = _select_visible_surface_prior_bands(obs.sensor_config)
    geometry = resample_geometry_for_surface_prior(
        obs,
        resolution=float(plan.config.solver.aerosol_resolution),
    )
    monthly_provider = resolve_monthly_composite_provider(plan.config, auth=plan.auth)
    source_bands = getattr(monthly_provider, "source_bands", None)
    if source_bands is None:
        monthly_composites = monthly_provider.get_monthly_composites(
            obs,
            float(plan.config.solver.aerosol_resolution),
        )
        source_bands = tuple(monthly_composites.source_bands)
    spectral_library, spectral_k_neighbors = _surface_prior_mapping_state(
        plan.config,
        source_bands=tuple(source_bands),
        target_bands=visible_bands,
        context="direct same-day MCD43 surface-prior comparison",
    )
    return {
        "visible_bands": visible_bands,
        "geometry": geometry,
        "spectral_library": spectral_library,
        "spectral_k_neighbors": spectral_k_neighbors,
    }


def _toa_simulation_metrics(
    plan: Any,
    obs: Any,
    atmo_state: AtmosphericState,
    surface_prior: SurfacePrior,
    *,
    label: str,
) -> tuple[dict[str, Any], xr.Dataset]:
    aerosol_resolution = float(plan.config.solver.aerosol_resolution)
    solver_inputs = _call_grid_assembler(
        plan.grid_assembler,
        obs,
        atmo_state,
        surface_prior,
        plan.rt_model,
        aerosol_resolution_m=aerosol_resolution,
    )
    valid_mask = build_solver_valid_mask(
        solver_inputs.cloud_mask,
        solver_inputs.toa,
        solver_inputs.surface_prior,
    ).values.astype(bool)

    metrics: dict[str, Any] = {}
    outputs = xr.Dataset()
    for band_index, band in enumerate(solver_inputs.bands):
        band_name = band.name
        toa_band = _select_band_slice(solver_inputs.toa, band_name=band_name, band_index=band_index)
        surface_band = _select_band_slice(
            solver_inputs.surface_prior.boa,
            band_name=band_name,
            band_index=band_index,
        )
        if toa_band is None or surface_band is None:
            continue
        coeffs = solver_inputs.rt_model.compute_coefficients(
            solver_inputs.geometry,
            solver_inputs.atmo_prior,
            band,
            False,
        )
        simulated = coeffs.simulate_toa(surface_band).astype(np.float32)
        band_valid = (
            valid_mask
            & np.isfinite(toa_band.values)
            & np.isfinite(simulated.values)
            & np.isfinite(surface_band.values)
        )
        count = int(np.count_nonzero(band_valid))
        if count == 0:
            metrics[band_name] = {"valid_count": 0}
            continue
        observed_values = np.asarray(toa_band.values[band_valid], dtype=np.float64)
        simulated_values = np.asarray(simulated.values[band_valid], dtype=np.float64)
        residual = simulated_values - observed_values
        metrics[band_name] = {
            "valid_count": count,
            "observed_mean": float(np.mean(observed_values)),
            "simulated_mean": float(np.mean(simulated_values)),
            "bias_simulated_minus_observed": float(np.mean(residual)),
            "mae_simulated_minus_observed": float(np.mean(np.abs(residual))),
            "rmse_simulated_minus_observed": float(np.sqrt(np.mean(residual**2))),
            "corr": _corrcoef(observed_values, simulated_values),
        }
        outputs[f"observed_toa_{band_name}"] = toa_band.astype(np.float32)
        outputs[f"simulated_toa_{label}_{band_name}"] = simulated
        outputs[f"residual_toa_{label}_{band_name}"] = (simulated - toa_band).astype(np.float32)
    return metrics, outputs


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def _finite_quantile_limits(values: np.ndarray, low: float = 0.02, high: float = 0.98) -> tuple[float, float]:
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return (0.0, 1.0)
    lower = float(np.quantile(finite, low))
    upper = float(np.quantile(finite, high))
    if np.isclose(lower, upper):
        upper = lower + 1.0e-6
    return lower, upper


def _plot_prior_maps(
    monthly_band: xr.DataArray,
    direct_band: xr.DataArray,
    *,
    band_name: str,
    output_path: Path,
) -> None:
    if plt is None:
        return
    delta = (direct_band - monthly_band).astype(np.float32)
    common_low, common_high = _finite_quantile_limits(
        np.concatenate(
            [
                np.asarray(monthly_band.values, dtype=np.float32).ravel(),
                np.asarray(direct_band.values, dtype=np.float32).ravel(),
            ]
        )
    )
    delta_low, delta_high = _finite_quantile_limits(np.asarray(delta.values, dtype=np.float32).ravel())
    delta_abs = max(abs(delta_low), abs(delta_high))

    fig, axes = plt.subplots(1, 3, figsize=(12, 4), constrained_layout=True)
    panels = (
        (monthly_band, "monthly_database", "viridis", common_low, common_high),
        (direct_band, "direct_mcd43", "viridis", common_low, common_high),
        (delta, "direct - monthly", "RdBu_r", -delta_abs, delta_abs),
    )
    for axis, (panel, title, cmap, vmin, vmax) in zip(axes, panels, strict=True):
        image = axis.imshow(panel.values, cmap=cmap, vmin=vmin, vmax=vmax)
        axis.set_title(f"{band_name}: {title}")
        axis.set_xticks([])
        axis.set_yticks([])
        fig.colorbar(image, ax=axis, shrink=0.85)
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def _plot_toa_residuals(
    observed: xr.DataArray,
    simulated_monthly: xr.DataArray,
    simulated_direct: xr.DataArray,
    *,
    band_name: str,
    output_path: Path,
) -> None:
    if plt is None:
        return
    valid_monthly = (
        np.isfinite(observed.values)
        & np.isfinite(simulated_monthly.values)
    )
    valid_direct = (
        np.isfinite(observed.values)
        & np.isfinite(simulated_direct.values)
    )
    obs_monthly = np.asarray(observed.values[valid_monthly], dtype=np.float32)
    sim_monthly = np.asarray(simulated_monthly.values[valid_monthly], dtype=np.float32)
    obs_direct = np.asarray(observed.values[valid_direct], dtype=np.float32)
    sim_direct = np.asarray(simulated_direct.values[valid_direct], dtype=np.float32)
    residual_monthly = sim_monthly - obs_monthly
    residual_direct = sim_direct - obs_direct

    all_obs = np.concatenate([obs_monthly, obs_direct]) if obs_monthly.size and obs_direct.size else obs_monthly
    all_sim = np.concatenate([sim_monthly, sim_direct]) if sim_monthly.size and sim_direct.size else sim_monthly
    if all_obs.size == 0 or all_sim.size == 0:
        return
    lower = float(min(np.nanmin(all_obs), np.nanmin(all_sim)))
    upper = float(max(np.nanmax(all_obs), np.nanmax(all_sim)))
    if np.isclose(lower, upper):
        upper = lower + 1.0e-6

    fig, axes = plt.subplots(1, 3, figsize=(14, 4), constrained_layout=True)
    axes[0].hexbin(obs_monthly, sim_monthly, gridsize=60, mincnt=1, cmap="Blues")
    axes[0].plot([lower, upper], [lower, upper], "k--", linewidth=1)
    axes[0].set_title(f"{band_name}: observed vs monthly")
    axes[0].set_xlabel("Observed TOA")
    axes[0].set_ylabel("Simulated TOA")

    axes[1].hexbin(obs_direct, sim_direct, gridsize=60, mincnt=1, cmap="Greens")
    axes[1].plot([lower, upper], [lower, upper], "k--", linewidth=1)
    axes[1].set_title(f"{band_name}: observed vs direct")
    axes[1].set_xlabel("Observed TOA")
    axes[1].set_ylabel("Simulated TOA")

    bins = 80
    axes[2].hist(residual_monthly, bins=bins, alpha=0.6, label="monthly", color="tab:blue")
    axes[2].hist(residual_direct, bins=bins, alpha=0.6, label="direct", color="tab:green")
    axes[2].axvline(0.0, color="k", linestyle="--", linewidth=1)
    axes[2].set_title(f"{band_name}: TOA residuals")
    axes[2].set_xlabel("Simulated - observed")
    axes[2].set_ylabel("Count")
    axes[2].legend()

    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def _granule_name(granule: Any) -> str:
    render_dict = getattr(granule, "render_dict", None)
    if render_dict is not None:
        umm = render_dict.get("umm", {})
        if isinstance(umm, dict):
            value = umm.get("GranuleUR")
            if isinstance(value, str):
                return value
    if isinstance(granule, dict):
        umm = granule.get("umm", {})
        if isinstance(umm, dict):
            value = umm.get("GranuleUR")
            if isinstance(value, str):
                return value
    return ""


def _resolve_exact_day_brdf_paths(
    brdf_provider: Any,
    obs: Any,
) -> list[Path]:
    obs_time = obs.metadata["observation_time"]
    granule_stamp = f"A{obs_time.strftime('%Y%j')}"
    short_name = brdf_provider._resolved_short_name()
    cache_dir = earthaccess_cache_dir(brdf_provider.cache_dir, short_name)
    exact_day_paths = sorted(cache_dir.glob(f"*{granule_stamp}*.hdf"))
    if exact_day_paths:
        LOGGER.info(
            "Using %d cached exact-day BRDF granule(s) for %s",
            len(exact_day_paths),
            granule_stamp,
        )
        return brdf_provider._select_candidate_paths(
            exact_day_paths,
            obs_time,
            obs.bounds,
            obs.crs,
            sample_dates=np.array([np.datetime64(obs_time.date(), "D")]),
        )

    LOGGER.info("Searching Earthdata for exact-day BRDF granules matching %s", granule_stamp)
    granules = brdf_provider.source.search_granules(
        short_name=short_name,
        bounds=obs.bounds,
        crs=obs.crs,
        temporal=(
            f"{obs_time.date().isoformat()}T00:00:00Z",
            f"{obs_time.date().isoformat()}T23:59:59Z",
        ),
        provider=brdf_provider.provider,
        count=-1,
    )
    exact_day_granules = [granule for granule in granules if granule_stamp in _granule_name(granule)]
    if not exact_day_granules:
        raise RuntimeError(f"No exact-day BRDF granules found for {granule_stamp}")

    LOGGER.info(
        "Downloading %d exact-day BRDF granule(s) for %s",
        len(exact_day_granules),
        granule_stamp,
    )
    downloaded = brdf_provider.source.download_granules(exact_day_granules, cache_dir)
    return brdf_provider._select_candidate_paths(
        downloaded,
        obs_time,
        obs.bounds,
        obs.crs,
        sample_dates=np.array([np.datetime64(obs_time.date(), "D")]),
    )


def _build_direct_prior(
    plan: Any,
    obs: Any,
    direct_runtime: dict[str, Any],
) -> SurfacePrior:
    brdf_provider = resolve_brdf_provider(plan.config, auth=plan.auth)
    request = _surface_prior_brdf_request(
        obs,
        brdf_provider=brdf_provider,
        target_resolution=float(plan.config.solver.aerosol_resolution),
        temporal_window=plan.config.brdf.temporal_window,
    )
    exact_day_paths = _resolve_exact_day_brdf_paths(brdf_provider, obs)
    requested = brdf_provider._resolve_requested_bands(list(request["bands"]))
    target_template = build_target_template(
        request["bounds"],
        request["crs"],
        request["target_resolution"],
    )
    payload = brdf_provider._merge_requested_payload(
        exact_day_paths,
        requested=requested,
        bounds=request["bounds"],
        crs=request["crs"],
        target_resolution=request["target_resolution"],
        target_template=target_template,
    )
    params, unc = brdf_provider._unpack_payload_stack(payload, requested=requested)
    if np.isfinite(params.values).any():
        brdf_weights = brdf_provider._weights_from_layers(
            brdf_provider._fill_parameter_defaults(params),
            unc.fillna(0.08),
        )
    else:
        brdf_weights = brdf_provider._weights_from_layers(
            params.astype(np.float32),
            unc.astype(np.float32),
        )
    deriver = KernelModelDeriver(
        psf_sigma_x=plan.config.surface_prior.psf_sigma_x,
        psf_sigma_y=plan.config.surface_prior.psf_sigma_y,
        apply_psf=plan.config.surface_prior.apply_psf,
    )
    return deriver.compute_surface_prior(
        brdf_weights,
        direct_runtime["geometry"],
        source_bands=brdf_provider.source_bands,
        target_bands=direct_runtime["visible_bands"],
        spectral_library=direct_runtime["spectral_library"],
        spectral_k_neighbors=direct_runtime["spectral_k_neighbors"],
    )


def main() -> None:
    args = _parse_args()
    _setup_logging(args.log_level)

    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    plots_dir = output_dir / "plots"
    plots_dir.mkdir(parents=True, exist_ok=True)

    config = SIACConfig.from_file(args.config).with_overrides(sensor=args.sensor)
    input_path = _infer_input_path(config, args.input)
    request = SceneProcessRequest(config=config, input_path=input_path)
    plan = build_execution_plan(request)

    LOGGER.info("Preprocessing observation from %s", input_path)
    obs = plan.preprocessor(plan.input_path, plan.runtime_aoi)
    aerosol_resolution = float(plan.config.solver.aerosol_resolution)

    LOGGER.info("Fetching atmospheric prior")
    atmo_prior = plan.atmo_provider(
        obs.bounds,
        obs.crs,
        obs.metadata["observation_time"],
        aerosol_resolution,
    )

    reference_output_dir = args.reference_output_dir.resolve() if args.reference_output_dir is not None else None
    stored_prior_path = reference_output_dir / "surface_prior.nc" if reference_output_dir is not None else None
    use_reference_monthly = (
        args.monthly_source in {"auto", "reference"}
        and stored_prior_path is not None
        and stored_prior_path.exists()
    )

    if use_reference_monthly:
        LOGGER.info("Using stored monthly_database prior from %s", stored_prior_path)
        monthly_prior = _surface_prior_from_reference(stored_prior_path)
        monthly_source_label = "reference_output.surface_prior.nc"
    else:
        LOGGER.info("Building monthly_database runtime and prior")
        monthly_provider = resolve_monthly_composite_provider(plan.config, auth=plan.auth)
        monthly_runtime = _prepare_monthly_surface_prior_runtime(
            plan.config,
            monthly_provider,
            observation=obs,
            resolution=aerosol_resolution,
        )
        monthly_prior = _query_monthly_surface_prior(obs, atmo_prior, plan.rt_model, monthly_runtime)
        monthly_source_label = "recomputed_monthly_database"

    LOGGER.info("Preparing direct same-day MCD43 runtime")
    direct_runtime = _prepare_direct_runtime(plan, obs)
    LOGGER.info("Building direct same-day MCD43/kernel-model prior")
    direct_prior = _build_direct_prior(plan, obs, direct_runtime)
    LOGGER.info("Aligning direct prior onto the monthly prior grid")
    direct_prior = _align_surface_prior_to_reference(direct_prior, monthly_prior)

    LOGGER.info("Computing prior comparison metrics")
    prior_metrics = _surface_prior_metrics(monthly_prior, direct_prior)
    difference_ds = _difference_dataset(monthly_prior, direct_prior)
    difference_path = output_dir / "surface_prior_comparison.nc"
    difference_ds.to_netcdf(difference_path)

    summary: dict[str, Any] = {
        "config_path": str(args.config.resolve()),
        "input_path": str(input_path.resolve()),
        "aerosol_resolution_m": aerosol_resolution,
        "observation_time": obs.metadata["observation_time"].isoformat(),
        "monthly_prior_source": monthly_source_label,
        "plots_enabled": plt is not None,
        "surface_prior_metrics": prior_metrics,
    }

    if reference_output_dir is not None:
        auxiliary_path = reference_output_dir / "auxiliary.nc"
        if stored_prior_path is not None and stored_prior_path.exists() and not use_reference_monthly:
            LOGGER.info("Comparing recomputed monthly prior against stored run prior")
            stored_prior = _surface_prior_from_reference(stored_prior_path)
            if _same_spatial_grid(monthly_prior, stored_prior):
                summary["stored_run_surface_prior_metrics"] = _surface_prior_metrics(monthly_prior, stored_prior)
            else:
                summary["stored_run_surface_prior_metrics"] = {
                    "skipped": True,
                    "reason": "stored run prior grid does not match recomputed aerosol-grid prior",
                    "stored_shape": [
                        int(stored_prior.boa.sizes.get("y", 0)),
                        int(stored_prior.boa.sizes.get("x", 0)),
                    ],
                    "recomputed_shape": [
                        int(monthly_prior.boa.sizes.get("y", 0)),
                        int(monthly_prior.boa.sizes.get("x", 0)),
                    ],
                }
        if auxiliary_path.exists():
            LOGGER.info("Reconstructing run atmosphere from %s", auxiliary_path)
            run_atmo = _reconstruct_run_atmosphere(reference_output_dir, atmo_prior)
            try:
                monthly_run_metrics, monthly_run_toa = _toa_simulation_metrics(
                    plan,
                    obs,
                    run_atmo,
                    monthly_prior,
                    label="monthly",
                )
                direct_run_metrics, direct_run_toa = _toa_simulation_metrics(
                    plan,
                    obs,
                    run_atmo,
                    direct_prior,
                    label="direct",
                )
                run_toa_ds = xr.merge([monthly_run_toa, direct_run_toa], compat="override")
                run_toa_path = output_dir / "toa_comparison_run_solution.nc"
                run_toa_ds.to_netcdf(run_toa_path)
                summary["toa_metrics_run_solution"] = {
                    "monthly": monthly_run_metrics,
                    "direct": direct_run_metrics,
                }

                for band_name in ("B02", "B04"):
                    observed_name = f"observed_toa_{band_name}"
                    monthly_name = f"simulated_toa_monthly_{band_name}"
                    direct_name = f"simulated_toa_direct_{band_name}"
                    if (
                        observed_name in run_toa_ds.data_vars
                        and monthly_name in run_toa_ds.data_vars
                        and direct_name in run_toa_ds.data_vars
                    ):
                        _plot_toa_residuals(
                            run_toa_ds[observed_name],
                            run_toa_ds[monthly_name],
                            run_toa_ds[direct_name],
                            band_name=band_name,
                            output_path=plots_dir / f"toa_run_solution_{band_name}.png",
                        )
            except Exception as exc:
                LOGGER.warning("Skipping run-solution TOA diagnostics: %s", exc)
                summary["toa_metrics_run_solution"] = {
                    "skipped": True,
                    "reason": str(exc),
                }

    LOGGER.info("Computing TOA comparison under atmospheric prior")
    try:
        monthly_prior_metrics_atmo, monthly_prior_toa = _toa_simulation_metrics(
            plan,
            obs,
            atmo_prior,
            monthly_prior,
            label="monthly",
        )
        direct_prior_metrics_atmo, direct_prior_toa = _toa_simulation_metrics(
            plan,
            obs,
            atmo_prior,
            direct_prior,
            label="direct",
        )
        atmo_toa_ds = xr.merge([monthly_prior_toa, direct_prior_toa], compat="override")
        atmo_toa_ds.to_netcdf(output_dir / "toa_comparison_atmo_prior.nc")
        summary["toa_metrics_atmo_prior"] = {
            "monthly": monthly_prior_metrics_atmo,
            "direct": direct_prior_metrics_atmo,
        }

        for band_name in ("B02", "B04"):
            observed_name = f"observed_toa_{band_name}"
            monthly_name = f"simulated_toa_monthly_{band_name}"
            direct_name = f"simulated_toa_direct_{band_name}"
            if (
                observed_name in atmo_toa_ds.data_vars
                and monthly_name in atmo_toa_ds.data_vars
                and direct_name in atmo_toa_ds.data_vars
            ):
                _plot_toa_residuals(
                    atmo_toa_ds[observed_name],
                    atmo_toa_ds[monthly_name],
                    atmo_toa_ds[direct_name],
                    band_name=band_name,
                    output_path=plots_dir / f"toa_atmo_prior_{band_name}.png",
                )
    except Exception as exc:
        LOGGER.warning("Skipping atmospheric-prior TOA diagnostics: %s", exc)
        summary["toa_metrics_atmo_prior"] = {
            "skipped": True,
            "reason": str(exc),
        }

    for band_name in ("B01", "B02", "B04"):
        if band_name not in [str(value) for value in np.asarray(monthly_prior.boa.coords["band"].values)]:
            continue
        _plot_prior_maps(
            monthly_prior.boa.sel(band=band_name, drop=True),
            direct_prior.boa.sel(band=band_name, drop=True),
            band_name=band_name,
            output_path=plots_dir / f"surface_prior_{band_name}.png",
        )

    _write_json(output_dir / "summary.json", summary)
    LOGGER.info("Experiment complete. Outputs written to %s", output_dir)


if __name__ == "__main__":
    main()
