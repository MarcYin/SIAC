"""
Pipeline orchestration for SIAC atmospheric correction.
"""

from __future__ import annotations

import inspect
import logging
import time
from collections.abc import Callable
from concurrent.futures import ThreadPoolExecutor
from concurrent.futures import TimeoutError as FuturesTimeoutError
from contextlib import nullcontext, suppress
from dataclasses import replace
from datetime import datetime
from pathlib import Path
from typing import TYPE_CHECKING, Any, TypeAlias, cast

import numpy as np
import xarray as xr
from numpy import typing as npt

from siac.algorithms.solver import build_solver_valid_mask
from siac.observability import (
    ExecutionObserver,
    bind_execution_observer,
    current_execution_observer,
    derive_summary_report_path,
    observe_stage,
    resolve_execution_observer,
)
from siac.runtime import (
    AOTScatterBandDiagnostics,
    AtmosphericState,
    CorrectionResult,
    GeometryAngles,
    MonthlyCompositeOutput,
    ObservationBundle,
    SolvedAtmosphere,
    SolverInputBundle,
    SurfacePrior,
)
from siac.runtime.models import copy_spatial_metadata_like
from siac.runtime.validation import (
    validate_atmospheric_state,
    validate_correction_result,
    validate_observation_bundle,
    validate_solved_atmosphere,
    validate_solver_input_bundle,
    validate_surface_prior,
)

if TYPE_CHECKING:
    from siac.domain.aoi import AOI
    from siac.domain.sensors import SensorConfig

logger = logging.getLogger(__name__)

Float32Array: TypeAlias = npt.NDArray[np.float32]

PreprocessorFn = Callable[[Path, Any], ObservationBundle]
AtmoPriorFn = Callable[
    [tuple[float, float, float, float], str, datetime, float],
    AtmosphericState,
]
SurfacePriorFn = Callable[[ObservationBundle, AtmosphericState | None, Any, float], SurfacePrior]
GridAssemblerFn = Callable[..., SolverInputBundle]
SolverFn = Callable[[SolverInputBundle, Any], SolvedAtmosphere]
CorrectorFn = Callable[..., CorrectionResult]
_OUTPUTS_WRITTEN_METADATA_KEY = "_siac_outputs_written"

# Module-level constants extracted from previously-hardcoded magic numbers.
_MAX_SCATTER_POINTS_PER_BAND = 4096
_DEFAULT_AUX_RESOLUTION_M = 500.0

_EXECUTION_KEYS = (
    "backend",
    "max_workers",
    "retries",
    "stage_timeout_s",
    "dashboard",
    "dashboard_address",
    "performance_report_path",
    "show_progress",
    "profiling_sample_interval_s",
    "progress_heartbeat_s",
)


def _execution_values(source: Any) -> dict[str, Any]:
    """Extract execution keys from dict/model-like objects."""
    if source is None:
        return {}
    if isinstance(source, dict):
        return {k: source[k] for k in _EXECUTION_KEYS if k in source}

    out: dict[str, Any] = {}
    for key in _EXECUTION_KEYS:
        if hasattr(source, key):
            out[key] = getattr(source, key)
    return out


def _resolve_execution_settings(
    config: Any,
    *,
    execution: Any | None,
    max_workers: int | None,
) -> dict[str, Any]:
    """Resolve execution settings from config + call overrides."""
    settings: dict[str, Any] = {
        "backend": "thread",
        "max_workers": 4,
        "retries": 2,
        "stage_timeout_s": None,
        "dashboard": False,
        "dashboard_address": None,
        "performance_report_path": None,
        "show_progress": False,
        "profiling_sample_interval_s": 5.0,
        "progress_heartbeat_s": 30.0,
    }

    settings.update(_execution_values(getattr(config, "execution", None)))
    settings.update(_execution_values(execution))
    if max_workers is not None:
        settings["max_workers"] = max_workers

    backend = str(settings.get("backend", "thread")).lower()
    if backend not in {"thread", "dask"}:
        raise ValueError(f"Unsupported execution backend: {backend!r}")
    settings["backend"] = backend

    workers = int(settings.get("max_workers", 4))
    if workers < 1:
        raise ValueError("max_workers must be >= 1")
    settings["max_workers"] = workers

    retries = int(settings.get("retries", 2))
    if retries < 0:
        raise ValueError("retries must be >= 0")
    settings["retries"] = retries

    timeout = settings.get("stage_timeout_s")
    if timeout is not None:
        timeout = float(timeout)
        if timeout <= 0:
            raise ValueError("stage_timeout_s must be > 0 when provided")
    settings["stage_timeout_s"] = timeout

    report_path = settings.get("performance_report_path")
    if report_path is not None:
        report_path = Path(report_path)
    settings["performance_report_path"] = report_path

    settings["dashboard"] = bool(settings.get("dashboard", False))
    settings["show_progress"] = bool(settings.get("show_progress", False))

    sample_interval = float(settings.get("profiling_sample_interval_s", 5.0))
    if sample_interval <= 0:
        raise ValueError("profiling_sample_interval_s must be > 0")
    settings["profiling_sample_interval_s"] = sample_interval

    heartbeat = float(settings.get("progress_heartbeat_s", 30.0))
    if heartbeat <= 0:
        raise ValueError("progress_heartbeat_s must be > 0")
    settings["progress_heartbeat_s"] = heartbeat
    return settings


def _aerosol_resolution(config: Any) -> float:
    solver_config = getattr(config, "solver", None)
    if solver_config is not None:
        return float(getattr(solver_config, "aerosol_resolution", 120.0))
    return 120.0


def _select_solver_bands_for_preload(sensor_config: SensorConfig) -> list[Any]:
    """Mirror M4 band-selection logic for LUT preloading hints."""
    default_selector = getattr(sensor_config, "default_aerosol_solver_bands", None)
    if callable(default_selector):
        return list(default_selector())

    range_selector = getattr(sensor_config, "select_bands_in_range", None)
    if callable(range_selector):
        bands = list(range_selector(400.0, 520.0))
        if bands:
            return bands

    return list(getattr(sensor_config, "bands", ())[:2])


def _should_skip_correction(config: Any) -> bool:
    output = getattr(config, "output", None)
    defaults = getattr(output, "defaults", None)
    if defaults is None:
        return False
    return bool(getattr(defaults, "skip_correction", False))


def _should_capture_aot_scatter(config: Any) -> bool:
    output = getattr(config, "output", None)
    defaults = getattr(output, "defaults", None)
    if defaults is None:
        return True
    return bool(getattr(defaults, "include_auxiliary", True))


def _sample_scatter_values(values: np.ndarray, *, max_points: int) -> Float32Array:
    if values.size <= max_points:
        return cast("Float32Array", values.astype(np.float32, copy=False))
    indices = np.linspace(0, values.size - 1, max_points, dtype=np.int64)
    return cast("Float32Array", values[indices].astype(np.float32, copy=False))


def _select_band_slice(
    data: xr.DataArray,
    *,
    band_name: str,
    band_index: int,
) -> xr.DataArray | None:
    if "band" not in data.dims:
        return data
    band_coord = data.coords.get("band")
    if band_coord is not None:
        band_values = [str(value) for value in np.asarray(band_coord.values).tolist()]
        if band_name in band_values:
            return cast("xr.DataArray", data.sel(band=band_name, drop=True))
        if np.asarray(band_coord.values).dtype.kind in {"U", "S", "O"}:
            return None
    return cast("xr.DataArray", data.isel(band=band_index, drop=True))


def _finite_diagnostic_field(
    values: xr.DataArray,
    fallback: xr.DataArray,
) -> xr.DataArray:
    source = np.asarray(values.values, dtype=np.float32)
    if np.all(np.isfinite(source)):
        return values

    filled = source.copy()
    missing = ~np.isfinite(filled)
    if fallback.shape == values.shape:
        fallback_values = np.asarray(fallback.values, dtype=np.float32)
        fallback_finite = missing & np.isfinite(fallback_values)
        filled[fallback_finite] = fallback_values[fallback_finite]
        missing = ~np.isfinite(filled)

    if np.any(missing):
        finite_values = filled[np.isfinite(filled)]
        fill_value = float(np.mean(finite_values)) if finite_values.size else 0.0
        filled[missing] = np.float32(fill_value)

    return xr.DataArray(
        filled,
        dims=values.dims,
        coords=values.coords,
        attrs=values.attrs,
        name=values.name,
    )


def _build_aot_scatter_diagnostics(
    solver_inputs: SolverInputBundle,
    solved: SolvedAtmosphere,
    *,
    max_points_per_band: int = _MAX_SCATTER_POINTS_PER_BAND,
) -> tuple[AOTScatterBandDiagnostics, ...]:
    valid_mask = build_solver_valid_mask(
        solver_inputs.cloud_mask,
        solver_inputs.toa,
        solver_inputs.surface_prior,
        sharp_transition_mask=solver_inputs.sharp_transition_mask,
    ).values.astype(bool)
    atmo_state = solved.atmo_state
    atmo_finite_mask = (
        np.isfinite(atmo_state.aot.values)
        & np.isfinite(atmo_state.tcwv.values)
        & np.isfinite(atmo_state.tco3.values)
        & np.isfinite(atmo_state.elevation.values)
    )
    valid_mask = valid_mask & atmo_finite_mask
    diagnostic_atmo_state = AtmosphericState(
        aot=_finite_diagnostic_field(atmo_state.aot, solver_inputs.atmo_prior.aot),
        tcwv=_finite_diagnostic_field(atmo_state.tcwv, solver_inputs.atmo_prior.tcwv),
        tco3=_finite_diagnostic_field(atmo_state.tco3, solver_inputs.atmo_prior.tco3),
        elevation=_finite_diagnostic_field(atmo_state.elevation, solver_inputs.atmo_prior.elevation),
        aot_unc=atmo_state.aot_unc,
        tcwv_unc=atmo_state.tcwv_unc,
        tco3_unc=atmo_state.tco3_unc,
    )
    diagnostics: list[AOTScatterBandDiagnostics] = []

    for band_index, band in enumerate(solver_inputs.bands):
        toa_band = _select_band_slice(solver_inputs.toa, band_name=band.name, band_index=band_index)
        surface_band = _select_band_slice(
            solver_inputs.surface_prior.boa,
            band_name=band.name,
            band_index=band_index,
        )
        if toa_band is None or surface_band is None:
            continue
        coeffs = solver_inputs.rt_model.compute_coefficients(
            solver_inputs.geometry,
            diagnostic_atmo_state,
            band,
            False,
        )
        simulated_toa = coeffs.simulate_toa(surface_band)

        band_valid = (
            valid_mask
            & np.isfinite(surface_band.values)
            & np.isfinite(toa_band.values)
            & np.isfinite(simulated_toa.values)
        )
        if not np.any(band_valid):
            continue

        surface_values = surface_band.values[band_valid].astype(np.float32, copy=False)
        observed_values = toa_band.values[band_valid].astype(np.float32, copy=False)
        simulated_values = simulated_toa.values[band_valid].astype(np.float32, copy=False)

        if surface_values.size > max_points_per_band:
            order = np.argsort(surface_values, kind="mergesort")
            surface_values = surface_values[order]
            observed_values = observed_values[order]
            simulated_values = simulated_values[order]
        diagnostics.append(
            AOTScatterBandDiagnostics(
                band_name=band.name,
                surface_reflectance=_sample_scatter_values(surface_values, max_points=max_points_per_band),
                observed_toa=_sample_scatter_values(observed_values, max_points=max_points_per_band),
                simulated_toa=_sample_scatter_values(simulated_values, max_points=max_points_per_band),
                total_valid_count=int(np.count_nonzero(band_valid)),
            )
        )

    return tuple(diagnostics)


def _geometry_for_atmo_grid(
    geometry: GeometryAngles,
    atmo: AtmosphericState,
) -> GeometryAngles:
    from siac.geo.resample import resample_field_to_template, shares_template_grid

    template = atmo.aot
    if all(
        shares_template_grid(field, template)
        for field in (geometry.sza, geometry.saa, geometry.vza, geometry.vaa)
    ):
        return geometry
    return GeometryAngles(
        sza=resample_field_to_template(geometry.sza, template),
        saa=resample_field_to_template(geometry.saa, template),
        vza=resample_field_to_template(geometry.vza, template),
        vaa=resample_field_to_template(geometry.vaa, template),
    )


def _maybe_submit_lut_preload(
    executor: Any,
    rt_model: Any,
    *,
    obs: ObservationBundle,
    atmo: AtmosphericState,
    retries: int,
    observer_id: str | None = None,
) -> Any | None:
    """Start LUT scene preload task when backend supports it."""
    preload_fn = getattr(rt_model, "preload_scene_subset", None)
    if not callable(preload_fn):
        return None

    bands = _select_solver_bands_for_preload(obs.sensor_config)
    preload_geometry = _geometry_for_atmo_grid(obs.geometry, atmo)
    logger.info(
        "Starting LUT preload in parallel on the atmospheric grid (scene subset + %d band grid%s).",
        len(bands),
        "" if len(bands) == 1 else "s",
    )
    observer = resolve_execution_observer(observer_id)
    if observer is not None:
        observer.increment_counter("lut_preload_started", stage="LUT.preload")
        observer.emit(
            "progress",
            stage="LUT.preload",
            message="Submitting LUT preload task.",
            band_count=len(bands),
            geometry_shape=tuple(preload_geometry.sza.shape),
        )
    return executor.submit(
        _call_with_retries,
        preload_fn,
        (preload_geometry, atmo, bands),
        retries=retries,
        stage_name="LUT.preload",
        observer_id=observer_id,
    )


def _surface_prior_requires_atmo(provider: SurfacePriorFn) -> bool:
    """Return whether the surface-prior provider depends on the M2 result."""
    return bool(getattr(provider, "requires_atmo_prior", False))


def _surface_template(data: xr.DataArray) -> xr.DataArray:
    if "band" in data.dims:
        band_coord = data.coords["band"].values[0] if "band" in data.coords else 0
        return data.sel(band=band_coord, drop=True)
    return data


def _band_name(value: object, index: int) -> str:
    if hasattr(value, "item"):
        with suppress(Exception):
            value = value.item()
    name = getattr(value, "name", None)
    if isinstance(name, str) and name:
        return name
    text = str(value)
    return text if text else f"band_{index + 1:02d}"


def _banded_dataarray_to_dataset(
    data: xr.DataArray,
    *,
    default_name: str,
    template: xr.DataArray,
) -> xr.Dataset:
    if "band" not in data.dims:
        return xr.Dataset({default_name: copy_spatial_metadata_like(data, template)})

    band_values = data.coords["band"].values if "band" in data.coords else np.arange(data.sizes["band"])
    return xr.Dataset(
        {
            _band_name(band, index): copy_spatial_metadata_like(
                data.sel(band=band, drop=True),
                template,
            )
            for index, band in enumerate(band_values)
        }
    )


def _monthly_composite_outputs(
    composites: tuple[Any, ...],
    *,
    template: xr.DataArray,
) -> dict[str, MonthlyCompositeOutput] | None:
    if not composites:
        return None

    outputs: dict[str, MonthlyCompositeOutput] = {}
    for composite in composites:
        label = f"{int(composite.year):04d}_{int(composite.month):02d}"
        outputs[label] = MonthlyCompositeOutput(
            reflectance=_banded_dataarray_to_dataset(
                composite.reflectance,
                default_name="reflectance",
                template=template,
            ),
            quality=copy_spatial_metadata_like(composite.quality.astype(np.float32), template),
            sample_index=copy_spatial_metadata_like(composite.sample_index.astype(np.int16), template),
        )
    return outputs


def _call_grid_assembler(
    grid_assembler: GridAssemblerFn,
    obs: ObservationBundle,
    atmo: AtmosphericState,
    surface: SurfacePrior,
    rt_model: Any,
    *,
    aerosol_resolution_m: float,
    sharp_transition_filter: Any | None = None,
) -> SolverInputBundle:
    """Call the grid assembler with a standardised interface.

    The assembler is expected to accept the signature::

        (obs, atmo, surface, rt_model, *, aerosol_resolution_m, sharp_transition_filter=None)

    For backwards compatibility with assemblers that do not accept the
    ``sharp_transition_filter`` keyword, a single fallback attempt is made
    without it.
    """
    try:
        return grid_assembler(
            obs, atmo, surface, rt_model,
            aerosol_resolution_m=aerosol_resolution_m,
            sharp_transition_filter=sharp_transition_filter,
        )
    except TypeError:
        logger.debug(
            "Grid assembler rejected sharp_transition_filter keyword; retrying without it",
            exc_info=True,
        )
        return grid_assembler(
            obs, atmo, surface, rt_model,
            aerosol_resolution_m=aerosol_resolution_m,
        )


def _call_with_retries(
    fn: Callable[..., Any],
    args: tuple[Any, ...],
    *,
    retries: int,
    stage_name: str,
    observer_id: str | None = None,
) -> Any:
    """Call a stage function with bounded retries."""
    observer = resolve_execution_observer(observer_id)
    max_attempts = retries + 1
    with bind_execution_observer(observer):
        for attempt in range(max_attempts):
            attempt_index = attempt + 1
            if observer is not None:
                observer.increment_counter("attempts_total", stage=stage_name)
            t0 = time.monotonic()
            try:
                with observe_stage(
                    stage_name,
                    details={"attempt": attempt_index, "max_attempts": max_attempts},
                ):
                    return fn(*args)
            except Exception as exc:
                duration = time.monotonic() - t0
                if observer is not None:
                    observer.increment_counter("attempts_failed", stage=stage_name)
                if attempt >= retries:
                    raise
                if observer is not None:
                    observer.record_retry(
                        stage=stage_name,
                        attempt=attempt_index,
                        max_attempts=max_attempts,
                        duration_s=duration,
                        error_type=type(exc).__name__,
                        error_message=str(exc),
                    )
                logger.warning(
                    "%s failed on attempt %d/%d after %.2fs (%s: %s); retrying",
                    stage_name,
                    attempt_index,
                    max_attempts,
                    duration,
                    type(exc).__name__,
                    exc,
                )

    raise RuntimeError(f"Unreachable retry state for {stage_name}")


def _prepare_observation(
    input_path: Path,
    aoi: AOI | None,
    config: Any,
    *,
    preprocessor: PreprocessorFn,
) -> tuple[ObservationBundle, tuple[float, float, float, float], str, datetime, float]:
    logger.info("M1: Preprocessing satellite data...")
    t0 = time.monotonic()
    with observe_stage("M1.preprocessing", details={"input_path": str(input_path)}):
        obs = preprocessor(input_path, aoi)
        validate_observation_bundle(obs)
    elapsed = time.monotonic() - t0
    band_count = len(list(obs.toa.data_vars))
    logger.info(
        "M1: Preprocessing complete (%.2fs), %d bands.",
        elapsed,
        band_count,
    )

    bounds = obs.bounds
    crs = obs.crs
    obs_time = obs.metadata["observation_time"]
    resolution = _aerosol_resolution(config)
    observer = current_execution_observer()
    if observer is not None:
        observer.emit(
            "progress",
            stage="M1.preprocessing",
            message="Observation bundle ready.",
            duration_s=elapsed,
            band_count=band_count,
            bounds=bounds,
            crs=crs,
            observation_time=obs_time,
            aerosol_resolution_m=resolution,
        )
    return obs, bounds, crs, obs_time, resolution


def _open_correction_output_stream(
    output_writer: Any | None,
    output_path: Path | str | None,
    *,
    metadata: dict[str, Any],
) -> Any | None:
    if output_writer is None or output_path is None:
        return None
    opener = getattr(output_writer, "open_correction_boa_stream", None)
    if not callable(opener):
        return None
    return opener(output_path, metadata=metadata)


def _call_corrector(
    corrector: CorrectorFn,
    obs: ObservationBundle,
    solved: SolvedAtmosphere,
    rt_model: Any,
    *,
    output_stream: Any | None,
) -> CorrectionResult:
    if output_stream is None:
        return corrector(obs, solved, rt_model)
    try:
        signature = inspect.signature(corrector)
    except (TypeError, ValueError):
        return corrector(obs, solved, rt_model)
    if "output_stream" not in signature.parameters and not any(
        param.kind == inspect.Parameter.VAR_KEYWORD
        for param in signature.parameters.values()
    ):
        return corrector(obs, solved, rt_model)
    return corrector(obs, solved, rt_model, output_stream=output_stream)


def _run_tail(
    obs: ObservationBundle,
    atmo: AtmosphericState,
    surface: SurfacePrior,
    config: Any,
    *,
    grid_assembler: GridAssemblerFn,
    solver: SolverFn,
    corrector: CorrectorFn,
    rt_model: Any,
    output_path: Path | str | None = None,
    output_writer: Any | None = None,
) -> CorrectionResult:
    validate_atmospheric_state(atmo)
    validate_surface_prior(surface)
    aerosol_resolution = _aerosol_resolution(config)
    solver_config = getattr(config, "solver", None)

    t0 = time.monotonic()
    logger.info("M4: Assembling solver grids...")
    with observe_stage("M4.grid_assembly"):
        solver_inputs = _call_grid_assembler(
            grid_assembler,
            obs,
            atmo,
            surface,
            rt_model,
            aerosol_resolution_m=aerosol_resolution,
            sharp_transition_filter=getattr(solver_config, "sharp_transition_filter", None),
        )
    validate_solver_input_bundle(solver_inputs)
    logger.info("M4: Grid assembly complete (%.2fs).", time.monotonic() - t0)

    # ------------------------------------------------------------------
    # skip_correction: write available auxiliary data without running the
    # solver (M5) or atmospheric correction (M6).
    # ------------------------------------------------------------------
    if _should_skip_correction(config):
        logger.info("skip_correction enabled — skipping M5 solver and M6 correction.")
        surface_template = _surface_template(solver_inputs.surface_prior.boa)
        result = CorrectionResult(
            boa=xr.Dataset(),
            boa_unc=None,
            aot=solver_inputs.atmo_prior.aot.astype(np.float32),
            tcwv=solver_inputs.atmo_prior.tcwv.astype(np.float32),
            cloud_mask=solver_inputs.cloud_mask,
            surface_prior=_banded_dataarray_to_dataset(
                solver_inputs.surface_prior.boa,
                default_name="surface_prior",
                template=surface_template,
            ),
            surface_prior_unc=_banded_dataarray_to_dataset(
                solver_inputs.surface_prior.boa_unc,
                default_name="surface_prior_unc",
                template=surface_template,
            ),
            solver_qa=None,
            monthly_composites=_monthly_composite_outputs(
                tuple(getattr(surface, "monthly_composites", ())),
                template=surface_template,
            ),
            metadata={
                **obs.metadata,
                "skip_correction": True,
                "sensor_config": obs.sensor_config,
                "geometry": obs.geometry,
                "crs": obs.crs,
                "bounds": obs.bounds,
            },
        )
        logger.info("Pipeline complete (auxiliary-only mode).")
        return result

    t0 = time.monotonic()
    logger.info("M5: Solving for aerosol parameters...")
    with observe_stage("M5.solver"):
        solved = solver(solver_inputs, config)
    validate_solved_atmosphere(solved)
    logger.info("M5: Solver complete (%.2fs).", time.monotonic() - t0)
    observer = current_execution_observer()
    if observer is not None:
        observer.emit(
            "progress",
            stage="M5.solver",
            message="Solver finished.",
            converged=solved.converged,
            iterations=solved.n_iterations,
            cost_final=solved.cost_final,
        )

    t0 = time.monotonic()
    logger.info("M6: Applying atmospheric correction...")
    with observe_stage("M6.correction"):
        output_stream = _open_correction_output_stream(
            output_writer,
            output_path,
            metadata=obs.metadata,
        )
        corrected = _call_corrector(
            corrector,
            obs,
            solved,
            rt_model,
            output_stream=output_stream,
        )
        diagnostics = corrected.diagnostics
        if _should_capture_aot_scatter(config):
            try:
                scatter = _build_aot_scatter_diagnostics(solver_inputs, solved)
            except Exception:
                logger.exception("Failed to build aerosol-solver scatter diagnostics.")
            else:
                diagnostics = replace(corrected.diagnostics, aot_scatter_plots=scatter)
        surface_template = _surface_template(solver_inputs.surface_prior.boa)
        result = CorrectionResult(
            boa=corrected.boa,
            boa_unc=corrected.boa_unc,
            aot=corrected.aot,
            tcwv=corrected.tcwv,
            cloud_mask=corrected.cloud_mask,
            surface_prior=_banded_dataarray_to_dataset(
                solver_inputs.surface_prior.boa,
                default_name="surface_prior",
                template=surface_template,
            ),
            surface_prior_unc=_banded_dataarray_to_dataset(
                solver_inputs.surface_prior.boa_unc,
                default_name="surface_prior_unc",
                template=surface_template,
            ),
            solver_qa=corrected.solver_qa,
            monthly_composites=_monthly_composite_outputs(
                tuple(getattr(surface, "monthly_composites", ())),
                template=surface_template,
            ),
            metadata={
                **obs.metadata,
                **corrected.metadata,
                "sensor_config": obs.sensor_config,
                "geometry": obs.geometry,
                "crs": obs.crs,
                "bounds": obs.bounds,
            },
            diagnostics=diagnostics,
        )
        validate_correction_result(result)
        if output_stream is not None and getattr(output_stream, "has_written", False):
            output_stream.finish(result)
            result = replace(
                result,
                metadata={
                    **result.metadata,
                    _OUTPUTS_WRITTEN_METADATA_KEY: True,
                },
            )
    logger.info("M6: Correction complete (%.2fs).", time.monotonic() - t0)

    if observer is not None:
        observer.emit(
            "progress",
            stage="M6.correction",
            message="Correction result ready.",
            boa_bands=len(list(result.boa.data_vars)),
        )
    logger.info("Pipeline complete.")
    return result


def _fetch_priors(
    *,
    submit_fn: Callable[..., Any],
    lut_submit_fn: Callable[..., Any] | None,
    obs: ObservationBundle,
    config: Any,
    atmo_provider: AtmoPriorFn,
    surface_prior_provider: SurfacePriorFn,
    rt_model: Any,
    settings: dict[str, Any],
    backend_label: str,
    timeout_errors: tuple[type[BaseException], ...] = (FuturesTimeoutError,),
) -> tuple[AtmosphericState, SurfacePrior]:
    """Fetch atmospheric and surface priors using the provided submit callable.

    ``submit_fn(fn, args, *, retries, stage_name, observer_id)`` must return a
    future whose ``.result(timeout=...)`` / ``.cancel()`` interface matches
    ``concurrent.futures.Future``.
    """
    observer = current_execution_observer()
    observer_id = observer.run_id if observer is not None else None
    timeout = settings["stage_timeout_s"]
    retries = settings["retries"]
    resolution = _aerosol_resolution(config)
    requires_atmo = _surface_prior_requires_atmo(surface_prior_provider)

    f_m2 = submit_fn(
        _call_with_retries,
        atmo_provider,
        (obs.bounds, obs.crs, obs.metadata["observation_time"], resolution),
        retries=retries,
        stage_name="M2.atmospheric_prior",
        observer_id=observer_id,
    )
    if observer is not None:
        observer.increment_counter("m2_started", stage="M2.atmospheric_prior")
        observer.emit("progress", stage="M2.atmospheric_prior", message="Atmospheric prior submitted.")

    f_m3 = None
    if not requires_atmo:
        logger.info("M2+M3: Fetching atmospheric & surface priors...")
        f_m3 = submit_fn(
            _call_with_retries,
            surface_prior_provider,
            (obs, None, rt_model, resolution),
            retries=retries,
            stage_name="M3.surface_prior",
            observer_id=observer_id,
        )
        if observer is not None:
            observer.increment_counter("m3_started", stage="M3.surface_prior")
            observer.emit(
                "progress",
                stage="M3.surface_prior",
                message="Surface prior submitted in parallel with atmospheric prior.",
            )
    else:
        logger.info("M2: Fetching atmospheric prior...")

    f_lut = None
    try:
        atmo = f_m2.result(timeout=timeout)
        if observer is not None:
            observer.increment_counter("m2_done", stage="M2.atmospheric_prior")
            observer.emit(
                "progress",
                stage="M2.atmospheric_prior",
                message="Atmospheric prior ready.",
            )
        if f_m3 is None:
            logger.info("M3: Deriving surface prior from observation + atmospheric prior...")
            f_m3 = submit_fn(
                _call_with_retries,
                surface_prior_provider,
                (obs, atmo, rt_model, resolution),
                retries=retries,
                stage_name="M3.surface_prior",
                observer_id=observer_id,
            )
            if observer is not None:
                observer.increment_counter("m3_started", stage="M3.surface_prior")
                observer.emit(
                    "progress",
                    stage="M3.surface_prior",
                    message="Surface prior submitted after atmospheric prior.",
                )
        _submit_fn = lut_submit_fn if lut_submit_fn is not None else submit_fn
        f_lut = _maybe_submit_lut_preload(
            _SubmitAdapter(_submit_fn),
            rt_model,
            obs=obs,
            atmo=atmo,
            retries=retries,
            observer_id=observer_id,
        )
        surface = f_m3.result(timeout=timeout)
        if observer is not None:
            observer.increment_counter("m3_done", stage="M3.surface_prior")
            observer.emit(
                "progress",
                stage="M3.surface_prior",
                message="Surface prior ready.",
            )
        if f_lut is not None:
            try:
                f_lut.result(timeout=timeout)
                if observer is not None:
                    observer.increment_counter("lut_preload_done", stage="LUT.preload")
                    observer.emit(
                        "progress",
                        stage="LUT.preload",
                        message="LUT preload complete.",
                    )
            except FuturesTimeoutError:
                if observer is not None and timeout is not None:
                    observer.record_timeout(
                        stage="LUT.preload",
                        timeout_s=float(timeout),
                        backend=backend_label,
                    )
                logger.warning(
                    "LUT preload timed out after %.1fs; proceeding with on-demand LUT reads.",
                    float(timeout),
                )
            except Exception as exc:
                if observer is not None:
                    observer.record_error(
                        stage="LUT.preload",
                        error_type=type(exc).__name__,
                        error_message=str(exc),
                    )
                logger.warning(
                    "LUT preload failed (%s); proceeding with on-demand LUT reads.",
                    exc,
                )
    except timeout_errors as exc:
        f_m2.cancel()
        if f_m3 is not None:
            f_m3.cancel()
        if f_lut is not None:
            f_lut.cancel()
        if observer is not None and timeout is not None:
            observer.record_timeout(
                stage="M2/M3",
                timeout_s=float(timeout),
                backend=backend_label,
            )
        raise TimeoutError(
            f"M2/M3 timed out after {float(timeout):.1f}s ({backend_label} backend)"
        ) from exc
    except Exception:
        f_m2.cancel()
        if f_m3 is not None:
            f_m3.cancel()
        if f_lut is not None:
            f_lut.cancel()
        raise

    return atmo, surface


class _SubmitAdapter:
    """Adapt a bare submit callable to look like an Executor for _maybe_submit_lut_preload."""

    def __init__(self, submit_fn: Callable[..., Any]) -> None:
        self._submit_fn = submit_fn

    def submit(self, fn: Callable[..., Any], *args: Any, **kwargs: Any) -> Any:
        return self._submit_fn(fn, *args, **kwargs)


def _run_pipeline_thread(
    input_path: Path,
    aoi: AOI | None,
    config: Any,
    *,
    preprocessor: PreprocessorFn,
    atmo_provider: AtmoPriorFn,
    surface_prior_provider: SurfacePriorFn,
    grid_assembler: GridAssemblerFn,
    solver: SolverFn,
    corrector: CorrectorFn,
    rt_model: Any,
    settings: dict[str, Any],
    output_path: Path | str | None = None,
    output_writer: Any | None = None,
) -> CorrectionResult:
    obs, bounds, crs, obs_time, resolution = _prepare_observation(
        input_path,
        aoi,
        config,
        preprocessor=preprocessor,
    )

    with ThreadPoolExecutor(max_workers=settings["max_workers"]) as executor:
        atmo, surface = _fetch_priors(
            submit_fn=executor.submit,
            lut_submit_fn=None,
            obs=obs,
            config=config,
            atmo_provider=atmo_provider,
            surface_prior_provider=surface_prior_provider,
            rt_model=rt_model,
            settings=settings,
            backend_label="thread",
        )

    return _run_tail(
        obs,
        atmo,
        surface,
        config,
        grid_assembler=grid_assembler,
        solver=solver,
        corrector=corrector,
        rt_model=rt_model,
        output_path=output_path,
        output_writer=output_writer,
    )


def _run_pipeline_dask(
    input_path: Path,
    aoi: AOI | None,
    config: Any,
    *,
    preprocessor: PreprocessorFn,
    atmo_provider: AtmoPriorFn,
    surface_prior_provider: SurfacePriorFn,
    grid_assembler: GridAssemblerFn,
    solver: SolverFn,
    corrector: CorrectorFn,
    rt_model: Any,
    settings: dict[str, Any],
    output_path: Path | str | None = None,
    output_writer: Any | None = None,
) -> CorrectionResult:
    try:
        from dask.distributed import (  # type: ignore[import-not-found]
            Client,
            LocalCluster,
            performance_report,
        )
        from dask.distributed import TimeoutError as DaskTimeoutError
    except Exception as exc:
        raise RuntimeError(
            "Dask backend requested but dask.distributed is not installed. "
            "Install dask/distributed or set execution.backend='thread'."
        ) from exc

    cluster_kwargs: dict[str, Any] = {
        "n_workers": settings["max_workers"],
        "threads_per_worker": 1,
        "processes": False,
        "dashboard_address": settings["dashboard_address"] if settings["dashboard"] else None,
    }
    observer = current_execution_observer()

    with LocalCluster(**cluster_kwargs) as cluster, Client(cluster) as client:
        if settings["show_progress"] and getattr(client, "dashboard_link", None):
            logger.info("Dask dashboard: %s", client.dashboard_link)
            if observer is not None:
                observer.emit(
                    "progress",
                    stage="dask.cluster",
                    message="Dask dashboard available.",
                    dashboard_link=client.dashboard_link,
                )

        report_ctx = nullcontext()
        report_path = settings["performance_report_path"]
        if report_path is not None:
            report_path.parent.mkdir(parents=True, exist_ok=True)
            report_ctx = performance_report(filename=str(report_path))

        with report_ctx:
            obs, bounds, crs, obs_time, resolution = _prepare_observation(
                input_path,
                aoi,
                config,
                preprocessor=preprocessor,
            )

            preload_executor = ThreadPoolExecutor(max_workers=1)
            try:
                atmo, surface = _fetch_priors(
                    submit_fn=client.submit,
                    lut_submit_fn=preload_executor.submit,
                    obs=obs,
                    config=config,
                    atmo_provider=atmo_provider,
                    surface_prior_provider=surface_prior_provider,
                    rt_model=rt_model,
                    settings=settings,
                    backend_label="dask",
                    timeout_errors=(FuturesTimeoutError, DaskTimeoutError),
                )
            finally:
                preload_executor.shutdown(wait=False, cancel_futures=False)

    return _run_tail(
        obs,
        atmo,
        surface,
        config,
        grid_assembler=grid_assembler,
        solver=solver,
        corrector=corrector,
        rt_model=rt_model,
        output_path=output_path,
        output_writer=output_writer,
    )


def run_pipeline(
    input_path: Path,
    aoi: AOI | None,
    config: Any,
    *,
    preprocessor: PreprocessorFn,
    atmo_provider: AtmoPriorFn,
    surface_prior_provider: SurfacePriorFn,
    grid_assembler: GridAssemblerFn,
    solver: SolverFn,
    corrector: CorrectorFn,
    rt_model: Any,
    max_workers: int | None = None,
    execution: Any | None = None,
    observer: ExecutionObserver | None = None,
    output_path: Path | str | None = None,
    output_writer: Any | None = None,
) -> CorrectionResult:
    """Orchestrate module execution with a selectable execution backend."""
    settings = _resolve_execution_settings(
        config,
        execution=execution,
        max_workers=max_workers,
    )
    summary_path = derive_summary_report_path(
        settings["performance_report_path"],
        backend=settings["backend"],
    )
    managed_observer = False
    if observer is None and (summary_path is not None or settings["show_progress"]):
        observer = ExecutionObserver(
            backend=settings["backend"],
            input_path=input_path,
            summary_path=summary_path,
            show_progress=settings["show_progress"],
            sample_interval_s=settings["profiling_sample_interval_s"],
            heartbeat_interval_s=settings["progress_heartbeat_s"],
            configured_report_path=settings["performance_report_path"],
        )
        managed_observer = True

    if observer is not None:
        observer.start()

    try:
        with bind_execution_observer(observer):
            if observer is not None:
                observer.emit(
                    "progress",
                    stage="pipeline",
                    message="Pipeline execution started.",
                    backend=settings["backend"],
                    max_workers=settings["max_workers"],
                )
            if settings["backend"] == "dask":
                result = _run_pipeline_dask(
                    input_path,
                    aoi,
                    config,
                    preprocessor=preprocessor,
                    atmo_provider=atmo_provider,
                    surface_prior_provider=surface_prior_provider,
                    grid_assembler=grid_assembler,
                    solver=solver,
                    corrector=corrector,
                    rt_model=rt_model,
                    settings=settings,
                    output_path=output_path,
                    output_writer=output_writer,
                )
            else:
                result = _run_pipeline_thread(
                    input_path,
                    aoi,
                    config,
                    preprocessor=preprocessor,
                    atmo_provider=atmo_provider,
                    surface_prior_provider=surface_prior_provider,
                    grid_assembler=grid_assembler,
                    solver=solver,
                    corrector=corrector,
                    rt_model=rt_model,
                    settings=settings,
                    output_path=output_path,
                    output_writer=output_writer,
                )
    except Exception as exc:
        if observer is not None and managed_observer:
            status = "timeout" if isinstance(exc, TimeoutError) else "error"
            observer.finish(status, error=exc)
        raise

    if observer is not None and managed_observer:
        observer.finish("success")
    return result
