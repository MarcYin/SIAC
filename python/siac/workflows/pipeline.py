"""
Pipeline orchestration for SIAC atmospheric correction.
"""

from __future__ import annotations

import inspect
import logging
import time
from collections.abc import Callable
from concurrent.futures import TimeoutError as FuturesTimeoutError
from dataclasses import replace
from datetime import datetime
from pathlib import Path
from typing import TYPE_CHECKING, Any, cast

import xarray as xr

from siac.adapters.data.water_mask import DEFAULT_WATER_MASK_VRT_URL
from siac.errors import ValidationError as SIACValidationError
from siac.observability import (
    ExecutionObserver,
    bind_execution_observer,
    current_execution_observer,
    derive_summary_report_path,
    observe_stage,
    resolve_execution_observer,
)
from siac.runtime import (
    AtmosphericState,
    CorrectionResult,
    GeometryAngles,
    ObservationBundle,
    SolvedAtmosphere,
    SolverInputBundle,
    SurfacePrior,
)
from siac.runtime.validation import (
    validate_atmospheric_state,
    validate_correction_result,
    validate_observation_bundle,
    validate_solved_atmosphere,
    validate_solver_input_bundle,
    validate_surface_prior,
)
from siac.workflows._pipeline_config import (
    _DEFAULT_AUX_RESOLUTION_M as _DEFAULT_AUX_RESOLUTION_M,
)
from siac.workflows._pipeline_config import (
    _EXECUTION_KEYS as _EXECUTION_KEYS,
)
from siac.workflows._pipeline_config import (
    PipelineExecutionSettings,
    _aerosol_resolution,
    _requested_solver_band_names,
    _resolve_execution_settings,
    _select_solver_bands_for_preload,
    _should_capture_aot_scatter,
    _should_skip_correction,
)
from siac.workflows._pipeline_config import (
    _config_aux_resolution as _config_aux_resolution,
)
from siac.workflows._pipeline_config import (
    _config_scatter_limit as _config_scatter_limit,
)
from siac.workflows._pipeline_config import (
    _execution_values as _execution_values,
)
from siac.workflows._pipeline_diagnostics import (
    build_aot_scatter_diagnostics as _build_aot_scatter_diagnostics,
)
from siac.workflows._pipeline_outputs import (
    banded_dataarray_to_dataset as _banded_dataarray_to_dataset,
)
from siac.workflows._pipeline_outputs import (
    monthly_composite_outputs as _monthly_composite_outputs,
)
from siac.workflows._pipeline_outputs import (
    surface_template as _surface_template,
)

if TYPE_CHECKING:
    from siac.domain.aoi import AOI

logger = logging.getLogger(__name__)

_NON_RETRYABLE_EXCEPTIONS = (
    SIACValidationError,
    ValueError,
    TypeError,
    AssertionError,
    AttributeError,
    KeyError,
    NotImplementedError,
)

PreprocessorFn = Callable[[Path, Any], ObservationBundle]
AtmoPriorFn = Callable[
    [tuple[float, float, float, float], str, datetime, float],
    AtmosphericState,
]
SurfacePriorFn = Callable[[ObservationBundle, AtmosphericState | None, Any, float], SurfacePrior]
GridAssemblerFn = Callable[..., SolverInputBundle]
SolverFn = Callable[[SolverInputBundle, Any], SolvedAtmosphere]
CorrectorFn = Callable[..., CorrectionResult]


def _is_retryable_exception(exc: BaseException) -> bool:
    return not isinstance(exc, _NON_RETRYABLE_EXCEPTIONS)


def _stage_timeout(settings: PipelineExecutionSettings, *stage_names: str) -> float | None:
    stage_timeouts = settings.get("stage_timeouts", {})
    if isinstance(stage_timeouts, dict):
        for stage_name in stage_names:
            timeout = stage_timeouts.get(stage_name)
            if timeout is not None:
                return float(timeout)
    timeout = settings.get("stage_timeout_s")
    return None if timeout is None else float(timeout)


_OUTPUTS_WRITTEN_METADATA_KEY = "_siac_outputs_written"


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
    requested_band_names: tuple[str, ...] | None,
    retries: int,
    observer_id: str | None = None,
) -> Any | None:
    """Start LUT scene preload task when backend supports it."""
    preload_fn = getattr(rt_model, "preload_scene_subset", None)
    if not callable(preload_fn):
        return None

    bands = _select_solver_bands_for_preload(
        obs.sensor_config,
        requested_band_names=requested_band_names,
    )
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


def _call_grid_assembler(
    grid_assembler: GridAssemblerFn,
    obs: ObservationBundle,
    atmo: AtmosphericState,
    surface: SurfacePrior,
    rt_model: Any,
    *,
    aerosol_resolution_m: float,
    sharp_transition_filter: Any | None = None,
    water_mask_path: str | Path | None = None,
    water_mask_cache_dir: str | Path | None = None,
    water_mask_buffer_pixels: int = 0,
    solver_band_names: tuple[str, ...] | None = None,
) -> SolverInputBundle:
    """Call the grid assembler with the current standardized interface."""
    return grid_assembler(
        obs,
        atmo,
        surface,
        rt_model,
        aerosol_resolution_m=aerosol_resolution_m,
        sharp_transition_filter=sharp_transition_filter,
        water_mask_path=water_mask_path,
        water_mask_cache_dir=water_mask_cache_dir,
        water_mask_buffer_pixels=water_mask_buffer_pixels,
        solver_band_names=solver_band_names,
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
                if attempt >= retries or not _is_retryable_exception(exc):
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


def _set_rt_observation_time(rt_model: Any, observation_time: datetime) -> None:
    setter = getattr(rt_model, "set_observation_time", None)
    if callable(setter):
        setter(observation_time)


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
        param.kind == inspect.Parameter.VAR_KEYWORD for param in signature.parameters.values()
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
    solver_config = getattr(getattr(config, "algorithms", None), "solver", None)
    paths_config = getattr(config, "paths", None)
    water_mask_path = getattr(paths_config, "water_mask", None) or DEFAULT_WATER_MASK_VRT_URL
    cache_root = getattr(paths_config, "cache_root", None)
    water_mask_cache_dir = (
        None if cache_root is None else Path(cache_root).expanduser() / "water-mask"
    )
    water_mask_buffer_pixels = int(getattr(solver_config, "water_mask_buffer_pixels", 0))
    solver_band_names = _requested_solver_band_names(config)

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
            water_mask_path=water_mask_path,
            water_mask_cache_dir=water_mask_cache_dir,
            water_mask_buffer_pixels=water_mask_buffer_pixels,
            solver_band_names=solver_band_names,
        )
    validate_solver_input_bundle(solver_inputs)
    logger.info("M4: Grid assembly complete (%.2fs).", time.monotonic() - t0)

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

    # ------------------------------------------------------------------
    # skip_correction: aerosol retrieval (M5) ran, but skip applying the
    # retrieved parameters for atmospheric correction (M6).
    # ------------------------------------------------------------------
    if _should_skip_correction(config):
        logger.info("skip_correction enabled — skipping M6 atmospheric correction.")
        surface_template = _surface_template(solver_inputs.surface_prior.boa)
        result = CorrectionResult(
            boa=xr.Dataset(),
            boa_unc=None,
            aot=solved.aot,
            tcwv=solved.tcwv,
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
            solver_qa=solved.qa,
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
        logger.info("Pipeline complete (retrieval-only mode, no BOA correction).")
        return result

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
            except (ValueError, KeyError, IndexError, TypeError) as exc:
                logger.warning(
                    "Failed to build aerosol-solver scatter diagnostics (%s: %s).",
                    type(exc).__name__,
                    exc,
                )
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
    settings: PipelineExecutionSettings,
    backend_label: str,
    timeout_errors: tuple[type[BaseException], ...] = (FuturesTimeoutError,),
) -> tuple[AtmosphericState, SurfacePrior]:
    from siac.workflows._pipeline_priors import fetch_priors

    return cast(
        "tuple[AtmosphericState, SurfacePrior]",
        fetch_priors(
            submit_fn=submit_fn,
            lut_submit_fn=lut_submit_fn,
            obs=obs,
            config=config,
            atmo_provider=atmo_provider,
            surface_prior_provider=surface_prior_provider,
            rt_model=rt_model,
            settings=settings,
            backend_label=backend_label,
            call_with_retries_fn=_call_with_retries,
            maybe_submit_lut_preload_fn=_maybe_submit_lut_preload,
            stage_timeout_fn=_stage_timeout,
            aerosol_resolution_fn=_aerosol_resolution,
            surface_prior_requires_atmo_fn=_surface_prior_requires_atmo,
            requested_solver_band_names_fn=_requested_solver_band_names,
            timeout_errors=timeout_errors,
        ),
    )


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
    settings: PipelineExecutionSettings,
    output_path: Path | str | None = None,
    output_writer: Any | None = None,
) -> CorrectionResult:
    from siac.workflows._pipeline_executors import PipelineExecutionContext, run_pipeline_thread

    return cast(
        "CorrectionResult",
        run_pipeline_thread(
            input_path,
            aoi,
            config,
            context=PipelineExecutionContext(
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
                prepare_observation=_prepare_observation,
                set_rt_observation_time=_set_rt_observation_time,
                fetch_priors=_fetch_priors,
                run_tail=_run_tail,
            ),
        ),
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
    settings: PipelineExecutionSettings,
    output_path: Path | str | None = None,
    output_writer: Any | None = None,
) -> CorrectionResult:
    from siac.workflows._pipeline_executors import PipelineExecutionContext, run_pipeline_dask

    return cast(
        "CorrectionResult",
        run_pipeline_dask(
            input_path,
            aoi,
            config,
            context=PipelineExecutionContext(
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
                prepare_observation=_prepare_observation,
                set_rt_observation_time=_set_rt_observation_time,
                fetch_priors=_fetch_priors,
                run_tail=_run_tail,
            ),
        ),
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
    """Run the SIAC pipeline with caller-supplied stage callables.

    This is the **advanced public extension point** for SIAC. Most users
    should call :func:`siac.api.process_sentinel2`, :func:`siac.api.siac_process`,
    or :class:`siac.api.SIAC` instead — those high-level entry points
    construct sensible defaults for every stage callable from the
    configuration object. Reach for ``run_pipeline`` only when you need to
    inject a custom RT backend, a non-standard surface prior, or otherwise
    swap a single stage of the M1-M6 pipeline without re-implementing
    the rest.

    REVIEW.md §3.6 originally flagged this function as having a 10+
    parameter signature with private-typed callables and no documented
    contract; it's now ratified as the canonical extension point. The
    callable types
    (:type:`PreprocessorFn`, :type:`AtmoPriorFn`,
    :type:`SurfacePriorFn`, :type:`GridAssemblerFn`, :type:`SolverFn`,
    :type:`CorrectorFn`) are exported from this module and form the
    public extension contract.

    Stages run in this order:

    1. ``preprocessor(input_path, aoi, config)``
       → :class:`siac.runtime.ObservationBundle` (M1)
    2. ``atmo_provider(observation, config)``
       → :class:`siac.runtime.AtmosphericState` (M2)
    3. ``surface_prior_provider(observation, atmo, rt_model, resolution)``
       → :class:`siac.runtime.SurfacePrior` (M3)
    4. ``grid_assembler(observation, atmo, surface_prior, ...)``
       → :class:`siac.runtime.SolverInputBundle` (M4)
    5. ``solver(inputs, config)``
       → :class:`siac.runtime.SolvedAtmosphere` (M5)
    6. ``corrector(observation, solved_atmo, rt_model, ...)``
       → :class:`siac.runtime.CorrectionResult` (M6)

    Args:
        input_path: SAFE directory for the scene to process.
        aoi: Optional clipping AOI; ``None`` means full-scene.
        config: A :class:`siac.config.SIACConfig` (or compatible) object.
        preprocessor: Stage M1 callable.
        atmo_provider: Stage M2 callable.
        surface_prior_provider: Stage M3 callable.
        grid_assembler: Stage M4 callable.
        solver: Stage M5 callable.
        corrector: Stage M6 callable.
        rt_model: RT model passed to M3/M5/M6 (typically built via
            :func:`siac.adapters.rt.build_rt_model`).
        max_workers: Optional override of ``config.runtime.execution.max_workers``.
        execution: Optional pre-resolved :class:`PipelineExecutionSettings`
            (skips re-resolution from config).
        observer: Optional :class:`ExecutionObserver` to receive
            stage-by-stage progress; one is constructed if a summary
            report path or progress display is configured.
        output_path: Optional directory where the configured output
            writer will materialise the BOA / AOT / TCWV artifacts.
        output_writer: Optional callable matching the
            :class:`OutputWriter` protocol; built from config when
            ``None``.

    Returns:
        :class:`siac.runtime.CorrectionResult` carrying BOA reflectance,
        atmospheric state, cloud mask, and any QA layers the corrector
        emitted.

    Raises:
        :class:`siac.errors.SIACError` (or a subclass) on configuration
        / data / RT / solver failures. Stage timeouts surface as
        :class:`TimeoutError`.
    """
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
