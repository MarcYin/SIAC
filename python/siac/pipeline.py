"""
Pipeline orchestration for SIAC atmospheric correction.

This module contains the ``run_pipeline()`` function that wires modules
M1-M6 together. Each module is passed in as a plain callable.
"""

from __future__ import annotations

import logging
from collections.abc import Callable
from concurrent.futures import ThreadPoolExecutor
from concurrent.futures import TimeoutError as FuturesTimeoutError
from contextlib import nullcontext
from datetime import datetime
from pathlib import Path
from typing import TYPE_CHECKING, Any

from siac.core.types import (
    AtmosphericState,
    CorrectionResult,
    GeometryAngles,
    ObservationBundle,
    SensorConfig,
    SolvedAtmosphere,
    SolverInputBundle,
    SurfacePrior,
)
from siac.core.validation import (
    _validate_atmospheric_state,
    _validate_observation_bundle,
    _validate_surface_prior,
)

if TYPE_CHECKING:
    from siac.core.aoi import AOI

logger = logging.getLogger(__name__)

# ── Type aliases for module callables ──────────────────────────────────

# M1 preprocessor
PreprocessorFn = Callable[[Path, Any], ObservationBundle]  # (path, aoi|None) -> ObservationBundle

# M2 atmospheric prior provider
AtmoPriorFn = Callable[
    [tuple[float, float, float, float], str, datetime, float],
    AtmosphericState,
]

# M3 surface prior provider
SurfacePriorFn = Callable[
    [tuple[float, float, float, float], str, datetime, SensorConfig, GeometryAngles, float],
    SurfacePrior,
]

# M4 grid assembler
GridAssemblerFn = Callable[
    [ObservationBundle, AtmosphericState, SurfacePrior, Any, float, float],
    SolverInputBundle,
]

# M5 aerosol solver
SolverFn = Callable[[SolverInputBundle, Any], SolvedAtmosphere]

# M6 atmospheric corrector
CorrectorFn = Callable[[ObservationBundle, SolvedAtmosphere, Any], CorrectionResult]


_EXECUTION_KEYS = (
    "backend",
    "max_workers",
    "retries",
    "stage_timeout_s",
    "dashboard",
    "dashboard_address",
    "performance_report_path",
    "show_progress",
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
    return settings


def _aerosol_resolution(config: Any) -> float:
    solver_cfg = getattr(config, "solver", None)
    if solver_cfg is not None:
        return float(getattr(solver_cfg, "aerosol_resolution", 1000.0))
    return 1000.0


def _select_solver_bands_for_preload(sensor_config: SensorConfig) -> list[Any]:
    """Mirror M4 band-selection logic for LUT preloading hints."""
    bands = sensor_config.select_bands_in_range(400.0, 520.0)
    if bands:
        return list(bands)
    return list(sensor_config.bands[:2])


def _maybe_submit_lut_preload(
    executor: Any,
    rt_model: Any,
    *,
    obs: ObservationBundle,
    atmo: AtmosphericState,
    retries: int,
) -> Any | None:
    """Start LUT scene preload task when backend supports it."""
    preload_fn = getattr(rt_model, "preload_scene_subset", None)
    if not callable(preload_fn):
        return None

    bands = _select_solver_bands_for_preload(obs.sensor_config)
    logger.info(
        "Starting LUT preload in parallel (scene subset + %d band grid%s).",
        len(bands),
        "" if len(bands) == 1 else "s",
    )
    return executor.submit(
        _call_with_retries,
        preload_fn,
        (obs.geometry, atmo, bands),
        retries=retries,
        stage_name="LUT preload",
    )


def _call_with_retries(
    fn: Callable[..., Any],
    args: tuple[Any, ...],
    *,
    retries: int,
    stage_name: str,
) -> Any:
    """Call a stage function with bounded retries."""
    for attempt in range(retries + 1):
        try:
            return fn(*args)
        except Exception:
            if attempt >= retries:
                raise
            logger.warning(
                "%s failed on attempt %d/%d; retrying",
                stage_name,
                attempt + 1,
                retries + 1,
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
    obs = preprocessor(input_path, aoi)
    _validate_observation_bundle(obs)

    bounds = obs.bounds
    crs = obs.crs
    obs_time = obs.metadata["observation_time"]
    resolution = _aerosol_resolution(config)
    return obs, bounds, crs, obs_time, resolution


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
) -> CorrectionResult:
    _validate_atmospheric_state(atmo)
    _validate_surface_prior(surface)

    logger.info("M4: Assembling solver grids...")
    solver_inputs = grid_assembler(obs, atmo, surface, rt_model)

    logger.info("M5: Solving for aerosol parameters...")
    solved = solver(solver_inputs, config)

    logger.info("M6: Applying atmospheric correction...")
    result = corrector(obs, solved, rt_model)

    logger.info("Pipeline complete.")
    return result


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
) -> CorrectionResult:
    obs, bounds, crs, obs_time, resolution = _prepare_observation(
        input_path,
        aoi,
        config,
        preprocessor=preprocessor,
    )

    timeout = settings["stage_timeout_s"]
    retries = settings["retries"]
    logger.info("M2+M3: Fetching atmospheric & surface priors...")
    with ThreadPoolExecutor(max_workers=settings["max_workers"]) as executor:
        f_m2 = executor.submit(
            _call_with_retries,
            atmo_provider,
            (bounds, crs, obs_time, resolution),
            retries=retries,
            stage_name="M2",
        )
        f_m3 = executor.submit(
            _call_with_retries,
            surface_prior_provider,
            (bounds, crs, obs_time, obs.sensor_config, obs.geometry, resolution),
            retries=retries,
            stage_name="M3",
        )
        f_lut = None
        try:
            atmo = f_m2.result(timeout=timeout)
            f_lut = _maybe_submit_lut_preload(
                executor,
                rt_model,
                obs=obs,
                atmo=atmo,
                retries=retries,
            )
            surface = f_m3.result(timeout=timeout)
            if f_lut is not None:
                try:
                    f_lut.result(timeout=timeout)
                except FuturesTimeoutError:
                    logger.warning(
                        "LUT preload timed out after %.1fs; proceeding with on-demand LUT reads.",
                        float(timeout),
                    )
                except Exception as exc:
                    logger.warning(
                        "LUT preload failed (%s); proceeding with on-demand LUT reads.",
                        exc,
                    )
        except FuturesTimeoutError as exc:
            f_m2.cancel()
            f_m3.cancel()
            if f_lut is not None:
                f_lut.cancel()
            raise TimeoutError(
                f"M2/M3 timed out after {timeout:.1f}s (thread backend)"
            ) from exc
        except Exception:
            f_m2.cancel()
            f_m3.cancel()
            if f_lut is not None:
                f_lut.cancel()
            raise

    return _run_tail(
        obs,
        atmo,
        surface,
        config,
        grid_assembler=grid_assembler,
        solver=solver,
        corrector=corrector,
        rt_model=rt_model,
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
) -> CorrectionResult:
    try:
        from dask.distributed import (  # type: ignore[import-not-found]
            Client,
            LocalCluster,
            performance_report,
        )
        from dask.distributed import (
            TimeoutError as DaskTimeoutError,
        )
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
    timeout = settings["stage_timeout_s"]
    retries = settings["retries"]

    with LocalCluster(**cluster_kwargs) as cluster, Client(cluster) as client:
        if settings["show_progress"] and getattr(client, "dashboard_link", None):
            logger.info("Dask dashboard: %s", client.dashboard_link)

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

            logger.info("M2+M3: Fetching atmospheric & surface priors...")
            f_m2 = client.submit(
                _call_with_retries,
                atmo_provider,
                (bounds, crs, obs_time, resolution),
                retries=retries,
                stage_name="M2",
            )
            f_m3 = client.submit(
                _call_with_retries,
                surface_prior_provider,
                (bounds, crs, obs_time, obs.sensor_config, obs.geometry, resolution),
                retries=retries,
                stage_name="M3",
            )
            preload_future = None
            preload_executor = ThreadPoolExecutor(max_workers=1)
            try:
                atmo = f_m2.result(timeout=timeout)
                preload_future = _maybe_submit_lut_preload(
                    preload_executor,
                    rt_model,
                    obs=obs,
                    atmo=atmo,
                    retries=retries,
                )
                surface = f_m3.result(timeout=timeout)
                if preload_future is not None:
                    try:
                        preload_future.result(timeout=timeout)
                    except FuturesTimeoutError:
                        logger.warning(
                            "LUT preload timed out after %.1fs; proceeding with on-demand LUT reads.",
                            float(timeout),
                        )
                    except Exception as exc:
                        logger.warning(
                            "LUT preload failed (%s); proceeding with on-demand LUT reads.",
                            exc,
                        )
            except DaskTimeoutError as exc:
                f_m2.cancel()
                f_m3.cancel()
                if preload_future is not None:
                    preload_future.cancel()
                raise TimeoutError(
                    f"M2/M3 timed out after {timeout:.1f}s (dask backend)"
                ) from exc
            except Exception:
                f_m2.cancel()
                f_m3.cancel()
                if preload_future is not None:
                    preload_future.cancel()
                raise
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
) -> CorrectionResult:
    """Orchestrate module execution with a selectable execution backend."""
    settings = _resolve_execution_settings(
        config,
        execution=execution,
        max_workers=max_workers,
    )

    if settings["backend"] == "dask":
        return _run_pipeline_dask(
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
        )

    return _run_pipeline_thread(
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
    )
