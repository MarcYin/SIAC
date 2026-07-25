"""Prior-fetch execution helpers for the SIAC pipeline."""

from __future__ import annotations

import logging
from dataclasses import replace
from typing import TYPE_CHECKING, Any

import numpy as np

from siac.observability import current_execution_observer

if TYPE_CHECKING:
    from collections.abc import Callable

    from siac.workflows._pipeline_config import PipelineExecutionSettings

logger = logging.getLogger(__name__)


class SubmitAdapter:
    """Adapt a bare submit callable to look like an executor."""

    def __init__(self, submit_fn: Callable[..., Any]) -> None:
        self._submit_fn = submit_fn

    def submit(self, fn: Callable[..., Any], *args: Any, **kwargs: Any) -> Any:
        return self._submit_fn(fn, *args, **kwargs)


def _with_dem_elevation(atmo: Any, config: Any) -> Any:
    """Return ``atmo`` with terrain elevation sampled from the configured DEM.

    The atmospheric provider reports elevation from its own coarse model grid (or
    not at all). M4 later replaces it from the DEM for the solver, but the
    surface prior is built before that — so without this the prior is derived at
    a different elevation than the retrieval that consumes it. Reading the DEM
    here makes the two agree.

    Degrades to the incoming state on any failure: a missing DEM must not abort
    a run, and ``read_elevation_km`` already falls back to sea level itself.
    """

    if atmo is None:
        return atmo
    dem_path = getattr(getattr(config, "paths", None), "dem", None)
    try:
        from siac.geo.dem import read_elevation_km, use_sea_level_elevation

        if use_sea_level_elevation(None if dem_path is None else str(dem_path)):
            return atmo
        elevation = read_elevation_km(atmo.aot, str(dem_path))
        logger.info(
            "M2->M3: terrain elevation from the DEM, median %.3f km "
            "(prior and solver now share one elevation).",
            float(np.nanmedian(np.asarray(elevation.values, dtype=float))),
        )
        return replace(atmo, elevation=elevation)
    except Exception:  # noqa: BLE001 - elevation is an refinement, never fatal
        logger.warning(
            "Could not sample the DEM for the surface prior; keeping the provider's elevation.",
            exc_info=True,
        )
        return atmo


def fetch_priors(
    *,
    submit_fn: Callable[..., Any],
    lut_submit_fn: Callable[..., Any] | None,
    obs: Any,
    config: Any,
    atmo_provider: Callable[..., Any],
    surface_prior_provider: Callable[..., Any],
    rt_model: Any,
    settings: PipelineExecutionSettings,
    backend_label: str,
    call_with_retries_fn: Callable[..., Any],
    maybe_submit_lut_preload_fn: Callable[..., Any],
    stage_timeout_fn: Callable[..., float | None],
    aerosol_resolution_fn: Callable[[Any], float],
    surface_prior_requires_atmo_fn: Callable[[Callable[..., Any]], bool],
    requested_solver_band_names_fn: Callable[[Any], tuple[str, ...] | None],
    timeout_errors: tuple[type[BaseException], ...],
) -> tuple[Any, Any]:
    """Fetch atmospheric and surface priors using the provided submit callable."""
    observer = current_execution_observer()
    observer_id = observer.run_id if observer is not None else None
    m2_timeout = stage_timeout_fn(settings, "M2.atmospheric_prior", "M2", "atmospheric_prior")
    m3_timeout = stage_timeout_fn(settings, "M3.surface_prior", "M3", "surface_prior")
    lut_timeout = stage_timeout_fn(settings, "LUT.preload", "LUT", "lut_preload")
    retries = settings["retries"]
    resolution = aerosol_resolution_fn(config)
    requires_atmo = surface_prior_requires_atmo_fn(surface_prior_provider)
    solver_band_names = requested_solver_band_names_fn(config)

    f_m2 = submit_fn(
        call_with_retries_fn,
        atmo_provider,
        (obs.bounds, obs.crs, obs.metadata["observation_time"], resolution),
        retries=retries,
        stage_name="M2.atmospheric_prior",
        observer_id=observer_id,
    )
    if observer is not None:
        observer.increment_counter("m2_started", stage="M2.atmospheric_prior")
        observer.emit(
            "progress", stage="M2.atmospheric_prior", message="Atmospheric prior submitted."
        )

    f_m3 = None
    if not requires_atmo:
        logger.info("M2+M3: Fetching atmospheric & surface priors...")
        f_m3 = submit_fn(
            call_with_retries_fn,
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
    timeout_stage = "M2/M3"
    active_timeout: float | None = None
    try:
        timeout_stage = "M2.atmospheric_prior"
        active_timeout = m2_timeout
        atmo = f_m2.result(timeout=m2_timeout)
        if observer is not None:
            observer.increment_counter("m2_done", stage="M2.atmospheric_prior")
            observer.emit(
                "progress",
                stage="M2.atmospheric_prior",
                message="Atmospheric prior ready.",
            )
        if f_m3 is None:
            # Give the surface prior the same terrain the solver will use. A prior
            # built at one elevation and solved at another disagree about surface
            # pressure, so molecular scattering is mis-attributed to aerosol; the
            # DEM is otherwise first read in M4, after this stage.
            atmo = _with_dem_elevation(atmo, config)
            logger.info("M3: Deriving surface prior from observation + atmospheric prior...")
            f_m3 = submit_fn(
                call_with_retries_fn,
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
        # Wave 18 (opt 4): pass the solver config through so the LUT
        # preload can build the joint grid-search LUT when the solver's
        # aot/tcwv axes are derivable from config alone — this overlaps
        # the joint LUT 6S build (~98s, the post-wave-17 dominant cost)
        # with the M3 surface-prior reprojection.
        solver_config = getattr(getattr(config, "algorithms", None), "solver", None)
        f_lut = maybe_submit_lut_preload_fn(
            SubmitAdapter(_submit_fn),
            rt_model,
            obs=obs,
            atmo=atmo,
            requested_band_names=solver_band_names,
            retries=retries,
            observer_id=observer_id,
            solver_config=solver_config,
        )
        timeout_stage = "M3.surface_prior"
        active_timeout = m3_timeout
        surface = f_m3.result(timeout=m3_timeout)
        if observer is not None:
            observer.increment_counter("m3_done", stage="M3.surface_prior")
            observer.emit(
                "progress",
                stage="M3.surface_prior",
                message="Surface prior ready.",
            )
        if f_lut is not None:
            try:
                f_lut.result(timeout=lut_timeout)
                if observer is not None:
                    observer.increment_counter("lut_preload_done", stage="LUT.preload")
                    observer.emit(
                        "progress",
                        stage="LUT.preload",
                        message="LUT preload complete.",
                    )
            except timeout_errors:
                if observer is not None and lut_timeout is not None:
                    observer.record_timeout(
                        stage="LUT.preload",
                        timeout_s=float(lut_timeout),
                        backend=backend_label,
                    )
                if lut_timeout is not None:
                    logger.warning(
                        "LUT preload timed out after %.1fs; proceeding with on-demand LUT reads.",
                        float(lut_timeout),
                    )
                else:
                    logger.warning("LUT preload timed out; proceeding with on-demand LUT reads.")
            except (OSError, RuntimeError, ValueError) as exc:
                if observer is not None:
                    observer.record_error(
                        stage="LUT.preload",
                        error_type=type(exc).__name__,
                        error_message=str(exc),
                    )
                logger.warning(
                    "LUT preload failed (%s: %s); proceeding with on-demand LUT reads.",
                    type(exc).__name__,
                    exc,
                )
    except timeout_errors as exc:
        f_m2.cancel()
        if f_m3 is not None:
            f_m3.cancel()
        if f_lut is not None:
            f_lut.cancel()
        if observer is not None and active_timeout is not None:
            observer.record_timeout(
                stage=timeout_stage,
                timeout_s=float(active_timeout),
                backend=backend_label,
            )
        timeout_description = (
            f" after {float(active_timeout):.1f}s" if active_timeout is not None else ""
        )
        raise TimeoutError(
            f"{timeout_stage} timed out{timeout_description} ({backend_label} backend)"
        ) from exc
    except Exception:
        f_m2.cancel()
        if f_m3 is not None:
            f_m3.cancel()
        if f_lut is not None:
            f_lut.cancel()
        raise

    return atmo, surface
