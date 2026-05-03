"""Thread and Dask execution backends for pipeline orchestration."""

from __future__ import annotations

import logging
from concurrent.futures import ThreadPoolExecutor
from concurrent.futures import TimeoutError as FuturesTimeoutError
from contextlib import nullcontext
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any

from siac.observability import current_execution_observer

if TYPE_CHECKING:
    from collections.abc import Callable
    from pathlib import Path

    from siac.workflows._pipeline_config import PipelineExecutionSettings

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class PipelineExecutionContext:
    preprocessor: Callable[..., Any]
    atmo_provider: Callable[..., Any]
    surface_prior_provider: Callable[..., Any]
    grid_assembler: Callable[..., Any]
    solver: Callable[..., Any]
    corrector: Callable[..., Any]
    rt_model: Any
    settings: PipelineExecutionSettings
    output_path: Path | str | None
    output_writer: Any | None
    prepare_observation: Callable[..., tuple[Any, Any, Any, Any, Any]]
    set_rt_observation_time: Callable[[Any, Any], None]
    fetch_priors: Callable[..., tuple[Any, Any]]
    run_tail: Callable[..., Any]


def run_pipeline_thread(
    input_path: Path,
    aoi: Any | None,
    config: Any,
    *,
    context: PipelineExecutionContext,
) -> Any:
    obs, _bounds, _crs, obs_time, _resolution = context.prepare_observation(
        input_path,
        aoi,
        config,
        preprocessor=context.preprocessor,
    )
    context.set_rt_observation_time(context.rt_model, obs_time)

    with ThreadPoolExecutor(max_workers=context.settings["max_workers"]) as executor:
        atmo, surface = context.fetch_priors(
            submit_fn=executor.submit,
            lut_submit_fn=None,
            obs=obs,
            config=config,
            atmo_provider=context.atmo_provider,
            surface_prior_provider=context.surface_prior_provider,
            rt_model=context.rt_model,
            settings=context.settings,
            backend_label="thread",
        )

    return context.run_tail(
        obs,
        atmo,
        surface,
        config,
        grid_assembler=context.grid_assembler,
        solver=context.solver,
        corrector=context.corrector,
        rt_model=context.rt_model,
        output_path=context.output_path,
        output_writer=context.output_writer,
    )


def run_pipeline_dask(
    input_path: Path,
    aoi: Any | None,
    config: Any,
    *,
    context: PipelineExecutionContext,
) -> Any:
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
        "n_workers": context.settings["max_workers"],
        "threads_per_worker": 1,
        "processes": False,
        "dashboard_address": (
            context.settings["dashboard_address"] if context.settings["dashboard"] else None
        ),
    }
    observer = current_execution_observer()

    with LocalCluster(**cluster_kwargs) as cluster, Client(cluster) as client:
        if context.settings["show_progress"] and getattr(client, "dashboard_link", None):
            logger.info("Dask dashboard: %s", client.dashboard_link)
            if observer is not None:
                observer.emit(
                    "progress",
                    stage="dask.cluster",
                    message="Dask dashboard available.",
                    dashboard_link=client.dashboard_link,
                )

        report_ctx = nullcontext()
        report_path = context.settings["performance_report_path"]
        if report_path is not None:
            report_path.parent.mkdir(parents=True, exist_ok=True)
            report_ctx = performance_report(filename=str(report_path))

        with report_ctx:
            obs, _bounds, _crs, obs_time, _resolution = context.prepare_observation(
                input_path,
                aoi,
                config,
                preprocessor=context.preprocessor,
            )
            context.set_rt_observation_time(context.rt_model, obs_time)

            preload_executor = ThreadPoolExecutor(max_workers=1)
            try:
                atmo, surface = context.fetch_priors(
                    submit_fn=client.submit,
                    lut_submit_fn=preload_executor.submit,
                    obs=obs,
                    config=config,
                    atmo_provider=context.atmo_provider,
                    surface_prior_provider=context.surface_prior_provider,
                    rt_model=context.rt_model,
                    settings=context.settings,
                    backend_label="dask",
                    timeout_errors=(FuturesTimeoutError, DaskTimeoutError),
                )
            finally:
                preload_executor.shutdown(wait=True, cancel_futures=False)

    return context.run_tail(
        obs,
        atmo,
        surface,
        config,
        grid_assembler=context.grid_assembler,
        solver=context.solver,
        corrector=context.corrector,
        rt_model=context.rt_model,
        output_path=context.output_path,
        output_writer=context.output_writer,
    )
