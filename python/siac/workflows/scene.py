"""Scene-processing workflows."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Any

from siac.app.planning import build_execution_plan
from siac.observability import derive_execution_report_path
from siac.workflows.pipeline import (
    _OUTPUTS_WRITTEN_METADATA_KEY,
    _resolve_execution_settings,
    run_pipeline,
)

if TYPE_CHECKING:
    from siac.app.planning import ExecutionPlan
    from siac.app.requests import SceneProcessRequest
    from siac.domain.protocols import OutputWriter
    from siac.runtime import CorrectionResult


def write_output(
    result: CorrectionResult,
    output_path: Path | str,
    *,
    output_writer: OutputWriter,
) -> None:
    """Persist SIAC outputs through the configured output adapter."""
    output_writer.write(result, Path(output_path))


def execute_plan(
    plan: ExecutionPlan,
    *,
    execution: Any | None = None,
) -> CorrectionResult:
    """Execute a fully resolved :class:`ExecutionPlan`."""
    return run_pipeline(
        plan.input_path,
        plan.runtime_aoi,
        plan.config,
        preprocessor=plan.preprocessor,
        atmo_provider=plan.atmo_provider,
        surface_prior_provider=plan.surface_prior_provider,
        grid_assembler=plan.grid_assembler,
        solver=plan.solver,
        corrector=plan.corrector,
        rt_model=plan.rt_model,
        execution=execution,
        output_path=getattr(plan, "output_path", None),
        output_writer=getattr(plan, "output_writer", None),
    )


def process_scene(
    request: SceneProcessRequest,
    *,
    preprocessor: Any | None = None,
    atmo_provider: Any | None = None,
    surface_prior_provider: Any | None = None,
    grid_assembler: Any | None = None,
    solver: Any | None = None,
    corrector: Any | None = None,
    rt_model: Any | None = None,
    aoi_resolver: Any | None = None,
    build_preprocessor_runtime_fn: Any | None = None,
    resolve_atmo_provider_fn: Any | None = None,
    resolve_surface_prior_provider_fn: Any | None = None,
    resolve_grid_assembler_fn: Any | None = None,
    resolve_solver_fn: Any | None = None,
    resolve_corrector_fn: Any | None = None,
    resolve_rt_model_fn: Any | None = None,
) -> CorrectionResult:
    """Build a runtime plan and execute it for one scene."""
    plan = build_execution_plan(
        request,
        preprocessor=preprocessor,
        atmo_provider=atmo_provider,
        surface_prior_provider=surface_prior_provider,
        grid_assembler=grid_assembler,
        solver=solver,
        corrector=corrector,
        rt_model=rt_model,
        aoi_resolver=aoi_resolver,
        build_preprocessor_runtime_fn=build_preprocessor_runtime_fn,
        resolve_atmo_provider_fn=resolve_atmo_provider_fn,
        resolve_surface_prior_provider_fn=resolve_surface_prior_provider_fn,
        resolve_grid_assembler_fn=resolve_grid_assembler_fn,
        resolve_solver_fn=resolve_solver_fn,
        resolve_corrector_fn=resolve_corrector_fn,
        resolve_rt_model_fn=resolve_rt_model_fn,
    )
    execution_override: dict[str, Any] | None = None
    settings_config = getattr(plan, "config", request.config)
    settings = _resolve_execution_settings(settings_config, execution=None, max_workers=None)
    if settings["show_progress"] and settings["performance_report_path"] is None:
        default_report_path = derive_execution_report_path(plan.output_path)
        if default_report_path is not None:
            execution_override = {"performance_report_path": default_report_path}

    result = execute_plan(plan, execution=execution_override)
    if (
        plan.output_path is not None
        and plan.output_writer is not None
        and not getattr(result, "metadata", {}).get(_OUTPUTS_WRITTEN_METADATA_KEY, False)
    ):
        write_output(result, plan.output_path, output_writer=plan.output_writer)
    return result
