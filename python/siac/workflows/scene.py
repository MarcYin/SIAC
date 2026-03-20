"""Scene-processing workflows."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Any

from siac.app.planning import build_execution_plan
from siac.workflows.pipeline import run_pipeline

if TYPE_CHECKING:
    from siac.app.requests import SceneProcessRequest
    from siac.runtime import CorrectionResult


def save_output(result: CorrectionResult, output_path: Path | str) -> None:
    """Persist SIAC outputs to disk."""
    from siac.storage import write_dataset

    resolved = Path(output_path)
    resolved.mkdir(parents=True, exist_ok=True)
    write_dataset(result.boa, resolved / "boa.nc")


def execute_plan(plan) -> CorrectionResult:
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
    )


def process_scene(
    request: SceneProcessRequest,
    *,
    preprocessor=None,
    atmo_provider=None,
    surface_prior_provider=None,
    grid_assembler=None,
    solver=None,
    corrector=None,
    rt_model: Any | None = None,
    aoi_resolver=None,
    build_preprocessor_runtime_fn=None,
    resolve_atmo_provider_fn=None,
    resolve_surface_prior_provider_fn=None,
    resolve_grid_assembler_fn=None,
    resolve_solver_fn=None,
    resolve_corrector_fn=None,
    resolve_rt_model_fn=None,
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
    result = execute_plan(plan)
    if request.output_path is not None:
        save_output(result, request.output_path)
    return result


__all__ = ["execute_plan", "process_scene", "save_output"]
