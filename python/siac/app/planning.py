"""Execution planning for SIAC runs."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Any

from siac.adapters.auth import CredentialManager
from siac.config import RunRequest
from siac.domain.aoi import AOI

if TYPE_CHECKING:
    from pathlib import Path

    from siac.workflows.pipeline import (
        AtmoPriorFn,
        CorrectorFn,
        GridAssemblerFn,
        PreprocessorFn,
        SolverFn,
        SurfacePriorFn,
    )


def resolve_run_config(
    config,
    *,
    input_path: Path | None = None,
    output_path: Path | str | None = None,
    sensor: str | None = None,
    aoi: AOI | Path | str | tuple[float, float, float, float] | list[float] | None = None,
    s2_query: str | Path | None = None,
):
    request = RunRequest(
        input_path=input_path,
        output_path=output_path,
        sensor=sensor,
        aoi=aoi,
        s2_query=s2_query,
    )
    return config.resolve(request)


def coerce_aoi_spec(
    aoi: AOI | Path | str | tuple[float, float, float, float] | list[float] | None,
) -> AOI | Any | None:
    if aoi is None or isinstance(aoi, AOI):
        return aoi
    if isinstance(aoi, (list, tuple)) and len(aoi) == 4:
        return AOI.from_bounds(tuple(aoi))
    return AOI.from_geojson(aoi)


@dataclass(frozen=True)
class ExecutionPlan:
    """Fully assembled runtime plan for one SIAC execution."""

    input_path: Path
    output_path: Path | str | None
    config: Any
    runtime_aoi: AOI | Any | None
    auth: CredentialManager
    preprocessor: PreprocessorFn
    atmo_provider: AtmoPriorFn
    surface_prior_provider: SurfacePriorFn
    grid_assembler: GridAssemblerFn
    solver: SolverFn
    corrector: CorrectorFn
    rt_model: Any


def build_execution_plan(
    config,
    input_path: Path,
    *,
    output_path: Path | str | None = None,
    aoi: AOI | Path | str | tuple[float, float, float, float] | list[float] | None = None,
    auth: CredentialManager | None = None,
    preprocessor: PreprocessorFn | None = None,
    atmo_provider: AtmoPriorFn | None = None,
    surface_prior_provider: SurfacePriorFn | None = None,
    grid_assembler: GridAssemblerFn | None = None,
    solver: SolverFn | None = None,
    corrector: CorrectorFn | None = None,
    rt_model: Any | None = None,
    aoi_resolver=None,
    build_preprocessor_runtime_fn=None,
    resolve_atmo_provider_fn=None,
    resolve_surface_prior_provider_fn=None,
    resolve_grid_assembler_fn=None,
    resolve_solver_fn=None,
    resolve_corrector_fn=None,
    resolve_rt_model_fn=None,
) -> ExecutionPlan:
    if build_preprocessor_runtime_fn is None:
        from siac.app.assembly import build_preprocessor_runtime as build_preprocessor_runtime_fn
    if resolve_atmo_provider_fn is None:
        from siac.app.assembly import resolve_atmo_provider as resolve_atmo_provider_fn
    if resolve_surface_prior_provider_fn is None:
        from siac.app.assembly import (
            resolve_surface_prior_provider as resolve_surface_prior_provider_fn,
        )
    if resolve_grid_assembler_fn is None:
        from siac.app.assembly import resolve_grid_assembler as resolve_grid_assembler_fn
    if resolve_solver_fn is None:
        from siac.app.assembly import resolve_solver as resolve_solver_fn
    if resolve_corrector_fn is None:
        from siac.app.assembly import resolve_corrector as resolve_corrector_fn
    if resolve_rt_model_fn is None:
        from siac.app.assembly import resolve_rt_model_for_pipeline as resolve_rt_model_fn

    resolved_config = resolve_run_config(
        config,
        input_path=input_path,
        output_path=output_path,
        sensor=config.sensor,
        aoi=aoi if aoi is not None else config.aoi,
    )
    runtime_aoi = coerce_aoi_spec(resolved_config.aoi)
    auth_obj = auth or CredentialManager.from_config(resolved_config)

    preprocessor_runtime = build_preprocessor_runtime_fn(
        resolved_config,
        input_path=input_path,
        sensor=resolved_config.sensor,
        default_aoi_resolver=lambda toa: aoi_resolver(toa) if callable(aoi_resolver) else AOI.from_raster(toa[list(toa.data_vars)[0]]),
    )

    return ExecutionPlan(
        input_path=input_path,
        output_path=output_path,
        config=resolved_config,
        runtime_aoi=runtime_aoi,
        auth=auth_obj,
        preprocessor=preprocessor or preprocessor_runtime.preprocessor,
        atmo_provider=atmo_provider or resolve_atmo_provider_fn(resolved_config, auth=auth_obj),
        surface_prior_provider=surface_prior_provider or resolve_surface_prior_provider_fn(resolved_config, auth=auth_obj),
        grid_assembler=grid_assembler or resolve_grid_assembler_fn(),
        solver=solver or resolve_solver_fn(resolved_config),
        corrector=corrector or resolve_corrector_fn(resolved_config),
        rt_model=rt_model or resolve_rt_model_fn(
            resolved_config,
            auth=auth_obj,
            sensor_config=preprocessor_runtime.sensor_config,
        ),
    )


__all__ = [
    "ExecutionPlan",
    "build_execution_plan",
    "coerce_aoi_spec",
    "resolve_run_config",
]
