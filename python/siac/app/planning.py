"""Execution planning for SIAC runs."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Any, cast

from siac.adapters.auth import CredentialManager
from siac.config import RunRequest
from siac.domain.aoi import AOI

if TYPE_CHECKING:
    from collections.abc import Callable

    import xarray as xr

    from siac.app.requests import AOISpec, PathLike, SceneProcessRequest
    from siac.domain.protocols import OutputWriter
    from siac.workflows.pipeline import (
        AtmoPriorFn,
        CorrectorFn,
        GridAssemblerFn,
        PreprocessorFn,
        SolverFn,
        SurfacePriorFn,
    )


def resolve_run_config(
    config: Any,
    *,
    input_path: Path | None = None,
    output_path: PathLike | None = None,
    sensor: Any | None = None,
    aoi: AOISpec = None,
    s2_query: str | Path | None = None,
) -> Any:
    request = RunRequest(
        input_path=input_path,
        output_path=Path(output_path) if isinstance(output_path, str) else output_path,
        sensor=cast("Any", sensor),
        aoi=cast("Any", aoi),
        s2_query=s2_query,
    )
    return config.resolve(request)


def coerce_aoi_spec(
    aoi: AOISpec,
) -> AOI | Any | None:
    if aoi is None or isinstance(aoi, AOI):
        return aoi
    if isinstance(aoi, (list, tuple)) and len(aoi) == 4:
        return AOI.from_bounds(cast("tuple[float, float, float, float]", tuple(aoi)))
    if isinstance(aoi, (dict, Path, str)):
        return AOI.from_geojson(aoi)
    raise ValueError(f"Could not parse AOI from {type(aoi).__name__}")


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
    output_writer: OutputWriter | None


def build_execution_plan(
    request: SceneProcessRequest,
    *,
    preprocessor: PreprocessorFn | None = None,
    atmo_provider: AtmoPriorFn | None = None,
    surface_prior_provider: SurfacePriorFn | None = None,
    grid_assembler: GridAssemblerFn | None = None,
    solver: SolverFn | None = None,
    corrector: CorrectorFn | None = None,
    rt_model: Any | None = None,
    aoi_resolver: Callable[[xr.Dataset], AOI] | None = None,
    build_preprocessor_runtime_fn: Any | None = None,
    resolve_atmo_provider_fn: Any | None = None,
    resolve_surface_prior_provider_fn: Any | None = None,
    resolve_grid_assembler_fn: Any | None = None,
    resolve_solver_fn: Any | None = None,
    resolve_corrector_fn: Any | None = None,
    resolve_rt_model_fn: Any | None = None,
    resolve_output_writer_fn: Any | None = None,
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
    if resolve_output_writer_fn is None:
        from siac.app.assembly import resolve_output_writer as resolve_output_writer_fn

    build_preprocessor_runtime = cast("Any", build_preprocessor_runtime_fn)
    resolve_atmo_provider = cast("Any", resolve_atmo_provider_fn)
    resolve_surface_prior_provider = cast("Any", resolve_surface_prior_provider_fn)
    resolve_grid_assembler = cast("Any", resolve_grid_assembler_fn)
    resolve_solver = cast("Any", resolve_solver_fn)
    resolve_corrector = cast("Any", resolve_corrector_fn)
    resolve_rt_model = cast("Any", resolve_rt_model_fn)
    resolve_output_writer = cast("Any", resolve_output_writer_fn)

    resolved_config = resolve_run_config(
        request.config,
        input_path=Path(request.input_path),
        output_path=request.output_path,
        sensor=request.config.sensor,
        aoi=request.aoi if request.aoi is not None else request.config.aoi,
    )
    runtime_aoi = coerce_aoi_spec(resolved_config.aoi)
    auth_obj = request.auth or CredentialManager.from_config(resolved_config)

    def _default_aoi_resolver(toa: xr.Dataset) -> AOI:
        if callable(aoi_resolver):
            return aoi_resolver(toa)
        return AOI.from_raster(toa[list(toa.data_vars)[0]])

    preprocessor_runtime = build_preprocessor_runtime(
        resolved_config,
        input_path=Path(request.input_path),
        sensor=resolved_config.sensor,
        default_aoi_resolver=_default_aoi_resolver,
    )

    return ExecutionPlan(
        input_path=Path(request.input_path),
        output_path=resolved_config.run.output_path,
        config=resolved_config,
        runtime_aoi=runtime_aoi,
        auth=auth_obj,
        preprocessor=preprocessor or preprocessor_runtime.preprocessor,
        atmo_provider=atmo_provider or resolve_atmo_provider(resolved_config, auth=auth_obj),
        surface_prior_provider=surface_prior_provider
        or resolve_surface_prior_provider(resolved_config, auth=auth_obj),
        grid_assembler=grid_assembler or resolve_grid_assembler(),
        solver=solver or resolve_solver(resolved_config),
        corrector=corrector or resolve_corrector(resolved_config),
        rt_model=rt_model
        or resolve_rt_model(
            resolved_config,
            auth=auth_obj,
            sensor_config=preprocessor_runtime.sensor_config,
        ),
        output_writer=(
            resolve_output_writer(resolved_config)
            if resolved_config.run.output_path is not None
            else None
        ),
    )


__all__ = [
    "ExecutionPlan",
    "build_execution_plan",
    "coerce_aoi_spec",
    "resolve_run_config",
]
