"""Application-layer planning and assembly."""

from __future__ import annotations

from importlib import import_module
from typing import Any

_EXPORTS: dict[str, tuple[str, str]] = {
    "build_preprocessor_runtime": ("siac.app.assembly", "build_preprocessor_runtime"),
    "resolve_atmo_provider": ("siac.app.assembly", "resolve_atmo_provider"),
    "resolve_brdf_provider": ("siac.app.assembly", "resolve_brdf_provider"),
    "resolve_corrector": ("siac.app.assembly", "resolve_corrector"),
    "resolve_grid_assembler": ("siac.app.assembly", "resolve_grid_assembler"),
    "resolve_output_writer": ("siac.app.assembly", "resolve_output_writer"),
    "resolve_preprocessor": ("siac.app.assembly", "resolve_preprocessor"),
    "resolve_rt_model_for_pipeline": ("siac.app.assembly", "resolve_rt_model_for_pipeline"),
    "resolve_s2_backend": ("siac.app.assembly", "resolve_s2_backend"),
    "resolve_solver": ("siac.app.assembly", "resolve_solver"),
    "resolve_surface_prior_provider": ("siac.app.assembly", "resolve_surface_prior_provider"),
    "ExecutionPlan": ("siac.app.planning", "ExecutionPlan"),
    "build_execution_plan": ("siac.app.planning", "build_execution_plan"),
    "coerce_aoi_spec": ("siac.app.planning", "coerce_aoi_spec"),
    "resolve_run_config": ("siac.app.planning", "resolve_run_config"),
    "AOISpec": ("siac.app.requests", "AOISpec"),
    "DateSpec": ("siac.app.requests", "DateSpec"),
    "PathLike": ("siac.app.requests", "PathLike"),
    "SceneProcessRequest": ("siac.app.requests", "SceneProcessRequest"),
    "Sentinel2ProcessRequest": ("siac.app.requests", "Sentinel2ProcessRequest"),
    "Sentinel2QuerySpec": ("siac.app.requests", "Sentinel2QuerySpec"),
    "Sentinel2ResolveRequest": ("siac.app.requests", "Sentinel2ResolveRequest"),
    "Sentinel2SearchRequest": ("siac.app.requests", "Sentinel2SearchRequest"),
    "apply_s2_query_defaults": ("siac.app.sentinel2", "apply_s2_query_defaults"),
    "coerce_date": ("siac.app.sentinel2", "coerce_date"),
    "coerce_s2_query": ("siac.app.sentinel2", "coerce_s2_query"),
    "resolve_s2_input": ("siac.app.sentinel2", "resolve_s2_input"),
    "search_sentinel2": ("siac.app.sentinel2", "search_sentinel2"),
}
_SUBMODULES = frozenset({"assembly", "planning", "requests", "sentinel2"})

__all__ = list(_EXPORTS)


def __getattr__(name: str) -> Any:
    if name in _SUBMODULES:
        module = import_module(f"{__name__}.{name}")
        globals()[name] = module
        return module

    try:
        module_name, attr_name = _EXPORTS[name]
    except KeyError as exc:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}") from exc

    value = getattr(import_module(module_name), attr_name)
    globals()[name] = value
    return value


def __dir__() -> list[str]:
    return sorted(list(globals()) + __all__ + list(_SUBMODULES))
