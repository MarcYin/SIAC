"""Application-layer planning and assembly."""

from __future__ import annotations

from importlib import import_module
from typing import Any

_EXPORTS: dict[str, tuple[str, str]] = {
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
    "resolve_s2_backend": ("siac.app.s2_backend", "resolve_s2_backend"),
    "search_sentinel2": ("siac.app.sentinel2", "search_sentinel2"),
}
_SUBMODULES = frozenset({"planning", "requests", "s2_backend", "sentinel2"})


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
    return sorted(list(globals()) + list(_EXPORTS) + list(_SUBMODULES))
