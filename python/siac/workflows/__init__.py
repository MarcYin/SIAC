"""Workflow-layer entrypoints."""

from __future__ import annotations

from importlib import import_module
from typing import Any

_EXPORTS: dict[str, tuple[str, str]] = {
    "run_pipeline": ("siac.workflows.pipeline", "run_pipeline"),
    "execute_plan": ("siac.workflows.scene", "execute_plan"),
    "process_scene": ("siac.workflows.scene", "process_scene"),
    "process_s2": ("siac.workflows.sentinel2", "process_s2"),
}
_SUBMODULES = frozenset({"pipeline", "scene", "sentinel2"})

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
