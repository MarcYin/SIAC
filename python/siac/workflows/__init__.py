"""Workflow-layer entrypoints."""

from __future__ import annotations

from importlib import import_module
from typing import Any

_EXPORTS: dict[str, tuple[str, str]] = {
    "run_pipeline": ("siac.workflows.pipeline", "run_pipeline"),
    "aod_scene_mean": ("siac.workflows.scene", "aod_scene_mean"),
    "aod_quality_score": ("siac.workflows.scene", "aod_quality_score"),
    "aod_rail_value": ("siac.workflows.scene", "aod_rail_value"),
    "execute_plan": ("siac.workflows.scene", "execute_plan"),
    "process_scene_with_aod_rail_fallback": (
        "siac.workflows.scene",
        "process_scene_with_aod_rail_fallback",
    ),
    "process_scene_with_aod_selector": (
        "siac.workflows.scene",
        "process_scene_with_aod_selector",
    ),
    "process_scene_with_aod_ensemble": (
        "siac.workflows.scene",
        "process_scene_with_aod_ensemble",
    ),
    "process_scene": ("siac.workflows.scene", "process_scene"),
    "process_s2": ("siac.workflows.sentinel2", "process_s2"),
}
_SUBMODULES = frozenset({"pipeline", "scene", "sentinel2"})


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
