"""Tests for lazy package-level re-exports."""

from __future__ import annotations

import importlib
import sys

import pytest


def _purge_modules(*module_names: str) -> None:
    for module_name in module_names:
        sys.modules.pop(module_name, None)


@pytest.mark.parametrize(
    ("package_name", "blocked_modules", "export_name", "target_module"),
    [
        (
            "siac.adapters",
            (
                "siac.adapters.atmo",
                "siac.adapters.auth",
                "siac.adapters.brdf",
                "siac.adapters.earthdata",
                "siac.adapters.earthdata_common",
                "siac.adapters.output",
                "siac.adapters.rsrf",
                "siac.adapters.rt",
                "siac.adapters.s2_backend",
                "siac.adapters.satellite",
            ),
            "CredentialManager",
            "siac.adapters.auth",
        ),
        (
            "siac.app",
            (
                "siac.app.assembly",
                "siac.app.planning",
                "siac.app.requests",
                "siac.app.sentinel2",
            ),
            "resolve_run_config",
            "siac.app.planning",
        ),
        (
            "siac.workflows",
            (
                "siac.workflows.pipeline",
                "siac.workflows.scene",
                "siac.workflows.sentinel2",
            ),
            "run_pipeline",
            "siac.workflows.pipeline",
        ),
    ],
)
def test_package_re_exports_are_lazy(
    package_name: str,
    blocked_modules: tuple[str, ...],
    export_name: str,
    target_module: str,
) -> None:
    original_modules = {
        module_name: sys.modules.get(module_name)
        for module_name in (package_name, *blocked_modules)
    }

    try:
        _purge_modules(package_name, *blocked_modules)

        package = importlib.import_module(package_name)

        for blocked_module in blocked_modules:
            assert blocked_module not in sys.modules

        exported = getattr(package, export_name)

        assert callable(exported) or exported is not None
        assert target_module in sys.modules
    finally:
        _purge_modules(package_name, *blocked_modules)
        for module_name, module in original_modules.items():
            if module is not None:
                sys.modules[module_name] = module


@pytest.mark.parametrize(
    "package_name",
    ["siac.adapters", "siac.app", "siac.workflows"],
)
def test_lazy_package_re_exports_unknown_attr(package_name: str) -> None:
    package = importlib.import_module(package_name)

    with pytest.raises(AttributeError, match="has no attribute"):
        _ = package.NOT_A_REAL_EXPORT


@pytest.mark.parametrize(
    ("package_name", "submodule_name"),
    [
        ("siac.adapters", "auth"),
        ("siac.app", "planning"),
        ("siac.workflows", "scene"),
    ],
)
def test_lazy_package_dir_lists_lazy_submodules(package_name: str, submodule_name: str) -> None:
    package = importlib.import_module(package_name)

    assert submodule_name in dir(package)
