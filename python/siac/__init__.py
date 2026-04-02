"""
SIAC - Sensor-Invariant Atmospheric Correction

A modular framework for atmospheric correction of satellite imagery.
Supports Sentinel-2, Landsat 8, and extensible to other sensors.

Example usage:
    >>> from siac import SIAC
    >>> from siac.config import SIACConfig
    >>>
    >>> config = SIACConfig.from_file("siac_config.toml")
    >>> siac = SIAC(config)
    >>> result = siac.process("/path/to/S2_SAFE/")
"""

from importlib import import_module
from importlib.metadata import PackageNotFoundError, version
from typing import Any

try:
    __version__ = version("siac")
except PackageNotFoundError:
    __version__ = "2.0.0-dev"


# Lazy imports for main API: maps name -> (module_path, attribute_name)
_LAZY_IMPORTS: dict[str, tuple[str, str]] = {
    "SIAC": ("siac.api", "SIAC"),
    "SIACConfig": ("siac.config", "SIACConfig"),
    "PreparedMonthlyCompositeBuildResult": ("siac.api", "PreparedMonthlyCompositeBuildResult"),
    "process_sentinel2": ("siac.api", "process_sentinel2"),
    "prepare_monthly_composites": ("siac.api", "prepare_monthly_composites"),
    "process_landsat8": ("siac.api", "process_landsat8"),
    "resolve_s2_input": ("siac.api", "resolve_s2_input"),
    "siac_process_s2": ("siac.api", "siac_process_s2"),
    "search_sentinel2": ("siac.api", "search_sentinel2"),
}


def __getattr__(name: str) -> Any:
    if name in _LAZY_IMPORTS:
        module_path, attr_name = _LAZY_IMPORTS[name]
        return getattr(import_module(module_path), attr_name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


__all__ = [
    "__version__",
    "PreparedMonthlyCompositeBuildResult",
    "SIAC",
    "SIACConfig",
    "prepare_monthly_composites",
    "process_sentinel2",
    "process_landsat8",
    "resolve_s2_input",
    "siac_process_s2",
    "search_sentinel2",
]
