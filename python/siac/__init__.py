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

from importlib.metadata import PackageNotFoundError, version
from typing import Any

try:
    __version__ = version("siac")
except PackageNotFoundError:
    __version__ = "2.0.0-dev"


# Lazy imports for main API
def __getattr__(name: str) -> Any:
    if name == "SIAC":
        from siac.api import SIAC

        return SIAC
    if name == "SIACConfig":
        from siac.config import SIACConfig

        return SIACConfig
    if name == "process_sentinel2":
        from siac.api import process_sentinel2

        return process_sentinel2
    if name == "process_landsat8":
        from siac.api import process_landsat8

        return process_landsat8
    if name == "resolve_s2_input":
        from siac.api import resolve_s2_input

        return resolve_s2_input
    if name == "siac_process_s2":
        from siac.api import siac_process_s2

        return siac_process_s2
    if name == "search_sentinel2":
        from siac.api import search_sentinel2

        return search_sentinel2
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


__all__ = [
    "__version__",
    "SIAC",
    "SIACConfig",
    "process_sentinel2",
    "process_landsat8",
    "resolve_s2_input",
    "siac_process_s2",
    "search_sentinel2",
]
