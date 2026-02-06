"""
SIAC - Sensor-Invariant Atmospheric Correction

A modular framework for atmospheric correction of satellite imagery.
Supports Sentinel-2, Landsat 8, and extensible to other sensors.

Example usage:
    >>> from siac import SIAC
    >>> from siac.core.config import SIACConfig
    >>>
    >>> config = SIACConfig.from_yaml("siac_config.yaml")
    >>> siac = SIAC(config)
    >>> result = siac.process("/path/to/S2_SAFE/")
"""

from importlib.metadata import version, PackageNotFoundError

try:
    __version__ = version("siac")
except PackageNotFoundError:
    __version__ = "2.0.0-dev"

# Lazy imports for main API
def __getattr__(name: str):
    if name == "SIAC":
        from siac.siac import SIAC
        return SIAC
    if name == "SIACConfig":
        from siac.core.config import SIACConfig
        return SIACConfig
    if name == "process_sentinel2":
        from siac.siac import process_sentinel2
        return process_sentinel2
    if name == "process_landsat8":
        from siac.siac import process_landsat8
        return process_landsat8
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")

__all__ = [
    "__version__",
    "SIAC",
    "SIACConfig",
    "process_sentinel2",
    "process_landsat8",
]
