"""
Satellite preprocessing module.

Provides sensor-specific data readers for TOA reflectance, viewing geometry,
cloud masks, and metadata extraction.

Example:
    >>> from siac.satellite import get_preprocessor, detect_sensor
    >>>
    >>> # Auto-detect sensor
    >>> sensor_id = detect_sensor("/path/to/data/")
    >>>
    >>> # Get preprocessor
    >>> preprocessor = get_preprocessor(sensor_id)
    >>>
    >>> # Run preprocessing
    >>> result = preprocessor.preprocess("/path/to/data/")
    >>> print(result["toa"])
"""

from siac.satellite.base import (
    BaseSatellitePreprocessor,
    apply_scale_offset,
    create_valid_mask,
    degrees_to_radians,
    detect_sensor,
    get_preprocessor,
    list_available_sensors,
    radians_to_degrees,
    register_preprocessor,
)

# Import sensor-specific preprocessors to trigger registration
from siac.satellite.sentinel2 import Sentinel2Preprocessor

__all__ = [
    # Base classes
    "BaseSatellitePreprocessor",
    # Registry functions
    "get_preprocessor",
    "detect_sensor",
    "register_preprocessor",
    "list_available_sensors",
    # Utilities
    "degrees_to_radians",
    "radians_to_degrees",
    "apply_scale_offset",
    "create_valid_mask",
    # Preprocessors
    "Sentinel2Preprocessor",
]
