"""
Satellite preprocessing module.

Provides sensor-specific data readers for TOA reflectance, viewing geometry,
cloud masks, and metadata extraction.

Example:
    >>> from siac.adapters.satellite import get_preprocessor, detect_sensor
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

from siac.adapters.satellite.base import (
    BaseSatellitePreprocessor as BaseSatellitePreprocessor,
)
from siac.adapters.satellite.base import (
    apply_scale_offset as apply_scale_offset,
)
from siac.adapters.satellite.base import (
    create_valid_mask as create_valid_mask,
)
from siac.adapters.satellite.base import (
    degrees_to_radians as degrees_to_radians,
)
from siac.adapters.satellite.base import (
    detect_sensor as detect_sensor,
)
from siac.adapters.satellite.base import (
    get_preprocessor as get_preprocessor,
)
from siac.adapters.satellite.base import (
    list_available_sensors as list_available_sensors,
)
from siac.adapters.satellite.base import (
    radians_to_degrees as radians_to_degrees,
)
from siac.adapters.satellite.base import (
    register_preprocessor as register_preprocessor,
)

# Import sensor-specific preprocessors to trigger registration
from siac.adapters.satellite.sentinel2 import Sentinel2Preprocessor as Sentinel2Preprocessor
