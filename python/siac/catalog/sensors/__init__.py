"""Built-in sensor catalog definitions."""

from siac.catalog.sensors.landsat import LANDSAT8_OLI_CONFIG, LANDSAT9_OLI2_CONFIG
from siac.catalog.sensors.registry import SENSOR_CONFIGS, get_sensor_config
from siac.catalog.sensors.sentinel2 import (
    SENTINEL2A_CONFIG,
    SENTINEL2B_CONFIG,
    SENTINEL2C_CONFIG,
)

__all__ = [
    "SENTINEL2A_CONFIG",
    "SENTINEL2B_CONFIG",
    "SENTINEL2C_CONFIG",
    "LANDSAT8_OLI_CONFIG",
    "LANDSAT9_OLI2_CONFIG",
    "SENSOR_CONFIGS",
    "get_sensor_config",
]
