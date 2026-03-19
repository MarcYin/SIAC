"""Static catalog data for built-in SIAC assets."""

from siac.catalog.sensors import (
    LANDSAT8_OLI_CONFIG,
    LANDSAT9_OLI2_CONFIG,
    SENSOR_CONFIGS,
    SENTINEL2A_CONFIG,
    SENTINEL2B_CONFIG,
    SENTINEL2C_CONFIG,
    get_sensor_config,
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
