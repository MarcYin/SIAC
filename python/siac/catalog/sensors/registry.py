"""Built-in sensor registry."""

from siac.catalog.sensors.landsat import LANDSAT8_OLI_CONFIG, LANDSAT9_OLI2_CONFIG
from siac.catalog.sensors.sentinel2 import (
    SENTINEL2A_CONFIG,
    SENTINEL2B_CONFIG,
    SENTINEL2C_CONFIG,
)
from siac.domain.sensors import SensorConfig

SENSOR_CONFIGS: dict[tuple[str, str], SensorConfig] = {
    ("MSI", "S2A"): SENTINEL2A_CONFIG,
    ("MSI", "S2B"): SENTINEL2B_CONFIG,
    ("MSI", "S2C"): SENTINEL2C_CONFIG,
    ("OLI", "L8"): LANDSAT8_OLI_CONFIG,
    ("OLI", "L9"): LANDSAT9_OLI2_CONFIG,
}


def get_sensor_config(sensor_id: str, satellite_id: str) -> SensorConfig:
    """Get sensor configuration by ID."""
    key = (sensor_id, satellite_id)
    if key not in SENSOR_CONFIGS:
        raise KeyError(
            f"Unknown sensor {sensor_id}/{satellite_id}. Available: {list(SENSOR_CONFIGS.keys())}"
        )
    return SENSOR_CONFIGS[key]
