"""Built-in sensor catalog definitions."""

from siac.catalog.sensors.landsat import LANDSAT8_OLI_CONFIG as LANDSAT8_OLI_CONFIG
from siac.catalog.sensors.landsat import LANDSAT9_OLI2_CONFIG as LANDSAT9_OLI2_CONFIG
from siac.catalog.sensors.registry import SENSOR_CONFIGS as SENSOR_CONFIGS
from siac.catalog.sensors.registry import get_sensor_config as get_sensor_config
from siac.catalog.sensors.registry import register as register
from siac.catalog.sensors.sentinel2 import (
    SENTINEL2A_CONFIG as SENTINEL2A_CONFIG,
)
from siac.catalog.sensors.sentinel2 import (
    SENTINEL2B_CONFIG as SENTINEL2B_CONFIG,
)
from siac.catalog.sensors.sentinel2 import (
    SENTINEL2C_CONFIG as SENTINEL2C_CONFIG,
)
