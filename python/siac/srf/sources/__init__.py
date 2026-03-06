"""Sensor-family SRF source providers."""

from siac.srf.sources.landsat import load_landsat_sensor_config
from siac.srf.sources.planet import load_planet_sensor_config
from siac.srf.sources.sentinel2 import load_sentinel2_sensor_config

__all__ = [
    "load_landsat_sensor_config",
    "load_planet_sensor_config",
    "load_sentinel2_sensor_config",
]
