"""Metadata-driven spectral characterization providers."""

from siac.srf.metadata.emit import load_emit_sensor_config_from_metadata
from siac.srf.metadata.enmap import load_enmap_sensor_config_from_metadata
from siac.srf.metadata.prisma import load_prisma_sensor_config_from_metadata

__all__ = [
    "load_emit_sensor_config_from_metadata",
    "load_enmap_sensor_config_from_metadata",
    "load_prisma_sensor_config_from_metadata",
]
