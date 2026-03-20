"""Pure domain-layer types and protocols."""

from siac.domain.aoi import AOI
from siac.domain.protocols import (
    AerosolSolver,
    AtmosphericPriorProvider,
    BRDFProductProvider,
    OutputWriter,
    RTModelBackend,
    SatellitePreprocessor,
    SurfacePriorDeriver,
)
from siac.domain.sensors import SensorBand, SensorConfig
from siac.domain.spectral import SpectralResponseFunction

__all__ = [
    "AOI",
    "AerosolSolver",
    "AtmosphericPriorProvider",
    "BRDFProductProvider",
    "OutputWriter",
    "RTModelBackend",
    "SatellitePreprocessor",
    "SurfacePriorDeriver",
    "SensorBand",
    "SensorConfig",
    "SpectralResponseFunction",
]
