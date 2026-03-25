"""Pure domain-layer types and protocols."""

from siac.domain.aoi import AOI
from siac.domain.protocols import (
    AerosolSolver,
    AtmosphericPriorProvider,
    BRDFProductProvider,
    MonthlyCompositeProvider,
    OutputWriter,
    RTModelBackend,
    SatellitePreprocessor,
    SurfacePriorDeriver,
)
from siac.domain.sensors import SensorBand, SensorConfig
from siac.domain.spectral import RelativeSpectralResponse, SpectralResponseFunction

__all__ = [
    "AOI",
    "AerosolSolver",
    "AtmosphericPriorProvider",
    "BRDFProductProvider",
    "MonthlyCompositeProvider",
    "OutputWriter",
    "RTModelBackend",
    "SatellitePreprocessor",
    "SurfacePriorDeriver",
    "SensorBand",
    "SensorConfig",
    "RelativeSpectralResponse",
    "SpectralResponseFunction",
]
