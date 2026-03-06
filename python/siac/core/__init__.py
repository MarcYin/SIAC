"""
Core module containing data types, protocols, and configuration.
"""

from siac.core.aoi import AOI
from siac.core.protocols import (
    AerosolSolver,
    AtmosphericPriorProvider,
    BRDFProductProvider,
    RTModelBackend,
    SatellitePreprocessor,
    SurfacePriorDeriver,
)
from siac.core.srf import SpectralResponseFunction
from siac.core.types import (
    AtmosphericState,
    BRDFKernelWeights,
    CorrectionResult,
    GeometryAngles,
    ObservationBundle,
    RTCoefficients,
    SensorBand,
    SensorConfig,
    SolvedAtmosphere,
    SolverInputBundle,
    SurfacePrior,
)


# Lazy import for config (requires pydantic_settings)
def __getattr__(name):
    if name == "SIACConfig":
        from siac.core.config import SIACConfig
        return SIACConfig
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")

__all__ = [
    # Data types
    "GeometryAngles",
    "AtmosphericState",
    "RTCoefficients",
    "BRDFKernelWeights",
    "SurfacePrior",
    "SensorBand",
    "SensorConfig",
    # Pipeline contracts
    "ObservationBundle",
    "SolverInputBundle",
    "SolvedAtmosphere",
    "CorrectionResult",
    "AOI",
    "SpectralResponseFunction",
    # Protocols
    "SatellitePreprocessor",
    "AtmosphericPriorProvider",
    "BRDFProductProvider",
    "SurfacePriorDeriver",
    "RTModelBackend",
    "AerosolSolver",
    # Configuration
    "SIACConfig",
]
