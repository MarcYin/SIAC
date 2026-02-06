"""
Core module containing data types, protocols, and configuration.
"""

from siac.core.types import (
    GeometryAngles,
    AtmosphericState,
    RTCoefficients,
    BRDFKernelWeights,
    SurfacePrior,
    SensorBand,
    SensorConfig,
)
from siac.core.protocols import (
    SatellitePreprocessor,
    AtmosphericPriorProvider,
    BRDFProductProvider,
    SurfacePriorDeriver,
    RTModelBackend,
    AerosolSolver,
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
