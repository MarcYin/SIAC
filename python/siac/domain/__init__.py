"""Domain-layer contracts and protocol types.

Built-in sensor catalog entries live in ``siac.catalog`` and are re-exported
here as a convenience import surface.
"""

from siac.catalog import (
    LANDSAT8_OLI_CONFIG,
    LANDSAT9_OLI2_CONFIG,
    SENSOR_CONFIGS,
    SENTINEL2A_CONFIG,
    SENTINEL2B_CONFIG,
    SENTINEL2C_CONFIG,
    get_sensor_config,
)
from siac.domain.aoi import AOI
from siac.domain.contracts import (
    AtmosphericState,
    BRDFKernelWeights,
    CorrectionResult,
    GeometryAngles,
    ObservationBundle,
    RTCoefficients,
    SolvedAtmosphere,
    SolverInputBundle,
    SurfacePrior,
)
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

__all__ = [
    "AOI",
    "AtmosphericState",
    "AerosolSolver",
    "AtmosphericPriorProvider",
    "BRDFProductProvider",
    "BRDFKernelWeights",
    "CorrectionResult",
    "GeometryAngles",
    "ObservationBundle",
    "OutputWriter",
    "RTCoefficients",
    "RTModelBackend",
    "SatellitePreprocessor",
    "SolvedAtmosphere",
    "SolverInputBundle",
    "SurfacePrior",
    "SurfacePriorDeriver",
    "SensorBand",
    "SensorConfig",
    "SENTINEL2A_CONFIG",
    "SENTINEL2B_CONFIG",
    "SENTINEL2C_CONFIG",
    "LANDSAT8_OLI_CONFIG",
    "LANDSAT9_OLI2_CONFIG",
    "SENSOR_CONFIGS",
    "get_sensor_config",
]
