"""Runtime payloads and validators for SIAC execution."""

from siac.runtime.models import (
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
from siac.runtime.validation import (
    _spatial_shape,
    _validate_atmospheric_state,
    _validate_correction_result,
    _validate_observation_bundle,
    _validate_solved_atmosphere,
    _validate_solver_input_bundle,
    _validate_surface_prior,
)

__all__ = [
    "AtmosphericState",
    "BRDFKernelWeights",
    "CorrectionResult",
    "GeometryAngles",
    "ObservationBundle",
    "RTCoefficients",
    "SolvedAtmosphere",
    "SolverInputBundle",
    "SurfacePrior",
    "_spatial_shape",
    "_validate_atmospheric_state",
    "_validate_correction_result",
    "_validate_observation_bundle",
    "_validate_solved_atmosphere",
    "_validate_solver_input_bundle",
    "_validate_surface_prior",
]
