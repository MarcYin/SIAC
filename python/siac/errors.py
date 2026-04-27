"""
Shared exceptions for SIAC.
"""

from __future__ import annotations


class SIACError(Exception):
    """Base exception for all SIAC errors.

    Attributes:
        stage: Optional pipeline stage name where the error originated
            (e.g. ``"M1_preprocessor"``, ``"M5_solver"``).
        is_recoverable: Hint to callers whether the pipeline can continue.
    """

    stage: str | None = None
    is_recoverable: bool = False

    def __init__(
        self,
        *args: object,
        stage: str | None = None,
        is_recoverable: bool | None = None,
    ) -> None:
        super().__init__(*args)
        if stage is not None:
            self.stage = stage
        if is_recoverable is not None:
            self.is_recoverable = is_recoverable


class ConfigurationError(SIACError):
    """Error in SIAC configuration."""


class DataNotFoundError(SIACError):
    """Required data not found."""

    is_recoverable = False


class SensorNotSupportedError(SIACError):
    """Sensor is not supported by the requested backend."""


class RTModelError(SIACError):
    """Error in radiative transfer model computation."""


class SolverConvergenceError(SIACError):
    """Aerosol solver failed to converge."""

    is_recoverable = True


class ValidationError(SIACError):
    """Data validation error."""


class AuthenticationError(SIACError):
    """Credentials missing or authentication failed for a data provider."""


class PreprocessingError(SIACError):
    """Error during satellite data preprocessing."""


class ResamplingError(SIACError):
    """Spatial resampling between grids failed."""

    is_recoverable = True


class EmulatorError(RTModelError):
    """Error during neural-network emulator evaluation."""


class LUTInterpolationError(RTModelError):
    """Error during look-up table interpolation."""

    is_recoverable = True
