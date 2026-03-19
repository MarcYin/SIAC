"""
Shared exceptions for SIAC.
"""


class SIACError(Exception):
    """Base exception for all SIAC errors."""


class ConfigurationError(SIACError):
    """Error in SIAC configuration."""


class DataNotFoundError(SIACError):
    """Required data not found."""


class SensorNotSupportedError(SIACError):
    """Sensor is not supported by the requested backend."""


class RTModelError(SIACError):
    """Error in radiative transfer model computation."""


class SolverConvergenceError(SIACError):
    """Aerosol solver failed to converge."""


class ValidationError(SIACError):
    """Data validation error."""


class AuthenticationError(SIACError):
    """Credentials missing or authentication failed for a data provider."""


class PreprocessingError(SIACError):
    """Error during satellite data preprocessing."""
