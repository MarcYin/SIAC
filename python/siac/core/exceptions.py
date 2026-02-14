"""
Custom exceptions for SIAC.
"""


class SIACError(Exception):
    """Base exception for all SIAC errors."""

    pass


class ConfigurationError(SIACError):
    """Error in SIAC configuration."""

    pass


class DataNotFoundError(SIACError):
    """Required data not found."""

    pass


class SensorNotSupportedError(SIACError):
    """Sensor is not supported by the requested backend."""

    pass


class RTModelError(SIACError):
    """Error in radiative transfer model computation."""

    pass


class SolverConvergenceError(SIACError):
    """Aerosol solver failed to converge."""

    pass


class ValidationError(SIACError):
    """Data validation error."""

    pass


class AuthenticationError(SIACError):
    """Credentials missing or authentication failed for a data provider."""

    pass
