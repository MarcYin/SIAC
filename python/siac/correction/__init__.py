"""
Atmospheric correction module.

Applies retrieved atmospheric parameters to convert TOA to BOA reflectance.
"""

from siac.core.types import CorrectionResult
from siac.correction.atmospheric import AtmosphericCorrector

__all__ = [
    "AtmosphericCorrector",
    "CorrectionResult",
]
