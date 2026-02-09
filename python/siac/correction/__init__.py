"""
Atmospheric correction module.

Applies retrieved atmospheric parameters to convert TOA to BOA reflectance.
"""

from siac.correction.atmospheric import AtmosphericCorrector
from siac.core.types import CorrectionResult

__all__ = [
    "AtmosphericCorrector",
    "CorrectionResult",
]
