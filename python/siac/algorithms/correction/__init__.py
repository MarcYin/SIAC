"""
Atmospheric correction module.

Applies retrieved atmospheric parameters to convert TOA to BOA reflectance.
"""

from siac.algorithms.correction.atmospheric import AtmosphericCorrector
from siac.runtime import CorrectionResult

__all__ = [
    "AtmosphericCorrector",
    "CorrectionResult",
]
