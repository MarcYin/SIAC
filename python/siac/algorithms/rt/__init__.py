"""
Radiative transfer model backends.

Provides pluggable RT backends:
- emulator: Pre-trained neural network emulators (fast, S2/L8 only)
- lut: Look-up table interpolation (medium speed, any sensor)
- direct: Direct native/Python RT backends (slow, any sensor)
"""

from siac.algorithms.rt.direct import SixSBackend
from siac.algorithms.rt.emulator import EmulatorRegistry, TwoLayerNNEmulator
from siac.algorithms.rt.lut import ZarrLUTBackend, create_lut_from_py6s

__all__ = [
    # Direct backends
    "SixSBackend",
    # Emulator backend
    "TwoLayerNNEmulator",
    "EmulatorRegistry",
    # LUT backend
    "ZarrLUTBackend",
    "create_lut_from_py6s",
]
