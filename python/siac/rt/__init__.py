"""
Radiative transfer model backends.

Provides pluggable RT backends:
- emulator: Pre-trained neural network emulators (fast, S2/L8 only)
- lut: Look-up table interpolation (medium speed, any sensor)
- direct: Direct Py6S simulation (slow, any sensor)
"""

from siac.rt.emulator import TwoLayerNNEmulator, EmulatorRegistry
from siac.rt.lut import ZarrLUTBackend, create_lut_from_py6s

__all__ = [
    # Emulator backend
    "TwoLayerNNEmulator",
    "EmulatorRegistry",
    # LUT backend
    "ZarrLUTBackend",
    "create_lut_from_py6s",
]
