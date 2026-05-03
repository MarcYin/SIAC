"""
Radiative transfer model backends.

Provides pluggable RT backends:
- emulator: Pre-trained neural network emulators (fast, S2/L8 only)
- lut: Look-up table interpolation (medium speed, any sensor)
- direct: Direct native RT backends (slow, any sensor)
"""

from siac.algorithms.rt.direct import SixSBackend as SixSBackend
from siac.algorithms.rt.emulator import EmulatorRegistry as EmulatorRegistry
from siac.algorithms.rt.emulator import TwoLayerNNEmulator as TwoLayerNNEmulator
from siac.algorithms.rt.lut import ZarrLUTBackend as ZarrLUTBackend
