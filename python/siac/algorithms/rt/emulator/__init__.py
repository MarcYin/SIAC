"""
Neural network emulator backend.

Provides fast RT coefficient computation using pre-trained neural networks
that approximate 6S radiative transfer model outputs.
"""

from siac.algorithms.rt.emulator.two_nn import (
    EmulatorRegistry as EmulatorRegistry,
)
from siac.algorithms.rt.emulator.two_nn import (
    TwoLayerNNEmulator as TwoLayerNNEmulator,
)
