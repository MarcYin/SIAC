"""
Look-up table (Zarr) backend for radiative transfer.

Provides RT coefficient computation using pre-computed LUTs stored in
Zarr format. Supports any sensor with known spectral response functions.
"""

from siac.algorithms.rt.lut.backend import ZarrLUTBackend as ZarrLUTBackend
from siac.algorithms.rt.lut.constants import DEFAULT_LUT_URL as DEFAULT_LUT_URL
