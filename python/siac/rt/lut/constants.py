"""Constants for LUT backend modules."""

from __future__ import annotations

DEFAULT_LUT_URL = (
    "https://gws-access.jasmin.ac.uk/public/nceo_isp/"
    "libradtran_continental_average_lut_1nm.zarr.zip"
)

LUT_COORD_DIMS = ("sza", "vza", "raa", "aot", "tcwv", "wavelength")
