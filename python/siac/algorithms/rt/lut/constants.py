"""Constants for LUT backend modules."""

from __future__ import annotations

DEFAULT_LUT_URL = (
    "https://gws-access.jasmin.ac.uk/public/nceo_isp/"
    "libradtran_continental_average_lut_1nm.zarr.zip"
)

LUT_COORD_DIMS = ("sza", "vza", "raa", "aot", "tcwv", "wavelength")

#: SIAC state units -> spectral-LUT schema units. The SIAC atmospheric state
#: carries tcwv in cm of precipitable water and tco3 in atm-cm, while the
#: spectral LUT schema (the remote ``*_lut_1nm`` zarr and the in-memory
#: libRadtran scene LUT that mirrors it) stores tcwv in mm and ozone in Dobson
#: units. Every state->LUT crossing must use these shared factors — scattered
#: bare literals are exactly how the cm-vs-mm / atm-cm-vs-DU lookup bugs
#: shipped (fixed in fd1f65b).
TCWV_CM_TO_LUT_MM = 10.0
TCO3_ATMCM_TO_DU = 1000.0
