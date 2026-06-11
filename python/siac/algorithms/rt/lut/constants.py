"""Constants and unit-boundary helpers for LUT backend modules."""

from __future__ import annotations

import numpy as np

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


def state_tcwv_to_lut(tcwv_cm: np.ndarray | float) -> np.ndarray | float:
    """SIAC state tcwv (cm precipitable water) -> spectral-LUT axis units (mm)."""
    return tcwv_cm * TCWV_CM_TO_LUT_MM


def state_tco3_to_lut(tco3_atmcm: np.ndarray | float) -> np.ndarray | float:
    """SIAC state tco3 (atm-cm) -> spectral-LUT ozone axis units (Dobson)."""
    return tco3_atmcm * TCO3_ATMCM_TO_DU


def normalize_compact_ozone(tco3: np.ndarray, ozone_axis: np.ndarray) -> np.ndarray:
    """Legacy COMPACT-LUT ozone normalisation by magnitude heuristic.

    The spectral-LUT schema declares its ozone axis in DU, so the spectral
    path converts unconditionally (:func:`state_tco3_to_lut`). Compact
    xap/xbp/xcp LUTs predate that schema and may index ozone in either
    atm-cm (~0.3) or DU (~300); with no declared units, the only safe rule
    is the magnitude heuristic below. New LUTs must declare units and use
    the explicit conversion instead of extending this fallback.
    """
    ozone = np.asarray(tco3, dtype=np.float32)
    finite_ozone = ozone[np.isfinite(ozone)]
    if (
        ozone_axis.size
        and finite_ozone.size
        and np.nanmax(np.abs(ozone_axis)) > 20
        and np.nanmax(np.abs(finite_ozone)) < 10
    ):
        return ozone * np.float32(TCO3_ATMCM_TO_DU)
    return ozone
