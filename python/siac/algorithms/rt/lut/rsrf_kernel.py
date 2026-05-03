"""LUT-aligned kernels derived from canonical relative spectral responses."""

from __future__ import annotations

import hashlib
from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from siac.domain.spectral import RelativeSpectralResponse


def _trapezoid(y: np.ndarray, x: np.ndarray) -> float:
    """Return trapezoid integration for a spectral curve."""
    return float(np.trapezoid(y, x))


@dataclass(frozen=True)
class AlignedRSRFKernel:
    """RSRF interpolated onto a specific LUT wavelength axis slice."""

    sensor_id: str
    satellite_id: str
    band_name: str
    lut_id: str
    wavelength_axis_hash: str
    start_index: int
    end_index: int
    wavelengths_nm: np.ndarray
    response_on_lut: np.ndarray
    solar_weighted_response_on_lut: np.ndarray | None = None


def build_aligned_rsrf_kernel(
    rsrf: RelativeSpectralResponse,
    *,
    lut_wavelengths_nm: np.ndarray,
    lut_id: str,
    solar_irradiance: np.ndarray | None = None,
    support_padding: int = 1,
) -> AlignedRSRFKernel:
    """Interpolate an RSRF onto the active LUT wavelength axis and trim to support."""
    wavelengths = np.asarray(lut_wavelengths_nm, dtype=np.float32).reshape(-1)
    if wavelengths.ndim != 1 or wavelengths.size < 2:
        raise ValueError("LUT wavelength axis must be 1-D with at least two samples")
    if np.any(np.diff(wavelengths) <= 0.0):
        raise ValueError("LUT wavelength axis must be strictly increasing")
    if support_padding < 0:
        raise ValueError("support_padding must be non-negative")

    support_mask = (wavelengths >= rsrf.support_min_nm) & (wavelengths <= rsrf.support_max_nm)
    if np.any(support_mask):
        support_idx = np.flatnonzero(support_mask)
        start_index = max(int(support_idx[0]) - support_padding, 0)
        end_index = min(int(support_idx[-1]) + support_padding + 1, wavelengths.size)
    else:
        left = max(
            int(np.searchsorted(wavelengths, rsrf.support_min_nm, side="left")) - support_padding,
            0,
        )
        right = min(
            int(np.searchsorted(wavelengths, rsrf.support_max_nm, side="right")) + support_padding,
            wavelengths.size,
        )
        if left == right:
            right = min(left + 1, wavelengths.size)
            left = max(right - 1, 0)
        start_index = left
        end_index = right

    wl_slice = wavelengths[start_index:end_index]
    response_interp = np.interp(
        wl_slice,
        rsrf.wavelengths_nm,
        rsrf.response,
        left=0.0,
        right=0.0,
    ).astype(np.float32)
    response_area = _trapezoid(response_interp, wl_slice)
    if response_area <= 0.0:
        raise ValueError("Aligned RSRF support has zero overlap with LUT wavelength axis")
    response_interp = response_interp / response_area

    solar_weighted = None
    if solar_irradiance is not None:
        solar = np.asarray(solar_irradiance, dtype=np.float32).reshape(-1)
        if solar.shape != wavelengths.shape:
            raise ValueError("Solar irradiance must match the full LUT wavelength axis")
        weighted = response_interp * solar[start_index:end_index]
        weighted_area = _trapezoid(weighted, wl_slice)
        if weighted_area > 0.0:
            solar_weighted = (weighted / weighted_area).astype(np.float32)

    axis_hash = hashlib.sha1(wavelengths.tobytes(), usedforsecurity=False).hexdigest()
    return AlignedRSRFKernel(
        sensor_id=rsrf.sensor_id,
        satellite_id=rsrf.satellite_id,
        band_name=rsrf.band_name,
        lut_id=lut_id,
        wavelength_axis_hash=axis_hash,
        start_index=start_index,
        end_index=end_index,
        wavelengths_nm=wl_slice,
        response_on_lut=response_interp.astype(np.float32),
        solar_weighted_response_on_lut=solar_weighted,
    )
