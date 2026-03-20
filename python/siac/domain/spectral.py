"""Canonical spectral response types."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from numpy.typing import NDArray


def _trapezoid(y: NDArray, x: NDArray) -> float:
    """Compatibility wrapper for NumPy 1.x/2.x integration API."""
    if hasattr(np, "trapezoid"):
        return float(np.trapezoid(y, x))
    return float(np.add.reduce((y[1:] + y[:-1]) * np.diff(x) * 0.5))


def _fwhm(wavelengths_nm: NDArray, response: NDArray) -> float | None:
    """Estimate full-width half-maximum from a sampled response curve."""
    if wavelengths_nm.size < 2:
        return None

    peak = float(np.nanmax(response))
    if not np.isfinite(peak) or peak <= 0.0:
        return None
    half = 0.5 * peak
    above = np.flatnonzero(response >= half)
    if above.size == 0:
        return None

    left_idx = int(above[0])
    right_idx = int(above[-1])

    left = float(wavelengths_nm[left_idx])
    if left_idx > 0:
        x0 = float(wavelengths_nm[left_idx - 1])
        x1 = float(wavelengths_nm[left_idx])
        y0 = float(response[left_idx - 1])
        y1 = float(response[left_idx])
        if y1 != y0:
            left = float(np.interp(half, np.array([y0, y1]), np.array([x0, x1])))

    right = float(wavelengths_nm[right_idx])
    if right_idx < wavelengths_nm.size - 1:
        x0 = float(wavelengths_nm[right_idx])
        x1 = float(wavelengths_nm[right_idx + 1])
        y0 = float(response[right_idx])
        y1 = float(response[right_idx + 1])
        if y0 != y1:
            right = float(np.interp(half, np.array([y1, y0]), np.array([x1, x0])))

    width = right - left
    return width if width >= 0.0 else None


@dataclass(frozen=True)
class SpectralResponseFunction:
    """Canonical, area-normalized spectral response function."""

    sensor_id: str
    satellite_id: str
    band_name: str
    wavelengths_nm: NDArray[np.floating]
    response: NDArray[np.floating]
    response_raw: NDArray[np.floating] | None = None
    source_id: str | None = None
    source_version: str | None = None
    source_url: str | None = None
    centre_wavelength_nm: float | None = None
    effective_wavelength_nm: float | None = None
    fwhm_nm: float | None = None

    def __post_init__(self) -> None:
        wavelengths = np.asarray(self.wavelengths_nm, dtype=np.float32)
        response = np.asarray(self.response, dtype=np.float32)
        raw = None if self.response_raw is None else np.asarray(self.response_raw, dtype=np.float32)

        if wavelengths.ndim != 1 or response.ndim != 1:
            raise ValueError("SRF wavelengths and response must be 1-D")
        if wavelengths.size != response.size:
            raise ValueError("SRF wavelengths and response must have equal length")
        if wavelengths.size < 2:
            raise ValueError("SRF must contain at least two samples")
        if not np.isfinite(wavelengths).all() or not np.isfinite(response).all():
            raise ValueError("SRF wavelengths and response must be finite")
        if np.any(np.diff(wavelengths) <= 0.0):
            raise ValueError("SRF wavelengths must be strictly increasing")
        if np.any(response < 0.0):
            raise ValueError("SRF response must be non-negative")
        if raw is not None and raw.shape != wavelengths.shape:
            raise ValueError("SRF raw response must match wavelengths")

        area = _trapezoid(response, wavelengths)
        if not np.isfinite(area) or area <= 0.0:
            raise ValueError("SRF response must integrate to a positive value")
        if not np.isclose(area, 1.0, rtol=1e-4, atol=1e-6):
            raise ValueError("Canonical SRF response must be area-normalized")

        object.__setattr__(self, "wavelengths_nm", wavelengths)
        object.__setattr__(self, "response", response)
        object.__setattr__(self, "response_raw", raw)

    @property
    def support_min_nm(self) -> float:
        """First wavelength stored in the canonical support."""
        return float(self.wavelengths_nm[0])

    @property
    def support_max_nm(self) -> float:
        """Last wavelength stored in the canonical support."""
        return float(self.wavelengths_nm[-1])

    @classmethod
    def from_tabulated(
        cls,
        *,
        sensor_id: str,
        satellite_id: str,
        band_name: str,
        wavelengths_nm: NDArray[np.floating],
        response: NDArray[np.floating],
        source_id: str | None = None,
        source_version: str | None = None,
        source_url: str | None = None,
    ) -> SpectralResponseFunction:
        """Build a canonical SRF from a raw tabulated response."""
        wavelengths = np.asarray(wavelengths_nm, dtype=np.float32).reshape(-1)
        raw_response = np.asarray(response, dtype=np.float32).reshape(-1)
        if wavelengths.size != raw_response.size:
            raise ValueError("SRF wavelengths and response must have equal length")
        if wavelengths.size < 2:
            raise ValueError("SRF must contain at least two samples")
        if not np.isfinite(wavelengths).all() or not np.isfinite(raw_response).all():
            raise ValueError("SRF wavelengths and response must be finite")

        order = np.argsort(wavelengths, kind="stable")
        wavelengths = wavelengths[order]
        raw_response = raw_response[order]

        unique_wavelengths, unique_idx = np.unique(wavelengths, return_index=True)
        wavelengths = unique_wavelengths.astype(np.float32, copy=False)
        raw_response = raw_response[unique_idx].astype(np.float32, copy=False)
        raw_response = np.clip(raw_response, 0.0, None)

        nonzero = np.flatnonzero(raw_response > 0.0)
        if nonzero.size == 0:
            raise ValueError("SRF response must include at least one non-zero sample")

        start = max(int(nonzero[0]) - 1, 0)
        stop = min(int(nonzero[-1]) + 2, raw_response.size)
        wavelengths = wavelengths[start:stop]
        raw_response = raw_response[start:stop]

        area = _trapezoid(raw_response, wavelengths)
        if not np.isfinite(area) or area <= 0.0:
            raise ValueError("SRF response must integrate to a positive value")
        response_norm = (raw_response / area).astype(np.float32)

        peak_idx = int(np.argmax(response_norm))
        centre = float(wavelengths[peak_idx])
        effective = _trapezoid(wavelengths * response_norm, wavelengths)
        width = _fwhm(wavelengths, response_norm)

        return cls(
            sensor_id=sensor_id,
            satellite_id=satellite_id,
            band_name=band_name,
            wavelengths_nm=wavelengths.astype(np.float32, copy=False),
            response=response_norm,
            response_raw=raw_response.astype(np.float32, copy=False),
            source_id=source_id,
            source_version=source_version,
            source_url=source_url,
            centre_wavelength_nm=centre,
            effective_wavelength_nm=float(effective),
            fwhm_nm=None if width is None else float(width),
        )


__all__ = ["SpectralResponseFunction"]
