"""Pure helpers for spectral-response curve preparation."""

from __future__ import annotations

from typing import TYPE_CHECKING, TypeAlias

import numpy as np
from numpy import typing as npt

from siac.adapters.rsrf import band_convolution_weights

if TYPE_CHECKING:
    from collections.abc import Sequence

    from siac.domain import SensorBand


Float32Array: TypeAlias = npt.NDArray[np.float32]

VISIBLE_UPPER_NM = 700.0
SEGMENT_RANGES_NM: dict[str, tuple[float, float]] = {
    "vnir": (400.0, 1000.0),
    "swir": (800.0, 2500.0),
}


def trapezoid(y: np.ndarray, x: np.ndarray) -> float:
    return float(np.trapezoid(y, x))


def normalized_band_response(band: SensorBand, wavelengths_nm: Float32Array) -> Float32Array:
    response = np.asarray(band_convolution_weights(band, wavelengths_nm), dtype=np.float32)
    if not np.all(np.isfinite(response)) or float(np.sum(response, dtype=np.float64)) <= 0.0:
        raise ValueError(f"Band {band.name!r} has zero support on the requested wavelength grid")
    return response


def classify_band_region(band: SensorBand) -> str:
    return "visible" if band.center_wavelength < VISIBLE_UPPER_NM else "infrared"


def segment_for_band(band: SensorBand) -> str:
    return "swir" if float(band.center_wavelength) >= 1000.0 else "vnir"


def segmentize_curve(
    wavelengths_nm: Float32Array,
    response: Float32Array,
    *,
    segment: str,
) -> tuple[Float32Array, Float32Array]:
    segment_min_nm, segment_max_nm = SEGMENT_RANGES_NM[segment]
    segment_wavelengths: Float32Array = np.arange(
        segment_min_nm,
        segment_max_nm + 1.0,
        1.0,
        dtype=np.float32,
    )
    segment_response = np.interp(
        segment_wavelengths,
        np.asarray(wavelengths_nm, dtype=np.float32),
        np.asarray(response, dtype=np.float32),
        left=0.0,
        right=0.0,
    ).astype(np.float32, copy=False)
    segment_response = np.clip(segment_response, 0.0, None)

    positive = np.flatnonzero(segment_response > 0.0)
    if positive.size == 0:
        raise ValueError(
            f"Band response does not overlap the supported {segment!r} segment range "
            f"{segment_min_nm:.0f}-{segment_max_nm:.0f} nm"
        )

    start = max(int(positive[0]) - 1, 0)
    stop = min(int(positive[-1]) + 2, segment_response.size)
    return segment_wavelengths[start:stop], segment_response[start:stop]


def primary_nir_band_index(bands: Sequence[SensorBand]) -> int | None:
    candidates = [
        (idx, abs(float(band.center_wavelength) - 865.0))
        for idx, band in enumerate(bands)
        if 750.0 <= float(band.center_wavelength) <= 1000.0
    ]
    if not candidates:
        return None
    candidates.sort(key=lambda item: (item[1], item[0]))
    return int(candidates[0][0])


def canonicalize_curve(
    wavelengths_nm: Float32Array,
    response: Float32Array,
) -> tuple[Float32Array, Float32Array]:
    wavelengths = np.asarray(wavelengths_nm, dtype=np.float32).reshape(-1)
    weights = np.asarray(response, dtype=np.float32).reshape(-1)
    if wavelengths.size != weights.size or wavelengths.size < 2:
        raise ValueError("Spectral response curves require at least two wavelength samples")
    order = np.argsort(wavelengths, kind="stable")
    wavelengths = wavelengths[order]
    weights = np.clip(weights[order], 0.0, None)
    unique_wavelengths, unique_idx = np.unique(wavelengths, return_index=True)
    wavelengths = unique_wavelengths.astype(np.float32, copy=False)
    weights = weights[unique_idx].astype(np.float32, copy=False)
    positive = np.flatnonzero(weights > 0.0)
    if positive.size == 0:
        raise ValueError("Spectral response curves must contain at least one positive response sample")
    start = max(int(positive[0]) - 1, 0)
    stop = min(int(positive[-1]) + 2, weights.size)
    return wavelengths[start:stop], weights[start:stop]
