"""Route-B monthly best-pixel BRDF composites."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import xarray as xr


@dataclass(frozen=True)
class MonthlyBestPixelComposite:
    """Monthly best-pixel composite of BRDF-simulated reflectance."""

    reflectance: xr.DataArray
    quality: xr.DataArray
    sample_index: xr.DataArray
    year: int
    month: int


def build_monthly_best_pixel_composite(
    reflectance: xr.DataArray,
    quality: xr.DataArray,
    *,
    year: int,
    month: int,
) -> MonthlyBestPixelComposite:
    """Select the best-quality BRDF-simulated sample for each pixel in the month."""
    if reflectance.ndim != 4 or reflectance.dims[:2] != ("time", "band"):
        raise ValueError("reflectance must have dims ('time', 'band', 'y', 'x')")
    if quality.ndim != 3 or quality.dims[0] != "time":
        raise ValueError("quality must have dims ('time', 'y', 'x')")
    if reflectance.sizes.get("time") != quality.sizes.get("time"):
        raise ValueError("reflectance and quality must share the same time axis")
    if reflectance.sizes.get("y") != quality.sizes.get("y") or reflectance.sizes.get("x") != quality.sizes.get("x"):
        raise ValueError("reflectance and quality must share the same spatial shape")

    refl_values = np.asarray(reflectance.values, dtype=np.float32)
    quality_values = np.asarray(quality.values, dtype=np.float32)

    valid_refl = np.all(np.isfinite(refl_values), axis=1)
    valid_quality = np.isfinite(quality_values)
    valid = valid_refl & valid_quality

    masked_quality = np.where(valid, quality_values, np.inf)
    sample_index = np.argmin(masked_quality, axis=0).astype(np.int16)
    has_valid = np.isfinite(np.min(masked_quality, axis=0))

    gather_index = sample_index[np.newaxis, np.newaxis, ...]
    selected = np.take_along_axis(refl_values, gather_index, axis=0)[0]
    selected_quality = np.take_along_axis(quality_values[:, np.newaxis, ...], sample_index[np.newaxis, np.newaxis, ...], axis=0)[0, 0]

    selected = np.where(has_valid[np.newaxis, ...], selected, np.nan)
    selected_quality = np.where(has_valid, selected_quality, np.nan)
    sample_index = np.where(has_valid, sample_index, -1).astype(np.int16)

    coords = {
        "band": reflectance.coords["band"],
        "y": reflectance.coords["y"],
        "x": reflectance.coords["x"],
    }
    return MonthlyBestPixelComposite(
        reflectance=xr.DataArray(selected, dims=["band", "y", "x"], coords=coords),
        quality=xr.DataArray(
            selected_quality,
            dims=["y", "x"],
            coords={"y": reflectance.coords["y"], "x": reflectance.coords["x"]},
        ),
        sample_index=xr.DataArray(
            sample_index,
            dims=["y", "x"],
            coords={"y": reflectance.coords["y"], "x": reflectance.coords["x"]},
        ),
        year=int(year),
        month=int(month),
    )
