"""Route-B monthly best-pixel BRDF composites."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Any, cast

import numpy as np
import xarray as xr

from siac.runtime import BRDFKernelWeights
from siac.runtime.models import copy_spatial_metadata_like

if TYPE_CHECKING:
    from siac.domain import SensorBand


@dataclass(frozen=True)
class MonthlyBestPixelComposite:
    """Monthly best-pixel composite of BRDF-simulated reflectance."""

    reflectance: xr.DataArray
    quality: xr.DataArray
    sample_index: xr.DataArray
    year: int
    month: int


@dataclass(frozen=True)
class MonthlyKernelWeightComposite:
    """Monthly best-pixel composite of BRDF kernel weights."""

    kernels: BRDFKernelWeights
    quality: xr.DataArray
    sample_index: xr.DataArray
    year: int
    month: int


@dataclass(frozen=True)
class MonthlyCompositeCollection:
    """Prepared monthly composites plus the source spectral basis metadata."""

    composites: tuple[MonthlyBestPixelComposite | MonthlyKernelWeightComposite, ...]
    source_bands: tuple[SensorBand, ...]
    source_name: str | None = None


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
    reflectance_reference = reflectance.isel(time=0, drop=True)
    quality_reference = quality.isel(time=0, drop=True)
    return MonthlyBestPixelComposite(
        reflectance=copy_spatial_metadata_like(
            xr.DataArray(selected, dims=["band", "y", "x"], coords=coords),
            reflectance_reference,
        ),
        quality=copy_spatial_metadata_like(
            xr.DataArray(
                selected_quality,
                dims=["y", "x"],
                coords={"y": reflectance.coords["y"], "x": reflectance.coords["x"]},
            ),
            quality_reference,
        ),
        sample_index=copy_spatial_metadata_like(
            xr.DataArray(
                sample_index,
                dims=["y", "x"],
                coords={"y": reflectance.coords["y"], "x": reflectance.coords["x"]},
            ),
            quality_reference,
        ),
        year=int(year),
        month=int(month),
    )


def build_monthly_best_pixel_kernel_composite(
    temporal_weights: BRDFKernelWeights,
    quality: xr.DataArray,
    *,
    year: int,
    month: int,
) -> MonthlyKernelWeightComposite:
    """Select the best-quality BRDF kernel-weight sample for each pixel in the month."""
    if quality.ndim != 3 or quality.dims[0] != "time":
        raise ValueError("quality must have dims ('time', 'y', 'x')")

    def _validate(name: str, data: xr.DataArray) -> None:
        if data.ndim != 4 or data.dims[:2] != ("time", "band"):
            raise ValueError(f"{name} must have dims ('time', 'band', 'y', 'x')")
        if data.sizes.get("time") != quality.sizes.get("time"):
            raise ValueError(f"{name} and quality must share the same time axis")
        if data.sizes.get("y") != quality.sizes.get("y") or data.sizes.get("x") != quality.sizes.get("x"):
            raise ValueError(f"{name} and quality must share the same spatial shape")

    _validate("f0", temporal_weights.f0)
    _validate("f1", temporal_weights.f1)
    _validate("f2", temporal_weights.f2)
    _validate("f0_unc", temporal_weights.f0_unc)
    _validate("f1_unc", temporal_weights.f1_unc)
    _validate("f2_unc", temporal_weights.f2_unc)
    if temporal_weights.reflectance_unc is not None:
        _validate("reflectance_unc", temporal_weights.reflectance_unc)

    quality_values = np.asarray(quality.values, dtype=np.float32)

    def _valid_weight_stack(data: xr.DataArray) -> np.ndarray:
        return cast(
            "np.ndarray[Any, np.dtype[np.bool_]]",
            np.all(np.isfinite(np.asarray(data.values, dtype=np.float32)), axis=1),
        )

    valid = _valid_weight_stack(temporal_weights.f0)
    valid &= _valid_weight_stack(temporal_weights.f1)
    valid &= _valid_weight_stack(temporal_weights.f2)
    valid &= _valid_weight_stack(temporal_weights.f0_unc)
    valid &= _valid_weight_stack(temporal_weights.f1_unc)
    valid &= _valid_weight_stack(temporal_weights.f2_unc)
    if temporal_weights.reflectance_unc is not None:
        valid &= _valid_weight_stack(temporal_weights.reflectance_unc)
    valid &= np.isfinite(quality_values)

    masked_quality = np.where(valid, quality_values, np.inf)
    sample_index = np.argmin(masked_quality, axis=0).astype(np.int16)
    has_valid = np.isfinite(np.min(masked_quality, axis=0))

    def _select(data: xr.DataArray) -> xr.DataArray:
        values = np.asarray(data.values, dtype=np.float32)
        gather_index = sample_index[np.newaxis, np.newaxis, ...]
        selected = np.take_along_axis(values, gather_index, axis=0)[0]
        selected = np.where(has_valid[np.newaxis, ...], selected, np.nan)
        return copy_spatial_metadata_like(
            xr.DataArray(
                selected,
                dims=["band", "y", "x"],
                coords={
                    "band": data.coords["band"],
                    "y": data.coords["y"],
                    "x": data.coords["x"],
                },
            ),
            data.isel(time=0, drop=True),
        )

    selected_quality = np.take_along_axis(
        quality_values[:, np.newaxis, ...],
        sample_index[np.newaxis, np.newaxis, ...],
        axis=0,
    )[0, 0]
    selected_quality = np.where(has_valid, selected_quality, np.nan)
    sample_index = np.where(has_valid, sample_index, -1).astype(np.int16)

    reflectance_unc = (
        _select(temporal_weights.reflectance_unc)
        if temporal_weights.reflectance_unc is not None
        else None
    )
    quality_reference = quality.isel(time=0, drop=True)
    return MonthlyKernelWeightComposite(
        kernels=BRDFKernelWeights(
            f0=_select(temporal_weights.f0),
            f1=_select(temporal_weights.f1),
            f2=_select(temporal_weights.f2),
            f0_unc=_select(temporal_weights.f0_unc),
            f1_unc=_select(temporal_weights.f1_unc),
            f2_unc=_select(temporal_weights.f2_unc),
            reflectance_unc=reflectance_unc,
        ),
        quality=copy_spatial_metadata_like(
            xr.DataArray(
                selected_quality,
                dims=["y", "x"],
                coords={"y": quality.coords["y"], "x": quality.coords["x"]},
            ),
            quality_reference,
        ),
        sample_index=copy_spatial_metadata_like(
            xr.DataArray(
                sample_index,
                dims=["y", "x"],
                coords={"y": quality.coords["y"], "x": quality.coords["x"]},
            ),
            quality_reference,
        ),
        year=int(year),
        month=int(month),
    )
