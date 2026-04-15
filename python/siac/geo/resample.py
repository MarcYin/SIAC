"""Shared resampling utilities for aligning xarray fields to a common grid.

These functions are used across the pipeline orchestrator, the atmospheric
corrector, and the assembly layer to resample 2-D xarray DataArrays between
different spatial grids (e.g. atmospheric prior grid -> solver grid ->
native pixel grid).
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING, Any

import numpy as np
import xarray as xr
from scipy import ndimage

from siac.runtime.models import copy_spatial_metadata_like

if TYPE_CHECKING:
    from siac.runtime import RTCoefficients

logger = logging.getLogger(__name__)

# xarray .interp() is precise but extremely slow for large arrays (minutes
# for 10980×10980).  Use scipy.ndimage.zoom instead when pixel count exceeds
# this threshold.
_INTERP_PIXEL_LIMIT = 1_000_000  # ~1000×1000


# ---------------------------------------------------------------------------
# Grid comparison
# ---------------------------------------------------------------------------

def shares_template_grid(field: xr.DataArray, template: xr.DataArray) -> bool:
    """Return True when *field* already sits on the same spatial grid as *template*.

    Checks shape, dims, and coordinate arrays for equality.
    """
    if tuple(field.shape) != tuple(template.shape):
        return False
    if tuple(field.dims) != tuple(template.dims):
        return False
    for dim in template.dims:
        template_has_coord = dim in template.coords
        field_has_coord = dim in field.coords
        if template_has_coord != field_has_coord:
            return False
        if not template_has_coord:
            continue
        if not np.array_equal(
            np.asarray(template.coords[dim].values),
            np.asarray(field.coords[dim].values),
            equal_nan=True,
        ):
            return False
    return True


# ---------------------------------------------------------------------------
# Axis resolution helper
# ---------------------------------------------------------------------------

def axis_resolution(values: np.ndarray) -> float | None:
    """Return the median absolute step size along a coordinate axis."""
    if values.size <= 1:
        return None
    diffs = np.abs(np.diff(np.asarray(values, dtype=np.float64)))
    finite = diffs[np.isfinite(diffs) & (diffs > 0.0)]
    if finite.size == 0:
        return None
    return float(np.median(finite))


def _can_coord_interp(source: xr.DataArray, template: xr.DataArray) -> bool:
    """Return True if both arrays are 2-D with matching named dims and coords."""
    return (
        len(source.dims) == 2
        and source.dims == template.dims
        and all(dim in source.coords for dim in source.dims)
        and all(dim in template.coords for dim in template.dims)
    )


# ---------------------------------------------------------------------------
# Field resampling (continuous data)
# ---------------------------------------------------------------------------

def resample_field_to_template(
    field: xr.DataArray,
    template: xr.DataArray,
    *,
    copy_metadata: bool = True,
    gap_fill: bool = True,
) -> xr.DataArray:
    """Resample a 2-D field to the spatial grid of *template*.

    Strategy:
    1. If the grids already match, return (optionally with metadata copy).
    2. If both have coordinate arrays, use xarray linear interpolation.
    3. Fall back to ``scipy.ndimage.zoom`` with bilinear interpolation.

    When *gap_fill* is True (default), any remaining NaN pixels are filled
    via nearest-neighbour interpolation then source-mean fallback.
    """
    if shares_template_grid(field, template):
        return copy_spatial_metadata_like(field, template) if copy_metadata else field

    resampled: xr.DataArray | None = None
    src_pixels = field.sizes[field.dims[0]] * field.sizes[field.dims[1]] if field.ndim == 2 else 0
    if _can_coord_interp(field, template) and src_pixels <= _INTERP_PIXEL_LIMIT:
        try:
            interpolated = field.interp(
                coords={dim: template.coords[dim] for dim in template.dims},
                method="linear",
            )
            resampled = copy_spatial_metadata_like(interpolated, template) if copy_metadata else interpolated
        except Exception as exc:
            logger.warning(
                "Linear interpolation failed for field %s (%s); falling back to zoom resampling",
                getattr(field, "name", "?"),
                exc,
            )

    if resampled is None:
        src = np.asarray(field.values, dtype=np.float32)
        if src.ndim != 2 or len(template.dims) != 2:
            return field

        h_out = int(template.sizes[template.dims[0]])
        w_out = int(template.sizes[template.dims[1]])
        if src.shape[0] == 0 or src.shape[1] == 0:
            out: np.ndarray[Any, Any] = np.full((h_out, w_out), np.nan, dtype=np.float32)
        else:
            out = ndimage.zoom(src, (h_out / src.shape[0], w_out / src.shape[1]), order=1)
            out = out[:h_out, :w_out]
            if out.shape != (h_out, w_out):
                padded: np.ndarray[Any, Any] = np.full((h_out, w_out), np.nan, dtype=np.float32)
                padded[: out.shape[0], : out.shape[1]] = out
                out = padded

        resampled = xr.DataArray(
            out,
            dims=template.dims,
            coords={dim: template.coords[dim] for dim in template.dims if dim in template.coords},
            attrs=field.attrs,
        )
        if copy_metadata:
            resampled = copy_spatial_metadata_like(resampled, template)

    if gap_fill:
        resampled = fill_nonfinite_like_template(resampled, field, template)

    return resampled


# ---------------------------------------------------------------------------
# Mask resampling (boolean data, conservative / dilated)
# ---------------------------------------------------------------------------

def resample_mask_to_template(mask: xr.DataArray, template: xr.DataArray) -> xr.DataArray:
    """Resample a boolean mask to *template* grid, dilating to avoid false negatives.

    When downsampling from a finer grid, a maximum-filter is applied before
    interpolation so that masked (True) pixels are not lost.
    """
    if shares_template_grid(mask, template):
        return copy_spatial_metadata_like(mask.astype(bool), template)

    src_pixels = mask.sizes[mask.dims[0]] * mask.sizes[mask.dims[1]] if mask.ndim == 2 else 0
    if _can_coord_interp(mask, template) and src_pixels <= _INTERP_PIXEL_LIMIT:
        source = np.asarray(mask.values, dtype=np.float32)
        factor_y = 1
        factor_x = 1
        src_y_res = axis_resolution(np.asarray(mask.coords[mask.dims[0]].values))
        src_x_res = axis_resolution(np.asarray(mask.coords[mask.dims[1]].values))
        dst_y_res = axis_resolution(np.asarray(template.coords[template.dims[0]].values))
        dst_x_res = axis_resolution(np.asarray(template.coords[template.dims[1]].values))
        if src_y_res is not None and dst_y_res is not None and dst_y_res > src_y_res:
            factor_y = max(1, int(np.ceil(dst_y_res / src_y_res)))
        if src_x_res is not None and dst_x_res is not None and dst_x_res > src_x_res:
            factor_x = max(1, int(np.ceil(dst_x_res / src_x_res)))
        if factor_y > 1 or factor_x > 1:
            source = ndimage.maximum_filter(source, size=(factor_y, factor_x), mode="nearest")

        remapped = xr.DataArray(
            source,
            dims=mask.dims,
            coords={dim: mask.coords[dim] for dim in mask.dims},
            attrs=mask.attrs,
        )
        try:
            interpolated = remapped.interp(
                coords={dim: template.coords[dim] for dim in template.dims},
                method="nearest",
            )
            values = np.asarray(interpolated.values, dtype=np.float32)
            values = np.where(np.isfinite(values), values, 1.0)
            return copy_spatial_metadata_like(
                xr.DataArray(
                    values > 0.5,
                    dims=template.dims,
                    coords={dim: template.coords[dim] for dim in template.dims if dim in template.coords},
                    attrs=mask.attrs,
                ),
                template,
            )
        except Exception as exc:
            logger.warning(
                "Mask coordinate-based resampling failed (%s); falling back to zoom resampling",
                exc,
            )

    # Zoom-based fallback
    src = np.asarray(mask.values, dtype=np.float32)
    h_out = int(template.sizes[template.dims[0]])
    w_out = int(template.sizes[template.dims[1]])
    factor_y = max(1, src.shape[0] // h_out) if src.ndim == 2 and h_out > 0 else 1
    factor_x = max(1, src.shape[1] // w_out) if src.ndim == 2 and w_out > 0 else 1
    if src.ndim != 2:
        out = np.ones((h_out, w_out), dtype=np.float32)
    else:
        dilated = ndimage.maximum_filter(src, size=(factor_y, factor_x), mode="nearest")
        out = ndimage.zoom(dilated, (h_out / dilated.shape[0], w_out / dilated.shape[1]), order=0)
        out = out[:h_out, :w_out]
        if out.shape != (h_out, w_out):
            padded: np.ndarray[Any, Any] = np.ones((h_out, w_out), dtype=np.float32)
            padded[: out.shape[0], : out.shape[1]] = out
            out = padded
    return copy_spatial_metadata_like(
        xr.DataArray(
            out > 0.5,
            dims=template.dims,
            coords={dim: template.coords[dim] for dim in template.dims if dim in template.coords},
            attrs=mask.attrs,
        ),
        template,
    )


# ---------------------------------------------------------------------------
# RT coefficient resampling
# ---------------------------------------------------------------------------

def resample_coefficients_to_template(
    coeffs: RTCoefficients,
    template: xr.DataArray,
) -> RTCoefficients:
    """Resample all fields of an RTCoefficients bundle to *template* grid."""
    from siac.runtime import RTCoefficients as RTCoefficientsClass

    def _resample_optional(field: xr.DataArray | None) -> xr.DataArray | None:
        if field is None:
            return None
        if len(field.dims) == 2:
            return resample_field_to_template(field, template)
        if len(field.dims) == 3 and "param" in field.dims:
            param_values = field.coords["param"].values if "param" in field.coords else np.arange(field.sizes["param"])
            stacked = xr.concat(
                [resample_field_to_template(field.sel(param=param, drop=True), template) for param in param_values],
                dim="param",
            )
            result: xr.DataArray = stacked.assign_coords(param=param_values).transpose("param", *template.dims)
            return result
        return field

    return RTCoefficientsClass(
        xap=resample_field_to_template(coeffs.xap, template),
        xbp=resample_field_to_template(coeffs.xbp, template),
        xcp=resample_field_to_template(coeffs.xcp, template),
        d_xap=_resample_optional(coeffs.d_xap),
        d_xbp=_resample_optional(coeffs.d_xbp),
        d_xcp=_resample_optional(coeffs.d_xcp),
    )


# ---------------------------------------------------------------------------
# NaN gap-fill helper
# ---------------------------------------------------------------------------

def fill_nonfinite_like_template(
    field: xr.DataArray,
    source: xr.DataArray,
    template: xr.DataArray,
    *,
    fallback_value: float = 0.0,
    field_name: str | None = None,
) -> xr.DataArray:
    """Fill NaN values in *field* using nearest-neighbour from *source*, then a constant."""
    values = np.asarray(field.values, dtype=np.float32)
    if np.all(np.isfinite(values)):
        return field

    filled = values.copy()
    src_pixels = source.sizes[source.dims[0]] * source.sizes[source.dims[1]] if source.ndim == 2 else 0
    if _can_coord_interp(source, template) and src_pixels <= _INTERP_PIXEL_LIMIT:
        try:
            nearest = source.interp(
                coords={dim: template.coords[dim] for dim in template.dims},
                method="nearest",
            )
            nearest_values = np.asarray(nearest.values, dtype=np.float32)
            filled = np.where(np.isfinite(filled), filled, nearest_values)
        except Exception:
            pass

    if not np.all(np.isfinite(filled)):
        source_values = np.asarray(source.values, dtype=np.float32)
        finite_source = source_values[np.isfinite(source_values)]
        fill_value = float(finite_source.mean()) if finite_source.size else float(fallback_value)
        n_nan = int(np.sum(~np.isfinite(filled)))
        logger.warning(
            "Filling %d remaining NaN pixels in field %s with source mean %.4g",
            n_nan,
            field_name or getattr(field, "name", "?"),
            fill_value,
        )
        filled = np.where(np.isfinite(filled), filled, np.float32(fill_value))

    result: xr.DataArray = field.copy(data=filled)
    return result


def resample_field_for_correction(
    field: xr.DataArray,
    template: xr.DataArray,
    *,
    fallback_value: float = 0.0,
    field_name: str | None = None,
) -> xr.DataArray:
    """Resample an atmosphere field to the correction grid and guarantee finiteness."""
    return fill_nonfinite_like_template(
        resample_field_to_template(field, template, copy_metadata=False, gap_fill=False),
        field,
        template,
        fallback_value=fallback_value,
        field_name=field_name,
    )


def should_resample_for_policy(
    current_resolution_m: float,
    target_resolution_m: float,
    *,
    policy: str = "auto",
    allow_upsample: bool = False,
) -> bool:
    """Decide whether resampling is needed given a resolution policy.

    Args:
        current_resolution_m: Current pixel resolution in metres.
        target_resolution_m: Desired resolution in metres.
        policy: ``"auto"`` (resample only when beneficial) or ``"force"``.
        allow_upsample: If True, upsampling from coarser to finer is allowed.

    Returns:
        True if the data should be resampled.
    """
    if target_resolution_m <= 0:
        raise ValueError(f"target_resolution_m must be > 0, got {target_resolution_m}")
    if policy == "force":
        return abs(current_resolution_m - target_resolution_m) > 1e-6
    if policy == "auto":
        if current_resolution_m < target_resolution_m - 1e-6:
            return True  # downsample finer to coarser
        return current_resolution_m > target_resolution_m + 1e-6 and allow_upsample
    raise ValueError(f"resolution_policy must be 'auto' or 'force', got {policy!r}")


__all__ = [
    "axis_resolution",
    "should_resample_for_policy",
    "fill_nonfinite_like_template",
    "resample_coefficients_to_template",
    "resample_field_for_correction",
    "resample_field_to_template",
    "resample_mask_to_template",
    "shares_template_grid",
]
