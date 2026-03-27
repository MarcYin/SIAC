"""Atmospheric correction: TOA to BOA conversion."""
from __future__ import annotations

import time
from typing import TYPE_CHECKING, Any

import numpy as np
import xarray as xr
from scipy import ndimage

from siac.domain.protocols import RTModelBackend
from siac.runtime import (
    AtmosphericState,
    CorrectionDiagnostics,
    CorrectionResult,
    GeometryAngles,
    RTCoefficients,
)
from siac.runtime.models import copy_spatial_metadata_like

if TYPE_CHECKING:
    from siac.domain import SensorConfig

_TOA_BAND_LOADER_ATTR = "_siac_toa_band_loader"


def _shares_template_grid(field: xr.DataArray, template: xr.DataArray) -> bool:
    if tuple(field.shape) != tuple(template.shape):
        return False
    if tuple(field.dims) != tuple(template.dims):
        return False
    for dim in template.dims:
        if dim not in field.coords or dim not in template.coords:
            return False
        if not np.array_equal(np.asarray(field.coords[dim].values), np.asarray(template.coords[dim].values), equal_nan=True):
            return False
    return True


def _axis_resolution(values: np.ndarray) -> float | None:
    if values.size <= 1:
        return None
    diffs = np.abs(np.diff(np.asarray(values, dtype=np.float64)))
    finite = diffs[np.isfinite(diffs) & (diffs > 0.0)]
    if finite.size == 0:
        return None
    return float(np.median(finite))


def _resample_mask_to_template(mask: xr.DataArray, template: xr.DataArray) -> xr.DataArray:
    if _shares_template_grid(mask, template):
        return copy_spatial_metadata_like(mask.astype(bool), template)

    if (
        len(mask.dims) == 2
        and mask.dims == template.dims
        and all(dim in mask.coords for dim in mask.dims)
        and all(dim in template.coords for dim in template.dims)
    ):
        source = np.asarray(mask.values, dtype=np.float32)
        factor_y = 1
        factor_x = 1
        src_y_res = _axis_resolution(np.asarray(mask.coords[mask.dims[0]].values))
        src_x_res = _axis_resolution(np.asarray(mask.coords[mask.dims[1]].values))
        dst_y_res = _axis_resolution(np.asarray(template.coords[template.dims[0]].values))
        dst_x_res = _axis_resolution(np.asarray(template.coords[template.dims[1]].values))
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
        except Exception:
            pass

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


def _resample_field_to_template(field: xr.DataArray, template: xr.DataArray) -> xr.DataArray:
    if _shares_template_grid(field, template):
        return copy_spatial_metadata_like(field, template)

    resampled: xr.DataArray | None = None
    if (
        len(field.dims) == 2
        and field.dims == template.dims
        and all(dim in field.coords for dim in field.dims)
        and all(dim in template.coords for dim in template.dims)
    ):
        try:
            interpolated = field.interp(
                coords={dim: template.coords[dim] for dim in template.dims},
                method="linear",
            )
            resampled = copy_spatial_metadata_like(interpolated, template)
        except Exception:
            pass

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
            out.astype(np.float32),
            dims=template.dims,
            coords={dim: template.coords[dim] for dim in template.dims if dim in template.coords},
            attrs=field.attrs,
        )
        resampled = copy_spatial_metadata_like(resampled, template)

    values = np.asarray(resampled.values, dtype=np.float32)
    if np.all(np.isfinite(values)):
        return resampled

    filled = values.copy()
    if (
        len(field.dims) == 2
        and field.dims == template.dims
        and all(dim in field.coords for dim in field.dims)
        and all(dim in template.coords for dim in template.dims)
    ):
        try:
            nearest = field.interp(
                coords={dim: template.coords[dim] for dim in template.dims},
                method="nearest",
            )
            filled = np.where(np.isfinite(filled), filled, np.asarray(nearest.values, dtype=np.float32))
        except Exception:
            pass

    if not np.all(np.isfinite(filled)):
        source_values = np.asarray(field.values, dtype=np.float32)
        finite_source = source_values[np.isfinite(source_values)]
        fill_value = float(finite_source.mean()) if finite_source.size else 0.0
        filled = np.where(np.isfinite(filled), filled, np.float32(fill_value))

    return resampled.copy(data=filled)


def _resample_coefficients_to_template(
    coeffs: RTCoefficients,
    template: xr.DataArray,
) -> RTCoefficients:
    def _resample_optional(field: xr.DataArray | None) -> xr.DataArray | None:
        if field is None:
            return None
        if len(field.dims) == 2:
            return _resample_field_to_template(field, template)
        if len(field.dims) == 3 and "param" in field.dims:
            param_values = field.coords["param"].values if "param" in field.coords else np.arange(field.sizes["param"])
            stacked = xr.concat(
                [_resample_field_to_template(field.sel(param=param, drop=True), template) for param in param_values],
                dim="param",
            )
            return stacked.assign_coords(param=param_values).transpose("param", *template.dims)
        return field

    return RTCoefficients(
        xap=_resample_field_to_template(coeffs.xap, template),
        xbp=_resample_field_to_template(coeffs.xbp, template),
        xcp=_resample_field_to_template(coeffs.xcp, template),
        d_xap=_resample_optional(coeffs.d_xap),
        d_xbp=_resample_optional(coeffs.d_xbp),
        d_xcp=_resample_optional(coeffs.d_xcp),
    )


class AtmosphericCorrector:
    def __init__(self, rt_model: Any, sensor_config: SensorConfig):
        if not isinstance(rt_model, RTModelBackend):
            raise TypeError(
                f"rt_model must implement RTModelBackend protocol, "
                f"got {type(rt_model).__name__}"
            )
        self.rt_model = rt_model
        self.sensor_config = sensor_config

    def correct(self, toa: xr.Dataset, geometry: GeometryAngles, atmo_state: AtmosphericState,
                cloud_mask: xr.DataArray | None = None) -> CorrectionResult:
        t0 = time.monotonic()
        boa_vars = {}
        invalid_boa_mask: xr.DataArray | None = None
        band_loader = toa.attrs.get(_TOA_BAND_LOADER_ATTR)
        for band_name in self.sensor_config.band_names:
            try:
                band_spec = self.sensor_config.get_band(band_name)
            except KeyError:
                continue
            band_data = toa.data_vars.get(band_name)
            if band_data is None and callable(band_loader):
                try:
                    band_data = band_loader(band_name)
                except KeyError:
                    band_data = None
            if band_data is None:
                continue
            coeffs = self.rt_model.compute_coefficients(geometry, atmo_state, band_spec, False)
            coeffs = _resample_coefficients_to_template(coeffs, band_data)
            boa = coeffs.apply_correction(band_data)
            band_valid = np.isfinite(boa) & (boa > -0.05) & (boa < 1.5)
            boa_vars[band_name] = boa.where(band_valid)
            band_invalid = ~band_valid
            invalid_boa_mask = (
                band_invalid
                if invalid_boa_mask is None
                else (invalid_boa_mask | band_invalid)
            )

        if not boa_vars:
            raise ValueError(
                "No bands in TOA dataset matched sensor_config bands. "
                f"TOA vars: {list(toa.data_vars)}, "
                f"sensor bands: {self.sensor_config.band_names}"
            )

        boa_ds = xr.Dataset(boa_vars)
        first = list(boa_vars.keys())[0]
        if invalid_boa_mask is None:
            invalid_boa_mask = xr.zeros_like(boa_ds[first], dtype=bool)

        # cloud_mask contract: True = cloudy/invalid, preserved from M1
        if cloud_mask is not None:
            final_cloud_mask = _resample_mask_to_template(cloud_mask, invalid_boa_mask) | invalid_boa_mask.astype(bool)
        else:
            final_cloud_mask = invalid_boa_mask.astype(bool)
        elapsed = time.monotonic() - t0
        return CorrectionResult(
            boa=boa_ds,
            boa_unc=None,
            aot=atmo_state.aot,
            tcwv=atmo_state.tcwv,
            cloud_mask=final_cloud_mask,
            diagnostics=CorrectionDiagnostics(processing_time_s=elapsed),
        )
