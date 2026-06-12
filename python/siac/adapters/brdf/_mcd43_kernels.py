"""Kernel-weight payload assembly for the Earthaccess BRDF providers.

Pure array/xarray logic moved out of ``siac.adapters.brdf.mcd43_earthaccess``:
allocating the temporal/spatial payload arrays, packing the per-band BRDF
parameter + uncertainty layers into a single mergeable stack (and unpacking it
again), and assembling ``BRDFKernelWeights`` from those arrays. The provider
classes bind these functions as staticmethods so the existing method seams
(``provider._pack_payload_stack`` etc.) keep resolving, and re-import the
shared error tuples and type aliases defined here. This module must not
import ``mcd43_earthaccess`` (one-way).
"""

from __future__ import annotations

from typing import TYPE_CHECKING, TypeAlias

import numpy as np
import xarray as xr

from siac.adapters.earthdata_common import ProductBandDefinition
from siac.domain import SensorBand
from siac.runtime import BRDFKernelWeights

if TYPE_CHECKING:
    from collections.abc import Sequence

# HDF4/HDF5 libraries raise their own exception types on I/O failures.
# Import them defensively so they can be included in except clauses.
_HDF4Error: type[BaseException]
try:
    from pyhdf.error import HDF4Error as _HDF4Error  # type: ignore[import-untyped,no-redef]
except ImportError:  # pragma: no cover
    _HDF4Error = OSError
_HDF5Error: type[BaseException]
try:
    from h5py import HDF5ExtError as _HDF5Error  # type: ignore[no-redef]
except (ImportError, AttributeError):  # pragma: no cover
    _HDF5Error = OSError

# Combined tuple used in except clauses throughout this module.
#
# REVIEW.md §2.1, §3.3 mcd43_earthaccess.py:62-69:
# This tuple is intentionally wide because it covers two distinct concerns:
#   (a) genuine I/O failures from HDF4/HDF5 reading (``OSError`` and the
#       library-specific ``_HDF4Error``/``_HDF5Error`` types), and
#   (b) downstream parsing/scaling failures (``KeyError`` for missing
#       dataset names, ``ValueError`` for bad scale-factor metadata,
#       ``RuntimeError`` for GDAL warp/translate failures, ``TypeError``
#       for shape mismatches in ``apply_scale_and_mask``).
#
# Each consumer of this tuple must call ``logger.warning(..., exc_info=True)``
# so that a typo in a dataset name doesn't silently downgrade to "use
# defaults". Use ``_TRANSIENT_DATA_READ_ERRORS`` instead for blocks where
# only category (a) is expected — that lets programming bugs in category (b)
# propagate.
_TRANSIENT_DATA_READ_ERRORS: tuple[type[BaseException], ...] = (
    OSError,
    _HDF4Error,
    _HDF5Error,
)
_DATA_READ_ERRORS: tuple[type[BaseException], ...] = (
    OSError,
    KeyError,
    ValueError,
    TypeError,
    RuntimeError,
    _HDF4Error,
    _HDF5Error,
)

_RequestedBand: TypeAlias = SensorBand
_RequestedBandCoord: TypeAlias = str
_RequestedBandSpec: TypeAlias = tuple[_RequestedBandCoord, ProductBandDefinition]
_PAYLOAD_FIELDS = ("f0", "f1", "f2", "unc")


def allocate_temporal_payload_arrays(
    time_axis: np.ndarray,
    requested: Sequence[_RequestedBandSpec],
    *,
    y_size: int,
    x_size: int,
) -> tuple[np.ndarray, np.ndarray]:
    params_values: np.ndarray = np.full(
        (len(time_axis), len(requested), 3, y_size, x_size),
        np.nan,
        dtype=np.float32,
    )
    unc_values: np.ndarray = np.full(
        (len(time_axis), len(requested), y_size, x_size),
        np.nan,
        dtype=np.float32,
    )
    return params_values, unc_values


def allocate_spatial_payload_arrays(
    requested: Sequence[_RequestedBandSpec],
    target_template: xr.DataArray,
) -> tuple[np.ndarray, np.ndarray]:
    y_size = int(target_template.sizes["y"])
    x_size = int(target_template.sizes["x"])
    params_values: np.ndarray = np.full(
        (len(requested), 3, y_size, x_size),
        np.nan,
        dtype=np.float32,
    )
    unc_values: np.ndarray = np.full(
        (len(requested), y_size, x_size),
        np.nan,
        dtype=np.float32,
    )
    return params_values, unc_values


def temporal_weights_from_arrays(
    params_values: np.ndarray,
    unc_values: np.ndarray,
    *,
    requested: Sequence[_RequestedBandSpec],
    time_axis: np.ndarray,
    target_template: xr.DataArray,
) -> BRDFKernelWeights:
    coords = {
        "time": xr.IndexVariable("time", time_axis),
        "band": xr.IndexVariable("band", [band_coord for band_coord, _band in requested]),
        "y": target_template.coords["y"],
        "x": target_template.coords["x"],
    }

    def _wrap(values: np.ndarray) -> xr.DataArray:
        return target_array_like(
            target_template,
            values,
            dims=("time", "band", "y", "x"),
            coords=coords,
        )

    scaled_unc = unc_values * np.float32(1.1)
    return BRDFKernelWeights(
        f0=_wrap(params_values[:, :, 0, :, :]),
        f1=_wrap(params_values[:, :, 1, :, :]),
        f2=_wrap(params_values[:, :, 2, :, :]),
        f0_unc=_wrap(unc_values),
        f1_unc=_wrap(scaled_unc),
        f2_unc=_wrap(scaled_unc),
        reflectance_unc=_wrap(unc_values),
    )


def native_array_like(
    reference: xr.DataArray,
    values: np.ndarray,
    *,
    dims: tuple[str, ...],
    coords: dict[str, object],
) -> xr.DataArray:
    out = xr.DataArray(np.asarray(values, dtype=np.float32), dims=dims, coords=coords)
    out = out.rio.set_spatial_dims(x_dim="x", y_dim="y")
    reference_crs = reference.rio.crs
    if reference_crs is not None:
        out = out.rio.write_crs(reference_crs)
    try:
        return out.rio.write_transform(reference.rio.transform(recalc=True))
    except _DATA_READ_ERRORS:
        return out


def target_array_like(
    target_template: xr.DataArray,
    values: np.ndarray,
    *,
    dims: tuple[str, ...],
    coords: dict[str, object],
) -> xr.DataArray:
    out = xr.DataArray(np.asarray(values, dtype=np.float32), dims=dims, coords=coords)
    out = out.rio.set_spatial_dims(x_dim="x", y_dim="y")
    out = out.rio.write_crs(target_template.rio.crs)
    return out.rio.write_transform(target_template.rio.transform(recalc=True))


def pack_payload_stack(
    params: xr.DataArray,
    unc: xr.DataArray,
) -> xr.DataArray:
    params = params.transpose("band", "parameter", "y", "x")
    unc = unc.transpose("band", "y", "x")
    band_coords = np.asarray(params.coords["band"].values, dtype=object)
    payload_values = np.concatenate(
        [
            np.asarray(params.values, dtype=np.float32),
            np.asarray(unc.values, dtype=np.float32)[:, np.newaxis, :, :],
        ],
        axis=1,
    )
    layer_count = band_coords.size * len(_PAYLOAD_FIELDS)
    layer_values = payload_values.reshape(layer_count, *payload_values.shape[-2:])
    return native_array_like(
        params,
        layer_values,
        dims=("layer", "y", "x"),
        coords={
            "layer": xr.IndexVariable("layer", np.arange(layer_count, dtype=np.int32)),
            "y": params.coords["y"],
            "x": params.coords["x"],
        },
    )


def unpack_payload_stack(
    payload: xr.DataArray,
    *,
    requested: Sequence[_RequestedBandSpec],
) -> tuple[xr.DataArray, xr.DataArray]:
    payload = payload.transpose("layer", "y", "x")
    band_coords = [band_coord for band_coord, _band in requested]
    expected_layers = len(band_coords) * len(_PAYLOAD_FIELDS)
    if payload.sizes["layer"] != expected_layers:
        raise ValueError(
            f"Expected {expected_layers} payload layers for {len(band_coords)} band(s), "
            f"got {payload.sizes['layer']}"
        )
    values = np.asarray(payload.values, dtype=np.float32).reshape(
        len(band_coords),
        len(_PAYLOAD_FIELDS),
        payload.sizes["y"],
        payload.sizes["x"],
    )
    coords = {
        "band": xr.IndexVariable("band", band_coords),
        "y": payload.coords["y"],
        "x": payload.coords["x"],
    }
    params = native_array_like(
        payload,
        values[:, :3, :, :],
        dims=("band", "parameter", "y", "x"),
        coords={
            **coords,
            "parameter": xr.IndexVariable("parameter", ["f0", "f1", "f2"]),
        },
    )
    unc = native_array_like(
        payload,
        values[:, 3, :, :],
        dims=("band", "y", "x"),
        coords=coords,
    )
    return params, unc


def stack_parameter_cube(
    params: tuple[xr.DataArray, xr.DataArray, xr.DataArray],
) -> xr.DataArray:
    return xr.concat(
        list(params),
        dim=xr.IndexVariable("parameter", ["f0", "f1", "f2"]),
    )


def fill_parameter_defaults(params: xr.DataArray) -> xr.DataArray:
    defaults = xr.DataArray(
        np.array([0.20, 0.05, 0.02], dtype=np.float32),
        dims=["parameter"],
        coords={"parameter": ["f0", "f1", "f2"]},
    )
    return params.fillna(defaults)


def weights_from_layers(
    params: xr.DataArray,
    unc: xr.DataArray,
) -> BRDFKernelWeights:
    unc = unc.transpose("band", "y", "x")
    return BRDFKernelWeights(
        f0=params.sel(parameter="f0", drop=True).transpose("band", "y", "x"),
        f1=params.sel(parameter="f1", drop=True).transpose("band", "y", "x"),
        f2=params.sel(parameter="f2", drop=True).transpose("band", "y", "x"),
        f0_unc=unc,
        f1_unc=(unc * np.float32(1.1)).transpose("band", "y", "x"),
        f2_unc=(unc * np.float32(1.1)).transpose("band", "y", "x"),
        reflectance_unc=unc,
    )
