"""Shared readers and grid helpers for Earthdata BRDF and MAIAC products."""

from __future__ import annotations

import re
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import TYPE_CHECKING, Any

import h5py
import numpy as np
import rioxarray  # noqa: F401
import xarray as xr
from pyhdf.SD import SD, SDC

from siac.io.reprojection import transform_bounds

if TYPE_CHECKING:
    from rasterio.enums import Resampling

MODLAND_SINUSOIDAL_CRS = (
    "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +R=6371007.181 +units=m +no_defs"
)
_MODLAND_TILE_SIZE_M = 1111950.5196666666
_MODLAND_X_MIN_M = -20015109.354
_MODLAND_Y_MAX_M = 10007554.677
_TILE_RE = re.compile(r"\.h(?P<h>\d{2})v(?P<v>\d{2})\.")
_DATE_RE = re.compile(r"\.A(?P<year>\d{4})(?P<doy>\d{3})\.")


@dataclass(frozen=True)
class ProductBandDefinition:
    """One spectral band exposed by a BRDF product."""

    label: str
    wavelength_nm: float
    parameter_dataset: str
    qa_dataset: str | None = None


def decode_attr(value: Any) -> Any:
    """Convert HDF/HDF5 attribute values into plain Python scalars where possible."""
    if isinstance(value, np.ndarray):
        if value.size == 1:
            return decode_attr(value.reshape(-1)[0])
        return [decode_attr(item) for item in value.tolist()]
    if isinstance(value, (np.floating, np.integer)):
        return value.item()
    if isinstance(value, (bytes, np.bytes_)):
        return value.decode("utf-8", errors="ignore")
    return value


def attr_scalar(attrs: dict[str, Any], key: str, default: float | int | None = None) -> float | int | None:
    """Return a numeric scalar attribute value."""
    if key not in attrs:
        return default
    value = decode_attr(attrs[key])
    if isinstance(value, list):
        if not value:
            return default
        value = value[0]
    if value is None:
        return default
    return float(value) if isinstance(value, float) else int(value) if isinstance(value, int) else float(value)


def parse_tile_indices(path: str | Path) -> tuple[int, int]:
    """Extract MODLAND h/v tile indices from a granule filename."""
    match = _TILE_RE.search(Path(path).name)
    if match is None:
        raise ValueError(f"Could not parse tile indices from {Path(path).name!r}")
    return int(match.group("h")), int(match.group("v"))


def parse_granule_date(path: str | Path) -> datetime:
    """Extract acquisition date from a granule filename."""
    match = _DATE_RE.search(Path(path).name)
    if match is None:
        raise ValueError(f"Could not parse acquisition date from {Path(path).name!r}")
    return datetime.strptime(
        f"{match.group('year')}{match.group('doy')}",
        "%Y%j",
    )


def modland_tile_coords(
    h_index: int,
    v_index: int,
    height: int,
    width: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Return 1-D x/y center coordinates for a MODLAND sinusoidal tile."""
    if height <= 0 or width <= 0:
        raise ValueError("height and width must be positive")

    ulx = _MODLAND_X_MIN_M + h_index * _MODLAND_TILE_SIZE_M
    uly = _MODLAND_Y_MAX_M - v_index * _MODLAND_TILE_SIZE_M
    xres = _MODLAND_TILE_SIZE_M / float(width)
    yres = _MODLAND_TILE_SIZE_M / float(height)

    x = ulx + (np.arange(width, dtype=np.float64) + 0.5) * xres
    y = uly - (np.arange(height, dtype=np.float64) + 0.5) * yres
    return x, y


def modland_tile_bounds(h_index: int, v_index: int) -> tuple[float, float, float, float]:
    """Return the full MODLAND tile bounds in the native sinusoidal CRS."""
    ulx = _MODLAND_X_MIN_M + h_index * _MODLAND_TILE_SIZE_M
    uly = _MODLAND_Y_MAX_M - v_index * _MODLAND_TILE_SIZE_M
    return (
        ulx,
        uly - _MODLAND_TILE_SIZE_M,
        ulx + _MODLAND_TILE_SIZE_M,
        uly,
    )


def make_native_grid_dataarray(
    values: np.ndarray,
    *,
    granule_path: str | Path,
    dims: tuple[str, ...] | None = None,
    coords: dict[str, Any] | None = None,
) -> xr.DataArray:
    """Attach MODLAND sinusoidal x/y coordinates and CRS to a native granule array."""
    h_index, v_index = parse_tile_indices(granule_path)
    if values.ndim < 2:
        raise ValueError(f"Expected at least 2 dimensions, got shape={values.shape}")

    height = int(values.shape[-2])
    width = int(values.shape[-1])
    x, y = modland_tile_coords(h_index, v_index, height, width)

    data_dims = dims or tuple(f"dim_{i}" for i in range(values.ndim - 2)) + ("y", "x")
    data_coords = dict(coords or {})
    data_coords.setdefault("x", x)
    data_coords.setdefault("y", y)

    da = xr.DataArray(np.asarray(values), dims=data_dims, coords=data_coords)
    da = da.rio.set_spatial_dims(x_dim="x", y_dim="y")
    return da.rio.write_crs(MODLAND_SINUSOIDAL_CRS)


def build_target_template(
    bounds: tuple[float, float, float, float],
    crs: str,
    resolution: float,
    *,
    fill_value: float = np.nan,
) -> xr.DataArray:
    """Build a target grid template for rioxarray reprojection."""
    xmin, ymin, xmax, ymax = bounds
    if resolution <= 0:
        raise ValueError(f"resolution must be > 0, got {resolution}")

    width = max(1, int(np.ceil((xmax - xmin) / resolution)))
    height = max(1, int(np.ceil((ymax - ymin) / resolution)))
    x = xmin + (np.arange(width, dtype=np.float64) + 0.5) * resolution
    y = ymax - (np.arange(height, dtype=np.float64) + 0.5) * resolution

    target = xr.DataArray(
        np.full((height, width), fill_value, dtype=np.float32),
        dims=("y", "x"),
        coords={"x": x, "y": y},
    )
    target = target.rio.set_spatial_dims(x_dim="x", y_dim="y")
    return target.rio.write_crs(crs)


def clip_native_to_target_bounds(
    data: xr.DataArray,
    *,
    target_bounds: tuple[float, float, float, float],
    target_crs: str,
    pad_pixels: int = 2,
) -> xr.DataArray:
    """Clip a native MODLAND array to the AOI transformed into source CRS."""
    source_bounds = transform_bounds(target_bounds, target_crs, data.rio.crs)
    xres, yres = data.rio.resolution()
    x_pad = abs(float(xres)) * max(0, int(pad_pixels))
    y_pad = abs(float(yres)) * max(0, int(pad_pixels))
    xmin, ymin, xmax, ymax = source_bounds

    mask_x = (data.coords["x"] >= xmin - x_pad) & (data.coords["x"] <= xmax + x_pad)
    mask_y = (data.coords["y"] >= ymin - y_pad) & (data.coords["y"] <= ymax + y_pad)
    return data.where(mask_x & mask_y, drop=True)


def reproject_native_to_target(
    data: xr.DataArray,
    *,
    target_bounds: tuple[float, float, float, float],
    target_crs: str,
    target_resolution: float,
    resampling: Resampling,
    nodata: float | None = None,
) -> xr.DataArray:
    """Clip a native array to the AOI and reproject it to the target grid."""
    clipped = clip_native_to_target_bounds(
        data,
        target_bounds=target_bounds,
        target_crs=target_crs,
    )
    if nodata is not None:
        clipped = clipped.rio.write_nodata(nodata)
    target = build_target_template(target_bounds, target_crs, target_resolution)
    return clipped.rio.reproject_match(target, resampling=resampling)


def read_hdf4_dataset(path: str | Path, dataset_name: str) -> tuple[np.ndarray, dict[str, Any]]:
    """Read an HDF4 SDS plus decoded attributes."""
    sd = SD(str(path), SDC.READ)
    sds = sd.select(dataset_name)
    return np.asarray(sds.get()), {key: decode_attr(value) for key, value in sds.attributes().items()}


def read_hdf5_dataset(path: str | Path, dataset_name: str) -> tuple[np.ndarray, dict[str, Any]]:
    """Read an HDF5 dataset plus decoded attributes."""
    with h5py.File(path, "r") as handle:
        dataset = handle[dataset_name]
        return np.asarray(dataset[...]), {key: decode_attr(value) for key, value in dataset.attrs.items()}


def apply_scale_and_mask(
    values: np.ndarray,
    attrs: dict[str, Any],
) -> np.ndarray:
    """Convert a raw product array to float32, applying fill value and scale/offset."""
    out = np.asarray(values, dtype=np.float32)
    fill_value = attr_scalar(attrs, "_FillValue", None)
    if fill_value is not None:
        out = np.where(out == float(fill_value), np.nan, out)

    valid_range = decode_attr(attrs.get("valid_range"))
    if isinstance(valid_range, list) and len(valid_range) == 2:
        lo, hi = float(valid_range[0]), float(valid_range[1])
        out = np.where((out < lo) | (out > hi), np.nan, out)

    scale_factor = float(attr_scalar(attrs, "scale_factor", 1.0) or 1.0)
    add_offset = float(attr_scalar(attrs, "add_offset", 0.0) or 0.0)
    out = out * scale_factor + add_offset
    return out.astype(np.float32, copy=False)


def reduce_orbit_stack(values: np.ndarray) -> np.ndarray:
    """Reduce a per-orbit stack to one 2-D field via nanmean."""
    arr = np.asarray(values, dtype=np.float32)
    if arr.ndim <= 2:
        return arr.astype(np.float32, copy=False)
    valid = np.isfinite(arr)
    count = np.sum(valid, axis=0, dtype=np.int32)
    total = np.nansum(arr, axis=0, dtype=np.float32)
    out = np.full(arr.shape[1:], np.nan, dtype=np.float32)
    np.divide(total, count, out=out, where=count > 0)
    return out
