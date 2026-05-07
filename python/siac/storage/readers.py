"""
Raster reading utilities using rioxarray.

This module provides functions for reading geospatial raster data into
xarray DataArrays and Datasets, with support for various formats including
GeoTIFF, COG, VRT, and remote URLs.

Example:
    >>> from siac.storage.readers import read_raster, read_multiband
    >>>
    >>> # Read single band
    >>> da = read_raster("/path/to/file.tif")
    >>>
    >>> # Read with window
    >>> da = read_raster_window("/path/to/file.tif", bounds=(0, 0, 100, 100))
    >>>
    >>> # Read multiple bands as dataset
    >>> ds = read_multiband(["/path/to/b1.tif", "/path/to/b2.tif"], ["B01", "B02"])
"""

from __future__ import annotations

import logging
import re
from pathlib import Path
from typing import TYPE_CHECKING, Any, cast

import numpy as np
import rioxarray  # noqa: F401 - needed for .rio accessor
import xarray as xr
from rasterio.enums import Resampling

from siac.geo._crs_compat import crs_equivalent as _crs_equivalent

if TYPE_CHECKING:
    from collections.abc import MutableMapping, Sequence

# RFC 3986 scheme prefix for fsspec/GDAL-style remote URLs.
# Matches "http://", "https://", "s3://", "gs://", "azure://", "abfs://",
# "file://", etc. — anything that looks like ``<scheme>://...`` with a
# lowercase-leading scheme. Local paths (``/tmp/foo`` or ``C:\foo``) do not
# match, which is the behaviour we want.
_REMOTE_SCHEME_RE = re.compile(r"^[a-z][a-z0-9+.\-]*://", re.IGNORECASE)


def _is_remote_uri(path: str) -> bool:
    """Return True if ``path`` looks like a URI that fsspec can open."""
    return bool(_REMOTE_SCHEME_RE.match(path))

logger = logging.getLogger(__name__)


def _get_remote_zarr_mapper(path: str) -> MutableMapping[str, bytes]:
    """Open a typed fsspec mapping for remote Zarr stores."""
    import fsspec  # type: ignore[import-untyped]

    return cast("MutableMapping[str, bytes]", fsspec.get_mapper(path))


# =============================================================================
# Single Raster Reading
# =============================================================================


def read_raster(
    path: str | Path,
    band: int | None = None,
    chunks: dict[str, int] | None = None,
    masked: bool = True,
    overview_level: int | None = None,
) -> xr.DataArray:
    """
    Read a raster file into an xarray DataArray.

    Uses rioxarray for reading, which supports GeoTIFF, COG, VRT,
    and remote URLs (http://, /vsicurl/, s3://).

    Args:
        path: Path to raster file or URL
        band: Specific band to read (1-indexed). If None, reads all bands.
        chunks: Chunk sizes for dask arrays, e.g., {"x": 2048, "y": 2048}.
                If None, loads entire array into memory.
        masked: Whether to mask nodata values as NaN
        overview_level: Overview level to read (0 = highest resolution).
                        If None, reads full resolution.

    Returns:
        xr.DataArray with spatial coordinates and CRS information.
        The array has dimensions (band, y, x) or (y, x) if single band.

    Example:
        >>> da = read_raster("/path/to/image.tif")
        >>> print(da.rio.crs)
        >>> print(da.rio.bounds())
    """
    path = str(path)

    # Build rioxarray open options
    open_kwargs: dict[str, Any] = {
        "masked": masked,
    }

    if chunks is not None:
        open_kwargs["chunks"] = chunks

    if overview_level is not None:
        open_kwargs["overview_level"] = overview_level

    # Open with rioxarray. The xr.open_dataarray("...", engine="rasterio")
    # path keeps the underlying rasterio dataset alive on the returned object,
    # which leaks file handles in long-running pipelines (REVIEW.md §2.8).
    # Eagerly load values into memory and close the dataset before returning.
    # All current callers consume the array via .values / .astype paths, so
    # this matches existing behaviour without changing the contract.
    with xr.open_dataarray(path, engine="rasterio", **open_kwargs) as opened:
        da = opened.load()

    # Select specific band if requested
    if band is not None and "band" in da.dims:
        da = da.sel(band=band)

    # Squeeze single-band dimension to return (y, x) instead of (1, y, x)
    if "band" in da.dims and da.sizes["band"] == 1:
        da = da.squeeze("band", drop=True)

    # Ensure we have CRS info
    if da.rio.crs is None:
        logger.warning(f"No CRS found in {path}")

    return da


def read_raster_window(
    path: str | Path,
    bounds: tuple[float, float, float, float],
    bounds_crs: str | None = None,
    chunks: dict[str, int] | None = None,
    masked: bool = True,
) -> xr.DataArray:
    """
    Read a subset of a raster defined by bounding box.

    Args:
        path: Path to raster file or URL
        bounds: Bounding box (xmin, ymin, xmax, ymax)
        bounds_crs: CRS of the bounds. If None, assumes same as raster.
        chunks: Chunk sizes for dask arrays
        masked: Whether to mask nodata values

    Returns:
        xr.DataArray clipped to the specified bounds
    """
    del chunks  # Windowed reads are eager for now.

    import rasterio
    from rasterio.transform import xy
    from rasterio.warp import transform_bounds
    from rasterio.windows import Window, from_bounds

    path = str(path)
    with rasterio.open(path) as src:
        target_bounds = bounds
        # Compare by authority/WKT semantics, not strings
        # (REVIEW.md §2.8, §3.7 readers.py).
        if bounds_crs is not None and not _crs_equivalent(src.crs, bounds_crs):
            target_bounds = transform_bounds(bounds_crs, src.crs, *bounds)

        window = from_bounds(*target_bounds, transform=src.transform)
        window = window.intersection(Window(0, 0, src.width, src.height))
        window = window.round_offsets().round_lengths()
        if int(window.width) <= 0 or int(window.height) <= 0:
            raise ValueError(f"Bounds {bounds!r} do not intersect raster {path}")

        values = src.read(window=window, masked=masked)
        transform = src.window_transform(window)
        height = int(window.height)
        width = int(window.width)

        x = np.asarray(
            xy(
                transform,
                np.zeros(width, dtype=np.int32),
                np.arange(width, dtype=np.int32),
                offset="center",
            )[0],
            dtype=np.float64,
        )
        y = np.asarray(
            xy(
                transform,
                np.arange(height, dtype=np.int32),
                np.zeros(height, dtype=np.int32),
                offset="center",
            )[1],
            dtype=np.float64,
        )

        if np.ma.isMaskedArray(values):
            values = np.asarray(
                np.ma.filled(np.ma.asarray(values, dtype=np.float32), np.nan), dtype=np.float32
            )

        if values.ndim == 3 and values.shape[0] == 1:
            values = values[0]

        if values.ndim == 3:
            da = xr.DataArray(
                np.asarray(values),
                dims=("band", "y", "x"),
                coords={"band": np.asarray(src.indexes, dtype=np.int32), "y": y, "x": x},
            )
        else:
            da = xr.DataArray(
                np.asarray(values),
                dims=("y", "x"),
                coords={"y": y, "x": x},
            )
        da = da.rio.write_crs(src.crs)
        return da.rio.write_transform(transform)


def read_raster_at_resolution(
    path: str | Path,
    target_resolution: float,
    resampling: Resampling = Resampling.bilinear,
    chunks: dict[str, int] | None = None,
) -> xr.DataArray:
    """
    Read a raster and resample to target resolution.

    Args:
        path: Path to raster file
        target_resolution: Target resolution in CRS units (usually meters)
        resampling: Resampling method (default: bilinear)
        chunks: Chunk sizes for dask arrays

    Returns:
        xr.DataArray resampled to target resolution
    """
    da = read_raster(path, chunks=chunks)

    # Get current resolution
    current_res = abs(da.rio.resolution()[0])

    if abs(current_res - target_resolution) < 0.01:
        # Already at target resolution
        return da

    # Calculate scale factor
    scale = current_res / target_resolution

    # Resample
    new_width = int(da.rio.width * scale)
    new_height = int(da.rio.height * scale)

    da = da.rio.reproject(
        da.rio.crs,
        shape=(new_height, new_width),
        resampling=resampling,
    )

    return da


# =============================================================================
# Multi-band / Multi-file Reading
# =============================================================================


def read_multiband(
    paths: Sequence[str | Path],
    band_names: Sequence[str] | None = None,
    chunks: dict[str, int] | None = None,
    masked: bool = True,
) -> xr.Dataset:
    """
    Read multiple raster files as a single Dataset.

    Each file becomes a variable in the Dataset, useful for reading
    satellite bands stored as separate files.

    Args:
        paths: List of paths to raster files
        band_names: Names for each band. If None, uses filenames.
        chunks: Chunk sizes for dask arrays
        masked: Whether to mask nodata values

    Returns:
        xr.Dataset with each band as a variable

    Example:
        >>> paths = ["B02.tif", "B03.tif", "B04.tif"]
        >>> ds = read_multiband(paths, band_names=["blue", "green", "red"])
    """
    if band_names is None:
        band_names = [Path(p).stem for p in paths]

    if len(paths) != len(band_names):
        raise ValueError(
            f"Number of paths ({len(paths)}) must match band_names ({len(band_names)})"
        )

    data_vars = {}
    ref_crs = None

    for path, name in zip(paths, band_names):
        da = read_raster(path, chunks=chunks, masked=masked)

        # Squeeze out band dimension if present and single
        if "band" in da.dims and da.sizes["band"] == 1:
            da = da.squeeze("band", drop=True)

        # Store CRS info from first band
        if ref_crs is None:
            ref_crs = da.rio.crs
            da.rio.transform()

        data_vars[name] = da

    ds = xr.Dataset(data_vars)

    # Write CRS to dataset
    if ref_crs is not None:
        ds = ds.rio.write_crs(ref_crs)

    return ds


def read_multiband_stack(
    paths: Sequence[str | Path],
    band_names: Sequence[str] | None = None,
    chunks: dict[str, int] | None = None,
    masked: bool = True,
) -> xr.DataArray:
    """
    Read multiple raster files and stack into a single DataArray.

    Unlike read_multiband which creates a Dataset, this creates a single
    DataArray with a 'band' dimension.

    Args:
        paths: List of paths to raster files
        band_names: Names for each band (used as coordinate values)
        chunks: Chunk sizes for dask arrays
        masked: Whether to mask nodata values

    Returns:
        xr.DataArray with dimensions (band, y, x)
    """
    if band_names is None:
        band_names = [Path(p).stem for p in paths]

    arrays = []
    for path in paths:
        da = read_raster(path, chunks=chunks, masked=masked)
        if "band" in da.dims and da.sizes["band"] == 1:
            da = da.squeeze("band", drop=True)
        arrays.append(da)

    # Stack along new band dimension
    stacked = xr.concat(arrays, dim="band")
    stacked = stacked.assign_coords(band=list(band_names))

    return stacked


# =============================================================================
# Specialized Readers
# =============================================================================


def read_jp2(
    path: str | Path,
    chunks: dict[str, int] | None = None,
) -> xr.DataArray:
    """
    Read JPEG2000 file (commonly used for Sentinel-2).

    Args:
        path: Path to JP2 file
        chunks: Chunk sizes for dask arrays

    Returns:
        xr.DataArray with image data
    """
    return read_raster(path, chunks=chunks, masked=True)


def read_hdf_subdataset(
    path: str | Path,
    subdataset: str,
    chunks: dict[str, int] | None = None,
) -> xr.DataArray:
    """
    Read a subdataset from an HDF4/HDF5 file.

    Commonly used for MODIS products.

    Args:
        path: Path to HDF file
        subdataset: Name of subdataset to read
        chunks: Chunk sizes for dask arrays

    Returns:
        xr.DataArray with subdataset data

    Example:
        >>> da = read_hdf_subdataset(
        ...     "MCD43A1.hdf",
        ...     "BRDF_Albedo_Parameters_Band1_iso"
        ... )
    """
    # Construct GDAL-style subdataset path
    hdf_path = f'HDF4_EOS:EOS_GRID:"{path}":{subdataset}'

    return read_raster(hdf_path, chunks=chunks)


def read_netcdf_variable(
    path: str | Path,
    variable: str,
    chunks: dict[str, int] | None = None,
) -> xr.DataArray:
    """
    Read a variable from a NetCDF file with spatial coordinates.

    Args:
        path: Path to NetCDF file
        variable: Name of variable to read
        chunks: Chunk sizes for dask arrays

    Returns:
        xr.DataArray with variable data
    """
    # The dataset must outlive the returned variable; previously the dataset
    # handle was leaked because we returned ds[variable] without closing
    # ``ds`` (REVIEW.md §2.8). Eagerly load the variable into memory, capture
    # any CRS metadata, then close the dataset on the way out.
    with xr.open_dataset(path, chunks=chunks) as ds:
        da = ds[variable].load()

        crs_wkt: str | None = None
        if "crs" in ds.attrs:
            crs_wkt = ds.attrs["crs"]
        elif "spatial_ref" in ds:
            crs_wkt = ds["spatial_ref"].attrs.get("crs_wkt")

    if crs_wkt is not None:
        da = da.rio.write_crs(crs_wkt)

    return da


def read_zarr_array(
    path: str | Path,
    variable: str | None = None,
    chunks: dict[str, int] | None = None,
) -> xr.DataArray | xr.Dataset:
    """
    Read data from a Zarr store.

    Args:
        path: Path or URL to Zarr store
        variable: Specific variable to read. If None, returns Dataset.
        chunks: Chunk sizes (if None, uses Zarr's native chunks)

    Returns:
        xr.DataArray if variable specified, else xr.Dataset
    """
    path = str(path)

    # Handle remote URLs. The previous explicit list missed `gs://`,
    # `azure://`, `abfs://`, `file://`, etc. (REVIEW.md §3.7 readers.py:443).
    if _is_remote_uri(path):
        mapper = _get_remote_zarr_mapper(path)
        ds = xr.open_zarr(mapper, chunks=chunks)
    else:
        ds = xr.open_zarr(path, chunks=chunks)

    if variable is not None:
        return ds[variable]

    return ds


# =============================================================================
# Utility Functions
# =============================================================================


def get_raster_info(path: str | Path) -> dict[str, Any]:
    """
    Get metadata information about a raster file without reading data.

    Args:
        path: Path to raster file

    Returns:
        Dictionary with raster metadata including:
        - crs: Coordinate reference system
        - bounds: Bounding box (xmin, ymin, xmax, ymax)
        - resolution: Pixel resolution (x, y)
        - shape: Array shape (bands, height, width)
        - dtype: Data type
        - nodata: NoData value
    """
    import rasterio

    path = str(path)

    with rasterio.open(path) as src:
        info = {
            "crs": str(src.crs) if src.crs else None,
            "bounds": src.bounds,
            "resolution": src.res,
            "shape": (src.count, src.height, src.width),
            "dtype": str(src.dtypes[0]),
            "nodata": src.nodata,
            "transform": src.transform,
            "driver": src.driver,
        }

    return info


def check_rasters_aligned(
    path1: str | Path,
    path2: str | Path,
    tolerance: float = 1e-6,
) -> bool:
    """
    Check if two rasters are spatially aligned (same grid).

    Args:
        path1: Path to first raster
        path2: Path to second raster
        tolerance: Tolerance for floating point comparison

    Returns:
        True if rasters have same CRS, bounds, and resolution
    """
    info1 = get_raster_info(path1)
    info2 = get_raster_info(path2)

    # Check CRS by authority/WKT semantics, not string equality.
    # ``"EPSG:4326"`` and the verbose WKT for the same authority must
    # compare equal (REVIEW.md §2.8, §3.7 readers.py:495-527).
    if not _crs_equivalent(info1["crs"], info2["crs"]):
        return False

    # Check resolution
    if not np.allclose(info1["resolution"], info2["resolution"], atol=tolerance):
        return False

    # Check bounds
    return np.allclose(
        [info1["bounds"].left, info1["bounds"].bottom, info1["bounds"].right, info1["bounds"].top],
        [info2["bounds"].left, info2["bounds"].bottom, info2["bounds"].right, info2["bounds"].top],
        atol=tolerance,
    )
