"""
Raster writing utilities using rioxarray.

This module provides functions for writing xarray DataArrays and Datasets
to various geospatial raster formats including GeoTIFF, Cloud-Optimized
GeoTIFF (COG), Zarr, and NetCDF.

Example:
    >>> from siac.io.writers import write_raster, write_cog
    >>>
    >>> # Write DataArray to GeoTIFF
    >>> write_raster(da, "/path/to/output.tif")
    >>>
    >>> # Write as Cloud-Optimized GeoTIFF
    >>> write_cog(da, "/path/to/output.tif", overviews=[2, 4, 8])
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Literal, Sequence

import numpy as np
import rioxarray  # noqa: F401
import xarray as xr

logger = logging.getLogger(__name__)

# Compression settings for different algorithms
COMPRESSION_SETTINGS = {
    "deflate": {"compress": "deflate", "zlevel": 6},
    "lzw": {"compress": "lzw"},
    "zstd": {"compress": "zstd", "zstd_level": 9},
    "none": {},
}

# COG-specific settings
COG_SETTINGS = {
    "driver": "COG",
    "blocksize": 512,
    "overview_resampling": "average",
}


# =============================================================================
# Basic Writers
# =============================================================================


def write_raster(
    data: xr.DataArray,
    path: str | Path,
    compression: Literal["deflate", "lzw", "zstd", "none"] = "deflate",
    dtype: str | None = None,
    nodata: float | int | None = None,
    tiled: bool = True,
    blockxsize: int = 512,
    blockysize: int = 512,
    **kwargs: Any,
) -> Path:
    """
    Write an xarray DataArray to a GeoTIFF file.

    Args:
        data: DataArray with spatial coordinates and CRS
        path: Output file path
        compression: Compression algorithm
        dtype: Output data type (e.g., "float32", "uint16"). If None, uses input dtype.
        nodata: NoData value. If None, uses existing or NaN for floats.
        tiled: Whether to write as tiled TIFF
        blockxsize: Tile width
        blockysize: Tile height
        **kwargs: Additional arguments passed to rioxarray.to_raster()

    Returns:
        Path to written file

    Example:
        >>> write_raster(da, "output.tif", compression="lzw", dtype="float32")
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    # Prepare data
    data = _prepare_for_write(data, dtype, nodata)

    # Build write options
    write_kwargs = {
        "driver": "GTiff",
        "tiled": tiled,
        **COMPRESSION_SETTINGS.get(compression, {}),
        **kwargs,
    }

    if tiled:
        write_kwargs["blockxsize"] = blockxsize
        write_kwargs["blockysize"] = blockysize

    # Write using rioxarray
    data.rio.to_raster(str(path), **write_kwargs)

    logger.info(f"Wrote raster to {path}")
    return path


def write_cog(
    data: xr.DataArray,
    path: str | Path,
    compression: Literal["deflate", "lzw", "zstd", "none"] = "deflate",
    dtype: str | None = None,
    nodata: float | int | None = None,
    overviews: Sequence[int] | None = None,
    overview_resampling: str = "average",
    blocksize: int = 512,
    **kwargs: Any,
) -> Path:
    """
    Write an xarray DataArray to a Cloud-Optimized GeoTIFF (COG).

    COGs are optimized for cloud storage and HTTP range requests,
    making them ideal for web-based visualization and analysis.

    Args:
        data: DataArray with spatial coordinates and CRS
        path: Output file path
        compression: Compression algorithm
        dtype: Output data type
        nodata: NoData value
        overviews: Overview levels (e.g., [2, 4, 8, 16]). If None, auto-computed.
        overview_resampling: Resampling method for overviews
        blocksize: Block size for tiling
        **kwargs: Additional arguments passed to rioxarray

    Returns:
        Path to written file
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    # Prepare data
    data = _prepare_for_write(data, dtype, nodata)

    # Auto-compute overview levels if not specified
    if overviews is None:
        overviews = _compute_overview_levels(data)

    # Build COG write options
    write_kwargs = {
        "driver": "COG",
        "blocksize": blocksize,
        "overview_resampling": overview_resampling,
        **COMPRESSION_SETTINGS.get(compression, {}),
        **kwargs,
    }

    if overviews:
        write_kwargs["overviews"] = list(overviews)

    # Write using rioxarray
    data.rio.to_raster(str(path), **write_kwargs)

    logger.info(f"Wrote COG to {path}")
    return path


def write_dataset(
    ds: xr.Dataset,
    output_dir: str | Path,
    prefix: str = "",
    compression: Literal["deflate", "lzw", "zstd", "none"] = "deflate",
    dtype: str | None = None,
    as_cog: bool = True,
    **kwargs: Any,
) -> dict[str, Path]:
    """
    Write an xarray Dataset to multiple raster files.

    Each variable in the Dataset is written as a separate file.

    Args:
        ds: Dataset with spatial coordinates and CRS
        output_dir: Output directory
        prefix: Prefix for output filenames
        compression: Compression algorithm
        dtype: Output data type
        as_cog: Whether to write as COG (True) or regular GeoTIFF (False)
        **kwargs: Additional arguments passed to write functions

    Returns:
        Dictionary mapping variable names to output paths

    Example:
        >>> paths = write_dataset(ds, "/output/", prefix="S2A_BOA_")
        >>> print(paths)  # {"B02": Path("/output/S2A_BOA_B02.tif"), ...}
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    write_fn = write_cog if as_cog else write_raster
    output_paths = {}

    for var_name in ds.data_vars:
        da = ds[var_name]

        # Skip non-spatial variables
        if not hasattr(da, "rio") or da.rio.crs is None:
            logger.warning(f"Skipping {var_name}: no spatial info")
            continue

        filename = f"{prefix}{var_name}.tif"
        output_path = output_dir / filename

        write_fn(da, output_path, compression=compression, dtype=dtype, **kwargs)
        output_paths[var_name] = output_path

    return output_paths


# =============================================================================
# Specialized Writers
# =============================================================================


def write_zarr(
    data: xr.DataArray | xr.Dataset,
    path: str | Path,
    chunks: dict[str, int] | None = None,
    mode: Literal["w", "a"] = "w",
    **kwargs: Any,
) -> Path:
    """
    Write data to a Zarr store.

    Zarr is ideal for large datasets and supports efficient
    chunked, compressed storage.

    Args:
        data: DataArray or Dataset to write
        path: Output Zarr store path
        chunks: Chunk sizes for Zarr. If None, uses existing chunks.
        mode: Write mode ('w' for overwrite, 'a' for append)
        **kwargs: Additional arguments passed to to_zarr()

    Returns:
        Path to Zarr store
    """
    path = Path(path)

    if chunks is not None:
        data = data.chunk(chunks)

    data.to_zarr(str(path), mode=mode, **kwargs)

    logger.info(f"Wrote Zarr store to {path}")
    return path


def write_netcdf(
    data: xr.DataArray | xr.Dataset,
    path: str | Path,
    compression: dict[str, Any] | None = None,
    **kwargs: Any,
) -> Path:
    """
    Write data to NetCDF-CF format.

    Args:
        data: DataArray or Dataset to write
        path: Output file path
        compression: Compression settings per variable
        **kwargs: Additional arguments passed to to_netcdf()

    Returns:
        Path to NetCDF file
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    # Default compression
    if compression is None:
        compression = {"zlib": True, "complevel": 4}

    # Apply encoding to all variables
    if isinstance(data, xr.Dataset):
        encoding = {var: compression for var in data.data_vars}
    else:
        encoding = {data.name or "data": compression}

    data.to_netcdf(str(path), encoding=encoding, **kwargs)

    logger.info(f"Wrote NetCDF to {path}")
    return path


# =============================================================================
# SIAC-Specific Writers
# =============================================================================


def write_boa_products(
    boa: xr.Dataset,
    output_dir: str | Path,
    sensor: str,
    tile_id: str,
    datetime_str: str,
    compression: Literal["deflate", "lzw", "zstd", "none"] = "deflate",
    include_uncertainty: bool = True,
) -> dict[str, Path]:
    """
    Write SIAC BOA (bottom-of-atmosphere) products.

    Creates standardized output filenames following SIAC conventions.

    Args:
        boa: Dataset with BOA reflectance bands
        output_dir: Output directory
        sensor: Sensor name (e.g., "S2A", "L8")
        tile_id: Tile identifier
        datetime_str: Datetime string for filename
        compression: Compression algorithm
        include_uncertainty: Whether to write uncertainty bands

    Returns:
        Dictionary of output paths
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    output_paths = {}
    prefix = f"{sensor}_{tile_id}_{datetime_str}"

    for var_name in boa.data_vars:
        # Skip uncertainty if not requested
        if "_unc" in var_name and not include_uncertainty:
            continue

        da = boa[var_name]
        filename = f"{prefix}_{var_name}_BOA.tif"
        output_path = output_dir / filename

        write_cog(da, output_path, compression=compression, dtype="float32")
        output_paths[var_name] = output_path

    return output_paths


def write_auxiliary_products(
    aot: xr.DataArray,
    tcwv: xr.DataArray,
    output_dir: str | Path,
    sensor: str,
    tile_id: str,
    datetime_str: str,
    compression: Literal["deflate", "lzw", "zstd", "none"] = "deflate",
) -> dict[str, Path]:
    """
    Write SIAC auxiliary products (AOT, TCWV).

    Args:
        aot: Aerosol optical thickness
        tcwv: Total column water vapor
        output_dir: Output directory
        sensor: Sensor name
        tile_id: Tile identifier
        datetime_str: Datetime string
        compression: Compression algorithm

    Returns:
        Dictionary with "aot" and "tcwv" paths
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    prefix = f"{sensor}_{tile_id}_{datetime_str}"

    aot_path = output_dir / f"{prefix}_AOT.tif"
    tcwv_path = output_dir / f"{prefix}_TCWV.tif"

    write_cog(aot, aot_path, compression=compression, dtype="float32")
    write_cog(tcwv, tcwv_path, compression=compression, dtype="float32")

    return {"aot": aot_path, "tcwv": tcwv_path}


def write_rgb_quicklook(
    boa: xr.Dataset,
    output_path: str | Path,
    red_band: str = "B04",
    green_band: str = "B03",
    blue_band: str = "B02",
    scale: tuple[float, float] = (0.0, 0.3),
    target_resolution: float | None = None,
) -> Path:
    """
    Write RGB quicklook image from BOA data.

    Args:
        boa: Dataset with BOA reflectance
        output_path: Output file path
        red_band: Name of red band
        green_band: Name of green band
        blue_band: Name of blue band
        scale: (min, max) reflectance values for scaling to 0-255
        target_resolution: Target resolution for quicklook (downsample if larger)

    Returns:
        Path to quicklook image
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Stack RGB bands
    rgb = xr.concat(
        [boa[red_band], boa[green_band], boa[blue_band]],
        dim="band",
    )
    rgb = rgb.assign_coords(band=["red", "green", "blue"])

    # Resample if needed
    if target_resolution is not None:
        current_res = abs(rgb.rio.resolution()[0])
        if current_res < target_resolution:
            scale_factor = current_res / target_resolution
            new_width = int(rgb.rio.width * scale_factor)
            new_height = int(rgb.rio.height * scale_factor)
            rgb = rgb.rio.reproject(
                rgb.rio.crs,
                shape=(new_height, new_width),
            )

    # Scale to 0-255
    rgb_scaled = (rgb - scale[0]) / (scale[1] - scale[0])
    rgb_scaled = rgb_scaled.clip(0, 1) * 255
    rgb_scaled = rgb_scaled.astype(np.uint8)

    # Write
    rgb_scaled.rio.to_raster(str(output_path), driver="GTiff")

    logger.info(f"Wrote RGB quicklook to {output_path}")
    return output_path


# =============================================================================
# Helper Functions
# =============================================================================


def _prepare_for_write(
    data: xr.DataArray,
    dtype: str | None,
    nodata: float | int | None,
) -> xr.DataArray:
    """Prepare DataArray for writing."""
    # Convert dtype if specified
    if dtype is not None:
        # Handle scaling for integer types
        if dtype in ("uint16", "int16", "uint8"):
            if nodata is None:
                nodata = 0 if dtype.startswith("u") else -9999

        data = data.astype(dtype)

    # Set nodata
    if nodata is not None:
        data = data.rio.write_nodata(nodata)
    elif data.rio.nodata is None and np.issubdtype(data.dtype, np.floating):
        # Default to NaN for floats
        data = data.rio.write_nodata(np.nan)

    return data


def _compute_overview_levels(data: xr.DataArray, max_overview_size: int = 256) -> list[int]:
    """Compute appropriate overview levels based on image size."""
    height, width = data.rio.height, data.rio.width
    min_dim = min(height, width)

    levels = []
    factor = 2
    while min_dim // factor > max_overview_size:
        levels.append(factor)
        factor *= 2

    return levels
