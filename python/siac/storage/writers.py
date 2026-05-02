"""
Raster writing utilities using rioxarray.

This module provides functions for writing xarray DataArrays and Datasets
to various geospatial raster formats including GeoTIFF, Cloud-Optimized
GeoTIFF (COG), Zarr, and NetCDF.

Example:
    >>> from siac.storage.writers import write_raster, write_cog
    >>>
    >>> # Write DataArray to GeoTIFF
    >>> write_raster(da, "/path/to/output.tif")
    >>>
    >>> # Write as Cloud-Optimized GeoTIFF
    >>> write_cog(da, "/path/to/output.tif")
"""

from __future__ import annotations

import importlib.util
import logging
import stat
from pathlib import Path
from typing import TYPE_CHECKING, Any, Literal, TypeAlias, cast

import numpy as np
import rioxarray  # noqa: F401
import xarray as xr
from numpy import typing as npt

__all__ = [
    "write_aot_scatter_plot",
    "write_auxiliary_products",
    "write_boa_products",
    "write_cloud_mask_preview",
    "write_cog",
    "write_dataset",
    "write_false_colour_preview",
    "write_field_preview",
    "write_netcdf",
    "write_raster",
    "write_rgb_quicklook",
    "write_zarr",
]

logger = logging.getLogger(__name__)

if TYPE_CHECKING:
    from siac.runtime import AOTScatterBandDiagnostics

CompressionSettings: TypeAlias = dict[str, str | int]
UInt8Array: TypeAlias = npt.NDArray[np.uint8]

# Compression settings for GTiff/GeoTIFF drivers.
GTIFF_COMPRESSION_SETTINGS: dict[str, CompressionSettings] = {
    "deflate": {"compress": "deflate", "zlevel": 6},
    "lzw": {"compress": "lzw"},
    "zstd": {"compress": "zstd", "zstd_level": 9},
    "none": {},
}

# Compression settings accepted by GDAL's COG driver.
COG_COMPRESSION_SETTINGS: dict[str, CompressionSettings] = {
    "deflate": {"compress": "deflate", "level": 6},
    "lzw": {"compress": "lzw"},
    "zstd": {"compress": "zstd", "level": 9},
    "none": {},
}

# COG-specific settings
COG_SETTINGS = {
    "driver": "COG",
    "blocksize": 512,
    "overview_resampling": "average",
}

_NETCDF_DEFAULT = object()


def _repair_directory_mode(directory: Path) -> None:
    current_mode = stat.S_IMODE(directory.stat().st_mode)
    required_owner_bits = stat.S_IRUSR | stat.S_IWUSR | stat.S_IXUSR
    updated_mode = current_mode | required_owner_bits
    if current_mode & stat.S_IRGRP:
        updated_mode |= stat.S_IXGRP
    if current_mode & stat.S_IROTH:
        updated_mode |= stat.S_IXOTH
    if updated_mode != current_mode:
        directory.chmod(updated_mode)


def ensure_writable_directory(path: str | Path) -> Path:
    """Create a directory tree and make sure each level is owner-traversable."""
    directory = Path(path).expanduser()
    if directory.exists():
        if not directory.is_dir():
            raise NotADirectoryError(directory)
        _repair_directory_mode(directory)
        return directory

    missing: list[Path] = []
    current = directory
    while not current.exists():
        missing.append(current)
        parent = current.parent
        if parent == current:
            break
        current = parent

    if not current.is_dir():
        raise NotADirectoryError(current)
    if current != Path(current.anchor) and current != Path():
        _repair_directory_mode(current)

    for current in reversed(missing):
        try:
            current.mkdir()
        except FileExistsError as exc:
            if not current.is_dir():
                raise NotADirectoryError(current) from exc
        _repair_directory_mode(current)
    return directory


def _netcdf_fill_value(data: xr.DataArray) -> Any:
    """Choose a NetCDF fill-value policy compatible with variable dtype."""
    dtype = np.dtype(data.dtype)
    if np.issubdtype(dtype, np.integer) or np.issubdtype(dtype, np.bool_):
        return None
    if "_FillValue" in data.encoding:
        return data.encoding["_FillValue"]
    return _NETCDF_DEFAULT


def _sanitize_netcdf_array(data: xr.DataArray) -> xr.DataArray:
    """Remove fill-value attrs that are invalid for integer/bool NetCDF variables."""
    out = data.copy(deep=False)
    dtype = np.dtype(out.dtype)
    if np.issubdtype(dtype, np.integer) or np.issubdtype(dtype, np.bool_):
        for key in ("_FillValue", "missing_value", "fill_value"):
            out.attrs.pop(key, None)
    return out


def _netcdf_base_encoding(data: xr.DataArray) -> dict[str, Any]:
    """Preserve NetCDF-relevant metadata from an xarray variable encoding."""
    encoding: dict[str, Any] = {}
    grid_mapping = data.encoding.get("grid_mapping")
    if grid_mapping is None and "spatial_ref" in data.coords:
        grid_mapping = "spatial_ref"
    if grid_mapping is not None:
        encoding["grid_mapping"] = grid_mapping
    return encoding


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
    ensure_writable_directory(path.parent)

    # Prepare data
    data = _prepare_for_write(data, dtype, nodata)

    # Build write options
    compression_options = GTIFF_COMPRESSION_SETTINGS.get(compression, {})
    write_kwargs: dict[str, Any] = {
        "driver": "GTiff",
        "tiled": tiled,
        **compression_options,
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
        overview_resampling: Resampling method for overviews
        blocksize: Block size for tiling
        **kwargs: Additional arguments passed to rioxarray

    Returns:
        Path to written file
    """
    path = Path(path)
    ensure_writable_directory(path.parent)

    # Prepare data
    data = _prepare_for_write(data, dtype, nodata)

    # Build COG write options
    compression_options = COG_COMPRESSION_SETTINGS.get(compression, {})
    write_kwargs: dict[str, Any] = {
        "driver": "COG",
        "blocksize": blocksize,
        "overview_resampling": overview_resampling,
        **compression_options,
        **kwargs,
    }

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
    output_dir = ensure_writable_directory(output_dir)

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
    ensure_writable_directory(path.parent)

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
    ensure_writable_directory(path.parent)
    payload: xr.DataArray | xr.Dataset

    # Default compression
    if compression is None:
        compression = {"zlib": True, "complevel": 4}
    effective_compression = dict(compression)

    engine = kwargs.get("engine")
    if engine is None:
        if importlib.util.find_spec("h5netcdf") is not None:
            kwargs["engine"] = "h5netcdf"
            engine = "h5netcdf"
        elif importlib.util.find_spec("netCDF4") is not None:
            kwargs["engine"] = "netcdf4"
            engine = "netcdf4"
        else:
            raise RuntimeError(
                "NetCDF output requires h5netcdf or netCDF4; scipy fallback is disabled. "
                "Install h5netcdf or netCDF4."
            )
    elif str(engine).lower() == "scipy":
        raise ValueError(
            "NetCDF scipy backend is not supported; use engine='h5netcdf' or engine='netcdf4'."
        )

    if isinstance(data, xr.Dataset):
        payload = data.copy(deep=False)
        for name in list(payload.data_vars):
            payload[name] = _sanitize_netcdf_array(payload[name])
        for name in list(payload.coords):
            payload.coords[name] = _sanitize_netcdf_array(payload.coords[name])
    else:
        payload = _sanitize_netcdf_array(data)

    # Apply encoding to all variables and coordinate variables.
    encoding: dict[str, dict[str, Any]] = {}
    if isinstance(payload, xr.Dataset):
        for name, var in payload.data_vars.items():
            var_encoding = _netcdf_base_encoding(var)
            var_encoding.update(effective_compression)
            fill_value = _netcdf_fill_value(var)
            if fill_value is not _NETCDF_DEFAULT:
                var_encoding["_FillValue"] = fill_value
            encoding[name] = var_encoding
        for name, _coord in payload.coords.items():
            encoding[name] = {"_FillValue": None}
    else:
        var_name = payload.name or "data"
        var_encoding = _netcdf_base_encoding(payload)
        var_encoding.update(effective_compression)
        fill_value = _netcdf_fill_value(payload)
        if fill_value is not _NETCDF_DEFAULT:
            var_encoding["_FillValue"] = fill_value
        encoding[var_name] = var_encoding
        for name, _coord in payload.coords.items():
            encoding[name] = {"_FillValue": None}

    payload.to_netcdf(str(path), encoding=encoding, **kwargs)

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
    output_dir = ensure_writable_directory(output_dir)

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
    output_dir = ensure_writable_directory(output_dir)

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
    ensure_writable_directory(output_path.parent)

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
    # Reproject/scaling paths can carry float NaN nodata metadata.
    # Replace NaN pixels and use an integer nodata before casting to uint8.
    rgb_scaled = rgb_scaled.fillna(0).astype(np.uint8)
    rgb_scaled = rgb_scaled.rio.write_nodata(0)

    # Write
    rgb_scaled.rio.to_raster(str(output_path), driver="GTiff")

    logger.info(f"Wrote RGB quicklook to {output_path}")
    return output_path


def write_aot_scatter_plot(
    scatter: AOTScatterBandDiagnostics,
    output_path: str | Path,
    *,
    width: int = 720,
    height: int = 720,
) -> Path:
    """Write a PNG scatter plot for one aerosol-solver band."""
    from PIL import Image, ImageDraw, ImageFont

    output_path = Path(output_path)
    ensure_writable_directory(output_path.parent)

    image = Image.new("RGBA", (width, height), (255, 255, 255, 255))
    draw = ImageDraw.Draw(image, "RGBA")
    font = ImageFont.load_default()

    margin_left = 76
    margin_right = 28
    margin_top = 56
    margin_bottom = 80
    plot_left = margin_left
    plot_top = margin_top
    plot_right = width - margin_right
    plot_bottom = height - margin_bottom

    draw.rectangle(
        (plot_left, plot_top, plot_right, plot_bottom), outline=(60, 60, 60, 255), width=2
    )

    x_values = np.asarray(scatter.surface_reflectance, dtype=np.float64)
    observed_values = np.asarray(scatter.observed_toa, dtype=np.float64)
    simulated_values = np.asarray(scatter.simulated_toa, dtype=np.float64)

    if x_values.size:
        x_min = float(np.nanmin(x_values))
        x_max = float(np.nanmax(x_values))
        y_min = float(np.nanmin(np.concatenate([observed_values, simulated_values])))
        y_max = float(np.nanmax(np.concatenate([observed_values, simulated_values])))
    else:
        x_min = 0.0
        x_max = 1.0
        y_min = 0.0
        y_max = 1.0

    if not np.isfinite(x_min) or not np.isfinite(x_max) or x_min == x_max:
        x_min, x_max = 0.0, 1.0
    if not np.isfinite(y_min) or not np.isfinite(y_max) or y_min == y_max:
        y_min, y_max = 0.0, 1.0

    x_pad = max((x_max - x_min) * 0.05, 1.0e-3)
    y_pad = max((y_max - y_min) * 0.05, 1.0e-3)
    x_min -= x_pad
    x_max += x_pad
    y_min -= y_pad
    y_max += y_pad

    def _sx(value: float) -> int:
        fraction = (value - x_min) / (x_max - x_min)
        return int(round(plot_left + fraction * (plot_right - plot_left)))

    def _sy(value: float) -> int:
        fraction = (value - y_min) / (y_max - y_min)
        return int(round(plot_bottom - fraction * (plot_bottom - plot_top)))

    for step in range(6):
        x_tick = x_min + (x_max - x_min) * step / 5.0
        y_tick = y_min + (y_max - y_min) * step / 5.0
        x_pos = _sx(x_tick)
        y_pos = _sy(y_tick)
        draw.line((x_pos, plot_bottom, x_pos, plot_bottom + 6), fill=(60, 60, 60, 255), width=1)
        draw.line((plot_left - 6, y_pos, plot_left, y_pos), fill=(60, 60, 60, 255), width=1)
        draw.text(
            (x_pos - 18, plot_bottom + 10), f"{x_tick:.2f}", fill=(40, 40, 40, 255), font=font
        )
        draw.text((8, y_pos - 6), f"{y_tick:.2f}", fill=(40, 40, 40, 255), font=font)

    observed_color = (70, 70, 70, 110)
    simulated_palette = {
        "B02": (20, 90, 220, 120),
        "B04": (210, 55, 40, 120),
    }
    simulated_color = simulated_palette.get(scatter.band_name, (0, 140, 120, 120))

    def _draw_points(
        x_series: np.ndarray, y_series: np.ndarray, color: tuple[int, int, int, int]
    ) -> None:
        for x_value, y_value in zip(x_series, y_series, strict=False):
            if not np.isfinite(x_value) or not np.isfinite(y_value):
                continue
            px = _sx(float(x_value))
            py = _sy(float(y_value))
            draw.ellipse((px - 1, py - 1, px + 1, py + 1), fill=color, outline=None)

    _draw_points(x_values, observed_values, observed_color)
    _draw_points(x_values, simulated_values, simulated_color)

    draw.text(
        (plot_left, 16),
        f"AOT Solver Scatter {scatter.band_name}  sampled={x_values.size} total={scatter.total_valid_count}",
        fill=(20, 20, 20, 255),
        font=font,
    )
    draw.text(
        (plot_left, height - 28), "Simulated surface reflectance", fill=(20, 20, 20, 255), font=font
    )
    draw.text((16, 24), "TOA", fill=(20, 20, 20, 255), font=font)

    legend_x = plot_right - 146
    legend_y = plot_top + 12
    draw.rectangle((legend_x, legend_y, legend_x + 12, legend_y + 12), fill=observed_color)
    draw.text((legend_x + 18, legend_y), "Observed TOA", fill=(20, 20, 20, 255), font=font)
    draw.rectangle((legend_x, legend_y + 20, legend_x + 12, legend_y + 32), fill=simulated_color)
    draw.text((legend_x + 18, legend_y + 20), "Simulated TOA", fill=(20, 20, 20, 255), font=font)

    image.save(output_path, format="PNG")
    logger.info("Wrote aerosol scatter plot to %s", output_path)
    return output_path


# =============================================================================
# Preview PNG Writers
# =============================================================================


def _field_to_uint8(
    data: np.ndarray,
    vmin: float,
    vmax: float,
    *,
    mask: np.ndarray | None = None,
) -> UInt8Array:
    """Scale a 2-D float field to 0-255 uint8, masking NaN/invalid to 0."""
    arr = np.asarray(data, dtype=np.float64)
    span = max(vmax - vmin, 1.0e-9)
    scaled = np.clip((arr - vmin) / span, 0.0, 1.0) * 255.0
    scaled = np.where(np.isfinite(scaled), scaled, 0.0)
    if mask is not None:
        scaled = np.where(mask, 0.0, scaled)
    return cast("UInt8Array", np.asarray(scaled, dtype=np.uint8))


def _apply_colourmap(
    data: np.ndarray,
    vmin: float,
    vmax: float,
    *,
    palette: str = "viridis",
    mask: np.ndarray | None = None,
) -> UInt8Array:
    """Map a 2-D float field to an (H, W, 3) uint8 RGB array via a colourmap.

    Supported palettes: ``viridis`` (blue-green-yellow), ``magma`` (black-magenta-yellow),
    ``coolwarm`` (blue-white-red).
    """
    idx = _field_to_uint8(data, vmin, vmax, mask=mask)

    # Build 256-entry lookup tables (no matplotlib dependency).
    lut = _LUT_CACHE.get(palette)
    if lut is None:
        lut = _build_lut(palette)
        _LUT_CACHE[palette] = lut

    rgb = np.asarray(lut[idx], dtype=np.uint8)
    # Mask → dark grey
    if mask is not None:
        rgb[mask] = 40
    return cast("UInt8Array", np.asarray(rgb, dtype=np.uint8))


_LUT_CACHE: dict[str, UInt8Array] = {}


def _build_lut(palette: str) -> UInt8Array:
    """Build a 256 x 3 uint8 lookup table for *palette*."""
    t = np.linspace(0.0, 1.0, 256)
    if palette == "viridis":
        # Simplified viridis:  dark-blue → teal → green → yellow
        r = np.clip(np.where(t < 0.5, 0.26 + 0.1 * t, -0.4 + 2.4 * t), 0, 1)
        g = np.clip(np.where(t < 0.3, 0.0 + 1.5 * t, 0.35 + 0.25 * t + 0.5 * t**2), 0, 1)
        b = np.clip(np.where(t < 0.6, 0.33 + 0.6 * t, 1.5 - 1.5 * t), 0, 1)
    elif palette == "magma":
        r = np.clip(t**0.6, 0, 1)
        g = np.clip(0.15 * t + 0.5 * t**3, 0, 1)
        b = np.clip(0.3 + 0.5 * t - 0.7 * t**2 + 0.5 * t**3, 0, 1)
    elif palette == "coolwarm":
        r = np.clip(np.where(t < 0.5, 0.2 + 1.6 * t, 1.0), 0, 1)
        g = np.clip(np.where(t < 0.5, 0.2 + 1.6 * t, 1.0 - 1.6 * (t - 0.5)), 0, 1)
        b = np.clip(np.where(t < 0.5, 1.0, 1.0 - 1.6 * (t - 0.5)), 0, 1)
    else:
        # Grey fallback
        r = g = b = t
    lut = np.stack([r, g, b], axis=-1)
    return cast("UInt8Array", np.asarray(lut * 255, dtype=np.uint8))


def write_false_colour_preview(
    boa: xr.Dataset,
    output_path: str | Path,
    *,
    red_band: str = "B08",
    green_band: str = "B04",
    blue_band: str = "B03",
    scale: tuple[float, float] = (0.0, 0.4),
) -> Path | None:
    """Write a false-colour composite PNG (default: NIR-Red-Green)."""
    from PIL import Image

    output_path = Path(output_path)
    bands_available = set(boa.data_vars)
    if not {red_band, green_band, blue_band} <= bands_available:
        # Try B8A if B08 not available
        if red_band == "B08" and "B8A" in bands_available:
            red_band = "B8A"
        else:
            return None

    ensure_writable_directory(output_path.parent)
    r = _field_to_uint8(boa[red_band].values, scale[0], scale[1])
    g = _field_to_uint8(boa[green_band].values, scale[0], scale[1])
    b = _field_to_uint8(boa[blue_band].values, scale[0], scale[1])
    rgb = np.stack([r, g, b], axis=-1)

    Image.fromarray(rgb, "RGB").save(output_path, format="PNG")
    logger.info("Wrote false-colour preview to %s", output_path)
    return output_path


def write_field_preview(
    field: xr.DataArray,
    output_path: str | Path,
    *,
    vmin: float | None = None,
    vmax: float | None = None,
    palette: str = "viridis",
    title: str = "",
    unit: str = "",
    cloud_mask: xr.DataArray | None = None,
) -> Path:
    """Write a colour-mapped 2-D field preview PNG with an embedded colour bar."""
    from PIL import Image, ImageDraw, ImageFont

    output_path = Path(output_path)
    ensure_writable_directory(output_path.parent)

    values = np.asarray(field.values, dtype=np.float64)
    mask = None
    if cloud_mask is not None:
        mask_vals = np.asarray(cloud_mask.values, dtype=bool)
        if mask_vals.shape == values.shape:
            mask = mask_vals

    finite = values[np.isfinite(values)]
    if mask is not None:
        finite = values[np.isfinite(values) & ~mask]
    if finite.size == 0:
        finite = np.array([0.0, 1.0])

    if vmin is None:
        vmin = float(np.percentile(finite, 2))
    if vmax is None:
        vmax = float(np.percentile(finite, 98))
    if vmin == vmax:
        vmax = vmin + 1.0

    rgb = _apply_colourmap(values, vmin, vmax, palette=palette, mask=mask)

    # Add colour bar strip (30px) and label area (30px) at the bottom
    h, w = rgb.shape[:2]
    bar_h, label_h = 20, 30
    canvas = np.full((h + bar_h + label_h, w, 3), 255, dtype=np.uint8)
    canvas[:h, :, :] = rgb

    # Colour bar gradient
    bar_t = np.linspace(0.0, 1.0, w)
    bar_idx = (bar_t * 255).astype(np.uint8)
    lut = _LUT_CACHE.get(palette) or _build_lut(palette)
    bar_row = lut[bar_idx]  # (w, 3)
    for row in range(h, h + bar_h):
        canvas[row, :, :] = bar_row

    image = Image.fromarray(canvas, "RGB")
    draw = ImageDraw.Draw(image)
    font = ImageFont.load_default()

    # Title at top
    if title:
        draw.text((4, 2), title, fill=(255, 255, 255), font=font)

    # Min / max labels under colour bar
    label_y = h + bar_h + 4
    draw.text((4, label_y), f"{vmin:.3f}", fill=(40, 40, 40), font=font)
    max_label = f"{vmax:.3f}"
    if unit:
        max_label = f"{vmax:.3f} {unit}"
    draw.text(
        (max(0, w - 8 * len(max_label) - 4), label_y), max_label, fill=(40, 40, 40), font=font
    )

    image.save(output_path, format="PNG")
    logger.info("Wrote field preview to %s", output_path)
    return output_path


def write_cloud_mask_preview(
    boa: xr.Dataset,
    cloud_mask: xr.DataArray,
    output_path: str | Path,
    *,
    red_band: str = "B04",
    green_band: str = "B03",
    blue_band: str = "B02",
    scale: tuple[float, float] = (0.0, 0.3),
    cloud_colour: tuple[int, int, int] = (255, 40, 40),
    cloud_alpha: float = 0.55,
) -> Path | None:
    """Write an RGB image with the cloud mask overlaid in a semi-transparent colour."""
    from PIL import Image

    output_path = Path(output_path)
    if not {red_band, green_band, blue_band} <= set(boa.data_vars):
        return None

    ensure_writable_directory(output_path.parent)
    r = _field_to_uint8(boa[red_band].values, scale[0], scale[1])
    g = _field_to_uint8(boa[green_band].values, scale[0], scale[1])
    b = _field_to_uint8(boa[blue_band].values, scale[0], scale[1])

    mask = np.asarray(cloud_mask.values, dtype=bool)
    # Resize mask to match BOA if needed
    if mask.shape != r.shape:
        from scipy.ndimage import zoom as _zoom

        zoom_y = r.shape[0] / max(mask.shape[0], 1)
        zoom_x = r.shape[1] / max(mask.shape[1], 1)
        mask = _zoom(mask.astype(np.float32), (zoom_y, zoom_x), order=0) > 0.5

    # Blend cloud colour into masked pixels
    a = cloud_alpha
    r = np.where(mask, (a * cloud_colour[0] + (1 - a) * r).astype(np.uint8), r)
    g = np.where(mask, (a * cloud_colour[1] + (1 - a) * g).astype(np.uint8), g)
    b = np.where(mask, (a * cloud_colour[2] + (1 - a) * b).astype(np.uint8), b)

    rgb = np.stack([r, g, b], axis=-1)
    Image.fromarray(rgb, "RGB").save(output_path, format="PNG")
    logger.info("Wrote cloud mask preview to %s", output_path)
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
        if dtype in ("uint16", "int16", "uint8") and nodata is None:
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
