"""
Reprojection and resampling utilities using rioxarray.

This module provides functions for coordinate transformations, reprojection,
and resampling of geospatial raster data. It replaces the GDAL-based
reproject.py from the original SIAC.

Example:
    >>> from siac.geo.reprojection import reproject_match, resample
    >>>
    >>> # Reproject to match another raster
    >>> aligned = reproject_match(source, target)
    >>>
    >>> # Resample to different resolution
    >>> resampled = resample(da, target_resolution=30.0)
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING, Any, Literal

import numpy as np
import rioxarray  # noqa: F401
import xarray as xr
from pyproj import CRS, Transformer
from rasterio.enums import Resampling

from siac.geo._crs_compat import (
    crs_equivalent as _crs_equivalent,
)
from siac.geo._crs_compat import (
    default_resampling_for_dtype as _default_resampling_for_dtype,
)

if TYPE_CHECKING:
    from rasterio.transform import Affine

logger = logging.getLogger(__name__)

# Mapping of string names to Resampling enums
RESAMPLING_METHODS = {
    "nearest": Resampling.nearest,
    "bilinear": Resampling.bilinear,
    "cubic": Resampling.cubic,
    "cubic_spline": Resampling.cubic_spline,
    "lanczos": Resampling.lanczos,
    "average": Resampling.average,
    "mode": Resampling.mode,
    "gauss": Resampling.gauss,
    "max": Resampling.max,
    "min": Resampling.min,
    "med": Resampling.med,
    "q1": Resampling.q1,
    "q3": Resampling.q3,
    "sum": Resampling.sum,
}


def _resolve_resampling(
    resampling: str | Resampling | None,
    *,
    dtype: Any = None,
) -> Resampling:
    """Resolve a user-supplied ``resampling`` argument to a :class:`Resampling`.

    * ``None`` → dtype-aware default (nearest for integer/bool, bilinear
      otherwise; REVIEW.md §3.7 reprojection.py:111-143, 170);
    * a :class:`Resampling` instance → returned as-is;
    * a known string → looked up in ``RESAMPLING_METHODS``;
    * an unknown string → ``ValueError`` instead of silently falling back to
      bilinear (REVIEW.md §3.7 reprojection.py:88).
    """
    if resampling is None:
        return _default_resampling_for_dtype(dtype)
    if isinstance(resampling, Resampling):
        return resampling
    if isinstance(resampling, str):
        method = RESAMPLING_METHODS.get(resampling)
        if method is None:
            raise ValueError(
                f"Unknown resampling method: {resampling!r}; "
                f"valid: {sorted(RESAMPLING_METHODS)}"
            )
        return method
    raise TypeError(
        f"resampling must be a str, Resampling enum, or None; got {type(resampling).__name__}"
    )


# =============================================================================
# Reprojection Functions
# =============================================================================


def reproject_to_crs(
    data: xr.DataArray,
    target_crs: str | CRS,
    resolution: float | tuple[float, float] | None = None,
    resampling: str | Resampling | None = None,
    nodata: float | None = None,
) -> xr.DataArray:
    """
    Reproject a DataArray to a target CRS.

    Args:
        data: Input DataArray with CRS information
        target_crs: Target CRS (EPSG code string or pyproj CRS)
        resolution: Target resolution. If tuple, (x_res, y_res).
                    If None, preserves approximate resolution.
        resampling: Resampling method (string or Resampling enum). When
                    ``None`` (the default), the kernel is selected from the
                    input dtype: nearest for integer/boolean rasters
                    (e.g. classification masks), bilinear for floating-point
                    rasters (e.g. reflectance).
        nodata: NoData value for output. If None, uses input nodata.

    Returns:
        Reprojected DataArray

    Raises:
        ValueError: If ``resampling`` is a string that is not in
            :data:`RESAMPLING_METHODS`.

    Example:
        >>> # Reproject to UTM zone 32N
        >>> reprojected = reproject_to_crs(da, "EPSG:32632")
        >>>
        >>> # Reproject with specific resolution
        >>> reprojected = reproject_to_crs(da, "EPSG:4326", resolution=0.0001)
    """
    resampling = _resolve_resampling(resampling, dtype=data.dtype)

    # Handle nodata
    if nodata is not None:
        data = data.rio.write_nodata(nodata)

    # Build reproject kwargs
    reproject_kwargs: dict[str, Any] = {
        "resampling": resampling,
    }

    if resolution is not None:
        if isinstance(resolution, (int, float)):
            reproject_kwargs["resolution"] = (resolution, resolution)
        else:
            reproject_kwargs["resolution"] = resolution

    # Reproject
    result = data.rio.reproject(target_crs, **reproject_kwargs)

    return result


def reproject_match(
    source: xr.DataArray,
    target: xr.DataArray,
    resampling: str | Resampling | None = None,
    nodata: float | None = None,
) -> xr.DataArray:
    """
    Reproject source DataArray to match the grid of target DataArray.

    This ensures pixel-perfect alignment between two rasters.

    Args:
        source: Source DataArray to reproject
        target: Target DataArray defining the output grid
        resampling: Resampling method. When ``None`` (the default), the kernel
            is selected from the source dtype: nearest for integer/boolean
            rasters (e.g. classification masks), bilinear for floating-point
            rasters.
        nodata: NoData value for output

    Returns:
        Source reprojected to match target's CRS, extent, and resolution

    Raises:
        ValueError: If ``resampling`` is a string not in
            :data:`RESAMPLING_METHODS`.

    Example:
        >>> # Align BRDF data to satellite imagery grid
        >>> aligned_brdf = reproject_match(brdf_500m, sentinel2_10m)
    """
    resampling = _resolve_resampling(resampling, dtype=source.dtype)

    if nodata is not None:
        source = source.rio.write_nodata(nodata)

    result = source.rio.reproject_match(target, resampling=resampling)

    return result


def reproject_dataset_match(
    source: xr.Dataset,
    target: xr.DataArray,
    resampling: str | Resampling | None = None,
) -> xr.Dataset:
    """
    Reproject all variables in a Dataset to match a target grid.

    Args:
        source: Source Dataset to reproject
        target: Target DataArray defining the output grid
        resampling: Resampling method. When ``None`` (the default), each
            variable's kernel is chosen from its dtype: nearest for
            integer/boolean (mask, class ID) variables, bilinear for floats.
            Pass an explicit string or :class:`Resampling` to use one method
            for the whole dataset.

    Returns:
        Dataset with all variables reprojected

    Raises:
        ValueError: If ``resampling`` is a string not in
            :data:`RESAMPLING_METHODS`.
    """
    # Resolve once for the whole dataset only when an explicit choice is
    # given; ``None`` means "decide per variable".
    explicit_resampling: Resampling | None
    if resampling is None:
        explicit_resampling = None
    else:
        explicit_resampling = _resolve_resampling(resampling)

    reprojected_vars = {}

    for var_name in source.data_vars:
        var = source[var_name]
        if hasattr(var, "rio") and var.rio.crs is not None:
            method = (
                explicit_resampling
                if explicit_resampling is not None
                else _default_resampling_for_dtype(var.dtype)
            )
            reprojected_vars[var_name] = var.rio.reproject_match(target, resampling=method)
        else:
            # Keep non-spatial variables as-is
            reprojected_vars[var_name] = var

    result = xr.Dataset(reprojected_vars)

    # Copy CRS from target
    if target.rio.crs is not None:
        result = result.rio.write_crs(target.rio.crs)

    return result


# =============================================================================
# Resampling Functions
# =============================================================================


def resample(
    data: xr.DataArray,
    target_resolution: float | tuple[float, float],
    resampling: str | Resampling | None = None,
) -> xr.DataArray:
    """
    Resample a DataArray to a different resolution.

    Args:
        data: Input DataArray
        target_resolution: Target resolution in CRS units (usually meters).
                          If tuple, (x_res, y_res).
        resampling: Resampling method. When ``None`` (the default), the kernel
            is selected from the input dtype: nearest for integer/boolean
            rasters, bilinear for floating-point rasters.

    Returns:
        Resampled DataArray at target resolution

    Raises:
        ValueError: If ``resampling`` is a string not in
            :data:`RESAMPLING_METHODS`.

    Example:
        >>> # Resample 10m data to 500m for MODIS matching
        >>> coarse = resample(sentinel2, target_resolution=500.0)
    """
    resampling = _resolve_resampling(resampling, dtype=data.dtype)

    if isinstance(target_resolution, (int, float)):
        target_resolution = (target_resolution, target_resolution)

    # Calculate new dimensions. ``int(...)`` truncates toward zero and was
    # losing the fractional border — use round() so the resampled raster
    # covers the original extent (REVIEW.md §3.7 reprojection.py:217-228).
    current_res = data.rio.resolution()
    scale_x = abs(current_res[0]) / target_resolution[0]
    scale_y = abs(current_res[1]) / target_resolution[1]

    new_width = max(1, int(round(data.rio.width * scale_x)))
    new_height = max(1, int(round(data.rio.height * scale_y)))

    result = data.rio.reproject(
        data.rio.crs,
        shape=(new_height, new_width),
        resampling=resampling,
    )

    return result


def resample_to_shape(
    data: xr.DataArray,
    target_shape: tuple[int, int],
    resampling: str | Resampling | None = None,
) -> xr.DataArray:
    """
    Resample a DataArray to a specific shape.

    Args:
        data: Input DataArray
        target_shape: Target (height, width)
        resampling: Resampling method. When ``None`` (the default), the kernel
            is selected from the input dtype: nearest for integer/boolean
            rasters, bilinear for floating-point rasters.

    Returns:
        Resampled DataArray with target shape

    Raises:
        ValueError: If ``resampling`` is a string not in
            :data:`RESAMPLING_METHODS`.
    """
    resampling = _resolve_resampling(resampling, dtype=data.dtype)

    result = data.rio.reproject(
        data.rio.crs,
        shape=target_shape,
        resampling=resampling,
    )

    return result


# =============================================================================
# Clipping Functions
# =============================================================================


def clip_to_bounds(
    data: xr.DataArray,
    bounds: tuple[float, float, float, float],
    bounds_crs: str | CRS | None = None,
) -> xr.DataArray:
    """
    Clip a DataArray to a bounding box.

    Args:
        data: Input DataArray
        bounds: Bounding box (xmin, ymin, xmax, ymax)
        bounds_crs: CRS of the bounds. If None, assumes same as data.

    Returns:
        Clipped DataArray

    Example:
        >>> # Clip to specific extent
        >>> clipped = clip_to_bounds(da, (500000, 4000000, 600000, 4100000))
    """
    return data.rio.clip_box(*bounds, crs=bounds_crs)


def clip_to_geometry(
    data: xr.DataArray,
    geometry: dict | list,
    geometry_crs: str | CRS | None = None,
    all_touched: bool = False,
    drop: bool = True,
) -> xr.DataArray:
    """
    Clip a DataArray to a geometry (polygon).

    Args:
        data: Input DataArray
        geometry: GeoJSON-like geometry dict or list of geometries
        geometry_crs: CRS of the geometry. If None, assumes same as data.
        all_touched: If True, all pixels touched by geometry are included
        drop: If True, drops rows/cols that are all nodata after clipping

    Returns:
        Clipped DataArray
    """
    geometries = geometry if isinstance(geometry, list) else [geometry]

    return data.rio.clip(geometries, crs=geometry_crs, all_touched=all_touched, drop=drop)


# =============================================================================
# Coordinate Utilities
# =============================================================================


def transform_bounds(
    bounds: tuple[float, float, float, float],
    source_crs: str | CRS,
    target_crs: str | CRS,
) -> tuple[float, float, float, float]:
    """
    Transform bounding box coordinates between CRS.

    Args:
        bounds: Source bounds (xmin, ymin, xmax, ymax)
        source_crs: Source CRS
        target_crs: Target CRS

    Returns:
        Transformed bounds in target CRS
    """
    from rasterio.warp import transform_bounds as rio_transform_bounds

    xmin, ymin, xmax, ymax = rio_transform_bounds(source_crs, target_crs, *bounds)
    return (float(xmin), float(ymin), float(xmax), float(ymax))


def transform_points(
    x: np.ndarray | list[float],
    y: np.ndarray | list[float],
    source_crs: str | CRS,
    target_crs: str | CRS,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Transform point coordinates between CRS.

    Args:
        x: X coordinates
        y: Y coordinates
        source_crs: Source CRS
        target_crs: Target CRS

    Returns:
        Tuple of (transformed_x, transformed_y)
    """
    transformer = Transformer.from_crs(source_crs, target_crs, always_xy=True)
    x_out, y_out = transformer.transform(x, y)

    return np.array(x_out), np.array(y_out)


def get_resolution(data: xr.DataArray) -> tuple[float, float]:
    """
    Get the resolution of a DataArray.

    Args:
        data: DataArray with spatial coordinates

    Returns:
        Tuple of (x_resolution, y_resolution) as positive values
    """
    res = data.rio.resolution()
    return (abs(res[0]), abs(res[1]))


def get_bounds(data: xr.DataArray) -> tuple[float, float, float, float]:
    """
    Get the bounding box of a DataArray.

    Args:
        data: DataArray with spatial coordinates

    Returns:
        Bounds as (xmin, ymin, xmax, ymax)
    """
    bounds = data.rio.bounds()
    return (float(bounds[0]), float(bounds[1]), float(bounds[2]), float(bounds[3]))


def get_transform(data: xr.DataArray) -> Affine:
    """
    Get the affine transform of a DataArray.

    Args:
        data: DataArray with spatial coordinates

    Returns:
        Affine transform
    """
    return data.rio.transform()


def get_crs(data: xr.DataArray) -> str | None:
    """
    Get the CRS of a DataArray as string.

    Args:
        data: DataArray with spatial coordinates

    Returns:
        CRS as string (e.g., "EPSG:32632") or None
    """
    crs = data.rio.crs
    if crs is None:
        return None
    return str(crs)


# =============================================================================
# Grid Alignment Utilities
# =============================================================================


def align_grids(
    *arrays: xr.DataArray,
    reference_idx: int = 0,
    resampling: str | Resampling | None = None,
) -> list[xr.DataArray]:
    """
    Align multiple DataArrays to the same grid.

    Args:
        *arrays: DataArrays to align
        reference_idx: Index of the reference array (others align to this)
        resampling: Resampling method. When ``None`` (the default), each
            array's kernel is chosen from its dtype (nearest for
            integer/boolean rasters, bilinear for floats).

    Returns:
        List of aligned DataArrays
    """
    if len(arrays) < 2:
        return list(arrays)

    reference = arrays[reference_idx]
    aligned = []

    for i, arr in enumerate(arrays):
        if i == reference_idx:
            aligned.append(arr)
        else:
            aligned.append(reproject_match(arr, reference, resampling=resampling))

    return aligned


def compute_common_bounds(
    *arrays: xr.DataArray,
    method: Literal["intersection", "union"] = "intersection",
) -> tuple[float, float, float, float]:
    """
    Compute common bounds from multiple DataArrays.

    Args:
        *arrays: DataArrays with spatial coordinates
        method: 'intersection' for overlap area, 'union' for combined extent

    Returns:
        Common bounds as (xmin, ymin, xmax, ymax)

    Raises:
        ValueError: If arrays don't overlap (for intersection method)
    """
    if len(arrays) == 0:
        raise ValueError("At least one array required")

    # Ensure all arrays are in same CRS (use first as reference)
    ref_crs = arrays[0].rio.crs

    all_bounds = []
    for arr in arrays:
        bounds = get_bounds(arr)
        # Compare CRS by authority/WKT semantics rather than by Python
        # object identity / string match (REVIEW.md §3.7 reprojection.py:484).
        if not _crs_equivalent(arr.rio.crs, ref_crs):
            bounds = transform_bounds(bounds, arr.rio.crs, ref_crs)
        all_bounds.append(bounds)

    if method == "intersection":
        xmin = max(b[0] for b in all_bounds)
        ymin = max(b[1] for b in all_bounds)
        xmax = min(b[2] for b in all_bounds)
        ymax = min(b[3] for b in all_bounds)

        if xmin >= xmax or ymin >= ymax:
            raise ValueError("Arrays do not overlap")

    else:  # union
        xmin = min(b[0] for b in all_bounds)
        ymin = min(b[1] for b in all_bounds)
        xmax = max(b[2] for b in all_bounds)
        ymax = max(b[3] for b in all_bounds)

    return (xmin, ymin, xmax, ymax)
