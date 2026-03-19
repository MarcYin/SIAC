"""
Geometry and vector utilities for SIAC.

This module provides functions for working with geometries, including
loading AOIs, converting between formats, and spatial operations.

Example:
    >>> from siac.geo.geometry import load_aoi, get_raster_footprint
    >>>
    >>> # Load AOI from GeoJSON
    >>> aoi = load_aoi("/path/to/aoi.geojson")
    >>>
    >>> # Get raster footprint
    >>> footprint = get_raster_footprint(da)
"""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Any

import numpy as np
import xarray as xr
from pyproj import CRS, Transformer

logger = logging.getLogger(__name__)


# =============================================================================
# AOI Loading
# =============================================================================


def load_aoi(
    aoi: str | Path | dict | tuple | list,
    target_crs: str | CRS | None = None,
) -> dict:
    """
    Load an Area of Interest from various formats.

    Supported formats:
    - GeoJSON file path
    - WKT string
    - GeoJSON dict
    - Bounding box tuple (xmin, ymin, xmax, ymax)
    - Bounding box list [xmin, ymin, xmax, ymax]

    Args:
        aoi: AOI specification in one of the supported formats
        target_crs: Target CRS to transform the AOI to

    Returns:
        GeoJSON-like dict representing the AOI polygon

    Example:
        >>> # From file
        >>> aoi = load_aoi("/path/to/aoi.geojson")
        >>>
        >>> # From bounds
        >>> aoi = load_aoi((500000, 4000000, 600000, 4100000))
        >>>
        >>> # From WKT
        >>> aoi = load_aoi("POLYGON ((0 0, 1 0, 1 1, 0 1, 0 0))")
    """
    if isinstance(aoi, (str, Path)):
        aoi_str = str(aoi)

        # Check if it's a file path
        if Path(aoi_str).exists():
            return _load_geojson_file(aoi_str, target_crs)

        # Check if it's WKT
        if aoi_str.upper().startswith(("POLYGON", "MULTIPOLYGON", "POINT", "LINESTRING")):
            return _parse_wkt(aoi_str, target_crs)

        raise ValueError(f"Could not parse AOI: {aoi_str}")

    elif isinstance(aoi, dict):
        # Already GeoJSON-like
        return _normalize_geojson(aoi, target_crs)

    elif isinstance(aoi, (tuple, list)) and len(aoi) == 4:
        # Bounding box
        return bounds_to_polygon(aoi, target_crs)

    else:
        raise ValueError(f"Unsupported AOI type: {type(aoi)}")


def _load_geojson_file(path: str, target_crs: str | CRS | None) -> dict:
    """Load GeoJSON from file."""
    with Path(path).open() as f:
        geojson = json.load(f)

    return _normalize_geojson(geojson, target_crs)


def _normalize_geojson(geojson: dict, target_crs: str | CRS | None) -> dict:
    """Normalize GeoJSON and optionally transform CRS."""
    # Handle FeatureCollection
    if geojson.get("type") == "FeatureCollection":
        if len(geojson.get("features", [])) == 0:
            raise ValueError("Empty FeatureCollection")
        # Use first feature's geometry
        geometry = geojson["features"][0]["geometry"]
    elif geojson.get("type") == "Feature":
        geometry = geojson["geometry"]
    else:
        # Assume it's a geometry directly
        geometry = geojson

    # Get source CRS if available
    source_crs = None
    if "crs" in geojson:
        crs_props = geojson["crs"].get("properties", {})
        if "name" in crs_props:
            source_crs = crs_props["name"]

    # Transform if needed
    if target_crs is not None and source_crs is not None:
        geometry = _transform_geometry(geometry, source_crs, target_crs)

    return geometry


def _parse_wkt(wkt: str, _target_crs: str | CRS | None) -> dict:
    """Parse WKT string to GeoJSON geometry."""
    try:
        from shapely import wkt as shapely_wkt
        from shapely.geometry import mapping

        geom = shapely_wkt.loads(wkt)
        return mapping(geom)
    except ImportError:
        # Fallback: simple WKT parsing for polygons
        return _simple_wkt_parse(wkt)


def _simple_wkt_parse(wkt: str) -> dict:
    """Simple WKT parser for POLYGON without shapely."""
    wkt = wkt.strip()

    if not wkt.upper().startswith("POLYGON"):
        raise ValueError(f"Only POLYGON WKT supported without shapely: {wkt[:50]}")

    # Extract coordinates
    # POLYGON ((x1 y1, x2 y2, ...))
    coords_str = wkt[wkt.find("((") + 2 : wkt.rfind("))")]
    coords = []

    for point_str in coords_str.split(","):
        parts = point_str.strip().split()
        if len(parts) >= 2:
            coords.append([float(parts[0]), float(parts[1])])

    return {
        "type": "Polygon",
        "coordinates": [coords],
    }


def _transform_geometry(
    geometry: dict,
    source_crs: str | CRS,
    target_crs: str | CRS,
) -> dict:
    """Transform geometry coordinates between CRS."""
    transformer = Transformer.from_crs(source_crs, target_crs, always_xy=True)

    def transform_coords(coords: list[Any]) -> list[Any]:
        if isinstance(coords[0], (list, tuple)):
            return [transform_coords(c) for c in coords]
        else:
            x, y = transformer.transform(coords[0], coords[1])
            return [x, y]

    new_geometry = geometry.copy()
    new_geometry["coordinates"] = transform_coords(geometry["coordinates"])

    return new_geometry


# =============================================================================
# Bounds and Polygon Utilities
# =============================================================================


def bounds_to_polygon(
    bounds: tuple[float, float, float, float],
    _crs: str | CRS | None = None,
) -> dict:
    """
    Convert bounding box to GeoJSON polygon.

    Args:
        bounds: (xmin, ymin, xmax, ymax)
        crs: CRS of the bounds (stored in properties)

    Returns:
        GeoJSON polygon dict
    """
    xmin, ymin, xmax, ymax = bounds

    polygon = {
        "type": "Polygon",
        "coordinates": [
            [
                [xmin, ymin],
                [xmax, ymin],
                [xmax, ymax],
                [xmin, ymax],
                [xmin, ymin],  # Close the ring
            ]
        ],
    }

    return polygon


def polygon_to_bounds(geometry: dict) -> tuple[float, float, float, float]:
    """
    Extract bounding box from a GeoJSON geometry.

    Args:
        geometry: GeoJSON geometry dict

    Returns:
        Bounds as (xmin, ymin, xmax, ymax)
    """
    coords = _flatten_coordinates(geometry["coordinates"])

    xs = [c[0] for c in coords]
    ys = [c[1] for c in coords]

    return (min(xs), min(ys), max(xs), max(ys))


def _flatten_coordinates(coords: list) -> list:
    """Recursively flatten nested coordinate lists."""
    if isinstance(coords[0], (int, float)):
        return [coords]

    result = []
    for item in coords:
        result.extend(_flatten_coordinates(item))
    return result


# =============================================================================
# Raster Geometry Functions
# =============================================================================


def get_raster_bounds(
    data: xr.DataArray,
    crs: str | CRS | None = None,
) -> tuple[float, float, float, float]:
    """
    Get bounding box of a raster DataArray.

    Args:
        data: DataArray with spatial coordinates
        crs: Target CRS for bounds. If None, uses raster's CRS.

    Returns:
        Bounds as (xmin, ymin, xmax, ymax)
    """
    bounds = data.rio.bounds()

    if crs is not None and data.rio.crs is not None:
        from siac.geo.reprojection import transform_bounds

        bounds = transform_bounds(bounds, data.rio.crs, crs)

    return bounds


def get_raster_footprint(
    data: xr.DataArray,
    crs: str | CRS | None = None,
) -> dict:
    """
    Get the footprint polygon of a raster.

    Args:
        data: DataArray with spatial coordinates
        crs: Target CRS for footprint

    Returns:
        GeoJSON polygon representing raster extent
    """
    bounds = get_raster_bounds(data, crs)
    return bounds_to_polygon(bounds, crs)


def get_raster_center(
    data: xr.DataArray,
    crs: str | CRS | None = None,
) -> tuple[float, float]:
    """
    Get the center point of a raster.

    Args:
        data: DataArray with spatial coordinates
        crs: Target CRS for center point

    Returns:
        Center as (x, y)
    """
    bounds = get_raster_bounds(data, crs)
    x = (bounds[0] + bounds[2]) / 2
    y = (bounds[1] + bounds[3]) / 2
    return (x, y)


def raster_to_geojson_feature(
    data: xr.DataArray,
    properties: dict | None = None,
) -> dict:
    """
    Convert raster extent to GeoJSON Feature.

    Args:
        data: DataArray with spatial coordinates
        properties: Optional properties dict

    Returns:
        GeoJSON Feature dict
    """
    footprint = get_raster_footprint(data, crs="EPSG:4326")

    feature = {
        "type": "Feature",
        "geometry": footprint,
        "properties": properties or {},
    }

    # Add CRS info
    if data.rio.crs is not None:
        feature["properties"]["crs"] = str(data.rio.crs)

    return feature


# =============================================================================
# Spatial Queries
# =============================================================================


def bounds_intersect(
    bounds1: tuple[float, float, float, float],
    bounds2: tuple[float, float, float, float],
) -> bool:
    """
    Check if two bounding boxes intersect (share interior area).

    Adjacent bounds that only touch at an edge are NOT considered intersecting.

    Args:
        bounds1: First bounds (xmin, ymin, xmax, ymax)
        bounds2: Second bounds (xmin, ymin, xmax, ymax)

    Returns:
        True if bounds have overlapping interior area
    """
    return not (
        bounds1[2] <= bounds2[0]  # bounds1 is left of or adjacent to bounds2
        or bounds1[0] >= bounds2[2]  # bounds1 is right of or adjacent to bounds2
        or bounds1[3] <= bounds2[1]  # bounds1 is below or adjacent to bounds2
        or bounds1[1] >= bounds2[3]  # bounds1 is above or adjacent to bounds2
    )


def bounds_contains(
    outer: tuple[float, float, float, float],
    inner: tuple[float, float, float, float],
) -> bool:
    """
    Check if outer bounds completely contain inner bounds.

    Args:
        outer: Outer bounds (xmin, ymin, xmax, ymax)
        inner: Inner bounds (xmin, ymin, xmax, ymax)

    Returns:
        True if outer contains inner
    """
    return (
        outer[0] <= inner[0]
        and outer[1] <= inner[1]
        and outer[2] >= inner[2]
        and outer[3] >= inner[3]
    )


def bounds_intersection(
    bounds1: tuple[float, float, float, float],
    bounds2: tuple[float, float, float, float],
) -> tuple[float, float, float, float] | None:
    """
    Compute intersection of two bounding boxes.

    Args:
        bounds1: First bounds (xmin, ymin, xmax, ymax)
        bounds2: Second bounds (xmin, ymin, xmax, ymax)

    Returns:
        Intersection bounds, or None if no intersection
    """
    if not bounds_intersect(bounds1, bounds2):
        return None

    return (
        max(bounds1[0], bounds2[0]),
        max(bounds1[1], bounds2[1]),
        min(bounds1[2], bounds2[2]),
        min(bounds1[3], bounds2[3]),
    )


def bounds_union(
    bounds1: tuple[float, float, float, float],
    bounds2: tuple[float, float, float, float],
) -> tuple[float, float, float, float]:
    """
    Compute union of two bounding boxes.

    Args:
        bounds1: First bounds (xmin, ymin, xmax, ymax)
        bounds2: Second bounds (xmin, ymin, xmax, ymax)

    Returns:
        Union bounds
    """
    return (
        min(bounds1[0], bounds2[0]),
        min(bounds1[1], bounds2[1]),
        max(bounds1[2], bounds2[2]),
        max(bounds1[3], bounds2[3]),
    )


def bounds_area(bounds: tuple[float, float, float, float]) -> float:
    """
    Compute area of bounding box.

    Args:
        bounds: Bounds (xmin, ymin, xmax, ymax)

    Returns:
        Area in square CRS units
    """
    width = bounds[2] - bounds[0]
    height = bounds[3] - bounds[1]
    return width * height


# =============================================================================
# Grid Utilities
# =============================================================================


def create_grid_from_bounds(
    bounds: tuple[float, float, float, float],
    resolution: float,
    crs: str | CRS,
) -> xr.DataArray:
    """
    Create an empty grid DataArray from bounds and resolution.

    Args:
        bounds: (xmin, ymin, xmax, ymax)
        resolution: Pixel size in CRS units
        crs: Coordinate reference system

    Returns:
        Empty DataArray with correct spatial coordinates
    """
    xmin, ymin, xmax, ymax = bounds

    # Calculate dimensions
    width = int(np.ceil((xmax - xmin) / resolution))
    height = int(np.ceil((ymax - ymin) / resolution))

    # Create coordinates (pixel centers)
    x = np.linspace(xmin + resolution / 2, xmin + (width - 0.5) * resolution, width)
    y = np.linspace(ymax - resolution / 2, ymax - (height - 0.5) * resolution, height)

    # Create DataArray
    data = np.zeros((height, width), dtype=np.float32)
    da = xr.DataArray(
        data,
        dims=["y", "x"],
        coords={"y": y, "x": x},
    )

    # Set CRS
    da = da.rio.write_crs(crs)

    # Set transform
    from rasterio.transform import from_bounds

    transform = from_bounds(xmin, ymin, xmax, ymax, width, height)
    da = da.rio.write_transform(transform)

    return da


def pixel_coords_to_geo(
    row: int | np.ndarray,
    col: int | np.ndarray,
    transform: tuple,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Convert pixel coordinates to geographic coordinates.

    Args:
        row: Row indices
        col: Column indices
        transform: Affine transform (from rasterio)

    Returns:
        Tuple of (x, y) geographic coordinates
    """
    from rasterio.transform import Affine

    if not isinstance(transform, Affine):
        transform = Affine(*transform[:6])

    col = np.asarray(col)
    row = np.asarray(row)

    x = transform.c + col * transform.a + row * transform.b
    y = transform.f + col * transform.d + row * transform.e

    return x, y


def geo_to_pixel_coords(
    x: float | np.ndarray,
    y: float | np.ndarray,
    transform: tuple,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Convert geographic coordinates to pixel coordinates.

    Args:
        x: X geographic coordinates
        y: Y geographic coordinates
        transform: Affine transform (from rasterio)

    Returns:
        Tuple of (row, col) pixel coordinates
    """
    from rasterio.transform import Affine

    if not isinstance(transform, Affine):
        transform = Affine(*transform[:6])

    inv_transform = ~transform

    x = np.asarray(x)
    y = np.asarray(y)

    col = inv_transform.a * x + inv_transform.b * y + inv_transform.c
    row = inv_transform.d * x + inv_transform.e * y + inv_transform.f

    return np.round(row).astype(int), np.round(col).astype(int)
