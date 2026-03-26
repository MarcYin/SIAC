"""
Area of Interest (AOI) container and constructors.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np

from siac.geo.geometry import bounds_to_polygon, load_aoi, polygon_to_bounds
from siac.geo.reprojection import transform_bounds

if TYPE_CHECKING:
    import xarray as xr
    from pyproj import CRS


@dataclass(frozen=True)
class AOI:
    """Area of Interest geometry + CRS."""

    geometry: dict[str, Any]
    crs: str = "EPSG:4326"

    @classmethod
    def from_bounds(
        cls,
        bounds: tuple[float, float, float, float],
        crs: str | CRS | None = None,
    ) -> AOI:
        """Create an AOI from bounds."""
        return cls(geometry=bounds_to_polygon(bounds), crs=_resolve_crs(bounds, crs))

    @classmethod
    def from_geojson(
        cls,
        aoi: str | Path | dict[str, Any],
        crs: str | CRS | None = None,
    ) -> AOI:
        """Create an AOI from GeoJSON-like input or a file path."""
        geometry = load_aoi(aoi)
        bounds = polygon_to_bounds(geometry)
        detected_crs = _detect_crs(aoi)
        return cls(geometry=geometry, crs=_resolve_crs(bounds, detected_crs or crs))

    @classmethod
    def from_raster(cls, raster: xr.DataArray) -> AOI:
        """Create an AOI from a raster extent."""
        if not hasattr(raster, "rio"):
            raise ValueError("Raster is missing rioxarray accessor; cannot derive AOI.")
        if raster.rio.crs is None:
            raise ValueError("Raster CRS is not set; cannot derive AOI.")

        bounds = raster.rio.bounds()
        return cls(geometry=bounds_to_polygon(bounds), crs=str(raster.rio.crs))

    def get_bounds(self, target_crs: str | CRS | None = None) -> tuple[float, float, float, float]:
        """Return AOI bounds, optionally transformed to `target_crs`."""
        bounds = polygon_to_bounds(self.geometry)
        if target_crs is None or str(target_crs) == self.crs:
            return bounds
        return transform_bounds(bounds, self.crs, target_crs)

    def to_geojson(self) -> dict[str, Any]:
        """Return AOI geometry as a GeoJSON geometry dict."""
        return self.geometry


def _detect_crs(aoi: str | Path | dict[str, Any]) -> str | None:
    """Best-effort CRS detection from a GeoJSON-like source."""
    payload: dict[str, Any] | None = None

    if isinstance(aoi, dict):
        payload = aoi
    elif isinstance(aoi, (str, Path)):
        path = Path(str(aoi))
        if path.exists():
            try:
                payload = json.loads(path.read_text())
            except (OSError, json.JSONDecodeError):
                return None

    if not isinstance(payload, dict):
        return None

    crs_value = payload.get("crs")
    if isinstance(crs_value, str):
        return crs_value
    if isinstance(crs_value, dict):
        props = crs_value.get("properties")
        if isinstance(props, dict):
            name = props.get("name")
            if isinstance(name, str):
                return name
    return None


def _resolve_crs(
    bounds: tuple[float, float, float, float],
    crs: str | CRS | None,
) -> str:
    if crs is None:
        _validate_default_wgs84_bounds(bounds)
        return "EPSG:4326"

    from pyproj import CRS as PyprojCRS

    return str(PyprojCRS.from_user_input(crs))


def _validate_default_wgs84_bounds(bounds: tuple[float, float, float, float]) -> None:
    xmin, ymin, xmax, ymax = (float(value) for value in bounds)
    if not all(np.isfinite([xmin, ymin, xmax, ymax])):
        raise ValueError("AOI bounds must be finite when defaulting CRS to WGS84")
    if xmin >= xmax or ymin >= ymax:
        raise ValueError("AOI bounds must satisfy min < max when defaulting CRS to WGS84")
    longitudes_are_signed = xmin >= -180.0 and xmax <= 180.0
    longitudes_are_wrapped = xmin >= 0.0 and xmax <= 360.0
    if not (longitudes_are_signed or longitudes_are_wrapped):
        raise ValueError(
            "AOI longitude bounds must look like degrees in [-180, 180] or [0, 360] when CRS is omitted"
        )
    if ymin < -90.0 or ymax > 90.0:
        raise ValueError("AOI latitude bounds must be within [-90, 90] when CRS is omitted")
