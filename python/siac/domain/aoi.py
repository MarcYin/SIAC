"""
Area of Interest (AOI) container and constructors.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Any

from siac.io.geometry import bounds_to_polygon, load_aoi, polygon_to_bounds
from siac.io.reprojection import transform_bounds

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
        crs: str | CRS = "EPSG:4326",
    ) -> AOI:
        """Create an AOI from bounds."""
        return cls(geometry=bounds_to_polygon(bounds), crs=str(crs))

    @classmethod
    def from_geojson(
        cls,
        aoi: str | Path | dict[str, Any],
        crs: str | CRS = "EPSG:4326",
    ) -> AOI:
        """Create an AOI from GeoJSON-like input or a file path."""
        geometry = load_aoi(aoi)
        detected_crs = _detect_crs(aoi) or str(crs)
        return cls(geometry=geometry, crs=detected_crs)

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
