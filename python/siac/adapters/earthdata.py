"""Earthdata adapter helpers."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, cast

import numpy as np
import xarray as xr

from siac.adapters.data.earthaccess_catalog import EarthAccessCatalog
from siac.adapters.data.earthaccess_source import EarthAccessSource
from siac.adapters.earthdata_common import (
    granule_intersects_bounds,
    parse_granule_date,
    parse_tile_indices,
    reproject_native_to_target,
)

if TYPE_CHECKING:
    from datetime import datetime

    from rasterio.enums import Resampling

    from siac.adapters.auth import CredentialManager


@dataclass(frozen=True)
class EarthAccessRuntime:
    """Normalized Earthaccess runtime dependencies for provider wrappers."""

    cache_dir: Path | None
    source: EarthAccessSource
    catalog: EarthAccessCatalog


def earthaccess_source_from_auth(
    auth: CredentialManager | None,
    *,
    provider: str | None = None,
) -> EarthAccessSource:
    """Build an ``EarthAccessSource`` using configured credentials when available."""
    if auth is None:
        return EarthAccessSource(provider=provider)
    return auth.earthdata().build_earthaccess_source(provider=provider)


def build_earthaccess_runtime(
    *,
    cache_dir: str | Path | None,
    source: EarthAccessSource | None,
    catalog: EarthAccessCatalog | None,
    provider: str | None,
) -> EarthAccessRuntime:
    """Build normalized source/catalog/cache dependencies for Earthaccess providers."""
    resolved_cache_dir = Path(cache_dir).expanduser() if cache_dir is not None else None
    resolved_source = source or EarthAccessSource(provider=provider)
    resolved_catalog = catalog or EarthAccessCatalog(source=resolved_source)
    return EarthAccessRuntime(
        cache_dir=resolved_cache_dir,
        source=resolved_source,
        catalog=resolved_catalog,
    )


def earthaccess_cache_dir(cache_dir: Path | None, short_name: str) -> Path:
    """Return the on-disk cache destination for a downloaded Earthaccess product."""
    return cache_dir or Path.home() / ".cache" / "siac" / "earthdata" / short_name


def select_candidate_paths(
    paths: list[Path],
    *,
    obs_time: datetime,
    bounds: tuple[float, float, float, float],
    crs: str,
    sample_dates: np.ndarray | None = None,
) -> list[Path]:
    """Rank downloaded granules by tile, temporal proximity, and filename."""
    if not paths:
        return []

    sample_day_set: set[np.datetime64] | None = None
    if sample_dates is not None:
        sample_day_set = {np.datetime64(day, "D") for day in sample_dates.tolist()}

    selected: list[tuple[tuple[int, int], float, str, Path]] = []
    for path in paths:
        try:
            granule_dt = parse_granule_date(path)
            granule_day = np.datetime64(granule_dt.date(), "D")
            if sample_day_set is not None and granule_day not in sample_day_set:
                continue
            delta = abs((granule_dt - obs_time).total_seconds())
            intersects = granule_intersects_bounds(path, bounds=bounds, crs=crs)
        except Exception:
            return paths

        if not intersects:
            continue

        try:
            tile = parse_tile_indices(path)
        except Exception:
            tile = (999, 999)
        selected.append((tile, delta, path.name, path))

    if not selected:
        return []
    return [item[3] for item in sorted(selected, key=lambda value: (value[0], value[1], value[2]))]


def target_grid_coords(
    bounds: tuple[float, float, float, float],
    resolution: float,
    *,
    resolution_name: str = "resolution",
) -> tuple[np.ndarray, np.ndarray]:
    """Build center-point y/x coordinates for a target raster grid."""
    xmin, ymin, xmax, ymax = bounds
    if resolution <= 0:
        raise ValueError(f"{resolution_name} must be > 0, got {resolution}")

    nx = max(1, int(np.ceil((xmax - xmin) / resolution)))
    ny = max(1, int(np.ceil((ymax - ymin) / resolution)))
    x = xmin + (np.arange(nx, dtype=np.float32) + 0.5) * resolution
    y = ymax - (np.arange(ny, dtype=np.float32) + 0.5) * resolution
    return y, x


def constant_target_array(
    bounds: tuple[float, float, float, float],
    resolution: float,
    value: float,
    *,
    resolution_name: str = "resolution",
) -> xr.DataArray:
    """Return a constant 2-D target-grid raster."""
    y, x = target_grid_coords(bounds, resolution, resolution_name=resolution_name)
    arr = cast("np.ndarray", np.full((y.size, x.size), value, dtype=np.float32))
    return xr.DataArray(arr, dims=["y", "x"], coords={"y": y, "x": x})


def constant_target_band_array(
    bands: list[int | str],
    bounds: tuple[float, float, float, float],
    resolution: float,
    value: float,
    *,
    resolution_name: str = "resolution",
) -> xr.DataArray:
    """Return a constant 3-D ``(band, y, x)`` target-grid raster."""
    y, x = target_grid_coords(bounds, resolution, resolution_name=resolution_name)
    arr = cast("np.ndarray", np.full((len(bands), y.size, x.size), value, dtype=np.float32))
    return xr.DataArray(
        arr,
        dims=["band", "y", "x"],
        coords={"band": np.array(bands, dtype=object), "y": y, "x": x},
    )


def merge_reprojected_tiles(
    arrays: list[xr.DataArray],
    *,
    bounds: tuple[float, float, float, float],
    crs: str,
    resolution: float,
    resampling: Resampling,
    nodata: float | None,
) -> xr.DataArray:
    """Reproject a list of native-tile rasters and merge them in priority order."""
    if not arrays:
        raise ValueError("Expected at least one array to merge")

    reprojected = [
        reproject_native_to_target(
            arr,
            target_bounds=bounds,
            target_crs=crs,
            target_resolution=resolution,
            resampling=resampling,
            nodata=nodata,
        )
        for arr in arrays
    ]
    merged = reprojected[0]
    for arr in reprojected[1:]:
        merged = merged.combine_first(arr)
    return cast("xr.DataArray", merged.astype(np.float32))


__all__ = [
    "EarthAccessRuntime",
    "build_earthaccess_runtime",
    "constant_target_array",
    "constant_target_band_array",
    "earthaccess_cache_dir",
    "earthaccess_source_from_auth",
    "merge_reprojected_tiles",
    "select_candidate_paths",
    "target_grid_coords",
]
