"""Water-mask loading with on-demand remote tile caching.

Network downloads share a single retry-enabled :class:`requests.Session` so
transient HTTP errors are retried with backoff and VRT companion ``.tif``
tiles avoid recreating a new Session per request (REVIEW.md §2.6, §3.3
water_mask.py:182-205).
"""

from __future__ import annotations

import logging
import shutil
from pathlib import Path
from tempfile import NamedTemporaryFile
from typing import TYPE_CHECKING
from urllib.parse import urlparse, urlunparse

import numpy as np
import xarray as xr

if TYPE_CHECKING:
    import requests

from siac.adapters._http import make_session
from siac.geo.reprojection import transform_bounds
from siac.runtime.models import copy_spatial_metadata_like
from siac.storage.readers import read_raster_window

logger = logging.getLogger(__name__)

DEFAULT_WATER_MASK_VRT_URL = (
    "https://zenodo.org/records/14899246/files/landWater2020.vrt?download=1"
)
DEFAULT_WATER_MASK_CACHE_DIR = Path.home() / ".cache" / "siac" / "water-mask"

_WGS84_CRS = "EPSG:4326"
_TILE_SIZE_DEG = 30
_LON_TILE_STARTS = tuple(range(-180, 180, _TILE_SIZE_DEG))
_LAT_TILE_STARTS = tuple(range(-90, 90, _TILE_SIZE_DEG))


def default_water_mask_cache_dir(cache_root: Path | None = None) -> Path:
    """Return the cache directory used for downloaded water-mask assets."""
    if cache_root is None:
        return DEFAULT_WATER_MASK_CACHE_DIR
    return Path(cache_root).expanduser() / "water-mask"


def load_water_mask_subset(
    bounds: tuple[float, float, float, float],
    crs: str,
    *,
    source: str | Path | None = None,
    cache_dir: Path | str | None = None,
    session: requests.Session | None = None,
) -> xr.DataArray:
    """Load a cropped boolean water mask for ``bounds``.

    The source raster stores ``1`` for land and ``0`` for water. The returned
    mask uses ``True`` for water pixels so callers can OR/aggregate it like any
    other exclusion mask.
    """
    local_source = ensure_local_water_mask_source(
        bounds,
        crs,
        source=source,
        cache_dir=cache_dir,
        session=session,
    )
    data = read_raster_window(local_source, bounds=bounds, bounds_crs=crs, masked=True)
    if "band" in data.dims and data.sizes["band"] == 1:
        data = data.squeeze("band", drop=True)
    water = xr.DataArray(
        np.asarray(data.values, dtype=np.float32) == 0.0,
        dims=data.dims,
        coords=data.coords,
        attrs=data.attrs,
        name="water_mask",
    )
    return copy_spatial_metadata_like(water, data)


def ensure_local_water_mask_source(
    bounds: tuple[float, float, float, float],
    crs: str,
    *,
    source: str | Path | None = None,
    cache_dir: Path | str | None = None,
    session: requests.Session | None = None,
) -> Path:
    """Resolve a water-mask source to a local raster/VRT path."""
    resolved_source = source or DEFAULT_WATER_MASK_VRT_URL
    explicit_cache_dir = Path(cache_dir).expanduser() if cache_dir is not None else None

    if isinstance(resolved_source, Path):
        local_source = resolved_source.expanduser()
    else:
        parsed = urlparse(str(resolved_source))
        local_source = (
            None
            if parsed.scheme.lower() in {"http", "https"}
            else Path(str(resolved_source)).expanduser()
        )

    if local_source is not None:
        if local_source.is_dir():
            candidate = local_source / "landWater2020.vrt"
            if candidate.exists():
                return candidate
            if explicit_cache_dir is None:
                explicit_cache_dir = local_source
            resolved_source = DEFAULT_WATER_MASK_VRT_URL
        elif local_source.exists():
            return local_source
        else:
            if explicit_cache_dir is None:
                explicit_cache_dir = local_source.parent if local_source.suffix else local_source
            resolved_source = DEFAULT_WATER_MASK_VRT_URL

    remote_source = str(resolved_source)
    remote_path = urlparse(remote_source).path
    remote_suffix = Path(remote_path).suffix.lower()
    if remote_suffix not in {".vrt", ".tif", ".tiff"}:
        raise ValueError(f"Unsupported water-mask source: {remote_source!r}")

    cache_root = explicit_cache_dir or default_water_mask_cache_dir()
    cache_root.mkdir(parents=True, exist_ok=True)
    local_path = cache_root / Path(remote_path).name

    # Reuse one retry-enabled session for the VRT plus all companion tiles
    # (REVIEW.md §3.3 water_mask.py:182-205).
    owns_session = session is None
    client = session if session is not None else make_session()
    try:
        if not local_path.exists():
            _download_file(remote_source, local_path, session=client)

        if remote_suffix != ".vrt":
            return local_path

        for tile_name in required_water_mask_tiles(bounds, crs):
            tile_path = cache_root / tile_name
            if tile_path.exists():
                continue
            _download_file(_tile_url(remote_source, tile_name), tile_path, session=client)
    finally:
        if owns_session:
            client.close()

    return local_path


def required_water_mask_tiles(
    bounds: tuple[float, float, float, float],
    crs: str,
) -> tuple[str, ...]:
    """Return the Zenodo tile names needed to cover ``bounds``."""
    xmin, ymin, xmax, ymax = transform_bounds(bounds, crs, _WGS84_CRS)
    lat_min = max(-90.0, min(float(ymin), float(ymax)))
    lat_max = min(90.0, max(float(ymin), float(ymax)))
    lon_intervals = _longitude_intervals(float(xmin), float(xmax))

    tiles: list[str] = []
    for lon_start in _LON_TILE_STARTS:
        lon_end = lon_start + _TILE_SIZE_DEG
        if not any(_intervals_intersect(l0, l1, lon_start, lon_end) for l0, l1 in lon_intervals):
            continue
        for lat_start in _LAT_TILE_STARTS:
            lat_end = lat_start + _TILE_SIZE_DEG
            if _intervals_intersect(lat_min, lat_max, lat_start, lat_end):
                tiles.append(_tile_filename(lon_start, lat_start))
    return tuple(tiles)


def _tile_filename(lon_start: int, lat_start: int) -> str:
    return f"landWater2020_{lon_start:04d}_{lat_start:03d}.tif"


def _longitude_intervals(xmin: float, xmax: float) -> tuple[tuple[float, float], ...]:
    if xmax - xmin >= 360.0:
        return ((-180.0, 180.0),)
    if xmax < xmin:
        return ((xmin, 180.0), (-180.0, xmax))
    if xmin < -180.0 or xmax > 180.0:
        shifted = (((xmin + 180.0) % 360.0) - 180.0, ((xmax + 180.0) % 360.0) - 180.0)
        if shifted[1] < shifted[0]:
            return ((shifted[0], 180.0), (-180.0, shifted[1]))
        return (shifted,)
    return ((xmin, xmax),)


def _intervals_intersect(a0: float, a1: float, b0: float, b1: float) -> bool:
    return not (a1 <= b0 or a0 >= b1)


def _tile_url(vrt_url: str, tile_name: str) -> str:
    parsed = urlparse(vrt_url)
    base_path = str(Path(parsed.path).parent / tile_name)
    return urlunparse(parsed._replace(path=base_path))


def _download_file(url: str, destination: Path, *, session: requests.Session | None = None) -> Path:
    """Stream *url* into *destination* via an atomic move.

    If *session* is ``None`` a retry-enabled session is created for this
    download only (REVIEW.md §2.6) and closed afterwards.
    """
    logger.info("Downloading water-mask asset %s -> %s", url, destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    owns_client = session is None
    client = session if session is not None else make_session()
    temp_path: Path | None = None
    try:
        with client.get(url, stream=True, timeout=120) as response:
            response.raise_for_status()
            with NamedTemporaryFile(delete=False, dir=destination.parent, suffix=".part") as handle:
                temp_path = Path(handle.name)
                for chunk in response.iter_content(chunk_size=1024 * 1024):
                    if chunk:
                        handle.write(chunk)
        assert temp_path is not None
        shutil.move(str(temp_path), str(destination))
        temp_path = None
        return destination
    finally:
        if temp_path is not None:
            temp_path.unlink(missing_ok=True)
        if owns_client:
            client.close()
