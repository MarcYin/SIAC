#!/usr/bin/env python3
"""Per-day MAIAC AOD over an AOI, reduced on Earth Engine instead of downloaded.

The earthaccess provider resolves the same quantity by downloading whole MCD19A2
HDF granules and reducing them locally. That is the right shape for one scene and
the wrong shape for a catalogue: measured at roughly 13 minutes per AOI, it was
the entire cost of the corpus index stage, dwarfing both the SCL survey and the
Cloud Score+ fetch it was meant to support.

The AOI is 7.68 km and MAIAC is 1 km, so the answer is a median over about 64
pixels per day. Asking Earth Engine to reduce server-side returns those medians
as numbers rather than as rasters, which is the difference between transferring
granules and transferring a list of floats.

Matches :class:`MAIACDayAODProvider`: ``Optical_Depth_055``, an area median per
day, at 1 km, returned in physical units.
"""

from __future__ import annotations

import datetime as dt
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from collections.abc import Sequence

COLLECTION = "MODIS/061/MCD19A2_GRANULES"
SERVER_URL = "https://earthengine-highvolume.googleapis.com"
BAND = "Optical_Depth_055"

#: MAIAC stores AOD as scaled integers; the catalogue documents 0.001 and Earth
#: Engine serves the raw values. The reference provider returns physical units,
#: so the scale is applied here rather than left for callers to remember.
AOD_SCALE = 0.001

#: Reduction resolution, matching the reference provider's AOI clip.
REDUCTION_SCALE_M = 1000.0

#: Days per Earth Engine request. A whole +/-45 day window in one call risks
#: "User memory limit exceeded" -- already hit once on this project by a large
#: getInfo -- so windows are split and each chunk retried smaller on failure.
CHUNK_DAYS = 46


def _geometry(bbox: tuple[float, float, float, float]) -> Any:
    import ee

    # edown's auth, which resolves ADC credentials and the quota project, rather
    # than a bare ee.Initialize. Note edown must be installed *into this
    # environment* (--no-deps): putting its bundled runtime on PYTHONPATH lets
    # its own pyproj and rasterio shadow this environment's and breaks the PROJ
    # database for every GDAL call in the process.
    from edown.auth import initialize_earth_engine

    initialize_earth_engine(SERVER_URL)
    west, south, east, north = bbox
    return ee.Geometry.Rectangle([float(west), float(south), float(east), float(north)])


def _reduce_window(geometry: Any, start: str, end: str) -> list[dict[str, Any]]:
    """One request: one AOI value per day in the range.

    Mirrors the reference reduction order exactly -- a per-pixel **mean** across
    the day's overpasses, then a spatial **median** over the AOI. The order is
    not interchangeable. Reducing each granule separately and combining the
    results afterwards gives a granule that merely clips the AOI corner the same
    weight as one covering it whole, which measured as a systematic +0.028 bias
    against the values the production pipeline stored.
    """

    import ee

    collection = (
        ee.ImageCollection(COLLECTION)
        .filterBounds(geometry)
        .filterDate(start, (dt.date.fromisoformat(end) + dt.timedelta(days=1)).isoformat())
        .select(BAND)
    )
    days = (
        collection.aggregate_array("system:time_start")
        .map(lambda millis: ee.Date(millis).format("YYYY-MM-dd"))
        .distinct()
    )

    def per_day(day: Any) -> Any:
        date = ee.Date(day)
        # ImageCollection.mean() ignores masked pixels, matching the reference's
        # nanmean over the orbit stack.
        daily = collection.filterDate(date, date.advance(1, "day")).mean()
        value = daily.reduceRegion(
            reducer=ee.Reducer.median(),
            geometry=geometry,
            scale=REDUCTION_SCALE_M,
            bestEffort=True,
            maxPixels=1e8,
        ).get(BAND)
        return ee.Feature(None, {"day": day, "aod": value})

    return list(ee.FeatureCollection(days.map(per_day)).getInfo().get("features", []) or [])


def _chunks(start: str, end: str, days: int = CHUNK_DAYS):
    first = dt.date.fromisoformat(start)
    last = dt.date.fromisoformat(end)
    while first <= last:
        stop = min(last, first + dt.timedelta(days=days - 1))
        yield first.isoformat(), stop.isoformat()
        first = stop + dt.timedelta(days=1)


def day_aod_map(
    bbox: tuple[float, float, float, float],
    windows: Sequence[tuple[str, str]],
) -> dict[str, float]:
    """``day -> physical AOD`` over the AOI across the seasonal windows.

    One value per day already, so the only combination left is across seasonal
    windows, which do not overlap.
    """

    from ee.ee_exception import EEException

    geometry = _geometry(bbox)
    by_day: dict[str, list[float]] = {}
    for start, end in windows:
        for chunk_start, chunk_end in _chunks(start, end):
            try:
                features = _reduce_window(geometry, chunk_start, chunk_end)
            except EEException:
                # Halve the chunk once rather than losing the whole window.
                features = []
                for half_start, half_end in _chunks(chunk_start, chunk_end, CHUNK_DAYS // 2):
                    features.extend(_reduce_window(geometry, half_start, half_end))
            for feature in features:
                properties = feature.get("properties") or {}
                value = properties.get("aod")
                day = properties.get("day")
                if value is None or day is None:
                    continue  # No retrieval over the AOI that overpass.
                by_day.setdefault(str(day), []).append(float(value) * AOD_SCALE)
    return {
        day: float(sorted(values)[len(values) // 2]) for day, values in by_day.items() if values
    }


def lookup_for(bbox: tuple[float, float, float, float], windows: Sequence[tuple[str, str]]):
    """An ``AODLookup`` over a pre-resolved window, for the index builder."""

    resolved = day_aod_map(bbox, windows)

    def lookup(days: Sequence[str]) -> dict[str, float | None]:
        return {day: resolved.get(day) for day in days}

    return lookup
