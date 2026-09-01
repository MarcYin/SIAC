#!/usr/bin/env python3
"""Generate the global land-point pool the catalogue is sampled from.

Land cover is read once per coarse tile through ``edown`` rather than per
candidate: the existing AERONET pipeline learned this the hard way, cutting
118,140 per-matchup lookups to 738 per-site ones. At stratification scale a
single class per 500 m is ample, so the whole globe is a few hundred requests.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

import numpy as np

if TYPE_CHECKING:
    from collections.abc import Callable, Iterable, Sequence
from tools.aeronet_validation.global_catalog_sampler import Candidate

# edown downloads at the collection's native resolution and exposes no scale
# control, so ESA WorldCover (10 m) cannot be tiled at stratification scale --
# one 30-degree tile would be 330,000 pixels square. MODIS MCD12Q1 is natively
# 500 m, which makes a 10-degree tile about 2,200 pixels square and a few MB.
# Class labels are mapped onto the same names the AERONET manifests use so the
# two corpora stay directly comparable.
LANDCOVER_COLLECTION = "MODIS/061/MCD12Q1"
LANDCOVER_BAND = "LC_Type1"

#: IGBP (``LC_Type1``) codes mapped to the manifest class names.
LANDCOVER_CLASSES: dict[int, str] = {
    1: "tree_cover",
    2: "tree_cover",
    3: "tree_cover",
    4: "tree_cover",
    5: "tree_cover",
    6: "shrubland",
    7: "shrubland",
    8: "savanna",
    9: "savanna",
    10: "grassland",
    11: "herbaceous_wetland",
    12: "cropland",
    13: "built_up",
    14: "cropland",
    15: "snow_and_ice",
    16: "bare_sparse",
    17: "permanent_water",
}

#: Classes that cannot support a surface-reflectance teacher: open water has no
#: stable surface to predict, and permanent snow/ice is the state the recurrent
#: snow policy already handles separately.
EXCLUDED_CLASSES: frozenset[str] = frozenset({"permanent_water", "snow_and_ice"})

#: Coarse longitude/latitude boxes. Only used to spread the draw across regions,
#: so approximate boundaries are adequate and a point landing in a neighbouring
#: box costs nothing.
_CONTINENT_BOXES: tuple[tuple[str, float, float, float, float], ...] = (
    ("africa", -20.0, -35.0, 55.0, 38.0),
    ("europe", -25.0, 38.0, 45.0, 72.0),
    ("asia", 45.0, -12.0, 180.0, 82.0),
    ("oceania", 110.0, -50.0, 180.0, -12.0),
    ("north_america", -170.0, 12.0, -50.0, 84.0),
    ("south_america", -82.0, -56.0, -34.0, 12.0),
)


def continent_of(longitude: float, latitude: float) -> str:
    for name, west, south, east, north in _CONTINENT_BOXES:
        if west <= longitude <= east and south <= latitude <= north:
            return name
    return "other"


def global_tiles(step_degrees: float = 30.0) -> tuple[tuple[float, float, float, float], ...]:
    """Tile the globe for WorldCover retrieval, skipping the polar caps.

    Above 82 degrees Sentinel-2 has no usable acquisitions, so those tiles would
    only waste requests.
    """

    if step_degrees <= 0:
        raise ValueError("step_degrees must be positive")
    tiles = []
    latitude = -60.0
    while latitude < 82.0:
        longitude = -180.0
        while longitude < 180.0:
            tiles.append(
                (
                    longitude,
                    latitude,
                    min(longitude + step_degrees, 180.0),
                    min(latitude + step_degrees, 82.0),
                )
            )
            longitude += step_degrees
        latitude += step_degrees
    return tuple(tiles)


#: Draw this multiple of ``per_tile`` before clipping to the nominal tile. A
#: sinusoidal raster covering a lat/lon box is much larger than the box, so most
#: drawn pixels fall outside it and would otherwise leave tiles short.
DEFAULT_OVERSAMPLE = 10


def inside(longitude: float, latitude: float, bbox: tuple[float, float, float, float]) -> bool:
    west, south, east, north = bbox
    return west <= longitude <= east and south <= latitude <= north


def candidates_from_class_grid(
    codes: np.ndarray,
    *,
    to_lonlat: Callable[[np.ndarray, np.ndarray], tuple[Any, Any]],
    tile_id: str,
    per_tile: int,
    seed: int = 0,
    excluded: Iterable[str] = EXCLUDED_CLASSES,
    within: tuple[float, float, float, float] | None = None,
    oversample: int = DEFAULT_OVERSAMPLE,
) -> tuple[Candidate, ...]:
    """Draw candidate points from one tile's land-cover class raster.

    The draw is *balanced across classes present in the tile* rather than
    uniform over pixels. A uniform draw would return the tile's dominant class
    almost exclusively, which is how a global sample ends up looking like the
    AERONET one -- mostly vegetation -- and misses the sparse classes that the
    stratified catalogue then has nothing to select from.

    Coordinates come from ``to_lonlat(rows, columns)`` rather than from a pair of
    1-D axis arrays. Indexing separate longitude and latitude vectors is only
    valid on an equirectangular grid; MCD12Q1 is sinusoidal, where longitude
    depends on the row as well as the column -- at the same easting, the top and
    bottom of one tile differ by more than 200 degrees. Resolving per point also
    keeps the cost proportional to the draw rather than to the 67-megapixel grid.

    ``within`` clips to the nominal tile so overlapping rasters do not
    over-represent the regions they share.
    """

    codes = np.asarray(codes)
    if codes.ndim != 2:
        raise ValueError(f"class grid must be 2-D, got shape {codes.shape}")
    excluded = set(excluded)
    rng = np.random.default_rng(seed)
    by_class: dict[str, list[tuple[int, int]]] = {}
    for row, column in zip(*np.nonzero(np.isin(codes, list(LANDCOVER_CLASSES)))):
        name = LANDCOVER_CLASSES[int(codes[row, column])]
        if name in excluded:
            continue
        by_class.setdefault(name, []).append((int(row), int(column)))
    if not by_class:
        return ()

    order = sorted(by_class)
    for name in order:
        rng.shuffle(by_class[name])

    # Round-robin across classes, so truncating after the clip keeps the balance.
    limit = int(per_tile) * max(1, int(oversample))
    picks: list[tuple[int, int, str]] = []
    while len(picks) < limit and any(by_class[name] for name in order):
        for name in order:
            if len(picks) >= limit:
                break
            if by_class[name]:
                row, column = by_class[name].pop()
                picks.append((row, column, name))
    if not picks:
        return ()

    rows = np.asarray([pick[0] for pick in picks])
    columns = np.asarray([pick[1] for pick in picks])
    longitudes, latitudes = to_lonlat(rows, columns)
    selected: list[Candidate] = []
    for (row, column, name), longitude, latitude in zip(picks, longitudes, latitudes):
        if len(selected) >= per_tile:
            break
        longitude, latitude = float(longitude), float(latitude)
        if not (-180.0 <= longitude <= 180.0 and -90.0 <= latitude <= 90.0):
            continue
        if within is not None and not inside(longitude, latitude, within):
            continue
        selected.append(
            Candidate(
                sample_id=f"{tile_id}_{name}_{row:05d}_{column:05d}",
                longitude=longitude,
                latitude=latitude,
                land_cover=name,
                continent=continent_of(longitude, latitude),
            )
        )
    return tuple(selected)


def pool_composition(pool: Sequence[Candidate]) -> dict[str, int]:
    counts: dict[str, int] = {}
    for candidate in pool:
        counts[candidate.land_cover] = counts.get(candidate.land_cover, 0) + 1
    return dict(sorted(counts.items(), key=lambda kv: -kv[1]))
