#!/usr/bin/env python3
"""MAIAC for the AERONET-lineage prep chain, from Earth Engine, reproduced exactly.

The lineage's two MAIAC inputs -- the per-day AOD that ranks winner-index days
and the M5 teacher's gridded AOD prior -- were built from MCD19A2 HDF granules
downloaded through Earthdata. This module rebuilds both without an Earthdata
login by splitting the job three ways, each part taken from where it is exact:

* WHICH granules: the lineage's own search replayed against NASA's public CMR
  catalogue (:mod:`tools.aeronet_validation.maiac_cmr_replay`). CMR's temporal
  semantics and ``count`` caps decide which days enter, and nothing else
  reproduces them.
* WHICH values: Earth Engine, for exactly those granules, matched down to the
  production version CMR named. Earth Engine keeps superseded versions beside
  their reprocessings; a missing version fails the scene instead of quietly
  substituting another.
* WHICH arithmetic: the lineage's own functions. Earth Engine returns each
  orbit's raw values on the granule's NATIVE MODLAND pixels, placed into a full
  1200 x 1200 tile with the coordinates the HDF path gives it, so GDAL's warp
  onto the AOI (``merge_reprojected_tiles``) sees the same geotransform it did.
  Letting Earth Engine resample instead flips the nearest native pixel whenever a
  target cell centre lies within float reach of a pixel edge (0.05-0.1% of 60 m
  cells, up to 0.24 AOD).

The day AOD follows the schema-1 ``MAIACDayAODProvider`` the winner indices
were built with: no QA filter, a float32 orbit mean per tile-day, a 1 km AOI
grid, the median over finite cells. The teacher prior follows
``MCD19AODProvider`` as committed on 2026-07-25: best-quality QA with the loose
fallback always on, the nearest valid orbit per pixel, first-valid merge in
(tile, day distance, name) order, and median gap fill.
"""

from __future__ import annotations

import datetime as dt
from dataclasses import dataclass
from typing import Any

import numpy as np

COLLECTION = "MODIS/061/MCD19A2_GRANULES"
SERVER_URL = "https://earthengine-highvolume.googleapis.com"
MODLAND_PIXELS = 1200
#: Sentinels written for masked or absent cells, chosen so the reference
#: decoding rejects them: outside the AOD and uncertainty valid ranges, and the
#: QA fill value.
_SENTINEL = {"Optical_Depth_055": -32768, "AOD_Uncertainty": -32768, "AOD_QA": 0}
#: MCD19A2 collection 6.1 HDF attributes, read from the granules themselves. The
#: sentinel stands in for the product fill value, which Earth Engine masks.
_ATTRS = {
    "Optical_Depth_055": {
        "_FillValue": -32768,
        "valid_range": [-100, 6000],
        "scale_factor": 0.001,
        "add_offset": 0.0,
    },
    "AOD_Uncertainty": {
        "_FillValue": -32768,
        "valid_range": [0, 30000],
        "scale_factor": 0.0001,
        "add_offset": 0.0,
    },
    "AOD_QA": {"_FillValue": 0, "valid_range": [1, 65535]},
}


class GranuleVersionUnavailable(RuntimeError):
    """Earth Engine lacks the exact granule version the catalogue named."""


class GranuleEmptyInEarthEngine(RuntimeError):
    """Earth Engine's copy of a granule holds no retrieval anywhere in the tile.

    Seen for MCD19A2.A2020245.h12v01 (version 2023140191025): the HDF the
    lineage downloaded has 836 valid AOD pixels across 11 orbits, while Earth
    Engine's copy of the same version is masked in every orbit over the whole
    tile. Read as "no retrieval", it would hand those days to the CAMS gap fill
    and publish a different index as if nothing were wrong, so it fails instead.
    """


@dataclass
class NativeGranule:
    """One tile-day's orbits on the native pixels around an AOI."""

    granule: Any  # maiac_cmr_replay.GranuleId
    times: tuple[dt.datetime, ...]
    values: dict[str, np.ndarray]  # band -> (orbits, rows, cols) raw, the block only
    window: tuple[int, int, int, int]  # r0, r1, c0, c1 inclusive, in the 1200 x 1200 tile
    x: np.ndarray  # the full tile's pixel-centre coordinates
    y: np.ndarray


_INITIALIZED = False


def _initialize() -> None:
    """Initialise Earth Engine once per process, through edown's credentials."""
    global _INITIALIZED
    if _INITIALIZED:
        return
    from edown.auth import initialize_earth_engine

    initialize_earth_engine(SERVER_URL)
    _INITIALIZED = True


def _tile_indices(tile: str) -> tuple[int, int]:
    return int(tile[1:3]), int(tile[4:6])


#: Earth Engine's code for the MODLAND sinusoidal projection MCD19A2 is served in.
_EE_SINUSOIDAL = "SR-ORG:6974"


def _native_window(tile: str, extent, crs: str, margin: int):
    """Tile pixel-centre coordinates and the native block covering ``extent``."""
    from pyproj import Transformer

    from siac.adapters.earthdata_common import MODLAND_SINUSOIDAL_CRS, modland_tile_coords

    h, v = _tile_indices(tile)
    x, y = modland_tile_coords(h, v, MODLAND_PIXELS, MODLAND_PIXELS)
    size = float(x[1] - x[0])
    left, top = float(x[0]) - size / 2.0, float(y[0]) + size / 2.0
    xmin, ymin, xmax, ymax = extent
    sin_x, sin_y = Transformer.from_crs(crs, MODLAND_SINUSOIDAL_CRS, always_xy=True).transform(
        [xmin, xmax, xmin, xmax], [ymin, ymin, ymax, ymax]
    )
    cols = [int(np.floor((value - left) / size)) for value in sin_x]
    rows = [int(np.floor((top - value) / size)) for value in sin_y]
    c0, c1 = max(0, min(cols) - margin), min(MODLAND_PIXELS - 1, max(cols) + margin)
    r0, r1 = max(0, min(rows) - margin), min(MODLAND_PIXELS - 1, max(rows) + margin)
    if c0 > c1 or r0 > r1:
        return None
    return x, y, size, left, top, (r0, r1, c0, c1)


def native_granules(
    granules: list[Any],
    extent: tuple[float, float, float, float],
    crs: str,
    bands: tuple[str, ...],
    *,
    margin: int = 2,
) -> dict[Any, NativeGranule]:
    """Raw orbit values of each granule on its native pixels around ``extent``.

    One Earth Engine request per tile covers every granule of that tile. Tiles
    that do not reach the AOI are left out: in a nearest-neighbour merge they
    contribute nothing but masked cells. A granule whose exact version Earth
    Engine lacks raises :class:`GranuleVersionUnavailable`.
    """
    from collections import defaultdict

    import ee

    _initialize()
    by_tile: dict[str, list[Any]] = defaultdict(list)
    for granule in granules:
        by_tile[granule.tile].append(granule)
    result: dict[Any, NativeGranule] = {}
    for tile, group in sorted(by_tile.items()):
        window = _native_window(tile, extent, crs, margin)
        if window is None:
            continue
        x, y, size, left, top, (r0, r1, c0, c1) = window
        inset = 0.25 * size
        # In the projection's CRS units (metres): a projection with its pixel
        # transform would make Earth Engine read these as pixel indices.
        region = ee.Geometry.Rectangle(
            [
                left + c0 * size + inset,
                top - (r1 + 1) * size + inset,
                left + (c1 + 1) * size - inset,
                top - r0 * size - inset,
            ],
            proj=ee.Projection(_EE_SINUSOIDAL),
            geodesic=False,
        )
        prefixes = {f"MCD19A2_{g.day}_{g.tile}_061_{g.production}_": g for g in group}
        dates = [g.date for g in group]
        collection = (
            ee.ImageCollection(COLLECTION)
            .filterDate(
                (min(dates) - dt.timedelta(days=1)).isoformat(),
                (max(dates) + dt.timedelta(days=2)).isoformat(),
            )
            .filter(
                ee.Filter.Or(*[ee.Filter.stringStartsWith("system:index", p) for p in prefixes])
            )
            .select(list(bands))
        )

        def sample(image: Any, region: Any = region) -> Any:
            out = {"index": image.get("system:index"), "time": image.get("system:time_start")}
            first = image.select(bands[0])
            # Retrievals anywhere in the tile, not just the AOI: an empty copy
            # of a granule cannot otherwise be told from a cloudy AOI.
            out["tile_valid"] = (
                first.mask()
                .gt(0)
                .reduceRegion(
                    reducer=ee.Reducer.sum(),
                    geometry=first.geometry(),
                    crs=first.projection(),
                    maxPixels=1e10,
                )
                .values()
                .get(0)
            )
            for band in bands:
                selected = image.select(band)
                out[f"{band}__projection"] = selected.projection()
                out[band] = (
                    selected.unmask(_SENTINEL[band])
                    .sampleRectangle(region=region, defaultValue=_SENTINEL[band])
                    .get(band)
                )
            return ee.Feature(None, out)

        features = ee.FeatureCollection(collection.map(sample)).getInfo().get("features", [])
        orbits: dict[str, list[dict[str, Any]]] = defaultdict(list)
        for feature in features or []:
            properties = feature.get("properties") or {}
            for prefix in prefixes:
                if str(properties.get("index", "")).startswith(prefix):
                    orbits[prefix].append(properties)
        for prefix, granule in prefixes.items():
            found = orbits.get(prefix)
            if not found:
                raise GranuleVersionUnavailable(f"Earth Engine has no {prefix}* images")
            # Earth Engine numbers a granule's orbits _01.._NN in the HDF's layer
            # order, which the float32 orbit mean depends on.
            found.sort(key=lambda item: str(item["index"]))
            if not any(float(item.get("tile_valid") or 0) > 0 for item in found):
                raise GranuleEmptyInEarthEngine(
                    f"{prefix}*: every orbit is masked over the whole tile in Earth Engine"
                )
            values: dict[str, list[np.ndarray]] = {band: [] for band in bands}
            for orbit in found:
                for band in bands:
                    info = orbit[f"{band}__projection"]
                    scale, _, tx, _, negative_scale, ty = info["transform"][:6]
                    if (
                        info.get("crs") != _EE_SINUSOIDAL
                        or abs(scale - size) > 1e-6 * size
                        or abs(negative_scale + size) > 1e-6 * size
                        or abs((tx - left) / size - round((tx - left) / size)) > 1e-6
                        or abs((top - ty) / size - round((top - ty) / size)) > 1e-6
                    ):
                        raise ValueError(f"{orbit['index']}/{band} is not on the MODLAND tile grid")
                    block = np.asarray(orbit[band], dtype=np.int32)
                    if block.shape != (r1 - r0 + 1, c1 - c0 + 1):
                        raise ValueError(f"{orbit['index']}: sampled {block.shape}, not the block")
                    values[band].append(block)
            result[granule] = NativeGranule(
                granule=granule,
                # The HDF's Orbit_time_stamp has minute resolution.
                times=tuple(
                    dt.datetime.fromtimestamp(
                        float(orbit["time"]) / 1000.0, tz=dt.timezone.utc
                    ).replace(second=0, microsecond=0)
                    for orbit in found
                ),
                values={band: np.stack(stack) for band, stack in values.items()},
                window=(r0, r1, c0, c1),
                x=x,
                y=y,
            )
    return result


def _tile_dataarray(block: np.ndarray, native: NativeGranule) -> Any:
    """A 2-D block result laid into its full tile, NaN elsewhere.

    Every per-pixel step runs on the block alone; only the warp needs the whole
    tile, so that GDAL sees the geotransform the HDF path gave it.
    """
    import xarray as xr

    from siac.adapters.earthdata_common import MODLAND_SINUSOIDAL_CRS

    r0, r1, c0, c1 = native.window
    full = np.full((MODLAND_PIXELS, MODLAND_PIXELS), np.nan, dtype=np.float32)
    full[r0 : r1 + 1, c0 : c1 + 1] = block
    array = xr.DataArray(full, dims=("y", "x"), coords={"x": native.x, "y": native.y})
    return array.rio.set_spatial_dims(x_dim="x", y_dim="y").rio.write_crs(MODLAND_SINUSOIDAL_CRS)


def native_granule_reaches(granule: Any, extent, crs: str) -> bool:
    """Whether the granule's tile reaches the AOI (cheap, no Earth Engine call)."""
    from pyproj import Transformer

    from siac.adapters.earthdata_common import MODLAND_SINUSOIDAL_CRS, modland_tile_coords

    h, v = _tile_indices(granule.tile)
    x, y = modland_tile_coords(h, v, MODLAND_PIXELS, MODLAND_PIXELS)
    size = float(x[1] - x[0])
    left, top = float(x[0]) - size / 2.0, float(y[0]) + size / 2.0
    xmin, ymin, xmax, ymax = extent
    sin_x, sin_y = Transformer.from_crs(crs, MODLAND_SINUSOIDAL_CRS, always_xy=True).transform(
        [xmin, xmax, xmin, xmax], [ymin, ymin, ymax, ymax]
    )
    right, bottom = left + MODLAND_PIXELS * size, top - MODLAND_PIXELS * size
    return not (max(sin_x) < left or min(sin_x) > right or max(sin_y) < bottom or min(sin_y) > top)


def lineage_day_aod(
    bounds: tuple[float, float, float, float],
    crs: str,
    year: int,
    month: int,
    *,
    resolution: float = 1000.0,
    max_granules: int = 64,
) -> dict[str, float]:
    """One calendar month of per-day MAIAC AOD, as the winner indices stored it."""
    from rasterio.enums import Resampling
    from tools.aeronet_validation.maiac_cmr_replay import month_granules

    from siac.adapters.earthdata import merge_reprojected_tiles
    from siac.adapters.earthdata_common import apply_scale_and_mask, reduce_orbit_stack

    _initialize()
    band = "Optical_Depth_055"
    granules = month_granules(bounds, crs, year, month, max_granules=max_granules)
    natives = native_granules(granules, bounds, crs, (band,))
    by_day: dict[str, list[Any]] = {}
    for granule in granules:
        native = natives.get(granule)
        if native is None:
            continue
        field = reduce_orbit_stack(apply_scale_and_mask(native.values[band], _ATTRS[band]))
        by_day.setdefault(granule.date.isoformat(), []).append(_tile_dataarray(field, native))
    day_aod: dict[str, float] = {}
    for day, tiles in by_day.items():
        merged = merge_reprojected_tiles(
            tiles,
            bounds=bounds,
            crs=crs,
            resolution=resolution,
            resampling=Resampling.nearest,
            nodata=np.nan,
        )
        values = np.asarray(merged.values, dtype=np.float64)
        finite = values[np.isfinite(values)]
        if finite.size:
            day_aod[day] = float(np.median(finite))
    return day_aod


def _whole_granule_has(granule: Any, loose: bool) -> bool:
    """Whether any pixel of any orbit, anywhere in the granule, passes the mask.

    The committed provider tests ``valid.any()`` over the whole 1200 x 1200
    tile before falling back to the loose mask, so this cannot be decided from
    the AOI block alone.
    """
    import ee

    prefix = f"MCD19A2_{granule.day}_{granule.tile}_061_{granule.production}_"
    collection = (
        ee.ImageCollection(COLLECTION)
        .filterDate(
            (granule.date - dt.timedelta(days=1)).isoformat(),
            (granule.date + dt.timedelta(days=2)).isoformat(),
        )
        .filter(ee.Filter.stringStartsWith("system:index", prefix))
    )

    def flag(image: Any) -> Any:
        aod, unc, qa = (image.select(b) for b in ("Optical_Depth_055", "AOD_Uncertainty", "AOD_QA"))
        finite = (
            aod.mask()
            .And(aod.gte(-100))
            .And(aod.lte(6000))
            .And(unc.mask())
            .And(unc.gte(0))
            .And(unc.lte(30000))
        )
        passing = (
            finite.And(qa.mask()).And(qa.gt(0))
            if loose
            else finite.And(
                qa.bitwiseAnd(7)
                .eq(1)
                .And(qa.rightShift(5).bitwiseAnd(7).eq(0))
                .And(qa.rightShift(8).bitwiseAnd(15).eq(0))
            )
        )
        value = (
            passing.unmask(0)
            .reduceRegion(
                reducer=ee.Reducer.max(),
                # Bands do not share one projection; these three are on the
                # 1 km MODLAND grid, the QA band stands for them.
                geometry=qa.geometry(),
                crs=qa.projection(),
                maxPixels=1e10,
            )
            .values()
            .get(0)
        )
        return ee.Feature(None, {"any": value})

    flags = ee.FeatureCollection(collection.map(flag)).aggregate_array("any").getInfo()
    return any(bool(value) for value in flags or [])


@dataclass
class TeacherMaiacPrior:
    aot: Any  # xr.DataArray on the AOI grid, gaps median-filled
    aot_unc: Any
    granules: list[str]
    fallback_granules: list[str]

    def atmospheric_state(self) -> Any:
        """The prior as ``MCD19AODProvider.get_prior`` returned it.

        MCD19A2's water vapour is not reproduced: the teacher keeps only the
        aerosol fields of this state and takes water vapour, ozone and
        elevation from the scene. The provider's no-TCWV defaults stand in.
        """
        import xarray as xr

        from siac.runtime import AtmosphericState

        aot = self.aot.astype(np.float32)
        return AtmosphericState(
            aot=aot,
            tcwv=xr.full_like(aot, 1.5).astype(np.float32),
            tco3=xr.full_like(aot, 0.30).astype(np.float32),
            aot_unc=self.aot_unc.astype(np.float32),
            tcwv_unc=xr.full_like(aot, 0.3).astype(np.float32),
            tco3_unc=xr.full_like(aot, 0.03).astype(np.float32),
            elevation=xr.zeros_like(aot).astype(np.float32),
        )


def lineage_teacher_prior(
    bounds: tuple[float, float, float, float],
    crs: str,
    obs_time: dt.datetime,
    *,
    resolution: float = 60.0,
    window_days: int = 2,
    max_granules: int = 8,
) -> TeacherMaiacPrior:
    """The teacher's MAIAC prior on its AOI grid, as the committed provider built it.

    Raises :class:`~siac.adapters.atmo.mcd19_earthaccess.NoAtmosphericDataError`
    exactly where the lineage found no MAIAC for the scene -- no granule reaching
    the AOI, or no QA-valid AOD on it -- which the teacher answers with its
    CAMS-only fallback. Every other failure (Earth Engine errors, a missing
    production version, an Earth Engine copy emptier than its HDF) propagates:
    the lineage had data there, so a CAMS-only prior would not reproduce it.
    """
    from rasterio.enums import Resampling
    from tools.aeronet_validation.maiac_cmr_replay import teacher_granules

    from siac.adapters.atmo.mcd19_earthaccess import (
        NoAtmosphericDataError,
        _maiac_best_quality_mask,
        _nearest_valid_orbit_indices,
        _select_orbit_values,
    )
    from siac.adapters.earthdata import merge_reprojected_tiles
    from siac.adapters.earthdata_common import apply_scale_and_mask

    _initialize()
    obs_naive = obs_time.replace(tzinfo=None)
    granules = [
        g
        for g in teacher_granules(
            bounds, crs, obs_naive, window_days=window_days, max_granules=max_granules
        )
        if native_granule_reaches(g, bounds, crs)
    ]
    # select_candidate_paths: tile (h, v), distance from the day's naive 00:00, name.
    granules.sort(
        key=lambda g: (
            _tile_indices(g.tile),
            abs((dt.datetime.combine(g.date, dt.time()) - obs_naive).total_seconds()),
            g.filename,
        )
    )
    bands = ("Optical_Depth_055", "AOD_Uncertainty", "AOD_QA")
    aot_tiles, unc_tiles, used, fallback = [], [], [], []
    natives = native_granules(granules, bounds, crs, bands)
    for granule in granules:
        native = natives.get(granule)
        if native is None:
            continue
        raw = native.values
        aod = apply_scale_and_mask(raw["Optical_Depth_055"], _ATTRS["Optical_Depth_055"])
        unc = apply_scale_and_mask(raw["AOD_Uncertainty"], _ATTRS["AOD_Uncertainty"])
        qa = apply_scale_and_mask(raw["AOD_QA"], _ATTRS["AOD_QA"])
        finite = np.isfinite(aod) & np.isfinite(unc)
        loose = finite & np.isfinite(qa) & (qa > 0)
        valid = finite & _maiac_best_quality_mask(raw["AOD_QA"])
        if (
            not valid.any()
            and not _whole_granule_has(granule, loose=False)
            and (loose.any() or _whole_granule_has(granule, loose=True))
        ):
            valid = loose
            fallback.append(granule.filename)
        aod = np.where(valid, aod, np.nan)
        unc = np.where(valid, unc, np.nan)
        indices = _nearest_valid_orbit_indices(valid, native.times, obs_naive)
        if indices is None:
            raise ValueError(f"{granule.filename}: orbit times do not match the orbit stack")
        aot_tiles.append(_tile_dataarray(_select_orbit_values(aod, indices), native))
        unc_tiles.append(_tile_dataarray(_select_orbit_values(unc, indices), native))
        used.append(granule.filename)
    if not aot_tiles:
        raise NoAtmosphericDataError("MCD19 has no granule reaching the requested AOI")

    def merge(tiles):
        return merge_reprojected_tiles(
            tiles,
            bounds=bounds,
            crs=crs,
            resolution=resolution,
            resampling=Resampling.nearest,
            nodata=np.nan,
        )

    aot, aot_unc = merge(aot_tiles), merge(unc_tiles)
    finite_aot = np.asarray(aot.values, dtype=np.float64)
    finite_aot = finite_aot[np.isfinite(finite_aot)]
    if finite_aot.size == 0:
        raise NoAtmosphericDataError(
            "MCD19 has no QA-valid AOD after reprojection to the requested AOI"
        )
    finite_unc = np.asarray(aot_unc.values, dtype=np.float64)
    finite_unc = finite_unc[np.isfinite(finite_unc)]
    aot_fill = float(np.median(finite_aot))
    unc_fill = max(float(np.median(finite_unc)) if finite_unc.size else 0.10, 0.05)
    return TeacherMaiacPrior(
        aot=aot.fillna(aot_fill).astype(np.float32),
        aot_unc=aot_unc.fillna(unc_fill).astype(np.float32),
        granules=used,
        fallback_granules=fallback,
    )
