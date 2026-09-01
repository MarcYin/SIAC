#!/usr/bin/env python3
"""Summarise Sentinel-2 L2A SCL over an AOI to feed the Cloud Score+ prefilter.

``scout_prefilter`` decides *which* acquisitions deserve a Cloud Score+ request;
this module supplies the free evidence it ranks on. Both signals come from the
public Planetary Computer L2A archive, which is an anonymous HTTP read under no
quota, so the whole seasonal window can be surveyed before any Earth Engine
request is spent.

Two things keep the survey cheap enough to run at catalogue scale.

*The STAC footprint gates before the raster does.* Overlap with the AOI is
computable from the item geometry alone, so a swath-edge acquisition is rejected
without opening its COG at all. The footprint is an upper bound on real coverage
-- it ignores interior nodata -- so the gate can only reject acquisitions that
would have failed the same gate on pixels.

*The read is small.* A prefilter needs two scalars over the AOI, not a mask, so
the read is a window, not a scene: a 7.7 km AOI is 384x384 native SCL pixels,
inside a single 512x512 COG block.

Decimating that window further was tried and rejected on measurement. Against 35
real Planetary Computer acquisitions, reading overview level 2 instead of the
native band saved 13% serially and nothing at all at 8 concurrent workers -- the
read is latency-bound at this AOI size, not bandwidth-bound -- while reporting
the AOI 0.016 *more* clear on average (max error 0.076, RMSE 0.033), because
subsampling a categorical field at 160 m drops the small cloud and shadow
features that the 0.98 clarity gate turns on. It changed the shortlist. So the
survey reads native and ``choose_overview`` stays available only for an AOI far
larger than a COG block, where the trade would actually be worth making.
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import time
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np
from tools.aeronet_validation.build_l2a_scl_index import (
    SAS_URL,
    SCL_CLEAR,
    STAC_URL,
    _coverage_ratio,
    _imports,
    _item_day,
    _item_idx,
    _numeric_baseline,
    _sensing_token,
    _session,
)
from tools.aeronet_validation.cloudscore_index_policy import LIBRARY_YEARS, seasonal_windows
from tools.aeronet_validation.scout_prefilter import MIN_COVERAGE, ScoutRecord, shortlist

if TYPE_CHECKING:
    from collections.abc import Sequence

#: Native SCL pixel size.
SCL_NATIVE_RESOLUTION_M = 20.0

#: Scout grid resolution. Native, by measurement: decimating gains nothing
#: concurrently at AOI sizes that fit one COG block and biases the clear
#: fraction upward. See the module docstring.
SCOUT_RESOLUTION_M = SCL_NATIVE_RESOLUTION_M

#: Target used when a caller opts into decimation for a large AOI. 240 m is 12x
#: native: coarse enough to cut a tile-sized read by two orders of magnitude,
#: fine enough that overview level 2 (160 m) is still the one selected.
DECIMATED_SCOUT_RESOLUTION_M = 240.0

#: Clear classes for scouting: vegetation, bare, water, unclassified and snow.
#: This is ``SCL_CLEAR["standard"]``; snow counts as clear because a bright
#: surface is a valid library target, not a defect.
SCOUT_CLEAR_CLASSES = frozenset(SCL_CLEAR["standard"])

#: SCL nodata. Distinct from "not clear": it is absence of observation, and it
#: is what the coverage fraction measures.
SCL_NODATA = 0

DEFAULT_WORKERS = 8

#: Sentinel for "probe the collection once and reuse", for a survey whose AOI is
#: large enough that decimating pays. ``None`` -- the default -- means read the
#: native band.
_AUTO_OVERVIEW = object()


@dataclass(frozen=True)
class ScoutGrid:
    """Coarse AOI raster definition that every scouted acquisition warps onto."""

    crs: str
    transform: Any
    width: int
    height: int
    bounds: tuple[float, float, float, float]
    wgs84_polygon: dict[str, Any]

    @property
    def pixels(self) -> int:
        return int(self.width) * int(self.height)


def scout_grid(
    crs: str,
    transform: Any,
    height: int,
    width: int,
    *,
    resolution_m: float = SCOUT_RESOLUTION_M,
    imports: dict[str, Any] | None = None,
) -> ScoutGrid:
    """Decimate a scene grid to the scout resolution, preserving its extent.

    The extent is preserved exactly and the pixel count rounded up, so the scout
    grid always covers the whole AOI. A fraction measured on it is therefore a
    fraction of the AOI, not of some inscribed subset.
    """

    imported = imports if imports is not None else _imports()
    affine = transform
    scale_x = abs(float(affine.a))
    scale_y = abs(float(affine.e))
    if scale_x <= 0.0 or scale_y <= 0.0:
        raise ValueError("scene transform has non-positive pixel size")
    west = float(affine.c)
    north = float(affine.f)
    east = west + scale_x * int(width)
    south = north - scale_y * int(height)
    out_width = max(1, int(np.ceil(scale_x * int(width) / float(resolution_m))))
    out_height = max(1, int(np.ceil(scale_y * int(height) / float(resolution_m))))
    out_transform = imported["Affine"](
        float(resolution_m), 0.0, west, 0.0, -float(resolution_m), north
    )
    corners = [(0.0, 0.0), (out_width, 0.0), (out_width, out_height), (0.0, out_height)]
    xy = [(out_transform * point) for point in corners]
    to_wgs84 = imported["Transformer"].from_crs(crs, "EPSG:4326", always_xy=True).transform
    ring = [list(to_wgs84(*point)) for point in xy]
    ring.append(ring[0])
    return ScoutGrid(
        crs=str(crs),
        transform=out_transform,
        width=out_width,
        height=out_height,
        bounds=(min(west, east), min(north, south), max(west, east), max(north, south)),
        wgs84_polygon={"type": "Polygon", "coordinates": [ring]},
    )


#: Keys that carry the scene raster, in preference order. Teacher archives store
#: ``comp`` as (realization, band, y, x); index archives store ``winners`` as
#: (month, y, x). Only the trailing two axes are spatial in either case.
GRID_KEYS = ("surface", "comp", "winners", "boa")


def archive_grid_spec(path: Path) -> tuple[str, tuple[float, ...], int, int]:
    """Read ``(crs, transform, height, width)`` from a teacher or index npz.

    Campaign archives are not uniform. Teacher archives store ``comp`` as
    (realization, band, y, x) and declare ``epsg``; index archives store
    ``winners`` as (month, y, x) and declare ``crs``. Reading ``shape[:2]`` off
    the wrong key sizes the grid from the realization and band axes, so the
    spatial axes are always taken from the trailing two.
    """

    archive = np.load(path, allow_pickle=False)
    transform = tuple(np.asarray(archive["transform"], dtype=np.float64).ravel()[:6])
    keys = [key for key in GRID_KEYS if key in archive.files]
    keys += [key for key in archive.files if key not in keys and np.asarray(archive[key]).ndim >= 2]
    if not keys:
        raise ValueError(f"{path}: no raster to take a grid shape from")
    height, width = np.asarray(archive[keys[0]]).shape[-2:]
    if "crs" in archive.files:
        crs = str(archive["crs"])
    elif "epsg" in archive.files:
        crs = f"EPSG:{int(np.asarray(archive['epsg']).item())}"
    else:
        raise ValueError(f"{path}: archive declares neither crs nor epsg")
    return crs, transform, int(height), int(width)


def grid_from_archive(path: Path, *, imports: dict[str, Any] | None = None) -> ScoutGrid:
    """Build a scout grid from a teacher or index npz."""

    imported = imports if imports is not None else _imports()
    crs, transform, height, width = archive_grid_spec(path)
    return scout_grid(
        crs,
        imported["Affine"](*transform),
        height,
        width,
        imports=imported,
    )


def choose_overview(
    decimations: Sequence[int],
    *,
    native_resolution_m: float = SCL_NATIVE_RESOLUTION_M,
    target_resolution_m: float = DECIMATED_SCOUT_RESOLUTION_M,
) -> int | None:
    """Index of the coarsest overview still at least as fine as the scout grid.

    Returning ``None`` means read the full-resolution band. Never returns an
    overview coarser than the target: warping up from coarser data would alias
    the class field, and the bytes saved past that point are negligible anyway.

    Only reached when a caller sets ``overview_level``; the default survey reads
    native. See the module docstring for why.
    """

    limit = float(target_resolution_m) / float(native_resolution_m)
    usable = [
        index for index, decimation in enumerate(decimations) if 0 < float(decimation) <= limit
    ]
    if not usable:
        return None
    return max(usable, key=lambda index: float(decimations[index]))


def summarize_scl(
    scl: np.ndarray,
    *,
    clear_classes: frozenset[int] = SCOUT_CLEAR_CLASSES,
) -> tuple[float, float]:
    """Return ``(clear_fraction, coverage_fraction)`` over the whole array.

    Clear fraction is against the AOI, not against the observed part of it: an
    acquisition covering half the AOI perfectly is half clear, because the
    missing half cannot contribute a composite pixel. The coverage gate makes
    that distinction moot for anything that survives, but the definition has to
    be right for the diagnostics.
    """

    values = np.asarray(scl)
    if values.size == 0:
        return 0.0, 0.0
    covered = values != SCL_NODATA
    clear = np.isin(values, sorted(clear_classes))
    return float(np.mean(clear)), float(np.mean(covered))


def _asset_url(item: dict[str, Any], sas: str) -> str:
    href = str(item["assets"]["SCL"]["href"])
    # An empty token means the href is already usable -- a public mirror or a
    # local file -- and appending a bare "?" would corrupt a filesystem path.
    if not sas:
        return href
    return f"{href}&{sas}" if "?" in href else f"{href}?{sas}"


def _gdal_env() -> Any:
    import rasterio

    return rasterio.Env(
        GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR",
        GDAL_HTTP_MAX_RETRY="3",
        GDAL_HTTP_RETRY_DELAY="2",
        CPL_VSIL_CURL_ALLOWED_EXTENSIONS=".tif",
    )


def probe_overview_level(
    item: dict[str, Any],
    *,
    sas: str,
    target_resolution_m: float = DECIMATED_SCOUT_RESOLUTION_M,
) -> int | None:
    """Overview level to scout at, decided once and reused across a survey.

    Every SCL COG in the collection is built the same way, so probing per item
    would double the HTTP opens -- the dominant cost of a small-AOI read -- to
    re-derive a constant.
    """

    import rasterio

    with _gdal_env(), rasterio.open(_asset_url(item, sas)) as probe:
        native = abs(float(probe.transform.a)) or SCL_NATIVE_RESOLUTION_M
        return choose_overview(
            probe.overviews(1),
            native_resolution_m=native,
            target_resolution_m=target_resolution_m,
        )


def read_scl(
    item: dict[str, Any],
    *,
    sas: str,
    grid: ScoutGrid,
    imports: dict[str, Any],
    overview_level: int | None = None,
) -> np.ndarray:
    """Warp one item's SCL onto the scout grid, reading the given overview."""

    import rasterio

    url = _asset_url(item, sas)
    with _gdal_env():
        opener = (
            rasterio.open(url)
            if overview_level is None
            else rasterio.open(url, overview_level=int(overview_level))
        )
        with (
            opener as source,
            imports["WarpedVRT"](
                source,
                crs=grid.crs,
                transform=grid.transform,
                width=grid.width,
                height=grid.height,
                resampling=imports["Resampling"].nearest,
                nodata=SCL_NODATA,
            ) as vrt,
        ):
            return vrt.read(1)


def query_items(
    session: Any,
    grid: ScoutGrid,
    windows: Sequence[tuple[str, str]],
) -> list[dict[str, Any]]:
    """Every L2A item intersecting the AOI in the seasonal windows.

    Reprocessed duplicates collapse to the highest baseline and newest
    generation, matching ``build_l2a_scl_index``: two products of one
    acquisition are one candidate, and shortlisting both would waste the Cloud
    Score+ request the prefilter exists to save.
    """

    items: dict[str, dict[str, Any]] = {}
    for start, end in windows:
        payload = {
            "collections": ["sentinel-2-l2a"],
            "datetime": f"{start}T00:00:00Z/{end}T23:59:59Z",
            "intersects": grid.wgs84_polygon,
            "limit": 1000,
        }
        response = session.post(f"{STAC_URL}/search", json=payload, timeout=120)
        response.raise_for_status()
        for feature in response.json().get("features", []) or []:
            if not isinstance(feature, dict) or not feature.get("id"):
                continue
            if "SCL" not in (feature.get("assets") or {}):
                continue
            items[str(feature["id"])] = feature
    selected: dict[tuple[str, str], dict[str, Any]] = {}
    for item in items.values():
        properties = item.get("properties") or {}
        key = (_sensing_token(item), str(properties.get("s2:mgrs_tile") or ""))
        rank = (
            _numeric_baseline(properties.get("s2:processing_baseline")),
            str(properties.get("s2:generation_time") or ""),
        )
        previous = selected.get(key)
        if previous is None:
            selected[key] = item
            continue
        previous_properties = previous.get("properties") or {}
        previous_rank = (
            _numeric_baseline(previous_properties.get("s2:processing_baseline")),
            str(previous_properties.get("s2:generation_time") or ""),
        )
        if rank > previous_rank:
            selected[key] = item
    return sorted(selected.values(), key=lambda item: (_item_day(item), _item_idx(item)))


def scout(
    items: Sequence[dict[str, Any]],
    *,
    sas: str,
    grid: ScoutGrid,
    imports: dict[str, Any],
    aod_by_day: dict[str, float | None] | None = None,
    min_footprint: float = MIN_COVERAGE,
    workers: int = DEFAULT_WORKERS,
    overview_level: int | None | object = None,
) -> tuple[tuple[ScoutRecord, ...], dict[str, Any]]:
    """Scout every item, skipping COG reads the footprint already rules out."""

    lookup = dict(aod_by_day or {})
    gated: list[dict[str, Any]] = []
    footprint_rejected = 0
    for item in items:
        ratio = _coverage_ratio(item, grid.bounds, grid.crs, imports)
        if ratio < float(min_footprint):
            footprint_rejected += 1
            continue
        gated.append(item)

    records: list[ScoutRecord] = []
    failures: dict[str, str] = {}

    level = overview_level
    if level is _AUTO_OVERVIEW:
        level = probe_overview_level(gated[0], sas=sas) if gated else None

    def one(item: dict[str, Any]) -> tuple[str, ScoutRecord | None, str | None]:
        item_id = str(item["id"])
        try:
            values = read_scl(item, sas=sas, grid=grid, imports=imports, overview_level=level)
        except Exception as error:  # noqa: BLE001 - one unreadable COG must not end the survey
            return item_id, None, f"{type(error).__name__}: {error}"
        clear, coverage = summarize_scl(values)
        day = _item_day(item)
        return (
            item_id,
            ScoutRecord(
                image_id=_item_idx(item),
                day=day,
                clear_fraction=clear,
                coverage_fraction=coverage,
                aod=lookup.get(day),
            ),
            None,
        )

    started = time.perf_counter()
    with ThreadPoolExecutor(max_workers=max(1, int(workers))) as pool:
        for item_id, record, error in pool.map(one, gated):
            if error is not None:
                failures[item_id] = error
                continue
            if record is not None:
                records.append(record)
    elapsed = time.perf_counter() - started

    diagnostics = {
        "items_found": len(items),
        "footprint_rejected": footprint_rejected,
        "scl_reads": len(gated),
        "scl_read_failures": len(failures),
        "read_seconds": elapsed,
        "seconds_per_read": elapsed / len(gated) if gated else 0.0,
        "overview_level": level,
        "failures": failures,
        "scout_grid": {"width": grid.width, "height": grid.height, "pixels": grid.pixels},
    }
    return tuple(sorted(records, key=lambda r: (r.day, r.image_id))), diagnostics


def maiac_by_day(grid: ScoutGrid, days: Sequence[str]) -> dict[str, float | None]:
    """Per-day MAIAC AOD over the AOI; days without a retrieval return None."""

    from siac.adapters.atmo.maiac_day_aod import MAIACDayAODProvider

    if not days:
        return {}
    west, south, east, north = grid.bounds
    periods = sorted({(int(day[:4]), int(day[5:7])) for day in days})
    resolved = MAIACDayAODProvider().day_aod_map((west, south, east, north), grid.crs, periods)
    return {day: resolved.get(day) for day in days}


def run(args: argparse.Namespace) -> dict[str, Any]:
    imports = _imports()
    grid = grid_from_archive(Path(args.grid_archive), imports=imports)
    scene_day = (
        dt.datetime.strptime(args.matchup_id.split("_")[-1][:8], "%Y%m%d").date().isoformat()
    )
    windows = seasonal_windows(scene_day, library_years=tuple(args.library_years))
    session = _session(imports["requests"])
    sas = str(session.get(SAS_URL, timeout=60).json()["token"])

    t0 = time.perf_counter()
    items = query_items(session, grid, windows)
    stac_seconds = time.perf_counter() - t0
    if args.limit:
        items = items[: int(args.limit)]

    records, diagnostics = scout(
        items,
        sas=sas,
        grid=grid,
        imports=imports,
        workers=args.workers,
    )
    if args.maiac:
        aod = maiac_by_day(grid, sorted({r.day for r in records}))
        records = tuple(
            ScoutRecord(
                image_id=r.image_id,
                day=r.day,
                clear_fraction=r.clear_fraction,
                coverage_fraction=r.coverage_fraction,
                aod=aod.get(r.day),
            )
            for r in records
        )

    selected = shortlist(records, top_k=args.top_k)
    cloud_score_requests = sum(len(v) for v in selected.values())
    payload = {
        "matchup_id": args.matchup_id,
        "scene_day": scene_day,
        "windows": len(windows),
        "library_years": list(args.library_years),
        "stac_seconds": stac_seconds,
        **diagnostics,
        "months_shortlisted": len(selected),
        "cloud_score_requests": cloud_score_requests,
        "requests_without_prefilter": diagnostics["items_found"],
        "shortlist": {
            month: [
                {
                    "image_id": r.image_id,
                    "day": r.day,
                    "clear_fraction": round(r.clear_fraction, 4),
                    "coverage_fraction": round(r.coverage_fraction, 4),
                    "aod": r.aod,
                }
                for r in entries
            ]
            for month, entries in selected.items()
        },
    }
    if args.output:
        Path(args.output).parent.mkdir(parents=True, exist_ok=True)
        Path(args.output).write_text(json.dumps(payload, indent=2), encoding="utf-8")
    return payload


def parser() -> argparse.ArgumentParser:
    value = argparse.ArgumentParser(description=__doc__)
    value.add_argument("--matchup-id", required=True)
    value.add_argument("--grid-archive", required=True, help="npz supplying the AOI grid.")
    value.add_argument("--library-years", type=int, nargs="+", default=list(LIBRARY_YEARS))
    value.add_argument("--top-k", type=int, default=5)
    value.add_argument("--workers", type=int, default=DEFAULT_WORKERS)
    value.add_argument("--limit", type=int, default=0, help="Scout only the first N items.")
    value.add_argument("--maiac", action="store_true", help="Attach MAIAC day AOD.")
    value.add_argument("--output")
    return value


def main() -> None:
    args = parser().parse_args()
    payload = run(args)
    payload.pop("shortlist", None)
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
