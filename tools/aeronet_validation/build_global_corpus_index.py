#!/usr/bin/env python3
"""Scout, fetch and compose a surface-library index for each catalogue AOI.

The three stages that were built and validated separately, run per AOI:

1. survey the seasonal window with free L2A SCL reads and shortlist it;
2. fetch Cloud Score+ for the shortlist alone, in one Earth Engine call;
3. compose the winner index locally from those planes.

The AOI grid is derived from the catalogue point rather than read from a teacher
archive, because no teacher exists yet -- this stage runs before capture. A
Sentinel-2 tile is always in the UTM zone its MGRS name encodes, so the CRS
follows from the tile id without opening a raster to ask.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import time
from pathlib import Path
from typing import TYPE_CHECKING, Any

from tools.aeronet_validation.build_cloudscore_index_local import (
    build_index,
    fetch_selected,
    grid_template,
    read_planes,
    write_index,
)
from tools.aeronet_validation.build_l2a_scl_index import _imports, _session
from tools.aeronet_validation.cloudscore_index_policy import LIBRARY_YEARS, seasonal_windows
from tools.aeronet_validation.maiac_gee_day_aod import day_aod_map
from tools.aeronet_validation.scout_prefilter import ScoutRecord, shortlist
from tools.aeronet_validation.scout_scl import SasToken, query_items, scout, scout_grid

if TYPE_CHECKING:
    from collections.abc import Sequence

#: AOI size on the 60 m solver grid, and the 20 m grid nested inside it. 128 is
#: a power of two at 60 m; 2^k is never divisible by 3, so a power of two at
#: 20 m cannot coexist with an integer 60 m grid.
TEMPLATE_SIZE = 128
TEMPLATE_RESOLUTION_M = 60.0
GRID_RESOLUTION_M = 20.0

#: MGRS latitude bands. N onwards is the northern hemisphere; the letters I and
#: O are unused, which is why this is a lookup rather than arithmetic.
MGRS_BANDS = "CDEFGHJKLMNPQRSTUVWX"


def crs_for_mgrs_tile(tile: str) -> str:
    """UTM CRS a Sentinel-2 tile is delivered in, from its MGRS name.

    Avoids opening a band just to read its CRS: the tile id already determines
    it, and at catalogue scale that is thousands of remote reads saved.
    """

    name = str(tile).upper().lstrip("T")
    if len(name) < 3 or not name[:2].isdigit():
        raise ValueError(f"not an MGRS tile id: {tile!r}")
    zone = int(name[:2])
    band = name[2]
    if band not in MGRS_BANDS:
        raise ValueError(f"invalid MGRS latitude band {band!r} in {tile!r}")
    if not 1 <= zone <= 60:
        raise ValueError(f"invalid UTM zone {zone} in {tile!r}")
    northern = band >= "N"
    return f"EPSG:{326 if northern else 327}{zone:02d}"


def grid_for_point(
    longitude: float,
    latitude: float,
    tile: str,
    *,
    imports: dict[str, Any],
    size: int = TEMPLATE_SIZE,
    template_resolution_m: float = TEMPLATE_RESOLUTION_M,
    resolution_m: float = GRID_RESOLUTION_M,
):
    """Grid spec for an AOI centred on a catalogue point.

    Snapped to the *60 m* template grid, not the 20 m one, so the two nest
    exactly -- the multiscale model's ``context60`` input depends on that
    nesting, and snapping to 20 m would let the coarse grid straddle it.
    """

    crs = crs_for_mgrs_tile(tile)
    transformer = imports["Transformer"].from_crs("EPSG:4326", crs, always_xy=True)
    center_x, center_y = transformer.transform(float(longitude), float(latitude))
    half = float(size) * float(template_resolution_m) / 2.0
    left = math.floor((center_x - half) / template_resolution_m) * template_resolution_m
    top = math.ceil((center_y + half) / template_resolution_m) * template_resolution_m
    span = int(round(float(size) * float(template_resolution_m) / float(resolution_m)))
    transform = imports["Affine"](float(resolution_m), 0.0, left, 0.0, -float(resolution_m), top)
    return crs, transform, span


def scene_day_of(datetime_text: str) -> str:
    return str(datetime_text)[:10]


def process_scene(
    row: dict[str, str],
    *,
    output_root: Path,
    cache_root: Path,
    edown: str,
    imports: dict[str, Any],
    sas: str,
    session: Any,
    library_years: Sequence[int],
    top_k: int,
    workers: int,
    use_maiac: bool,
) -> dict[str, Any]:
    """Run all three stages for one AOI, returning timings and counts."""

    matchup = row["matchup_id"]
    crs, transform, span = grid_for_point(
        float(row["longitude"]), float(row["latitude"]), row["mgrs_tile"], imports=imports
    )
    grid = scout_grid(crs, transform, span, span, resolution_m=GRID_RESOLUTION_M, imports=imports)
    template, bbox = grid_template(crs, tuple(transform)[:6], span, span)
    day = scene_day_of(row["datetime"])
    windows = seasonal_windows(day, library_years=library_years)

    started = time.perf_counter()
    items = query_items(session, grid, windows)
    t_stac = time.perf_counter() - started

    started = time.perf_counter()
    records, diagnostics = scout(items, sas=sas, grid=grid, imports=imports, workers=workers)
    t_scout = time.perf_counter() - started

    # One MAIAC resolution per scene, covering the whole window, reused for both
    # the scout's aerosol ranking and the index's day weights. Reduced on Earth
    # Engine rather than downloaded: the earthaccess path measured ~13 minutes
    # per AOI, which was the entire cost of this stage.
    started = time.perf_counter()
    aod = day_aod_map(bbox, windows) if use_maiac else {}
    t_aod = time.perf_counter() - started
    if aod:
        records = tuple(
            ScoutRecord(
                image_id=record.image_id,
                day=record.day,
                clear_fraction=record.clear_fraction,
                coverage_fraction=record.coverage_fraction,
                aod=aod.get(record.day),
            )
            for record in records
        )

    chosen = shortlist(records, top_k=top_k)
    image_ids = sorted({record.image_id for entries in chosen.values() for record in entries})
    if not image_ids:
        raise RuntimeError("scouting shortlisted no acquisitions")

    cache = cache_root / matchup
    started = time.perf_counter()
    fetch_selected(bbox, image_ids, cache, edown)
    t_fetch = time.perf_counter() - started

    started = time.perf_counter()
    planes = read_planes(cache, template)

    def lookup(days: Sequence[str]) -> dict[str, float | None]:
        return {day: aod.get(day) for day in days}

    payload = build_index(planes, lookup, library_years=tuple(library_years))
    write_index(payload, output_root / f"{matchup}.npz", matchup)
    t_index = time.perf_counter() - started

    return {
        "matchup_id": matchup,
        "acquisitions": diagnostics["items_found"],
        "footprint_rejected": diagnostics["footprint_rejected"],
        "scl_reads": diagnostics["scl_reads"],
        "shortlisted": len(image_ids),
        "planes": len(planes),
        "maiac_days": len(aod),
        "seconds": {
            "stac": round(t_stac, 2),
            "scout": round(t_scout, 2),
            "aod": round(t_aod, 2),
            "fetch": round(t_fetch, 2),
            "index": round(t_index, 2),
        },
        **payload["counts"],
    }


def run(args: argparse.Namespace) -> dict[str, Any]:
    with Path(args.scene_catalog).open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    if args.shard_count > 1:
        rows = [row for index, row in enumerate(rows) if index % args.shard_count == args.shard]
    if args.limit:
        rows = rows[: int(args.limit)]

    output_root = Path(args.output_root)
    cache_root = Path(args.cache_root)
    output_root.mkdir(parents=True, exist_ok=True)

    imports = _imports()
    session = _session(imports["requests"])
    # Renewed as the run proceeds: a token fetched once lapses ~45 minutes in
    # and takes every remaining AOI with it.
    token = SasToken(session)

    done: list[dict[str, Any]] = []
    failures: dict[str, str] = {}
    skipped = 0
    started = time.perf_counter()
    for position, row in enumerate(rows):
        target = output_root / f"{row['matchup_id']}.npz"
        if target.is_file() and not args.overwrite:
            skipped += 1
            continue
        try:
            done.append(
                process_scene(
                    row,
                    output_root=output_root,
                    cache_root=cache_root,
                    edown=args.edown,
                    imports=imports,
                    sas=token.value,
                    session=session,
                    library_years=tuple(args.library_years),
                    top_k=args.top_k,
                    workers=args.workers,
                    use_maiac=args.maiac,
                )
            )
        except Exception as error:  # noqa: BLE001 - one bad AOI must not end the shard
            failures[row["matchup_id"]] = f"{type(error).__name__}: {error}"
        if (position + 1) % 5 == 0:
            print(
                f"  {position + 1}/{len(rows)} AOIs, {len(done)} built, {len(failures)} failed",
                flush=True,
            )
    elapsed = time.perf_counter() - started

    def mean(key: str) -> float:
        values = [entry["seconds"][key] for entry in done]
        return round(sum(values) / len(values), 2) if values else 0.0

    summary = {
        "aois": len(rows),
        "built": len(done),
        "skipped_existing": skipped,
        "failed": len(failures),
        "failures": dict(list(failures.items())[:10]),
        "wall_seconds": round(elapsed, 1),
        "sas_renewals": token.renewals,
        "seconds_per_scene": round(elapsed / max(1, len(done)), 2),
        "mean_seconds": {
            stage: mean(stage) for stage in ("stac", "scout", "aod", "fetch", "index")
        },
        "mean_shortlisted": (
            round(sum(entry["shortlisted"] for entry in done) / len(done), 1) if done else 0.0
        ),
        "mean_acquisitions": (
            round(sum(entry["acquisitions"] for entry in done) / len(done), 1) if done else 0.0
        ),
    }
    if args.summary:
        Path(args.summary).parent.mkdir(parents=True, exist_ok=True)
        Path(args.summary).write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return summary


def parser() -> argparse.ArgumentParser:
    value = argparse.ArgumentParser(description=__doc__)
    value.add_argument("--scene-catalog", required=True)
    value.add_argument("--output-root", required=True)
    value.add_argument("--cache-root", required=True)
    value.add_argument("--edown", default="/home/users/marcyin/SIAC/tools/edown_runtime")
    value.add_argument("--library-years", type=int, nargs="+", default=list(LIBRARY_YEARS))
    value.add_argument("--top-k", type=int, default=5)
    value.add_argument("--workers", type=int, default=8)
    value.add_argument("--maiac", action="store_true", default=True)
    value.add_argument("--no-maiac", dest="maiac", action="store_false")
    value.add_argument("--overwrite", action="store_true")
    value.add_argument("--limit", type=int, default=0)
    value.add_argument("--shard", type=int, default=0)
    value.add_argument("--shard-count", type=int, default=1)
    value.add_argument("--summary")
    return value


def main() -> None:
    print(json.dumps(run(parser().parse_args()), indent=2))


if __name__ == "__main__":
    main()
