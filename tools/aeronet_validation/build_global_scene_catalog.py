#!/usr/bin/env python3
"""Pick one current-date Sentinel-2 acquisition per catalogue AOI.

The AOI catalogue is a set of land points; the corpus needs a specific L1C
acquisition at each. Selection runs against Earth Engine's
``COPERNICUS/S2_HARMONIZED`` rather than a STAC catalogue because it is the only
source that names the **L1C** product directly: ``PRODUCT_ID`` is the SAFE name
the capture stage fetches from GCS, and no reliable mapping exists from an L2A
item to its L1C sibling. Only image listings are requested here, never pixels.

Years are visited in a per-AOI shuffled order and the search stops as soon as
enough usable acquisitions have been seen. Visiting them in calendar order would
make the corpus's scene years pile up on 2018 wherever the first year already
offers a clear scene, which quietly correlates acquisition year with cloud
climatology.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
import math
from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np
from tools.aeronet_validation.global_scene_selection import (
    ELIGIBLE_YEARS,
    MAX_SCENE_CLOUD_COVER,
    SceneOption,
    assign_target_months,
    choose_scene,
    matchup_id,
    season_balance,
)

if TYPE_CHECKING:
    from collections.abc import Sequence

COLLECTION = "COPERNICUS/S2_HARMONIZED"

#: Half-width of the AOI, matching ``--template-size 128`` at 60 m in
#: ``capture_mixed_resolution_l1c``. Search only has to intersect; full coverage
#: is enforced later, at capture.
AOI_HALF_WIDTH_M = 128 * 60.0 / 2.0

#: Days either side of the target month's midpoint to search in each year.
WINDOW_HALF_DAYS = 20

#: Stop once this many usable acquisitions have been found; more years only add
#: choice that ``choose_scene`` rarely needs, at one Earth Engine listing each.
MIN_OPTIONS = 6

#: S2C joined the constellation in 2025 and carries the same L1C product naming.
SUPPORTED_SPACECRAFT = ("Sentinel-2A", "Sentinel-2B", "Sentinel-2C")


def aoi_bbox(longitude: float, latitude: float, *, half_width_m: float = AOI_HALF_WIDTH_M):
    """WGS84 bbox of the AOI footprint around a catalogue point."""

    latitude = float(latitude)
    # A projected coordinate reaching here would be silently accepted by
    # Earth Engine as a geometry spanning the planet many times over, which it
    # answers with "User memory limit exceeded" after a long wait rather than a
    # useful error. Refuse it here instead.
    if not -180.0 <= float(longitude) <= 180.0 or not -90.0 <= latitude <= 90.0:
        raise ValueError(
            f"({longitude}, {latitude}) is not a WGS84 coordinate; "
            "the catalogue may hold projected metres"
        )
    degrees_lat = half_width_m / 111_320.0
    scale = max(math.cos(math.radians(latitude)), 1e-6)
    degrees_lon = half_width_m / (111_320.0 * scale)
    return (
        float(longitude) - degrees_lon,
        latitude - degrees_lat,
        float(longitude) + degrees_lon,
        latitude + degrees_lat,
    )


def year_order(sample_id: str, years: Sequence[int]) -> tuple[int, ...]:
    """Deterministic per-AOI shuffle of the eligible years."""

    generator = np.random.default_rng(abs(hash(sample_id)) % (2**32))
    ordered = list(years)
    generator.shuffle(ordered)
    return tuple(int(year) for year in ordered)


def month_window(year: int, month: int, *, half_days: int = WINDOW_HALF_DAYS):
    midpoint = dt.date(int(year), int(month), 15)
    delta = dt.timedelta(days=int(half_days))
    return (midpoint - delta).isoformat(), (midpoint + delta).isoformat()


def options_from_images(images: Sequence[Any]) -> tuple[SceneOption, ...]:
    """Turn Earth Engine image records into selection candidates."""

    options = []
    for image in images:
        properties = dict(getattr(image, "properties", {}) or {})
        product_id = str(properties.get("PRODUCT_ID") or "")
        tile = str(properties.get("MGRS_TILE") or "")
        cloud = properties.get("CLOUDY_PIXEL_PERCENTAGE")
        spacecraft = str(properties.get("SPACECRAFT_NAME") or "")
        if not product_id or not tile or cloud is None:
            continue
        if spacecraft and spacecraft not in SUPPORTED_SPACECRAFT:
            continue
        options.append(
            SceneOption(
                product_id=product_id,
                datetime=image.acquisition_time_utc.isoformat(),
                mgrs_tile=tile if tile.startswith("T") else f"T{tile}",
                cloud_cover=float(cloud),
            )
        )
    return tuple(options)


def search_options(
    longitude: float,
    latitude: float,
    sample_id: str,
    target_month: int,
    *,
    years: Sequence[int] = ELIGIBLE_YEARS,
    min_options: int = MIN_OPTIONS,
    max_cloud: float = MAX_SCENE_CLOUD_COVER,
) -> tuple[tuple[SceneOption, ...], int]:
    """Collect candidate acquisitions, stopping once there are enough."""

    from edown import AOI, SearchConfig
    from edown.discovery import search_images
    from edown.errors import DiscoveryError

    aoi = AOI.from_bbox(aoi_bbox(longitude, latitude))
    collected: list[SceneOption] = []
    searches = 0
    for year in year_order(sample_id, years):
        start, end = month_window(year, target_month)
        searches += 1
        try:
            result = search_images(
                SearchConfig(
                    collection_id=COLLECTION,
                    start_date=start,
                    end_date=end,
                    aoi=aoi,
                    bands=("B2",),
                )
            )
        except DiscoveryError:
            continue  # No acquisitions that month; a polar winter, most likely.
        collected.extend(options_from_images(result.images))
        usable = sum(1 for option in collected if option.cloud_cover <= float(max_cloud))
        if usable >= int(min_options):
            break
    return tuple(collected), searches


def read_catalog(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def run(args: argparse.Namespace) -> dict[str, Any]:
    rows = read_catalog(Path(args.catalog))
    if args.shard_count > 1:
        rows = [row for index, row in enumerate(rows) if index % args.shard_count == args.shard]
    if args.limit:
        rows = rows[: int(args.limit)]
    months = assign_target_months([row["sample_id"] for row in rows], seed=args.seed)

    selected: dict[str, SceneOption] = {}
    records: list[dict[str, Any]] = []
    unresolved: list[str] = []
    searches = 0
    for position, row in enumerate(rows):
        sample_id = row["sample_id"]
        target = months[sample_id]
        options, used = search_options(
            float(row["longitude"]),
            float(row["latitude"]),
            sample_id,
            target,
            years=tuple(args.eligible_years),
            min_options=args.min_options,
        )
        searches += used
        scene = choose_scene(options, target_month=target)
        if scene is None:
            unresolved.append(sample_id)
            continue
        selected[sample_id] = scene
        records.append(
            {
                "sample_id": sample_id,
                "matchup_id": matchup_id(sample_id, scene),
                "product_id": scene.product_id,
                "datetime": scene.datetime,
                "mgrs_tile": scene.mgrs_tile,
                "cloud_cover": scene.cloud_cover,
                "target_month": target,
                "longitude": row["longitude"],
                "latitude": row["latitude"],
                "land_cover": row.get("land_cover", ""),
                "continent": row.get("continent", ""),
            }
        )
        if (position + 1) % 25 == 0:
            print(f"  {position + 1}/{len(rows)} AOIs, {len(records)} scenes", flush=True)

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    if records:
        with output.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(records[0]))
            writer.writeheader()
            writer.writerows(records)
    summary = {
        "aois": len(rows),
        "scenes": len(records),
        "unresolved": len(unresolved),
        "unresolved_ids": unresolved[:20],
        "earth_engine_searches": searches,
        "searches_per_aoi": searches / len(rows) if rows else 0.0,
        "season_balance": season_balance(selected),
        "catalog": str(output),
    }
    output.with_suffix(".summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return summary


def parser() -> argparse.ArgumentParser:
    value = argparse.ArgumentParser(description=__doc__)
    value.add_argument("--catalog", required=True, help="AOI catalogue CSV.")
    value.add_argument("--output", required=True)
    value.add_argument("--eligible-years", type=int, nargs="+", default=list(ELIGIBLE_YEARS))
    value.add_argument("--min-options", type=int, default=MIN_OPTIONS)
    value.add_argument("--seed", type=int, default=20260901)
    value.add_argument("--limit", type=int, default=0)
    value.add_argument("--shard", type=int, default=0)
    value.add_argument("--shard-count", type=int, default=1)
    return value


def main() -> None:
    print(json.dumps(run(parser().parse_args()), indent=2))


if __name__ == "__main__":
    main()
