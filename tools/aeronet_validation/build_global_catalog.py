#!/usr/bin/env python3
"""Build the global AOI catalogue: WorldCover tiles -> pool -> stratified draw."""

from __future__ import annotations

import argparse
import csv
import json
import subprocess
from dataclasses import asdict
from pathlib import Path

from tools.aeronet_validation.global_candidate_pool import (
    LANDCOVER_BAND,
    LANDCOVER_COLLECTION,
    candidates_from_class_grid,
    global_tiles,
    pool_composition,
)
from tools.aeronet_validation.global_catalog_sampler import (
    DEFAULT_TARGETS,
    composition,
    latitude_band,
    sample_catalog,
)

# MCD12Q1 is annual; the window selects one epoch.
LANDCOVER_WINDOW = ("2021-01-01", "2021-12-31")


def fetch_tile(
    bbox: tuple[float, float, float, float], cache: Path, edown: str, timeout_s: int = 1800
) -> None:
    cache.mkdir(parents=True, exist_ok=True)
    result = subprocess.run(
        [
            edown,
            "download",
            "--collection-id",
            LANDCOVER_COLLECTION,
            "--band",
            LANDCOVER_BAND,
            "--bbox",
            *[str(c) for c in bbox],
            "--start-date",
            LANDCOVER_WINDOW[0],
            "--end-date",
            LANDCOVER_WINDOW[1],
            "--output-root",
            str(cache),
        ],
        capture_output=True,
        text=True,
        timeout=timeout_s,
    )
    if result.returncode != 0:
        raise RuntimeError(f"edown failed for {bbox}: {result.stderr[-400:]}")


def read_tile(cache: Path):
    """Return the class grid and a resolver from pixel indices to lon/lat.

    The raster is delivered on the collection's native grid, which for MCD12Q1
    is sinusoidal. Longitude there is a function of both axes, so the resolver
    transforms points rather than returning separable coordinate vectors.
    """

    import rasterio
    from pyproj import Transformer

    rasters = sorted(cache.glob("images/*/*.tif"))
    if not rasters:
        return None
    with rasterio.open(rasters[-1]) as source:
        codes = source.read(1).astype(int)
        transform = source.transform
        to_wgs84 = Transformer.from_crs(source.crs, "EPSG:4326", always_xy=True)

    def to_lonlat(rows, columns):
        # Pixel centres through the affine, then out of the native projection.
        x = transform.c + transform.a * (columns + 0.5) + transform.b * (rows + 0.5)
        y = transform.f + transform.d * (columns + 0.5) + transform.e * (rows + 0.5)
        return to_wgs84.transform(x, y)

    return codes, to_lonlat


def check_geography(selected: list) -> None:
    """Refuse a catalogue whose coordinates are not plausible lon/lat.

    Native-grid rasters arrive in projected metres, and treating those as
    degrees produces a catalogue that looks complete -- right row count, right
    land-cover quotas -- while every point is nonsense and every latitude band
    collapses to one value. That failure is invisible in the row count, so it is
    asserted here instead.
    """

    longitudes = [float(candidate.longitude) for candidate in selected]
    latitudes = [float(candidate.latitude) for candidate in selected]
    if not all(-180.0 <= value <= 180.0 for value in longitudes):
        raise RuntimeError("catalogue longitudes are outside [-180, 180]; check the raster CRS")
    if not all(-90.0 <= value <= 90.0 for value in latitudes):
        raise RuntimeError("catalogue latitudes are outside [-90, 90]; check the raster CRS")
    bands = {latitude_band(value) for value in latitudes}
    if len(bands) < 2:
        raise RuntimeError(f"every catalogue point landed in one latitude band ({bands})")


def run(args: argparse.Namespace) -> dict:
    tiles = global_tiles(float(args.tile_degrees))
    pool = []
    fetched, empty = 0, 0
    for index, bbox in enumerate(tiles):
        cache = Path(args.cache) / f"tile_{index:04d}"
        if not args.skip_fetch:
            try:
                fetch_tile(bbox, cache, args.edown)
            except RuntimeError as error:
                print(f"tile {index} fetch failed, skipping: {error}", flush=True)
                continue
        payload = read_tile(cache)
        if payload is None:
            empty += 1
            continue
        fetched += 1
        codes, to_lonlat = payload
        pool.extend(
            candidates_from_class_grid(
                codes,
                to_lonlat=to_lonlat,
                tile_id=f"t{index:04d}",
                per_tile=int(args.per_tile),
                seed=index,
                within=bbox,
            )
        )
        if (index + 1) % 10 == 0:
            print(f"  {index + 1}/{len(tiles)} tiles, pool {len(pool)}", flush=True)

    if not pool:
        raise RuntimeError("WorldCover produced no candidates; check edown and the cache")
    selected = sample_catalog(pool, total=int(args.total), seed=int(args.seed))
    # Checked before writing, so a catalogue that fails never lands on disk to
    # be picked up by a downstream stage.
    check_geography(selected)
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(asdict(selected[0])))
        writer.writeheader()
        writer.writerows(asdict(c) for c in selected)
    summary = {
        "tiles_total": len(tiles),
        "tiles_with_data": fetched,
        "tiles_empty": empty,
        "pool_size": len(pool),
        "pool_composition": pool_composition(pool),
        "selected": len(selected),
        "targets": dict(DEFAULT_TARGETS),
        "composition": composition(selected),
        "catalog": str(output),
    }
    output.with_suffix(".summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return summary


def parser() -> argparse.ArgumentParser:
    value = argparse.ArgumentParser(description=__doc__)
    value.add_argument("--output", required=True)
    value.add_argument("--cache", required=True)
    value.add_argument("--total", type=int, default=5000)
    # 10 degrees at MCD12Q1's native 500 m is ~2,200 pixels square, well inside
    # edown's 48 MB request limit; 30 degrees would not be.
    value.add_argument("--tile-degrees", type=float, default=10.0)
    value.add_argument("--per-tile", type=int, default=400)
    value.add_argument("--seed", type=int, default=20260901)
    value.add_argument("--edown", default="/home/users/marcyin/SIAC/tools/edown_runtime")
    value.add_argument("--skip-fetch", action="store_true")
    return value


def main() -> None:
    print(json.dumps(run(parser().parse_args()), indent=2))


if __name__ == "__main__":
    main()
