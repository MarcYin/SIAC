"""Stage 2: find Sentinel-2 acquisitions coincident with AERONET measurements.

For every AERONET site with data in the period, search the Sentinel-2 catalog
(anonymous CDSE OData by default) over a point bounding box, then pair each
acquisition with AERONET measurements inside the temporal window. Per-site
search responses are cached under ``matchups/search/`` so the stage resumes.
"""

from __future__ import annotations

import argparse
import logging
import time
from dataclasses import dataclass
from datetime import date, datetime, timedelta
from pathlib import Path

import pandas as pd
from tools.aeronet_validation.common import (
    ExperimentPaths,
    matchup_id_for,
    read_json,
    write_json,
)

logger = logging.getLogger("aeronet_validation.matchup")

#: Half-width of the point bbox handed to the catalog search, in degrees.
SEARCH_BBOX_HALF_WIDTH_DEG = 0.001


@dataclass(frozen=True)
class MatchupSettings:
    start: date
    end: date
    window_minutes: float
    max_cloud_cover: float
    search_backend: str
    min_aeronet_points: int


def _search_site_products(
    site: str,
    longitude: float,
    latitude: float,
    settings: MatchupSettings,
    cache_path: Path,
    *,
    retries: int = 3,
) -> list[dict]:
    if cache_path.exists():
        return read_json(cache_path)["products"]

    from siac.api import search_sentinel2

    eps = SEARCH_BBOX_HALF_WIDTH_DEG
    bbox = (longitude - eps, latitude - eps, longitude + eps, latitude + eps)
    last_error: Exception | None = None
    for attempt in range(retries):
        try:
            products = search_sentinel2(
                bbox=bbox,
                start_date=settings.start.isoformat(),
                end_date=settings.end.isoformat(),
                max_cloud_cover=settings.max_cloud_cover,
                backend=settings.search_backend,
            )
            payload = [
                {
                    "product_id": p.product_id,
                    "mgrs_tile": p.mgrs_tile,
                    "sensing_date": p.sensing_date.isoformat(),
                    "processing_baseline": p.processing_baseline,
                    "cloud_cover": p.cloud_cover,
                    "satellite": p.satellite,
                }
                for p in products
            ]
            write_json(cache_path, {"site": site, "products": payload})
            return payload
        except Exception as error:  # noqa: BLE001 - retried network search
            last_error = error
            time.sleep(5.0 * (attempt + 1))
    logger.warning("S2 search failed for %s after %d attempts: %s", site, retries, last_error)
    return []


def _dedupe_products(products: list[dict]) -> list[dict]:
    """One product per acquisition: drop tile duplicates and older baselines."""
    best: dict[str, dict] = {}
    for product in products:
        sensing_minute = product["sensing_date"][:16]
        current = best.get(sensing_minute)
        if current is None:
            best[sensing_minute] = product
            continue
        baseline = int(product["processing_baseline"].replace("N", "") or 0)
        current_baseline = int(current["processing_baseline"].replace("N", "") or 0)
        if (baseline, product["mgrs_tile"]) > (current_baseline, current["mgrs_tile"]):
            best[sensing_minute] = product
    return [best[key] for key in sorted(best)]


def _match_measurements(
    measurements: pd.DataFrame, sensing_time: datetime, window_minutes: float
) -> pd.DataFrame:
    window = timedelta(minutes=window_minutes)
    in_window = measurements[
        (measurements["datetime_utc"] >= sensing_time - window)
        & (measurements["datetime_utc"] <= sensing_time + window)
    ]
    return in_window


def build_matchups(paths: ExperimentPaths, settings: MatchupSettings) -> pd.DataFrame:
    paths.ensure()
    sites = pd.read_csv(paths.sites_with_data_file)
    if sites.empty:
        raise SystemExit("No sites with AERONET data; run the fetch-aeronet stage first.")
    logger.info("Building matchups for %d sites", len(sites))

    rows: list[dict] = []
    for i, site_row in enumerate(sites.itertuples(index=False), start=1):
        site = str(site_row.site)
        measurements = pd.read_csv(
            paths.aeronet_parsed_dir / f"{site}.csv", parse_dates=["datetime_utc"]
        )
        cache_path = paths.matchup_search_dir / f"{site}.json"
        products = _search_site_products(
            site,
            float(site_row.longitude),
            float(site_row.latitude),
            settings,
            cache_path,
        )
        for product in _dedupe_products(products):
            sensing_time = datetime.fromisoformat(product["sensing_date"]).replace(tzinfo=None)
            matched = _match_measurements(measurements, sensing_time, settings.window_minutes)
            if len(matched) < settings.min_aeronet_points:
                continue
            offsets_min = (matched["datetime_utc"] - sensing_time).dt.total_seconds().abs() / 60.0
            rows.append(
                {
                    "matchup_id": matchup_id_for(
                        site,
                        product["mgrs_tile"],
                        sensing_time.strftime("%Y%m%dT%H%M%S"),
                    ),
                    "site": site,
                    "longitude": float(site_row.longitude),
                    "latitude": float(site_row.latitude),
                    "elevation_m": float(site_row.elevation_m),
                    "product_id": product["product_id"],
                    "mgrs_tile": product["mgrs_tile"],
                    "satellite": product["satellite"],
                    "sensing_time_utc": sensing_time.isoformat(),
                    "scene_cloud_cover": float(product["cloud_cover"]),
                    "n_aeronet": int(len(matched)),
                    "mean_abs_time_offset_min": float(offsets_min.mean()),
                    "aeronet_aod550_mean": float(matched["aod_550"].mean()),
                    "aeronet_aod550_std": float(matched["aod_550"].std(ddof=0)),
                    "aeronet_aod440_mean": float(matched["aod_440"].mean()),
                    "aeronet_aod500_mean": float(matched["aod_500"].mean()),
                    "aeronet_angstrom_mean": float(matched["angstrom_440_870"].mean()),
                }
            )
        if i % 25 == 0:
            logger.info("Progress: %d/%d sites searched, %d matchups", i, len(sites), len(rows))

    matchups = pd.DataFrame(rows)
    if not matchups.empty:
        matchups = matchups.sort_values(["site", "sensing_time_utc"]).reset_index(drop=True)
    matchups.to_csv(paths.matchups_file, index=False)
    logger.info("Wrote %d matchups to %s", len(matchups), paths.matchups_file)
    return matchups


def add_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--start-date", type=date.fromisoformat, default=date(2024, 1, 1))
    parser.add_argument("--end-date", type=date.fromisoformat, default=date(2024, 12, 31))
    parser.add_argument(
        "--window-minutes",
        type=float,
        default=30.0,
        help="Half-width of the AERONET temporal matchup window around the S2 overpass.",
    )
    parser.add_argument("--max-cloud-cover", type=float, default=80.0)
    parser.add_argument("--search-backend", choices=("cdse", "gcs"), default="cdse")
    parser.add_argument(
        "--min-aeronet-points",
        type=int,
        default=1,
        help="Minimum AERONET measurements inside the window to accept a matchup.",
    )


def run(args: argparse.Namespace, paths: ExperimentPaths) -> None:
    settings = MatchupSettings(
        start=args.start_date,
        end=args.end_date,
        window_minutes=args.window_minutes,
        max_cloud_cover=args.max_cloud_cover,
        search_backend=args.search_backend,
        min_aeronet_points=args.min_aeronet_points,
    )
    build_matchups(paths, settings)
