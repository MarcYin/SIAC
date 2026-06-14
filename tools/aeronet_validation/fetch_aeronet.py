"""Stage 1: download global AERONET V3 direct-sun AOD measurements.

Uses the AERONET web service (https://aeronet.gsfc.nasa.gov/print_web_data_help_v3.html)
per site over the experiment period. Raw responses are cached so the stage is
resumable; sites already fetched are skipped unless ``--refetch`` is passed.
"""

from __future__ import annotations

import argparse
import io
import logging
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from datetime import date

import pandas as pd
import requests
from tools.aeronet_validation.common import (
    ExperimentPaths,
    aod_550_from_angstrom,
    aod_550_from_channels,
)

logger = logging.getLogger("aeronet_validation.fetch")

SITE_LIST_URL = "https://aeronet.gsfc.nasa.gov/aeronet_locations_v3.txt"
WEB_SERVICE_URL = "https://aeronet.gsfc.nasa.gov/cgi-bin/print_web_data_v3"
MISSING_VALUE = -999.0
_CHANNEL_COLUMNS = {
    440.0: "AOD_440nm",
    500.0: "AOD_500nm",
    675.0: "AOD_675nm",
    870.0: "AOD_870nm",
}
_ANGSTROM_COLUMN = "440-870_Angstrom_Exponent"


@dataclass(frozen=True)
class FetchSettings:
    start: date
    end: date
    level: str  # "AOD20" or "AOD15"
    workers: int
    timeout_s: float
    retries: int


def fetch_site_list(paths: ExperimentPaths, *, timeout_s: float = 60.0) -> pd.DataFrame:
    """Download (or load cached) AERONET V3 site list."""
    if paths.sites_file.exists():
        return pd.read_csv(paths.sites_file)
    response = requests.get(SITE_LIST_URL, timeout=timeout_s)
    response.raise_for_status()
    # First line is a generation-date banner; the CSV header follows.
    sites = pd.read_csv(io.StringIO(response.text), skiprows=1)
    sites = sites.rename(
        columns={
            "Site_Name": "site",
            "Longitude(decimal_degrees)": "longitude",
            "Latitude(decimal_degrees)": "latitude",
            "Elevation(meters)": "elevation_m",
        }
    )
    sites = sites.dropna(subset=["site", "longitude", "latitude"])
    sites.to_csv(paths.sites_file, index=False)
    logger.info("Fetched AERONET site list: %d sites", len(sites))
    return sites


def _fetch_site_raw(
    site: str, settings: FetchSettings, raw_path, session: requests.Session
) -> str | None:
    params = {
        "site": site,
        "year": settings.start.year,
        "month": settings.start.month,
        "day": settings.start.day,
        "year2": settings.end.year,
        "month2": settings.end.month,
        "day2": settings.end.day,
        settings.level: 1,
        "AVG": 10,  # all points (no daily averaging)
        "if_no_html": 1,
    }
    last_error: Exception | None = None
    for attempt in range(settings.retries):
        try:
            response = session.get(WEB_SERVICE_URL, params=params, timeout=settings.timeout_s)
            response.raise_for_status()
            raw_path.write_text(response.text)
            return response.text
        except Exception as error:  # noqa: BLE001 - retried network fetch
            last_error = error
            time.sleep(2.0 * (attempt + 1))
    logger.warning("Failed to fetch %s after %d attempts: %s", site, settings.retries, last_error)
    return None


def parse_site_measurements(raw_text: str) -> pd.DataFrame | None:
    """Parse an AERONET web-service response into canonical measurement rows."""
    lines = raw_text.splitlines()
    header_index = next(
        (i for i, line in enumerate(lines) if line.startswith("AERONET_Site,")), None
    )
    if header_index is None:
        return None
    frame = pd.read_csv(io.StringIO("\n".join(lines[header_index:])))
    if frame.empty:
        return None
    frame = frame.replace(MISSING_VALUE, pd.NA)

    timestamp = pd.to_datetime(
        frame["Date(dd:mm:yyyy)"] + " " + frame["Time(hh:mm:ss)"],
        format="%d:%m:%Y %H:%M:%S",
        utc=True,
    ).dt.tz_localize(None)

    parsed = pd.DataFrame({"datetime_utc": timestamp})
    for wavelength, column in _CHANNEL_COLUMNS.items():
        parsed[f"aod_{int(wavelength)}"] = pd.to_numeric(frame.get(column), errors="coerce")
    parsed["angstrom_440_870"] = pd.to_numeric(frame.get(_ANGSTROM_COLUMN), errors="coerce")
    parsed["aod_550"] = [_row_aod_550(row) for row in parsed.itertuples(index=False)]
    parsed = parsed.dropna(subset=["aod_550"]).reset_index(drop=True)
    return parsed if not parsed.empty else None


def _row_aod_550(row) -> float | None:
    aod_500 = getattr(row, "aod_500", None)
    angstrom = getattr(row, "angstrom_440_870", None)
    if pd.notna(aod_500) and pd.notna(angstrom) and float(aod_500) > 0.0:
        return aod_550_from_angstrom(float(aod_500), 500.0, float(angstrom))
    channels = {
        wavelength: float(value)
        for wavelength in (440.0, 500.0, 675.0, 870.0)
        if pd.notna(value := getattr(row, f"aod_{int(wavelength)}", None))
    }
    return aod_550_from_channels(channels)


def fetch_all_sites(
    paths: ExperimentPaths,
    settings: FetchSettings,
    *,
    site_filter: list[str] | None = None,
    max_sites: int | None = None,
    refetch: bool = False,
) -> pd.DataFrame:
    """Fetch AOD for every site; return the index of sites with data in period."""
    paths.ensure()
    sites = fetch_site_list(paths)
    if site_filter:
        wanted = {name.lower() for name in site_filter}
        sites = sites[sites["site"].str.lower().isin(wanted)]
    if max_sites is not None:
        sites = sites.head(max_sites)
    logger.info(
        "Fetching %s AOD for %d sites, %s to %s",
        settings.level,
        len(sites),
        settings.start,
        settings.end,
    )

    records: list[dict] = []
    session = requests.Session()

    def _process(site_row) -> dict | None:
        site = str(site_row.site)
        raw_path = paths.aeronet_raw_dir / f"{site}.csv"
        parsed_path = paths.aeronet_parsed_dir / f"{site}.csv"
        if parsed_path.exists() and not refetch:
            parsed = pd.read_csv(parsed_path, parse_dates=["datetime_utc"])
        else:
            if raw_path.exists() and not refetch:
                raw_text = raw_path.read_text()
            else:
                raw_text = _fetch_site_raw(site, settings, raw_path, session)
            if raw_text is None:
                return None
            parsed = parse_site_measurements(raw_text)
            if parsed is None:
                return None
            parsed.to_csv(parsed_path, index=False)
        if parsed.empty:
            return None
        return {
            "site": site,
            "longitude": float(site_row.longitude),
            "latitude": float(site_row.latitude),
            "elevation_m": float(site_row.elevation_m),
            "n_measurements": int(len(parsed)),
            "first_measurement": parsed["datetime_utc"].min(),
            "last_measurement": parsed["datetime_utc"].max(),
        }

    with ThreadPoolExecutor(max_workers=settings.workers) as pool:
        futures = {pool.submit(_process, row): row.site for row in sites.itertuples(index=False)}
        for i, future in enumerate(as_completed(futures), start=1):
            record = future.result()
            if record is not None:
                records.append(record)
            if i % 100 == 0:
                logger.info(
                    "Progress: %d/%d sites checked, %d with data", i, len(futures), len(records)
                )

    index = (
        pd.DataFrame(records).sort_values("site")
        if records
        else pd.DataFrame(
            columns=[
                "site",
                "longitude",
                "latitude",
                "elevation_m",
                "n_measurements",
                "first_measurement",
                "last_measurement",
            ]
        )
    )
    index.to_csv(paths.sites_with_data_file, index=False)
    logger.info(
        "Done: %d/%d sites have %s data in period; index at %s",
        len(index),
        len(sites),
        settings.level,
        paths.sites_with_data_file,
    )
    return index


def add_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--start-date", type=date.fromisoformat, default=date(2024, 1, 1))
    parser.add_argument("--end-date", type=date.fromisoformat, default=date(2024, 12, 31))
    parser.add_argument(
        "--level",
        choices=("AOD20", "AOD15"),
        default="AOD20",
        help="AERONET quality level (AOD20 = quality assured, AOD15 = near-realtime).",
    )
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--timeout", type=float, default=120.0)
    parser.add_argument("--retries", type=int, default=3)
    parser.add_argument("--sites", nargs="*", default=None, help="Optional site-name filter.")
    parser.add_argument("--max-sites", type=int, default=None)
    parser.add_argument("--refetch", action="store_true")


def run(args: argparse.Namespace, paths: ExperimentPaths) -> None:
    settings = FetchSettings(
        start=args.start_date,
        end=args.end_date,
        level=args.level,
        workers=args.workers,
        timeout_s=args.timeout,
        retries=args.retries,
    )
    fetch_all_sites(
        paths,
        settings,
        site_filter=args.sites,
        max_sites=args.max_sites,
        refetch=args.refetch,
    )
