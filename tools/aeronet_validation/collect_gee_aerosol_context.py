"""Collect coarse aerosol-composition context for campaign AOD diagnostics.

The scalar CAMS/MERRA AOD baselines show that independent aerosol products add
oracle coverage but need a trust signal. This collector stores composition and
spectral-AOD fields from the same coarse products so selectors can be tested
without using AERONET labels as features.

Run campaign collection via SLURM; one task per matchup.
"""

from __future__ import annotations

import csv
import datetime as dt
import json
import math
import os
import sys
import time
from pathlib import Path

import ee

ROOT = Path("/gws/ssde/j25a/nceo_isp/public/siac_refactor")
MATCHUPS = ROOT / "matchups" / "matchups.csv"
OUT = ROOT / os.environ.get("GEE_AEROSOL_CONTEXT_OUTDIR", "gee_aerosol_context_campaign250")
HALF = float(os.environ.get("GEE_AEROSOL_CONTEXT_HALF_DEG", "0.06"))

CAMS_BANDS = (
    "total_aerosol_optical_depth_at_550nm_surface",
    "sea_salt_aerosol_optical_depth_at_550nm_surface",
    "dust_aerosol_optical_depth_at_550nm_surface",
    "organic_matter_aerosol_optical_depth_at_550nm_surface",
    "black_carbon_aerosol_optical_depth_at_550nm_surface",
    "sulphate_aerosol_optical_depth_at_550nm_surface",
    "total_aerosol_optical_depth_at_469nm_surface",
    "total_aerosol_optical_depth_at_670nm_surface",
    "total_aerosol_optical_depth_at_865nm_surface",
)

MERRA_BANDS = (
    "TOTEXTTAU",
    "TOTSCATAU",
    "TOTANGSTR",
    "DUEXTTAU",
    "OCEXTTAU",
    "BCEXTTAU",
    "SUEXTTAU",
    "SSEXTTAU",
)


def _init_ee() -> None:
    service_account = os.environ.get("GEE_SERVICE_ACCOUNT", "python-gee@gee-marc.iam.gserviceaccount.com")
    key = os.environ.get("GEE_SERVICE_ACCOUNT_KEY", "/home/users/marcyin/gee-service-account.json")
    ee.Initialize(ee.ServiceAccountCredentials(service_account, key))


def _safe_float(value: object) -> float | None:
    if value is None:
        return None
    try:
        out = float(value)
    except (TypeError, ValueError):
        return None
    return out if math.isfinite(out) else None


def _with_ratios(prefix: str, values: dict[str, float | None], total_key: str) -> dict[str, float | None]:
    total = values.get(f"{prefix}_{total_key}")
    out = dict(values)
    if total is None or abs(total) < 1e-12:
        return out
    for key, value in values.items():
        if key.startswith(prefix) and key != f"{prefix}_{total_key}" and value is not None:
            out[f"{key}_frac"] = value / total
    return out


def _aerosol_context(lon: float, lat: float, sensing_time: str) -> dict[str, float | None]:
    obs_time = dt.datetime.fromisoformat(sensing_time.replace("Z", "+00:00"))
    geom = ee.Geometry.Rectangle([lon - HALF, lat - HALF, lon + HALF, lat + HALF])
    t = ee.Date(obs_time.isoformat())

    cams_img = (
        ee.ImageCollection("ECMWF/CAMS/NRT")
        .select(CAMS_BANDS)
        .filterDate(t.advance(-3, "hour"), t.advance(3, "hour"))
        .mean()
    )
    merra_img = (
        ee.ImageCollection("NASA/GSFC/MERRA/aer/2")
        .select(MERRA_BANDS)
        .filterDate(t.advance(-1, "hour"), t.advance(1, "hour"))
        .mean()
    )
    raw = ee.Dictionary(
        {
            "cams": cams_img.reduceRegion(
                ee.Reducer.mean(),
                geom,
                40000,
                bestEffort=True,
                maxPixels=1_000_000,
            ),
            "merra": merra_img.reduceRegion(
                ee.Reducer.mean(),
                geom,
                50000,
                bestEffort=True,
                maxPixels=1_000_000,
            ),
        }
    ).getInfo()

    out: dict[str, float | None] = {}
    for band in CAMS_BANDS:
        out[f"cams_{band}"] = _safe_float(raw.get("cams", {}).get(band))
    for band in MERRA_BANDS:
        out[f"merra_{band}"] = _safe_float(raw.get("merra", {}).get(band))

    out = _with_ratios("cams", out, "total_aerosol_optical_depth_at_550nm_surface")
    out = _with_ratios("merra", out, "TOTEXTTAU")
    merra_total = out.get("merra_TOTEXTTAU")
    merra_scatter = out.get("merra_TOTSCATAU")
    if merra_total is not None and merra_scatter is not None and abs(merra_total) > 1e-12:
        out["merra_scatter_fraction"] = merra_scatter / merra_total
    return out


def collect_one(matchup_id: str, rows: dict[str, dict[str, str]]) -> None:
    row = rows[matchup_id]
    t0 = time.monotonic()
    status = "OK"
    error = None
    try:
        values = _aerosol_context(
            float(row["longitude"]),
            float(row["latitude"]),
            row["sensing_time_utc"],
        )
    except Exception as exc:  # noqa: BLE001 - resilient array jobs
        values = {}
        status = "FAILED"
        error = f"{type(exc).__name__}: {exc}"
    record: dict[str, object] = {
        "matchup_id": matchup_id,
        "site": row.get("site") or matchup_id.split("__")[0],
        "status": status,
        "source": "gee_aerosol_context",
        "runtime_s": time.monotonic() - t0,
        "sensing_time_utc": row["sensing_time_utc"],
        "lon": float(row["longitude"]),
        "lat": float(row["latitude"]),
        "values": values,
    }
    if error:
        record["reason"] = error
    OUT.mkdir(parents=True, exist_ok=True)
    out_path = OUT / f"{matchup_id}.json"
    if not out_path.exists() or os.environ.get("FORCE"):
        out_path.write_text(json.dumps(record, indent=2, sort_keys=True), encoding="utf-8")
    print(
        "GEE_AEROSOL_CONTEXT_DONE "
        + json.dumps({"matchup_id": matchup_id, "status": status, "runtime_s": record["runtime_s"]}),
        flush=True,
    )


def main(argv: list[str]) -> int:
    if not argv:
        raise SystemExit("Usage: collect_gee_aerosol_context.py <matchup_id> [<matchup_id> ...]")
    _init_ee()
    rows = {row["matchup_id"]: row for row in csv.DictReader(MATCHUPS.open())}
    for matchup_id in argv:
        collect_one(matchup_id, rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
