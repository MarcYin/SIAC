"""Collect lightweight GEE model AOD baselines for AERONET campaign sites.

This reads only coarse model AOD fields around the matchup time/AOI:

* ECMWF/CAMS/NRT ``total_aerosol_optical_depth_at_550nm_surface``
* NASA/GSFC/MERRA/aer/2 ``TOTEXTTAU``
* their simple mean when both exist

Each source is written as Phase-D-like validation JSON records under separate
directories so ``summarize_saved_results.py`` can score them directly. Run via
SLURM for the 250-site campaign.
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
HALF = float(os.environ.get("GEE_MODEL_AOD_HALF_DEG", "0.06"))
CAMS_OUT = ROOT / os.environ.get("GEE_CAMS_AOD_OUTDIR", "gee_cams_aod_campaign250")
MERRA_OUT = ROOT / os.environ.get("GEE_MERRA_AOD_OUTDIR", "gee_merra_aod_campaign250")
MEAN_OUT = ROOT / os.environ.get("GEE_MODEL_MEAN_AOD_OUTDIR", "gee_model_mean_aod_campaign250")
CAMS_BAND = "total_aerosol_optical_depth_at_550nm_surface"
MERRA_BAND = "TOTEXTTAU"


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


def _within_ee(retrieved: float, truth: float) -> bool:
    return abs(retrieved - truth) <= 0.05 + 0.15 * truth


def _model_values(lon: float, lat: float, sensing_time: str) -> dict[str, float | None]:
    obs_time = dt.datetime.fromisoformat(sensing_time.replace("Z", "+00:00"))
    geom = ee.Geometry.Rectangle([lon - HALF, lat - HALF, lon + HALF, lat + HALF])
    t = ee.Date(obs_time.isoformat())
    cams = (
        ee.ImageCollection("ECMWF/CAMS/NRT")
        .select(CAMS_BAND)
        .filterDate(t.advance(-3, "hour"), t.advance(3, "hour"))
        .mean()
    )
    merra = (
        ee.ImageCollection("NASA/GSFC/MERRA/aer/2")
        .select(MERRA_BAND)
        .filterDate(t.advance(-1, "hour"), t.advance(1, "hour"))
        .mean()
    )
    feature = ee.Dictionary(
        {
            "cams": cams.reduceRegion(
                ee.Reducer.mean(), geom, 40000, bestEffort=True, maxPixels=1_000_000
            ).get(CAMS_BAND),
            "merra": merra.reduceRegion(
                ee.Reducer.mean(), geom, 50000, bestEffort=True, maxPixels=1_000_000
            ).get(MERRA_BAND),
        }
    ).getInfo()
    cams_value = _safe_float(feature.get("cams"))
    merra_value = _safe_float(feature.get("merra"))
    mean_value = None
    if cams_value is not None and merra_value is not None:
        mean_value = (cams_value + merra_value) / 2.0
    return {"cams": cams_value, "merra": merra_value, "model_mean": mean_value}


def _record(
    matchup_id: str,
    row: dict[str, str],
    source: str,
    value: float | None,
    runtime_s: float,
    error: str | None = None,
) -> dict[str, object]:
    truth = float(row["aeronet_aod550_mean"])
    base: dict[str, object] = {
        "matchup_id": matchup_id,
        "site": row.get("site") or matchup_id.split("__")[0],
        "source": source,
        "truth": truth,
        "runtime_s": runtime_s,
    }
    if value is None:
        base.update(
            {
                "status": "FAILED",
                "error_type": "MissingAOD",
                "reason": error or f"{source} AOD is null",
            }
        )
        return base
    base.update(
        {
            "status": "OK",
            "retrieved": value,
            "err": value - truth,
            "within_ee": _within_ee(value, truth),
        }
    )
    return base


def collect_one(matchup_id: str, rows: dict[str, dict[str, str]]) -> None:
    row = rows[matchup_id]
    lon = float(row["longitude"])
    lat = float(row["latitude"])
    t0 = time.monotonic()
    values: dict[str, float | None]
    error = None
    try:
        values = _model_values(lon, lat, row["sensing_time_utc"])
    except Exception as exc:  # noqa: BLE001 - per-site array resilience
        values = {"cams": None, "merra": None, "model_mean": None}
        error = f"{type(exc).__name__}: {exc}"
    runtime_s = time.monotonic() - t0

    outputs = {
        CAMS_OUT: ("gee_cams", values["cams"]),
        MERRA_OUT: ("gee_merra", values["merra"]),
        MEAN_OUT: ("gee_cams_merra_mean", values["model_mean"]),
    }
    for out_dir, (source, value) in outputs.items():
        out_dir.mkdir(parents=True, exist_ok=True)
        out_path = out_dir / f"{matchup_id}.json"
        if out_path.exists() and not os.environ.get("FORCE"):
            continue
        record = _record(matchup_id, row, source, value, runtime_s, error)
        out_path.write_text(json.dumps(record, indent=2), encoding="utf-8")
    print(
        "GEE_MODEL_AOD_DONE "
        + json.dumps({"matchup_id": matchup_id, **values, "runtime_s": runtime_s}),
        flush=True,
    )


def main(argv: list[str]) -> int:
    if not argv:
        raise SystemExit("Usage: collect_gee_model_aod.py <matchup_id> [<matchup_id> ...]")
    _init_ee()
    rows = {row["matchup_id"]: row for row in csv.DictReader(MATCHUPS.open())}
    for matchup_id in argv:
        collect_one(matchup_id, rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
