"""Collect NOAA VIIRS AOD EDR diagnostics for campaign sites.

The Earth Engine collection ``NOAA/VIIRS/AOD_EDR/V3`` provides daily 750 m AOD
at 550 nm plus land aerosol-model alternatives. This collector writes several
Phase-D-like validation JSON directories so each variant can be scored directly.
Run the 250-site campaign through SLURM.
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
HALF = float(os.environ.get("VIIRS_AOD_HALF_DEG", "0.06"))
COLLECTION = "NOAA/VIIRS/AOD_EDR/V3"
VARIANTS = {
    "viirs_aod550_qc0": ("AOD550", "qc0"),
    "viirs_aod550_qc2": ("AOD550", "qc2"),
    "viirs_dust_qc2": ("AOD550LndMdl_Dust", "qc2"),
    "viirs_generic_qc2": ("AOD550LndMdl_Generic", "qc2"),
    "viirs_urban_qc2": ("AOD550LndMdl_Urban", "qc2"),
    "viirs_smoke_qc2": ("AOD550LndMdl_Smoke", "qc2"),
}
OUT_SUFFIX = os.environ.get("VIIRS_AOD_OUTDIR_SUFFIX", "campaign250")


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
    if not math.isfinite(out) or out < 0.0 or out > 5.0:
        return None
    return out


def _within_ee(retrieved: float, truth: float) -> bool:
    return abs(retrieved - truth) <= 0.05 + 0.15 * truth


def _mask_image(image: ee.Image, band: str, qc_mode: str) -> ee.Image:
    value = image.select(band)
    qc = image.select("QCAll")
    mask = qc.eq(0) if qc_mode == "qc0" else qc.lte(2)
    return value.updateMask(mask).updateMask(value.gte(0)).updateMask(value.lte(5))


def _values(lon: float, lat: float, sensing_time: str) -> dict[str, float | None]:
    day = sensing_time[:10]
    next_day = (dt.date.fromisoformat(day) + dt.timedelta(days=1)).isoformat()
    geom = ee.Geometry.Rectangle([lon - HALF, lat - HALF, lon + HALF, lat + HALF])
    col = ee.ImageCollection(COLLECTION).filterBounds(geom).filterDate(day, next_day)
    count = int(col.size().getInfo())
    if count == 0:
        return {name: None for name in VARIANTS}
    out: dict[str, float | None] = {}
    for name, (band, qc_mode) in VARIANTS.items():
        image = col.map(lambda img, b=band, q=qc_mode: _mask_image(img, b, q)).median()
        raw = image.reduceRegion(
            ee.Reducer.median(),
            geom,
            750,
            bestEffort=True,
            maxPixels=1_000_000,
        ).get(band).getInfo()
        out[name] = _safe_float(raw)
    return out


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
    error = None
    try:
        values = _values(lon, lat, row["sensing_time_utc"])
    except Exception as exc:  # noqa: BLE001 - per-site array resilience
        values = {name: None for name in VARIANTS}
        error = f"{type(exc).__name__}: {exc}"
    runtime_s = time.monotonic() - t0
    for source, value in values.items():
        out_dir = ROOT / f"{source}_{OUT_SUFFIX}"
        out_dir.mkdir(parents=True, exist_ok=True)
        out_path = out_dir / f"{matchup_id}.json"
        if out_path.exists() and not os.environ.get("FORCE"):
            continue
        record = _record(matchup_id, row, source, value, runtime_s, error)
        out_path.write_text(json.dumps(record, indent=2), encoding="utf-8")
    print(
        "VIIRS_AOD_DONE "
        + json.dumps({"matchup_id": matchup_id, **values, "runtime_s": runtime_s}),
        flush=True,
    )


def main(argv: list[str]) -> int:
    if not argv:
        raise SystemExit("Usage: collect_viirs_aod.py <matchup_id> [<matchup_id> ...]")
    _init_ee()
    rows = {row["matchup_id"]: row for row in csv.DictReader(MATCHUPS.open())}
    for matchup_id in argv:
        collect_one(matchup_id, rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
