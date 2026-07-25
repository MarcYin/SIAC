"""Collect same-scene Sentinel-2 L2A/Sen2Cor AOT for campaign diagnostics.

This is a lightweight independent AOD baseline/fallback candidate. It reads the
GEE COPERNICUS/S2_SR_HARMONIZED ``AOT`` band for the matchup day/tile over a
small station-centred buffer and writes Phase-D-like validation JSON records.
Run through SLURM; do not call this for the 250-site campaign on a login node.
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
OUT = ROOT / os.environ.get("S2_L2A_AOT_OUTDIR", "s2_l2a_aot_campaign250")
BUFFER_M = float(os.environ.get("S2_L2A_AOT_BUFFER_M", "3000"))


def _init_ee() -> None:
    service_account = os.environ.get("GEE_SERVICE_ACCOUNT", "python-gee@gee-marc.iam.gserviceaccount.com")
    key = os.environ.get("GEE_SERVICE_ACCOUNT_KEY", "/home/users/marcyin/gee-service-account.json")
    ee.Initialize(ee.ServiceAccountCredentials(service_account, key))


def _ee_safe_number(value: object) -> float | None:
    if value is None:
        return None
    try:
        out = float(value)
    except (TypeError, ValueError):
        return None
    return out if math.isfinite(out) else None


def _normalise_s2_aot(raw: float | None) -> tuple[float | None, str | None]:
    if raw is None:
        return None, None
    # GEE S2 SR AOT is commonly delivered as scaled integer. Keep already
    # physical values, but convert DN-like values such as 100 -> 0.1.
    if abs(raw) > 5.0:
        return raw / 1000.0, "dn_div_1000"
    return raw, "physical"


def _within_ee(retrieved: float, truth: float) -> bool:
    return abs(retrieved - truth) <= 0.05 + 0.15 * truth


def collect_one(matchup_id: str, rows: dict[str, dict[str, str]]) -> dict[str, object]:
    row = rows[matchup_id]
    site = row.get("site") or matchup_id.split("__")[0]
    lon = float(row["longitude"])
    lat = float(row["latitude"])
    truth = float(row["aeronet_aod550_mean"])
    sensing_time = row["sensing_time_utc"]
    day = sensing_time[:10]
    tile = matchup_id.split("__", 1)[1].split("_", 1)[0].removeprefix("T")
    next_day = (dt.date.fromisoformat(day) + dt.timedelta(days=1)).isoformat()
    geom = ee.Geometry.Point([lon, lat]).buffer(BUFFER_M)

    col = (
        ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
        .filterBounds(geom)
        .filterDate(day, next_day)
        .filter(ee.Filter.eq("MGRS_TILE", tile))
    )
    count = int(col.size().getInfo())
    if count == 0:
        # Tile-edge fallback: same day and AOI, no MGRS filter.
        col = (
            ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
            .filterBounds(geom)
            .filterDate(day, next_day)
        )
        count = int(col.size().getInfo())
    if count == 0:
        raise RuntimeError(f"no S2_SR_HARMONIZED image for {matchup_id} on {day}")

    image = ee.Image(col.sort("system:time_start").first())
    raw = image.select("AOT").reduceRegion(
        ee.Reducer.median(),
        geom,
        60,
        bestEffort=True,
        maxPixels=1_000_000,
    ).get("AOT").getInfo()
    raw_aot = _ee_safe_number(raw)
    aot, scale = _normalise_s2_aot(raw_aot)
    if aot is None:
        raise RuntimeError(f"S2_SR_HARMONIZED AOT is null for {matchup_id} on {day}")

    return {
        "matchup_id": matchup_id,
        "site": site,
        "status": "OK",
        "source": "COPERNICUS/S2_SR_HARMONIZED.AOT",
        "day": day,
        "tile": tile,
        "n_images": count,
        "aot_raw": raw_aot,
        "aot_scale": scale,
        "truth": truth,
        "retrieved": aot,
        "err": aot - truth,
        "within_ee": _within_ee(aot, truth),
        "buffer_m": BUFFER_M,
    }


def main(argv: list[str]) -> int:
    if not argv:
        raise SystemExit("Usage: collect_s2_l2a_aot.py <matchup_id> [<matchup_id> ...]")
    _init_ee()
    rows = {row["matchup_id"]: row for row in csv.DictReader(MATCHUPS.open())}
    OUT.mkdir(parents=True, exist_ok=True)
    status = 0
    for matchup_id in argv:
        out_path = OUT / f"{matchup_id}.json"
        if out_path.exists() and not os.environ.get("FORCE"):
            print(f"S2AOT {matchup_id} SKIP existing", flush=True)
            continue
        t0 = time.monotonic()
        try:
            record = collect_one(matchup_id, rows)
        except Exception as exc:  # noqa: BLE001 - per-site resilience for arrays
            row = rows.get(matchup_id, {})
            record = {
                "matchup_id": matchup_id,
                "site": row.get("site") or matchup_id.split("__")[0],
                "status": "FAILED",
                "truth": float(row["aeronet_aod550_mean"]) if row.get("aeronet_aod550_mean") else None,
                "error_type": type(exc).__name__,
                "reason": str(exc)[:500],
            }
            status = 1
        record["runtime_s"] = time.monotonic() - t0
        out_path.write_text(json.dumps(record, indent=2), encoding="utf-8")
        print(
            "S2AOT_DONE "
            + json.dumps(
                {
                    key: record.get(key)
                    for key in ("matchup_id", "status", "truth", "retrieved", "err", "within_ee")
                }
            ),
            flush=True,
        )
    return status


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
