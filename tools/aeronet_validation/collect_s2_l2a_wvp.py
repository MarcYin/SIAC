"""Collect same-scene Sentinel-2 L2A/Sen2Cor water vapour.

The selected L2A image is the closest acquisition on the matchup day and MGRS
tile. The station-centred median WVP is written in centimetres for use as one
fixed TCWV throughout a controlled 6S aerosol retrieval. Run campaign-scale
collection through Slurm, not on a login node.
"""

from __future__ import annotations

import csv
import datetime as dt
import json
import os
import sys
import time
from pathlib import Path

import ee
from tools.aeronet_validation.target_tcwv import (
    normalise_s2_l2a_wvp,
    product_acquisition_time,
)

ROOT = Path("/gws/ssde/j25a/nceo_isp/public/siac_refactor")
MATCHUPS = ROOT / "matchups" / "matchups.csv"
OUT = ROOT / os.environ.get("S2_L2A_WVP_OUTDIR", "s2_l2a_wvp_campaign250")
BUFFER_M = float(os.environ.get("S2_L2A_WVP_BUFFER_M", "3000"))


def _init_ee() -> None:
    service_account = os.environ.get(
        "GEE_SERVICE_ACCOUNT", "python-gee@gee-marc.iam.gserviceaccount.com"
    )
    key = os.environ.get("GEE_SERVICE_ACCOUNT_KEY", "/home/users/marcyin/gee-service-account.json")
    ee.Initialize(ee.ServiceAccountCredentials(service_account, key))


def collect_one(matchup_id: str, rows: dict[str, dict[str, str]]) -> dict[str, object]:
    row = rows[matchup_id]
    lon = float(row["longitude"])
    lat = float(row["latitude"])
    sensing_time = row["sensing_time_utc"]
    day = sensing_time[:10]
    next_day = (dt.date.fromisoformat(day) + dt.timedelta(days=1)).isoformat()
    tile = row.get("mgrs_tile") or matchup_id.split("__", 1)[1].split("_", 1)[0]
    tile = tile.removeprefix("T")
    geom = ee.Geometry.Point([lon, lat]).buffer(BUFFER_M)
    target_ms = ee.Date(sensing_time).millis()

    collection = (
        ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
        .filterBounds(geom)
        .filterDate(day, next_day)
        .filter(ee.Filter.eq("MGRS_TILE", tile))
    )
    count = int(collection.size().getInfo())
    tile_fallback = False
    if count == 0:
        collection = (
            ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
            .filterBounds(geom)
            .filterDate(day, next_day)
        )
        count = int(collection.size().getInfo())
        tile_fallback = True
    if count == 0:
        raise RuntimeError(f"no S2_SR_HARMONIZED image for {matchup_id} on {day}")

    def _with_time_delta(image: ee.Image) -> ee.Image:
        delta = ee.Number(image.get("system:time_start")).subtract(target_ms).abs()
        return image.set("_siac_time_delta_ms", delta)

    image = ee.Image(collection.map(_with_time_delta).sort("_siac_time_delta_ms").first())
    wvp = (
        image.select("WVP")
        .reduceRegion(
            ee.Reducer.median(),
            geom,
            60,
            bestEffort=True,
            maxPixels=1_000_000,
        )
        .get("WVP")
    )
    info = ee.Dictionary(
        {
            "wvp": wvp,
            "system_index": image.get("system:index"),
            "product_id": image.get("PRODUCT_ID"),
            "system_time_start": image.get("system:time_start"),
            "time_delta_ms": image.get("_siac_time_delta_ms"),
        }
    ).getInfo()
    tcwv_cm, scale = normalise_s2_l2a_wvp(info.get("wvp"))
    if tcwv_cm is None:
        raise RuntimeError(f"invalid S2 L2A WVP for {matchup_id}: {info.get('wvp')!r}")
    target_acquisition = product_acquisition_time(row.get("product_id"))
    selected_acquisition = product_acquisition_time(info.get("product_id"))
    if target_acquisition is None or selected_acquisition is None:
        raise RuntimeError(f"cannot verify L1C/L2A product acquisition time for {matchup_id}")
    acquisition_delta_s = abs((selected_acquisition - target_acquisition).total_seconds())
    if acquisition_delta_s > 60.0:
        raise RuntimeError(
            f"selected L2A acquisition differs from target by {acquisition_delta_s:.1f}s "
            f"for {matchup_id}"
        )

    return {
        "matchup_id": matchup_id,
        "site": row.get("site") or matchup_id.split("__")[0],
        "status": "OK",
        "source": "COPERNICUS/S2_SR_HARMONIZED.WVP",
        "target_product_id": row.get("product_id"),
        "target_sensing_time": sensing_time,
        "selected_product_id": info.get("product_id"),
        "selected_system_index": info.get("system_index"),
        "selected_system_time_start": info.get("system_time_start"),
        "acquisition_delta_s": acquisition_delta_s,
        "system_time_delta_s": float(info["time_delta_ms"]) / 1000.0,
        "day": day,
        "tile": tile,
        "tile_fallback": tile_fallback,
        "n_images": count,
        "wvp_raw": float(info["wvp"]),
        "wvp_scale": scale,
        "tcwv_cm": tcwv_cm,
        "buffer_m": BUFFER_M,
    }


def main(argv: list[str]) -> int:
    if not argv:
        raise SystemExit("Usage: collect_s2_l2a_wvp.py <matchup_id> [<matchup_id> ...]")
    _init_ee()
    rows = {row["matchup_id"]: row for row in csv.DictReader(MATCHUPS.open())}
    OUT.mkdir(parents=True, exist_ok=True)
    status = 0
    for matchup_id in argv:
        out_path = OUT / f"{matchup_id}.json"
        if out_path.exists() and not os.environ.get("FORCE"):
            try:
                if json.loads(out_path.read_text(encoding="utf-8")).get("status") == "OK":
                    print(f"S2WVP {matchup_id} SKIP existing OK", flush=True)
                    continue
            except (OSError, ValueError):
                pass
        started = time.monotonic()
        try:
            record = collect_one(matchup_id, rows)
        except Exception as exc:  # noqa: BLE001 - preserve per-site array progress
            row = rows.get(matchup_id, {})
            record = {
                "matchup_id": matchup_id,
                "site": row.get("site") or matchup_id.split("__")[0],
                "status": "FAILED",
                "error_type": type(exc).__name__,
                "reason": str(exc)[:500],
            }
            status = 1
        record["runtime_s"] = time.monotonic() - started
        out_path.write_text(json.dumps(record, indent=2), encoding="utf-8")
        print(
            "S2WVP_DONE "
            + json.dumps(
                {
                    key: record.get(key)
                    for key in (
                        "matchup_id",
                        "status",
                        "tcwv_cm",
                        "acquisition_delta_s",
                        "reason",
                    )
                }
            ),
            flush=True,
        )
    return status


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
