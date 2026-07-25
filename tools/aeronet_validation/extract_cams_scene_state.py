"""Extract per-scene CAMS state (real ozone, AOD, water vapour) for the exact pairs.

Sen2Cor takes its ozone column from ECMWF per scene; the Van Heuklon climatology
proxy was falsified, so this pulls the actual CAMS analysis value at each pair
scene's site and local overpass hour from the daily global archive. The output
table keys on (matchup_id, day) and is consumed by the nonlinear trainer.
"""

from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path

import numpy as np
import xarray as xr
from tools.aeronet_validation.train_l2a_l1c_harmonizer import _case_ids, load_pairs

ROOT = Path("/gws/ssde/j25a/nceo_isp/public/siac_refactor")
CAMS_DIR = Path("/gws/ssde/j25b/nceo_ard/public/cams")
DEFAULT_PAIRS = ROOT / "analysis/l2a_l1c_physical_pairs_lowcloud152_20260716"
DEFAULT_CASES = ROOT / "reports/aod-final-performance-dashboard-20260713/data/all-cases.csv"
DEFAULT_MATCHUPS = ROOT / "matchups/matchups.csv"
DEFAULT_OUTPUT = DEFAULT_PAIRS / "cams_scene_state.csv"
KG_M2_PER_DU = 2.1415e-5


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pairs", type=Path, default=DEFAULT_PAIRS)
    parser.add_argument("--cases", type=Path, default=DEFAULT_CASES)
    parser.add_argument("--matchups", type=Path, default=DEFAULT_MATCHUPS)
    parser.add_argument("--cams", type=Path, default=CAMS_DIR)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    coordinates = {
        row["matchup_id"]: (float(row["latitude"]), float(row["longitude"]))
        for row in csv.DictReader(args.matchups.open())
    }
    dataset = load_pairs(
        args.pairs,
        _case_ids(args.cases),
        scene_day_max="2023-12-31",
        max_samples_per_scene=10,
        allow_missing_matchups=True,
    )
    requests: dict[str, list[tuple[str, str, float, float]]] = defaultdict(list)
    seen: set[tuple[str, str]] = set()
    for scene in dataset.scene_metadata:
        matchup_id = str(scene["matchup_id"])
        day = str(scene["day"])
        key = (matchup_id, day)
        if key in seen:
            continue
        seen.add(key)
        lat, lon = coordinates[matchup_id]
        requests[day].append((matchup_id, day, lat, lon))

    rows: list[dict[str, object]] = []
    missing_days: list[str] = []
    for day in sorted(requests):
        path = args.cams / f"{day}.nc"
        if not path.exists():
            missing_days.append(day)
            continue
        with xr.open_dataset(path) as data:
            # Two archive eras: newer files carry (forecast_period,
            # forecast_reference_time, ...), older ones a plain hourly time axis.
            field = data[["gtco3", "aod550", "tcwv"]]
            if "forecast_reference_time" in field.dims:
                field = field.squeeze("forecast_reference_time", drop=True)
            for matchup_id, scene_day, lat, lon in requests[day]:
                # Sentinel-2 descending node crosses ~10:30 local solar time.
                hour = int(round(10.5 - lon / 15.0)) % 24
                if "forecast_period" in field.dims:
                    time_selector = {"forecast_period": float(hour)}
                else:
                    time_selector = {"time": np.datetime64(f"{day}T{hour:02d}:00")}
                point = field.sel(
                    latitude=lat,
                    longitude=lon,
                    **time_selector,
                    method="nearest",
                )
                rows.append(
                    {
                        "matchup_id": matchup_id,
                        "day": scene_day,
                        "utc_hour": hour,
                        "cams_tco3_du": float(point["gtco3"].item()) / KG_M2_PER_DU,
                        "cams_aod550": float(point["aod550"].item()),
                        "cams_tcwv_cm": float(point["tcwv"].item()) / 10.0,
                    }
                )
        if len(rows) % 500 < len(requests[day]):
            print(f"processed through {day} ({len(rows)} scene states)", flush=True)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    values = np.asarray([row["cams_tco3_du"] for row in rows], dtype=np.float64)
    print(
        json.dumps(
            {
                "output": str(args.output),
                "scene_states": len(rows),
                "days": len(requests),
                "missing_days": missing_days,
                "tco3_du": {
                    "p05": float(np.percentile(values, 5)),
                    "median": float(np.median(values)),
                    "p95": float(np.percentile(values, 95)),
                },
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
