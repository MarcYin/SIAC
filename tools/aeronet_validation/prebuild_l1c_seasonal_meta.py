"""Prebuild L1C seasonal metadata needed by L2A monthly composites.

This runs only stage 1 from the Phase-D L1C seasonal path:
earth-search L1C enumeration plus GEE MAIAC/WVP/geometry sidecar metadata.
It intentionally does not run custom AC or build L1C composites.
"""

from __future__ import annotations

import calendar
import csv
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path("/gws/ssde/j25a/nceo_isp/public/siac_refactor")
MATCHUPS = ROOT / "matchups" / "matchups.csv"
WORK = ROOT / "l1c_test" / "seasonal"
HALF = 0.06
MONTH_OFFSETS = (-1, 0, 1)
POOL_YEARS = range(2020, 2026)


def _windows(scene_year: int, month: int):
    for year in POOL_YEARS:
        if year == scene_year:
            continue
        for offset in MONTH_OFFSETS:
            mm = month + offset
            yy = year
            if mm < 1:
                mm += 12
                yy -= 1
            elif mm > 12:
                mm -= 12
                yy += 1
            if yy == scene_year:
                continue
            yield yy, mm


def _run_meta_builder(bbox: tuple[float, float, float, float], year: int, month: int, base: Path) -> bool:
    meta_path = Path(f"{base}_meta.json")
    if meta_path.exists() and not os.environ.get("FORCE"):
        print(f"    {year}-{month:02d} SKIP meta exists", flush=True)
        return True

    last = calendar.monthrange(year, month)[1]
    start = f"{year}-{month:02d}-01"
    end = f"{year}-{month:02d}-{last:02d}"
    cmd = [
        sys.executable,
        str(ROOT / "build_atmo_sidecar_robust.py"),
        "--bbox",
        *[str(value) for value in bbox],
        "--month",
        start,
        end,
        "--out",
        str(base),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True, env=os.environ.copy(), check=False)
    if result.stdout:
        print(result.stdout, end="" if result.stdout.endswith("\n") else "\n", flush=True)
    if result.stderr:
        print(result.stderr[-2000:], file=sys.stderr, flush=True)
    ok = meta_path.exists()
    if not ok:
        print(f"    {year}-{month:02d} FAILED meta stage rc={result.returncode}", flush=True)
    return ok


def prebuild_one(matchup_id: str, rows: dict[str, dict[str, str]]) -> int:
    row = rows[matchup_id]
    lon = float(row["longitude"])
    lat = float(row["latitude"])
    sensing_time = row["sensing_time_utc"]
    scene_year = int(sensing_time[:4])
    month = int(sensing_time[5:7])
    bbox = (lon - HALF, lat - HALF, lon + HALF, lat + HALF)
    out_dir = WORK / matchup_id
    out_dir.mkdir(parents=True, exist_ok=True)

    ok = 0
    total = 0
    print(f"META {matchup_id} bbox={bbox}", flush=True)
    for year, mm in _windows(scene_year, month):
        total += 1
        base = out_dir / f"{year}-{mm:02d}"
        if _run_meta_builder(bbox, year, mm, base):
            ok += 1
    print(f"META_DONE {matchup_id} ok={ok}/{total}", flush=True)
    return 0 if ok > 0 else 1


def main(argv: list[str]) -> int:
    if not argv:
        raise SystemExit("Usage: prebuild_l1c_seasonal_meta.py <matchup_id> [<matchup_id> ...]")
    rows = {row["matchup_id"]: row for row in csv.DictReader(MATCHUPS.open())}
    status = 0
    for matchup_id in argv:
        status = max(status, prebuild_one(matchup_id, rows))
    return status


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
