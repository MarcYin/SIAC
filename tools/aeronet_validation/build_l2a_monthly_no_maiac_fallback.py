"""Build L2A monthly composites when MAIAC day selection has no valid days.

The normal Phase-D monthly L2A builder gates candidate days with MAIAC AOD from
the L1C sidecar metadata. A few high-latitude campaign sites have valid scene
metadata files but zero MAIAC-valid selected days. For those coverage holes, use
the same monthly Sen2Cor L2A composite path without an AOD day gate.
"""

from __future__ import annotations

import csv
import glob
import os
import re
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path("/gws/ssde/j25a/nceo_isp/public/siac_refactor")
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import bestpixel as bp  # noqa: E402
from build_l2a_monthly import BANDS, CACHE, ENDPOINT, HALF, OUT, TOP_K, slab_from_period  # noqa: E402
from build_l2a_seasonal import make_template  # noqa: E402

META = ROOT / "l1c_test" / "seasonal"


def build_one(matchup_id: str, rows: dict[str, dict[str, str]]) -> None:
    out_path = OUT / f"{matchup_id}.npz"
    if out_path.exists() and not os.environ.get("FORCE"):
        print(f"L2AmonFallback {matchup_id.split('__')[0]:18s} SKIP existing", flush=True)
        return

    row = rows[matchup_id]
    lon = float(row["longitude"])
    lat = float(row["latitude"])
    sensing_time = row["sensing_time_utc"]
    scene_year = int(sensing_time[:4])
    month0 = int(sensing_time[5:7])
    bbox = (lon - HALF, lat - HALF, lon + HALF, lat + HALF)
    tmpl, epsg, transform = make_template(lon, lat)

    metas = sorted(glob.glob(str(META / matchup_id / "*_meta.json")))
    win_order: list[tuple[int, int]] = []
    for meta_path in metas:
        match = re.search(r"(\d{4})-(\d{2})_meta", meta_path)
        if match:
            win_order.append((int(match.group(1)), int(match.group(2))))
    if not win_order:
        print(f"L2AmonFallback {matchup_id.split('__')[0]:18s} NO META", flush=True)
        return

    years = sorted({year for year, _ in win_order})
    months = sorted({month for _, month in win_order})
    t0 = time.perf_counter()
    periods = bp.build_monthly_composites(
        bbox,
        years,
        months,
        resolution=60.0,
        top_k=TOP_K,
        endpoint=ENDPOINT,
        bands=BANDS,
        disk_cache=CACHE,
        max_cloud_cover=90.0,
        output_crs="utm",
    )
    build_s = time.perf_counter() - t0
    print(
        f"BUILD_TIME fallback mid={matchup_id} seconds={build_s:.1f} periods={len(periods)}",
        flush=True,
    )

    slabs = []
    used = []
    for period in periods:
        try:
            slab = slab_from_period(period, tmpl)
        except Exception as exc:  # noqa: BLE001 - per-period resilience mirrors Phase-D builder
            print(f"    {period['year']}-{period['month']:02d} EXC {exc!r}", flush=True)
            continue
        if np.isfinite(slab[1]).mean() > 0.2:
            slabs.append(slab.astype(np.float32))
            used.append(f"{period['year']}-{period['month']:02d}")

    OUT.mkdir(parents=True, exist_ok=True)
    if not slabs:
        print(f"L2AmonFallback {matchup_id.split('__')[0]:18s} FAILED", flush=True)
        return
    comp = np.stack(slabs, 0)
    np.savez(
        out_path,
        comp=comp,
        epsg=epsg,
        transform=np.array(transform, float),
        realizations=np.array(used),
        scene_year=scene_year,
        month=month0,
        fallback=np.array("no_maiac_day_gate"),
    )
    print(
        f"L2AmonFallback {matchup_id.split('__')[0]:18s} n={len(used)} "
        f"grid={comp.shape[2]}x{comp.shape[3]} build={build_s:.1f}s -> {out_path.name}",
        flush=True,
    )


def main(argv: list[str]) -> int:
    if not argv:
        raise SystemExit("Usage: build_l2a_monthly_no_maiac_fallback.py <matchup_id> [...]")
    rows = {row["matchup_id"]: row for row in csv.DictReader((ROOT / "matchups" / "matchups.csv").open())}
    for matchup_id in argv:
        build_one(matchup_id, rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
