#!/usr/bin/env python3
"""Measure how the winner index changes when the ordering band changes.

Candidate days and daily scalars come from the reference index, so the only
thing that varies is which Cloud Score+ band feeds the per-pixel argmax. The
candidate rule is relative and therefore band-insensitive (measured: cs mean
coverage 0.751 versus cs_cdf 0.798, identical day set); the ordering is an
absolute argmax and is not.
"""

from __future__ import annotations

import json
from collections import defaultdict
from pathlib import Path

import numpy as np
from tools.aeronet_validation.build_cloudscore_index_local import (
    read_planes,
    scene_grid,
)
from tools.aeronet_validation.cloudscore_local_mosaic import (
    Candidate,
    agreement,
    compose_monthly_winners,
)

ROOT = Path("/gws/ssde/j25a/nceo_isp/public/siac_refactor")
CAMPAIGN = ROOT / "current_date_toa_surface_allera_stac_20260809"
SCENE = "AAU_ET__T37PDL_20201224T074229"


def main() -> None:
    reference = np.load(
        CAMPAIGN / f"historical_l1c_cloudscore_index20_ablation/{SCENE}.npz", allow_pickle=False
    )
    images = [json.loads(str(r)) for r in reference["image_table"]]
    scalars = {json.loads(str(r))["day"]: json.loads(str(r)) for r in reference["day_scalars"]}
    template, _ = scene_grid(
        CAMPAIGN / f"optimal_m5_teacher20_cci_exact_control_20260831/{SCENE}.npz"
    )
    banks = {
        "cs": read_planes(ROOT / "cloudscore_cache_cs", template),
        "cs_cdf": read_planes(ROOT / "cloudscore_cache_full", template),
    }
    shared = set(banks["cs"]) & set(banks["cs_cdf"])
    month_index = {month: position for position, month in enumerate(reference["months"])}
    results = {}
    for band, planes in banks.items():
        by_month: dict[str, list[Candidate]] = defaultdict(list)
        for record in images:
            key = str(record["idx"])
            if key not in shared:
                continue
            day = str(record["day"])
            by_month[day[:7]].append(
                Candidate(
                    day=day,
                    score=planes[key],
                    day_weight=float(scalars[day]["weight"]),
                    coverage_ratio=float(record["ratio"]),
                )
            )
        by_month = {m: c for m, c in by_month.items() if m in month_index}
        months, winners, _ = compose_monthly_winners(by_month)
        scores = [
            agreement(reference["winners"][month_index[month]], winners[position])[
                "identical_fraction"
            ]
            for position, month in enumerate(months)
        ]
        results[band] = {"months": len(months), "mean_identical": float(np.mean(scores))}
        print(
            f"{band:7} months={len(months):>2} mean identical vs reference = {np.mean(scores):.4f}",
            flush=True,
        )
    if len(results) == 2:
        delta = results["cs"]["mean_identical"] - results["cs_cdf"]["mean_identical"]
        print(
            f"\nswitching the ordering band cs_cdf -> cs changes agreement by {delta:+.4f}",
            flush=True,
        )
        print(f"shared images compared: {len(shared)}", flush=True)


if __name__ == "__main__":
    main()
