"""Quantify how spatially variable Sen2Cor's own AC inputs are within the AOIs.

Sen2Cor produces per-pixel AOT (DDV) and WVP (APDA) rasters that our features
already consume per pixel. This report measures how much they actually vary
inside each retained tile-scene component, so we know whether per-pixel inputs
are binding or whether Sen2Cor effectively used one value per scene.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from tools.aeronet_validation.train_l2a_l1c_harmonizer import _case_ids, load_pairs

ROOT = Path("/gws/ssde/j25a/nceo_isp/public/siac_refactor")
DEFAULT_PAIRS = ROOT / "analysis/l2a_l1c_physical_pairs_lowcloud152_20260716"
DEFAULT_CASES = ROOT / "reports/aod-final-performance-dashboard-20260713/data/all-cases.csv"


def per_scene_spread(values: np.ndarray, scene_code: np.ndarray) -> dict[str, float]:
    order = np.argsort(scene_code, kind="stable")
    values = values[order]
    scenes = scene_code[order]
    boundaries = np.flatnonzero(scenes[1:] != scenes[:-1]) + 1
    starts = np.concatenate(([0], boundaries))
    stops = np.concatenate((boundaries, [values.size]))
    stds, spans, medians = [], [], []
    for start, stop in zip(starts, stops, strict=True):
        chunk = values[start:stop]
        if chunk.size < 20:
            continue
        stds.append(float(np.std(chunk)))
        spans.append(float(np.percentile(chunk, 95) - np.percentile(chunk, 5)))
        medians.append(float(np.median(chunk)))
    stds_array = np.asarray(stds)
    return {
        "components": int(stds_array.size),
        "within_component_std_median": float(np.median(stds_array)),
        "within_component_std_p90": float(np.percentile(stds_array, 90)),
        "within_component_p5_p95_span_median": float(np.median(spans)),
        "fraction_effectively_constant": float(np.mean(stds_array < 1.0e-3)),
        "between_component_std_of_medians": float(np.std(np.asarray(medians))),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pairs", type=Path, default=DEFAULT_PAIRS)
    parser.add_argument("--cases", type=Path, default=DEFAULT_CASES)
    args = parser.parse_args()
    dataset = load_pairs(
        args.pairs,
        _case_ids(args.cases),
        scene_day_max="2023-12-31",
        max_samples_per_scene=300,
        allow_missing_matchups=True,
    )
    scene_code = np.asarray(dataset.scene_code, dtype=np.int64)
    print(
        json.dumps(
            {
                "sen2cor_aot_within_aoi": per_scene_spread(
                    np.asarray(dataset.l2a_aot, dtype=np.float64), scene_code
                ),
                "sen2cor_wvp_cm_within_aoi": per_scene_spread(
                    np.asarray(dataset.l2a_tcwv, dtype=np.float64), scene_code
                ),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
