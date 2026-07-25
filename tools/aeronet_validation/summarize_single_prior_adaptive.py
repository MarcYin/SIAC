"""Summarize single-S2-prior adaptive AOD experiments on the low-cloud cohort."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from tools.aeronet_validation.summarize_low_cloud_aod_screen import (  # noqa: E402
    _load_records,
    _summary,
)

DEFAULT_ROOT = Path("/gws/ssde/j25a/nceo_isp/public/siac_refactor")
DEFAULT_ANALYSIS = Path("reports/aod-low-cloud-20260711/screen-analysis.json")
DEFAULT_OUTPUT = Path("reports/aod-low-cloud-20260711/single-prior-adaptive-analysis.json")
OPTIONS = (
    "b03_chi2",
    "b03_shape",
    "b03_auto2",
    "b03_profile_s05",
    "b03_profile_s10",
    "b03_profile_s20",
    "b03_loo_s05",
    "b03_loo_s10",
    "b03_loo_s20",
    "b03_trimmed",
    "b03_trimmed_loose",
    "b03_trimmed_bs3",
    "b03_trimmed_bs10",
    "b01b03_chi2",
    "b01b03_profile_s10",
    "b01b03_loo_s10",
    "b01b03_trimmed",
    "b03_anchor_weighted",
)


def _hit_set(summary: dict[str, Any]) -> set[str]:
    return set(summary.get("hit_matchup_ids", []))


def summarize(root: Path, analysis_path: Path) -> dict[str, Any]:
    analysis = json.loads(analysis_path.read_text(encoding="utf-8"))
    cohort = [
        line.strip()
        for line in (root / "campaign250_lowcloud20_mids.txt")
        .read_text(encoding="utf-8")
        .splitlines()
        if line.strip()
    ]
    if len(cohort) != len(set(cohort)) or len(cohort) != 152:
        raise ValueError("Expected 152 unique low-cloud matchup IDs")
    screen = list(analysis["screen_matchup_ids"])
    baseline = analysis["benchmarks"]["historical_r2"]["cohort"]
    baseline_hits = _hit_set(baseline)

    variants = {}
    for option in OPTIONS:
        directory = root / f"phaseD_results_lowcloud20_singleprior_{option}_20260711"
        records = _load_records(directory)
        cohort_summary = _summary(records, cohort)
        screen_summary = _summary(records, screen)
        candidate_hits = _hit_set(cohort_summary)
        variants[option] = {
            "directory": str(directory),
            "cohort": cohort_summary,
            "screen": screen_summary,
            "gains_vs_historical_r2": len(candidate_hits - baseline_hits),
            "losses_vs_historical_r2": len(baseline_hits - candidate_hits)
            if cohort_summary["present"] == len(cohort)
            else None,
            "net_vs_historical_r2": len(candidate_hits) - len(baseline_hits)
            if cohort_summary["present"] == len(cohort)
            else None,
        }

    return {
        "cohort_count": len(cohort),
        "screen_count": len(screen),
        "target_hits": 133,
        "surface_prior_type": "S2 monthly SWIR/NIR-anchored ExtraTree predictor",
        "baseline": baseline,
        "variants": variants,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument("--analysis", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    analysis_path = args.analysis or args.root / DEFAULT_ANALYSIS
    output_path = args.output or args.root / DEFAULT_OUTPUT
    result = summarize(args.root, analysis_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    for option, variant in result["variants"].items():
        cohort = variant["cohort"]
        screen = variant["screen"]
        print(
            f"{option}: present={cohort['present']}/152 hits={cohort['hits']} "
            f"hard={screen['hits']}/{screen['present']} statuses={cohort['statuses']}"
        )
    print(output_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
