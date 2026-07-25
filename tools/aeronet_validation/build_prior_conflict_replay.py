"""Build an AERONET-free adaptive atmospheric-prior replay.

The baseline calibrated backstop is retained unless the surface-only tile
optimum is significantly above the atmospheric prior.  At that point the
record from the otherwise identical product-uncertainty solve is selected.
The fixed default z=2.576 is the standard 99% two-sided Gaussian threshold.
Truth fields are carried through for later scoring but never enter selection.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

ROOT = Path("/gws/ssde/j25a/nceo_isp/public/siac_refactor")
DEFAULT_BASELINE = ROOT / "phaseD_results_lowcloud20_native_maiac_adaptive_b03_chi2_20260713"
DEFAULT_PRODUCT_UNC = ROOT / "phaseD_results_lowcloud20_mediumphysics_product_unc_20260713"
DEFAULT_OUTPUT = ROOT / "analysis/lowcloud20_prior_conflict_z2576_development_20260713"


def calibrated_backstop_sigma(aot: float) -> float:
    loose = max(0.5 * aot, 0.02)
    mid = max(0.07, 0.5 * aot / (1.0 + math.exp(-(aot - 0.5) / 0.15)))
    return loose if aot < 0.15 else mid


def replay_record(
    baseline: dict[str, Any],
    product_unc: dict[str, Any],
    *,
    z_threshold: float,
) -> dict[str, Any]:
    atmo_aot = float(baseline["atmo_prior"]["aot_median"])
    sigma = calibrated_backstop_sigma(atmo_aot)
    # Evaluate conflict on the loose-prior solve's own surface likelihood. The
    # harness performs one anchor fixed-point iteration, so its final surface
    # cube can differ from the calibrated-prior pass even though the source and
    # RT model are identical.
    surface_min_raw = product_unc.get("solver", {}).get("surface_cost_curve_min_aot")
    surface_min = float(surface_min_raw) if surface_min_raw is not None else math.nan
    standardized_conflict = (
        (surface_min - atmo_aot) / sigma if math.isfinite(surface_min) else math.nan
    )
    use_product_unc = bool(
        str(product_unc.get("status", "")).upper() == "OK"
        and math.isfinite(standardized_conflict)
        and standardized_conflict > z_threshold
    )

    selected = dict(product_unc if use_product_unc else baseline)
    selected["prior_conflict_replay"] = {
        "uses_aeronet_in_retrieval": False,
        "decision": "product_uncertainty" if use_product_unc else "calibrated_uncertainty",
        "surface_curve_min_aot": surface_min if math.isfinite(surface_min) else None,
        "surface_curve_source": "product_uncertainty_pass",
        "atmospheric_prior_aot": atmo_aot,
        "calibrated_prior_sigma": sigma,
        "standardized_positive_conflict": (
            standardized_conflict if math.isfinite(standardized_conflict) else None
        ),
        "z_threshold": z_threshold,
    }
    return selected


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline-dir", type=Path, default=DEFAULT_BASELINE)
    parser.add_argument("--product-unc-dir", type=Path, default=DEFAULT_PRODUCT_UNC)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--z-threshold", type=float, default=2.576)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.z_threshold <= 0.0:
        raise SystemExit("--z-threshold must be positive")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    selected_loose = 0
    completed = 0
    skipped = 0
    for baseline_path in sorted(args.baseline_dir.glob("*.json")):
        candidate_path = args.product_unc_dir / baseline_path.name
        if not candidate_path.exists():
            skipped += 1
            continue
        baseline = json.loads(baseline_path.read_text())
        candidate = json.loads(candidate_path.read_text())
        record = replay_record(baseline, candidate, z_threshold=args.z_threshold)
        selected_loose += int(
            record["prior_conflict_replay"]["decision"] == "product_uncertainty"
        )
        (args.output_dir / baseline_path.name).write_text(json.dumps(record, indent=2) + "\n")
        completed += 1
    summary = {
        "uses_aeronet_in_retrieval": False,
        "z_threshold": args.z_threshold,
        "completed": completed,
        "selected_product_uncertainty": selected_loose,
        "selected_calibrated_uncertainty": completed - selected_loose,
        "skipped": skipped,
    }
    (args.output_dir / "replay_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, sort_keys=True))


if __name__ == "__main__":
    main()
