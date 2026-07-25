"""Diagnose every outside-EE result in the frozen low-cloud AOD cohort."""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from collections import Counter, defaultdict
from pathlib import Path
from statistics import median
from typing import TYPE_CHECKING, Any

import numpy as np

if TYPE_CHECKING:
    from collections.abc import Callable

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from tools.aeronet_validation.analyze_single_prior_consensus import (  # noqa: E402, I001
    _classification_inputs,
    _load_rows,
)
from tools.aeronet_validation.summarize_single_prior_adaptive import OPTIONS  # noqa: E402


ROOT = Path("/gws/ssde/j25a/nceo_isp/public/siac_refactor")
COHORT_FILE = "campaign250_lowcloud20_mids.txt"
CURRENT_DIR = "phaseD_results_lowcloud20_singleprior_b03_chi2_20260711"
BASELINE_DIR = "phaseD_results_campaign250_R2_full_localdiag_20260705"
PRIOR_DIR = "prior_quality"
OUTPUT_DIR = "reports/aod-low-cloud-20260711"
FIXED_WEIGHTS = (0.175, 0.20, 0.225, 0.25, 0.275, 0.30)

MECHANISM_METADATA = {
    "A": {
        "mechanism": "Global curve inside EE; scene output outside",
        "evidence": "The median scene cost-curve minimum is inside EE, but the final scene-mean retrieval is outside.",
        "next_test": "Compare local station-window, scene-median, and scene-mean extraction with a fixed aggregation rule.",
    },
    "B": {
        "mechanism": "All CAMS and visible anchors below EE",
        "evidence": "CAMS and all B02/B03/B04 AOD minima lie below the lower EE bound.",
        "next_test": "Add independent high-AOD information or improve aerosol optics; selector tuning cannot create a missing signal.",
    },
    "C": {
        "mechanism": "Visible-band disagreement >=0.2 AOD",
        "evidence": "The spread between B02/B03/B04 AOD minima is at least 0.2.",
        "next_test": "Use a fixed robust band-consensus likelihood and validate its break rate on all 152 cases.",
    },
    "D": {
        "mechanism": "Coherent visible minima below EE",
        "evidence": "All visible-band minima are below EE while the CAMS anchor is not also below EE.",
        "next_test": "Audit surface-prior brightness and uncertainty calibration against the atmospheric anchor.",
    },
    "E": {
        "mechanism": "Coherent visible minima above EE",
        "evidence": "All visible-band minima are above the upper EE bound.",
        "next_test": "Audit dark or contaminated surface priors and adjacency before changing the atmospheric anchor.",
    },
    "F": {
        "mechanism": "Residual underread",
        "evidence": "The retrieval is below EE without one of the stronger diagnostic signatures.",
        "next_test": "Inspect local cost fields and aerosol-family sensitivity case by case.",
    },
    "G": {
        "mechanism": "Residual overread",
        "evidence": "The retrieval is above EE without one of the stronger diagnostic signatures.",
        "next_test": "Inspect cloud adjacency, local surface support, and aerosol-family sensitivity case by case.",
    },
}

TRUTH_BINS = (
    ("<0.1", -math.inf, 0.1),
    ("0.1-0.2", 0.1, 0.2),
    ("0.2-0.4", 0.2, 0.4),
    ("0.4-0.6", 0.4, 0.6),
    ("0.6-1.0", 0.6, 1.0),
    (">=1.0", 1.0, math.inf),
)

SEVERITY_BINS = (
    ("1.00-1.25 x EE", 1.0, 1.25),
    ("1.25-1.50 x EE", 1.25, 1.5),
    ("1.50-2.00 x EE", 1.5, 2.0),
    ("2.00-3.00 x EE", 2.0, 3.0),
    (">=3.00 x EE", 3.0, math.inf),
)


def _finite(value: object) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def _load_records(directory: Path) -> dict[str, dict[str, Any]]:
    records: dict[str, dict[str, Any]] = {}
    for path in sorted(directory.glob("*.json")):
        record = json.loads(path.read_text(encoding="utf-8"))
        records[str(record.get("matchup_id") or path.stem)] = record
    return records


def _ee(truth: float) -> float:
    return 0.05 + 0.15 * truth


def _hit_value(value: float, truth: float) -> bool:
    return abs(value - truth) <= _ee(truth) + 1.0e-12


def _record_hit(record: dict[str, Any] | None) -> bool:
    if record is None or str(record.get("status", "")).upper() != "OK":
        return False
    truth = _finite(record.get("truth"))
    retrieved = _finite(record.get("retrieved"))
    return truth is not None and retrieved is not None and _hit_value(retrieved, truth)


def _bin(value: float, bins: tuple[tuple[str, float, float], ...]) -> str:
    for label, lower, upper in bins:
        if lower <= value < upper:
            return label
    raise ValueError(f"Value {value!r} did not match a configured bin")


def _truth_bin(value: float) -> str:
    return _bin(value, TRUTH_BINS)


def _severity_bin(value: float) -> str:
    return _bin(value, SEVERITY_BINS)


def _cloud_bin(value: float) -> str:
    if value < 0.05:
        return "0-5%"
    if value < 0.10:
        return "5-10%"
    return "10-20%"


def _band_minima(record: dict[str, Any]) -> list[float]:
    solver = record.get("solver") or {}
    values = [
        _finite(solver.get(f"surface_band_{band}_argmin_aot")) for band in ("B02", "B03", "B04")
    ]
    if any(value is None for value in values):
        raise ValueError(f"Missing band-minimum diagnostic for {record.get('matchup_id')}")
    return [float(value) for value in values if value is not None]


def _mechanism(record: dict[str, Any], cams_aod: float) -> str:
    truth = float(record["truth"])
    retrieved = float(record["retrieved"])
    lower = truth - _ee(truth)
    upper = truth + _ee(truth)
    solver = record.get("solver") or {}
    bands = _band_minima(record)
    curve_min = float(solver["surface_cost_curve_min_aot"])
    spread = float(solver["surface_band_argmin_spread"])

    if lower <= curve_min <= upper:
        return "A"
    if retrieved < lower and cams_aod < lower and all(value < lower for value in bands):
        return "B"
    if spread >= 0.2:
        return "C"
    if retrieved < lower and all(value < lower for value in bands):
        return "D"
    if retrieved > upper and all(value > upper for value in bands):
        return "E"
    return "F" if retrieved < lower else "G"


def _median(rows: list[dict[str, Any]], field: str) -> float | None:
    values = [float(row[field]) for row in rows if _finite(row.get(field)) is not None]
    return float(median(values)) if values else None


def _extra_trees_oof_predictions(root: Path, cohort: list[str]) -> dict[str, float]:
    from sklearn.ensemble import ExtraTreesClassifier
    from sklearn.impute import SimpleImputer
    from sklearn.pipeline import make_pipeline

    rows = _load_rows(root, cohort)
    if [row["matchup_id"] for row in rows] != cohort:
        raise ValueError("ExtraTrees audit rows do not match the frozen cohort")
    _, matrix, baseline, posterior, _truth, label, folds = _classification_inputs(
        rows,
        weight=0.275,
        aggregation="median",
    )
    prediction = np.full(len(rows), np.nan, dtype=np.float64)
    for fold in range(5):
        test = np.flatnonzero(folds == fold)
        train = np.flatnonzero(folds != fold)
        model = make_pipeline(
            SimpleImputer(strategy="median"),
            ExtraTreesClassifier(
                n_estimators=300,
                max_depth=4,
                min_samples_leaf=6,
                class_weight="balanced",
                random_state=0,
                n_jobs=1,
            ),
        )
        model.fit(matrix[train], label[train])
        choose_posterior = model.predict(matrix[test]).astype(bool)
        prediction[test] = np.where(choose_posterior, posterior[test], baseline[test])
    if not np.isfinite(prediction).all():
        raise ValueError("ExtraTrees audit did not produce a prediction for every matchup")
    return {row["matchup_id"]: float(prediction[index]) for index, row in enumerate(rows)}


def _fixed_posterior_sources(
    record: dict[str, Any],
    cams_aod: float,
) -> list[str]:
    truth = float(record["truth"])
    bands = np.asarray(_band_minima(record), dtype=np.float64)
    sources = []
    for aggregation, surface in (
        ("mean", float(np.mean(bands))),
        ("median", float(np.median(bands))),
    ):
        for weight in FIXED_WEIGHTS:
            value = weight * surface + (1.0 - weight) * cams_aod
            if _hit_value(value, truth):
                sources.append(f"fixed_{aggregation}_{weight:.3f}")
    return sources


def _transition(baseline_hit: bool, current_hit: bool) -> str:
    if baseline_hit and current_hit:
        return "Retained hit"
    if baseline_hit:
        return "Regression: R2 hit to current miss"
    if current_hit:
        return "Fix: R2 miss to current hit"
    return "Retained miss"


def _risk_flag_rows(all_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    definitions: tuple[tuple[str, Callable[[dict[str, Any]], bool]], ...] = (
        ("Truth AOD >= 0.4", lambda row: row["truth_aod"] >= 0.4),
        ("Band spread >= 0.2", lambda row: row["band_spread"] >= 0.2),
        ("All visible minima below EE", lambda row: bool(row["all_bands_low"])),
        ("CAMS below EE", lambda row: bool(row["cams_low"])),
        ("Global curve inside EE", lambda row: bool(row["curve_min_within_ee"])),
        ("Cloud fraction >= 5%", lambda row: row["cloud_fraction"] >= 0.05),
    )
    rows = []
    for indicator, predicate in definitions:
        flagged = [row for row in all_rows if predicate(row)]
        misses = sum(not row["within_ee"] for row in flagged)
        hit_flag_rate = sum(predicate(row) for row in all_rows if row["within_ee"]) / sum(
            row["within_ee"] for row in all_rows
        )
        miss_flag_rate = sum(predicate(row) for row in all_rows if not row["within_ee"]) / sum(
            not row["within_ee"] for row in all_rows
        )
        rows.append(
            {
                "indicator": indicator,
                "flagged_count": len(flagged),
                "flagged_misses": misses,
                "miss_rate": misses / len(flagged) if flagged else None,
                "hit_flag_rate": hit_flag_rate,
                "miss_flag_rate": miss_flag_rate,
                "miss_to_hit_prevalence_ratio": (
                    miss_flag_rate / hit_flag_rate if hit_flag_rate > 0.0 else None
                ),
            }
        )
    return rows


def analyze(root: Path) -> dict[str, Any]:
    cohort = [
        line.strip()
        for line in (root / COHORT_FILE).read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    if len(cohort) != 152 or len(set(cohort)) != 152:
        raise ValueError("Expected 152 unique low-cloud matchup IDs")

    current = _load_records(root / CURRENT_DIR)
    baseline = _load_records(root / BASELINE_DIR)
    priors = _load_records(root / PRIOR_DIR)
    missing = {
        "current": sorted(set(cohort) - set(current)),
        "baseline": sorted(set(cohort) - set(baseline)),
        "prior": sorted(set(cohort) - set(priors)),
    }
    if any(missing.values()):
        raise ValueError(f"Missing required records: {missing}")
    non_ok = [mid for mid in cohort if str(current[mid].get("status", "")).upper() != "OK"]
    if non_ok:
        raise ValueError(f"Current retrieval has non-OK records: {non_ok}")

    variants = {
        option: _load_records(root / f"phaseD_results_lowcloud20_singleprior_{option}_20260711")
        for option in OPTIONS
        if option != "b03_chi2"
    }
    oof = _extra_trees_oof_predictions(root, cohort)
    all_rows: list[dict[str, Any]] = []
    failure_rows: list[dict[str, Any]] = []
    unresolved_rows: list[dict[str, Any]] = []

    for matchup_id in cohort:
        record = current[matchup_id]
        baseline_record = baseline[matchup_id]
        prior = priors[matchup_id]
        solver = record.get("solver") or {}
        truth = float(record["truth"])
        retrieved = float(record["retrieved"])
        error = retrieved - truth
        threshold = _ee(truth)
        lower = truth - threshold
        upper = truth + threshold
        bands = _band_minima(record)
        cams_aod = float(prior["cams_aot"])
        current_hit = _hit_value(retrieved, truth)
        baseline_hit = _record_hit(baseline_record)
        fixed_sources = _fixed_posterior_sources(record, cams_aod)
        variant_sources = [
            option for option, records in variants.items() if _record_hit(records.get(matchup_id))
        ]
        hit_sources = []
        if current_hit:
            hit_sources.append("current_b03_chi2")
        if baseline_hit:
            hit_sources.append("historical_r2")
        hit_sources.extend(fixed_sources)
        hit_sources.extend(variant_sources)
        curve_min = float(solver["surface_cost_curve_min_aot"])
        band_spread = float(solver["surface_band_argmin_spread"])
        row = {
            "matchup_id": matchup_id,
            "site": str(record.get("site") or matchup_id.split("__", 1)[0]),
            "product_id": str(record.get("product_id") or ""),
            "truth_aod": truth,
            "retrieved_aod": retrieved,
            "error": error,
            "abs_error": abs(error),
            "ee_threshold": threshold,
            "ee_ratio": abs(error) / threshold,
            "within_ee": current_hit,
            "direction": "Underread" if error < 0.0 else "Overread",
            "truth_bin": _truth_bin(truth),
            "cloud_bin": _cloud_bin(float(record["cloud_frac"])),
            "cloud_fraction": float(record["cloud_frac"]),
            "valid_support_fraction": float(
                solver.get("surface_valid_observation_fraction", math.nan)
            ),
            "cost_per_band": float(solver["cost_final_per_band"]),
            "curve_min_aod": curve_min,
            "curve_min_within_ee": lower <= curve_min <= upper,
            "curve_relative_second_delta": float(
                solver["surface_cost_curve_relative_second_delta"]
            ),
            "curve_curvature": _finite(solver.get("surface_cost_curve_curvature")),
            "curve_min_at_edge": bool(solver.get("surface_cost_curve_min_at_edge")),
            "band_b02_min_aod": bands[0],
            "band_b03_min_aod": bands[1],
            "band_b04_min_aod": bands[2],
            "band_spread": band_spread,
            "all_bands_low": all(value < lower for value in bands),
            "all_bands_high": all(value > upper for value in bands),
            "any_band_within_ee": any(lower <= value <= upper for value in bands),
            "cams_aod": cams_aod,
            "cams_low": cams_aod < lower,
            "cams_high": cams_aod > upper,
            "baseline_aod": float(baseline_record["retrieved"]),
            "baseline_within_ee": baseline_hit,
            "baseline_transition": _transition(baseline_hit, current_hit),
            "oof_aod": oof[matchup_id],
            "oof_within_ee": _hit_value(oof[matchup_id], truth),
            "tau_gate_fired": bool(solver.get("surface_tau_gate_fired")),
            "backstop_escape_applied": bool(solver.get("surface_backstop_escape_applied")),
            "cloud_mask_bypassed": bool(solver.get("surface_cloud_mask_bypassed")),
            "water_mask_bypassed": bool(solver.get("surface_water_mask_bypassed")),
            "fixed_posterior_hit_count": len(fixed_sources),
            "hard_variant_fixes": ", ".join(variant_sources),
            "tested_one_prior_hit_sources": ", ".join(hit_sources),
            "tested_one_prior_recoverable": bool(hit_sources),
        }
        all_rows.append(row)
        if not current_hit:
            code = _mechanism(record, cams_aod)
            row.update(
                {
                    "mechanism_code": code,
                    "mechanism": MECHANISM_METADATA[code]["mechanism"],
                    "severity": _severity_bin(row["ee_ratio"]),
                }
            )
            failure_rows.append(row)
            if not hit_sources:
                unresolved_rows.append(row)

    failure_rows.sort(key=lambda row: (-row["ee_ratio"], row["matchup_id"]))
    for rank, row in enumerate(failure_rows, start=1):
        row["severity_rank"] = rank
    unresolved_rows.sort(key=lambda row: (-row["ee_ratio"], row["matchup_id"]))

    mechanisms: list[dict[str, Any]] = []
    grouped: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in failure_rows:
        grouped[row["mechanism_code"]].append(row)
    for code in MECHANISM_METADATA:
        rows = grouped.get(code, [])
        if not rows:
            continue
        metadata = MECHANISM_METADATA[code]
        mechanisms.append(
            {
                "mechanism_code": code,
                "mechanism": metadata["mechanism"],
                "count": len(rows),
                "share_of_misses": len(rows) / len(failure_rows),
                "underreads": sum(row["direction"] == "Underread" for row in rows),
                "overreads": sum(row["direction"] == "Overread" for row in rows),
                "median_abs_error": _median(rows, "abs_error"),
                "median_ee_ratio": _median(rows, "ee_ratio"),
                "diagnostic_evidence": metadata["evidence"],
                "next_test": metadata["next_test"],
            }
        )

    truth_bins = []
    for label, _lower, _upper in TRUTH_BINS:
        rows = [row for row in all_rows if row["truth_bin"] == label]
        misses = [row for row in rows if not row["within_ee"]]
        truth_bins.append(
            {
                "truth_bin": label,
                "expected": len(rows),
                "hits": len(rows) - len(misses),
                "misses": len(misses),
                "miss_rate": len(misses) / len(rows) if rows else None,
                "underreads": sum(row["direction"] == "Underread" for row in misses),
                "overreads": sum(row["direction"] == "Overread" for row in misses),
                "median_abs_error_misses": _median(misses, "abs_error"),
            }
        )

    cloud_bins = []
    for label in ("0-5%", "5-10%", "10-20%"):
        rows = [row for row in all_rows if row["cloud_bin"] == label]
        misses = [row for row in rows if not row["within_ee"]]
        cloud_bins.append(
            {
                "cloud_bin": label,
                "expected": len(rows),
                "misses": len(misses),
                "miss_rate": len(misses) / len(rows) if rows else None,
            }
        )

    severity = []
    for label, _lower, _upper in SEVERITY_BINS:
        rows = [row for row in failure_rows if row["severity"] == label]
        severity.append(
            {
                "severity": label,
                "count": len(rows),
                "underreads": sum(row["direction"] == "Underread" for row in rows),
                "overreads": sum(row["direction"] == "Overread" for row in rows),
            }
        )

    transition_counts = Counter(row["baseline_transition"] for row in all_rows)
    transitions = [
        {"outcome": label, "count": transition_counts[label]}
        for label in (
            "Retained hit",
            "Fix: R2 miss to current hit",
            "Regression: R2 hit to current miss",
            "Retained miss",
        )
    ]

    current_misses = len(failure_rows)
    oof_hits = sum(row["oof_within_ee"] for row in all_rows)
    tested_oracle_hits = sum(row["tested_one_prior_recoverable"] for row in all_rows)
    product_counts = Counter(row["product_id"] for row in failure_rows)
    return {
        "summary": {
            "cohort_count": len(cohort),
            "current_hits": len(cohort) - current_misses,
            "current_misses": current_misses,
            "underreads": sum(row["direction"] == "Underread" for row in failure_rows),
            "overreads": sum(row["direction"] == "Overread" for row in failure_rows),
            "median_abs_error": _median(failure_rows, "abs_error"),
            "median_ee_ratio": _median(failure_rows, "ee_ratio"),
            "misses_ge_1_5_ee": sum(row["ee_ratio"] >= 1.5 for row in failure_rows),
            "misses_ge_2_ee": sum(row["ee_ratio"] >= 2.0 for row in failure_rows),
            "unique_failure_sites": len({row["site"] for row in failure_rows}),
            "unique_failure_scenes": len(product_counts),
            "repeated_failure_scenes": sum(count > 1 for count in product_counts.values()),
            "baseline_fixes": transition_counts["Fix: R2 miss to current hit"],
            "baseline_regressions": transition_counts["Regression: R2 hit to current miss"],
            "oof_hits": oof_hits,
            "oof_misses": len(cohort) - oof_hits,
            "oof_fixes_current_misses": sum(
                not row["within_ee"] and row["oof_within_ee"] for row in all_rows
            ),
            "oof_breaks_current_hits": sum(
                row["within_ee"] and not row["oof_within_ee"] for row in all_rows
            ),
            "tested_one_prior_oracle_hits": tested_oracle_hits,
            "tested_one_prior_oracle_rate": tested_oracle_hits / len(cohort),
            "tested_one_prior_unresolved": len(unresolved_rows),
        },
        "mechanisms": mechanisms,
        "truth_bins": truth_bins,
        "cloud_bins": cloud_bins,
        "severity": severity,
        "baseline_transitions": transitions,
        "risk_flags": _risk_flag_rows(all_rows),
        "failure_cases": failure_rows,
        "unresolved_cases": unresolved_rows,
        "data_quality": {
            "cohort_unique": len(set(cohort)),
            "current_present": len(current),
            "current_non_ok": len(non_ok),
            "failure_rows": len(failure_rows),
            "required_diagnostic_complete": sum(
                all(
                    key in (current[mid].get("solver") or {})
                    for key in (
                        "surface_cost_curve_min_aot",
                        "surface_band_argmin_spread",
                        "surface_band_B02_argmin_aot",
                        "surface_band_B03_argmin_aot",
                        "surface_band_B04_argmin_aot",
                    )
                )
                for mid in cohort
            ),
        },
        "definitions": {
            "within_ee": "abs(retrieved - AERONET) <= 0.05 + 0.15 * AERONET",
            "failure_cohort": "Current b03_chi2 status-OK retrievals outside expected error on the fixed 152-matchup cloud<20% cohort",
            "tested_one_prior_oracle": "Truth-selected union of current b03_chi2, historical R2, fixed CAMS/band posteriors, and all completed one-prior hard-screen variants",
            "mechanism_priority": "A global-curve/output disagreement; B all anchors low; C band spread >=0.2; D/E coherent visible minima; F/G residual direction",
            "oof_caveat": "ExtraTrees predictions are site-grouped out of fold, but the model family was selected after comparing alternatives on this campaign",
        },
    }


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    fields = list(rows[0])
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=ROOT)
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args()

    output_dir = args.output_dir or args.root / OUTPUT_DIR
    result = analyze(args.root)
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "low-cloud-failure-analysis.json").write_text(
        json.dumps(result, indent=2) + "\n",
        encoding="utf-8",
    )
    _write_csv(output_dir / "low-cloud-failure-cases.csv", result["failure_cases"])
    _write_csv(output_dir / "low-cloud-failure-mechanisms.csv", result["mechanisms"])
    (output_dir / "low-cloud-failure-mids.txt").write_text(
        "".join(f"{row['matchup_id']}\n" for row in result["failure_cases"]),
        encoding="utf-8",
    )
    (output_dir / "low-cloud-one-prior-unresolved-mids.txt").write_text(
        "".join(f"{row['matchup_id']}\n" for row in result["unresolved_cases"]),
        encoding="utf-8",
    )
    print(json.dumps(result["summary"], indent=2))
    print(output_dir / "low-cloud-failure-analysis.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
