"""Build the focused cloud-below-20% AOD experiment performance report."""

from __future__ import annotations

import argparse
import csv
import json
import math
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

DEFAULT_ROOT = Path("/gws/ssde/j25a/nceo_isp/public/siac_refactor")
DEFAULT_ANALYSIS = Path("reports/aod-low-cloud-20260711/screen-analysis.json")
DEFAULT_OUTPUT = Path("reports/aod-low-cloud-20260711")
EXPECTED_ERROR = "abs(retrieved - truth) <= 0.05 + 0.15 * truth"

BENCHMARK_LABELS = {
    "historical_r2": "Historical R2 (S2 anchored)",
    "masked_r2": "Masked R2",
    "multisource_tree_v1": "Multi-source selector v1",
    "multisource_tree_v2": "Multi-source selector v2",
}
CANDIDATE_LABELS = {
    "modis_single_day_followup": "MODIS single-day prior",
    "modis_timeseries_smooth_followup": "MODIS time-series smoothed prior",
    "s2_swir_nir_anchor_followup": "S2 SWIR/NIR anchored prior",
    "modis_monthly_anchor_calibrated": "MODIS SWIR/NIR anchor",
    "modis_monthly_anchor_loose_backstop": "MODIS anchor, loose backstop",
    "modis_monthly_anchor_spread_only": "MODIS anchor, spread-only uncertainty",
    "modis_monthly_multigrid_current": "MODIS monthly DB, full uncertainty",
    "modis_monthly_multigrid_june": "MODIS monthly DB, June control",
    "modis_monthly_multigrid_spread_only": "MODIS monthly DB, spread-only",
    "modis_monthly_multigrid_spread_only_aot25": "MODIS monthly DB, spread-only, AOT max 2.5",
    "modis_monthly_multigrid_legacy_resample_control": "MODIS monthly DB, legacy resampling control",
}
REQUIRED_CANDIDATES = {
    "modis_monthly_anchor_calibrated": 24,
    "modis_monthly_multigrid_current": 152,
    "modis_monthly_multigrid_spread_only": 152,
}
PLANNED_CANDIDATE_SCOPE = {
    "modis_single_day_followup": (29, "29-case historical low-cloud overlap"),
    "modis_timeseries_smooth_followup": (29, "29-case historical low-cloud overlap"),
    "s2_swir_nir_anchor_followup": (29, "29-case historical low-cloud overlap"),
    "modis_monthly_anchor_calibrated": (24, "24-case hard screen"),
    "modis_monthly_anchor_loose_backstop": (2, "2-case backstop control"),
    "modis_monthly_anchor_spread_only": (2, "2-case uncertainty control"),
    "modis_monthly_multigrid_current": (152, "152-case low-cloud cohort"),
    "modis_monthly_multigrid_june": (54, "54-case historical overlap"),
    "modis_monthly_multigrid_spread_only": (152, "152-case low-cloud cohort"),
    "modis_monthly_multigrid_spread_only_aot25": (2, "2-case AOT-bound control"),
    "modis_monthly_multigrid_legacy_resample_control": (
        1,
        "1-case resampling control",
    ),
}
SINGLE_PRIOR_LABELS = {
    "b03_chi2": "B02/B03/B04 chi-square",
    "b03_shape": "B02/B03/B04 spectral shape",
    "b03_auto2": "B02/B03/B04 auto2",
    "b03_profile_s05": "Three-band profiled scale, sigma 0.05",
    "b03_profile_s10": "Three-band profiled scale, sigma 0.10",
    "b03_profile_s20": "Three-band profiled scale, sigma 0.20",
    "b03_loo_s05": "Three-band leave-one-out, sigma 0.05",
    "b03_loo_s10": "Three-band leave-one-out, sigma 0.10",
    "b03_loo_s20": "Three-band leave-one-out, sigma 0.20",
    "b03_trimmed": "Three-band trimmed chi-square",
    "b03_trimmed_loose": "Trimmed chi-square, loose backstop",
    "b03_trimmed_bs3": "Trimmed chi-square, 3x backstop width",
    "b03_trimmed_bs10": "Trimmed chi-square, 10x backstop width",
    "b01b03_chi2": "B01/B02/B03/B04 chi-square",
    "b01b03_profile_s10": "Four-band profiled scale",
    "b01b03_loo_s10": "Four-band leave-one-out",
    "b01b03_trimmed": "Four-band trimmed chi-square",
    "b03_anchor_weighted": "NIR/SWIR realization-weighted prior",
}


def _read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _finite(value: object) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def _hit(retrieved: float, truth: float) -> bool:
    return abs(retrieved - truth) <= 0.05 + 0.15 * truth + 1.0e-12


def _surface_label(directory: str) -> str:
    label = directory.removeprefix("phaseD_results_campaign250_")
    replacements = {
        "K2": "K2 S2 L2A monthly + RT anchor",
        "R2_full_localdiag_20260705": "R2 S2 anchored + gated tau (local diagnostics)",
        "R2_full": "R2 S2 anchored + gated tau",
        "masked_r2_l2awvp_6s_20260710": "Masked R2",
        "masked_pc_direct_l2awvp_6s_20260711": "Masked PC direct",
        "l2a_monthly_median3_scene_mean": "L2A monthly median-3",
        "l2a_pc_production_mean3_scene_mean": "L2A PC production mean-3",
        "l2a_pc_direct_acix_20260705": "L2A PC direct ACIX",
        "l1clut_floor015_scene_mean_20260705": "L1C LUT floor 0.015",
    }
    return replacements.get(label, label.replace("_", " "))


def _candidate_scope(name: str, candidate: dict[str, Any]) -> dict[str, Any]:
    if name.startswith("modis_monthly_anchor"):
        return candidate["screen"]
    return candidate["cohort"]


def _failures(summary: dict[str, Any]) -> int:
    return sum(
        count
        for status, count in summary.get("statuses", {}).items()
        if status not in {"OK", "NO_VALID_OBSERVATION"}
    )


def _required_completion_issues(analysis: dict[str, Any]) -> list[str]:
    issues = []
    for name, expected in REQUIRED_CANDIDATES.items():
        candidate = analysis["candidates"][name]
        summary = _candidate_scope(name, candidate)
        if int(summary["present"]) != expected:
            issues.append(f"{name}: {summary['present']}/{expected} records")
        if _failures(summary):
            issues.append(f"{name}: {_failures(summary)} failed records")
    return issues


def _single_prior_completion_issues(analysis: dict[str, Any]) -> list[str]:
    expected_screen = int(analysis.get("screen_count", 24))
    variants = analysis.get("variants", {})
    issues = []
    for name, variant in variants.items():
        screen = variant.get("screen", {})
        present = int(screen.get("present", 0))
        if present != expected_screen:
            issues.append(f"single-prior {name} screen: {present}/{expected_screen} records")
        non_ok = {
            status: int(count)
            for status, count in screen.get("statuses", {}).items()
            if status != "OK" and int(count) > 0
        }
        if non_ok:
            issues.append(f"single-prior {name} screen: non-OK statuses {non_ok}")
    promoted = variants.get("b03_chi2", {}).get("cohort", {})
    promoted_present = int(promoted.get("present", 0))
    if promoted_present != 152:
        issues.append(f"single-prior b03_chi2 cohort: {promoted_present}/152 records")
    return issues


def _failure_analysis_completion_issues(analysis: dict[str, Any]) -> list[str]:
    summary = analysis.get("summary", {})
    data_quality = analysis.get("data_quality", {})
    cases = analysis.get("failure_cases", [])
    unresolved = analysis.get("unresolved_cases", [])
    issues = []
    cohort_count = int(summary.get("cohort_count", 0))
    hits = int(summary.get("current_hits", 0))
    misses = int(summary.get("current_misses", 0))
    if cohort_count != 152 or hits + misses != 152:
        issues.append(f"failure analysis score: {hits}+{misses}/{cohort_count}")
    if len(cases) != misses:
        issues.append(f"failure analysis rows: {len(cases)}/{misses}")
    if len(unresolved) != int(summary.get("tested_one_prior_unresolved", -1)):
        issues.append(
            "failure analysis unresolved rows: "
            f"{len(unresolved)}/{summary.get('tested_one_prior_unresolved')}"
        )
    if int(data_quality.get("current_non_ok", -1)) != 0:
        issues.append(f"failure analysis non-OK: {data_quality.get('current_non_ok')}")
    if int(data_quality.get("required_diagnostic_complete", 0)) != 152:
        issues.append(
            f"failure analysis diagnostics: {data_quality.get('required_diagnostic_complete')}/152"
        )
    return issues


def _load_regression_controls(root: Path) -> list[dict[str, Any]]:
    experiment_root = root / "experiments/lowcloud20_modis_monthly_database_20260711"
    bambey_id = "Bambey-ISRA__T28PCB_20240610T113321"
    aaq2_id = "AAQ2_SK_Suwon__T52SCG_20240401T021529"
    specs = (
        (
            "Bambey",
            bambey_id,
            "June baseline",
            root / "runs/monthly_database" / bambey_id / "result.json",
            0.0,
            2.5,
            6,
            "Pre-July result",
        ),
        (
            "Bambey",
            bambey_id,
            "Current full uncertainty",
            experiment_root / "runs/monthly_database" / bambey_id / "result.json",
            1.0,
            4.0,
            1,
            "July uncertainty propagation plus widened AOT grid",
        ),
        (
            "Bambey",
            bambey_id,
            "Spread-only, AOT max 4.0",
            experiment_root / "runs/monthly_database_spread_only" / bambey_id / "result.json",
            0.0,
            4.0,
            1,
            "Candidate used for full-cohort validation",
        ),
        (
            "Bambey",
            bambey_id,
            "Spread-only, AOT max 2.5",
            experiment_root / "runs/monthly_database_spread_only_aot25" / bambey_id / "result.json",
            0.0,
            2.5,
            1,
            "Controlled June reproduction",
        ),
        (
            "AAQ2 Suwon",
            aaq2_id,
            "June baseline",
            root / "runs/monthly_database" / aaq2_id / "result.json",
            0.0,
            2.5,
            6,
            "Historical hit",
        ),
        (
            "AAQ2 Suwon",
            aaq2_id,
            "Spread-only, AOT max 4.0",
            experiment_root / "runs/monthly_database_spread_only" / aaq2_id / "result.json",
            0.0,
            4.0,
            1,
            "Current deterministic candidate",
        ),
        (
            "AAQ2 Suwon",
            aaq2_id,
            "Spread-only, AOT max 2.5",
            experiment_root / "runs/monthly_database_spread_only_aot25" / aaq2_id / "result.json",
            0.0,
            2.5,
            1,
            "Current deterministic reproduction control",
        ),
        (
            "AAQ2 Suwon",
            aaq2_id,
            "AOT max 2.5, legacy resampling",
            experiment_root
            / "runs/monthly_database_spread_only_aot25_legacy_resample"
            / aaq2_id
            / "result.json",
            0.0,
            2.5,
            6,
            "Thread-count control",
        ),
    )
    rows = []
    for scene, matchup_id, label, path, uncertainty_scale, aot_max, workers, note in specs:
        truth_record = _read_json(
            root / "phaseD_results_campaign250_R2_full_localdiag_20260705" / f"{matchup_id}.json"
        )
        truth = float(truth_record["truth"])
        record = _read_json(path)
        retrieved = float(record["aot_window_mean"])
        rows.append(
            {
                "scene": scene,
                "configuration": label,
                "uncertainty_scale": uncertainty_scale,
                "aot_max": aot_max,
                "resample_workers": workers,
                "truth": truth,
                "retrieved": retrieved,
                "error": retrieved - truth,
                "within_ee": _hit(retrieved, truth),
                "aot_window_std": _finite(record.get("aot_window_std")),
                "note": note,
            }
        )
    return rows


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    fields = sorted({key for row in rows for key in row})
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def _build_datasets(
    analysis: dict[str, Any],
    regression_controls: list[dict[str, Any]],
    single_consensus: dict[str, Any] | None = None,
    single_adaptive: dict[str, Any] | None = None,
    failure_analysis: dict[str, Any] | None = None,
) -> dict[str, list[dict[str, Any]]]:
    cohort_count = int(analysis["cohort_count"])
    required_hits = int(analysis["required_hits"])
    target_rate = required_hits / cohort_count
    benchmarks = analysis["benchmarks"]
    candidates = analysis["candidates"]
    surface_methods = analysis["surface_method_summaries"]

    benchmark_rows = []
    for name, benchmark in benchmarks.items():
        summary = benchmark["cohort"]
        benchmark_rows.append(
            {
                "method": BENCHMARK_LABELS.get(name, name),
                "key": name,
                "hits": int(summary["hits"]),
                "valid": int(summary["valid"]),
                "expected": int(summary["expected"]),
                "strict_rate": float(summary["strict_rate"]),
                "coverage": float(summary["valid"]) / int(summary["expected"]),
                "failed": _failures(summary),
                "scope": "complete low-cloud cohort",
            }
        )
    for name in (
        "modis_monthly_multigrid_current",
        "modis_monthly_multigrid_spread_only",
    ):
        summary = candidates[name]["cohort"]
        benchmark_rows.append(
            {
                "method": CANDIDATE_LABELS[name],
                "key": name,
                "hits": int(summary["hits"]),
                "valid": int(summary["valid"]),
                "expected": int(summary["expected"]),
                "strict_rate": float(summary["strict_rate"]),
                "coverage": float(summary["valid"]) / int(summary["expected"]),
                "failed": _failures(summary),
                "scope": "complete low-cloud cohort",
            }
        )
    benchmark_rows.sort(key=lambda row: (-row["strict_rate"], row["method"]))

    full_hit_ids = set(candidates["modis_monthly_multigrid_current"]["cohort"]["hit_matchup_ids"])
    spread_hit_ids = set(
        candidates["modis_monthly_multigrid_spread_only"]["cohort"]["hit_matchup_ids"]
    )
    uncertainty_ab = [
        {
            "comparison": "Spread-only versus full uncertainty",
            "expected": cohort_count,
            "full_hits": len(full_hit_ids),
            "spread_hits": len(spread_hit_ids),
            "gains": len(spread_hit_ids - full_hit_ids),
            "losses": len(full_hit_ids - spread_hit_ids),
            "net_hits": len(spread_hit_ids) - len(full_hit_ids),
        }
    ]

    surface_rows = []
    for directory, summary in surface_methods.items():
        surface_rows.append(
            {
                "method": _surface_label(directory),
                "directory": directory,
                "hits": int(summary["hits"]),
                "valid": int(summary["valid"]),
                "expected": int(summary["expected"]),
                "strict_rate": float(summary["strict_rate"]),
                "coverage": float(summary["valid"]) / int(summary["expected"]),
                "rmse": _finite(summary.get("rmse")),
                "bias": _finite(summary.get("bias")),
                "failed": _failures(summary),
            }
        )
    surface_rows.sort(key=lambda row: (-row["hits"], row["method"]))

    candidate_rows = []
    for name, candidate in candidates.items():
        summary = _candidate_scope(name, candidate)
        planned_expected, scope = PLANNED_CANDIDATE_SCOPE.get(
            name, (int(summary["expected"]), "candidate run")
        )
        candidate_rows.append(
            {
                "method": CANDIDATE_LABELS.get(name, name),
                "key": name,
                "present": int(summary["present"]),
                "expected": planned_expected,
                "valid": int(summary["valid"]),
                "hits": int(summary["hits"]),
                "strict_rate": int(summary["hits"]) / planned_expected,
                "unique_hits": int(candidate["unique_vs_existing_surface_oracle"]),
                "rmse": _finite(summary.get("rmse")),
                "bias": _finite(summary.get("bias")),
                "coverage": int(summary["present"]) / planned_expected,
                "failed": _failures(summary),
                "scope": scope,
            }
        )
    candidate_rows.sort(key=lambda row: (-row["unique_hits"], row["method"]))

    hard_screen_rows = []
    screen_count = int(analysis["screen_count"])
    for name, candidate in candidates.items():
        summary = candidate["screen"]
        hard_screen_rows.append(
            {
                "method": CANDIDATE_LABELS.get(name, name),
                "key": name,
                "present": int(summary["present"]),
                "expected": screen_count,
                "valid": int(summary["valid"]),
                "hits": int(summary["hits"]),
                "unique_hits": int(candidate["unique_vs_existing_surface_oracle"]),
                "coverage": int(summary["present"]) / screen_count,
                "failed": _failures(summary),
            }
        )
    hard_screen_rows.sort(key=lambda row: (-row["unique_hits"], row["method"]))

    base_hits = int(analysis["existing_surface_oracle_hits"])
    # The summaries expose counts, while unique candidate IDs expose the exact
    # additions over the existing oracle. Count the staged union from those additions.
    june_unique = set(candidates["modis_monthly_multigrid_june"]["unique_hit_matchup_ids"])
    spread_unique = set(candidates["modis_monthly_multigrid_spread_only"]["unique_hit_matchup_ids"])
    current_unique = set(candidates["modis_monthly_multigrid_current"]["unique_hit_matchup_ids"])
    version_unique = june_unique | spread_unique | current_unique
    all_candidate_unique = {
        matchup_id
        for candidate in candidates.values()
        for matchup_id in candidate["unique_hit_matchup_ids"]
    }
    version_union_hits = base_hits + len(version_unique)
    expanded_hits = base_hits + len(all_candidate_unique)
    if expanded_hits != int(analysis["expanded_surface_oracle_hits"]):
        raise ValueError("Candidate unique-hit union does not match expanded oracle")
    progression = [
        {
            "stage": "Existing complete surface arms",
            "hits": base_hits,
            "rate": base_hits / cohort_count,
            "added_hits": 0,
        },
        {
            "stage": "Plus June MODIS monthly DB",
            "hits": base_hits + len(june_unique),
            "rate": (base_hits + len(june_unique)) / cohort_count,
            "added_hits": len(june_unique),
        },
        {
            "stage": "Plus spread-only AOT max 4.0 MODIS DB",
            "hits": base_hits + len(june_unique | spread_unique),
            "rate": (base_hits + len(june_unique | spread_unique)) / cohort_count,
            "added_hits": len(spread_unique - june_unique),
        },
        {
            "stage": "Plus full-uncertainty AOT max 4.0 MODIS DB",
            "hits": version_union_hits,
            "rate": version_union_hits / cohort_count,
            "added_hits": len(current_unique - (june_unique | spread_unique)),
        },
    ]
    remaining_unique = all_candidate_unique - version_unique
    if remaining_unique:
        progression.append(
            {
                "stage": "Plus remaining controlled candidates",
                "hits": expanded_hits,
                "rate": expanded_hits / cohort_count,
                "added_hits": len(remaining_unique),
            }
        )

    required_failures = sum(
        _failures(_candidate_scope(name, candidates[name])) for name in REQUIRED_CANDIDATES
    )
    best_surface_hits = max(row["hits"] for row in surface_rows)
    selector_v2 = next(row for row in benchmark_rows if row["key"] == "multisource_tree_v2")
    headline = [
        {
            "cohort_count": cohort_count,
            "required_hits": required_hits,
            "target_rate": target_rate,
            "historical_hits": int(benchmarks["historical_r2"]["cohort"]["hits"]),
            "historical_rate": float(benchmarks["historical_r2"]["cohort"]["strict_rate"]),
            "best_surface_hits": best_surface_hits,
            "best_surface_rate": best_surface_hits / cohort_count,
            "selector_hits": selector_v2["hits"],
            "selector_rate": selector_v2["strict_rate"],
            "expanded_oracle_hits": int(analysis["expanded_surface_oracle_hits"]),
            "expanded_oracle_rate": float(analysis["expanded_surface_oracle_rate"]),
            "required_experiment_failures": required_failures,
        }
    ]
    datasets = {
        "headline": headline,
        "benchmark_scores": benchmark_rows,
        "surface_methods": surface_rows,
        "candidate_completion": candidate_rows,
        "hard_screen_candidates": hard_screen_rows,
        "uncertainty_ab": uncertainty_ab,
        "oracle_progression": progression,
        "regression_controls": regression_controls,
    }
    if single_consensus is not None:
        available = int(single_consensus["available"])
        score_specs = (
            ("Historical R2", "Measured baseline", single_consensus["baseline"], True),
            (
                "Three-band chi-square",
                "Measured retrieval",
                single_consensus["b03_retrieval"],
                True,
            ),
            (
                f"Fixed consensus ({single_consensus['best_fixed_posterior']['candidate']})",
                "Validation-screened fixed posterior",
                single_consensus["best_fixed_posterior"],
                True,
            ),
            (
                "Depth-2 site-grouped gate",
                "Out-of-fold diagnostic",
                single_consensus["tree_oof_gate"]["summary"],
                True,
            ),
            (
                "Single-threshold site-grouped gate",
                "Out-of-fold diagnostic",
                single_consensus["oof_gate"],
                True,
            ),
            (
                "ExtraTrees site-grouped model audit",
                "Out-of-fold model audit selected after comparison",
                single_consensus["extra_trees_oof_gate"]["summary"],
                False,
            ),
        )
        score_rows = [
            {
                "method": method,
                "evidence": evidence,
                "hits": int(summary["hits"]),
                "count": int(summary["count"]),
                "strict_rate": float(summary["rate"]),
                "rmse": _finite(summary.get("rmse")),
                "bias": _finite(summary.get("bias")),
                "operational": operational,
            }
            for method, evidence, summary, operational in score_specs
            if summary is not None
        ]
        oracle_hits = int(single_consensus["baseline_plus_fixed_posterior_oracle_hits"])
        score_rows.append(
            {
                "method": "R2 plus consensus hindsight oracle",
                "evidence": "Truth-selected ceiling",
                "hits": oracle_hits,
                "count": available,
                "strict_rate": oracle_hits / available,
                "rmse": None,
                "bias": None,
                "operational": False,
            }
        )
        expanded_oracle_hits = int(single_consensus["expanded_single_prior_oracle_hits"])
        score_rows.append(
            {
                "method": "Expanded one-prior hindsight oracle",
                "evidence": "Truth-selected ceiling including measured retrieval and widened-backstop screen",
                "hits": expanded_oracle_hits,
                "count": available,
                "strict_rate": expanded_oracle_hits / available,
                "rmse": None,
                "bias": None,
                "operational": False,
            }
        )
        datasets["single_prior_scores"] = score_rows
        headline[0].update(
            {
                "single_prior_available": available,
                "single_prior_b03_hits": int(single_consensus["b03_retrieval"]["hits"]),
                "single_prior_fixed_hits": int(single_consensus["best_fixed_posterior"]["hits"]),
                "single_prior_tree_hits": int(single_consensus["tree_oof_gate"]["summary"]["hits"]),
                "single_prior_oracle_hits": expanded_oracle_hits,
                "single_prior_oracle_rate": expanded_oracle_hits / available,
                "single_prior_model_audit_hits": int(
                    single_consensus["extra_trees_oof_gate"]["summary"]["hits"]
                ),
            }
        )
    if single_adaptive is not None:
        hard_rows = []
        for key, variant in single_adaptive["variants"].items():
            summary = variant["screen"]
            if int(summary["present"]) != int(single_adaptive["screen_count"]):
                continue
            hard_rows.append(
                {
                    "method": SINGLE_PRIOR_LABELS.get(key, key),
                    "key": key,
                    "present": int(summary["present"]),
                    "expected": int(single_adaptive["screen_count"]),
                    "hits": int(summary["hits"]),
                    "valid": int(summary["valid"]),
                    "failed": _failures(summary),
                    "coverage": int(summary["present"]) / int(single_adaptive["screen_count"]),
                }
            )
        hard_rows.sort(key=lambda row: (-row["hits"], row["method"]))
        datasets["single_prior_hard_screen"] = hard_rows
    if failure_analysis is not None:
        failure_summary = dict(failure_analysis["summary"])
        datasets["failure_headline"] = [failure_summary]
        datasets["failure_mechanisms"] = [
            {
                **row,
                "mechanism_label": f"{row['mechanism_code']}. {row['mechanism']}",
            }
            for row in failure_analysis["mechanisms"]
        ]
        datasets["failure_truth_bins"] = list(failure_analysis["truth_bins"])
        datasets["failure_cloud_bins"] = list(failure_analysis["cloud_bins"])
        datasets["failure_severity"] = list(failure_analysis["severity"])
        datasets["failure_transitions"] = list(failure_analysis["baseline_transitions"])
        datasets["failure_risk_flags"] = list(failure_analysis["risk_flags"])
        datasets["failure_cases"] = [
            {
                **row,
                "mechanism_label": f"{row['mechanism_code']}. {row['mechanism']}",
            }
            for row in failure_analysis["failure_cases"]
        ]
        datasets["failure_unresolved"] = [
            {
                **row,
                "mechanism_label": f"{row['mechanism_code']}. {row['mechanism']}",
            }
            for row in failure_analysis["unresolved_cases"]
        ]
        headline[0].update(
            {
                "current_failure_misses": int(failure_summary["current_misses"]),
                "current_failure_underreads": int(failure_summary["underreads"]),
                "current_failure_overreads": int(failure_summary["overreads"]),
                "current_failure_unresolved": int(failure_summary["tested_one_prior_unresolved"]),
            }
        )
    return datasets


def _source(generated_at: str) -> dict[str, Any]:
    return {
        "id": "low-cloud-aod-analysis",
        "label": "Low-cloud AOD experiment analysis",
        "path": "report-data.csv",
        "query": {
            "engine": "DuckDB",
            "language": "sql",
            "sql": "SELECT * FROM read_csv_auto('report-data.csv', header = true);",
            "description": "Load the bounded reviewed rows generated from terminal experiment JSON records.",
            "executed_at": generated_at,
            "tables_used": ["report-data.csv"],
            "filters": [
                "Fixed cohort contains 152 campaign matchups with masked-R2 cloud_frac < 0.20",
                "A valid retrieval requires terminal status OK and finite truth and retrieved AOD",
                "Strict scores divide by all 152 expected matchups; no-valid and failed outcomes are misses",
                "The hard screen contains the 24 cohort matchups missed by every complete campaign surface arm",
                "The failure suite contains the 42 status-OK b03_chi2 retrievals outside EE",
            ],
            "metric_definitions": [
                f"Within expected error: {EXPECTED_ERROR}",
                "Strict within-EE rate: within-EE hits divided by the fixed expected denominator",
                "Coverage: terminal finite retrievals divided by the fixed expected denominator",
                "Surface oracle: a matchup is a hit if any included surface arm is within EE; this is a ceiling, not a runtime selector",
                "Unique hit: a candidate hit on a matchup missed by every complete campaign surface arm",
                "Failure severity: absolute AOD error divided by the matchup-specific expected-error threshold",
                "Failure mechanisms are mutually exclusive observational signatures, not causal labels",
            ],
        },
    }


def _artifact(
    *,
    generated_at: str,
    datasets: dict[str, list[dict[str, Any]]],
) -> dict[str, Any]:
    headline = datasets["headline"][0]
    source = _source(generated_at)
    title = "Low-cloud AOD experiment performance"
    current_row = next(
        row
        for row in datasets["benchmark_scores"]
        if row["key"] == "modis_monthly_multigrid_current"
    )
    spread_row = next(
        row
        for row in datasets["benchmark_scores"]
        if row["key"] == "modis_monthly_multigrid_spread_only"
    )
    anchor_row = next(
        row
        for row in datasets["candidate_completion"]
        if row["key"] == "modis_monthly_anchor_calibrated"
    )
    june_row = next(
        row
        for row in datasets["candidate_completion"]
        if row["key"] == "modis_monthly_multigrid_june"
    )
    spread_screen_row = next(
        row
        for row in datasets["candidate_completion"]
        if row["key"] == "modis_monthly_multigrid_spread_only"
    )
    current_candidate_row = next(
        row
        for row in datasets["candidate_completion"]
        if row["key"] == "modis_monthly_multigrid_current"
    )
    daily_row = next(
        row for row in datasets["candidate_completion"] if row["key"] == "modis_single_day_followup"
    )
    smooth_row = next(
        row
        for row in datasets["candidate_completion"]
        if row["key"] == "modis_timeseries_smooth_followup"
    )
    s2_anchor_row = next(
        row
        for row in datasets["candidate_completion"]
        if row["key"] == "s2_swir_nir_anchor_followup"
    )
    base_oracle_hits = datasets["oracle_progression"][0]["hits"]
    version_union_hits = datasets["oracle_progression"][-1]["hits"]
    uncertainty_ab = datasets["uncertainty_ab"][0]
    has_single_prior = bool(datasets.get("single_prior_scores"))
    has_failure_analysis = bool(datasets.get("failure_headline"))
    if uncertainty_ab["net_hits"] <= 0:
        uncertainty_recommendation = (
            "Keep full composite-uncertainty propagation as the default and retain "
            "spread-only as an explicit diagnostic alternative."
        )
    else:
        uncertainty_recommendation = (
            "Keep the uncertainty scale explicit and validate the spread-only gain on "
            "unseen sites before changing the default."
        )

    cards = [
        {
            "id": "cohort-card",
            "dataset": "headline",
            "sourceId": source["id"],
            "description": "Frozen masked-R2 cloud-below-20% denominator.",
            "metrics": [
                {"label": "Low-cloud matchups", "field": "cohort_count", "format": "number"}
            ],
        },
        {
            "id": "target-card",
            "dataset": "headline",
            "sourceId": source["id"],
            "description": "The first attainable strict score above 87% on 152 matchups.",
            "metrics": [
                {"label": "Required hits", "field": "required_hits", "format": "number"},
                {"label": "Required rate", "field": "target_rate", "format": "percent"},
            ],
        },
        {
            "id": "selector-card",
            "dataset": "headline",
            "sourceId": source["id"],
            "description": "Observed zero-failure score of the rounded multi-source selector; validation-fitted, not an independent holdout result.",
            "metrics": [
                {"label": "Selector v2", "field": "selector_rate", "format": "percent"},
                {"label": "Hits", "field": "selector_hits", "format": "number"},
            ],
        },
        {
            "id": "oracle-card",
            "dataset": "headline",
            "sourceId": source["id"],
            "description": "Best-case union across surface outputs and controlled MODIS versions; not a deployable single method.",
            "metrics": [
                {
                    "label": "Expanded surface oracle",
                    "field": "expanded_oracle_rate",
                    "format": "percent",
                },
                {"label": "Hits", "field": "expanded_oracle_hits", "format": "number"},
            ],
        },
        {
            "id": "failure-card",
            "dataset": "headline",
            "sourceId": source["id"],
            "description": "Terminal FAILED records in the three required new experiment sets.",
            "metrics": [
                {
                    "label": "Execution failures",
                    "field": "required_experiment_failures",
                    "format": "number",
                }
            ],
        },
    ]
    if has_single_prior:
        cards.append(
            {
                "id": "single-prior-card",
                "dataset": "headline",
                "sourceId": source["id"],
                "description": "Best site-grouped result and hindsight ceiling using one S2 SWIR/NIR-anchored surface-prior type.",
                "metrics": [
                    {
                        "label": "Best OOF audit",
                        "field": "single_prior_model_audit_hits",
                        "format": "number",
                    },
                    {
                        "label": "Hindsight oracle",
                        "field": "single_prior_oracle_hits",
                        "format": "number",
                    },
                ],
            }
        )
    if has_failure_analysis:
        cards.append(
            {
                "id": "failure-diagnostic-card",
                "dataset": "failure_headline",
                "sourceId": source["id"],
                "description": "Measured b03_chi2 misses and the subset unresolved by every tested one-prior output.",
                "metrics": [
                    {"label": "Measured misses", "field": "current_misses", "format": "number"},
                    {"label": "Underreads", "field": "underreads", "format": "number"},
                    {
                        "label": "One-prior unresolved",
                        "field": "tested_one_prior_unresolved",
                        "format": "number",
                    },
                ],
            }
        )

    charts = [
        {
            "id": "benchmark-chart",
            "title": "Strict within-EE score by complete method",
            "subtitle": "Fixed 152-matchup cloud-below-20% denominator; target is 133 hits (87.5%)",
            "type": "horizontalBar",
            "intent": "comparison",
            "question": "Which complete methods meet the low-cloud target?",
            "rationale": "A ranked bar keeps the fixed denominator and method labels visible.",
            "comparisonContext": {
                "baseline": "Historical R2",
                "denominator": "152 low-cloud matchups",
                "grain": "method",
                "unit": "strict within-EE fraction",
            },
            "dataset": "benchmark_scores",
            "sourceId": source["id"],
            "encodings": {
                "x": {"field": "method", "type": "nominal", "label": "Method"},
                "y": {
                    "field": "strict_rate",
                    "type": "quantitative",
                    "label": "Strict within EE",
                },
                "tooltip": [
                    {"field": "hits", "type": "quantitative"},
                    {"field": "valid", "type": "quantitative"},
                    {"field": "failed", "type": "quantitative"},
                ],
            },
            "valueFormat": "percent",
            "layout": "full",
            "palette": {"kind": "sequential"},
            "referenceLines": [
                {
                    "axis": "y",
                    "value": headline["target_rate"],
                    "label": "87% target (133/152)",
                    "color": "neutral",
                }
            ],
            "settings": {
                "sort": "descending",
                "showValues": True,
                "categoryLabelPolicy": "wrap",
            },
        },
        {
            "id": "oracle-chart",
            "title": "Surface-method oracle progression",
            "subtitle": "Best-case union of saved outputs; useful as a ceiling, not as a runtime score",
            "type": "horizontalBar",
            "intent": "comparison",
            "question": "How many previously unsolved matchups do the MODIS controls add?",
            "rationale": "Staged bars separate the existing surface ceiling from version-specific additions.",
            "comparisonContext": {
                "baseline": "Existing complete campaign surface arms",
                "denominator": "152 low-cloud matchups",
                "grain": "oracle stage",
                "unit": "within-EE hits",
            },
            "dataset": "oracle_progression",
            "sourceId": source["id"],
            "encodings": {
                "x": {"field": "stage", "type": "nominal", "label": "Included outputs"},
                "y": {"field": "hits", "type": "quantitative", "label": "Oracle hits"},
                "tooltip": [
                    {"field": "rate", "type": "quantitative", "format": "percent"},
                    {"field": "added_hits", "type": "quantitative"},
                ],
            },
            "valueFormat": "number",
            "unit": "hits",
            "layout": "full",
            "palette": {"kind": "sequential"},
            "referenceLines": [
                {
                    "axis": "y",
                    "value": headline["required_hits"],
                    "label": "Target",
                    "color": "neutral",
                }
            ],
            "settings": {"showValues": True, "categoryLabelPolicy": "wrap"},
        },
        {
            "id": "screen-chart",
            "title": "Unique fixes on the 24-case hard screen",
            "subtitle": "A unique fix is within EE where every complete campaign surface arm missed",
            "type": "horizontalBar",
            "intent": "comparison",
            "question": "Which new MODIS variants solve genuinely new low-cloud cases?",
            "rationale": "Unique-hit counts directly measure added support without mixing full and partial denominators.",
            "comparisonContext": {
                "baseline": "Existing surface oracle",
                "denominator": "24 hard-screen matchups",
                "grain": "candidate",
                "unit": "unique within-EE hits",
            },
            "dataset": "hard_screen_candidates",
            "sourceId": source["id"],
            "encodings": {
                "x": {"field": "method", "type": "nominal", "label": "Candidate"},
                "y": {"field": "unique_hits", "type": "quantitative", "label": "Unique hits"},
                "tooltip": [
                    {"field": "present", "type": "quantitative"},
                    {"field": "expected", "type": "quantitative"},
                    {"field": "hits", "type": "quantitative"},
                    {"field": "failed", "type": "quantitative"},
                ],
            },
            "valueFormat": "number",
            "unit": "hits",
            "layout": "full",
            "palette": {"kind": "sequential"},
            "settings": {
                "sort": "descending",
                "showValues": True,
                "categoryLabelPolicy": "wrap",
            },
        },
    ]
    if has_single_prior:
        charts.append(
            {
                "id": "single-prior-chart",
                "title": "One-prior operational results and diagnostic ceiling",
                "subtitle": "One S2 monthly SWIR/NIR-anchored ExtraTree prior; fixed 152-matchup denominator",
                "type": "horizontalBar",
                "intent": "comparison",
                "question": "Can a single surface-prior type exceed 87% without truth-selected routing?",
                "rationale": "Ranked bars separate measured and out-of-fold methods from the hindsight oracle.",
                "comparisonContext": {
                    "baseline": "Historical R2",
                    "denominator": "152 low-cloud matchups",
                    "grain": "single-prior method",
                    "unit": "strict within-EE fraction",
                },
                "dataset": "single_prior_scores",
                "sourceId": source["id"],
                "encodings": {
                    "x": {"field": "method", "type": "nominal", "label": "Method"},
                    "y": {
                        "field": "strict_rate",
                        "type": "quantitative",
                        "label": "Strict within EE",
                    },
                    "tooltip": [
                        {"field": "hits", "type": "quantitative"},
                        {"field": "count", "type": "quantitative"},
                        {"field": "evidence", "type": "nominal"},
                    ],
                },
                "valueFormat": "percent",
                "layout": "full",
                "palette": {"kind": "sequential"},
                "referenceLines": [
                    {
                        "axis": "y",
                        "value": headline["target_rate"],
                        "label": "Target (133/152)",
                        "color": "neutral",
                    }
                ],
                "settings": {
                    "sort": "descending",
                    "showValues": True,
                    "categoryLabelPolicy": "wrap",
                },
            }
        )
    if has_failure_analysis:
        failure_summary = datasets["failure_headline"][0]
        charts.extend(
            [
                {
                    "id": "failure-mechanism-chart",
                    "title": "Measured misses by diagnostic mechanism",
                    "subtitle": "42 outside-EE b03_chi2 retrievals; categories are mutually exclusive observational signatures",
                    "type": "horizontalBar",
                    "intent": "composition",
                    "question": "Which diagnostic signatures account for the measured failures?",
                    "rationale": "A ranked bar makes the mutually exclusive failure decomposition and long mechanism labels readable.",
                    "comparisonContext": {
                        "baseline": "All measured b03_chi2 misses",
                        "denominator": "42 outside-EE retrievals",
                        "grain": "diagnostic mechanism",
                        "unit": "misses",
                    },
                    "dataset": "failure_mechanisms",
                    "sourceId": source["id"],
                    "encodings": {
                        "x": {
                            "field": "mechanism_label",
                            "type": "nominal",
                            "label": "Mechanism",
                        },
                        "y": {"field": "count", "type": "quantitative", "label": "Misses"},
                        "tooltip": [
                            {
                                "field": "share_of_misses",
                                "type": "quantitative",
                                "format": "percent",
                            },
                            {"field": "underreads", "type": "quantitative"},
                            {"field": "overreads", "type": "quantitative"},
                            {"field": "median_ee_ratio", "type": "quantitative"},
                        ],
                    },
                    "valueFormat": "number",
                    "unit": "misses",
                    "layout": "full",
                    "palette": {"kind": "sequential"},
                    "settings": {
                        "sort": "descending",
                        "showValues": True,
                        "categoryLabelPolicy": "wrap",
                    },
                },
                {
                    "id": "failure-risk-chart",
                    "title": "Miss rate when diagnostic indicators are present",
                    "subtitle": "All 152 matchups; cloud fraction >=5% is included as a negative control",
                    "type": "horizontalBar",
                    "intent": "comparison",
                    "question": "Which observable flags concentrate failures rather than merely occurring often?",
                    "rationale": "Flag-conditioned miss rates compare candidate risk indicators on one fixed outcome definition.",
                    "comparisonContext": {
                        "baseline": "Overall measured miss rate",
                        "denominator": "Flagged matchups among 152",
                        "grain": "diagnostic indicator",
                        "unit": "outside-EE fraction",
                    },
                    "dataset": "failure_risk_flags",
                    "sourceId": source["id"],
                    "encodings": {
                        "x": {"field": "indicator", "type": "nominal", "label": "Indicator"},
                        "y": {
                            "field": "miss_rate",
                            "type": "quantitative",
                            "label": "Miss rate",
                        },
                        "tooltip": [
                            {"field": "flagged_count", "type": "quantitative"},
                            {"field": "flagged_misses", "type": "quantitative"},
                            {
                                "field": "miss_to_hit_prevalence_ratio",
                                "type": "quantitative",
                            },
                        ],
                    },
                    "valueFormat": "percent",
                    "layout": "full",
                    "palette": {"kind": "sequential"},
                    "referenceLines": [
                        {
                            "axis": "y",
                            "value": failure_summary["current_misses"]
                            / failure_summary["cohort_count"],
                            "label": "Overall miss rate",
                            "color": "neutral",
                        }
                    ],
                    "settings": {
                        "sort": "descending",
                        "showValues": True,
                        "categoryLabelPolicy": "wrap",
                    },
                },
                {
                    "id": "failure-truth-chart",
                    "title": "Measured miss rate by AERONET AOD",
                    "subtitle": "Fixed 152-matchup low-cloud cohort; bin denominators range from 7 to 35",
                    "type": "bar",
                    "intent": "comparison",
                    "question": "Where in the truth-AOD range do failures concentrate?",
                    "rationale": "Ordered bins preserve the physical AOD progression while showing denominator-aware miss rates.",
                    "comparisonContext": {
                        "baseline": "All truth-AOD bins",
                        "denominator": "Matchups in each AERONET AOD bin",
                        "grain": "truth-AOD bin",
                        "unit": "outside-EE fraction",
                    },
                    "dataset": "failure_truth_bins",
                    "sourceId": source["id"],
                    "encodings": {
                        "x": {"field": "truth_bin", "type": "ordinal", "label": "AERONET AOD"},
                        "y": {
                            "field": "miss_rate",
                            "type": "quantitative",
                            "label": "Miss rate",
                        },
                        "tooltip": [
                            {"field": "expected", "type": "quantitative"},
                            {"field": "misses", "type": "quantitative"},
                            {"field": "underreads", "type": "quantitative"},
                            {"field": "overreads", "type": "quantitative"},
                        ],
                    },
                    "valueFormat": "percent",
                    "layout": "full",
                    "palette": {"kind": "sequential"},
                    "settings": {"showValues": True},
                },
                {
                    "id": "failure-severity-chart",
                    "title": "How far failures lie outside expected error",
                    "subtitle": "Absolute error divided by each matchup's EE threshold; all 42 measured misses",
                    "type": "bar",
                    "intent": "distribution",
                    "question": "Are most failures marginal or materially outside EE?",
                    "rationale": "Ordered normalized-error bins separate near-boundary misses from severe retrieval failures.",
                    "comparisonContext": {
                        "baseline": "EE boundary at 1.0",
                        "denominator": "42 outside-EE retrievals",
                        "grain": "normalized-error bin",
                        "unit": "misses",
                    },
                    "dataset": "failure_severity",
                    "sourceId": source["id"],
                    "encodings": {
                        "x": {
                            "field": "severity",
                            "type": "ordinal",
                            "label": "Absolute error / EE",
                        },
                        "y": {"field": "count", "type": "quantitative", "label": "Misses"},
                        "tooltip": [
                            {"field": "underreads", "type": "quantitative"},
                            {"field": "overreads", "type": "quantitative"},
                        ],
                    },
                    "valueFormat": "number",
                    "unit": "misses",
                    "layout": "full",
                    "palette": {"kind": "sequential"},
                    "settings": {"showValues": True},
                },
            ]
        )

    tables = [
        {
            "id": "surface-table",
            "title": "Complete campaign surface arms",
            "subtitle": "All saved surface methods with at least 244 campaign records, rescored on the same 152 low-cloud matchups",
            "dataset": "surface_methods",
            "sourceId": source["id"],
            "layout": "full",
            "density": "dense",
            "defaultSort": {"field": "hits", "direction": "desc"},
            "columns": [
                {"field": "method", "label": "Surface arm", "type": "text"},
                {"field": "hits", "label": "EE hits", "format": "number"},
                {"field": "strict_rate", "label": "Strict EE", "format": "percent"},
                {"field": "valid", "label": "Valid", "format": "number"},
                {"field": "coverage", "label": "Coverage", "format": "percent"},
                {"field": "rmse", "label": "RMSE", "format": "number"},
                {"field": "bias", "label": "Bias", "format": "number", "signed": True},
                {"field": "failed", "label": "Failed", "format": "number"},
            ],
        },
        {
            "id": "candidate-table",
            "title": "MODIS and anchored-prior experiment coverage",
            "subtitle": "Full-cohort methods, historical overlaps, and intentional screens retain their actual planned denominators",
            "dataset": "candidate_completion",
            "sourceId": source["id"],
            "layout": "full",
            "density": "dense",
            "defaultSort": {"field": "unique_hits", "direction": "desc"},
            "columns": [
                {"field": "method", "label": "Experiment", "type": "text"},
                {"field": "scope", "label": "Scope", "type": "text"},
                {"field": "present", "label": "Present", "format": "number"},
                {"field": "expected", "label": "Expected", "format": "number"},
                {"field": "valid", "label": "Valid", "format": "number"},
                {"field": "hits", "label": "EE hits", "format": "number"},
                {"field": "strict_rate", "label": "Strict EE", "format": "percent"},
                {"field": "unique_hits", "label": "Unique", "format": "number"},
                {"field": "rmse", "label": "RMSE", "format": "number"},
                {"field": "bias", "label": "Bias", "format": "number", "signed": True},
                {"field": "failed", "label": "Failed", "format": "number"},
            ],
        },
        {
            "id": "regression-table",
            "title": "MODIS monthly-database regression controls",
            "subtitle": "Bambey reproduces the June configuration; AAQ2 remains a non-reproduced historical hit",
            "dataset": "regression_controls",
            "sourceId": source["id"],
            "layout": "full",
            "density": "spacious",
            "defaultSort": {"field": "aot_max", "direction": "asc"},
            "columns": [
                {"field": "scene", "label": "Scene", "type": "text"},
                {"field": "configuration", "label": "Configuration", "type": "text"},
                {"field": "uncertainty_scale", "label": "Quality scale", "format": "number"},
                {"field": "aot_max", "label": "AOT max", "format": "number"},
                {"field": "resample_workers", "label": "Grid workers", "format": "number"},
                {"field": "truth", "label": "AERONET", "format": "number"},
                {"field": "retrieved", "label": "Retrieved", "format": "number"},
                {"field": "error", "label": "Error", "format": "number", "signed": True},
                {"field": "within_ee", "label": "Within EE", "format": "boolean"},
                {"field": "note", "label": "Role", "type": "text"},
            ],
        },
    ]
    if datasets.get("single_prior_hard_screen"):
        tables.append(
            {
                "id": "single-prior-hard-table",
                "title": "Single-prior hard-case objective screen",
                "subtitle": "Same 24 matchups for every complete screen; all variants use one S2 SWIR/NIR-anchored prior type",
                "dataset": "single_prior_hard_screen",
                "sourceId": source["id"],
                "layout": "full",
                "density": "dense",
                "defaultSort": {"field": "hits", "direction": "desc"},
                "columns": [
                    {"field": "method", "label": "Variant", "type": "text"},
                    {"field": "hits", "label": "EE hits", "format": "number"},
                    {"field": "present", "label": "Present", "format": "number"},
                    {"field": "expected", "label": "Expected", "format": "number"},
                    {"field": "valid", "label": "Valid", "format": "number"},
                    {"field": "coverage", "label": "Coverage", "format": "percent"},
                    {"field": "failed", "label": "Failed", "format": "number"},
                ],
            }
        )
    if has_failure_analysis:
        tables.extend(
            [
                {
                    "id": "failure-mechanism-table",
                    "title": "Failure mechanisms and next discriminating tests",
                    "subtitle": "Mutually exclusive evidence groups over all 42 measured misses",
                    "dataset": "failure_mechanisms",
                    "sourceId": source["id"],
                    "layout": "full",
                    "density": "spacious",
                    "defaultSort": {"field": "count", "direction": "desc"},
                    "columns": [
                        {"field": "mechanism_label", "label": "Mechanism", "type": "text"},
                        {"field": "count", "label": "Misses", "format": "number"},
                        {"field": "share_of_misses", "label": "Share", "format": "percent"},
                        {"field": "underreads", "label": "Under", "format": "number"},
                        {"field": "overreads", "label": "Over", "format": "number"},
                        {
                            "field": "median_ee_ratio",
                            "label": "Median error / EE",
                            "format": "number",
                        },
                        {
                            "field": "diagnostic_evidence",
                            "label": "Observed evidence",
                            "type": "text",
                        },
                        {"field": "next_test", "label": "Next discriminating test", "type": "text"},
                    ],
                },
                {
                    "id": "failure-unresolved-table",
                    "title": "Twelve cases unsolved by every tested one-prior output",
                    "subtitle": "Truth-selected union of measured R2, b03_chi2, fixed posteriors, and all completed hard-screen variants",
                    "dataset": "failure_unresolved",
                    "sourceId": source["id"],
                    "layout": "full",
                    "density": "dense",
                    "defaultSort": {"field": "ee_ratio", "direction": "desc"},
                    "columns": [
                        {"field": "matchup_id", "label": "Matchup", "type": "text"},
                        {"field": "truth_aod", "label": "AERONET", "format": "number"},
                        {"field": "retrieved_aod", "label": "Retrieved", "format": "number"},
                        {
                            "field": "error",
                            "label": "Error",
                            "format": "number",
                            "movement": True,
                        },
                        {"field": "ee_ratio", "label": "Error / EE", "format": "number"},
                        {"field": "mechanism_label", "label": "Mechanism", "type": "text"},
                        {"field": "cams_aod", "label": "CAMS", "format": "number"},
                        {"field": "band_b02_min_aod", "label": "B02 min", "format": "number"},
                        {"field": "band_b03_min_aod", "label": "B03 min", "format": "number"},
                        {"field": "band_b04_min_aod", "label": "B04 min", "format": "number"},
                        {"field": "curve_min_aod", "label": "Global curve min", "format": "number"},
                        {"field": "cloud_fraction", "label": "Cloud", "format": "percent"},
                    ],
                },
                {
                    "id": "failure-case-table",
                    "title": "All 42 measured failures ranked by normalized error",
                    "subtitle": "Status-OK b03_chi2 retrievals outside EE; error / EE makes severity comparable across truth AOD",
                    "dataset": "failure_cases",
                    "sourceId": source["id"],
                    "layout": "full",
                    "density": "dense",
                    "defaultSort": {"field": "severity_rank", "direction": "asc"},
                    "columns": [
                        {"field": "severity_rank", "label": "Rank", "format": "number"},
                        {"field": "matchup_id", "label": "Matchup", "type": "text"},
                        {"field": "truth_aod", "label": "AERONET", "format": "number"},
                        {"field": "retrieved_aod", "label": "Retrieved", "format": "number"},
                        {
                            "field": "error",
                            "label": "Error",
                            "format": "number",
                            "movement": True,
                        },
                        {"field": "ee_ratio", "label": "Error / EE", "format": "number"},
                        {"field": "direction", "label": "Direction", "type": "text"},
                        {"field": "mechanism_label", "label": "Mechanism", "type": "text"},
                        {
                            "field": "baseline_transition",
                            "label": "Versus historical R2",
                            "type": "text",
                        },
                        {"field": "oof_within_ee", "label": "OOF fixes", "format": "boolean"},
                        {
                            "field": "tested_one_prior_recoverable",
                            "label": "Oracle recoverable",
                            "format": "boolean",
                        },
                        {"field": "cloud_fraction", "label": "Cloud", "format": "percent"},
                        {"field": "band_spread", "label": "Band spread", "format": "number"},
                        {"field": "cams_aod", "label": "CAMS", "format": "number"},
                    ],
                },
            ]
        )

    blocks = [
        {"id": "title", "type": "markdown", "body": f"# {title}"},
        {
            "id": "technical-summary",
            "type": "markdown",
            "sourceId": source["id"],
            "body": (
                "## Technical summary\n\n"
                f"**The fixed low-cloud cohort contains {headline['cohort_count']} matchups, so exceeding 87% requires at least "
                f"{headline['required_hits']}/{headline['cohort_count']} = {100 * headline['target_rate']:.1f}%.** Historical R2 scores "
                f"{headline['historical_hits']}/{headline['cohort_count']} = {100 * headline['historical_rate']:.1f}%, and the best complete single surface arm remains "
                f"{headline['best_surface_hits']}/{headline['cohort_count']} = {100 * headline['best_surface_rate']:.1f}%.\n\n"
                f"**The observed multi-source selector clears the numerical target with {headline['selector_hits']}/{headline['cohort_count']} = "
                f"{100 * headline['selector_rate']:.1f}% and no missing retrievals.** It was designed after inspecting this campaign, so this is a validation-fitted result rather than independent evidence.\n\n"
                f"**The expanded surface-output oracle reaches {headline['expanded_oracle_hits']}/{headline['cohort_count']} = "
                f"{100 * headline['expanded_oracle_rate']:.1f}%.** That proves the saved surface methods contain enough complementary hits to meet the target, but no deployable selector currently identifies all of them without AERONET truth.\n\n"
                f"**The required new experiment sets contain {headline['required_experiment_failures']} terminal FAILED records.** The actual MODIS SWIR/NIR anchored surface solver fixes "
                f"{anchor_row['unique_hits']} of its {anchor_row['present']} completed hard cases, with bias {anchor_row['bias']:.3f}; the standard MODIS monthly-database path is the useful branch."
            ),
        },
        {
            "id": "headline-strip",
            "type": "metric-strip",
            "cardIds": [
                "cohort-card",
                "target-card",
                "selector-card",
                "oracle-card",
                "failure-card",
            ],
        },
        {
            "id": "method-finding",
            "type": "markdown",
            "sourceId": source["id"],
            "body": (
                "## Only the multi-source selectors exceed 87% as complete methods\n\n"
                f"The rounded selector v2 records {headline['selector_hits']} strict hits. By comparison, the full-uncertainty MODIS monthly database records "
                f"{current_row['hits']}, and the corrected spread-only arm records {spread_row['hits']}. Spread-only gains {uncertainty_ab['gains']} paired hits and loses {uncertainty_ab['losses']} versus full uncertainty, for a net {uncertainty_ab['net_hits']:+d}. The chart uses one fixed denominator for every bar; failures and no-valid outcomes remain misses."
            ),
        },
        {"id": "benchmark-chart-block", "type": "chart", "chartId": "benchmark-chart"},
        {
            "id": "oracle-finding",
            "type": "markdown",
            "sourceId": source["id"],
            "body": (
                "## Complementary MODIS versions close the surface-method gap\n\n"
                f"The existing complete surface arms solve {base_oracle_hits} cases in an oracle union. "
                f"The June MODIS monthly-database control adds {june_row['unique_hits']} genuinely new hits, while the current spread-only version with AOT maximum 4.0 adds {spread_screen_row['unique_hits']}. "
                f"The full-uncertainty version adds {current_candidate_row['unique_hits']}; after overlap, the three versions contribute {version_union_hits - base_oracle_hits} distinct additions and their union reaches {version_union_hits} hits. This staged result is a diagnostic ceiling, not permission to select whichever version matches truth after the fact."
            ),
        },
        {"id": "oracle-chart-block", "type": "chart", "chartId": "oracle-chart"},
        {
            "id": "screen-finding",
            "type": "markdown",
            "sourceId": source["id"],
            "body": (
                "## The MODIS monthly database adds hard-case hits; the anchored surface solver does not\n\n"
                f"The 24-case screen contains only matchups missed by every complete campaign surface arm. The full-uncertainty monthly-database route adds {current_candidate_row['unique_hits']} unique hits and spread-only adds {spread_screen_row['unique_hits']}; the formulations are complementary rather than one globally dominating the other. In the anchored route, the database is built from MODIS MCD43 monthly BRDF surfaces and queried with the same-scene Sentinel-2 NIR/SWIR bands; it is not an S2 L2A surface dictionary. That actual MODIS anchored route was repaired so it runs, but its completed screen is 0/24 with bias {anchor_row['bias']:.3f} and RMSE {anchor_row['rmse']:.3f}."
            ),
        },
        {"id": "screen-chart-block", "type": "chart", "chartId": "screen-chart"},
        {
            "id": "temporal-prior-finding",
            "type": "markdown",
            "sourceId": source["id"],
            "body": (
                "## Daily and time-series MODIS priors remain partial controls\n\n"
                f"The earlier controlled surface-prior screen intersects this low-cloud cohort in 29 cases. MODIS single-day retrieves {daily_row['hits']}/29 within EE, Whittaker time-series smoothing retrieves {smooth_row['hits']}/29, and the S2 L2A SWIR/NIR-anchored control retrieves {s2_anchor_row['hits']}/29. All 29 planned records are present for each arm and none is a complete 152-case score. The one-hit smoothing improvement is too small and too selectively sampled to support promotion."
            ),
        },
        {"id": "candidate-table-block", "type": "table", "tableId": "candidate-table"},
        {
            "id": "regression-finding",
            "type": "markdown",
            "sourceId": source["id"],
            "body": (
                "## Bambey isolates two scientific regressions\n\n"
                "Full propagation of the matched-composite quality term changes Bambey from a within-EE June retrieval near 0.522 to 0.332. Setting the quality scale to zero restores the historical spread-plus-floor model. Returning the AOT upper bound from 4.0 to 2.5 then reproduces the June result to about 0.0001 AOD, showing that widening a logarithmic bound changes the full grid, not only the high-AOD tail. AAQ2 does not reproduce: both one-worker and six-worker current controls remain near 0.346 rather than June's 0.239, so thread count is not the remaining cause. "
                f"Across the full cohort, spread-only has a net {uncertainty_ab['net_hits']:+d} paired hits versus full propagation, so the Bambey fix is not a global improvement."
            ),
        },
        {"id": "regression-table-block", "type": "table", "tableId": "regression-table"},
        {
            "id": "surface-finding",
            "type": "markdown",
            "sourceId": source["id"],
            "body": (
                "## Complete single-surface arms remain at or below 73%\n\n"
                "Every complete campaign surface output is rescored below on the same 152 IDs. Multiple arms tie historical R2 at 111 hits; changing cloud handling, predictor variants, tau behavior, and L2A surface construction does not supply the 22 additional hits required for a single surface arm to reach 133."
            ),
        },
        {"id": "surface-table-block", "type": "table", "tableId": "surface-table"},
        {
            "id": "scope-method",
            "type": "markdown",
            "sourceId": source["id"],
            "body": (
                "## Scope, metric, and method\n\n"
                "The cohort is frozen by matchup ID from the masked-R2 diagnostic output: `cloud_frac < 0.20`, strictly less than 20%. These are 152 unique matchups. The earlier count of 528 was 8 experiment arms times 66 diagnostic matchups and is not used as a site or score denominator here. Within expected error is `abs(retrieved - AERONET) <= 0.05 + 0.15 * AERONET`. Strict rates divide by all 152 IDs. The 24-case hard screen is the complement of the oracle union over complete campaign surface directories. All metrics are recomputed from saved JSON rather than copied from job logs."
            ),
        },
        {
            "id": "limitations",
            "type": "markdown",
            "sourceId": source["id"],
            "body": (
                "## Limitations and robustness\n\n"
                f"The campaign has been inspected repeatedly, so the {100 * headline['selector_rate']:.1f}% selector score and the {100 * headline['expanded_oracle_rate']:.1f}% expanded oracle are not independent validation. The oracle uses AERONET truth to recognize a hit and cannot run in production. The June MODIS arm covers only 54 of the 152 low-cloud cases; the daily, smoothed, and S2-anchor follow-ups cover 29. They are retained as controlled historical subsets, not scored as complete methods. Cluster runtime and cache state were not controlled, so runtime is excluded from method ranking."
            ),
        },
        {
            "id": "next-steps",
            "type": "markdown",
            "body": (
                "## Recommended next steps\n\n"
                f"1. {uncertainty_recommendation}\n"
                "2. Freeze the SWIR/NIR-anchored prior definition and develop any stability rule using only truth-free diagnostics from a separate development cohort.\n"
                "3. Evaluate that fixed one-prior algorithm on unseen sites or a later-year campaign; do not treat the hindsight oracle as production performance.\n"
                "4. Retain terminal-status audits in every batch so missing or failed jobs cannot inflate strict accuracy."
            ),
        },
        {
            "id": "further-questions",
            "type": "markdown",
            "body": (
                "## Further questions\n\n"
                "- Can grid-bound sensitivity itself identify unstable moderate-AOD solutions without external truth?\n"
                "- Which independent atmospheric product is reliable enough to arbitrate between complementary surface solutions?\n"
                f"- Do the {spread_screen_row['unique_hits']} spread-only MODIS hard-case fixes persist on unseen sites and dates?"
            ),
        },
    ]
    if has_single_prior:
        single_scores = {row["method"]: row for row in datasets["single_prior_scores"]}
        tree_row = single_scores["Depth-2 site-grouped gate"]
        model_audit_row = single_scores["ExtraTrees site-grouped model audit"]
        fixed_row = next(
            row
            for row in datasets["single_prior_scores"]
            if row["method"].startswith("Fixed consensus")
        )
        b03_row = single_scores["Three-band chi-square"]
        oracle_row = single_scores["R2 plus consensus hindsight oracle"]
        expanded_oracle_row = single_scores["Expanded one-prior hindsight oracle"]
        hard_best = max(
            (row["hits"] for row in datasets.get("single_prior_hard_screen", [])),
            default=0,
        )
        blocks[1]["body"] += (
            "\n\n"
            f"**Restricting the retrieval to one S2 SWIR/NIR-anchored surface-prior type does not reach the target.** "
            f"The measured three-band solve records {b03_row['hits']}/{b03_row['count']}; the best fixed consensus records "
            f"{fixed_row['hits']}/{fixed_row['count']}; and the fixed depth-2 site-grouped gate records "
            f"{tree_row['hits']}/{tree_row['count']}. A richer out-of-fold model audit reaches "
            f"{model_audit_row['hits']}/{model_audit_row['count']}, still below target. Only the AERONET-selected hindsight union reaches "
            f"{expanded_oracle_row['hits']}/{expanded_oracle_row['count']}, which is not operational."
        )
        next(block for block in blocks if block["id"] == "headline-strip")["cardIds"].append(
            "single-prior-card"
        )
        insertion = next(
            index for index, block in enumerate(blocks) if block["id"] == "temporal-prior-finding"
        )
        single_blocks = [
            {
                "id": "single-prior-finding",
                "type": "markdown",
                "sourceId": source["id"],
                "body": (
                    "## One anchored surface prior remains below 87%\n\n"
                    f"The three-band chi-square retrieval scores {b03_row['hits']}/{b03_row['count']}. "
                    f"A fixed CAMS-plus-band-consensus posterior scores {fixed_row['hits']}/{fixed_row['count']}; "
                    f"the fixed depth-2 gate is evaluated out of fold by site and scores {tree_row['hits']}/{tree_row['count']}. "
                    f"A less-simple ExtraTrees gate reaches {model_audit_row['hits']}/{model_audit_row['count']} in a model-family audit but was selected after comparing alternatives. "
                    f"The strongest widened-backstop objective recovers {hard_best}/24 cases on the frozen hard screen. "
                    f"The fixed-posterior union reaches {oracle_row['hits']}/{oracle_row['count']}; adding the measured three-band output and widened-backstop screen raises the hindsight ceiling to {expanded_oracle_row['hits']}/{expanded_oracle_row['count']}. Both unions are selected with AERONET truth and therefore measure information availability, not deployable accuracy."
                ),
            },
            {"id": "single-prior-chart-block", "type": "chart", "chartId": "single-prior-chart"},
        ]
        if datasets.get("single_prior_hard_screen"):
            single_blocks.extend(
                [
                    {
                        "id": "single-prior-hard-finding",
                        "type": "markdown",
                        "sourceId": source["id"],
                        "body": (
                            "## Wider backstops expose real surface evidence but do not generalize\n\n"
                            f"Dropping the worst of B02/B03/B04 and widening the atmospheric backstop recovers at most {hard_best}/24 hard cases. "
                            "It fixes scenes such as Bambey where green and red agree against blue, but fails where all visible minima and CAMS are low. "
                            "Profiled brightness, leave-one-band-out, four-band, and realization-weighted variants add little or no hard-case coverage."
                        ),
                    },
                    {
                        "id": "single-prior-hard-table-block",
                        "type": "table",
                        "tableId": "single-prior-hard-table",
                    },
                ]
            )
        blocks[insertion:insertion] = single_blocks
        next(block for block in blocks if block["id"] == "limitations")["body"] += (
            " The single-prior fixed and tree gates were developed after inspecting this campaign; the tree score is site-grouped out of fold, but model-family selection remains development-set informed."
        )
    if has_failure_analysis:
        failure_summary = datasets["failure_headline"][0]
        mechanism_rows = {row["mechanism_code"]: row for row in datasets["failure_mechanisms"]}
        risk_rows = {row["indicator"]: row for row in datasets["failure_risk_flags"]}
        truth_rows = {row["truth_bin"]: row for row in datasets["failure_truth_bins"]}
        moderate_high_misses = truth_rows["0.4-0.6"]["misses"] + truth_rows["0.6-1.0"]["misses"]
        moderate_high_count = truth_rows["0.4-0.6"]["expected"] + truth_rows["0.6-1.0"]["expected"]
        blocks[1]["body"] += (
            "\n\n"
            f"**The measured three-band run has {failure_summary['current_misses']} outside-EE cases: "
            f"{failure_summary['underreads']} underreads and {failure_summary['overreads']} overreads.** "
            f"The misses span {failure_summary['unique_failure_sites']} sites and "
            f"{failure_summary['unique_failure_scenes']} scenes, so this is a mechanism problem rather than one bad site. "
            f"A truth-selected union of every tested one-prior output reaches {failure_summary['tested_one_prior_oracle_hits']}/"
            f"{failure_summary['cohort_count']}, leaving {failure_summary['tested_one_prior_unresolved']} cases with no tested one-prior solution."
        )
        next(block for block in blocks if block["id"] == "headline-strip")["cardIds"].append(
            "failure-diagnostic-card"
        )
        insertion = next(
            index for index, block in enumerate(blocks) if block["id"] == "temporal-prior-finding"
        )
        failure_blocks = [
            {
                "id": "failure-overview",
                "type": "markdown",
                "sourceId": source["id"],
                "body": (
                    "## The failures split into extraction disagreement, band conflict, and missing AOD signal\n\n"
                    f"The two largest mutually exclusive groups each contain 12 cases: a global median cost-curve minimum inside EE while the final scene output is outside, and visible-band minima spread by at least 0.2 AOD. Another {mechanism_rows['B']['count']} cases place CAMS and all three visible-band minima below EE; these are the clearest information-gap failures. The remaining groups are coherent visible underreads ({mechanism_rows['D']['count']}), coherent visible overreads ({mechanism_rows['E']['count']}), and residual underreads ({mechanism_rows['F']['count']}). These labels describe observed solver evidence, not proven physical causes."
                ),
            },
            {
                "id": "failure-mechanism-chart-block",
                "type": "chart",
                "chartId": "failure-mechanism-chart",
            },
            {
                "id": "failure-risk-finding",
                "type": "markdown",
                "sourceId": source["id"],
                "body": (
                    "## Low visible-band minima are the strongest risk flag; cloud fraction is not\n\n"
                    f"When all visible minima lie below EE, {risk_rows['All visible minima below EE']['flagged_misses']}/{risk_rows['All visible minima below EE']['flagged_count']} = {100 * risk_rows['All visible minima below EE']['miss_rate']:.1f}% of cases fail. Truth AOD at least 0.4, CAMS below EE, and band spread at least 0.2 each roughly double failure prevalence relative to hits. By contrast, cloud fraction at least 5% has a {100 * risk_rows['Cloud fraction >= 5%']['miss_rate']:.1f}% miss rate, essentially the overall {100 * failure_summary['current_misses'] / failure_summary['cohort_count']:.1f}%. Within this already-low-cloud cohort, residual cloud fraction does not explain the failures."
                ),
            },
            {
                "id": "failure-risk-chart-block",
                "type": "chart",
                "chartId": "failure-risk-chart",
            },
            {
                "id": "failure-truth-finding",
                "type": "markdown",
                "sourceId": source["id"],
                "body": (
                    "## Moderate and high AOD dominate the miss burden\n\n"
                    f"AERONET AOD 0.4-1.0 contributes {moderate_high_misses}/{moderate_high_count} misses ({100 * moderate_high_misses / moderate_high_count:.1f}%) and 22 of all 42 failures. In those two bins, 20 of 22 misses are underreads. The cleanest bin below 0.1 misses only {truth_rows['<0.1']['misses']}/{truth_rows['<0.1']['expected']} cases. The >=1.0 bin is small (n={truth_rows['>=1.0']['expected']}), but its two failures are the two largest absolute errors."
                ),
            },
            {
                "id": "failure-truth-chart-block",
                "type": "chart",
                "chartId": "failure-truth-chart",
            },
            {
                "id": "failure-severity-finding",
                "type": "markdown",
                "sourceId": source["id"],
                "body": (
                    "## Sixteen misses are at least 1.5 times outside the EE boundary\n\n"
                    f"The median failure is {failure_summary['median_ee_ratio']:.2f} times the matchup-specific EE threshold. "
                    f"{failure_summary['misses_ge_1_5_ee']} cases are at least 1.5 times EE and "
                    f"{failure_summary['misses_ge_2_ee']} are at least twice EE, so the gap is not only a collection of boundary-rounding cases."
                ),
            },
            {
                "id": "failure-severity-chart-block",
                "type": "chart",
                "chartId": "failure-severity-chart",
            },
            {
                "id": "failure-mechanism-detail",
                "type": "markdown",
                "sourceId": source["id"],
                "body": (
                    "## Each failure group requires a different discriminating test\n\n"
                    f"The 12 global-curve/output disagreements motivate a fixed extraction comparison, but replacing every output with the global curve minimum would score only {risk_rows['Global curve inside EE']['flagged_count']}/152 and is not a solution. Band-conflict cases require a full-cohort robust-likelihood test with break accounting. The all-anchor-low group requires new information or better aerosol optics rather than another selector over the same signals."
                ),
            },
            {
                "id": "failure-mechanism-table-block",
                "type": "table",
                "tableId": "failure-mechanism-table",
            },
            {
                "id": "failure-unresolved-finding",
                "type": "markdown",
                "sourceId": source["id"],
                "body": (
                    "## Twelve cases remain unsolved by every tested one-prior output\n\n"
                    f"The corrected hindsight union includes historical R2, the measured three-band solve, every fixed CAMS/band posterior, and all 17 completed alternative hard-screen formulations. It reaches {failure_summary['tested_one_prior_oracle_hits']}/152 = {100 * failure_summary['tested_one_prior_oracle_rate']:.1f}% and leaves the 12 cases below. Because the union uses AERONET truth to choose a winner, 92.1% is an information ceiling, not operational accuracy."
                ),
            },
            {
                "id": "failure-unresolved-table-block",
                "type": "table",
                "tableId": "failure-unresolved-table",
            },
            {
                "id": "failure-case-finding",
                "type": "markdown",
                "sourceId": source["id"],
                "body": (
                    "## The complete bad suite is auditable case by case\n\n"
                    f"Historical R2 fixes {failure_summary['baseline_regressions']} current misses but the new solve fixes {failure_summary['baseline_fixes']} historical misses, a net loss of one hit. The ExtraTrees OOF audit fixes {failure_summary['oof_fixes_current_misses']} of the 42 measured misses but breaks {failure_summary['oof_breaks_current_hits']} of 110 measured hits, explaining its net score of {failure_summary['oof_hits']}/152. The table ranks all failures by absolute error divided by EE and retains the transition, mechanism, cloud, CAMS, and band-spread diagnostics."
                ),
            },
            {
                "id": "failure-case-table-block",
                "type": "table",
                "tableId": "failure-case-table",
            },
        ]
        blocks[insertion:insertion] = failure_blocks
        next(block for block in blocks if block["id"] == "scope-method")["body"] += (
            " The failure suite is the 42 status-OK three-band retrievals outside EE. Mechanisms are assigned hierarchically from global-curve agreement, joint CAMS/visible evidence, band spread, coherent visible direction, then residual error direction."
        )
        next(block for block in blocks if block["id"] == "limitations")["body"] += (
            " Failure mechanisms are threshold-based observational groups derived after inspecting this campaign. They prioritize experiments but do not establish physical causality. The 140-hit one-prior union and every recoverability label use AERONET truth."
        )
        next(block for block in blocks if block["id"] == "next-steps")["body"] += (
            "\n5. Run the extraction comparison for group A and one fixed robust band-likelihood for group C across all 152 cases, reporting fixes and breaks.\n"
            "6. Treat group B as an information-gap cohort for aerosol-optics or independent-AOD experiments, not further selector tuning."
        )
        next(block for block in blocks if block["id"] == "further-questions")["body"] += (
            "\n- Can a truth-free local-versus-global consistency statistic isolate group A without breaking the 68 current hits whose global curve is also inside EE?\n"
            "- Which additional observation supplies high-AOD evidence for the six all-anchor-low failures?"
        )
    return {
        "surface": "report",
        "manifest": {
            "version": 1,
            "surface": "report",
            "title": title,
            "description": "Technical analysis of AOD retrieval performance for a fixed cloud-below-20% Sentinel-2/AERONET cohort.",
            "generatedAt": generated_at,
            "blocks": blocks,
            "cards": cards,
            "charts": charts,
            "tables": tables,
            "sources": [source],
        },
        "snapshot": {
            "version": 1,
            "generatedAt": generated_at,
            "status": "ready",
            "datasets": datasets,
        },
        "sources": [source],
    }


def build(
    root: Path,
    analysis_path: Path,
    output_dir: Path,
    *,
    allow_partial: bool,
) -> dict[str, Any]:
    analysis = _read_json(analysis_path)
    if int(analysis.get("cohort_count", 0)) != 152:
        raise ValueError("Expected the frozen 152-matchup low-cloud cohort")
    issues = _required_completion_issues(analysis)
    single_consensus_path = output_dir / "single-prior-consensus-analysis.json"
    single_adaptive_path = output_dir / "single-prior-adaptive-analysis.json"
    failure_analysis_path = output_dir / "low-cloud-failure-analysis.json"
    single_consensus = _read_json(single_consensus_path) if single_consensus_path.exists() else None
    single_adaptive = _read_json(single_adaptive_path) if single_adaptive_path.exists() else None
    failure_analysis = _read_json(failure_analysis_path) if failure_analysis_path.exists() else None
    if single_consensus is not None and int(single_consensus.get("available", 0)) != 152:
        issues.append(f"single-prior consensus: {single_consensus.get('available', 0)}/152 records")
    if single_adaptive is not None:
        issues.extend(_single_prior_completion_issues(single_adaptive))
    if failure_analysis is not None:
        issues.extend(_failure_analysis_completion_issues(failure_analysis))
    if issues and not allow_partial:
        raise ValueError("Required experiments are incomplete: " + "; ".join(issues))

    regression_controls = _load_regression_controls(root)
    datasets = _build_datasets(
        analysis,
        regression_controls,
        single_consensus=single_consensus,
        single_adaptive=single_adaptive,
        failure_analysis=failure_analysis,
    )
    generated_at = datetime.now(UTC).replace(microsecond=0).isoformat().replace("+00:00", "Z")
    artifact = _artifact(generated_at=generated_at, datasets=datasets)
    chart_map = [
        {
            "section": "Complete methods",
            "question": "Which complete methods meet the target?",
            "family": "Comparison and ranking",
            "type": "horizontal bar with target reference",
            "fields": ["method", "strict_rate", "hits", "valid"],
            "claim": "Only the validation-fitted multi-source selectors exceed 87% as complete methods.",
            "palette": "Single-root sequential blue with a neutral target line.",
        },
        {
            "section": "Surface oracle",
            "question": "How do MODIS versions add complementary hits?",
            "family": "Comparison",
            "type": "staged horizontal bar",
            "fields": ["stage", "hits", "added_hits"],
            "claim": (
                f"The saved-output union reaches "
                f"{datasets['headline'][0]['expanded_oracle_hits']} hits, but remains an oracle ceiling."
            ),
            "palette": "Single-root sequential blue with a neutral target line.",
        },
        {
            "section": "Hard-case screen",
            "question": "Which candidates add genuinely new hits?",
            "family": "Comparison and ranking",
            "type": "horizontal bar",
            "fields": ["method", "unique_hits", "present", "expected"],
            "claim": "The standard MODIS monthly database adds hard-case hits; the anchored surface route does not.",
            "palette": "Single-root sequential blue.",
        },
    ]
    if datasets.get("single_prior_scores"):
        chart_map.append(
            {
                "section": "One-prior performance",
                "question": "Can one S2 SWIR/NIR-anchored prior exceed the target operationally?",
                "family": "Comparison and benchmark",
                "type": "horizontal bar with target reference",
                "fields": ["method", "strict_rate", "hits", "evidence"],
                "claim": "Only the truth-selected one-prior oracle exceeds the target; measured and out-of-fold methods remain below it.",
                "palette": "Single-root sequential blue with a neutral target line.",
            }
        )
    if datasets.get("failure_headline"):
        chart_map.extend(
            [
                {
                    "section": "Measured failure mechanisms",
                    "question": "Which mutually exclusive diagnostic signatures describe the 42 misses?",
                    "family": "Composition",
                    "type": "ranked horizontal bar",
                    "fields": ["mechanism_label", "count", "share_of_misses"],
                    "claim": "Extraction disagreement and visible-band conflict are the two largest groups.",
                    "palette": "Single-root sequential blue.",
                },
                {
                    "section": "Failure indicators",
                    "question": "Which observable flags concentrate misses?",
                    "family": "Comparison",
                    "type": "horizontal bar with overall-rate reference",
                    "fields": ["indicator", "miss_rate", "flagged_count"],
                    "claim": "All visible minima below EE is the strongest risk flag; residual cloud is not.",
                    "palette": "Single-root sequential blue with a neutral reference line.",
                },
                {
                    "section": "Failure distribution",
                    "question": "Where in truth AOD do misses concentrate?",
                    "family": "Ordered comparison",
                    "type": "bar",
                    "fields": ["truth_bin", "miss_rate", "expected", "misses"],
                    "claim": "AOD 0.4-1.0 has the largest miss burden and is dominated by underreads.",
                    "palette": "Single-root sequential blue.",
                },
                {
                    "section": "Failure severity",
                    "question": "How far outside EE are the measured misses?",
                    "family": "Distribution",
                    "type": "ordered binned bar",
                    "fields": ["severity", "count", "underreads", "overreads"],
                    "claim": "Sixteen misses are at least 1.5 times the EE threshold.",
                    "palette": "Single-root sequential blue.",
                },
            ]
        )
    source_notes = {
        "audience": "technical",
        "delivery_mode": "html",
        "required_structure": [
            "title",
            "technical summary",
            "key findings with visual evidence",
            "scope, data, and metric definitions",
            "methodology",
            "limitations, uncertainty, and robustness checks",
            "recommended next steps",
            "further questions",
        ],
        "structure_mapping": {
            "key findings with visual evidence": [
                "complete-method comparison",
                "surface-oracle progression",
                "hard-case candidate screen",
                "single-prior operational comparison",
                "measured failure bad-suite diagnostics",
                "Bambey regression controls",
                "complete surface-arm audit",
            ],
            "methodology": "Merged into Scope, metric, and method.",
        },
        "chart_map": chart_map,
        "completion_issues": issues,
        "validation_assessment": (
            "Ready to share with explicit validation-fitting and oracle caveats."
            if not issues
            else "Partial preview only; required experiment arrays are incomplete."
        ),
        "omissions": [
            "No independent production-performance claim for the multi-source selector.",
            "No deployable performance claim for the method oracle.",
            "No runtime ranking because cache and cluster placement were uncontrolled.",
        ],
    }

    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "artifact.json").write_text(json.dumps(artifact, indent=2), encoding="utf-8")
    (output_dir / "analysis.json").write_text(json.dumps(analysis, indent=2), encoding="utf-8")
    (output_dir / "source-notes.json").write_text(
        json.dumps(source_notes, indent=2), encoding="utf-8"
    )
    report_rows = [{"dataset": name, **row} for name, rows in datasets.items() for row in rows]
    _write_csv(output_dir / "report-data.csv", report_rows)
    for name, rows in datasets.items():
        _write_csv(output_dir / f"{name.replace('_', '-')}.csv", rows)
    return {
        "output_dir": str(output_dir),
        "artifact": str(output_dir / "artifact.json"),
        "completion_issues": issues,
        "allow_partial": allow_partial,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument("--analysis", type=Path)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--allow-partial", action="store_true")
    args = parser.parse_args()
    analysis_path = args.analysis or args.root / DEFAULT_ANALYSIS
    output_dir = args.output_dir or args.root / DEFAULT_OUTPUT
    receipt = build(
        args.root,
        analysis_path,
        output_dir,
        allow_partial=args.allow_partial,
    )
    print(json.dumps(receipt, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
