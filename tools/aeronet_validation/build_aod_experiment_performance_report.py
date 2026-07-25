"""Build the completed 66-scene AOD experiment performance report.

The report keeps execution completeness, retrieval coverage, strict expected-error
score, valid-only accuracy, and paired common-case error separate.  This avoids
rewarding a method merely for abstaining on difficult scenes or, conversely,
rewarding an unsafe mask bypass solely for returning more values.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from dataclasses import asdict, dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np
from scipy.stats import binomtest

if TYPE_CHECKING:
    from collections.abc import Iterable

DEFAULT_ROOT = Path("/gws/ssde/j25a/nceo_isp/public/siac_refactor")
DEFAULT_REPORT_DIR = Path("reports/aod-experiment-performance-20260711")
EXPECTED_ERROR_FORMULA = "abs(retrieved - truth) <= 0.05 + 0.15 * truth"
FULL_OLD_DIR = "phaseD_results_campaign250_R2_full_localdiag_20260705"
FULL_MASKED_DIR = "phaseD_results_campaign250_masked_r2_l2awvp_6s_20260710"


@dataclass(frozen=True)
class MethodSpec:
    key: str
    label: str
    short_label: str
    family: str
    directory: str
    surface_prior: str
    rt_backend: str
    tcwv: str
    masking: str
    proper_mask: bool
    controlled_note: str


METHODS = (
    MethodSpec(
        key="modis_single_day",
        label="MODIS single-day prior",
        short_label="MODIS daily",
        family="Surface-prior follow-up",
        directory="phaseD_results_masked_modis_prior_modis_single_day_20260711",
        surface_prior="MCD43 nearest-valid daily kernel model",
        rt_backend="libRadtran LUT",
        tcwv="Atmospheric prior",
        masking="Cloud and water masked",
        proper_mask=True,
        controlled_note="Direct temporal-prior comparison with MODIS time-series smooth.",
    ),
    MethodSpec(
        key="modis_timeseries_smooth",
        label="MODIS time-series smoothed prior",
        short_label="MODIS smooth",
        family="Surface-prior follow-up",
        directory="phaseD_results_masked_modis_prior_modis_timeseries_smooth_20260711",
        surface_prior="MCD43 Whittaker time-series smoothing",
        rt_backend="libRadtran LUT",
        tcwv="Atmospheric prior",
        masking="Cloud and water masked",
        proper_mask=True,
        controlled_note="Direct temporal-prior comparison with MODIS single-day.",
    ),
    MethodSpec(
        key="swir_nir_anchored",
        label="S2 SWIR+NIR-anchored prior",
        short_label="SWIR+NIR anchor",
        family="Surface-prior follow-up",
        directory="phaseD_results_masked_modis_prior_swir_nir_anchored_20260711",
        surface_prior="S2 L2A monthly predictor with SWIR+NIR anchor iteration",
        rt_backend="libRadtran LUT",
        tcwv="Atmospheric prior",
        masking="Cloud and water masked",
        proper_mask=True,
        controlled_note="Uses the refactored anchored surface-prior path; not a one-factor change versus MODIS priors.",
    ),
    MethodSpec(
        key="lut_prior_masked",
        label="LUT + prior TCWV, masked",
        short_label="LUT masked",
        family="Controlled factorial",
        directory="phaseD_results_masked_r2_factorial_lut_prior_masked_20260711_v2",
        surface_prior="S2 L2A monthly predictor with SWIR+NIR anchor iteration and tau prior",
        rt_backend="libRadtran LUT",
        tcwv="Atmospheric prior",
        masking="Cloud and water masked",
        proper_mask=True,
        controlled_note="Factorial reference arm.",
    ),
    MethodSpec(
        key="lut_prior_cloud",
        label="LUT + prior TCWV, cloud allowed",
        short_label="LUT cloud",
        family="Controlled factorial",
        directory="phaseD_results_masked_r2_factorial_lut_prior_cloud_20260711_v2",
        surface_prior="S2 L2A monthly predictor with SWIR+NIR anchor iteration and tau prior",
        rt_backend="libRadtran LUT",
        tcwv="Atmospheric prior",
        masking="Cloud retrieval allowed; water masked",
        proper_mask=False,
        controlled_note="Cloud-mask ablation versus the factorial reference.",
    ),
    MethodSpec(
        key="lut_prior_legacy",
        label="LUT + prior TCWV, legacy mask bypass",
        short_label="LUT legacy",
        family="Controlled factorial",
        directory="phaseD_results_masked_r2_factorial_lut_prior_legacy_20260711_v2",
        surface_prior="S2 L2A monthly predictor with SWIR+NIR anchor iteration and tau prior",
        rt_backend="libRadtran LUT",
        tcwv="Atmospheric prior",
        masking="Cloud and water masks bypassed",
        proper_mask=False,
        controlled_note="Legacy cloud-and-water mask ablation versus the factorial reference.",
    ),
    MethodSpec(
        key="lut_l2awvp_masked",
        label="LUT + S2 L2A TCWV, masked",
        short_label="LUT L2A WVP",
        family="Controlled factorial",
        directory="phaseD_results_masked_r2_factorial_lut_l2awvp_masked_20260711_v2",
        surface_prior="S2 L2A monthly predictor with SWIR+NIR anchor iteration and tau prior",
        rt_backend="libRadtran LUT",
        tcwv="Same-scene S2 L2A WVP",
        masking="Cloud and water masked",
        proper_mask=True,
        controlled_note="TCWV-source ablation versus the factorial reference.",
    ),
    MethodSpec(
        key="sixs_prior_masked",
        label="6S + prior TCWV, masked",
        short_label="6S masked",
        family="Controlled factorial",
        directory="phaseD_results_masked_r2_factorial_sixs_prior_masked_20260711_v2",
        surface_prior="S2 L2A monthly predictor with SWIR+NIR anchor iteration and tau prior",
        rt_backend="Native 6S, continental aerosol",
        tcwv="Atmospheric prior",
        masking="Cloud and water masked",
        proper_mask=True,
        controlled_note="RT-backend ablation versus the factorial reference.",
    ),
)


def _finite(value: object) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return math.nan
    return number if math.isfinite(number) else math.nan


def _valid(record: dict[str, Any] | None) -> bool:
    return bool(
        record
        and str(record.get("status", "")).upper() == "OK"
        and math.isfinite(_finite(record.get("truth")))
        and math.isfinite(_finite(record.get("retrieved")))
    )


def _hit(record: dict[str, Any] | None) -> bool:
    if not _valid(record):
        return False
    truth = _finite(record.get("truth"))
    error = _finite(record.get("err"))
    if not math.isfinite(error):
        error = _finite(record.get("retrieved")) - truth
    threshold = _finite(record.get("ee_threshold"))
    if not math.isfinite(threshold):
        threshold = 0.05 + 0.15 * truth
    return abs(error) <= threshold + 1e-12


def _error(record: dict[str, Any]) -> float:
    error = _finite(record.get("err"))
    if math.isfinite(error):
        return error
    return _finite(record.get("retrieved")) - _finite(record.get("truth"))


def _wilson(
    successes: int, count: int, z: float = 1.959963984540054
) -> tuple[float | None, float | None]:
    if count <= 0:
        return None, None
    rate = successes / count
    denominator = 1.0 + z * z / count
    centre = (rate + z * z / (2.0 * count)) / denominator
    half = z * math.sqrt(rate * (1.0 - rate) / count + z * z / (4.0 * count * count)) / denominator
    return max(0.0, centre - half), min(1.0, centre + half)


def _quantile(values: np.ndarray, q: float) -> float | None:
    return float(np.quantile(values, q)) if values.size else None


def _correlation(truth: np.ndarray, retrieved: np.ndarray) -> float | None:
    if truth.size < 3 or np.std(truth) == 0.0 or np.std(retrieved) == 0.0:
        return None
    return float(np.corrcoef(truth, retrieved)[0, 1])


def _summary(
    spec: MethodSpec,
    records: dict[str, dict[str, Any]],
    mids: list[str],
) -> dict[str, Any]:
    selected = [records[mid] for mid in mids if mid in records]
    valid = [row for row in selected if _valid(row)]
    errors = np.asarray([_error(row) for row in valid], dtype=np.float64)
    absolute = np.abs(errors)
    truth = np.asarray([_finite(row.get("truth")) for row in valid], dtype=np.float64)
    retrieved = np.asarray([_finite(row.get("retrieved")) for row in valid], dtype=np.float64)
    hits = sum(_hit(row) for row in selected)
    runtimes = np.asarray(
        [
            value
            for row in selected
            if math.isfinite(value := _finite(row.get("total_s", row.get("runtime_s"))))
        ],
        dtype=np.float64,
    )
    strict_ci = _wilson(hits, len(mids))
    coverage_ci = _wilson(len(valid), len(mids))
    valid_ci = _wilson(hits, len(valid))
    slope = intercept = None
    if truth.size >= 2 and np.std(truth) > 0.0:
        slope_value, intercept_value = np.polyfit(truth, retrieved, 1)
        slope, intercept = float(slope_value), float(intercept_value)
    statuses: dict[str, int] = {}
    for mid in mids:
        status = (
            "MISSING"
            if mid not in records
            else str(records[mid].get("status", "MISSING_STATUS")).upper()
        )
        statuses[status] = statuses.get(status, 0) + 1
    return {
        "key": spec.key,
        "method": spec.label,
        "short_method": spec.short_label,
        "family": spec.family,
        "expected": len(mids),
        "present": len(selected),
        "valid": len(valid),
        "hits": hits,
        "failed": statuses.get("FAILED", 0),
        "no_valid": statuses.get("NO_VALID_OBSERVATION", 0),
        "missing": statuses.get("MISSING", 0),
        "coverage": len(valid) / len(mids) if mids else None,
        "coverage_ci_low": coverage_ci[0],
        "coverage_ci_high": coverage_ci[1],
        "strict_rate": hits / len(mids) if mids else None,
        "strict_ci_low": strict_ci[0],
        "strict_ci_high": strict_ci[1],
        "valid_hit_rate": hits / len(valid) if valid else None,
        "valid_hit_ci_low": valid_ci[0],
        "valid_hit_ci_high": valid_ci[1],
        "bias": float(np.mean(errors)) if errors.size else None,
        "mae": float(np.mean(absolute)) if absolute.size else None,
        "rmse": float(np.sqrt(np.mean(np.square(errors)))) if errors.size else None,
        "median_abs_error": float(np.median(absolute)) if absolute.size else None,
        "p90_abs_error": _quantile(absolute, 0.9),
        "p95_abs_error": _quantile(absolute, 0.95),
        "underestimate_fraction": float(np.mean(errors < 0.0)) if errors.size else None,
        "pearson_r": _correlation(truth, retrieved),
        "regression_slope": slope,
        "regression_intercept": intercept,
        "median_runtime_s": float(np.median(runtimes)) if runtimes.size else None,
        "p90_runtime_s": _quantile(runtimes, 0.9),
        "statuses": statuses,
    }


def _bootstrap_mean_ci(
    values: np.ndarray, seed: int, draws: int = 20_000
) -> tuple[float | None, float | None]:
    if values.size == 0:
        return None, None
    rng = np.random.default_rng(seed)
    sample_indices = rng.integers(0, values.size, size=(draws, values.size))
    means = values[sample_indices].mean(axis=1)
    low, high = np.quantile(means, [0.025, 0.975])
    return float(low), float(high)


def _paired_summary(
    group: str,
    baseline: MethodSpec,
    candidate: MethodSpec,
    records: dict[str, dict[str, dict[str, Any]]],
    mids: list[str],
    comparison: str,
) -> dict[str, Any]:
    baseline_records = records[baseline.key]
    candidate_records = records[candidate.key]
    gains = [
        mid
        for mid in mids
        if not _hit(baseline_records.get(mid)) and _hit(candidate_records.get(mid))
    ]
    losses = [
        mid
        for mid in mids
        if _hit(baseline_records.get(mid)) and not _hit(candidate_records.get(mid))
    ]
    common = [
        mid
        for mid in mids
        if _valid(baseline_records.get(mid)) and _valid(candidate_records.get(mid))
    ]
    baseline_abs = np.asarray(
        [abs(_error(baseline_records[mid])) for mid in common], dtype=np.float64
    )
    candidate_abs = np.asarray(
        [abs(_error(candidate_records[mid])) for mid in common], dtype=np.float64
    )
    delta_abs = candidate_abs - baseline_abs
    seed = sum(ord(char) for char in f"{baseline.key}:{candidate.key}") + 20260711
    ci_low, ci_high = _bootstrap_mean_ci(delta_abs, seed)
    discordant = len(gains) + len(losses)
    mcnemar_p = (
        float(binomtest(min(len(gains), len(losses)), discordant, 0.5).pvalue)
        if discordant
        else 1.0
    )
    baseline_hits = sum(_hit(baseline_records.get(mid)) for mid in mids)
    candidate_hits = sum(_hit(candidate_records.get(mid)) for mid in mids)
    baseline_valid = sum(_valid(baseline_records.get(mid)) for mid in mids)
    candidate_valid = sum(_valid(candidate_records.get(mid)) for mid in mids)
    return {
        "group": group,
        "comparison": comparison,
        "baseline_key": baseline.key,
        "baseline": baseline.short_label,
        "candidate_key": candidate.key,
        "candidate": candidate.short_label,
        "expected": len(mids),
        "baseline_hits": baseline_hits,
        "candidate_hits": candidate_hits,
        "strict_hit_delta": candidate_hits - baseline_hits,
        "baseline_valid": baseline_valid,
        "candidate_valid": candidate_valid,
        "coverage_delta": (candidate_valid - baseline_valid) / len(mids),
        "gains": len(gains),
        "losses": len(losses),
        "net_hits": len(gains) - len(losses),
        "mcnemar_exact_p": mcnemar_p,
        "common_valid": len(common),
        "baseline_common_hits": sum(_hit(baseline_records[mid]) for mid in common),
        "candidate_common_hits": sum(_hit(candidate_records[mid]) for mid in common),
        "baseline_common_mae": float(np.mean(baseline_abs)) if baseline_abs.size else None,
        "candidate_common_mae": float(np.mean(candidate_abs)) if candidate_abs.size else None,
        "mae_delta": float(np.mean(delta_abs)) if delta_abs.size else None,
        "mae_delta_ci_low": ci_low,
        "mae_delta_ci_high": ci_high,
        "gain_ids": gains,
        "loss_ids": losses,
    }


def _load_directory(directory: Path) -> dict[str, dict[str, Any]]:
    records: dict[str, dict[str, Any]] = {}
    for path in sorted(directory.glob("*.json")):
        record = json.loads(path.read_text(encoding="utf-8"))
        records[str(record.get("matchup_id") or path.stem)] = record
    return records


def _truth_bin(value: float) -> str:
    if value < 0.2:
        return "clean"
    if value < 0.6:
        return "moderate"
    if value < 1.0:
        return "high"
    return "very_high"


def _comparison_outcome(
    old_record: dict[str, Any] | None, new_record: dict[str, Any] | None
) -> str:
    old_valid = _valid(old_record)
    new_valid = _valid(new_record)
    old_hit = _hit(old_record)
    new_hit = _hit(new_record)
    if old_valid and new_valid:
        if not old_hit and new_hit:
            return "fix"
        if old_hit and not new_hit:
            return "break"
        return "retained_hit" if old_hit else "retained_miss"
    if old_valid:
        return "new_abstention"
    if new_valid:
        return "new_only"
    return "neither_valid"


def _rebuild_diagnostic_selection(
    full_mids: list[str],
    old_records: dict[str, dict[str, Any]],
    new_records: dict[str, dict[str, Any]],
) -> tuple[list[str], dict[str, str]]:
    outcomes = {
        mid: _comparison_outcome(old_records.get(mid), new_records.get(mid)) for mid in full_mids
    }
    selected = [mid for mid in full_mids if outcomes[mid] in {"fix", "break", "new_abstention"}]
    selected_set = set(selected)
    for outcome in ("retained_hit", "retained_miss"):
        for truth_bin in ("clean", "moderate", "high", "very_high"):
            controls = [
                mid
                for mid in full_mids
                if outcomes[mid] == outcome
                and _truth_bin(_finite((new_records.get(mid) or old_records[mid]).get("truth")))
                == truth_bin
                and mid not in selected_set
            ][:3]
            selected.extend(controls)
            selected_set.update(controls)
    return selected, outcomes


def _cohort_score_row(
    *,
    cohort: str,
    series: str,
    ids: list[str],
    records: dict[str, dict[str, Any]],
) -> dict[str, Any]:
    hits = sum(_hit(records.get(mid)) for mid in ids)
    valid = sum(_valid(records.get(mid)) for mid in ids)
    return {
        "cohort": cohort,
        "series": series,
        "expected": len(ids),
        "valid": valid,
        "hits": hits,
        "strict_rate": hits / len(ids) if ids else None,
        "coverage": valid / len(ids) if ids else None,
        "valid_hit_rate": hits / valid if valid else None,
    }


def _full_campaign_context(
    *,
    full_mids: list[str],
    diagnostic_mids: list[str],
    old_records: dict[str, dict[str, Any]],
    new_records: dict[str, dict[str, Any]],
) -> tuple[list[dict[str, Any]], dict[str, Any], list[str]]:
    rebuilt, outcomes = _rebuild_diagnostic_selection(full_mids, old_records, new_records)
    full_low_cloud = [
        mid
        for mid in full_mids
        if math.isfinite(cloud := _finite((new_records.get(mid) or {}).get("cloud_frac")))
        and cloud < 0.2
    ]
    diagnostic_low_cloud = [mid for mid in diagnostic_mids if mid in set(full_low_cloud)]
    contexts = [
        ("Full 250, all", full_mids),
        ("Full 250, cloud <20%", full_low_cloud),
        ("Diagnostic 66, all", diagnostic_mids),
        ("Diagnostic 66, cloud <20%", diagnostic_low_cloud),
    ]
    rows = [
        _cohort_score_row(cohort=cohort, series=series, ids=ids, records=records)
        for cohort, ids in contexts
        for series, records in (
            ("Historical R2", old_records),
            ("Masked R2", new_records),
        )
    ]
    selected_outcomes: dict[str, int] = {}
    for mid in diagnostic_mids:
        outcome = outcomes[mid]
        selected_outcomes[outcome] = selected_outcomes.get(outcome, 0) + 1
    full_truth_bins: dict[str, int] = {}
    diagnostic_truth_bins: dict[str, int] = {}
    for ids, target in (
        (full_mids, full_truth_bins),
        (diagnostic_mids, diagnostic_truth_bins),
    ):
        for mid in ids:
            record = new_records.get(mid) or old_records[mid]
            label = _truth_bin(_finite(record.get("truth")))
            target[label] = target.get(label, 0) + 1
    selection = {
        "definition": (
            "All fixes, breaks, and new abstentions, plus the first three retained hits "
            "and retained misses in each of four AOD truth bins."
        ),
        "full_campaign_count": len(full_mids),
        "diagnostic_count": len(diagnostic_mids),
        "affected_count": sum(
            selected_outcomes.get(name, 0) for name in ("fix", "break", "new_abstention")
        ),
        "control_count": sum(
            selected_outcomes.get(name, 0) for name in ("retained_hit", "retained_miss")
        ),
        "outcomes": selected_outcomes,
        "full_low_cloud_count": len(full_low_cloud),
        "diagnostic_low_cloud_count": len(diagnostic_low_cloud),
        "full_truth_bins": full_truth_bins,
        "diagnostic_truth_bins": diagnostic_truth_bins,
        "very_high_fraction_full": full_truth_bins.get("very_high", 0) / len(full_mids),
        "very_high_fraction_diagnostic": diagnostic_truth_bins.get("very_high", 0)
        / len(diagnostic_mids),
        "rebuilt_selection_matches_manifest": rebuilt == diagnostic_mids,
    }
    return rows, selection, diagnostic_low_cloud


def _load_records(
    root: Path,
    mids: list[str],
) -> tuple[dict[str, dict[str, dict[str, Any]]], dict[str, Any]]:
    expected = set(mids)
    all_records: dict[str, dict[str, dict[str, Any]]] = {}
    audit_methods: list[dict[str, Any]] = []
    reference_truth: dict[str, float] = {}
    reference_regime: dict[str, str] = {}
    issues: list[str] = []
    for method_index, spec in enumerate(METHODS):
        directory = root / spec.directory
        records: dict[str, dict[str, Any]] = {}
        malformed: list[str] = []
        duplicate_ids: list[str] = []
        for path in sorted(directory.glob("*.json")):
            try:
                record = json.loads(path.read_text(encoding="utf-8"))
            except (OSError, ValueError):
                malformed.append(path.name)
                continue
            matchup_id = str(record.get("matchup_id") or path.stem)
            if matchup_id in records:
                duplicate_ids.append(matchup_id)
            records[matchup_id] = record
        actual = set(records)
        missing = sorted(expected - actual)
        extra = sorted(actual - expected)
        status_counts: dict[str, int] = {}
        formula_mismatches: list[str] = []
        error_mismatches: list[str] = []
        nonfinite_ok: list[str] = []
        truth_mismatches: list[str] = []
        regime_mismatches: list[str] = []
        for mid in mids:
            record = records.get(mid)
            if record is None:
                continue
            status = str(record.get("status", "MISSING_STATUS")).upper()
            status_counts[status] = status_counts.get(status, 0) + 1
            truth = _finite(record.get("truth"))
            regime = str(record.get("regime", ""))
            if method_index == 0:
                reference_truth[mid] = truth
                reference_regime[mid] = regime
            else:
                if not math.isclose(truth, reference_truth[mid], rel_tol=0.0, abs_tol=1e-12):
                    truth_mismatches.append(mid)
                if regime != reference_regime[mid]:
                    regime_mismatches.append(mid)
            if status == "OK":
                retrieved = _finite(record.get("retrieved"))
                error = _finite(record.get("err"))
                if not (math.isfinite(truth) and math.isfinite(retrieved) and math.isfinite(error)):
                    nonfinite_ok.append(mid)
                    continue
                if not math.isclose(error, retrieved - truth, rel_tol=0.0, abs_tol=1e-10):
                    error_mismatches.append(mid)
                explicit = record.get("within_ee")
                if isinstance(explicit, bool) and explicit != _hit(record):
                    formula_mismatches.append(mid)
        method_audit = {
            "key": spec.key,
            "directory": spec.directory,
            "json_files": len(list(directory.glob("*.json"))),
            "records": len(records),
            "missing": missing,
            "extra": extra,
            "malformed": malformed,
            "duplicates": duplicate_ids,
            "status_counts": status_counts,
            "nonfinite_ok": nonfinite_ok,
            "error_mismatches": error_mismatches,
            "within_ee_mismatches": formula_mismatches,
            "truth_mismatches": truth_mismatches,
            "regime_mismatches": regime_mismatches,
        }
        audit_methods.append(method_audit)
        for field in (
            "missing",
            "extra",
            "malformed",
            "duplicates",
            "nonfinite_ok",
            "error_mismatches",
            "within_ee_mismatches",
            "truth_mismatches",
            "regime_mismatches",
        ):
            if method_audit[field]:
                issues.append(f"{spec.key}:{field}={len(method_audit[field])}")
        all_records[spec.key] = records
    return all_records, {
        "expected_matchups": len(mids),
        "methods": audit_methods,
        "issues": issues,
        "ready": not issues,
    }


def _write_csv(path: Path, rows: list[dict[str, Any]], fields: Iterable[str] | None = None) -> None:
    if fields is None:
        fields = rows[0].keys() if rows else ()
    fieldnames = list(fields)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _widget_rows(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Keep the bounded artifact snapshot to scalar reader cell values."""
    return [
        {
            key: value
            for key, value in row.items()
            if value is None or isinstance(value, (str, int, float, bool))
        }
        for row in rows
    ]


def _format_pct(value: float | None) -> str:
    return "-" if value is None else f"{100.0 * value:.1f}%"


def _format_num(value: float | None, digits: int = 3) -> str:
    return "-" if value is None else f"{value:.{digits}f}"


def _artifact(
    *,
    generated_at: str,
    audit: dict[str, Any],
    method_summary: list[dict[str, Any]],
    common_summary: list[dict[str, Any]],
    regime_summary: list[dict[str, Any]],
    paired_summary: list[dict[str, Any]],
    outliers: list[dict[str, Any]],
    method_configs: list[dict[str, Any]],
    common_count: int,
    cohort_context: list[dict[str, Any]],
    selection_context: dict[str, Any],
    low_cloud_count: int,
) -> dict[str, Any]:
    title = "AOD diagnostic ablation: 66 selected matchups"
    by_key = {row["key"]: row for row in method_summary}
    best_strict = max(method_summary, key=lambda row: row["strict_rate"])
    best_masked = max(
        (
            row
            for row in method_summary
            if next(spec for spec in METHODS if spec.key == row["key"]).proper_mask
        ),
        key=lambda row: row["strict_rate"],
    )
    best_low_cloud = max(method_summary, key=lambda row: row["low_cloud_strict_rate"])
    best_low_cloud_masked = max(
        (
            row
            for row in method_summary
            if next(spec for spec in METHODS if spec.key == row["key"]).proper_mask
        ),
        key=lambda row: row["low_cloud_strict_rate"],
    )
    lowest_mae = min(method_summary, key=lambda row: row["mae"])
    common_best = max(common_summary, key=lambda row: row["valid_hit_rate"])
    high_rows = [row for row in regime_summary if row["regime"] == "high"]
    high_best = max(high_rows, key=lambda row: row["strict_rate"])
    total_records = sum(row["present"] for row in method_summary)
    total_valid = sum(row["valid"] for row in method_summary)
    total_no_valid = sum(row["no_valid"] for row in method_summary)
    total_failed = sum(row["failed"] for row in method_summary)
    total_missing = sum(row["missing"] for row in method_summary)
    context = {(row["cohort"], row["series"]): row for row in cohort_context}
    full_old = context[("Full 250, all", "Historical R2")]
    full_old_low = context[("Full 250, cloud <20%", "Historical R2")]
    diagnostic_old = context[("Diagnostic 66, all", "Historical R2")]
    diagnostic_old_low = context[("Diagnostic 66, cloud <20%", "Historical R2")]
    masked_gap = selection_context["historical_vs_masked_reference"]

    factorial = [row for row in paired_summary if row["group"] == "Factorial"]
    surface = [row for row in paired_summary if row["group"] == "Surface prior"]
    cloud = next(row for row in factorial if row["candidate_key"] == "lut_prior_cloud")
    legacy = next(row for row in factorial if row["candidate_key"] == "lut_prior_legacy")
    l2awvp = next(row for row in factorial if row["candidate_key"] == "lut_l2awvp_masked")
    sixs = next(row for row in factorial if row["candidate_key"] == "sixs_prior_masked")
    smooth_vs_daily = next(
        row
        for row in surface
        if row["baseline_key"] == "modis_single_day"
        and row["candidate_key"] == "modis_timeseries_smooth"
    )
    anchor_vs_smooth = next(
        row
        for row in surface
        if row["baseline_key"] == "modis_timeseries_smooth"
        and row["candidate_key"] == "swir_nir_anchored"
    )

    source = {
        "id": "aod-performance-analysis",
        "label": "Completed AOD experiment performance analysis",
        "path": "report-data.csv",
        "query": {
            "engine": "DuckDB",
            "language": "sql",
            "sql": "SELECT * FROM read_csv_auto('report-data.csv', header = true);",
            "description": "Load the bounded reviewed rows produced by the deterministic matchup-level Python analysis.",
            "executed_at": generated_at,
            "tables_used": ["report-data.csv"],
            "filters": [
                "masked_r2_regression_diagnostic_mids.txt contains 66 selected diagnostic matchups, not 528 unique sites",
                "The 66 contain every full-campaign fix, break, and new abstention plus 24 AOD-stratified retained-hit/retained-miss controls",
                "A valid retrieval requires status=OK with finite truth and retrieved AOD",
                "Strict score counts no-valid, failed, malformed, and missing outcomes as misses",
                "Common-case metrics use the intersection of status-OK finite retrievals across all eight methods",
            ],
            "metric_definitions": [
                "Method-scene evaluation: one experiment arm applied to one matchup; 8 arms x 66 matchups = 528 evaluations",
                f"Within expected error: {EXPECTED_ERROR_FORMULA}",
                "Coverage: finite status-OK retrievals divided by 66 expected matchups",
                "Strict within-EE rate: within-EE hits divided by all 66 expected matchups",
                "Valid-only within-EE rate: within-EE hits divided by finite status-OK retrievals",
                "MAE/RMSE/bias: computed from retrieved minus AERONET truth on the stated valid cohort",
                "Rate intervals: two-sided 95% Wilson intervals",
                "Paired MAE intervals: 20,000 deterministic matchup bootstrap resamples",
                "McNemar p-values: two-sided exact binomial tests over discordant strict hit outcomes",
            ],
        },
    }

    headline = [
        {
            "unique_matchups": audit["expected_matchups"],
            "experiment_arms": len(METHODS),
            "method_scene_evaluations": total_records,
            "execution_failures": total_failed,
            "records_complete": total_records,
            "records_expected": len(METHODS) * audit["expected_matchups"],
            "missing_records": total_missing,
            "best_strict_rate": best_strict["strict_rate"],
            "best_masked_rate": best_masked["strict_rate"],
            "common_count": common_count,
        }
    ]
    coverage_accuracy = [
        {
            "method": row["short_method"],
            "metric": metric,
            "value": row[field],
            "valid_hit_rate": row["valid_hit_rate"],
            "hits": row["hits"],
            "valid": row["valid"],
            "expected": row["expected"],
        }
        for row in method_summary
        for metric, field in (
            ("Strict within EE", "strict_rate"),
            ("Retrieval coverage", "coverage"),
        )
    ]
    error_summary = [
        {
            "method": row["short_method"],
            "metric": metric,
            "value": row[field],
            "rmse": row["rmse"],
            "valid": row["valid"],
        }
        for row in method_summary
        for metric, field in (
            ("Mean absolute error", "mae"),
            ("Median absolute error", "median_abs_error"),
        )
    ]
    common_chart = [
        {
            "method": row["short_method"],
            "within_ee_rate": row["valid_hit_rate"],
            "hits": row["hits"],
            "common_valid": row["valid"],
            "mae": row["mae"],
        }
        for row in sorted(common_summary, key=lambda item: item["valid_hit_rate"], reverse=True)
    ]
    paired_chart = [
        {
            "comparison": row["candidate"],
            "group": row["group"],
            "net_hits": row["net_hits"],
            "gains": row["gains"],
            "losses": row["losses"],
            "common_valid": row["common_valid"],
            "mae_delta": row["mae_delta"],
        }
        for row in factorial
    ]
    runtime_summary = [
        {
            "method": row["short_method"],
            "median_runtime_s": row["median_runtime_s"],
            "p90_runtime_s": row["p90_runtime_s"],
        }
        for row in method_summary
    ]
    regime_chart = [
        {
            "method": row["short_method"],
            "regime": {"clean": "Clean", "mod": "Moderate", "high": "High"}.get(
                row["regime"], row["regime"]
            ),
            "strict_rate": row["strict_rate"],
            "coverage": row["coverage"],
            "hits": row["hits"],
            "expected": row["expected"],
        }
        for row in regime_summary
    ]

    cards = [
        {
            "id": "completion-card",
            "dataset": "headline",
            "sourceId": source["id"],
            "description": "The same 66 selected matchups were evaluated by eight experiment arms; 528 is not a unique-site count.",
            "metrics": [
                {"label": "Unique matchups", "field": "unique_matchups", "format": "number"},
                {"label": "Experiment arms", "field": "experiment_arms", "format": "number"},
                {
                    "label": "Method-scene evaluations",
                    "field": "method_scene_evaluations",
                    "format": "number",
                },
                {"label": "Execution failures", "field": "execution_failures", "format": "number"},
            ],
        },
        {
            "id": "best-strict-card",
            "dataset": "headline",
            "sourceId": source["id"],
            "description": f"Highest raw 66-scene score: {best_strict['method']}; this arm bypasses cloud and water masking.",
            "metrics": [
                {"label": "Best strict score", "field": "best_strict_rate", "format": "percent"},
            ],
        },
        {
            "id": "best-masked-card",
            "dataset": "headline",
            "sourceId": source["id"],
            "description": f"Highest strict score among arms retaining cloud and water masking: {best_masked['method']}.",
            "metrics": [
                {"label": "Best masked score", "field": "best_masked_rate", "format": "percent"},
            ],
        },
        {
            "id": "common-card",
            "dataset": "headline",
            "sourceId": source["id"],
            "description": "Matchups with finite status-OK retrievals from every method, used for coverage-neutral comparisons.",
            "metrics": [
                {"label": "All-method common cases", "field": "common_count", "format": "number"},
            ],
        },
    ]

    charts = [
        {
            "id": "cohort-context-chart",
            "title": "Full-campaign and diagnostic-cohort within-EE rates",
            "subtitle": "Historical and masked R2 on all scenes and on a fixed cloud-below-20% subset",
            "type": "bar",
            "intent": "comparison",
            "question": "Why are the 66-scene experiment rates lower than the previous full-250 validation?",
            "rationale": "A grouped comparison shows the effect of diagnostic selection separately from the retrieval configuration.",
            "comparisonContext": {
                "baseline": "Historical R2",
                "denominator": "Shown separately for each cohort",
                "grain": "cohort and R2 configuration",
                "unit": "strict within-EE fraction",
            },
            "dataset": "cohort_context",
            "sourceId": source["id"],
            "encodings": {
                "x": {"field": "cohort", "type": "nominal", "label": "Evaluation cohort"},
                "y": {
                    "field": "strict_rate",
                    "type": "quantitative",
                    "label": "Strict within EE",
                },
                "color": {"field": "series", "type": "nominal", "label": "Configuration"},
                "tooltip": [
                    {"field": "hits", "type": "quantitative"},
                    {"field": "expected", "type": "quantitative"},
                    {"field": "coverage", "type": "quantitative", "format": "percent"},
                ],
            },
            "valueFormat": "percent",
            "layout": "full",
            "palette": {"kind": "categorical"},
            "legend": {"position": "bottom", "sort": "spec"},
            "settings": {
                "groupMode": "grouped",
                "categoryLabelPolicy": "wrap",
                "showValues": True,
            },
        },
        {
            "id": "coverage-accuracy-chart",
            "title": "Strict within-EE score and retrieval coverage",
            "subtitle": "Fixed 66-scene denominator; strict score counts every abstention as a miss",
            "type": "bar",
            "intent": "comparison",
            "question": "How much of each method's strict score is influenced by retrieval coverage?",
            "rationale": "The two same-scale rates expose the quality-coverage trade-off without hiding abstentions.",
            "comparisonContext": {
                "baseline": "66 expected matchups",
                "denominator": "66 scenes per method",
                "grain": "method",
                "unit": "fraction",
            },
            "dataset": "coverage_accuracy",
            "sourceId": source["id"],
            "encodings": {
                "x": {"field": "method", "type": "nominal", "label": "Method"},
                "y": {"field": "value", "type": "quantitative", "label": "Rate"},
                "color": {"field": "metric", "type": "nominal", "label": "Metric"},
                "tooltip": [
                    {"field": "hits", "type": "quantitative"},
                    {"field": "valid", "type": "quantitative"},
                    {"field": "expected", "type": "quantitative"},
                    {
                        "field": "valid_hit_rate",
                        "type": "quantitative",
                        "format": "percent",
                    },
                ],
            },
            "valueFormat": "percent",
            "layout": "full",
            "palette": {"kind": "categorical"},
            "legend": {"position": "bottom", "sort": "spec"},
            "settings": {"groupMode": "grouped", "categoryLabelPolicy": "wrap", "showValues": True},
        },
        {
            "id": "common-score-chart",
            "title": "Within-EE rate on the all-method common-valid cohort",
            "subtitle": f"Same {common_count} matchups for every method; coverage differences removed",
            "type": "horizontalBar",
            "intent": "comparison",
            "question": "Which method scores best when every method retrieved the same matchups?",
            "rationale": "A ranked horizontal bar supports long method labels and a common denominator.",
            "comparisonContext": {
                "baseline": "All-method common-valid cohort",
                "denominator": f"{common_count} common matchups",
                "grain": "method",
                "unit": "within-EE fraction",
            },
            "dataset": "common_summary",
            "sourceId": source["id"],
            "encodings": {
                "x": {"field": "method", "type": "nominal", "label": "Method"},
                "y": {
                    "field": "within_ee_rate",
                    "type": "quantitative",
                    "label": "Within EE",
                },
                "tooltip": [
                    {"field": "hits", "type": "quantitative"},
                    {"field": "common_valid", "type": "quantitative"},
                    {"field": "mae", "type": "quantitative"},
                ],
            },
            "valueFormat": "percent",
            "layout": "full",
            "palette": {"kind": "sequential"},
            "settings": {"sort": "descending", "showValues": True, "categoryLabelPolicy": "wrap"},
        },
        {
            "id": "error-chart",
            "title": "Mean and median absolute AOD error",
            "subtitle": "Available status-OK retrievals; the mean reveals sensitivity to high-error outliers",
            "type": "bar",
            "intent": "comparison",
            "question": "How do typical error and outlier-sensitive error differ across methods?",
            "rationale": "Mean and median absolute error share units and expose skew without a misleading aggregate rank.",
            "comparisonContext": {
                "baseline": "AERONET AOD550",
                "denominator": "Method-specific status-OK retrievals",
                "grain": "method",
                "unit": "AOD",
            },
            "dataset": "error_summary",
            "sourceId": source["id"],
            "encodings": {
                "x": {"field": "method", "type": "nominal", "label": "Method"},
                "y": {
                    "field": "value",
                    "type": "quantitative",
                    "label": "Absolute AOD error",
                },
                "color": {"field": "metric", "type": "nominal", "label": "Metric"},
                "tooltip": [
                    {"field": "valid", "type": "quantitative"},
                    {"field": "rmse", "type": "quantitative"},
                ],
            },
            "valueFormat": "number",
            "unit": "AOD",
            "layout": "full",
            "palette": {"kind": "categorical"},
            "legend": {"position": "bottom", "sort": "spec"},
            "settings": {"groupMode": "grouped", "categoryLabelPolicy": "wrap", "showValues": True},
        },
        {
            "id": "factorial-net-chart",
            "title": "Factorial strict-hit changes versus LUT masked",
            "subtitle": "Paired gains minus losses over all 66 matchups; positive values can include recovered coverage",
            "type": "horizontalBar",
            "intent": "comparison",
            "question": "Which controlled factorial changes add or lose strict within-EE hits?",
            "rationale": "Signed paired net hits summarize each ablation while the adjacent table retains gains, losses, and common-case MAE.",
            "comparisonContext": {
                "baseline": "LUT + prior TCWV, masked",
                "denominator": "66 paired matchups",
                "grain": "factorial arm",
                "unit": "net within-EE hits",
            },
            "dataset": "factorial_net",
            "sourceId": source["id"],
            "encodings": {
                "x": {"field": "comparison", "type": "nominal", "label": "Candidate"},
                "y": {
                    "field": "net_hits",
                    "type": "quantitative",
                    "label": "Net strict hits",
                },
                "tooltip": [
                    {"field": "gains", "type": "quantitative"},
                    {"field": "losses", "type": "quantitative"},
                    {"field": "mae_delta", "type": "quantitative"},
                ],
            },
            "valueFormat": "number",
            "unit": "hits",
            "layout": "full",
            "palette": {"kind": "diverging", "midpoint": 0},
            "referenceLines": [
                {"axis": "y", "value": 0, "label": "No net change", "color": "neutral"}
            ],
            "settings": {"sort": "descending", "showValues": True, "categoryLabelPolicy": "wrap"},
        },
        {
            "id": "regime-chart",
            "title": "Strict within-EE rate by AERONET AOD regime",
            "subtitle": "Fixed regime denominators: 11 clean, 16 moderate, and 39 high-AOD matchups",
            "type": "heatmap",
            "intent": "relationship",
            "question": "Where in the AOD range does each method succeed or fail?",
            "rationale": "A compact method-by-regime matrix exposes systematic regime weakness without eight crowded grouped series.",
            "comparisonContext": {
                "baseline": "AERONET regime labels",
                "denominator": "Expected matchups in each regime",
                "grain": "method by regime",
                "unit": "strict within-EE fraction",
            },
            "dataset": "regime_summary",
            "sourceId": source["id"],
            "encodings": {
                "x": {"field": "method", "type": "nominal", "label": "Method"},
                "y": {
                    "field": "strict_rate",
                    "type": "quantitative",
                    "label": "Strict within EE",
                },
                "color": {"field": "regime", "type": "nominal", "label": "Regime"},
                "tooltip": [
                    {"field": "hits", "type": "quantitative"},
                    {"field": "expected", "type": "quantitative"},
                    {"field": "coverage", "type": "quantitative", "format": "percent"},
                ],
            },
            "valueFormat": "percent",
            "layout": "full",
            "palette": {"kind": "sequential"},
        },
        {
            "id": "runtime-chart",
            "title": "Observed median runtime per matchup",
            "subtitle": "Wall-clock seconds recorded by the harness; cache state and cluster placement were not controlled",
            "type": "horizontalBar",
            "intent": "comparison",
            "question": "What runtime did each completed experiment record?",
            "rationale": "A horizontal bar supports exact method comparison while keeping the confounding caveat visible.",
            "comparisonContext": {
                "baseline": "Recorded harness total_s",
                "denominator": "All records with finite runtime",
                "grain": "method",
                "unit": "seconds per matchup",
            },
            "dataset": "runtime_summary",
            "sourceId": source["id"],
            "encodings": {
                "x": {"field": "method", "type": "nominal", "label": "Method"},
                "y": {
                    "field": "median_runtime_s",
                    "type": "quantitative",
                    "label": "Median seconds",
                },
                "tooltip": [{"field": "p90_runtime_s", "type": "quantitative"}],
            },
            "valueFormat": "number",
            "unit": "seconds",
            "layout": "full",
            "palette": {"kind": "sequential"},
            "settings": {"sort": "ascending", "showValues": True, "categoryLabelPolicy": "wrap"},
        },
    ]

    tables = [
        {
            "id": "method-summary-table",
            "title": "Complete method performance",
            "subtitle": f"Exact metrics for all eight arms; low-cloud columns use the same {low_cloud_count} selected matchups",
            "dataset": "method_summary",
            "sourceId": source["id"],
            "layout": "full",
            "density": "dense",
            "defaultSort": {"field": "strict_rate", "direction": "desc"},
            "columns": [
                {"field": "short_method", "label": "Method", "type": "text"},
                {"field": "hits", "label": "EE hits", "format": "number"},
                {"field": "strict_rate", "label": "Strict EE", "format": "percent"},
                {"field": "valid", "label": "Valid", "format": "number"},
                {"field": "coverage", "label": "Coverage", "format": "percent"},
                {"field": "valid_hit_rate", "label": "EE of valid", "format": "percent"},
                {"field": "low_cloud_hits", "label": "Low-cloud hits", "format": "number"},
                {
                    "field": "low_cloud_strict_rate",
                    "label": "Low-cloud EE",
                    "format": "percent",
                },
                {"field": "mae", "label": "MAE", "format": "number"},
                {"field": "rmse", "label": "RMSE", "format": "number"},
                {"field": "bias", "label": "Bias", "format": "number"},
                {"field": "pearson_r", "label": "Pearson r", "format": "number"},
                {"field": "median_runtime_s", "label": "Median s", "format": "number"},
            ],
        },
        {
            "id": "paired-table",
            "title": "Paired surface-prior and factorial comparisons",
            "subtitle": "Strict gains/losses include abstentions; paired MAE uses common-valid matchups only",
            "dataset": "paired_summary",
            "sourceId": source["id"],
            "layout": "full",
            "density": "dense",
            "defaultSort": {"field": "net_hits", "direction": "desc"},
            "columns": [
                {"field": "group", "label": "Group", "type": "text"},
                {"field": "comparison", "label": "Comparison", "type": "text"},
                {"field": "gains", "label": "Gains", "format": "number"},
                {"field": "losses", "label": "Losses", "format": "number"},
                {"field": "net_hits", "label": "Net hits", "format": "number", "movement": True},
                {
                    "field": "coverage_delta",
                    "label": "Coverage delta",
                    "format": "percent",
                    "movement": True,
                },
                {"field": "common_valid", "label": "Common valid", "format": "number"},
                {"field": "mae_delta", "label": "MAE delta", "format": "number", "movement": True},
                {"field": "mae_delta_ci_low", "label": "MAE CI low", "format": "number"},
                {"field": "mae_delta_ci_high", "label": "MAE CI high", "format": "number"},
                {"field": "mcnemar_exact_p", "label": "Exact p", "format": "number"},
            ],
        },
        {
            "id": "method-config-table",
            "title": "Experiment configuration map",
            "subtitle": "The factorial arms are controlled against LUT masked; the surface-prior family is analyzed separately",
            "dataset": "method_configs",
            "sourceId": source["id"],
            "layout": "full",
            "density": "spacious",
            "defaultSort": {"field": "short_label", "direction": "asc"},
            "columns": [
                {"field": "short_label", "label": "Method", "type": "text"},
                {"field": "family", "label": "Family", "type": "text"},
                {"field": "surface_prior", "label": "Surface prior", "type": "text"},
                {"field": "rt_backend", "label": "RT", "type": "text"},
                {"field": "tcwv", "label": "TCWV", "type": "text"},
                {"field": "masking", "label": "Masking", "type": "text"},
            ],
        },
        {
            "id": "outlier-table",
            "title": "Largest absolute matchup errors",
            "subtitle": "Top 16 method-scene pairs; SSE share shows each row's contribution to that method's RMSE",
            "dataset": "outliers",
            "sourceId": source["id"],
            "layout": "full",
            "density": "dense",
            "defaultSort": {"field": "abs_error", "direction": "desc"},
            "columns": [
                {"field": "short_method", "label": "Method", "type": "text"},
                {"field": "site", "label": "Site", "type": "text"},
                {"field": "regime", "label": "Regime", "type": "text"},
                {"field": "truth", "label": "Truth", "format": "number"},
                {"field": "retrieved", "label": "Retrieved", "format": "number"},
                {"field": "error", "label": "Error", "format": "number", "movement": True},
                {"field": "abs_error", "label": "Absolute error", "format": "number"},
                {"field": "sse_share", "label": "Method SSE share", "format": "percent"},
            ],
        },
    ]

    blocks = [
        {"id": "title", "type": "markdown", "body": f"# {title}"},
        {
            "id": "technical-summary",
            "type": "markdown",
            "sourceId": source["id"],
            "body": (
                "## Technical summary\n\n"
                f"**This report contains 66 unique matchups, not 528 sites.** The same 66 selected Sentinel-2/AERONET matchups were run through eight experiment arms, giving "
                f"{total_records} method-scene evaluations. All {total_records}/{len(METHODS) * audit['expected_matchups']} JSON records are present, "
                f"with {total_failed} `FAILED`, {total_missing} missing, and no malformed or formula-inconsistent records. "
                f"There are {total_valid} finite retrievals and {total_no_valid} valid abstention records across the eight arms.\n\n"
                f"**The 66-scene rates are not comparable to the previous full-250 validation rate.** Historical R2 scores {full_old['hits']}/250 "
                f"({_format_pct(full_old['strict_rate'])}) on the full campaign and {full_old_low['hits']}/{full_old_low['expected']} "
                f"({_format_pct(full_old_low['strict_rate'])}) for cloud below 20%. On the deliberately selected 66, the same historical R2 falls to "
                f"{diagnostic_old['hits']}/66 ({_format_pct(diagnostic_old['strict_rate'])}), and its low-cloud subset falls to "
                f"{diagnostic_old_low['hits']}/{diagnostic_old_low['expected']} ({_format_pct(diagnostic_old_low['strict_rate'])}).\n\n"
                f"**Within this diagnostic ablation only,** {best_strict['method']} has the highest strict score "
                f"at {best_strict['hits']}/66 ({_format_pct(best_strict['strict_rate'])}), but it achieves full coverage by bypassing cloud and water masks. "
                f"Among properly masked arms, {best_masked['method']} leads at {best_masked['hits']}/66 ({_format_pct(best_masked['strict_rate'])}) "
                f"and also has the lowest available-case MAE ({_format_num(lowest_mae['mae'])}).\n\n"
                f"**Coverage and matchup difficulty materially affect the ranking.** Only {common_count} scenes are valid for every method. "
                f"On that fixed subset, {common_best['method']} has the highest within-EE rate ({_format_pct(common_best['valid_hit_rate'])}); "
                "the report therefore keeps strict, valid-only, and common-case results separate."
            ),
        },
        {
            "id": "headline-metrics",
            "type": "metric-strip",
            "cardIds": ["completion-card", "best-strict-card", "best-masked-card", "common-card"],
        },
        {
            "id": "selection-finding",
            "type": "markdown",
            "sourceId": source["id"],
            "body": (
                "## The diagnostic selection, not a broken denominator, explains the low headline rate\n\n"
                f"The 66 include every changed/problem case from the full comparison: {selection_context['outcomes'].get('fix', 0)} fixes, "
                f"{selection_context['outcomes'].get('break', 0)} regressions, and {selection_context['outcomes'].get('new_abstention', 0)} new abstentions. "
                f"Only {selection_context['control_count']} are controls: 12 retained hits and 12 retained misses, stratified across four truth-AOD bins. "
                f"Very-high AOD (truth >=1) therefore rises from {_format_pct(selection_context['very_high_fraction_full'])} of the full 250 to "
                f"{_format_pct(selection_context['very_high_fraction_diagnostic'])} of the diagnostic 66. Only {low_cloud_count}/66 selected scenes have cloud below 20%. "
                "This cohort was designed to discriminate failure mechanisms, so its raw score is expected to be much lower than a campaign prevalence score."
            ),
        },
        {"id": "cohort-context-chart-block", "type": "chart", "chartId": "cohort-context-chart"},
        {
            "id": "coverage-finding",
            "type": "markdown",
            "sourceId": source["id"],
            "body": (
                "## Selection explains the baseline shift, but current arms still lag historical R2 on the same scenes\n\n"
                f"On the identical 66 matchups, historical R2 has {diagnostic_old['hits']} strict hits. The best current arm has {best_strict['hits']}, "
                f"and the best properly masked arm has {best_masked['hits']}. For the fixed {low_cloud_count}-scene low-cloud subset, historical R2 has "
                f"{diagnostic_old_low['hits']} hits; the best current arm has {best_low_cloud['low_cloud_hits']} and the best properly masked arm, "
                f"{best_low_cloud_masked['method']}, has {best_low_cloud_masked['low_cloud_hits']}. This remaining same-scene gap is real, but it compares "
                "different retrieval configurations rather than indicating missing or failed jobs. For the LUT masked reference, "
                f"{masked_gap['historical_hits_lost_to_masked_abstention']} of the {masked_gap['historical_hits'] - masked_gap['masked_reference_hits']} lost hits come from scenes where historical R2 hit but the masked arm abstains; "
                f"among {masked_gap['shared_valid']} shared-valid scenes, hits move from {masked_gap['historical_shared_hits']} to {masked_gap['masked_reference_shared_hits']} (net "
                f"{masked_gap['masked_reference_shared_hits'] - masked_gap['historical_shared_hits']:+d}).\n\n"
                "**Mask bypass raises raw hits, but does not establish better retrieval quality.** "
                f"The factorial reference returns {by_key['lut_prior_masked']['valid']}/66 retrievals and {by_key['lut_prior_masked']['hits']} strict hits. "
                f"Allowing cloud retrieval raises coverage to {by_key['lut_prior_cloud']['valid']}/66 and adds {cloud['net_hits']:+d} net hits; "
                f"bypassing both cloud and water masks reaches 66/66 coverage and adds {legacy['net_hits']:+d}. "
                f"However, available-case MAE rises from {_format_num(by_key['lut_prior_masked']['mae'])} to {_format_num(by_key['lut_prior_cloud']['mae'])} and "
                f"{_format_num(by_key['lut_prior_legacy']['mae'])}, respectively. The extra strict hits therefore mix broader coverage with noisier estimates."
            ),
        },
        {"id": "coverage-chart-block", "type": "chart", "chartId": "coverage-accuracy-chart"},
        {
            "id": "common-finding",
            "type": "markdown",
            "sourceId": source["id"],
            "body": (
                "## A common-case comparison narrows the apparent differences\n\n"
                f"Exactly {common_count} matchups have finite AOD from all eight methods. Holding that cohort fixed removes the advantage from retrieving more scenes. "
                f"{common_best['method']} ranks first on expected-error hits in this subset at {common_best['hits']}/{common_count} "
                f"({_format_pct(common_best['valid_hit_rate'])}), but the sample is small and failure-enriched. This is a robustness check, not a population score."
            ),
        },
        {"id": "common-chart-block", "type": "chart", "chartId": "common-score-chart"},
        {
            "id": "error-finding",
            "type": "markdown",
            "sourceId": source["id"],
            "body": (
                "## The masked LUT reference has the best error magnitude, while all methods under-retrieve on average\n\n"
                f"{lowest_mae['method']} has the lowest available-case MAE ({_format_num(lowest_mae['mae'])}) and RMSE ({_format_num(lowest_mae['rmse'])}). "
                "Every arm has negative mean bias, and the gap between mean and median absolute error shows that a small set of high-AOD misses dominates RMSE. "
                "Consequently, expected-error hit rate, MAE, median error, and outlier concentration should be read together."
            ),
        },
        {"id": "error-chart-block", "type": "chart", "chartId": "error-chart"},
        {
            "id": "surface-finding",
            "type": "markdown",
            "sourceId": source["id"],
            "body": (
                "## MODIS smoothing is only a one-hit change; the anchored path trades hits for lower mean error\n\n"
                f"Against the MODIS single-day prior, Whittaker smoothing produces {smooth_vs_daily['gains']} gains and {smooth_vs_daily['losses']} losses "
                f"for a {smooth_vs_daily['net_hits']:+d} strict-hit change. Relative to the smoothed MODIS prior, the S2 SWIR+NIR-anchored path changes "
                f"strict hits by {anchor_vs_smooth['net_hits']:+d} and common-case MAE by {_format_num(anchor_vs_smooth['mae_delta'])}. "
                "The anchored arm is not a one-factor comparison because it also switches to the L2A monthly predictor and anchor iteration; the result is descriptive, not causal."
            ),
        },
        {
            "id": "factorial-finding",
            "type": "markdown",
            "sourceId": source["id"],
            "body": (
                "## Neither L2A water vapour nor 6S improves the controlled masked reference\n\n"
                f"Changing only target TCWV to same-scene S2 L2A WVP yields {l2awvp['gains']} gains and {l2awvp['losses']} losses "
                f"({l2awvp['net_hits']:+d} net), while changing LUT RT to native 6S yields {sixs['gains']} gains and {sixs['losses']} losses "
                f"({sixs['net_hits']:+d} net). Their common-valid MAE changes are {_format_num(l2awvp['mae_delta'])} and {_format_num(sixs['mae_delta'])}. "
                "The paired bootstrap intervals and exact discordant-hit tests in the table quantify the uncertainty; small net changes should not be treated as decisive."
            ),
        },
        {"id": "factorial-chart-block", "type": "chart", "chartId": "factorial-net-chart"},
        {"id": "paired-table-block", "type": "table", "tableId": "paired-table"},
        {
            "id": "regime-finding",
            "type": "markdown",
            "sourceId": source["id"],
            "body": (
                "## High-AOD scenes remain the main weakness\n\n"
                f"The cohort contains 39 high-AOD scenes, compared with 11 clean and 16 moderate scenes. Even the strongest strict high-regime result is "
                f"{high_best['hits']}/{high_best['expected']} ({_format_pct(high_best['strict_rate'])}) from {high_best['method']}. "
                "This concentration is expected in a regression-diagnostic cohort, but it also shows that aggregate improvement depends on resolving high-AOD under-retrieval rather than only recovering low-risk abstentions."
            ),
        },
        {"id": "regime-chart-block", "type": "chart", "chartId": "regime-chart"},
        {
            "id": "runtime-finding",
            "type": "markdown",
            "sourceId": source["id"],
            "body": (
                "## Recorded runtime is operational context, not a controlled benchmark\n\n"
                "The result files contain end-to-end wall time, including data access, cache reuse, LUT setup, and cluster variability. The anchored method also performs two passes. "
                "These medians are useful for planning reruns, but differences must not be attributed to algorithmic cost until every arm is run under a shared cold-cache or warm-cache protocol."
            ),
        },
        {"id": "runtime-chart-block", "type": "chart", "chartId": "runtime-chart"},
        {
            "id": "scope-definitions",
            "type": "markdown",
            "body": (
                "## Scope, data, and metric definitions\n\n"
                "The analysis covers eight completed experiment arms over the same 66 unique Sentinel-2/AERONET matchups listed in `masked_r2_regression_diagnostic_mids.txt`; 528 is the resulting count of method-scene evaluations. "
                "The cohort contains all 42 fixes, regressions, and new abstentions from the historical-to-masked R2 comparison, plus three retained-hit and three retained-miss controls in each of four truth-AOD bins (24 controls). It is not a random sample of the full campaign. "
                f"A within-EE hit satisfies `{EXPECTED_ERROR_FORMULA}`. Strict score uses all 66 expected scenes. Coverage is the status-OK finite-retrieval share. "
                f"The fixed low-cloud diagnostic subset contains {low_cloud_count} scenes whose masked-R2 cloud fraction is below 20%. Valid-only metrics use each method's own available retrievals. "
                "Common-case metrics use only scenes retrieved by every method. AERONET truth and saved scene-mean retrievals are in AOD550 units."
            ),
        },
        {"id": "method-config-block", "type": "table", "tableId": "method-config-table"},
        {
            "id": "methodology",
            "type": "markdown",
            "sourceId": source["id"],
            "body": (
                "## Paired methodology and validation\n\n"
                "Records are joined by exact matchup identifier. The audit independently checks record completeness, JSON parsing, duplicate IDs, truth/regime agreement, finite status-OK values, `err = retrieved - truth`, and recomputation of the expected-error flag. "
                "The factorial uses LUT + atmospheric-prior TCWV + proper masking as its reference. Strict paired gains and losses include coverage changes; paired MAE uses only matchups valid in both arms. "
                "Two-sided Wilson intervals describe rates, deterministic 20,000-resample matchup bootstraps describe paired MAE differences, and exact McNemar/binomial tests describe discordant hit outcomes. "
                "No multiplicity-adjusted hypothesis claim is made."
            ),
        },
        {"id": "method-summary-block", "type": "table", "tableId": "method-summary-table"},
        {
            "id": "limitations",
            "type": "markdown",
            "body": (
                "## Limitations, uncertainty, and robustness checks\n\n"
                "The 66-scene set is outcome-selected and deliberately failure-enriched, so its rates are ablation diagnostics rather than expected production prevalence. The full-250 and cloud-below-20% historical R2 rates are shown only as cohort context, not as directly interchangeable denominators. "
                "The common-valid subset removes coverage bias but introduces its own selection toward easier shared scenes. "
                "The MODIS single-day versus smoothed comparison is close to one-factor; comparisons involving the SWIR+NIR-anchored arm change additional surface-prior configuration. "
                "The cloud and legacy bypass arms return values where the proper mask abstains, so higher strict hit counts do not by themselves establish physically trustworthy retrievals. "
                "Rate intervals are wide at n=66, paired tests are exploratory across several arms, and runtime was not benchmarked under controlled cache and host conditions."
            ),
        },
        {
            "id": "outlier-finding",
            "type": "markdown",
            "sourceId": source["id"],
            "body": (
                "## A few high-AOD misses dominate squared error\n\n"
                "The largest method-scene errors account for a disproportionate share of each arm's sum of squared error. The table identifies those cases explicitly so subsequent surface-prior and cost-curve diagnosis can target the errors that control RMSE rather than optimizing only the median scene."
            ),
        },
        {"id": "outlier-table-block", "type": "table", "tableId": "outlier-table"},
        {
            "id": "next-steps",
            "type": "markdown",
            "body": (
                "## Recommended next steps\n\n"
                "1. Run the leading properly masked arms over the full 250 before comparing them with the 64.4% campaign score or the approximately 73% low-cloud score.\n"
                "2. Keep **LUT + prior TCWV, masked** as the controlled diagnostic reference; it currently has the best error magnitude among properly masked arms.\n"
                "3. Do not promote the legacy mask bypass from its raw hit lead. First inspect its added cloudy/water cases and require spatial/physical quality checks.\n"
                "4. Treat MODIS time-series smoothing as neutral until it produces a larger paired gain; one net hit in this cohort is not enough.\n"
                "5. Focus the next surface-prior work on the high-AOD outliers listed here and on sparse-support abstentions, while preserving fully cloudy abstention.\n"
                "6. Validate any winner on an independent year/site split and benchmark runtime with explicit cold-cache and warm-cache protocols."
            ),
        },
        {
            "id": "further-questions",
            "type": "markdown",
            "body": (
                "## Further questions\n\n"
                "- Which high-AOD misses are caused by surface-prior bias versus RT/cost-curve limitations?\n"
                "- Can sparse valid support be pooled without admitting cloud- or water-contaminated pixels?\n"
                "- Does the anchored prior improve after controlling the tau-prior and predictor differences against the MODIS arms?\n"
                "- Do the paired rankings persist on the full 250-matchup campaign and on unseen sites?"
            ),
        },
    ]

    datasets = {
        "headline": headline,
        "cohort_context": _widget_rows(cohort_context),
        "coverage_accuracy": coverage_accuracy,
        "common_summary": common_chart,
        "error_summary": error_summary,
        "factorial_net": paired_chart,
        "regime_summary": regime_chart,
        "runtime_summary": runtime_summary,
        "method_summary": _widget_rows(method_summary),
        "paired_summary": _widget_rows(paired_summary),
        "method_configs": _widget_rows(method_configs),
        "outliers": _widget_rows(outliers[:16]),
    }
    return {
        "surface": "report",
        "manifest": {
            "version": 1,
            "surface": "report",
            "title": title,
            "description": "Technical performance analysis of eight completed AOD retrieval experiments on a fixed 66-scene diagnostic cohort.",
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


def build(root: Path, output_dir: Path) -> dict[str, Any]:
    mids_path = root / "masked_r2_regression_diagnostic_mids.txt"
    mids = [
        line.strip() for line in mids_path.read_text(encoding="utf-8").splitlines() if line.strip()
    ]
    if len(mids) != len(set(mids)):
        raise ValueError("The matchup manifest contains duplicate IDs.")
    records, audit = _load_records(root, mids)
    if not audit["ready"]:
        raise ValueError(f"Result audit failed: {audit['issues']}")

    full_mids = [
        line.strip()
        for line in (root / "campaign250_mids.txt").read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    full_old_records = _load_directory(root / FULL_OLD_DIR)
    full_masked_records = _load_directory(root / FULL_MASKED_DIR)
    for label, full_records in (
        ("historical R2", full_old_records),
        ("masked R2", full_masked_records),
    ):
        missing_full = sorted(set(full_mids) - set(full_records))
        if missing_full:
            raise ValueError(f"{label} is missing {len(missing_full)} campaign records.")
    cohort_context, selection_context, low_cloud_ids = _full_campaign_context(
        full_mids=full_mids,
        diagnostic_mids=mids,
        old_records=full_old_records,
        new_records=full_masked_records,
    )
    if not selection_context["rebuilt_selection_matches_manifest"]:
        raise ValueError(
            "The 66-scene manifest does not match the documented diagnostic selection."
        )

    masked_reference_records = records["lut_prior_masked"]
    masked_reference_common = [
        mid
        for mid in mids
        if _valid(full_old_records.get(mid)) and _valid(masked_reference_records.get(mid))
    ]
    selection_context["historical_vs_masked_reference"] = {
        "historical_hits": sum(_hit(full_old_records.get(mid)) for mid in mids),
        "masked_reference_hits": sum(_hit(masked_reference_records.get(mid)) for mid in mids),
        "shared_valid": len(masked_reference_common),
        "historical_shared_hits": sum(
            _hit(full_old_records.get(mid)) for mid in masked_reference_common
        ),
        "masked_reference_shared_hits": sum(
            _hit(masked_reference_records.get(mid)) for mid in masked_reference_common
        ),
        "historical_hits_lost_to_masked_abstention": sum(
            _hit(full_old_records.get(mid)) and not _valid(masked_reference_records.get(mid))
            for mid in mids
        ),
    }

    method_summary = [_summary(spec, records[spec.key], mids) for spec in METHODS]
    for spec, row in zip(METHODS, method_summary, strict=True):
        low_cloud = _summary(spec, records[spec.key], low_cloud_ids)
        row.update(
            {
                "low_cloud_expected": low_cloud["expected"],
                "low_cloud_valid": low_cloud["valid"],
                "low_cloud_hits": low_cloud["hits"],
                "low_cloud_coverage": low_cloud["coverage"],
                "low_cloud_strict_rate": low_cloud["strict_rate"],
                "low_cloud_valid_hit_rate": low_cloud["valid_hit_rate"],
            }
        )
    common_ids = [
        mid for mid in mids if all(_valid(records[spec.key].get(mid)) for spec in METHODS)
    ]
    common_summary = [_summary(spec, records[spec.key], common_ids) for spec in METHODS]

    reference = records[METHODS[0].key]
    regimes = ["clean", "mod", "high"]
    regime_summary: list[dict[str, Any]] = []
    for spec in METHODS:
        for regime in regimes:
            regime_ids = [mid for mid in mids if str(reference[mid].get("regime")) == regime]
            row = _summary(spec, records[spec.key], regime_ids)
            row["regime"] = regime
            regime_summary.append(row)

    by_key = {spec.key: spec for spec in METHODS}
    pair_requests = [
        ("Surface prior", "modis_single_day", "modis_timeseries_smooth", "MODIS smooth vs daily"),
        (
            "Surface prior",
            "modis_single_day",
            "swir_nir_anchored",
            "SWIR+NIR anchor vs MODIS daily",
        ),
        (
            "Surface prior",
            "modis_timeseries_smooth",
            "swir_nir_anchored",
            "SWIR+NIR anchor vs MODIS smooth",
        ),
        ("Factorial", "lut_prior_masked", "lut_prior_cloud", "Cloud allowed vs LUT masked"),
        ("Factorial", "lut_prior_masked", "lut_prior_legacy", "Legacy bypass vs LUT masked"),
        ("Factorial", "lut_prior_masked", "lut_l2awvp_masked", "L2A WVP vs prior TCWV"),
        ("Factorial", "lut_prior_masked", "sixs_prior_masked", "6S vs LUT RT"),
    ]
    paired_summary = [
        _paired_summary(group, by_key[baseline], by_key[candidate], records, mids, comparison)
        for group, baseline, candidate, comparison in pair_requests
    ]

    detail_rows: list[dict[str, Any]] = []
    for spec in METHODS:
        valid_records = [record for record in records[spec.key].values() if _valid(record)]
        squared_sum = sum(_error(record) ** 2 for record in valid_records)
        for record in valid_records:
            error = _error(record)
            detail_rows.append(
                {
                    "key": spec.key,
                    "method": spec.label,
                    "short_method": spec.short_label,
                    "matchup_id": str(record.get("matchup_id")),
                    "site": str(record.get("site", "")),
                    "regime": str(record.get("regime", "")),
                    "truth": _finite(record.get("truth")),
                    "retrieved": _finite(record.get("retrieved")),
                    "error": error,
                    "abs_error": abs(error),
                    "within_ee": _hit(record),
                    "sse_share": error * error / squared_sum if squared_sum else 0.0,
                    "cloud_frac": _finite(record.get("cloud_frac")),
                    "valid_aot_fraction": _finite(record.get("valid_aot_fraction")),
                }
            )
    outliers = sorted(detail_rows, key=lambda row: row["abs_error"], reverse=True)
    method_configs = [{"order": index + 1, **asdict(spec)} for index, spec in enumerate(METHODS)]

    generated_at = datetime.now(UTC).replace(microsecond=0).isoformat().replace("+00:00", "Z")
    analysis = {
        "generated_at": generated_at,
        "root": str(root),
        "manifest": mids_path.name,
        "expected_error_formula": EXPECTED_ERROR_FORMULA,
        "audit": audit,
        "selection_context": selection_context,
        "cohort_context": cohort_context,
        "common_valid_ids": common_ids,
        "common_valid_count": len(common_ids),
        "method_summary": method_summary,
        "common_summary": common_summary,
        "regime_summary": regime_summary,
        "paired_summary": paired_summary,
        "outliers": outliers,
        "method_configs": method_configs,
    }
    artifact = _artifact(
        generated_at=generated_at,
        audit=audit,
        method_summary=method_summary,
        common_summary=common_summary,
        regime_summary=regime_summary,
        paired_summary=paired_summary,
        outliers=outliers,
        method_configs=method_configs,
        common_count=len(common_ids),
        cohort_context=cohort_context,
        selection_context=selection_context,
        low_cloud_count=len(low_cloud_ids),
    )
    chart_map = [
        {
            "section": "Cohort context",
            "question": "Why are diagnostic rates lower than the full-campaign validation?",
            "family": "Comparison",
            "type": "grouped bar",
            "fields": ["cohort", "series", "strict_rate"],
            "claim": "Outcome-selected diagnostic scenes are much harder than the full 250, including after a cloud-below-20% filter.",
            "palette": "Hard two-root cap: blue and orange, plus neutrals.",
        },
        {
            "section": "Coverage and strict score",
            "question": "How much does coverage influence strict ranking?",
            "family": "Comparison",
            "type": "grouped bar",
            "fields": ["method", "strict_rate", "coverage"],
            "claim": "Mask bypass raises coverage and raw hits while increasing available-case error.",
            "palette": "Hard two-root cap: blue and orange, plus neutrals.",
        },
        {
            "section": "Common-case robustness",
            "question": "Which method leads on identical valid scenes?",
            "family": "Comparison and ranking",
            "type": "ranked horizontal bar",
            "fields": ["method", "within_ee_rate"],
            "claim": "Common-case ranking removes the coverage advantage.",
            "palette": "Single-root sequential blue.",
        },
        {
            "section": "Error magnitude",
            "question": "Are mean errors driven by outliers?",
            "family": "Comparison",
            "type": "grouped bar",
            "fields": ["method", "mae", "median_abs_error"],
            "claim": "Mean error materially exceeds median error for outlier-sensitive arms.",
            "palette": "Hard two-root cap: blue and orange, plus neutrals.",
        },
        {
            "section": "Factorial effects",
            "question": "Which controlled changes add strict hits?",
            "family": "Comparison",
            "type": "signed horizontal bar",
            "fields": ["candidate", "net_hits"],
            "claim": "Mask bypass adds raw hits; L2A WVP and 6S do not improve the masked reference.",
            "palette": "Diverging around zero without red/green semantics.",
        },
        {
            "section": "AOD regime",
            "question": "Where in the AOD range do methods fail?",
            "family": "Matrix",
            "type": "heatmap",
            "fields": ["method", "regime", "strict_rate"],
            "claim": "High-AOD performance remains weak across all methods.",
            "palette": "Single-root sequential blue.",
        },
        {
            "section": "Runtime",
            "question": "What wall time did the harness record?",
            "family": "Comparison",
            "type": "horizontal bar",
            "fields": ["method", "median_runtime_s"],
            "claim": "Observed runtime varies, but cache and host effects prevent algorithmic attribution.",
            "palette": "Single-root sequential blue.",
        },
    ]
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
                "cohort context and selection",
                "coverage",
                "common-case ranking",
                "error magnitude",
                "surface-prior comparison",
                "factorial effects",
                "regime performance",
                "runtime",
            ],
            "validation details": "Merged into paired methodology and validation.",
        },
        "chart_map": chart_map,
        "data_quality": audit,
        "selection_context": selection_context,
        "validation_assessment": "Ready to share with explicit cohort, mask-bypass, multiplicity, and runtime caveats.",
        "omissions": [
            "No production-population estimate: the 66-scene cohort is failure-enriched.",
            "No causal claim for SWIR+NIR versus MODIS: more than one configuration dimension changes.",
            "No algorithmic runtime claim: cache and cluster placement were not controlled.",
        ],
    }

    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "analysis.json").write_text(json.dumps(analysis, indent=2), encoding="utf-8")
    (output_dir / "artifact.json").write_text(json.dumps(artifact, indent=2), encoding="utf-8")
    (output_dir / "source-notes.json").write_text(
        json.dumps(source_notes, indent=2), encoding="utf-8"
    )
    report_rows = [
        {"dataset": dataset, **row}
        for dataset, dataset_rows in artifact["snapshot"]["datasets"].items()
        for row in dataset_rows
    ]
    report_fields = ["dataset"] + sorted(
        {key for row in report_rows for key in row if key != "dataset"}
    )
    _write_csv(output_dir / "report-data.csv", report_rows, report_fields)
    _write_csv(output_dir / "method-summary.csv", method_summary)
    _write_csv(output_dir / "common-summary.csv", common_summary)
    _write_csv(output_dir / "regime-summary.csv", regime_summary)
    _write_csv(output_dir / "paired-summary.csv", paired_summary)
    _write_csv(output_dir / "outliers.csv", outliers)
    return {
        "output_dir": str(output_dir),
        "artifact": str(output_dir / "artifact.json"),
        "analysis": str(output_dir / "analysis.json"),
        "source_notes": str(output_dir / "source-notes.json"),
        "audit_ready": audit["ready"],
        "common_valid": len(common_ids),
        "methods": len(METHODS),
        "records": sum(row["present"] for row in method_summary),
        "failed": sum(row["failed"] for row in method_summary),
        "missing": sum(row["missing"] for row in method_summary),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args()
    output_dir = args.output_dir or args.root / DEFAULT_REPORT_DIR
    receipt = build(args.root, output_dir)
    print(json.dumps(receipt, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
