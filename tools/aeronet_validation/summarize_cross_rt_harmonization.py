"""Summarize the cross-RT surface harmonisation experiment without plots.

Surface-pair metrics use only same-day L2A/L1C pairs and grouped out-of-fold
predictions. AERONET appears only in the final retrieval evaluation. The output
is a compact HTML, JSON, CSV and Markdown package intended for experiment review.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
from collections import Counter
from pathlib import Path
from typing import Any

ROOT = Path("/gws/ssde/j25a/nceo_isp/public/siac_refactor")
DEFAULT_MODEL_ROOT = ROOT / "analysis/cross_rt_harmonizer_lowcloud152_20260716"
DEFAULT_CASES = ROOT / "reports/aod-final-performance-dashboard-20260713/data/all-cases.csv"
DEFAULT_SPLIT_MANIFEST = ROOT / "analysis/l2a_harmonization_retrieval_score_20260714/summary.json"
DEFAULT_COHORT_MANIFEST = ROOT / "campaign250_lowcloud20_mids.json"
VISIBLE_BANDS = ("blue", "green", "red")
ALL_BANDS = ("coastal", "blue", "green", "red", "nir08", "swir16", "swir22")
MODE_BANDS = {
    "blue": ("blue",),
    "solver": ("blue", "green", "red"),
    "visible": ("coastal", "blue", "green", "red"),
    "all": ALL_BANDS,
}
DELTA_AOT_BIN_LABELS = (
    "< -0.10",
    "-0.10 to -0.02",
    "-0.02 to 0.02",
    "0.02 to 0.10",
    "> 0.10",
)
HISTORY_TAGS = (
    "identity_daily",
    "cross_rt_baseline_a1_solver_cap0p030",
    "cross_rt_aod_a1_solver_cap0p030",
    "cross_rt_atmosphere_a1_solver_cap0p030",
    "cross_rt_terrain_a1_solver_cap0p015",
    "cross_rt_terrain_a1_solver_cap0p030",
    "cross_rt_terrain_a1_all_cap0p015",
    "cross_rt_terrain_a1_all_cap0p030",
    "cross_rt_terrain_a10_solver_cap0p030",
)
TARGET_WITHIN_EE_RATE = 0.87
CONTROL_DIRS = {
    "previous_best": "phaseD_results_lowcloud20_geometry_backstop05_b03_chi2_20260712",
    "current_physical": (
        "phaseD_results_lowcloud20_physical_anchor_maiac_l2awvp_dem_scenegeom_20260715"
    ),
}


def _finite(value: Any) -> float | None:
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None


def _retrieval_metrics(rows: list[dict[str, Any]]) -> dict[str, Any]:
    valid = [
        row
        for row in rows
        if row.get("status") == "OK"
        and _finite(row.get("truth")) is not None
        and _finite(row.get("retrieved")) is not None
    ]
    errors = [float(row["retrieved"]) - float(row["truth"]) for row in valid]
    hits = sum(bool(row.get("within_ee")) for row in valid)
    return {
        "expected": len(rows),
        "n": len(valid),
        "within_ee": hits,
        # The benchmark is strict: missing and non-OK outcomes remain in the
        # frozen denominator and therefore count as misses.
        "within_ee_rate": hits / len(rows) if rows else None,
        "valid_within_ee_rate": hits / len(valid) if valid else None,
        "bias": sum(errors) / len(errors) if errors else None,
        "mae": sum(abs(error) for error in errors) / len(errors) if errors else None,
        "rmse": (
            math.sqrt(sum(error * error for error in errors) / len(errors)) if errors else None
        ),
        "cloud_mask_bypass_count": sum(
            bool((row.get("solver") or {}).get("surface_cloud_mask_bypassed")) for row in valid
        ),
        "water_mask_bypass_count": sum(
            bool((row.get("solver") or {}).get("surface_water_mask_bypassed")) for row in valid
        ),
        "status_counts": dict(Counter(str(row.get("status", "MISSING")) for row in rows)),
    }


def _regime(truth: float) -> str:
    if truth < 0.25:
        return "low_lt_0p25"
    if truth <= 0.85:
        return "medium_0p25_0p85"
    return "high_gt_0p85"


def _load_results(
    path: Path,
    matchup_ids: list[str],
    truth_by_matchup_id: dict[str, float] | None = None,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for matchup_id in matchup_ids:
        result_path = path / f"{matchup_id}.json"
        if not result_path.exists():
            rows.append({"matchup_id": matchup_id, "status": "MISSING"})
            continue
        try:
            rows.append(json.loads(result_path.read_text(encoding="utf-8")))
        except (OSError, json.JSONDecodeError) as exc:
            rows.append(
                {
                    "matchup_id": matchup_id,
                    "status": "MALFORMED",
                    "reason": str(exc),
                }
            )
    if truth_by_matchup_id is not None:
        for row in rows:
            matchup_id = str(row["matchup_id"])
            if matchup_id not in truth_by_matchup_id:
                continue
            expected_truth = truth_by_matchup_id[matchup_id]
            result_truth = _finite(row.get("truth"))
            if result_truth is None:
                row["truth"] = expected_truth
            elif not math.isclose(result_truth, expected_truth, rel_tol=1.0e-9, abs_tol=1.0e-9):
                raise ValueError(f"retrieval truth differs from frozen case table: {matchup_id}")
    return rows


def _evaluation_split(path: Path, matchup_ids: list[str]) -> tuple[dict[str, str], dict[str, Any]]:
    """Load the locked development/holdout assignment used by earlier harmonisation tests."""
    payload = json.loads(path.read_text(encoding="utf-8"))
    cohort = payload["cohort"]
    fold_by_matchup_id = {key: int(value) for key, value in cohort["fold_by_matchup_id"].items()}
    holdout_folds = {int(value) for value in cohort["holdout_folds"]}
    missing = sorted(set(matchup_ids) - set(fold_by_matchup_id))
    if missing:
        raise ValueError(f"locked split is missing {len(missing)} matchup(s)")
    assignment = {
        matchup_id: (
            "holdout" if fold_by_matchup_id[matchup_id] in holdout_folds else "development"
        )
        for matchup_id in matchup_ids
    }
    return assignment, {
        "manifest": str(path),
        "seed": cohort.get("split_seed"),
        "holdout_folds": sorted(holdout_folds),
        "development_n": sum(value == "development" for value in assignment.values()),
        "holdout_n": sum(value == "holdout" for value in assignment.values()),
    }


def _cohort_definition(
    path: Path,
    matchup_ids: list[str],
    case_rows: list[dict[str, str]],
) -> dict[str, Any]:
    """Verify the frozen low-cloud IDs and retain the controlling metric definition."""
    manifest = json.loads(path.read_text(encoding="utf-8"))
    mids_path = Path(str(manifest["selected_mids_file"]))
    serialized = mids_path.read_text(encoding="utf-8")
    selected_ids = [line.strip() for line in serialized.splitlines() if line.strip()]
    if selected_ids != matchup_ids:
        raise ValueError("case CSV does not match the ordered frozen low-cloud cohort")
    digest = hashlib.sha256(serialized.encode("utf-8")).hexdigest()
    if digest != str(manifest.get("selected_mids_sha256")):
        raise ValueError("frozen low-cloud cohort checksum does not match its manifest")
    metadata_cloud = [_finite(row.get("cloud_cover")) for row in case_rows]
    return {
        "manifest": str(path),
        "selected_mids_file": str(mids_path),
        "selected_mids_sha256": digest,
        "cloud_source_directory": manifest.get("cloud_source_directory"),
        "cloud_field": manifest.get("cloud_field"),
        "comparison": manifest.get("comparison"),
        "threshold_fraction": manifest.get("threshold_fraction"),
        "campaign_count": manifest.get("campaign_count"),
        "selected_count": len(selected_ids),
        "metadata_scene_cloud_cover_below_20_count": sum(
            value is not None and value < 20.0 for value in metadata_cloud
        ),
        "metadata_scene_cloud_cover_note": (
            "Context field only; it did not define the frozen low-cloud cohort."
        ),
    }


def _pair_archive_health(pair_dir: Path | None, matchup_ids: list[str]) -> dict[str, Any]:
    """Summarize terminal cases separately from rejected tile-scene components."""
    audits: list[dict[str, Any]] = []
    missing: list[str] = []
    for matchup_id in matchup_ids:
        if pair_dir is None:
            missing.append(matchup_id)
            continue
        path = pair_dir / f"{matchup_id}.json"
        if not path.exists():
            missing.append(matchup_id)
            continue
        try:
            audits.append(json.loads(path.read_text(encoding="utf-8")))
        except (OSError, json.JSONDecodeError):
            missing.append(matchup_id)
    rejected = [error for audit in audits for error in audit.get("errors", [])]
    sparse_rejections = sum(
        str(error.get("reason", "")).startswith("only ")
        and str(error.get("reason", "")).endswith(" paired clear-land pixels")
        for error in rejected
    )
    status_counts = Counter(str(audit.get("status", "MISSING")) for audit in audits)
    attempted_components = sum(int(audit.get("attempted_scenes", 0)) for audit in audits)
    retained_components = sum(int(audit.get("successful_scenes", 0)) for audit in audits)
    return {
        "pair_archive": str(pair_dir) if pair_dir is not None else None,
        "audits": len(audits),
        "missing_or_malformed": len(missing),
        "missing_matchup_ids": missing,
        "status_counts": dict(status_counts),
        "terminal_cases": int(status_counts["OK"] + status_counts["DATA_UNAVAILABLE"]),
        "mapped_cases": int(status_counts["OK"]),
        "identity_fallback_cases": int(status_counts["DATA_UNAVAILABLE"]),
        "identity_fallback_matchup_ids": sorted(
            str(audit.get("matchup_id"))
            for audit in audits
            if audit.get("status") == "DATA_UNAVAILABLE"
        ),
        "attempted_tile_scene_components": attempted_components,
        "retained_tile_scene_components": retained_components,
        "rejected_tile_scene_components": len(rejected),
        # Kept as compatibility aliases for earlier machine-readable reports.
        "attempted_acquisitions": attempted_components,
        "retained_acquisitions": retained_components,
        "rejected_acquisitions": len(rejected),
        "sparse_clear_land_rejections": sparse_rejections,
        "other_rejections": len(rejected) - sparse_rejections,
        "rejection_types": dict(
            Counter(str(error.get("error_type", "unspecified")) for error in rejected)
        ),
    }


def _history_archive_health(audit_dir: Path, matchup_ids: list[str]) -> dict[str, Any]:
    """Audit per-case history outcomes and every expected candidate artifact."""
    audits: list[dict[str, Any]] = []
    missing_or_malformed: list[str] = []
    for matchup_id in matchup_ids:
        path = audit_dir / f"{matchup_id}.json"
        if not path.exists():
            missing_or_malformed.append(matchup_id)
            continue
        try:
            audits.append(json.loads(path.read_text(encoding="utf-8")))
        except (OSError, json.JSONDecodeError):
            missing_or_malformed.append(matchup_id)

    status_counts = Counter(str(audit.get("status", "MISSING")) for audit in audits)
    fallback = [
        audit
        for audit in audits
        if audit.get("status") == "OK" and audit.get("mapping_applied") is False
    ]
    expected_outputs = [
        Path(path)
        for audit in audits
        if audit.get("status") == "OK"
        for path in (audit.get("outputs") or {}).values()
    ]
    missing_outputs = [str(path) for path in expected_outputs if not path.exists()]
    return {
        "audit_directory": str(audit_dir),
        "audits": len(audits),
        "missing_or_malformed": len(missing_or_malformed),
        "missing_matchup_ids": missing_or_malformed,
        "status_counts": dict(status_counts),
        "terminal_cases": int(status_counts["OK"] + status_counts["DATA_UNAVAILABLE"]),
        "per_acquisition_harmonized_cases": sum(
            audit.get("status") == "OK"
            and audit.get("application") == "per acquisition before tile mosaic and temporal median"
            for audit in audits
        ),
        "uncorrected_single_source_fallback_cases": len(fallback),
        "uncorrected_single_source_fallback_matchup_ids": sorted(
            str(audit.get("matchup_id")) for audit in fallback
        ),
        "fallback_reason_counts": dict(
            Counter(str(audit.get("skip_reason", "unspecified")) for audit in fallback)
        ),
        "low_temporal_support_cases": sum(
            bool(audit.get("low_temporal_support")) for audit in audits
        ),
        "expected_candidate_outputs": len(expected_outputs),
        "present_candidate_outputs": len(expected_outputs) - len(missing_outputs),
        "missing_candidate_outputs": missing_outputs,
        "nonfatal_scene_errors": sum(len(audit.get("errors") or []) for audit in audits),
    }


def _surface_summary(metrics: dict[str, Any]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    references = {
        "identity": metrics["identity"],
        "current_aod_offset": metrics["current_aod_offset"],
    }
    for name, band_metrics in references.items():
        rows.append(_aggregate_surface_candidate(name, "uncapped", band_metrics))
    for name, caps in metrics["candidates"].items():
        for cap_name, band_metrics in caps.items():
            rows.append(_aggregate_surface_candidate(name, cap_name, band_metrics))
    return rows


def _aggregate_surface_candidate(
    model: str,
    cap: str,
    band_metrics: dict[str, dict[str, float]],
) -> dict[str, Any]:
    visible = [band_metrics[band] for band in VISIBLE_BANDS]
    return {
        "model": model,
        "cap": cap,
        "visible_scene_mae": sum(row["scene_mae"] for row in visible) / len(visible),
        "visible_scene_rmse": sum(row["scene_rmse"] for row in visible) / len(visible),
        "visible_abs_scene_bias": sum(abs(row["scene_bias"]) for row in visible) / len(visible),
        "visible_pixel_mae": sum(row["mae"] for row in visible) / len(visible),
    }


def _paired_transitions(
    baseline: list[dict[str, Any]],
    candidate: list[dict[str, Any]],
) -> dict[str, int]:
    before = {row.get("matchup_id"): row for row in baseline if row.get("status") == "OK"}
    after = {row.get("matchup_id"): row for row in candidate if row.get("status") == "OK"}
    counts: Counter[str] = Counter()
    for matchup_id in sorted(set(before) | set(after)):
        if matchup_id not in before or matchup_id not in after:
            counts["unpaired"] += 1
            continue
        old_hit = bool(before[matchup_id].get("within_ee"))
        new_hit = bool(after[matchup_id].get("within_ee"))
        if old_hit and new_hit:
            counts["stable_hit"] += 1
        elif old_hit:
            counts["loss"] += 1
        elif new_hit:
            counts["gain"] += 1
        else:
            counts["stable_miss"] += 1
    return dict(counts)


def _parse_variant_tag(tag: str) -> tuple[str, str, str, float] | None:
    """Split an experiment history tag into (model, mode, cap_key, cap_value)."""
    match = re.fullmatch(
        r"(?P<model>.+)_(?P<mode>blue|solver|visible|all)_cap(?P<cap>\d+p\d+)", tag
    )
    if match is None:
        return None
    cap_value = float(match.group("cap").replace("p", "."))
    return match.group("model"), match.group("mode"), f"cap_{cap_value:.3f}", cap_value


def _percentile(sorted_values: list[float], fraction: float) -> float:
    """Linear-interpolation percentile of an ascending list (numpy default rule)."""
    if not sorted_values:
        return float("nan")
    position = fraction * (len(sorted_values) - 1)
    lower = math.floor(position)
    upper = math.ceil(position)
    weight = position - lower
    return sorted_values[lower] * (1.0 - weight) + sorted_values[upper] * weight


def _delta_aot_bin(delta_aot: float) -> int:
    if delta_aot < -0.10:
        return 0
    if delta_aot < -0.02:
        return 1
    if delta_aot <= 0.02:
        return 2
    if delta_aot <= 0.10:
        return 3
    return 4


def _scene_correction_stats(
    scene_metrics_path: Path,
    model: str,
    cap_key: str,
    bands: tuple[str, ...],
) -> dict[str, Any]:
    """Summarize the applied per-component correction and its MAIAC-minus-Sen2Cor AOD response.

    The scene-mean applied correction per band equals the corrected-minus-identity
    component bias, because both errors share the same target reflectance.
    """
    prefix = f"{model}_{cap_key}".replace(".", "p")
    corrections: dict[str, list[float]] = {band: [] for band in bands}
    bins: list[dict[str, Any]] = [
        {
            "label": label,
            "delta_aot": [],
            "corrections": {band: [] for band in bands},
            "identity_blue_bias": [],
            "corrected_blue_bias": [],
        }
        for label in DELTA_AOT_BIN_LABELS
    ]
    with scene_metrics_path.open() as stream:
        for row in csv.DictReader(stream):
            delta_aot = _finite(row.get("delta_aot_maiac_minus_sen2cor"))
            component_bin = bins[_delta_aot_bin(delta_aot)] if delta_aot is not None else None
            if component_bin is not None:
                component_bin["delta_aot"].append(delta_aot)
            for band in bands:
                identity_bias = _finite(row.get(f"identity_{band}_bias"))
                corrected_bias = _finite(row.get(f"{prefix}_{band}_bias"))
                if identity_bias is None or corrected_bias is None:
                    continue
                correction = corrected_bias - identity_bias
                corrections[band].append(correction)
                if component_bin is not None:
                    component_bin["corrections"][band].append(correction)
                    if band == "blue":
                        component_bin["identity_blue_bias"].append(identity_bias)
                        component_bin["corrected_blue_bias"].append(corrected_bias)

    def _mean(values: list[float]) -> float | None:
        return sum(values) / len(values) if values else None

    per_band: dict[str, Any] = {}
    for band, values in corrections.items():
        ordered = sorted(values)
        per_band[band] = {
            "n": len(ordered),
            "mean": _mean(ordered),
            "p05": _percentile(ordered, 0.05),
            "median": _percentile(ordered, 0.50),
            "p95": _percentile(ordered, 0.95),
            "median_abs": _percentile(sorted(abs(value) for value in ordered), 0.50),
        }
    bin_rows = [
        {
            "label": component_bin["label"],
            "n": len(component_bin["delta_aot"]),
            "median_delta_aot": _percentile(sorted(component_bin["delta_aot"]), 0.50),
            "mean_correction": {
                band: _mean(values) for band, values in component_bin["corrections"].items()
            },
            "identity_blue_bias_mean": _mean(component_bin["identity_blue_bias"]),
            "corrected_blue_bias_mean": _mean(component_bin["corrected_blue_bias"]),
        }
        for component_bin in bins
    ]
    return {
        "unit": "tile-scene component mean (pixel corrections averaged per component)",
        "per_band": per_band,
        "delta_aot_bins": bin_rows,
    }


def _coefficient_table(
    artifact: dict[str, Any],
    model_name: str,
    bands: tuple[str, ...],
) -> dict[str, Any]:
    """Export the applied standardized ridge coefficients, band-aligned by feature."""
    model = artifact["models"][model_name]
    protocol = artifact.get("crossfit_protocol") or {}
    holdout_fold = str(protocol.get("holdout_model_fold", ""))
    fold_models = model.get("folds") or {}
    if holdout_fold and holdout_fold in fold_models:
        band_models = fold_models[holdout_fold]
        source = "development-trained holdout-fold model (per-fold coefficients differ slightly)"
    else:
        band_models = model["full"]
        source = "full-data model"
    aligned: dict[str, dict[str, float]] = {}
    intercepts: dict[str, float] = {}
    for band in bands:
        band_model = band_models[band]
        intercepts[band] = float(band_model["intercept"])
        for name, coefficient in zip(band_model["feature_names"], band_model["coef"], strict=True):
            generic = name.replace(f"l2a_{band}", "l2a_band")
            aligned.setdefault(generic, {})[band] = float(coefficient)
    rows = sorted(
        (
            {"feature": feature, "coefficients": coefficients}
            for feature, coefficients in aligned.items()
        ),
        key=lambda row: -max(abs(value) for value in row["coefficients"].values()),
    )
    return {"source": source, "intercepts": intercepts, "rows": rows}


def _correction_model_report(
    artifact: dict[str, Any],
    surface_metrics: dict[str, Any],
    scene_metrics_path: Path,
    variant_tag: str,
    *,
    selection_basis: str,
) -> dict[str, Any] | None:
    """Describe how the selected variant's correction is derived and how accurate it is."""
    parsed = _parse_variant_tag(variant_tag)
    if parsed is None:
        return None
    model, mode, cap_key, cap_value = parsed
    if model not in (artifact.get("models") or {}) or model not in surface_metrics["candidates"]:
        return None
    bands = MODE_BANDS[mode]
    model_payload = artifact["models"][model]
    oof_scale = model_payload["oof_residual_scale"]
    band_rows: list[dict[str, Any]] = []
    for band in bands:
        identity = surface_metrics["identity"][band]
        corrected = surface_metrics["candidates"][model][cap_key][band]
        band_rows.append(
            {
                "band": band,
                "identity_scene_bias": identity["scene_bias"],
                "corrected_scene_bias": corrected["scene_bias"],
                "identity_scene_mae": identity["scene_mae"],
                "corrected_scene_mae": corrected["scene_mae"],
                "scene_mae_change_pct": 100.0
                * (corrected["scene_mae"] / identity["scene_mae"] - 1.0),
                "identity_pixel_p95_abs": identity["p95_abs"],
                "corrected_pixel_p95_abs": corrected["p95_abs"],
                "oof_residual_median": oof_scale[band]["median"],
                "oof_residual_sigma_mad": oof_scale[band]["mad_to_sigma"],
                "oof_residual_rmse": oof_scale[band]["rmse"],
            }
        )
    sample_band = model_payload["full"][bands[0]]
    return {
        "variant": variant_tag,
        "selection_basis": selection_basis,
        "model": model,
        "mode": mode,
        "cap_key": cap_key,
        "correction_cap": cap_value,
        "corrected_bands": list(bands),
        "regression": {
            "family": "per-band ridge on standardized features",
            "alpha": model_payload.get("alpha"),
            "weighting_group": model_payload.get("weighting_group"),
            "feature_count": len(sample_band["feature_names"]),
            "target_residual": "target reflectance minus operational Sen2Cor L2A reflectance",
        },
        "target": artifact.get("target"),
        "effect": {"bands": band_rows},
        "scene_corrections": _scene_correction_stats(scene_metrics_path, model, cap_key, bands),
        "coefficients": _coefficient_table(artifact, model, bands),
    }


def _rate_text(value: Any) -> str:
    parsed = _finite(value)
    return f"{100.0 * parsed:.1f}%" if parsed is not None else "n/a"


def _non_ok_status_text(status_counts: dict[str, int]) -> str:
    non_ok = [
        f"{status}:{count}" for status, count in sorted(status_counts.items()) if status != "OK"
    ]
    return ", ".join(non_ok) if non_ok else "none"


def _best_complete_candidate(summary: dict[str, Any]) -> dict[str, Any] | None:
    """Rank variants only after the entire experimental matrix is complete."""
    candidates = [row for row in summary["retrieval"] if row["variant"] in HISTORY_TAGS]
    if len(candidates) != len(HISTORY_TAGS) or any(
        row["metrics"]["n"] != summary["cohort_size"]
        or row["metrics"]["within_ee_rate"] is None
        or row["development_metrics"]["n"] != summary["evaluation_split"]["development_n"]
        or row["development_metrics"]["within_ee_rate"] is None
        for row in candidates
    ):
        return None
    return max(
        candidates,
        key=lambda row: row["development_metrics"]["within_ee_rate"],
        default=None,
    )


def _result_commentary(summary: dict[str, Any]) -> list[str]:
    """Produce factual comments while refusing to rank incomplete experiments."""
    surface_best = min(summary["surface"], key=lambda row: row["visible_scene_mae"])
    surface_identity = next(row for row in summary["surface"] if row["model"] == "identity")
    surface_reduction = 100.0 * (
        1.0 - surface_best["visible_scene_mae"] / surface_identity["visible_scene_mae"]
    )
    comments = [
        f"Lowest grouped surface-pair error: {surface_best['model']} "
        f"({surface_best['cap']}) at {surface_best['visible_scene_mae']:.5f} visible "
        f"acquisition MAE, {surface_reduction:.1f}% below the unmodified L2A history."
    ]
    complete_candidate = _best_complete_candidate(summary)
    if complete_candidate is None:
        incomplete = [
            row
            for row in summary["retrieval"]
            if row["variant"] in HISTORY_TAGS and row["metrics"]["n"] != summary["cohort_size"]
        ]
        comments.extend(
            [
                "No experimental variant has all cohort results, so the retrieval "
                "variants must not yet be ranked.",
                f"Incomplete experimental variants: {len(incomplete)} of {len(HISTORY_TAGS)}.",
            ]
        )
        return comments

    identity = next(row for row in summary["retrieval"] if row["variant"] == "identity_daily")
    previous_best = next(
        (row for row in summary["retrieval"] if row["variant"] == "previous_best"),
        None,
    )
    best_rate = float(complete_candidate["metrics"]["within_ee_rate"])
    development_rate = float(complete_candidate["development_metrics"]["within_ee_rate"])
    holdout_rate = float(complete_candidate["holdout_metrics"]["within_ee_rate"])
    identity_rate = _finite(identity["development_metrics"]["within_ee_rate"])
    medium = complete_candidate["regimes"]["medium_0p25_0p85"]
    transition = complete_candidate.get("transitions_vs_identity", {})
    target_gap = 100.0 * (best_rate - TARGET_WITHIN_EE_RATE)
    medium_bias = _finite(medium["bias"])
    medium_mae = _finite(medium["mae"])
    comments.extend(
        [
            f"Development-selected variant: {complete_candidate['variant']} at "
            f"{_rate_text(development_rate)} within EE on the locked development split.",
            f"Its locked holdout score is {_rate_text(holdout_rate)}; the complete-cohort "
            f"score is {_rate_text(best_rate)} ({complete_candidate['metrics']['within_ee']}/"
            f"{complete_candidate['metrics']['n']}).",
            f"Its medium-AOD (0.25-0.85) score is {_rate_text(medium['within_ee_rate'])} "
            f"with bias {medium_bias:.4f} and MAE {medium_mae:.4f}."
            if medium_bias is not None and medium_mae is not None
            else "Its medium-AOD (0.25-0.85) regime has no valid score.",
            f"Paired against the unmodified history: {transition.get('gain', 0)} gains and "
            f"{transition.get('loss', 0)} losses.",
            f"The result is {abs(target_gap):.1f} percentage points "
            f"{'above' if target_gap >= 0.0 else 'below'} the 87% target.",
        ]
    )
    previous_rate = (
        _finite(previous_best["metrics"]["within_ee_rate"]) if previous_best is not None else None
    )
    if (
        previous_best is not None
        and previous_rate is not None
        and previous_best["metrics"]["n"] == summary["cohort_size"]
    ):
        comments.insert(
            3,
            f"Change from the previous-best control: "
            f"{100.0 * (best_rate - previous_rate):+.1f} percentage points "
            f"({complete_candidate['metrics']['within_ee']} versus "
            f"{previous_best['metrics']['within_ee']} hits).",
        )
    mapped_metrics = complete_candidate.get("mapped_history_metrics")
    fallback_metrics = complete_candidate.get("fallback_history_metrics")
    if mapped_metrics and fallback_metrics:
        comments.insert(
            4,
            f"Mapped-history cases score {_rate_text(mapped_metrics['within_ee_rate'])} "
            f"({mapped_metrics['within_ee']}/{mapped_metrics['expected']}); unchanged "
            f"same-source fallbacks score {_rate_text(fallback_metrics['within_ee_rate'])} "
            f"({fallback_metrics['within_ee']}/{fallback_metrics['expected']}).",
        )
    if (
        identity_rate is not None
        and identity["development_metrics"]["n"] == summary["evaluation_split"]["development_n"]
    ):
        comments.insert(
            3,
            f"Development change from the identity-history control: "
            f"{100.0 * (development_rate - identity_rate):+.1f} percentage points.",
        )
    return comments


def _reflectance_text(value: Any) -> str:
    parsed = _finite(value)
    return f"{parsed:+.5f}" if parsed is not None else "n/a"


def _correction_markdown(correction: dict[str, Any]) -> list[str]:
    """Render derivation, applied effect, and standalone accuracy of the correction."""
    bands = correction["corrected_bands"]
    regression = correction["regression"]
    scene_corrections = correction["scene_corrections"]
    lines = [
        "",
        "## Correction model: derivation",
        "",
        f"- Reported variant: {correction['variant']} ({correction['selection_basis']}).",
        f"- Target frame: {correction['target']} — the mapping moves operational Sen2Cor "
        "L2A reflectance onto the MAIAC-AOD/current-RT surface.",
        f"- Form: corrected = clip(L2A + clip(r_band(x), ±{correction['correction_cap']:.3f}), "
        "0.001, 0.8). Each r_band is a ridge regression "
        f"(alpha={regression['alpha']}, {regression['feature_count']} standardized features, "
        f"equal weight per {regression['weighting_group']}) fit to per-pixel residuals "
        f"({regression['target_residual']}) on the same-day exact pairs.",
        f"- Corrected bands: {', '.join(bands)}; other bands pass through unchanged.",
        "- Application uses the site's cross-fit fold model — locked-holdout sites receive "
        "the model trained only on development sites — with each historical acquisition's "
        "own same-day MAIAC and Sen2Cor AOD/WVP, scene geometry, sensor, processing "
        "baseline, month, and per-pixel GLO-30 terrain.",
        f"- Coefficients below are from the {correction['coefficients']['source']}, in "
        "reflectance units per +1 standard deviation of each feature.",
        "",
        "| Feature | " + " | ".join(bands) + " |",
        "|---|" + "---:|" * len(bands),
        "| (intercept) | "
        + " | ".join(
            _reflectance_text(correction["coefficients"]["intercepts"].get(band)) for band in bands
        )
        + " |",
    ]
    for row in correction["coefficients"]["rows"]:
        lines.append(
            f"| {row['feature']} | "
            + " | ".join(_reflectance_text(row["coefficients"].get(band)) for band in bands)
            + " |"
        )
    lines.extend(
        [
            "",
            "## Correction effect on the L2A history",
            "",
            "Out-of-fold, grouped by acquisition. The applied correction is summarized per "
            f"{scene_corrections['unit']}.",
            "",
            "| Band | Bias before | Bias after | MAE before | MAE after | MAE change | "
            "Correction p05/median/p95 | Median abs correction |",
            "|---|---:|---:|---:|---:|---:|---:|---:|",
        ]
    )
    for row in correction["effect"]["bands"]:
        applied = scene_corrections["per_band"][row["band"]]
        lines.append(
            f"| {row['band']} | {_reflectance_text(row['identity_scene_bias'])} | "
            f"{_reflectance_text(row['corrected_scene_bias'])} | "
            f"{row['identity_scene_mae']:.5f} | {row['corrected_scene_mae']:.5f} | "
            f"{row['scene_mae_change_pct']:+.1f}% | "
            f"{_reflectance_text(applied['p05'])} / {_reflectance_text(applied['median'])} / "
            f"{_reflectance_text(applied['p95'])} | {applied['median_abs']:.5f} |"
        )
    lines.extend(
        [
            "",
            "Response to the MAIAC-minus-Sen2Cor AOD difference (mean applied correction per "
            "retained tile-scene component):",
            "",
            "| Delta-AOT bin | Components | Median delta | "
            + " | ".join(f"Mean {band}" for band in bands)
            + " | Blue bias before -> after |",
            "|---|---:|---:|" + "---:|" * len(bands) + "---:|",
        ]
    )
    for bin_row in scene_corrections["delta_aot_bins"]:
        before = _reflectance_text(bin_row["identity_blue_bias_mean"])
        after = _reflectance_text(bin_row["corrected_blue_bias_mean"])
        lines.append(
            f"| {bin_row['label']} | {bin_row['n']} | "
            f"{_reflectance_text(bin_row['median_delta_aot'])} | "
            + " | ".join(_reflectance_text(bin_row["mean_correction"].get(band)) for band in bands)
            + f" | {before} -> {after} |"
        )
    lines.extend(
        [
            "",
            "## Correction accuracy",
            "",
            "Residuals are corrected-minus-target reflectance on grouped out-of-fold pairs; "
            "the spread columns come from uncapped out-of-fold predictions.",
            "",
            "| Band | Acquisition MAE after (before) | Residual median | Robust sigma (MAD) | "
            "Residual RMSE | Pixel p95 abs error before -> after |",
            "|---|---:|---:|---:|---:|---:|",
        ]
    )
    for row in correction["effect"]["bands"]:
        lines.append(
            f"| {row['band']} | {row['corrected_scene_mae']:.5f} "
            f"({row['identity_scene_mae']:.5f}) | "
            f"{_reflectance_text(row['oof_residual_median'])} | "
            f"{row['oof_residual_sigma_mad']:.5f} | {row['oof_residual_rmse']:.5f} | "
            f"{row['identity_pixel_p95_abs']:.5f} -> {row['corrected_pixel_p95_abs']:.5f} |"
        )
    return lines


def _markdown(summary: dict[str, Any]) -> str:
    def number(value: Any) -> str:
        parsed = _finite(value)
        return f"{parsed:.4f}" if parsed is not None else "n/a"

    result_comments = _result_commentary(summary)
    lines = [
        "# Cross-RT L2A harmonisation experiment",
        "",
        "The surface model uses no AERONET values. AERONET is used only for the final retrieval score.",
        "",
        "## Technical summary",
        "",
        *[f"- {comment}" for comment in result_comments],
        "",
        "## Scope and method",
        "",
        f"- Cohort: {summary['cohort_size']} Sentinel-2 matchups frozen by masked-R2 "
        f"`{summary['cohort']['cloud_field']} < "
        f"{summary['cohort']['threshold_fraction']:.2f}`. The separate scene-metadata "
        f"`cloud_cover` field did not select the cohort; "
        f"{summary['cohort']['metadata_scene_cloud_cover_below_20_count']}/"
        f"{summary['cohort_size']} rows are below 20% in that contextual field.",
        f"- Pair support: {summary['experiment']['sample_count']} sampled pixels from "
        f"{summary['experiment']['acquisition_count']} independent acquisitions "
        f"({summary['experiment']['scene_count']} tile-scene components) and "
        f"{summary['experiment']['site_count']} sites.",
        f"- Pair archive health before temporal filtering: "
        f"{summary['experiment']['pair_health']['terminal_cases']}/"
        f"{summary['cohort_size']} terminal case archives; "
        f"{summary['experiment']['pair_health']['mapped_cases']} mapped and "
        f"{summary['experiment']['pair_health']['identity_fallback_cases']} uncorrected "
        "single-source fallbacks; "
        f"{summary['experiment']['pair_health']['retained_tile_scene_components']}/"
        f"{summary['experiment']['pair_health']['attempted_tile_scene_components']} "
        "tile-scene components "
        "retained. Rejected components: "
        f"{summary['experiment']['pair_health']['sparse_clear_land_rejections']} sparse clear-land, "
        f"{summary['experiment']['pair_health']['other_rejections']} other.",
        f"- History health: {summary['experiment']['history_health']['terminal_cases']}/"
        f"{summary['cohort_size']} terminal cases and "
        f"{summary['experiment']['history_health']['present_candidate_outputs']}/"
        f"{summary['experiment']['history_health']['expected_candidate_outputs']} candidate "
        "files present; "
        f"{summary['experiment']['history_health']['per_acquisition_harmonized_cases']} cases "
        "use per-acquisition mapping and "
        f"{summary['experiment']['history_health']['uncorrected_single_source_fallback_cases']} "
        "retain the unchanged operational L2A prior. Fallback reasons: "
        f"{summary['experiment']['history_health']['fallback_reason_counts']}.",
        f"- AOD support by retained tile-scene component: MAIAC-minus-Sen2Cor delta "
        "p05/median/p95 = "
        f"{summary['experiment']['scene_distribution']['delta_aot_maiac_minus_sen2cor']['p05']:.3f}/"
        f"{summary['experiment']['scene_distribution']['delta_aot_maiac_minus_sen2cor']['median']:.3f}/"
        f"{summary['experiment']['scene_distribution']['delta_aot_maiac_minus_sen2cor']['p95']:.3f}; "
        f"{summary['experiment']['scene_distribution']['near_equal_aod_abs_delta_le_0p02']} "
        "near-equal-AOD components (absolute delta <= 0.02).",
        f"- Training cutoff: {summary['experiment']['training_cutoff']}; target-date imagery "
        "is excluded from history construction.",
        "- Validation: five-fold grouped cross-fitting holds out entire sites. Each independent "
        "overpass has equal regression weight, including overpasses split across adjacent tiles.",
        f"- Surface holdout isolation: development folds "
        f"{summary['experiment']['crossfit_protocol'].get('development_folds')} cross-fit only "
        f"against other development folds; holdout folds "
        f"{summary['experiment']['crossfit_protocol'].get('holdout_folds')} use one surface "
        "model trained on all development sites and no holdout sites.",
        f"- Retrieval selection: the locked development split has "
        f"{summary['evaluation_split']['development_n']} cases; the locked holdout has "
        f"{summary['evaluation_split']['holdout_n']} cases. Candidate ranking uses development "
        "only.",
        "- Retrieval scoring: within-EE rates use the frozen expected denominator; missing and "
        "non-OK outcomes count as misses.",
        "- Target: same-day L1C TOA corrected with MAIAC AOD, operational L2A per-pixel WVP, "
        "GLO-30 elevation, scene-mean angles, and the current libRadtran LUT.",
        "- Baseline stage learns the non-zero L2A/current-RT difference that can remain when "
        "MAIAC and Sen2Cor AOD are equal.",
        "- AOD, atmosphere, and terrain stages add delta-AOD interactions, water-vapour and "
        "geometry terms, then local elevation/slope/solar-incidence terms.",
        "- Band-scope ablation: solver variants change B02/B03/B04 only; all-band variants "
        "also place the historical B8A/B11/B12 anchors in the learned current-RT frame.",
        "- Application: correct each historical L2A acquisition before tile mosaicking and "
        "temporal compositing; retrieval physics and solver settings are fixed across variants.",
        "- Unsupported cases retain the same operational L2A prior; the experiment never "
        "switches those cases to another surface source.",
        "",
        "## Surface-pair validation",
        "",
        "| Model | Cap | Visible acquisition MAE | Visible acquisition RMSE | "
        "Absolute acquisition bias |",
        "|---|---:|---:|---:|---:|",
    ]
    for row in sorted(summary["surface"], key=lambda value: value["visible_scene_mae"]):
        lines.append(
            f"| {row['model']} | {row['cap']} | {row['visible_scene_mae']:.5f} | "
            f"{row['visible_scene_rmse']:.5f} | {row['visible_abs_scene_bias']:.5f} |"
        )
    if summary.get("correction_model"):
        lines.extend(_correction_markdown(summary["correction_model"]))
    lines.extend(
        [
            "",
            "## Retrieval validation",
            "",
            "| Variant | Valid | Non-OK | Mask bypass | Within EE | Bias | MAE | RMSE | Gains | Losses |",
            "|---|---:|---|---:|---:|---:|---:|---:|---:|---:|",
        ]
    )
    for row in summary["retrieval"]:
        metric = row["metrics"]
        transition = row.get("transitions_vs_identity", {})
        lines.append(
            f"| {row['variant']} | {metric['n']}/{summary['cohort_size']} | "
            f"{_non_ok_status_text(metric['status_counts'])} | "
            f"{metric['cloud_mask_bypass_count']} | "
            f"{_rate_text(metric.get('within_ee_rate'))} | "
            f"{number(metric.get('bias'))} | {number(metric.get('mae'))} | "
            f"{number(metric.get('rmse'))} | {transition.get('gain', 0)} | "
            f"{transition.get('loss', 0)} |"
        )
    lines.extend(
        [
            "",
            "## Generalisation check",
            "",
            "| Variant | Development | Holdout | Holdout minus development |",
            "|---|---:|---:|---:|",
        ]
    )
    for row in summary["retrieval"]:
        development_rate = _finite(row["development_metrics"]["within_ee_rate"])
        holdout_rate = _finite(row["holdout_metrics"]["within_ee_rate"])
        delta = (
            f"{100.0 * (holdout_rate - development_rate):+.1f} pp"
            if development_rate is not None and holdout_rate is not None
            else "n/a"
        )
        lines.append(
            f"| {row['variant']} | {_rate_text(development_rate)} "
            f"({row['development_metrics']['n']}) | {_rate_text(holdout_rate)} "
            f"({row['holdout_metrics']['n']}) | {delta} |"
        )
    lines.extend(
        [
            "",
            "## AOD regimes",
            "",
            "| Variant | Low <0.25 | Medium 0.25-0.85 | Medium bias | Medium MAE | High >0.85 |",
            "|---|---:|---:|---:|---:|---:|",
        ]
    )
    for row in summary["retrieval"]:
        regimes = row["regimes"]
        low = regimes["low_lt_0p25"]
        medium = regimes["medium_0p25_0p85"]
        high = regimes["high_gt_0p85"]
        lines.append(
            f"| {row['variant']} | {_rate_text(low['within_ee_rate'])} ({low['n']}) | "
            f"{_rate_text(medium['within_ee_rate'])} ({medium['n']}) | "
            f"{number(medium['bias'])} | {number(medium['mae'])} | "
            f"{_rate_text(high['within_ee_rate'])} ({high['n']}) |"
        )
    lines.extend(
        [
            "",
            "## Interpretation guardrails",
            "",
            "- Compare surface stages using grouped out-of-fold acquisition metrics, not pixel "
            "count.",
            "- Compare retrieval variants only when all 152 terminal results are present.",
            "- Select candidates on the locked development split; use holdout only to assess "
            "generalisation.",
            "- The holdout is isolated from fitting and ranking in this experiment, but the "
            "campaign has been inspected previously; it is not a fresh external confirmation.",
            "- AOD performance is an external evaluation and did not select pair-training targets.",
            "- The report identifies the development-selected complete variant and reports its "
            "locked holdout result; it does not promote it to production automatically.",
            "",
        ]
    )
    return "\n".join(lines)


def _html(markdown_text: str) -> str:
    """Render the compact report as a dependency-free static webpage."""
    import markdown

    body = markdown.markdown(markdown_text, extensions=("tables", "sane_lists"))
    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <meta name="color-scheme" content="light">
  <title>Cross-RT L2A harmonisation experiment</title>
  <style>
    :root {{
      color: #20252b;
      background: #f5f6f7;
      font-family: Inter, ui-sans-serif, system-ui, -apple-system, BlinkMacSystemFont,
        "Segoe UI", sans-serif;
      font-synthesis: none;
      letter-spacing: 0;
    }}
    * {{ box-sizing: border-box; }}
    body {{ margin: 0; line-height: 1.55; }}
    header {{
      background: #20252b;
      color: #fff;
      border-bottom: 4px solid #36a269;
    }}
    .header-inner, main, footer {{ width: min(1180px, calc(100% - 32px)); margin: 0 auto; }}
    .header-inner {{ padding: 20px 0 18px; }}
    .eyebrow {{ margin: 0 0 4px; color: #b9c2cb; font-size: 0.78rem; font-weight: 700; }}
    header h1 {{ margin: 0; font-size: 1.45rem; line-height: 1.25; }}
    nav {{ display: flex; flex-wrap: wrap; gap: 16px; margin-top: 12px; }}
    nav a {{ color: #dbe9ff; font-size: 0.86rem; font-weight: 650; }}
    main {{ padding: 28px 0 44px; }}
    main > h1 {{ display: none; }}
    h2 {{
      margin: 34px 0 12px;
      padding-top: 10px;
      border-top: 1px solid #c8cdd2;
      color: #17202a;
      font-size: 1.15rem;
    }}
    p, li {{ max-width: 96ch; }}
    ul {{ margin: 8px 0 18px; padding-left: 22px; }}
    li + li {{ margin-top: 6px; }}
    code {{
      padding: 1px 4px;
      border: 1px solid #d4d8dc;
      border-radius: 3px;
      background: #fff;
      font-size: 0.9em;
    }}
    table {{
      display: block;
      width: 100%;
      margin: 12px 0 24px;
      overflow-x: auto;
      border-collapse: collapse;
      background: #fff;
      border: 1px solid #c8cdd2;
    }}
    th, td {{
      min-width: 98px;
      padding: 9px 11px;
      border-bottom: 1px solid #dfe2e5;
      text-align: right;
      vertical-align: top;
      white-space: nowrap;
      font-variant-numeric: tabular-nums;
    }}
    th {{ background: #e9edf0; color: #17202a; font-size: 0.78rem; }}
    th:first-child, td:first-child {{ min-width: 255px; text-align: left; }}
    tbody tr:nth-child(even) {{ background: #f8f9fa; }}
    footer {{
      padding: 18px 0 32px;
      border-top: 1px solid #c8cdd2;
      color: #57616b;
      font-size: 0.82rem;
    }}
    @media (max-width: 620px) {{
      .header-inner, main, footer {{ width: min(100% - 20px, 1180px); }}
      .header-inner {{ padding-top: 16px; }}
      main {{ padding-top: 18px; }}
      th, td {{ padding: 8px; }}
    }}
  </style>
</head>
<body>
  <header>
    <div class="header-inner">
      <p class="eyebrow">Sentinel-2 low-cloud validation | 152 frozen matchups</p>
      <h1>Cross-RT L2A harmonisation experiment</h1>
      <nav aria-label="Report files">
        <a href="visual.html">Visual report</a>
        <a href="summary.json">Summary JSON</a>
        <a href="retrieval_metrics.csv">Retrieval CSV</a>
        <a href="surface_ablation.csv">Surface CSV</a>
        <a href="correction_effect.csv">Correction CSV</a>
        <a href="summary.md">Markdown</a>
      </nav>
    </div>
  </header>
  <main>{body}</main>
  <footer>Static experiment report. No AERONET observations were used to train the surface mapping.</footer>
</body>
</html>
"""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model-root", type=Path, default=DEFAULT_MODEL_ROOT)
    parser.add_argument("--cases", type=Path, default=DEFAULT_CASES)
    parser.add_argument("--cohort-manifest", type=Path, default=DEFAULT_COHORT_MANIFEST)
    parser.add_argument("--split-manifest", type=Path, default=DEFAULT_SPLIT_MANIFEST)
    parser.add_argument("--root", type=Path, default=ROOT)
    parser.add_argument("--output", type=Path, default=None)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    output = args.output or args.model_root / "summary"
    output.mkdir(parents=True, exist_ok=True)
    case_rows = list(csv.DictReader(args.cases.open()))
    matchup_ids = [row["matchup_id"] for row in case_rows]
    truth_by_matchup_id = {
        row["matchup_id"]: float(row["truth"])
        for row in case_rows
        if _finite(row.get("truth")) is not None
    }
    cohort_definition = _cohort_definition(args.cohort_manifest, matchup_ids, case_rows)
    evaluation_split, split_provenance = _evaluation_split(args.split_manifest, matchup_ids)
    surface_metrics = json.loads(
        (args.model_root / "surface_metrics.json").read_text(encoding="utf-8")
    )
    artifact = json.loads((args.model_root / "harmonizer.json").read_text(encoding="utf-8"))
    artifact_split = artifact.get("evaluation_split")
    if artifact_split is not None and (
        artifact_split.get("split_seed") != split_provenance.get("seed")
        or artifact_split.get("holdout_folds") != split_provenance.get("holdout_folds")
    ):
        raise ValueError("harmonizer and retrieval report use different locked splits")
    pair_archive = artifact.get("pair_archive") or surface_metrics.get("pair_archive")
    pair_health = _pair_archive_health(Path(pair_archive) if pair_archive else None, matchup_ids)
    history_health = _history_archive_health(
        args.model_root / "daily_histories/audits", matchup_ids
    )
    scene_distribution = surface_metrics.get("scene_distribution") or {
        "delta_aot_maiac_minus_sen2cor": {
            "n": 0,
            "p05": float("nan"),
            "median": float("nan"),
            "p95": float("nan"),
        },
        "near_equal_aod_abs_delta_le_0p02": 0,
    }
    surface = _surface_summary(surface_metrics)

    retrieval_rows: list[dict[str, Any]] = []
    loaded: dict[str, list[dict[str, Any]]] = {}
    fallback_matchup_ids = set(history_health["uncorrected_single_source_fallback_matchup_ids"])
    variant_dirs = {
        **CONTROL_DIRS,
        **{
            tag: f"phaseD_results_lowcloud20_crossrt_{tag}_physical_20260716"
            for tag in HISTORY_TAGS
        },
    }
    for variant, directory in variant_dirs.items():
        rows = _load_results(args.root / directory, matchup_ids, truth_by_matchup_id)
        for row in rows:
            row["evaluation_split"] = evaluation_split[str(row["matchup_id"])]
        loaded[variant] = rows
        by_regime: dict[str, dict[str, Any]] = {}
        for regime in ("low_lt_0p25", "medium_0p25_0p85", "high_gt_0p85"):
            by_regime[regime] = _retrieval_metrics(
                [
                    row
                    for row in rows
                    if _finite(row.get("truth")) is not None
                    and _regime(float(row["truth"])) == regime
                ]
            )
        retrieval_rows.append(
            {
                "variant": variant,
                "directory": directory,
                "metrics": _retrieval_metrics(rows),
                "development_metrics": _retrieval_metrics(
                    [row for row in rows if row["evaluation_split"] == "development"]
                ),
                "holdout_metrics": _retrieval_metrics(
                    [row for row in rows if row["evaluation_split"] == "holdout"]
                ),
                "mapped_history_metrics": _retrieval_metrics(
                    [row for row in rows if row["matchup_id"] not in fallback_matchup_ids]
                ),
                "fallback_history_metrics": _retrieval_metrics(
                    [row for row in rows if row["matchup_id"] in fallback_matchup_ids]
                ),
                "regimes": by_regime,
            }
        )
    identity = loaded["identity_daily"]
    for row in retrieval_rows:
        if row["variant"] not in {"identity_daily", "previous_best", "current_physical"}:
            row["transitions_vs_identity"] = _paired_transitions(identity, loaded[row["variant"]])

    summary = {
        "schema_version": 2,
        "uses_aeronet_for_surface_training": False,
        "cohort_size": len(matchup_ids),
        "cohort": cohort_definition,
        "evaluation_split": split_provenance,
        "target_within_ee_rate": TARGET_WITHIN_EE_RATE,
        "experiment": {
            "sample_count": surface_metrics["sample_count"],
            "scene_count": surface_metrics["scene_count"],
            "acquisition_count": surface_metrics.get(
                "acquisition_count", surface_metrics["scene_count"]
            ),
            "site_count": surface_metrics["site_count"],
            "training_cutoff": surface_metrics["training_cutoff"],
            "max_samples_per_scene": surface_metrics.get("max_samples_per_scene"),
            "model_family": surface_metrics.get("model_family"),
            "scene_distribution": scene_distribution,
            "target": artifact.get("target"),
            "target_rt": artifact.get("target_rt"),
            "fold_by_matchup_id": artifact.get("fold_by_matchup_id"),
            "crossfit_protocol": artifact.get("crossfit_protocol", {}),
            "pair_health": pair_health,
            "history_health": history_health,
        },
        "surface": surface,
        "retrieval": retrieval_rows,
    }
    complete_candidate = _best_complete_candidate(summary)
    if complete_candidate is not None:
        selected_tag = complete_candidate["variant"]
        selection_basis = "development-selected complete variant"
    else:
        selected_tag = next(
            (tag for tag in HISTORY_TAGS if _parse_variant_tag(tag) is not None), None
        )
        selection_basis = "first experimental variant; the matrix is incomplete, do not rank"
    summary["correction_model"] = (
        _correction_model_report(
            artifact,
            surface_metrics,
            args.model_root / "surface_scene_metrics.csv",
            selected_tag,
            selection_basis=selection_basis,
        )
        if selected_tag is not None
        else None
    )
    markdown_text = _markdown(summary)
    (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    (output / "summary.md").write_text(markdown_text, encoding="utf-8")
    (output / "index.html").write_text(_html(markdown_text), encoding="utf-8")
    with (output / "surface_ablation.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(surface[0]))
        writer.writeheader()
        writer.writerows(surface)
    if summary["correction_model"] is not None:
        correction = summary["correction_model"]
        effect_rows = [
            {
                **row,
                **{
                    f"applied_correction_{key}": value
                    for key, value in correction["scene_corrections"]["per_band"][
                        row["band"]
                    ].items()
                },
            }
            for row in correction["effect"]["bands"]
        ]
        with (output / "correction_effect.csv").open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(effect_rows[0]))
            writer.writeheader()
            writer.writerows(effect_rows)
    with (output / "retrieval_metrics.csv").open("w", newline="") as stream:
        fields = (
            "variant",
            "expected",
            "n",
            "within_ee",
            "within_ee_rate",
            "valid_within_ee_rate",
            "cloud_mask_bypass_count",
            "water_mask_bypass_count",
            "bias",
            "mae",
            "rmse",
            "development_expected",
            "development_n",
            "development_within_ee_rate",
            "holdout_expected",
            "holdout_n",
            "holdout_within_ee_rate",
        )
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for row in retrieval_rows:
            values = {
                "variant": row["variant"],
                **row["metrics"],
                "development_expected": row["development_metrics"]["expected"],
                "development_n": row["development_metrics"]["n"],
                "development_within_ee_rate": row["development_metrics"]["within_ee_rate"],
                "holdout_expected": row["holdout_metrics"]["expected"],
                "holdout_n": row["holdout_metrics"]["n"],
                "holdout_within_ee_rate": row["holdout_metrics"]["within_ee_rate"],
            }
            writer.writerow({field: values.get(field) for field in fields})
    print(json.dumps({"output": str(output), "cohort_size": len(matchup_ids)}, indent=2))


if __name__ == "__main__":
    main()
