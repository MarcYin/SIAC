"""Build the L2A-to-current-RT harmonization investigation webpage."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import shutil
import subprocess
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np

if TYPE_CHECKING:
    from collections.abc import Iterable

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from tools.aeronet_validation.analyze_medium_aod_physics import (  # noqa: E402
    compute_metrics,
    within_ee,
)
from tools.aeronet_validation.build_low_cloud_failure_explorer import (  # noqa: E402
    BAND_COLORS,
    BANDS,
    WAVELENGTHS,
    _median_curve,
    _normalise_curve,
    _pooled_maps,
    _save_spatial_figure,
    _window_mask,
)

ROOT = Path("/gws/ssde/j25a/nceo_isp/public/siac_refactor")
BASE_REPORT = ROOT / "reports/aod-medium-physics-20260713"
DEFAULT_OUTPUT = ROOT / "reports/l2a-harmonization-20260714"
PAIR_DIR = ROOT / "analysis/l2a_l1c_exact_pairs_mediumdev_20260713"
ALL_MODEL_DIR = ROOT / "analysis/l2a_l1c_harmonizer_mediumdev_20260713"
CLEAN_MODEL_DIR = ROOT / "analysis/l2a_l1c_harmonizer_clean015_mediumdev_20260714"
TERRAIN_PAIR_DIR = ROOT / "analysis/l2a_l1c_exact_pairs_terrain_mediumdev_20260714"
TERRAIN_ALL_MODEL_DIR = ROOT / "analysis/l2a_l1c_harmonizer_terrain_mediumdev_20260714"
TERRAIN_CLEAN_MODEL_DIR = ROOT / (
    "analysis/l2a_l1c_harmonizer_terrain_clean015_mediumdev_20260714"
)
BASELINE_DIR = ROOT / "phaseD_results_lowcloud20_mediumphysics_baseline_costcube_mediumdev_20260713"
BASELINE_CUBES = ROOT / "phaseD_cost_cubes_lowcloud20_mediumphysics_baseline_mediumdev_20260713"
CLEAN_RESULTS = ROOT / (
    "phaseD_results_lowcloud20_l2aharm_"
    "clean015_full_a100_all_cap0p030_direct_static_costcube_mediumdev_20260713"
)
CLEAN_CUBES = ROOT / (
    "phaseD_cost_cubes_lowcloud20_l2aharm_"
    "clean015_full_a100_all_cap0p030_direct_static_costcube_mediumdev_20260713"
)
TERRAIN_RESULTS_015 = ROOT / (
    "phaseD_results_lowcloud20_l2aharm_"
    "terrain_a100_all_cap0p015_direct_static_mediumdev_20260713"
)
TERRAIN_RESULTS_030 = ROOT / (
    "phaseD_results_lowcloud20_l2aharm_"
    "terrain_a100_all_cap0p030_direct_static_mediumdev_20260713"
)
TERRAIN_CUBES_015 = ROOT / (
    "phaseD_cost_cubes_lowcloud20_l2aharm_"
    "terrain_a100_all_cap0p015_direct_static_mediumdev_20260713"
)
TERRAIN_CUBES_030 = ROOT / (
    "phaseD_cost_cubes_lowcloud20_l2aharm_"
    "terrain_a100_all_cap0p030_direct_static_mediumdev_20260713"
)
TERRAIN_BLUE_RESULTS_015 = ROOT / (
    "phaseD_results_lowcloud20_l2aharm_"
    "terrain_a100_blue_cap0p015_direct_static_mediumdev_20260713"
)
TERRAIN_BLUE_RESULTS_030 = ROOT / (
    "phaseD_results_lowcloud20_l2aharm_"
    "terrain_a100_blue_cap0p030_direct_static_mediumdev_20260713"
)
TERRAIN_CLEAN_RESULTS_030 = ROOT / (
    "phaseD_results_lowcloud20_l2aharm_"
    "terrain_clean_a100_all_cap0p030_direct_static_mediumdev_20260713"
)
TERRAIN_CLEAN_CUBES_030 = ROOT / (
    "phaseD_cost_cubes_lowcloud20_l2aharm_"
    "terrain_clean_a100_all_cap0p030_direct_static_mediumdev_20260713"
)
TERRAIN_CLEAN_BLUE_RESULTS_030 = ROOT / (
    "phaseD_results_lowcloud20_l2aharm_"
    "terrain_clean_a100_blue_cap0p030_direct_static_mediumdev_20260713"
)
TERRAIN_CLEAN_BLUE_CUBES_030 = ROOT / (
    "phaseD_cost_cubes_lowcloud20_l2aharm_"
    "terrain_clean_a100_blue_cap0p030_direct_static_mediumdev_20260713"
)
TERRAIN_ARTIFACTS: dict[str, dict[str, Any]] = {
    "terrain_a100_cap015": {
        "results": TERRAIN_RESULTS_015,
        "cubes": TERRAIN_CUBES_015,
        "history_audits": TERRAIN_ALL_MODEL_DIR / "daily_histories/audits",
        "history_tag": "terrain_a100_all_cap0p015",
    },
    "terrain_a100_cap030": {
        "results": TERRAIN_RESULTS_030,
        "cubes": TERRAIN_CUBES_030,
        "history_audits": TERRAIN_ALL_MODEL_DIR / "daily_histories/audits",
        "history_tag": "terrain_a100_all_cap0p030",
    },
    "terrain_blue_cap015": {
        "results": TERRAIN_BLUE_RESULTS_015,
        "cubes": ROOT
        / "phaseD_cost_cubes_lowcloud20_l2aharm_"
        "terrain_a100_blue_cap0p015_direct_static_mediumdev_20260713",
        "history_audits": TERRAIN_ALL_MODEL_DIR / "daily_histories_blue/audits",
        "history_tag": "terrain_a100_blue_cap0p015",
    },
    "terrain_blue_cap030": {
        "results": TERRAIN_BLUE_RESULTS_030,
        "cubes": ROOT
        / "phaseD_cost_cubes_lowcloud20_l2aharm_"
        "terrain_a100_blue_cap0p030_direct_static_mediumdev_20260713",
        "history_audits": TERRAIN_ALL_MODEL_DIR / "daily_histories_blue/audits",
        "history_tag": "terrain_a100_blue_cap0p030",
    },
    "terrain_clean_a100_cap030": {
        "results": TERRAIN_CLEAN_RESULTS_030,
        "cubes": TERRAIN_CLEAN_CUBES_030,
        "history_audits": TERRAIN_CLEAN_MODEL_DIR / "daily_histories/audits",
        "history_tag": "terrain_a100_all_cap0p030",
    },
    "terrain_clean_blue_cap030": {
        "results": TERRAIN_CLEAN_BLUE_RESULTS_030,
        "cubes": TERRAIN_CLEAN_BLUE_CUBES_030,
        "history_audits": TERRAIN_CLEAN_MODEL_DIR / "daily_histories/audits",
        "history_tag": "terrain_a100_blue_cap0p030",
    },
}
WEB_ASSETS = Path(__file__).with_name("medium_aod_physics_report")
TARGET_HITS = 32
TARGET_RATE = 0.87
BAND_NAMES = ("coastal", "blue", "green", "red", "nir08", "swir16", "swir22")
DISPLAY_BANDS = ("B01", "B02", "B03", "B04", "B8A", "B11", "B12")

VARIANT_SPECS: tuple[dict[str, Any], ...] = (
    {
        "id": "fresh",
        "label": "Fresh current-code replay",
        "path": BASELINE_DIR,
        "family": "Reference",
        "description": "Current SIAC surface prior and solver with full cost-cube capture.",
    },
    {
        "id": "current_best",
        "label": "Current prior-conflict rule",
        "path": ROOT
        / "analysis/medium_aod_current_end_to_end_prior_conflict_z2576_development_20260713",
        "family": "Atmospheric prior",
        "description": "Best complete generic current-code candidate; no source routing.",
    },
    {
        "id": "identity",
        "label": "Operational L2A, no mapping",
        "path": ROOT
        / "phaseD_results_lowcloud20_l2aharm_identity_daily_off_mediumdev_20260713",
        "family": "Surface harmonization",
        "description": "Daily operational L2A histories without the legacy offset or learned mapping.",
    },
    {
        "id": "legacy_offset",
        "label": "Operational L2A + legacy offset",
        "path": ROOT
        / "phaseD_results_lowcloud20_l2aharm_identity_daily_on_mediumdev_20260713",
        "family": "Surface harmonization",
        "description": "Same histories with the existing three fixed AOD-dependent visible offsets.",
    },
    {
        "id": "all_aod_live",
        "label": "All-AOD pixel mapping, full path",
        "path": ROOT
        / "phaseD_results_lowcloud20_l2aharm_full_a100_all_cap0p030_off_mediumdev_20260713",
        "family": "Surface harmonization",
        "description": "Held-site all-AOD ridge mapping applied before the standard seasonal predictor.",
    },
    {
        "id": "all_aod",
        "label": "All-AOD pixel mapping, direct",
        "path": ROOT
        / "phaseD_results_lowcloud20_l2aharm_full_a100_all_cap0p030_direct_static_mediumdev_20260713",
        "family": "Surface harmonization",
        "description": "Same mapped histories through the direct prebuilt-history retrieval path.",
    },
    {
        "id": "all_aod_mapunc",
        "label": "All-AOD mapping + OOF floor",
        "path": ROOT
        / "phaseD_results_lowcloud20_l2aharm_full_a100_all_cap0p030_direct_mapunc_mediumdev_20260713",
        "family": "Surface uncertainty",
        "description": "All-AOD mapping with a visible uncertainty floor from held-site surface error.",
    },
    {
        "id": "all_aod_product_unc",
        "label": "All-AOD mapping + product uncertainty",
        "path": ROOT
        / "phaseD_results_lowcloud20_l2aharm_full_a100_all_cap0p030_direct_product_unc_mediumdev_20260713",
        "family": "Surface uncertainty",
        "description": "All-AOD mapping with the operational product uncertainty treatment.",
    },
    {
        "id": "all_aod_tau_gated",
        "label": "All-AOD mapping + gated tau surface",
        "path": ROOT
        / "phaseD_results_lowcloud20_l2aharm_full_a100_all_cap0p030_direct_gated8_mediumdev_20260713",
        "family": "Surface iteration",
        "description": "AOD-dependent surface evaluation enabled only above the fixed gate.",
    },
    {
        "id": "all_aod_tau_always",
        "label": "All-AOD mapping + tau surface",
        "path": ROOT
        / "phaseD_results_lowcloud20_l2aharm_full_a100_all_cap0p030_direct_always_mediumdev_20260713",
        "family": "Surface iteration",
        "description": "AOD-dependent surface evaluation at every candidate AOD node.",
    },
    {
        "id": "terrain_a100_cap015",
        "label": "Terrain mapping, cap 0.015",
        "path": TERRAIN_RESULTS_015,
        "family": "Terrain-conditioned surface mapping",
        "description": "One L2A source with GLO-30 elevation, local slope and acquisition-specific solar incidence; held-site ridge map with a conservative cap.",
    },
    {
        "id": "terrain_a100_cap030",
        "label": "Terrain mapping, cap 0.030",
        "path": TERRAIN_RESULTS_030,
        "family": "Terrain-conditioned surface mapping",
        "description": "Same single-source terrain-conditioned held-site map with the established 0.03 reflectance cap.",
    },
    {
        "id": "terrain_blue_cap015",
        "label": "Terrain B02-only mapping, cap 0.015",
        "path": TERRAIN_BLUE_RESULTS_015,
        "family": "Terrain-conditioned band ablation",
        "description": "The terrain-conditioned ridge correction is applied only to B02; B01/B03/B04 and the NIR/SWIR anchors remain operational L2A, with a conservative cap.",
    },
    {
        "id": "terrain_blue_cap030",
        "label": "Terrain B02-only mapping, cap 0.030",
        "path": TERRAIN_BLUE_RESULTS_030,
        "family": "Terrain-conditioned band ablation",
        "description": "The same fixed B02-only terrain mapping with the established 0.03 reflectance cap.",
    },
    {
        "id": "terrain_clean_a100_cap030",
        "label": "Clean terrain mapping, cap 0.030",
        "path": TERRAIN_CLEAN_RESULTS_030,
        "family": "Clean-day terrain-conditioned surface mapping",
        "description": "Terrain-conditioned ridge map trained only on clean MAIAC AOD <=0.15 pairs, then applied to the one operational L2A prior source.",
    },
    {
        "id": "terrain_clean_blue_cap030",
        "label": "Clean terrain B02-only mapping, cap 0.030",
        "path": TERRAIN_CLEAN_BLUE_RESULTS_030,
        "family": "Clean-day terrain-conditioned band ablation",
        "description": "The clean-domain terrain correction is applied only to B02; all other prior bands remain operational L2A.",
    },
    {
        "id": "paired_exact",
        "label": "Exact paired current-RT history",
        "path": ROOT
        / "phaseD_results_lowcloud20_l2aharm_paired_siac_daily_direct_static_mediumdev_20260713",
        "family": "Diagnostic upper bound",
        "description": "Historical L1C corrected directly with saved MAIAC/current-RT coefficients.",
    },
    {
        "id": "scene_mapping",
        "label": "Scene ExtraTrees mapping",
        "path": ROOT
        / "phaseD_results_lowcloud20_l2aharm_scene_et_all_cap0p030_direct_static_mediumdev_20260713",
        "family": "Surface harmonization",
        "description": "Held-site nonlinear scene-level residual model applied to all seven bands.",
    },
    {
        "id": "scene_mapping_mapunc",
        "label": "Scene mapping + uncertainty map",
        "path": ROOT
        / "phaseD_results_lowcloud20_l2aharm_scene_et_all_cap0p030_direct_mapunc_mediumdev_20260713",
        "family": "Surface uncertainty",
        "description": "Scene ExtraTrees mapping with its held-site prediction uncertainty.",
    },
    {
        "id": "scene_solver",
        "label": "Scene mapping, solver bands only",
        "path": ROOT
        / "phaseD_results_lowcloud20_l2aharm_scene_et_solver_cap0p030_direct_static_mediumdev_20260713",
        "family": "Surface harmonization",
        "description": "Scene mapping restricted to B02/B03/B04.",
    },
    {
        "id": "harmonized",
        "label": "Clean-day held-site harmonization",
        "path": CLEAN_RESULTS,
        "family": "Surface harmonization",
        "description": "MAIAC AOD <=0.15 training, full ridge features, alpha 100, all bands, +/-0.03 cap.",
    },
    {
        "id": "clean_visible",
        "label": "Clean-day mapping, B01-B04 only",
        "path": ROOT
        / "phaseD_results_lowcloud20_l2aharm_clean015_full_a100_visible_cap0p030_direct_static_mediumdev_20260713",
        "family": "Band ablation",
        "description": "Preserves raw NIR/SWIR anchors while mapping coastal and visible bands.",
    },
    {
        "id": "clean_solver",
        "label": "Clean-day mapping, B02-B04 only",
        "path": ROOT
        / "phaseD_results_lowcloud20_l2aharm_clean015_full_a100_solver_cap0p030_direct_static_mediumdev_20260713",
        "family": "Band ablation",
        "description": "Maps only the three bands used directly in the AOD surface cost.",
    },
    {
        "id": "clean_domain",
        "label": "Clean-day mapping, in-domain scenes only",
        "path": ROOT
        / "phaseD_results_lowcloud20_l2aharm_clean015_domain015_full_a100_all_cap0p030_direct_static_mediumdev_20260713",
        "family": "Domain guard",
        "description": "Applies the all-band map only when scene MAIAC AOD is within the <=0.15 training domain.",
    },
)

IMPLEMENTATION_FILES = (
    REPO_ROOT / "tools/aeronet_validation/build_l2a_l1c_harmonization_pairs.py",
    REPO_ROOT / "tools/aeronet_validation/train_l2a_l1c_harmonizer.py",
    REPO_ROOT / "tools/aeronet_validation/build_harmonized_l2a_histories.py",
    REPO_ROOT / "tools/aeronet_validation/terrain_features.py",
    REPO_ROOT / "tools/aeronet_validation/derive_harmonized_history_modes.py",
    REPO_ROOT / "tools/aeronet_validation/run_harmonized_retrieval.py",
    REPO_ROOT / "tools/aeronet_validation/l2a_l1c_harmonization_pairs_submit.sbatch",
    REPO_ROOT / "tools/aeronet_validation/clean015_harmonized_histories_submit.sbatch",
    REPO_ROOT / "tools/aeronet_validation/harmonized_l2a_retrieval_submit.sbatch",
    REPO_ROOT / "tools/aeronet_validation/l2a_l1c_terrain_pairs_submit.sbatch",
    REPO_ROOT / "tools/aeronet_validation/train_l2a_l1c_terrain_harmonizer.sbatch",
    REPO_ROOT / "tools/aeronet_validation/terrain_harmonized_histories_submit.sbatch",
    REPO_ROOT / "tools/aeronet_validation/terrain_blue_harmonized_histories_submit.sbatch",
    REPO_ROOT / "tools/aeronet_validation/terrain_harmonized_retrieval_submit.sbatch",
    REPO_ROOT / "tools/aeronet_validation/terrain_blue_harmonized_retrieval_submit.sbatch",
    REPO_ROOT / "tools/aeronet_validation/terrain_clean_harmonized_histories_submit.sbatch",
    REPO_ROOT / "tools/aeronet_validation/terrain_clean_harmonized_retrieval_submit.sbatch",
    REPO_ROOT / "tools/aeronet_validation/build_l2a_terrain_report.sbatch",
    REPO_ROOT / "python/siac/algorithms/surface/seasonal_predictor.py",
    REPO_ROOT / "python/siac/algorithms/solver/surface_driven.py",
    REPO_ROOT / "python/siac/algorithms/rt/lut/backend.py",
    REPO_ROOT / "src/siac_rs/src/surface_driven.rs",
)

JOB_IDS = (
    "39222612",
    "39222616",
    "39222717",
    "39222620",
    "39223049",
    "39223073",
    "39254023",
    "39254024",
    "39254071",
    "39254072",
    "39286948",
    "39286961",
    "39287410",
    "39287411",
    "39304054",
    "39301441",
    "39301442",
    "39306008",
    "39306010",
)
JOB_REPAIRS: dict[str, tuple[str, ...]] = {
    # Slurm retains the submission-time script.  The first terrain retrieval
    # array therefore retained a pre-fix nested executable invocation.
    "39287411": ("39304054",),
}


def _finite(value: Any) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def _json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _load_result(path: Path) -> dict[str, Any] | None:
    try:
        result = _json(path)
    except (OSError, json.JSONDecodeError):
        return None
    return result


def _load_records(
    directory: Path, matchup_ids: Iterable[str]
) -> tuple[dict[str, dict[str, Any]], Counter[str]]:
    records: dict[str, dict[str, Any]] = {}
    statuses: Counter[str] = Counter()
    for matchup_id in matchup_ids:
        record = _load_result(directory / f"{matchup_id}.json")
        if record is None:
            statuses["MISSING"] += 1
            continue
        status = str(record.get("status", "UNKNOWN")).upper()
        statuses[status] += 1
        if status == "OK" and _finite(record.get("retrieved")) is not None:
            records[matchup_id] = record
    return records, statuses


def _metrics(records: dict[str, dict[str, Any]], matchup_ids: list[str]) -> dict[str, Any]:
    result = compute_metrics(records[mid] for mid in matchup_ids if mid in records).as_dict()
    result["within_ee_percent"] = (
        None
        if result["within_ee_rate"] is None
        else 100.0 * result["within_ee_rate"]
    )
    return result


def _transitions(
    reference: dict[str, dict[str, Any]],
    candidate: dict[str, dict[str, Any]],
    matchup_ids: list[str],
) -> dict[str, int]:
    counts: Counter[str] = Counter()
    for matchup_id in matchup_ids:
        if matchup_id not in reference or matchup_id not in candidate:
            counts["unpaired"] += 1
            continue
        truth = float(reference[matchup_id]["truth"])
        old = within_ee(float(reference[matchup_id]["retrieved"]), truth)
        new = within_ee(float(candidate[matchup_id]["retrieved"]), truth)
        key = (
            "gain"
            if new and not old
            else "loss"
            if old and not new
            else "stable_hit"
            if old
            else "stable_miss"
        )
        counts[key] += 1
    return dict(counts)


def _candidate(label: str, value: Any, truth: float, source: str) -> dict[str, Any]:
    number = _finite(value)
    ee = 0.05 + 0.15 * truth
    return {
        "label": label,
        "source": source,
        "value": number,
        "error": None if number is None else number - truth,
        "error_over_ee": None if number is None else (number - truth) / ee,
        "within_ee": None if number is None else within_ee(number, truth),
    }


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _surface_metric_rows() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for training_set, directory in (
        ("All retained AOD", ALL_MODEL_DIR),
        ("Clean MAIAC AOD <=0.15", CLEAN_MODEL_DIR),
    ):
        metrics = _json(directory / "surface_metrics.json")
        methods = (
            ("Raw operational L2A", metrics["identity"]),
            ("Legacy AOD offset", metrics["current_aod_offset"]),
            (
                "Held-site full ridge, alpha 100, cap 0.03",
                metrics["candidates"]["full_a100"]["cap_0.030"],
            ),
        )
        for method, values in methods:
            for band_name, display_name in zip(BAND_NAMES, DISPLAY_BANDS, strict=True):
                item = values[band_name]
                rows.append(
                    {
                        "training_set": training_set,
                        "method": method,
                        "band": display_name,
                        "scene_bias": item["scene_bias"],
                        "scene_mae": item["scene_mae"],
                        "scene_rmse": item["scene_rmse"],
                        "pixel_bias": item["bias"],
                        "pixel_mae": item["mae"],
                        "pixel_rmse": item["rmse"],
                        "pixel_p95_abs": item["p95_abs"],
                        "samples": item["n"],
                    }
                )
    return rows


def _scene_key(row: dict[str, str]) -> tuple[str, str, str, str]:
    return (row["matchup_id"], row["window"], row["day"], row["scene_id"])


def _surface_site_rows(clean_rows: list[dict[str, str]]) -> list[dict[str, Any]]:
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in clean_rows:
        grouped[row["site"]].append(row)
    output: list[dict[str, Any]] = []
    for site, rows in sorted(grouped.items()):
        weights = np.asarray([float(row["sample_count"]) for row in rows])
        identity = np.asarray(
            [
                np.mean(
                    [
                        float(row[f"identity_{band}_mae"])
                        for band in ("blue", "green", "red")
                    ]
                )
                for row in rows
            ]
        )
        harmonized = np.asarray(
            [
                np.mean(
                    [
                        float(row[f"full_a100_cap_0p030_{band}_mae"])
                        for band in ("blue", "green", "red")
                    ]
                )
                for row in rows
            ]
        )
        output.append(
            {
                "site": site,
                "scene_count": len(rows),
                "sample_count": int(weights.sum()),
                "maiac_aot_median": float(
                    np.median([float(row["maiac_aot"]) for row in rows])
                ),
                "delta_aot_median": float(
                    np.median(
                        [
                            float(row["delta_aot_maiac_minus_sen2cor"])
                            for row in rows
                        ]
                    )
                ),
                "identity_visible_mae": float(np.average(identity, weights=weights)),
                "harmonized_visible_mae": float(
                    np.average(harmonized, weights=weights)
                ),
            }
        )
    return output


def _pair_audit() -> dict[str, Any]:
    audits = [_json(path) for path in sorted(PAIR_DIR.glob("*.json"))]
    errors: Counter[str] = Counter()
    for audit in audits:
        for error in audit.get("errors", []):
            errors[str(error.get("error_type") or error.get("reason") or "unknown")] += 1
    return {
        "matchups": len(audits),
        "ok_matchups": sum(str(item.get("status")).upper() == "OK" for item in audits),
        "sample_count": sum(int(item.get("sample_count", 0)) for item in audits),
        "attempted_scenes": sum(int(item.get("attempted_scenes", 0)) for item in audits),
        "successful_scenes": sum(int(item.get("successful_scenes", 0)) for item in audits),
        "rejected_scenes": sum(len(item.get("errors", [])) for item in audits),
        "uses_aeronet": False,
        "error_types": dict(errors),
    }


def _history_application_audit(audit_dir: Path) -> dict[str, Any]:
    scene_aot: list[float] = []
    mapping_applied_scenes = 0
    correction: dict[str, dict[str, list[float]]] = {
        band: defaultdict(list) for band in BAND_NAMES
    }
    for path in sorted(audit_dir.glob("*.json")):
        audit = _json(path)
        for scene in audit.get("scenes", []):
            scene_aot.append(float(scene["maiac_aot"]))
            stats = scene["candidates"]["full_a100_all_cap0p030"]
            mapping_applied = stats.get("mapping_applied")
            mapping_applied_scenes += int(
                bool(stats.get("bands"))
                if mapping_applied is None
                else bool(mapping_applied)
            )
            for band, values in stats.get("bands", {}).items():
                for key in ("median", "median_abs", "p95_abs", "cap_fraction"):
                    correction[band][key].append(float(values[key]))
    aot = np.asarray(scene_aot, dtype=np.float64)
    return {
        "scene_count": len(scene_aot),
        "inside_training_aod_domain": int(np.count_nonzero(aot <= 0.15)),
        "outside_training_aod_domain": int(np.count_nonzero(aot > 0.15)),
        "outside_training_aod_fraction": float(np.mean(aot > 0.15)),
        "mapping_applied_scenes": mapping_applied_scenes,
        "mapping_skipped_scenes": len(scene_aot) - mapping_applied_scenes,
        "maiac_aot_median": float(np.median(aot)),
        "maiac_aot_p95": float(np.percentile(aot, 95)),
        "per_band_median_scene_statistics": {
            band: {
                key: float(np.median(values))
                for key, values in band_values.items()
            }
            for band, band_values in correction.items()
        },
    }


def _history_case_audits() -> dict[str, dict[str, Any]]:
    output: dict[str, dict[str, Any]] = {}
    audit_dir = CLEAN_MODEL_DIR / "daily_histories/audits"
    for path in sorted(audit_dir.glob("*.json")):
        audit = _json(path)
        scene_aot: list[float] = []
        correction: dict[str, dict[str, list[float]]] = {
            band: defaultdict(list) for band in BAND_NAMES
        }
        for scene in audit.get("scenes", []):
            scene_aot.append(float(scene["maiac_aot"]))
            stats = scene["candidates"]["full_a100_all_cap0p030"]
            for band, values in stats.get("bands", {}).items():
                for key in ("median", "median_abs", "p95_abs", "cap_fraction"):
                    correction[band][key].append(float(values[key]))
        aot = np.asarray(scene_aot, dtype=np.float64)
        output[audit["matchup_id"]] = {
            "scene_count": len(scene_aot),
            "inside_training_aod_domain": int(np.count_nonzero(aot <= 0.15)),
            "outside_training_aod_domain": int(np.count_nonzero(aot > 0.15)),
            "outside_training_aod_fraction": float(np.mean(aot > 0.15)),
            "maiac_aot_median": float(np.median(aot)),
            "maiac_aot_p95": float(np.percentile(aot, 95)),
            "per_band_median_scene_statistics": {
                band: {
                    key: float(np.median(values))
                    for key, values in band_values.items()
                }
                for band, band_values in correction.items()
            },
        }
    return output


def _terrain_history_case_audits(
    audit_dir: Path, candidate_tag: str
) -> dict[str, dict[str, Any]]:
    """Summarize the per-acquisition terrain mapping used by one replay arm."""
    output: dict[str, dict[str, Any]] = {}
    for path in sorted(audit_dir.glob("*.json")):
        audit = _json(path)
        slopes: list[float] = []
        incidence: list[float] = []
        aot: list[float] = []
        correction: dict[str, list[float]] = defaultdict(list)
        for scene in audit.get("scenes", []):
            stats = (scene.get("candidates") or {}).get(candidate_tag) or {}
            terrain = stats.get("terrain") or {}
            slope = _finite(terrain.get("median_slope_deg"))
            local_incidence = _finite(terrain.get("median_incidence_cos"))
            scene_aot = _finite(scene.get("maiac_aot"))
            if None in (slope, local_incidence, scene_aot):
                continue
            slopes.append(float(slope))
            incidence.append(float(local_incidence))
            aot.append(float(scene_aot))
            for band, values in (stats.get("bands") or {}).items():
                value = _finite(values.get("median"))
                if value is not None:
                    correction[band].append(float(value))
        if not slopes:
            continue
        output[str(audit["matchup_id"])] = {
            "scene_count": len(slopes),
            "median_slope_deg": float(np.median(slopes)),
            "median_incidence_cos": float(np.median(incidence)),
            "median_maiac_aot": float(np.median(aot)),
            "per_band_median_correction": {
                band: float(np.median(values)) for band, values in correction.items()
            },
        }
    return output


def _copy_implementation(output: Path) -> list[dict[str, Any]]:
    copied: list[dict[str, Any]] = []
    for source in IMPLEMENTATION_FILES:
        relative = source.relative_to(REPO_ROOT)
        destination = output / "downloads/implementation" / relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, destination)
        payload = destination.read_bytes()
        copied.append(
            {
                "path": str(relative),
                "url": f"downloads/implementation/{relative}",
                "sha256": hashlib.sha256(payload).hexdigest(),
                "bytes": len(payload),
            }
        )
    (output / "downloads/implementation/manifest.json").write_text(
        json.dumps(copied, indent=2) + "\n", encoding="utf-8"
    )
    return copied


def _job_audit() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for job_id in JOB_IDS:
        result = subprocess.run(
            [
                "sacct",
                "-j",
                job_id,
                "--format=JobIDRaw,State,ExitCode,Elapsed",
                "-n",
                "-X",
                "-P",
            ],
            check=False,
            capture_output=True,
            text=True,
        )
        states: Counter[str] = Counter()
        exits: Counter[str] = Counter()
        for line in result.stdout.splitlines():
            fields = line.split("|")
            if len(fields) >= 3 and fields[1]:
                states[fields[1].split("+")[0]] += 1
                exits[fields[2]] += 1
        terminal_bad = sum(
            count
            for state, count in states.items()
            if state
            in {
                "BOOT_FAIL",
                "CANCELLED",
                "DEADLINE",
                "FAILED",
                "NODE_FAIL",
                "OUT_OF_MEMORY",
                "PREEMPTED",
                "TIMEOUT",
            }
        )
        active = sum(
            count
            for state, count in states.items()
            if state in {"CONFIGURING", "PENDING", "RUNNING"}
        )
        if result.returncode != 0 or not states:
            final_state = "UNRESOLVED"
        elif terminal_bad:
            final_state = "FAILED"
        elif active:
            final_state = "ACTIVE"
        elif set(states) == {"COMPLETED"}:
            final_state = "COMPLETE"
        else:
            final_state = "UNRESOLVED"
        repair_job_ids = list(JOB_REPAIRS.get(job_id, ()))
        if terminal_bad and repair_job_ids:
            repaired = True
            for repair_job_id in repair_job_ids:
                repair_result = subprocess.run(
                    [
                        "sacct",
                        "-j",
                        repair_job_id,
                        "--format=State",
                        "-n",
                        "-X",
                        "-P",
                    ],
                    check=False,
                    capture_output=True,
                    text=True,
                )
                repair_states = {
                    line.strip().split("|")[0].split("+")[0]
                    for line in repair_result.stdout.splitlines()
                    if line.strip()
                }
                if repair_result.returncode != 0 or repair_states != {"COMPLETED"}:
                    repaired = False
                    break
            if repaired:
                final_state = "REPAIRED"
        rows.append(
            {
                "job_id": job_id,
                "states": dict(states),
                "exit_codes": dict(exits),
                "repair_job_ids": repair_job_ids,
                "query_ok": result.returncode == 0,
                "final_state": final_state,
            }
        )
    return rows


def _exact_pair_scatter(output: Path) -> list[dict[str, Any]]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm

    pair_paths = sorted(PAIR_DIR.glob("*.npz"))
    if len(pair_paths) != 36:
        raise ValueError(f"expected 36 exact-pair archives, got {len(pair_paths)}")

    count = np.zeros(len(BAND_NAMES), dtype=np.int64)
    sum_x = np.zeros(len(BAND_NAMES), dtype=np.float64)
    sum_y = np.zeros(len(BAND_NAMES), dtype=np.float64)
    sum_xx = np.zeros(len(BAND_NAMES), dtype=np.float64)
    sum_yy = np.zeros(len(BAND_NAMES), dtype=np.float64)
    sum_xy = np.zeros(len(BAND_NAMES), dtype=np.float64)
    sum_error = np.zeros(len(BAND_NAMES), dtype=np.float64)
    sum_abs_error = np.zeros(len(BAND_NAMES), dtype=np.float64)
    sum_squared_error = np.zeros(len(BAND_NAMES), dtype=np.float64)
    absolute_errors: list[list[np.ndarray]] = [[] for _ in BAND_NAMES]
    plot_x: list[list[np.ndarray]] = [[] for _ in BAND_NAMES]
    plot_y: list[list[np.ndarray]] = [[] for _ in BAND_NAMES]

    for path in pair_paths:
        with np.load(path, allow_pickle=False) as pair:
            names = tuple(str(value) for value in pair["band_names"].tolist())
            if names != BAND_NAMES:
                raise ValueError(f"unexpected band order in {path}: {names}")
            l2a = np.asarray(pair["l2a"], dtype=np.float32)
            target = np.asarray(pair["siac"], dtype=np.float32)
            sample_indices = np.linspace(
                0,
                l2a.shape[0] - 1,
                num=min(3000, l2a.shape[0]),
                dtype=np.int64,
            )
            for band_index in range(len(BAND_NAMES)):
                valid = np.isfinite(l2a[:, band_index]) & np.isfinite(
                    target[:, band_index]
                )
                x = np.asarray(l2a[valid, band_index], dtype=np.float64)
                y = np.asarray(target[valid, band_index], dtype=np.float64)
                error = x - y
                count[band_index] += x.size
                sum_x[band_index] += np.sum(x, dtype=np.float64)
                sum_y[band_index] += np.sum(y, dtype=np.float64)
                sum_xx[band_index] += np.dot(x, x)
                sum_yy[band_index] += np.dot(y, y)
                sum_xy[band_index] += np.dot(x, y)
                sum_error[band_index] += np.sum(error, dtype=np.float64)
                sum_abs_error[band_index] += np.sum(
                    np.abs(error), dtype=np.float64
                )
                sum_squared_error[band_index] += np.dot(error, error)
                absolute_errors[band_index].append(
                    np.asarray(np.abs(error), dtype=np.float32)
                )

                sample_valid = np.isfinite(l2a[sample_indices, band_index]) & np.isfinite(
                    target[sample_indices, band_index]
                )
                selected = sample_indices[sample_valid]
                plot_x[band_index].append(l2a[selected, band_index])
                plot_y[band_index].append(target[selected, band_index])

    metrics: list[dict[str, Any]] = []
    for band_index, band in enumerate(DISPLAY_BANDS):
        n = int(count[band_index])
        centered_xx = sum_xx[band_index] - sum_x[band_index] ** 2 / n
        centered_yy = sum_yy[band_index] - sum_y[band_index] ** 2 / n
        centered_xy = sum_xy[band_index] - sum_x[band_index] * sum_y[band_index] / n
        slope = centered_xy / centered_xx
        intercept = sum_y[band_index] / n - slope * sum_x[band_index] / n
        correlation = centered_xy / math.sqrt(centered_xx * centered_yy)
        errors = np.concatenate(absolute_errors[band_index])
        metrics.append(
            {
                "band": band,
                "pair_count": n,
                "target_on_l2a_slope": float(slope),
                "target_on_l2a_intercept": float(intercept),
                "pearson_r": float(correlation),
                "r_squared": float(correlation**2),
                "raw_bias_l2a_minus_target": float(sum_error[band_index] / n),
                "raw_mae": float(sum_abs_error[band_index] / n),
                "raw_rmse": float(math.sqrt(sum_squared_error[band_index] / n)),
                "raw_p95_abs": float(np.percentile(errors, 95)),
                "plot_sample_count": int(
                    sum(values.size for values in plot_x[band_index])
                ),
            }
        )

    fig, axes = plt.subplots(2, 4, figsize=(16.0, 8.8), facecolor="white")
    for band_index, (ax, row) in enumerate(
        zip(axes.flat[:7], metrics, strict=True)
    ):
        x = np.concatenate(plot_x[band_index]).astype(np.float64)
        y = np.concatenate(plot_y[band_index]).astype(np.float64)
        upper = min(
            0.8,
            max(0.08, 1.06 * float(np.percentile(np.concatenate((x, y)), 99.7))),
        )
        ax.hexbin(
            x,
            y,
            gridsize=72,
            extent=(0.0, upper, 0.0, upper),
            mincnt=1,
            norm=LogNorm(vmin=1, vmax=500),
            cmap="viridis",
            linewidths=0,
        )
        axis = np.asarray([0.0, upper])
        ax.plot(axis, axis, color="#222222", linewidth=1.0, label="1:1")
        ax.plot(
            axis,
            row["target_on_l2a_intercept"] + row["target_on_l2a_slope"] * axis,
            color="#a46b19",
            linewidth=1.1,
            label="OLS",
        )
        ax.set(
            title=row["band"],
            xlabel="Operational L2A reflectance",
            ylabel="Current-RT L1C-corrected reflectance",
            xlim=(0.0, upper),
            ylim=(0.0, upper),
        )
        ax.set_aspect("equal", adjustable="box")
        ax.grid(alpha=0.12)
        ax.tick_params(labelsize=8)
        ax.xaxis.label.set_size(8)
        ax.yaxis.label.set_size(8)
        ax.text(
            0.035,
            0.965,
            (
                f"slope {row['target_on_l2a_slope']:.3f}  "
                f"r {row['pearson_r']:.3f}\n"
                f"bias {row['raw_bias_l2a_minus_target']:+.4f}  "
                f"MAE {row['raw_mae']:.4f}"
            ),
            transform=ax.transAxes,
            va="top",
            fontsize=7.5,
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.78},
        )

    audit = _pair_audit()
    axes.flat[7].axis("off")
    axes.flat[7].text(
        0.04,
        0.94,
        "Exact-pair audit",
        transform=axes.flat[7].transAxes,
        va="top",
        fontsize=12,
        fontweight="bold",
    )
    axes.flat[7].text(
        0.04,
        0.82,
        (
            f"{audit['sample_count']:,} paired pixels\n"
            f"{audit['successful_scenes']:,} historical scenes\n"
            f"{audit['ok_matchups']}/36 target matchups\n\n"
            "Metrics use every finite pair.\n"
            "Density uses <=3,000 rows per case.\n"
            "Purple to yellow: increasing log density.\n"
            "Bias is L2A minus current-RT target.\n"
            "Black: identity; gold: OLS target on L2A."
        ),
        transform=axes.flat[7].transAxes,
        va="top",
        fontsize=9,
        linespacing=1.55,
    )
    fig.suptitle(
        "Exact same-day L2A and MAIAC-conditioned current-RT surface pairs",
        x=0.035,
        y=0.99,
        ha="left",
    )
    fig.subplots_adjust(
        left=0.055,
        right=0.98,
        top=0.91,
        bottom=0.075,
        hspace=0.38,
        wspace=0.30,
    )
    fig.savefig(output / "assets/exact-pair-scatter.png", dpi=145, bbox_inches="tight")
    plt.close(fig)
    return metrics


def _surface_figures(
    output: Path,
    metric_rows: list[dict[str, Any]],
    all_scenes: list[dict[str, str]],
    clean_scenes: list[dict[str, str]],
    site_rows: list[dict[str, Any]],
    coordinates: dict[str, tuple[float, float]],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    pair_scatter_metrics = _exact_pair_scatter(output)
    figures: list[dict[str, Any]] = [
        {
            "path": "assets/exact-pair-scatter.png",
            "title": "Direct L2A-to-current-RT pair scatter",
            "caption": "Exact same-day surface-reflectance pairs for all seven bands. Metrics use every finite pair; density is a deterministic equal-case sample. The identity and OLS lines expose correlation separately from bias and scale mismatch.",
            "wide": True,
        }
    ]
    colors = {
        "Raw operational L2A": "#6c777a",
        "Legacy AOD offset": "#a46b19",
        "Held-site full ridge, alpha 100, cap 0.03": "#176d68",
    }
    methods = tuple(colors)
    x = np.arange(len(DISPLAY_BANDS))
    fig, axes = plt.subplots(2, 1, figsize=(14.5, 9.2), facecolor="white", sharex=True)
    width = 0.23
    for row_index, (training, ax) in enumerate(
        zip(("All retained AOD", "Clean MAIAC AOD <=0.15"), axes, strict=True)
    ):
        for index, method in enumerate(methods):
            selected = {
                row["band"]: row
                for row in metric_rows
                if row["training_set"] == training and row["method"] == method
            }
            values = [selected[band]["scene_mae"] for band in DISPLAY_BANDS]
            ax.bar(
                x + (index - 1) * width,
                values,
                width=width,
                color=colors[method],
                label=method,
            )
        ax.set_ylabel("Held-site scene MAE")
        ax.set_title(training, loc="left", fontsize=11)
        ax.grid(axis="y", alpha=0.2)
        if row_index == 0:
            ax.legend(fontsize=8, ncol=3)
    axes[-1].set_xticks(x, DISPLAY_BANDS)
    fig.suptitle("Surface-reflectance closure by Sentinel-2 band", x=0.02, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    path = output / "assets/surface-band-performance.png"
    fig.savefig(path, dpi=145, bbox_inches="tight")
    plt.close(fig)
    figures.append(
        {
            "path": "assets/surface-band-performance.png",
            "title": "Per-band held-site surface error",
            "caption": "Scene-median error for raw L2A, the legacy AOD offset, and the capped conditional mapping in both training domains.",
        }
    )

    fig, axes = plt.subplots(2, 1, figsize=(14.5, 8.8), facecolor="white", sharex=True)
    for training, ax in zip(
        ("All retained AOD", "Clean MAIAC AOD <=0.15"), axes, strict=True
    ):
        for method in methods:
            selected = {
                row["band"]: row
                for row in metric_rows
                if row["training_set"] == training and row["method"] == method
            }
            ax.plot(
                x,
                [selected[band]["scene_bias"] for band in DISPLAY_BANDS],
                "o-",
                color=colors[method],
                label=method,
            )
        ax.axhline(0.0, color="#222222", linewidth=0.8)
        ax.set_ylabel("Held-site scene bias")
        ax.set_title(training, loc="left", fontsize=11)
        ax.grid(alpha=0.2)
    axes[0].legend(fontsize=8, ncol=3)
    axes[-1].set_xticks(x, DISPLAY_BANDS)
    fig.suptitle("Signed surface-reflectance closure", x=0.02, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    path = output / "assets/surface-band-bias.png"
    fig.savefig(path, dpi=145, bbox_inches="tight")
    plt.close(fig)
    figures.append(
        {
            "path": "assets/surface-band-bias.png",
            "title": "Signed residual by band",
            "caption": "Positive means the tested prior is brighter than same-day L1C corrected with MAIAC and the current libRadTran LUT.",
        }
    )

    all_by_key = {_scene_key(row): row for row in all_scenes}
    merged = [(all_by_key[_scene_key(row)], row) for row in clean_scenes if _scene_key(row) in all_by_key]
    fig, axes = plt.subplots(1, 3, figsize=(15.2, 5.3), facecolor="white", sharey=True)
    for ax, band, label in zip(axes, ("blue", "green", "red"), BANDS, strict=True):
        distributions = [
            [float(clean[f"identity_{band}_bias"]) for _all, clean in merged],
            [float(all_row[f"full_a100_cap_0p030_{band}_bias"]) for all_row, _clean in merged],
            [float(clean[f"full_a100_cap_0p030_{band}_bias"]) for _all, clean in merged],
        ]
        box = ax.boxplot(
            distributions,
            tick_labels=("Raw L2A", "All-AOD map", "Clean-day map"),
            showfliers=False,
            patch_artist=True,
        )
        for patch, color in zip(box["boxes"], ("#6c777a", "#28689b", "#176d68"), strict=True):
            patch.set_facecolor(color)
            patch.set_alpha(0.75)
        ax.axhline(0.0, color="#222222", linewidth=0.8)
        ax.set_title(label)
        ax.tick_params(axis="x", rotation=20)
        ax.grid(axis="y", alpha=0.2)
    axes[0].set_ylabel("Scene bias: prior - current-RT target")
    fig.suptitle("Like-for-like clean-domain out-of-fold residuals", x=0.02, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    path = output / "assets/clean-domain-residuals.png"
    fig.savefig(path, dpi=145, bbox_inches="tight")
    plt.close(fig)
    figures.append(
        {
            "path": "assets/clean-domain-residuals.png",
            "title": "Like-for-like clean-scene residuals",
            "caption": "All-AOD and clean-day models evaluated on the same retained clean scenes, separated by solver band.",
        }
    )

    fig, axes = plt.subplots(2, 2, figsize=(14.5, 9.2), facecolor="white")
    delta = np.asarray(
        [float(row["delta_aot_maiac_minus_sen2cor"]) for row in clean_scenes]
    )
    maiac = np.asarray([float(row["maiac_aot"]) for row in clean_scenes])
    raw = np.asarray([float(row["identity_blue_bias"]) for row in clean_scenes])
    mapped = np.asarray(
        [float(row["full_a100_cap_0p030_blue_bias"]) for row in clean_scenes]
    )
    scatter = axes[0, 0].scatter(delta, raw, c=maiac, cmap="viridis", s=9, alpha=0.45)
    axes[0, 0].axhline(0, color="#222222", linewidth=0.8)
    axes[0, 0].set(title="Raw B02 residual", xlabel="MAIAC - Sen2Cor AOD", ylabel="Reflectance residual")
    fig.colorbar(scatter, ax=axes[0, 0], label="MAIAC AOD550")
    axes[0, 1].scatter(delta, mapped, c=maiac, cmap="viridis", s=9, alpha=0.45)
    axes[0, 1].axhline(0, color="#222222", linewidth=0.8)
    axes[0, 1].set(title="Harmonized B02 residual", xlabel="MAIAC - Sen2Cor AOD", ylabel="Reflectance residual")
    axes[1, 0].hist(maiac, bins=24, color="#28689b", alpha=0.85)
    axes[1, 0].axvline(0.15, color="#a34436", linestyle="--")
    axes[1, 0].set(title="Clean training AOD support", xlabel="MAIAC AOD550", ylabel="Scenes")
    axes[1, 1].hist(delta, bins=28, color="#176d68", alpha=0.85)
    axes[1, 1].axvline(0.0, color="#222222", linewidth=0.8)
    axes[1, 1].set(title="Atmospheric-product disagreement", xlabel="MAIAC - Sen2Cor AOD", ylabel="Scenes")
    for ax in axes.ravel():
        ax.grid(alpha=0.15)
    fig.suptitle("Atmospheric conditioning domain", x=0.02, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    path = output / "assets/atmospheric-conditioning.png"
    fig.savefig(path, dpi=145, bbox_inches="tight")
    plt.close(fig)
    figures.append(
        {
            "path": "assets/atmospheric-conditioning.png",
            "title": "Atmospheric conditioning and residual response",
            "caption": "The clean model is conditioned on the disagreement between MAIAC and Sen2Cor AOD/WV, geometry, elevation, month, sensor, and processing baseline.",
        }
    )

    fig, axes = plt.subplots(1, 2, figsize=(15.2, 5.7), facecolor="white")
    points = [row for row in site_rows if row["site"] in coordinates]
    scatter = axes[0].scatter(
        [coordinates[row["site"]][0] for row in points],
        [coordinates[row["site"]][1] for row in points],
        s=[30 + 2.3 * row["scene_count"] for row in points],
        c=[row["harmonized_visible_mae"] for row in points],
        cmap="magma",
        edgecolor="white",
        linewidth=0.5,
    )
    axes[0].set(title="Clean-pair site coverage", xlabel="Longitude", ylabel="Latitude")
    axes[0].grid(alpha=0.15)
    fig.colorbar(scatter, ax=axes[0], label="Visible scene MAE")
    ordered = sorted(points, key=lambda row: row["harmonized_visible_mae"], reverse=True)
    axes[1].barh(
        np.arange(len(ordered)),
        [row["harmonized_visible_mae"] for row in ordered],
        color="#176d68",
    )
    axes[1].set_yticks(np.arange(len(ordered)), [row["site"] for row in ordered], fontsize=7)
    axes[1].invert_yaxis()
    axes[1].set(title="Site-level harmonized visible error", xlabel="Weighted visible MAE")
    axes[1].grid(axis="x", alpha=0.2)
    fig.tight_layout()
    path = output / "assets/training-coverage.png"
    fig.savefig(path, dpi=145, bbox_inches="tight")
    plt.close(fig)
    figures.append(
        {
            "path": "assets/training-coverage.png",
            "title": "Spatial and site coverage",
            "caption": "Marker size is retained clean-scene count; colour and bars show held-site visible-band closure.",
        }
    )
    return figures, pair_scatter_metrics


def _terrain_surface_ablation(output: Path) -> dict[str, Any] | None:
    """Compare held-site closure with and without local terrain features.

    The two candidates are fitted from the same terrain-enabled pair archive;
    this deliberately keeps the sampled pixels and atmospheric inputs fixed.
    """
    metrics_path = TERRAIN_ALL_MODEL_DIR / "surface_metrics.json"
    scenes_path = TERRAIN_ALL_MODEL_DIR / "surface_scene_metrics.csv"
    model_path = TERRAIN_ALL_MODEL_DIR / "harmonizer.json"
    if not all(path.exists() for path in (metrics_path, scenes_path, model_path)):
        return None

    metrics = _json(metrics_path)
    candidates = metrics.get("candidates") or {}
    if "full_a100" not in candidates or "terrain_a100" not in candidates:
        return None
    baseline = candidates["full_a100"]["cap_0.030"]
    terrain = candidates["terrain_a100"]["cap_0.030"]
    band_rows: list[dict[str, Any]] = []
    display_by_band = dict(zip(BAND_NAMES, DISPLAY_BANDS, strict=True))
    for band in BAND_NAMES:
        base_values = baseline[band]
        terrain_values = terrain[band]
        band_rows.append(
            {
                "band": display_by_band[band],
                "baseline_scene_bias": base_values["scene_bias"],
                "terrain_scene_bias": terrain_values["scene_bias"],
                "scene_bias_delta": terrain_values["scene_bias"]
                - base_values["scene_bias"],
                "baseline_scene_mae": base_values["scene_mae"],
                "terrain_scene_mae": terrain_values["scene_mae"],
                "scene_mae_delta": terrain_values["scene_mae"]
                - base_values["scene_mae"],
                "baseline_pixel_mae": base_values["mae"],
                "terrain_pixel_mae": terrain_values["mae"],
                "pixel_mae_delta": terrain_values["mae"] - base_values["mae"],
            }
        )

    scenes: list[dict[str, Any]] = []
    for row in _read_csv(scenes_path):
        values: dict[str, float] = {}
        for band in ("blue", "green", "red"):
            baseline_mae = _finite(row.get(f"full_a100_cap_0p030_{band}_mae"))
            terrain_mae = _finite(row.get(f"terrain_a100_cap_0p030_{band}_mae"))
            baseline_bias = _finite(row.get(f"full_a100_cap_0p030_{band}_bias"))
            terrain_bias = _finite(row.get(f"terrain_a100_cap_0p030_{band}_bias"))
            if None in (baseline_mae, terrain_mae, baseline_bias, terrain_bias):
                break
            values[f"baseline_{band}_mae"] = float(baseline_mae)
            values[f"terrain_{band}_mae"] = float(terrain_mae)
            values[f"{band}_mae_delta"] = float(terrain_mae - baseline_mae)
            values[f"baseline_{band}_bias"] = float(baseline_bias)
            values[f"terrain_{band}_bias"] = float(terrain_bias)
            values[f"{band}_bias_delta"] = float(terrain_bias - baseline_bias)
        else:
            sza = _finite(row.get("sza_deg"))
            incidence = _finite(row.get("terrain_incidence_cos_mean"))
            slope = _finite(row.get("terrain_slope_deg_mean"))
            elevation = _finite(row.get("terrain_elevation_km_mean"))
            maiac_aot = _finite(row.get("maiac_aot"))
            if None in (sza, incidence, slope, elevation, maiac_aot):
                continue
            baseline_visible = float(
                np.mean([values[f"baseline_{band}_mae"] for band in ("blue", "green", "red")])
            )
            terrain_visible = float(
                np.mean([values[f"terrain_{band}_mae"] for band in ("blue", "green", "red")])
            )
            scenes.append(
                {
                    "matchup_id": row["matchup_id"],
                    "site": row["site"],
                    "window": row["window"],
                    "day": row["day"],
                    "scene_id": row["scene_id"],
                    "sample_count": int(float(row["sample_count"])),
                    "maiac_aot": float(maiac_aot),
                    "terrain_elevation_km_mean": float(elevation),
                    "terrain_slope_deg_mean": float(slope),
                    "terrain_incidence_cos_mean": float(incidence),
                    "flat_incidence_cos": float(np.cos(np.radians(sza))),
                    "terrain_incidence_delta": float(
                        incidence - np.cos(np.radians(sza))
                    ),
                    "baseline_visible_mae": baseline_visible,
                    "terrain_visible_mae": terrain_visible,
                    "visible_mae_delta": terrain_visible - baseline_visible,
                    **values,
                }
            )
    if not scenes:
        return None

    grouped: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in scenes:
        grouped[row["site"]].append(row)
    site_rows: list[dict[str, Any]] = []
    for site, rows in sorted(grouped.items()):
        weights = np.asarray([row["sample_count"] for row in rows], dtype=np.float64)
        base_visible = np.asarray(
            [row["baseline_visible_mae"] for row in rows], dtype=np.float64
        )
        terrain_visible = np.asarray(
            [row["terrain_visible_mae"] for row in rows], dtype=np.float64
        )
        site_rows.append(
            {
                "site": site,
                "scene_count": len(rows),
                "sample_count": int(weights.sum()),
                "baseline_visible_mae": float(np.average(base_visible, weights=weights)),
                "terrain_visible_mae": float(
                    np.average(terrain_visible, weights=weights)
                ),
                "visible_mae_delta": float(
                    np.average(terrain_visible - base_visible, weights=weights)
                ),
            }
        )

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(15.6, 10.4), facecolor="white")
    solver_rows = [row for row in band_rows if row["band"] in ("B02", "B03", "B04")]
    x = np.arange(len(solver_rows))
    width = 0.36
    axes[0, 0].bar(
        x - width / 2,
        [row["baseline_scene_mae"] for row in solver_rows],
        width=width,
        color="#6c777a",
        label="Atmosphere/geometry ridge",
    )
    axes[0, 0].bar(
        x + width / 2,
        [row["terrain_scene_mae"] for row in solver_rows],
        width=width,
        color="#176d68",
        label="Ridge + GLO-30 terrain",
    )
    axes[0, 0].set(
        title="Held-site solver-band scene MAE",
        xlabel="Sentinel-2 band",
        ylabel="Surface-reflectance MAE",
    )
    axes[0, 0].set_xticks(x, [row["band"] for row in solver_rows])
    axes[0, 0].grid(axis="y", alpha=0.18)
    axes[0, 0].legend(fontsize=8)

    slope = np.asarray([row["terrain_slope_deg_mean"] for row in scenes])
    visible_delta = np.asarray([row["visible_mae_delta"] for row in scenes])
    incidence_delta = np.asarray(
        [row["terrain_incidence_delta"] for row in scenes]
    )
    scatter = axes[0, 1].scatter(
        slope,
        visible_delta,
        c=incidence_delta,
        cmap="coolwarm",
        vmin=-0.16,
        vmax=0.16,
        s=22,
        alpha=0.78,
    )
    axes[0, 1].axhline(0.0, color="#222222", linewidth=0.8)
    axes[0, 1].set(
        title="Terrain response by local slope",
        xlabel="Mean local slope (degrees)",
        ylabel="Terrain - baseline visible MAE",
    )
    axes[0, 1].grid(alpha=0.18)
    fig.colorbar(scatter, ax=axes[0, 1], label="Local - flat solar incidence cosine")

    aot = np.asarray([row["maiac_aot"] for row in scenes])
    blue_bias_delta = np.asarray([row["blue_bias_delta"] for row in scenes])
    scatter = axes[1, 0].scatter(
        incidence_delta,
        blue_bias_delta,
        c=aot,
        cmap="viridis",
        s=22,
        alpha=0.78,
    )
    axes[1, 0].axhline(0.0, color="#222222", linewidth=0.8)
    axes[1, 0].axvline(0.0, color="#222222", linewidth=0.8)
    axes[1, 0].set(
        title="B02 bias response by local illumination",
        xlabel="Local - flat solar incidence cosine",
        ylabel="Terrain - baseline B02 scene bias",
    )
    axes[1, 0].grid(alpha=0.18)
    fig.colorbar(scatter, ax=axes[1, 0], label="MAIAC AOD550")

    ordered_sites = sorted(site_rows, key=lambda row: row["visible_mae_delta"])
    site_delta = [row["visible_mae_delta"] for row in ordered_sites]
    axes[1, 1].barh(
        np.arange(len(ordered_sites)),
        site_delta,
        color=["#176d68" if value <= 0.0 else "#a34436" for value in site_delta],
    )
    axes[1, 1].axvline(0.0, color="#222222", linewidth=0.8)
    axes[1, 1].set(
        title="Weighted visible-MAE change by site",
        xlabel="Terrain - baseline visible MAE",
    )
    axes[1, 1].set_yticks(
        np.arange(len(ordered_sites)),
        [row["site"] for row in ordered_sites],
        fontsize=7,
    )
    axes[1, 1].grid(axis="x", alpha=0.18)
    fig.suptitle(
        "Terrain-conditioned surface-prior ablation", x=0.035, ha="left"
    )
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    figure_path = output / "assets/terrain-surface-ablation.png"
    fig.savefig(figure_path, dpi=145, bbox_inches="tight")
    plt.close(fig)

    model = _json(model_path)
    return {
        "pair_audit": {
            "matchups": len(list(TERRAIN_PAIR_DIR.glob("*.json"))),
            "source": model.get("terrain_features", {}).get("source"),
            "fields": model.get("terrain_features", {}).get("fields", []),
        },
        "model": {
            "baseline_specification": "full_a100",
            "terrain_specification": "terrain_a100",
            "feature_set": model["models"]["terrain_a100"]["feature_set"],
            "ridge_alpha": model["models"]["terrain_a100"]["alpha"],
            "correction_cap": 0.03,
            "uses_aeronet": model["uses_aeronet"],
        },
        "band_rows": band_rows,
        "scene_count": len(scenes),
        "site_count": len(site_rows),
        "figure": "assets/terrain-surface-ablation.png",
        "scene_rows": scenes,
        "site_rows": site_rows,
    }


def _candidate_cost_figure(
    path: Path,
    *,
    cost: np.lib.npyio.NpzFile,
    record: dict[str, Any],
    local_mask: np.ndarray,
    title_prefix: str = "Clean-day harmonized",
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    axis = np.asarray(cost["aot_axis"], dtype=np.float64)
    valid = np.asarray(cost["solve_valid"], dtype=bool)
    scope = local_mask & valid
    if np.count_nonzero(scope) < 20:
        scope = valid
    total = _median_curve(np.asarray(cost["cube"], dtype=np.float64), scope)
    band_cost = np.asarray(cost["band_cost_cube"], dtype=np.float64)
    band_signed = np.asarray(cost["band_signed_residual_cube"], dtype=np.float64)
    per_band = [_median_curve(band_cost[index], scope) for index in range(3)]
    signed = [_median_curve(band_signed[index], scope) for index in range(3)]
    truth = float(record["truth"])
    retrieved = float(record["retrieved"])
    ee = 0.05 + 0.15 * truth
    atmo = float(np.nanmedian(np.asarray(cost["aot_prior"])[scope]))
    toa = np.asarray(cost["toa"], dtype=np.float64)
    prior = np.asarray(cost["boa_prior"], dtype=np.float64)
    prior_unc = np.asarray(cost["boa_unc"], dtype=np.float64)
    toa_median = np.asarray([np.nanmedian(toa[index][scope]) for index in range(3)])
    prior_median = np.asarray([np.nanmedian(prior[index][scope]) for index in range(3)])
    uncertainty = np.asarray([np.nanmedian(prior_unc[index][scope]) for index in range(3)])
    solved, _uncertainty, _cost = _pooled_maps(cost)

    def add_context(ax: Any) -> None:
        ax.axvspan(truth - ee, truth + ee, color="#e5eaeb", alpha=0.9)
        ax.axvline(truth, color="#161d20", linewidth=1.4, label="AERONET")
        ax.axvline(retrieved, color="#28689b", linestyle="--", label="Harmonized retrieval")
        ax.axvline(atmo, color="#a34436", linestyle=":", label="MAIAC")

    fig, axes = plt.subplots(2, 3, figsize=(16.2, 9.1), facecolor="white")
    axes = axes.ravel()
    axes[0].plot(axis, _normalise_curve(total), color="#28689b")
    add_context(axes[0])
    axes[0].set(title="Total station-window likelihood", xlabel="AOD550", ylabel="log10(1 + delta cost)")
    axes[0].legend(fontsize=7)
    for index, band in enumerate(BANDS):
        axes[1].plot(axis, _normalise_curve(per_band[index]), color=BAND_COLORS[band], label=band)
    add_context(axes[1])
    axes[1].set(title="Per-band likelihood", xlabel="AOD550", ylabel="log10(1 + delta cost)")
    axes[1].legend(fontsize=7)
    for index, band in enumerate(BANDS):
        axes[2].plot(axis, signed[index], color=BAND_COLORS[band], label=band)
    add_context(axes[2])
    axes[2].axhline(0.0, color="#777777", linewidth=0.8)
    axes[2].set(title="Signed BOA - prior residual", xlabel="AOD550", ylabel="Reflectance")
    axes[2].legend(fontsize=7)
    wavelengths = np.asarray([WAVELENGTHS[band] for band in BANDS])
    axes[3].plot(wavelengths, toa_median, "o-", color="#161d20", label="TOA")
    axes[3].errorbar(wavelengths, prior_median, yerr=uncertainty, fmt="o-", color="#176d68", capsize=4, label="Harmonized surface prior")
    axes[3].set_xticks(wavelengths, BANDS)
    axes[3].set(title="Station-window spectrum", ylabel="Reflectance")
    axes[3].legend(fontsize=7)
    image = axes[4].imshow(solved, cmap="viridis", vmin=0.0, vmax=max(1.0, truth + ee))
    axes[4].set_title("Solved AOD map")
    axes[4].set_xticks([])
    axes[4].set_yticks([])
    fig.colorbar(image, ax=axes[4], fraction=0.046, pad=0.03)
    difference = solved - np.asarray(cost["aot_prior"], dtype=np.float64)
    image = axes[5].imshow(difference, cmap="RdBu_r", vmin=-0.4, vmax=0.4)
    axes[5].set_title("Solved AOD - MAIAC")
    axes[5].set_xticks([])
    axes[5].set_yticks([])
    fig.colorbar(image, ax=axes[5], fraction=0.046, pad=0.03)
    for ax in axes[:3]:
        ax.set_xlim(float(axis[0]), min(float(axis[-1]), max(1.25, 1.5 * truth)))
        ax.grid(alpha=0.18)
    fig.suptitle(
        f"{record['matchup_id']} | {title_prefix} spectral and cost evidence",
        x=0.02,
        ha="left",
    )
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=125, bbox_inches="tight")
    plt.close(fig)


def _retrieval_figures(
    output: Path,
    cases: list[dict[str, Any]],
    variants: list[dict[str, Any]],
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    complete = sorted(
        (row for row in variants if row["complete"]),
        key=lambda row: row["metrics"]["within_ee_rate"],
        reverse=True,
    )
    fig, ax = plt.subplots(figsize=(13.0, 8.0), facecolor="white")
    labels = [row["label"] for row in complete]
    rates = [100.0 * row["metrics"]["within_ee_rate"] for row in complete]
    ax.barh(
        np.arange(len(labels)),
        rates,
        color=["#176d68" if value >= 87 else "#28689b" for value in rates],
    )
    ax.axvline(87.0, color="#a34436", linewidth=1.5, label="Target")
    ax.set_yticks(np.arange(len(labels)), labels, fontsize=8)
    ax.invert_yaxis()
    ax.set(xlabel="Within expected error (%)", xlim=(0, 100), title="Complete generic development experiments")
    ax.grid(axis="x", alpha=0.2)
    ax.legend()
    fig.tight_layout()
    fig.savefig(output / "assets/variant-ranking.png", dpi=145, bbox_inches="tight")
    plt.close(fig)

    complete_ids = {row["id"] for row in variants if row["complete"]}
    completed_terrain = [
        row
        for row in variants
        if row["id"] in TERRAIN_ARTIFACTS and row["complete"]
    ]
    terrain_id = (
        max(
            completed_terrain,
            key=lambda row: row["metrics"]["within_ee_rate"],
        )["id"]
        if completed_terrain
        else "harmonized"
    )
    selected_ids = ("fresh", "current_best", "all_aod_live", terrain_id)
    fig, axes = plt.subplots(1, 4, figsize=(19.0, 5.2), facecolor="white", sharex=True, sharey=True)
    truth = np.asarray([case["truth"] for case in cases])
    for ax, variant_id in zip(axes, selected_ids, strict=True):
        values = np.asarray([case["candidate_by_id"][variant_id]["value"] for case in cases], dtype=float)
        hit = np.asarray([case["candidate_by_id"][variant_id]["within_ee"] for case in cases], dtype=bool)
        grid = np.linspace(0.2, 0.9, 200)
        ax.fill_between(grid, grid - (0.05 + 0.15 * grid), grid + (0.05 + 0.15 * grid), color="#e5eaeb")
        ax.plot(grid, grid, color="#161d20", linewidth=1)
        ax.scatter(truth[hit], values[hit], s=28, color="#176d68")
        ax.scatter(truth[~hit], values[~hit], s=28, color="#a34436")
        label = next(row["label"] for row in variants if row["id"] == variant_id)
        ax.set(title=label, xlabel="AERONET AOD550", xlim=(0.2, 0.9), ylim=(0.0, 1.05))
        ax.grid(alpha=0.15)
    axes[0].set_ylabel("Retrieved AOD550")
    fig.tight_layout()
    fig.savefig(output / "assets/retrieval-scatter.png", dpi=145, bbox_inches="tight")
    plt.close(fig)

    matrix_ids = tuple(
        variant_id
        for variant_id in (
            "fresh",
            "current_best",
            "identity",
            "legacy_offset",
            "all_aod_live",
            "all_aod",
            "terrain_a100_cap015",
            "terrain_a100_cap030",
            "terrain_blue_cap015",
            "terrain_blue_cap030",
            "terrain_clean_a100_cap030",
            "terrain_clean_blue_cap030",
            "paired_exact",
            "scene_mapping",
            "harmonized",
        )
        if variant_id in complete_ids
    )
    ordered = sorted(cases, key=lambda case: case["fresh_error_over_ee"])
    matrix = np.asarray(
        [
            [case["candidate_by_id"][variant_id]["error_over_ee"] for variant_id in matrix_ids]
            for case in ordered
        ],
        dtype=float,
    )
    fig, ax = plt.subplots(figsize=(15.0, 11.0), facecolor="white")
    image = ax.imshow(matrix, cmap="RdBu_r", vmin=-3, vmax=3, aspect="auto")
    ax.set_xticks(
        np.arange(len(matrix_ids)),
        [next(row["label"] for row in variants if row["id"] == item) for item in matrix_ids],
        rotation=30,
        ha="right",
        fontsize=8,
    )
    ax.set_yticks(np.arange(len(ordered)), [case["site"] for case in ordered], fontsize=7)
    ax.set_title("Signed AOD error in expected-error units")
    fig.colorbar(image, ax=ax, label="(retrieved - AERONET) / EE")
    fig.tight_layout()
    fig.savefig(output / "assets/evidence-matrix.png", dpi=145, bbox_inches="tight")
    plt.close(fig)

    reference = np.asarray([case["candidate_by_id"]["fresh"]["error"] for case in cases])
    candidate_error = np.asarray(
        [case["candidate_by_id"][terrain_id]["error"] for case in cases]
    )
    delta = candidate_error - reference
    candidate_label = next(row["label"] for row in variants if row["id"] == terrain_id)
    fig, axes = plt.subplots(1, 2, figsize=(15.5, 5.8), facecolor="white")
    points = axes[0].scatter(
        [case["lon"] for case in cases],
        [case["lat"] for case in cases],
        c=delta,
        cmap="RdBu_r",
        vmin=-0.25,
        vmax=0.25,
        s=75,
        edgecolor="white",
        linewidth=0.5,
    )
    axes[0].set(title="Spatial change from current code", xlabel="Longitude", ylabel="Latitude")
    axes[0].grid(alpha=0.15)
    fig.colorbar(
        points, ax=axes[0], label=f"{candidate_label} error - current error"
    )
    order = np.argsort(delta)
    colors = ["#176d68" if value < 0 else "#a34436" for value in delta[order]]
    axes[1].barh(np.arange(len(cases)), delta[order], color=colors)
    axes[1].set_yticks(np.arange(len(cases)), [cases[index]["site"] for index in order], fontsize=7)
    axes[1].axvline(0, color="#161d20", linewidth=0.8)
    axes[1].set(
        title="Per-case signed change",
        xlabel=f"{candidate_label} error - current error",
    )
    axes[1].grid(axis="x", alpha=0.2)
    fig.tight_layout()
    fig.savefig(output / "assets/spatial-error.png", dpi=145, bbox_inches="tight")
    plt.close(fig)

    colors = {
        "gain": "#176d68",
        "loss": "#a34436",
        "stable_hit": "#28689b",
        "stable_miss": "#7b8588",
    }
    fig, axes = plt.subplots(2, 3, figsize=(16.0, 8.8), facecolor="white")
    retrieval_delta = np.asarray(
        [
            case["candidate_by_id"]["harmonized"]["value"]
            - case["candidate_by_id"]["fresh"]["value"]
            for case in cases
        ]
    )
    absolute_error_change = np.asarray(
        [
            abs(case["candidate_by_id"]["harmonized"]["error"])
            - abs(case["candidate_by_id"]["fresh"]["error"])
            for case in cases
        ]
    )
    point_colors = [colors[case["transition"]] for case in cases]
    for ax, band in zip(axes[0], BANDS, strict=True):
        surface_delta = np.asarray(
            [
                case["harmonized_prior_boa"].get(band, np.nan)
                - case["prior_boa"].get(band, np.nan)
                for case in cases
            ]
        )
        ax.scatter(surface_delta, retrieval_delta, c=point_colors, s=38, alpha=0.85)
        ax.axhline(0, color="#222222", linewidth=0.8)
        ax.axvline(0, color="#222222", linewidth=0.8)
        ax.set(
            title=f"{band} prior response",
            xlabel=f"Delta {band} prior reflectance",
        )
        ax.grid(alpha=0.18)
    axes[0, 0].set_ylabel("Delta retrieved AOD")
    nir_correction = np.asarray(
        [
            case["harmonization_history"]["per_band_median_scene_statistics"]["nir08"]["median"]
            for case in cases
        ]
    )
    axes[1, 0].scatter(nir_correction, retrieval_delta, c=point_colors, s=38, alpha=0.85)
    axes[1, 0].axhline(0, color="#222222", linewidth=0.8)
    axes[1, 0].set(
        title="NIR anchor correction",
        xlabel="Median B8A correction",
        ylabel="Delta retrieved AOD",
    )
    outside_fraction = np.asarray(
        [case["harmonization_history"]["outside_training_aod_fraction"] for case in cases]
    )
    axes[1, 1].scatter(outside_fraction, absolute_error_change, c=point_colors, s=38, alpha=0.85)
    axes[1, 1].axhline(0, color="#222222", linewidth=0.8)
    axes[1, 1].set(
        title="Application-domain exposure",
        xlabel="Fraction of history with AOD > 0.15",
        ylabel="Delta absolute AOD error",
    )
    history_aot = np.asarray(
        [case["harmonization_history"]["maiac_aot_median"] for case in cases]
    )
    axes[1, 2].scatter(history_aot, retrieval_delta, c=point_colors, s=38, alpha=0.85)
    axes[1, 2].axhline(0, color="#222222", linewidth=0.8)
    axes[1, 2].axvline(0.15, color="#a34436", linestyle="--", linewidth=0.9)
    axes[1, 2].set(
        title="Historical aerosol domain",
        xlabel="Median historical MAIAC AOD",
        ylabel="Delta retrieved AOD",
    )
    for ax in axes[1]:
        ax.grid(alpha=0.18)
    for ax in axes.flat:
        ax.tick_params(labelsize=8)
        ax.title.set_fontsize(10)
        ax.xaxis.label.set_size(9)
        ax.yaxis.label.set_size(9)
    handles = [
        plt.Line2D([], [], marker="o", linestyle="", color=color, label=label.replace("_", " "))
        for label, color in colors.items()
    ]
    fig.legend(
        handles=handles,
        loc="lower center",
        bbox_to_anchor=(0.5, 0.015),
        ncol=4,
        fontsize=8,
    )
    fig.suptitle(
        "How surface-prior changes propagate into AOD", x=0.06, y=0.98, ha="left"
    )
    fig.subplots_adjust(
        left=0.07, right=0.985, top=0.90, bottom=0.14, hspace=0.42, wspace=0.25
    )
    fig.savefig(output / "assets/surface-to-aod-response.png", dpi=145, bbox_inches="tight")
    plt.close(fig)


def _terrain_retrieval_figure(
    output: Path,
    cases: list[dict[str, Any]],
    *,
    candidate_id: str,
    candidate_label: str,
) -> dict[str, Any] | None:
    """Visualize terrain-conditioned AOD gains and losses without routing cases."""
    rows: list[dict[str, Any]] = []
    for case in cases:
        terrain = case["candidate_by_id"].get(candidate_id) or {}
        history = case.get("terrain_history") or {}
        value = _finite(terrain.get("value"))
        slope = _finite(history.get("median_slope_deg"))
        incidence = _finite(history.get("median_incidence_cos"))
        b02_correction = _finite(
            (history.get("per_band_median_correction") or {}).get("blue")
        )
        b02_prior_delta = _finite(
            case.get("terrain_prior_boa", {}).get("B02", np.nan)
        )
        reference_b02 = _finite(case.get("prior_boa", {}).get("B02", np.nan))
        if None in (value, slope, incidence, b02_correction, b02_prior_delta, reference_b02):
            continue
        truth = float(case["truth"])
        fresh_value = float(case["candidate_by_id"]["fresh"]["value"])
        fresh_error = fresh_value - truth
        terrain_error = float(value) - truth
        fresh_hit = within_ee(fresh_value, truth)
        terrain_hit = within_ee(float(value), truth)
        transition = (
            "gain"
            if terrain_hit and not fresh_hit
            else "loss"
            if fresh_hit and not terrain_hit
            else "stable_hit"
            if fresh_hit
            else "stable_miss"
        )
        rows.append(
            {
                "matchup_id": case["matchup_id"],
                "site": case["site"],
                "truth": truth,
                "fresh": fresh_value,
                "terrain": float(value),
                "fresh_error": fresh_error,
                "terrain_error": terrain_error,
                "retrieved_aod_delta": float(value) - fresh_value,
                "absolute_error_delta": abs(terrain_error) - abs(fresh_error),
                "fresh_within_ee": fresh_hit,
                "terrain_within_ee": terrain_hit,
                "transition": transition,
                "median_slope_deg": float(slope),
                "median_incidence_cos": float(incidence),
                "median_b02_correction": float(b02_correction),
                "terrain_b02_prior_delta": float(b02_prior_delta - reference_b02),
                "median_history_maiac_aot": _finite(history.get("median_maiac_aot")),
            }
        )
    if len(rows) < 3:
        return None

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    truth = np.asarray([row["truth"] for row in rows])
    terrain = np.asarray([row["terrain"] for row in rows])
    slope = np.asarray([row["median_slope_deg"] for row in rows])
    incidence = np.asarray([row["median_incidence_cos"] for row in rows])
    b02_correction = np.asarray([row["median_b02_correction"] for row in rows])
    prior_delta = np.asarray([row["terrain_b02_prior_delta"] for row in rows])
    aod_delta = np.asarray([row["retrieved_aod_delta"] for row in rows])
    absolute_delta = np.asarray([row["absolute_error_delta"] for row in rows])
    transitions = [row["transition"] for row in rows]
    colors = {
        "gain": "#176d68",
        "loss": "#a34436",
        "stable_hit": "#28689b",
        "stable_miss": "#7b8588",
    }
    point_colors = [colors[value] for value in transitions]
    fig, axes = plt.subplots(2, 3, figsize=(16.5, 10.2), facecolor="white")
    grid = np.linspace(0.2, 0.9, 200)
    hit = np.asarray([row["terrain_within_ee"] for row in rows], dtype=bool)
    axes[0, 0].fill_between(
        grid,
        grid - (0.05 + 0.15 * grid),
        grid + (0.05 + 0.15 * grid),
        color="#e5eaeb",
    )
    axes[0, 0].plot(grid, grid, color="#161d20", linewidth=1)
    axes[0, 0].scatter(truth[hit], terrain[hit], color="#176d68", s=34)
    axes[0, 0].scatter(truth[~hit], terrain[~hit], color="#a34436", s=34)
    axes[0, 0].set(
        title=f"{candidate_label} against AERONET",
        xlabel="AERONET AOD550",
        ylabel="Retrieved AOD550",
        xlim=(0.2, 0.9),
        ylim=(0.0, 1.05),
    )
    axes[0, 0].grid(alpha=0.18)

    axes[0, 1].scatter(slope, absolute_delta, c=point_colors, s=38, alpha=0.88)
    axes[0, 1].axhline(0.0, color="#222222", linewidth=0.8)
    axes[0, 1].set(
        title="AOD absolute-error response by slope",
        xlabel="Median history slope (degrees)",
        ylabel="Terrain - current absolute AOD error",
    )
    axes[0, 1].grid(alpha=0.18)

    axes[0, 2].scatter(incidence, aod_delta, c=point_colors, s=38, alpha=0.88)
    axes[0, 2].axhline(0.0, color="#222222", linewidth=0.8)
    axes[0, 2].set(
        title="AOD response by local illumination",
        xlabel="Median local solar-incidence cosine",
        ylabel="Terrain - current retrieved AOD",
    )
    axes[0, 2].grid(alpha=0.18)

    axes[1, 0].scatter(b02_correction, aod_delta, c=point_colors, s=38, alpha=0.88)
    axes[1, 0].axhline(0.0, color="#222222", linewidth=0.8)
    axes[1, 0].axvline(0.0, color="#222222", linewidth=0.8)
    axes[1, 0].set(
        title="AOD response to mapped B02 correction",
        xlabel="Median terrain B02 correction",
        ylabel="Terrain - current retrieved AOD",
    )
    axes[1, 0].grid(alpha=0.18)

    axes[1, 1].scatter(prior_delta, aod_delta, c=point_colors, s=38, alpha=0.88)
    axes[1, 1].axhline(0.0, color="#222222", linewidth=0.8)
    axes[1, 1].axvline(0.0, color="#222222", linewidth=0.8)
    axes[1, 1].set(
        title="AOD response to B02 prior change",
        xlabel="Terrain - current B02 prior",
        ylabel="Terrain - current retrieved AOD",
    )
    axes[1, 1].grid(alpha=0.18)

    order = np.argsort(absolute_delta)
    axes[1, 2].barh(
        np.arange(len(rows)),
        absolute_delta[order],
        color=["#176d68" if value <= 0.0 else "#a34436" for value in absolute_delta[order]],
    )
    axes[1, 2].axvline(0.0, color="#222222", linewidth=0.8)
    axes[1, 2].set(
        title="Per-case absolute AOD-error change",
        xlabel="Terrain - current absolute AOD error",
    )
    axes[1, 2].set_yticks(
        np.arange(len(rows)), [rows[index]["site"] for index in order], fontsize=7
    )
    axes[1, 2].grid(axis="x", alpha=0.18)
    handles = [
        plt.Line2D([], [], marker="o", linestyle="", color=color, label=label.replace("_", " "))
        for label, color in colors.items()
    ]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, 0.0), ncol=4, fontsize=8)
    fig.suptitle(
        f"{candidate_label}: retrieval gains and losses",
        x=0.04,
        ha="left",
    )
    fig.subplots_adjust(
        left=0.07, right=0.985, top=0.91, bottom=0.10, hspace=0.38, wspace=0.28
    )
    figure_path = output / "assets/terrain-to-aod-response.png"
    fig.savefig(figure_path, dpi=145, bbox_inches="tight")
    plt.close(fig)
    counts = Counter(transitions)
    return {
        "figure": "assets/terrain-to-aod-response.png",
        "candidate_id": candidate_id,
        "candidate_label": candidate_label,
        "case_count": len(rows),
        "transitions": dict(counts),
        "rows": rows,
    }


def _inventory(
    matchup_ids: list[str],
    reference: dict[str, dict[str, Any]],
    curated: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    rows = [{**row, "inventory_kind": "curated"} for row in curated]
    known = {Path(row["path"]).resolve() for row in curated}
    paths = sorted(ROOT.glob("phaseD_results_lowcloud20_l2aharm_*mediumdev_20260713"))
    for path in paths:
        if path.resolve() in known:
            continue
        records, statuses = _load_records(path, matchup_ids)
        if len(records) < 6:
            continue
        rows.append(
            {
                "id": "inventory_" + hashlib.sha1(path.name.encode()).hexdigest()[:10],
                "label": path.name.replace("phaseD_results_lowcloud20_l2aharm_", ""),
                "path": str(path),
                "family": "Surface harmonization diagnostic",
                "description": "Discovered saved harmonization arm; retained without production-equivalence assumptions.",
                "status_counts": dict(statuses),
                "metrics": _metrics(records, matchup_ids),
                "complete": len(records) == len(matchup_ids),
                "transitions_vs_fresh": _transitions(reference, records, matchup_ids),
                "inventory_kind": "discovered",
            }
        )
    return rows


def build(output: Path, *, reuse_figures: bool = False) -> dict[str, Any]:
    base = _json(BASE_REPORT / "data/report.json")
    base_cases = {row["matchup_id"]: row for row in base["cases"]}
    matchup_ids = list(base_cases)
    if len(matchup_ids) != 36:
        raise ValueError(f"expected 36 locked medium-development cases, got {len(matchup_ids)}")

    output.mkdir(parents=True, exist_ok=True)
    for directory in (
        output / "assets",
        output / "assets/harmonized-spatial",
        output / "assets/harmonized-diagnostic",
        output / "assets/terrain-spatial",
        output / "assets/terrain-diagnostic",
        output / "data",
        output / "downloads",
        output / "downloads/implementation",
    ):
        directory.mkdir(parents=True, exist_ok=True)
    shutil.copy2(WEB_ASSETS / "app.css", output / "app.css")
    shutil.copy2(WEB_ASSETS / "app.js", output / "app.js")

    variant_records: dict[str, dict[str, dict[str, Any]]] = {}
    variants: list[dict[str, Any]] = []
    for spec in VARIANT_SPECS:
        records, statuses = _load_records(Path(spec["path"]), matchup_ids)
        variant_records[spec["id"]] = records
        variants.append(
            {
                **spec,
                "path": str(spec["path"]),
                "status_counts": dict(statuses),
                "metrics": _metrics(records, matchup_ids),
                "complete": len(records) == len(matchup_ids),
            }
        )
    reference = variant_records["fresh"]
    harmonized_records = variant_records["harmonized"]
    terrain_variants = [
        row
        for row in variants
        if row["id"] in TERRAIN_ARTIFACTS
    ]
    terrain_cap030 = next(
        row for row in terrain_variants if row["id"] == "terrain_a100_cap030"
    )
    completed_terrain_variants = [
        row for row in terrain_variants if row["complete"]
    ]
    terrain_best = (
        max(
            completed_terrain_variants,
            key=lambda row: row["metrics"]["within_ee_rate"],
        )
        if completed_terrain_variants
        else terrain_cap030
    )
    terrain_complete = bool(completed_terrain_variants)
    terrain_artifacts = TERRAIN_ARTIFACTS[terrain_best["id"]]
    terrain_results = Path(terrain_artifacts["results"])
    terrain_cubes = Path(terrain_artifacts["cubes"])
    terrain_records = variant_records[terrain_best["id"]]
    if len(reference) != 36 or len(harmonized_records) != 36:
        raise ValueError(
            f"reference/harmonized output incomplete: {len(reference)}/36 and {len(harmonized_records)}/36"
        )
    terrain_evidence_available = terrain_complete and len(terrain_records) == len(
        matchup_ids
    ) and all(
        (terrain_cubes / f"{matchup_id}.npz").exists()
        for matchup_id in matchup_ids
    )
    for row in variants:
        row["transitions_vs_fresh"] = _transitions(
            reference, variant_records[row["id"]], matchup_ids
        )

    metric_rows = _surface_metric_rows()
    all_scene_rows = _read_csv(ALL_MODEL_DIR / "surface_scene_metrics.csv")
    clean_scene_rows = _read_csv(CLEAN_MODEL_DIR / "surface_scene_metrics.csv")
    site_rows = _surface_site_rows(clean_scene_rows)
    coordinates: dict[str, tuple[float, float]] = {}
    for row in base["cases"]:
        coordinates.setdefault(row["site"], (float(row["lon"]), float(row["lat"])))
    surface_figures, pair_scatter_metrics = _surface_figures(
        output,
        metric_rows,
        all_scene_rows,
        clean_scene_rows,
        site_rows,
        coordinates,
    )
    terrain_surface_ablation = _terrain_surface_ablation(output)
    if terrain_surface_ablation is not None:
        surface_figures.append(
            {
                "path": terrain_surface_ablation["figure"],
                "title": "Terrain-conditioned held-site ablation",
                "caption": "The baseline and terrain models use the same GLO-30-enabled exact-pair archive and held-site folds. Negative MAE differences favour terrain; the scene and site distributions remain visible instead of being reduced to one aggregate score.",
                "wide": True,
            }
        )
    pair_gallery_path = output / "assets/pair-examples/pair-gallery-metadata.json"
    if not pair_gallery_path.exists():
        raise ValueError(
            "missing spatial-pair gallery; run build_l2a_l1c_pair_gallery.py first"
        )
    pair_examples = _json(pair_gallery_path)
    surface_figures.insert(
        0,
        {
            "path": pair_examples["gallery_image"],
            "title": "Four selected spatial L2A/current-RT image pairs",
            "caption": "Actual same-day 60 m operational L2A and MAIAC-conditioned L1C/current-RT surface pairs. The examples span arid, vegetated and urban settings and were selected for visual diversity, not retrieval outcome. Each source pair shares a contrast stretch; the visible-band difference uses a fixed +/-0.04 scale. Blank pixels fail the exact clear-land mask used for the pair archive.",
            "wide": True,
        },
    )
    terrain_wvp_path = output / "assets/pair-examples/terrain-wvp-diagnostic.json"
    if not terrain_wvp_path.exists():
        raise ValueError(
            "missing terrain/WVP diagnostic; run "
            "build_l2a_l1c_terrain_wvp_diagnostic.py first"
        )
    terrain_wvp_diagnostic = _json(terrain_wvp_path)
    surface_figures.insert(
        1,
        {
            "path": terrain_wvp_diagnostic["figure"],
            "title": "Terrain and water-vapour attribution diagnostic",
            "caption": "For the same four exact image pairs, DEM elevation, slope, local solar incidence, L2A-minus-scalar WVP, the visible L2A-current-RT residual, and the WVP-only counterfactual target shift are shown side by side. The final column changes only the L1C correction WVP; it tests sensitivity rather than providing a production correction. Correlations are within-scene associations, not causal attribution.",
            "wide": True,
        },
    )
    history_case_audits = _history_case_audits()
    missing_history_audits = sorted(set(matchup_ids) - set(history_case_audits))
    if missing_history_audits:
        raise ValueError(
            "missing clean-day history audits for " + ", ".join(missing_history_audits)
        )
    terrain_history_case_audits = (
        _terrain_history_case_audits(
            Path(terrain_artifacts["history_audits"]),
            str(terrain_artifacts["history_tag"]),
        )
        if terrain_evidence_available
        else {}
    )

    cases: list[dict[str, Any]] = []
    for matchup_id in matchup_ids:
        old = base_cases[matchup_id]
        truth = float(old["truth"])
        candidate_list = [
            _candidate(
                next(item["label"] for item in variants if item["id"] == variant_id),
                variant_records[variant_id].get(matchup_id, {}).get("retrieved"),
                truth,
                variant_id,
            )
            for variant_id in (
                "fresh",
                "current_best",
                "identity",
                "legacy_offset",
                "all_aod_live",
                "all_aod",
                "all_aod_mapunc",
                "all_aod_product_unc",
                "all_aod_tau_gated",
                "all_aod_tau_always",
                "terrain_a100_cap015",
                "terrain_a100_cap030",
                "terrain_blue_cap015",
                "terrain_blue_cap030",
                "terrain_clean_a100_cap030",
                "terrain_clean_blue_cap030",
                "paired_exact",
                "scene_mapping",
                "scene_mapping_mapunc",
                "scene_solver",
                "harmonized",
                "clean_visible",
                "clean_solver",
                "clean_domain",
            )
        ]
        candidate_by_id = {row["source"]: row for row in candidate_list}
        fresh = reference[matchup_id]
        clean = harmonized_records[matchup_id]
        fresh_hit = bool(candidate_by_id["fresh"]["within_ee"])
        clean_hit = bool(candidate_by_id["harmonized"]["within_ee"])
        transition = (
            "gain"
            if clean_hit and not fresh_hit
            else "loss"
            if fresh_hit and not clean_hit
            else "stable_hit"
            if fresh_hit
            else "stable_miss"
        )
        cube_path = CLEAN_CUBES / f"{matchup_id}.npz"
        if not cube_path.exists():
            raise ValueError(f"missing harmonized cost cube: {cube_path}")
        with np.load(cube_path, allow_pickle=False) as cost:
            local_mask, iy, ix = _window_mask(
                np.asarray(cost["x"], dtype=np.float64),
                np.asarray(cost["y"], dtype=np.float64),
                lon=float(clean["lon"]),
                lat=float(clean["lat"]),
                crs=str(clean["scene_crs"]),
            )
            spatial_path = output / "assets/harmonized-spatial" / f"{matchup_id}.png"
            diagnostic_path = output / "assets/harmonized-diagnostic" / f"{matchup_id}.png"
            if not (reuse_figures and spatial_path.exists()):
                solved, solved_unc, _pooled_cost = _pooled_maps(cost)
                no_backstop, _unused_unc, _unused_cost = _pooled_maps(
                    cost,
                    prior_unc=np.full(np.asarray(cost["solve_valid"]).shape, np.inf),
                )
                _save_spatial_figure(
                    spatial_path,
                    cost=cost,
                    record=clean,
                    pooled_aod=solved,
                    pooled_unc=solved_unc,
                    no_backstop_aod=no_backstop,
                    iy=iy,
                    ix=ix,
                )
            if not (reuse_figures and diagnostic_path.exists()):
                _candidate_cost_figure(
                    diagnostic_path,
                    cost=cost,
                    record=clean,
                    local_mask=local_mask,
                )
        terrain_record = terrain_records.get(matchup_id)
        terrain_spatial_image: str | None = None
        terrain_diagnostic_image: str | None = None
        if terrain_evidence_available and terrain_record is not None:
            terrain_cube_path = terrain_cubes / f"{matchup_id}.npz"
            with np.load(terrain_cube_path, allow_pickle=False) as terrain_cost:
                terrain_mask, terrain_iy, terrain_ix = _window_mask(
                    np.asarray(terrain_cost["x"], dtype=np.float64),
                    np.asarray(terrain_cost["y"], dtype=np.float64),
                    lon=float(terrain_record["lon"]),
                    lat=float(terrain_record["lat"]),
                    crs=str(terrain_record["scene_crs"]),
                )
                terrain_spatial_path = (
                    output / "assets/terrain-spatial" / f"{matchup_id}.png"
                )
                terrain_diagnostic_path = (
                    output / "assets/terrain-diagnostic" / f"{matchup_id}.png"
                )
                if not (reuse_figures and terrain_spatial_path.exists()):
                    terrain_solved, terrain_unc, _terrain_cost = _pooled_maps(
                        terrain_cost
                    )
                    terrain_no_backstop, _unused_unc, _unused_cost = _pooled_maps(
                        terrain_cost,
                        prior_unc=np.full(
                            np.asarray(terrain_cost["solve_valid"]).shape, np.inf
                        ),
                    )
                    _save_spatial_figure(
                        terrain_spatial_path,
                        cost=terrain_cost,
                        record=terrain_record,
                        pooled_aod=terrain_solved,
                        pooled_unc=terrain_unc,
                        no_backstop_aod=terrain_no_backstop,
                        iy=terrain_iy,
                        ix=terrain_ix,
                    )
                if not (reuse_figures and terrain_diagnostic_path.exists()):
                    _candidate_cost_figure(
                        terrain_diagnostic_path,
                        cost=terrain_cost,
                        record=terrain_record,
                        local_mask=terrain_mask,
                        title_prefix=terrain_best["label"],
                    )
            terrain_spatial_image = f"assets/terrain-spatial/{matchup_id}.png"
            terrain_diagnostic_image = f"assets/terrain-diagnostic/{matchup_id}.png"
        solver = fresh.get("solver") or {}
        clean_solver = clean.get("solver") or {}
        terrain_solver = (
            terrain_record.get("solver") or {} if terrain_record else {}
        )
        bands = [
            _finite(solver.get(f"surface_band_{band}_argmin_aot")) for band in BANDS
        ]
        finite_bands = [value for value in bands if value is not None]
        cases.append(
            {
                **{
                    key: old.get(key)
                    for key in (
                        "case_index",
                        "matchup_id",
                        "site",
                        "truth",
                        "ee",
                        "ee_low",
                        "ee_high",
                        "lon",
                        "lat",
                        "cloud_fraction",
                        "scene_cloud_cover",
                        "aeronet_count",
                        "aeronet_std",
                        "angstrom",
                        "elevation_m",
                        "maiac",
                        "maiac_unc",
                        "surface_anchor_aod",
                        "surface_min",
                        "cost_per_band",
                        "curvature",
                        "valid_fraction",
                    )
                },
                "fresh": float(fresh["retrieved"]),
                "fresh_error": float(fresh["retrieved"]) - truth,
                "fresh_error_over_ee": (float(fresh["retrieved"]) - truth)
                / (0.05 + 0.15 * truth),
                "direction": "under" if float(fresh["retrieved"]) < truth else "over",
                "transition": transition,
                "band_min_spread": max(finite_bands) - min(finite_bands)
                if finite_bands
                else None,
                "candidate_list": candidate_list,
                "candidate_by_id": candidate_by_id,
                "solver": solver,
                "harmonized_solver": clean_solver,
                "terrain_solver": terrain_solver,
                "terrain_evidence_id": terrain_best["id"]
                if terrain_evidence_available
                else None,
                "terrain_evidence_label": terrain_best["label"]
                if terrain_evidence_available
                else None,
                "atmo_prior": fresh.get("atmo_prior") or {},
                "prior_boa": fresh.get("prior_boa") or {},
                "prior_unc": fresh.get("prior_unc") or {},
                "harmonized_prior_boa": clean.get("prior_boa") or {},
                "harmonized_prior_unc": clean.get("prior_unc") or {},
                "terrain_prior_boa": terrain_record.get("prior_boa") or {}
                if terrain_record
                else {},
                "terrain_prior_unc": terrain_record.get("prior_unc") or {}
                if terrain_record
                else {},
                "harmonization_history": history_case_audits[matchup_id],
                "terrain_history": terrain_history_case_audits.get(matchup_id),
                "anchor_iterate": fresh.get("anchor_iterate") or {},
                "diagnostic": old.get("diagnostic") or {},
                "spatial_image": f"../{BASE_REPORT.name}/assets/spatial/{matchup_id}.png",
                "diagnostic_image": f"../{BASE_REPORT.name}/assets/diagnostic/{matchup_id}.png",
                "surface_image": f"../{BASE_REPORT.name}/assets/surface/{matchup_id}.png",
                "harmonized_spatial_image": f"assets/harmonized-spatial/{matchup_id}.png",
                "harmonized_diagnostic_image": f"assets/harmonized-diagnostic/{matchup_id}.png",
                "terrain_spatial_image": terrain_spatial_image,
                "terrain_diagnostic_image": terrain_diagnostic_image,
                "result_url": f"../../{BASELINE_DIR.name}/{matchup_id}.json",
                "cost_url": f"../../{BASELINE_CUBES.name}/{matchup_id}.npz",
                "harmonized_result_url": f"../../{CLEAN_RESULTS.name}/{matchup_id}.json",
                "harmonized_cost_url": f"../../{CLEAN_CUBES.name}/{matchup_id}.npz",
                "terrain_result_url": f"../../{terrain_results.name}/{matchup_id}.json"
                if terrain_record
                else None,
                "terrain_cost_url": f"../../{terrain_cubes.name}/{matchup_id}.npz"
                if terrain_evidence_available
                else None,
            }
        )

    _retrieval_figures(output, cases, variants)
    terrain_retrieval_ablation = (
        _terrain_retrieval_figure(
            output,
            cases,
            candidate_id=terrain_best["id"],
            candidate_label=terrain_best["label"],
        )
        if terrain_evidence_available
        else None
    )
    if terrain_retrieval_ablation is not None:
        surface_figures.append(
            {
                "path": terrain_retrieval_ablation["figure"],
                "title": "Terrain-conditioned AOD gains and losses",
                "caption": terrain_best["label"]
                + " is shown for every replay case. Axes retain the individual slope, local illumination, surface-correction, and AOD outcomes rather than inferring a single cause from an aggregate score.",
                "wide": True,
            }
        )
    surface_figures.append(
        {
            "path": "assets/surface-to-aod-response.png",
            "title": "Surface-to-AOD propagation",
            "caption": "Per-case prior changes, historical correction magnitude, training-domain exposure, and the resulting AOD change; colors preserve gains, losses, and stable outcomes without assigning a cause.",
            "wide": True,
        }
    )
    implementation = _copy_implementation(output)
    inventory = _inventory(matchup_ids, reference, variants)
    pair_audit = _pair_audit()
    application_audit = _history_application_audit(
        CLEAN_MODEL_DIR / "daily_histories/audits"
    )
    domain_guard_audit = _history_application_audit(
        CLEAN_MODEL_DIR / "daily_histories_domain015/audits"
    )
    clean_metrics = _json(CLEAN_MODEL_DIR / "surface_metrics.json")
    clean_model = _json(CLEAN_MODEL_DIR / "harmonizer.json")
    all_sites = {mid.split("__", 1)[0] for mid in clean_model["fold_by_matchup_id"]}
    clean_sites = {row["site"] for row in site_rows}
    missing_sites = sorted(all_sites - clean_sites)
    best = max(
        (row for row in variants if row["complete"]),
        key=lambda row: row["metrics"]["within_ee_rate"],
    )
    harmonized = next(row for row in variants if row["id"] == "harmonized")
    report = {
        "schema_version": 1,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "title": "L2A surface-prior harmonization investigation",
        "subtitle": "Operational Sentinel-2 L2A mapped to MAIAC-conditioned current SIAC RT",
        "comparison": {
            "id": "harmonized",
            "label": "Clean-day harmonized",
            "short_label": "harmonized",
            "reference_id": "fresh",
            "reference_label": "Current SIAC",
            "reference_short_label": "current",
        },
        "target": {
            "within_ee_rate": TARGET_RATE,
            "development_required_hits": TARGET_HITS,
            "development_count": 36,
            "best_variant": best["label"],
            "best_hits": best["metrics"]["within_ee"],
            "best_rate": best["metrics"]["within_ee_rate"],
            "gap_cases": TARGET_HITS - best["metrics"]["within_ee"],
            "met": best["metrics"]["within_ee"] >= TARGET_HITS,
            "harmonized_hits": harmonized["metrics"]["within_ee"],
            "harmonized_rate": harmonized["metrics"]["within_ee_rate"],
            "terrain_cap030_hits": terrain_cap030["metrics"]["within_ee"],
            "terrain_cap030_rate": terrain_cap030["metrics"]["within_ee_rate"],
            "terrain_best_variant": terrain_best["label"],
            "terrain_best_hits": terrain_best["metrics"]["within_ee"],
            "terrain_best_rate": terrain_best["metrics"]["within_ee_rate"],
        },
        "cohort": base["cohort"],
        "overview_captions": {
            "retrieval_scatter": (
                "Current SIAC, the current prior-conflict candidate, the all-AOD mapping, and "
                + terrain_best["label"]
                + " on the same 36 locked cases."
                if terrain_complete
                else "Current SIAC, the current prior-conflict candidate, the all-AOD mapping, and the clean-day mapping on the same 36 locked cases."
            ),
            "evidence_matrix": "Every case and principal harmonization arm, expressed in signed AERONET expected-error units.",
            "spatial_error": (
                "Geographic and per-case change in signed AOD error from current SIAC to "
                + terrain_best["label"]
                + "."
                if terrain_complete
                else "Geographic and per-case change in signed AOD error from current SIAC to clean-day harmonization."
            ),
        },
        "surface": {
            "pair_audit": pair_audit,
            "application_audit": application_audit,
            "domain_guard_audit": domain_guard_audit,
            "clean_training": {
                "sample_count": clean_metrics["sample_count"],
                "scene_count": clean_metrics["scene_count"],
                "site_count": clean_metrics["site_count"],
                "maiac_aot_max": clean_metrics["training_maiac_aot_max"],
                "fold_count": len(set(clean_model["fold_by_matchup_id"].values())),
                "missing_site_count": len(missing_sites),
                "missing_sites": missing_sites,
            },
            "target_rt": clean_model["target_rt"],
            "model": {
                "target": clean_model["target"],
                "training_cutoff": clean_model["training_cutoff"],
                "uses_aeronet": clean_model["uses_aeronet"],
                "selected_specification": "full_a100",
                "feature_set": clean_model["models"]["full_a100"]["feature_set"],
                "ridge_alpha": clean_model["models"]["full_a100"]["alpha"],
                "correction_cap": 0.03,
                "application_bands": list(DISPLAY_BANDS),
                "feature_schema_version": clean_model["feature_schema_version"],
            },
            "band_metrics": metric_rows,
            "pair_scatter_metrics": pair_scatter_metrics,
            "pair_examples": pair_examples,
            "terrain_wvp_diagnostic": terrain_wvp_diagnostic,
            "terrain_ablation": terrain_surface_ablation,
            "terrain_retrieval_ablation": terrain_retrieval_ablation,
            "site_rows": site_rows,
            "figures": surface_figures,
        },
        "variants": variants,
        "experiment_inventory": inventory,
        "cases": cases,
        "jobs": _job_audit(),
        "audit_counts": {
            "case_figures": (7 if terrain_evidence_available else 5) * len(cases),
            "cost_cubes": len(list(CLEAN_CUBES.glob("*.npz"))),
            "terrain_evidence_variant": terrain_best["id"]
            if terrain_evidence_available
            else None,
            "terrain_cost_cubes": len(list(terrain_cubes.glob("*.npz"))),
        },
        "limitations": [
            "The same-day MAIAC-conditioned L1C target is a consistency target, not independent surface truth.",
            "The saved historical RT sidecar has AOD, water-vapour and geometry support but no per-pixel elevation axis; terrain is therefore a held-site surface-residual feature, not an approximated elevation RT correction.",
            "The exact paired-history diagnostic did not improve AOD retrieval, so L2A RT inconsistency is not the dominant remaining error in this cohort.",
            "Four development sites have no historical MAIAC AOD <=0.15 pair and are genuine unseen-site applications of the fixed folds.",
            "Surface-reflectance MAE and AERONET AOD within-EE measure different stages and must not be substituted for one another.",
            "The locked holdout was not opened because no candidate reached 32/36 on development.",
        ],
        "reproduction": {
            "constraints": [
                "One operational Sentinel-2 L2A surface-prior source type",
                "No AERONET values or labels in the harmonizer",
                "No case-level source selection",
                "One fixed tile-wide continental-average libRadTran aerosol model",
                "Held-site folds for development; full model only for production",
                "Holdout remains locked below 32/36 development hits",
            ],
            "pipeline": [
                {"stage": "Inputs", "detail": "Sentinel-2 L1C target scene; operational COPERNICUS/S2_SR_HARMONIZED history; per-acquisition L2A AOT/WVP; matched MAIAC AOD/TCWV; sun/view/relative azimuth; elevation; month; sensor; processing baseline."},
                {"stage": "Pair target", "detail": "For each historical acquisition, cloud-score and SCL clear-land pixels are paired exactly with the same-day L1C corrected using saved MAIAC-conditioned libRadTran coefficients in the current SIAC RT space."},
                {"stage": "Clean gate", "detail": "Training uses dates through 2023-12-31 with MAIAC AOD550 <=0.15. The gate limits aerosol contamination in the target while preserving all available sites with clean scenes."},
                {"stage": "Features", "detail": "Per band: all seven L2A reflectances; MAIAC, Sen2Cor and difference AOD/WV; two-way airmass; cos(relative azimuth); elevation; month sin/cos; S2B/S2C flags; processing baseline; band x AOD/WV/airmass interactions; delta-AOD x airmass."},
                {"stage": "Terrain ablation", "detail": "The terrain-conditioned arm adds GLO-30 local elevation, slope, and local solar-incidence cosine from the acquisition sun azimuth and zenith. It uses the same operational L2A source, exact-pair target, site folds, and ridge solver as the non-terrain arm."},
                {"stage": "Training", "detail": "Standardize each feature and fit one ridge residual model per band with alpha=100. Five fixed site-group folds produce leakage-safe development coefficients; each exported band model stores mean, scale, coefficient vector and intercept in JSON."},
                {"stage": "Per-acquisition mapping", "detail": "Predict current-RT minus L2A reflectance independently for B01/B02/B03/B04/B8A/B11/B12; clip each correction to +/-0.03 and output reflectance to [0.001,0.8]. Apply before same-day tile mosaicking or temporal compositing."},
                {"stage": "Composite", "detail": "Median-mosaic same-day tiles, median across days within each historical year-month realization, and retain at least two realizations. Development uses the held-site fold; production uses the full pre-2024 model."},
                {"stage": "Surface prior", "detail": "Fit the same ET20 seasonal predictor to each harmonized historical realization. Predict B02/B03/B04 from B8A/B11/B12 and historical localizers; aggregate realization medians and robust spread with a 0.006 visible uncertainty floor."},
                {"stage": "Target radiometry and masks", "detail": "Apply L1C quantification and radiometric offsets. OmniCloudMask, water and no-data exclusions are resampled to the 60 m solve grid; cloud bypass is disabled except explicit fallback after no valid observation."},
                {"stage": "Atmospheric prior", "detail": "Use staged MCD19A2 MAIAC AOD550 as one tile-wide centre with the fixed calibrated uncertainty. TCWV is fixed at 2.0 cm and elevation at sea level in this reproduction branch."},
                {"stage": "Radiative transfer", "detail": "Use the packaged libRadTran continental-average 1 nm Zarr LUT, Sentinel-2 spectral convolution, scene-mean geometry, and the fixed 68-node AOD axis from 0.01 to 4.0."},
                {"stage": "Solve and output", "detail": "Sum diagonal B02/B03/B04 chi-square surface costs, median-pool each AOD-node cost in a 20x20-pixel window, add the Gaussian MAIAC penalty, take per-pixel argmin, and report the mean finite solved AOD over the AOI."},
            ],
            "formula": "rho_target,b = rho_L2A,b + clip(beta_0,b + sum_j beta_j,b * (x_j - mu_j) / s_j, -0.03, 0.03)\n\nJ_surface(tau,p) = sum_b [(rho_BOA,b(tau,p) - rho_prior,b(p)) / sigma_b(p)]^2\nJ_total(tau,p) = median_window(J_surface) + [(tau - tau_MAIAC) / sigma_MAIAC]^2\ntau_hat(p) = argmin_tau J_total(tau,p)\nAOD_scene = mean of finite tau_hat over the retrieval AOI\n\nEE(A) = 0.05 + 0.15 A; within_EE = |AOD_scene - AERONET| <= EE(AERONET)",
            "calibration_title": "Exported residual mapping",
            "calibration_description": "Each band is an explicit standardized linear residual model. Development loads the held-site fold; operational deployment loads the full coefficient set. AERONET is absent from both.",
            "calibration_formula": "z_j = (x_j - mean_j) / scale_j\ndelta_rho_b = clip(intercept_b + sum_j coefficient_bj * z_j, -0.03, 0.03)\nrho_harmonized,b = clip(rho_L2A,b + delta_rho_b, 0.001, 0.8)",
            "parameters": {
                "training_cutoff": "2023-12-31",
                "clean_maiac_aod_max": 0.15,
                "cross_validation": "five fixed site-group folds",
                "model": "per-band Ridge",
                "ridge_alpha": 100.0,
                "correction_cap_reflectance": 0.03,
                "surface_source": "COPERNICUS/S2_SR_HARMONIZED only",
                "application_order": "per acquisition, then day mosaic, then temporal median",
                "surface_predictor": "ET20 seasonal ExtraTrees",
                "solve_bands": ["B02", "B03", "B04"],
                "surface_uncertainty_floor": 0.006,
                "rt_lut": str(clean_model["target_rt"]["lut"]),
                "aerosol_profile": "continental_average",
                "aod_axis_nodes": 68,
                "solver_resolution_m": 60.0,
                "pool_window_pixels": 20,
                "uses_aeronet": False,
            },
            "acceptance": "Freeze only at 32/36 or better with 36/36 terminal outputs; otherwise retain all evidence, do not score the holdout, and do not describe the mapping as a production improvement.",
            "conflict_rule": "Not used by the harmonization candidate; retained only as a current-code comparison.",
            "implementation_files": implementation,
        },
        "sources": {
            "cases_csv": "downloads/cases.csv",
            "variants_csv": "downloads/variants.csv",
            "surface_band_metrics": "downloads/surface-band-metrics.csv",
            "pair_scatter_metrics": "downloads/exact-pair-scatter-metrics.csv",
            "pair_image_gallery": "assets/pair-examples/pair-gallery-metadata.json",
            "terrain_wvp_diagnostic": "assets/pair-examples/terrain-wvp-diagnostic.json",
            "terrain_wvp_metrics": "downloads/terrain-wvp-diagnostic.csv",
            "terrain_surface_band_metrics": "downloads/terrain-surface-band-metrics.csv",
            "terrain_surface_scene_metrics": "downloads/terrain-surface-scenes.csv",
            "terrain_surface_site_metrics": "downloads/terrain-surface-sites.csv",
            "terrain_retrieval_cases": "downloads/terrain-retrieval-cases.csv",
            "surface_scene_metrics": "downloads/surface-scenes.csv",
            "surface_site_metrics": "downloads/surface-sites.csv",
            "pair_audits": f"../../analysis/{PAIR_DIR.name}/",
            "terrain_pair_audits": f"../../analysis/{TERRAIN_PAIR_DIR.name}/",
            "all_aod_model": f"../../analysis/{ALL_MODEL_DIR.name}/harmonizer.json",
            "clean_model": f"../../analysis/{CLEAN_MODEL_DIR.name}/harmonizer.json",
            "terrain_all_aod_model": f"../../analysis/{TERRAIN_ALL_MODEL_DIR.name}/harmonizer.json",
            "terrain_clean_model": f"../../analysis/{TERRAIN_CLEAN_MODEL_DIR.name}/harmonizer.json",
            "terrain_results_cap015": f"../../{TERRAIN_RESULTS_015.name}/",
            "terrain_results_cap030": f"../../{TERRAIN_RESULTS_030.name}/",
            "terrain_cost_cubes_cap015": f"../../{TERRAIN_CUBES_015.name}/",
            "terrain_cost_cubes_cap030": f"../../{TERRAIN_CUBES_030.name}/",
            "terrain_clean_results_cap030": f"../../{TERRAIN_CLEAN_RESULTS_030.name}/",
            "terrain_clean_cost_cubes_cap030": f"../../{TERRAIN_CLEAN_CUBES_030.name}/",
            "terrain_clean_blue_results_cap030": f"../../{TERRAIN_CLEAN_BLUE_RESULTS_030.name}/",
            "terrain_clean_blue_cost_cubes_cap030": f"../../{TERRAIN_CLEAN_BLUE_CUBES_030.name}/",
            "harmonized_results": f"../../{CLEAN_RESULTS.name}/",
            "harmonized_cost_cubes": f"../../{CLEAN_CUBES.name}/",
            "prior_full_investigation": f"../{BASE_REPORT.name}/",
            "report_json": "data/report.json",
            "implementation_manifest": "downloads/implementation/manifest.json",
        },
    }

    (output / "data/report.json").write_text(
        json.dumps(report, separators=(",", ":")) + "\n", encoding="utf-8"
    )
    _write_csv(output / "downloads/surface-band-metrics.csv", metric_rows)
    _write_csv(
        output / "downloads/exact-pair-scatter-metrics.csv", pair_scatter_metrics
    )
    shutil.copy2(
        output / "assets/pair-examples/terrain-wvp-diagnostic.csv",
        output / "downloads/terrain-wvp-diagnostic.csv",
    )
    _write_csv(output / "downloads/surface-scenes.csv", clean_scene_rows)
    _write_csv(output / "downloads/surface-sites.csv", site_rows)
    if terrain_surface_ablation is not None:
        _write_csv(
            output / "downloads/terrain-surface-band-metrics.csv",
            terrain_surface_ablation["band_rows"],
        )
        _write_csv(
            output / "downloads/terrain-surface-scenes.csv",
            terrain_surface_ablation["scene_rows"],
        )
        _write_csv(
            output / "downloads/terrain-surface-sites.csv",
            terrain_surface_ablation["site_rows"],
        )
    if terrain_retrieval_ablation is not None:
        _write_csv(
            output / "downloads/terrain-retrieval-cases.csv",
            terrain_retrieval_ablation["rows"],
        )
    _write_csv(
        output / "downloads/variants.csv",
        [
            {
                "id": row["id"],
                "label": row["label"],
                "family": row["family"],
                "ok": row["metrics"]["n"],
                "hits": row["metrics"]["within_ee"],
                "within_ee_percent": row["metrics"]["within_ee_percent"],
                "bias": row["metrics"]["bias"],
                "mae": row["metrics"]["mae"],
                "rmse": row["metrics"]["rmse"],
                "gains": row["transitions_vs_fresh"].get("gain", 0),
                "losses": row["transitions_vs_fresh"].get("loss", 0),
                "description": row["description"],
                "path": row["path"],
            }
            for row in variants
        ],
    )
    _write_csv(
        output / "downloads/experiment-inventory.csv",
        [
            {
                "id": row["id"],
                "label": row["label"],
                "family": row["family"],
                "inventory_kind": row["inventory_kind"],
                "ok": row["metrics"]["n"],
                "hits": row["metrics"]["within_ee"],
                "within_ee_percent": row["metrics"]["within_ee_percent"],
                "bias": row["metrics"]["bias"],
                "mae": row["metrics"]["mae"],
                "rmse": row["metrics"]["rmse"],
                "gains": row["transitions_vs_fresh"].get("gain", 0),
                "losses": row["transitions_vs_fresh"].get("loss", 0),
                "status_counts": json.dumps(row["status_counts"], sort_keys=True),
                "description": row["description"],
                "path": row["path"],
            }
            for row in inventory
        ],
    )
    _write_csv(
        output / "downloads/cases.csv",
        [
            {
                "matchup_id": row["matchup_id"],
                "site": row["site"],
                "truth": row["truth"],
                "current": row["candidate_by_id"]["fresh"]["value"],
                "current_within_ee": row["candidate_by_id"]["fresh"]["within_ee"],
                "terrain_cap015": row["candidate_by_id"]["terrain_a100_cap015"]["value"],
                "terrain_cap015_within_ee": row["candidate_by_id"]["terrain_a100_cap015"]["within_ee"],
                "terrain_cap030": row["candidate_by_id"]["terrain_a100_cap030"]["value"],
                "terrain_cap030_within_ee": row["candidate_by_id"]["terrain_a100_cap030"]["within_ee"],
                "terrain_cap030_error_change": row["candidate_by_id"]["terrain_a100_cap030"]["error"]
                - row["candidate_by_id"]["fresh"]["error"]
                if row["candidate_by_id"]["terrain_a100_cap030"]["error"] is not None
                else None,
                "terrain_blue_cap015": row["candidate_by_id"]["terrain_blue_cap015"]["value"],
                "terrain_blue_cap015_within_ee": row["candidate_by_id"]["terrain_blue_cap015"]["within_ee"],
                "terrain_blue_cap030": row["candidate_by_id"]["terrain_blue_cap030"]["value"],
                "terrain_blue_cap030_within_ee": row["candidate_by_id"]["terrain_blue_cap030"]["within_ee"],
                "terrain_blue_cap030_error_change": row["candidate_by_id"]["terrain_blue_cap030"]["error"]
                - row["candidate_by_id"]["fresh"]["error"]
                if row["candidate_by_id"]["terrain_blue_cap030"]["error"] is not None
                else None,
                "terrain_clean_cap030": row["candidate_by_id"]["terrain_clean_a100_cap030"]["value"],
                "terrain_clean_cap030_within_ee": row["candidate_by_id"]["terrain_clean_a100_cap030"]["within_ee"],
                "terrain_clean_cap030_error_change": row["candidate_by_id"]["terrain_clean_a100_cap030"]["error"]
                - row["candidate_by_id"]["fresh"]["error"]
                if row["candidate_by_id"]["terrain_clean_a100_cap030"]["error"] is not None
                else None,
                "terrain_clean_blue_cap030": row["candidate_by_id"]["terrain_clean_blue_cap030"]["value"],
                "terrain_clean_blue_cap030_within_ee": row["candidate_by_id"]["terrain_clean_blue_cap030"]["within_ee"],
                "terrain_clean_blue_cap030_error_change": row["candidate_by_id"]["terrain_clean_blue_cap030"]["error"]
                - row["candidate_by_id"]["fresh"]["error"]
                if row["candidate_by_id"]["terrain_clean_blue_cap030"]["error"] is not None
                else None,
                "harmonized": row["candidate_by_id"]["harmonized"]["value"],
                "harmonized_within_ee": row["candidate_by_id"]["harmonized"]["within_ee"],
                "harmonized_error_change": row["candidate_by_id"]["harmonized"]["error"]
                - row["candidate_by_id"]["fresh"]["error"],
                "harmonized_prior_delta_b02": row["harmonized_prior_boa"].get("B02")
                - row["prior_boa"].get("B02"),
                "harmonized_prior_delta_b03": row["harmonized_prior_boa"].get("B03")
                - row["prior_boa"].get("B03"),
                "harmonized_prior_delta_b04": row["harmonized_prior_boa"].get("B04")
                - row["prior_boa"].get("B04"),
                "history_scene_count": row["harmonization_history"]["scene_count"],
                "history_outside_training_aod_fraction": row["harmonization_history"][
                    "outside_training_aod_fraction"
                ],
                "history_maiac_aot_median": row["harmonization_history"][
                    "maiac_aot_median"
                ],
                "history_b8a_median_correction": row["harmonization_history"][
                    "per_band_median_scene_statistics"
                ]["nir08"]["median"],
                "clean_visible": row["candidate_by_id"]["clean_visible"]["value"],
                "clean_visible_within_ee": row["candidate_by_id"]["clean_visible"][
                    "within_ee"
                ],
                "clean_solver": row["candidate_by_id"]["clean_solver"]["value"],
                "clean_solver_within_ee": row["candidate_by_id"]["clean_solver"][
                    "within_ee"
                ],
                "clean_domain": row["candidate_by_id"]["clean_domain"]["value"],
                "clean_domain_within_ee": row["candidate_by_id"]["clean_domain"][
                    "within_ee"
                ],
                "transition": row["transition"],
                "cloud_fraction": row["cloud_fraction"],
            }
            for row in cases
        ],
    )
    (output / "index.html").write_text(
        """<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width,initial-scale=1">
  <meta name="color-scheme" content="light">
  <title>L2A surface-prior harmonization investigation</title>
  <link rel="icon" href="data:,">
  <link rel="stylesheet" href="app.css">
</head>
<body>
  <div id="app"><div class="loading">Loading harmonization evidence...</div></div>
  <noscript>This report requires JavaScript. CSV and JSON downloads remain available.</noscript>
  <script src="app.js"></script>
</body>
</html>
""",
        encoding="utf-8",
    )
    receipt = {
        "output": str(output),
        "cases": len(cases),
        "variants": len(variants),
        "harmonized_hits": harmonized["metrics"]["within_ee"],
        "best_hits": best["metrics"]["within_ee"],
        "target_met": report["target"]["met"],
        "holdout_scored": False,
        "cost_cubes": len(list(CLEAN_CUBES.glob("*.npz"))),
        "harmonized_spatial_images": len(list((output / "assets/harmonized-spatial").glob("*.png"))),
        "harmonized_diagnostic_images": len(list((output / "assets/harmonized-diagnostic").glob("*.png"))),
        "terrain_spatial_images": len(list((output / "assets/terrain-spatial").glob("*.png"))),
        "terrain_diagnostic_images": len(list((output / "assets/terrain-diagnostic").glob("*.png"))),
    }
    (output / "build-receipt.json").write_text(
        json.dumps(receipt, indent=2) + "\n", encoding="utf-8"
    )
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--reuse-figures", action="store_true")
    args = parser.parse_args()
    print(json.dumps(build(args.output, reuse_figures=args.reuse_figures), indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
