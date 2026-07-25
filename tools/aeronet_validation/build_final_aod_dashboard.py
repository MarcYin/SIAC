"""Build the final low-cloud AOD performance and evidence dashboard."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import shutil
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np

if TYPE_CHECKING:
    from collections.abc import Iterable, Sequence

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from tools.aeronet_validation.aod_residual_calibration import (  # noqa: E402
    CalibrationSample,
    load_samples,
    metrics,
)
from tools.aeronet_validation.build_low_cloud_failure_explorer import (  # noqa: E402
    _jsonable,
    _pooled_maps,
    _save_diagnostic_figure,
    _save_spatial_figure,
    _window_mask,
)
from tools.aeronet_validation.select_generic_aod_calibrator import _ids  # noqa: E402

ROOT = Path("/gws/ssde/j25a/nceo_isp/public/siac_refactor")
DEFAULT_OUTPUT = ROOT / "reports/aod-final-performance-dashboard-20260713"
WEB_ASSETS = Path(__file__).with_name("final_aod_dashboard")
FINAL_REPORT = ROOT / "reports/aod-global-offset-validation-20260713.json"
SELECTION_REPORT = ROOT / "reports/aod-target-domain-recipe-selection-20260713.json"
SEED_REPORT = ROOT / "reports/aod-selected-et-seed-robustness-20260713.json"
MODEL_MANIFEST = ROOT / "reports/aod-et35-global-offset-model-manifest-20260713.json"
ABLATION_REPORT = ROOT / "reports/aod-generic-feature-ablation-20260713.json"
LEGACY_EXPLORER = ROOT / "reports/aod-low-cloud-failure-explorer-20260712"
TARGET_RESULTS = ROOT / "phaseD_results_lowcloud20_native_maiac_adaptive_b03_chi2_20260713"
TARGET_CONTEXT = ROOT / "gee_aerosol_context_campaign250"
TARGET_ATMO = ROOT / "maiac_qa_lowcloud20_native_adaptive"
TARGET_MIDS = ROOT / "campaign250_lowcloud20_mids.txt"
PRIOR_RESULTS = ROOT / "prior_quality"

NEW_DIAGNOSTICS = (
    (
        ROOT / "phaseD_results_lowcloud20_finaloffset_lossdiag_b03_chi2_20260713",
        ROOT / "phaseD_cost_cubes_lowcloud20_finaloffset_lossdiag_b03_chi2_20260713",
    ),
    (
        ROOT / "phaseD_results_lowcloud20_finaloffset_gaindiag_b03_chi2_20260713",
        ROOT / "phaseD_cost_cubes_lowcloud20_finaloffset_gaindiag_b03_chi2_20260713",
    ),
)


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _load_records(directory: Path) -> dict[str, dict[str, Any]]:
    records: dict[str, dict[str, Any]] = {}
    for path in sorted(directory.glob("*.json")):
        record = _load_json(path)
        records[str(record.get("matchup_id") or path.stem)] = record
    return records


def _within_ee(value: float, truth: float) -> bool:
    return abs(value - truth) <= 0.05 + 0.15 * truth + 1e-12


def summarize_rows(rows: Sequence[dict[str, Any]], field: str) -> dict[str, Any]:
    """Summarize one prediction field for a bounded case subset."""
    truth = np.asarray([row["truth"] for row in rows], dtype=float)
    prediction = np.asarray([row[field] for row in rows], dtype=float)
    return metrics(truth, prediction)


def binned_summaries(
    rows: Sequence[dict[str, Any]],
    field: str,
    edges: Sequence[float],
    labels: Sequence[str],
) -> list[dict[str, Any]]:
    """Return reconciled baseline/model summaries for fixed numeric bins."""
    if len(labels) != len(edges) - 1:
        raise ValueError("Bin labels and edges do not agree")
    output: list[dict[str, Any]] = []
    for lower, upper, label in zip(edges[:-1], edges[1:], labels, strict=True):
        selected = [
            row
            for row in rows
            if row.get(field) is not None and lower <= float(row[field]) < upper
        ]
        if not selected:
            continue
        output.append(
            {
                "label": label,
                "lower": lower,
                "upper": upper if math.isfinite(upper) else None,
                "count": len(selected),
                "baseline": summarize_rows(selected, "baseline"),
                "unadjusted": summarize_rows(selected, "unadjusted"),
                "final": summarize_rows(selected, "adjusted"),
            }
        )
    return output


def _metric_block(rows: Sequence[dict[str, Any]]) -> dict[str, Any]:
    return {
        "count": len(rows),
        "baseline": summarize_rows(rows, "baseline"),
        "unadjusted": summarize_rows(rows, "unadjusted"),
        "final": summarize_rows(rows, "adjusted"),
    }


def _offset_curve(rows: Sequence[dict[str, Any]], offsets: Iterable[float]) -> list[dict[str, Any]]:
    output: list[dict[str, Any]] = []
    for offset in offsets:
        adjusted = []
        for row in rows:
            value = np.clip(
                (float(row["unadjusted"]) + 1 / 3) * np.exp(float(offset)) - 1 / 3,
                0.0,
                4.0,
            )
            adjusted.append({**row, "offset_prediction": float(value)})
        development = [row for row in adjusted if row["site_fold"] != 4]
        confirmation = [row for row in adjusted if row["site_fold"] == 4]
        output.append(
            {
                "offset": float(offset),
                "development": summarize_rows(development, "offset_prediction"),
                "confirmation": summarize_rows(confirmation, "offset_prediction"),
                "all": summarize_rows(adjusted, "offset_prediction"),
            }
        )
    return output


def _candidate(label: str, kind: str, value: float | None, truth: float) -> dict[str, Any]:
    return {
        "label": label,
        "kind": kind,
        "value": value,
        "within_ee": _within_ee(value, truth) if value is not None else None,
    }


def _new_diagnostic_source(matchup_id: str) -> tuple[Path, Path] | None:
    for result_dir, cube_dir in NEW_DIAGNOSTICS:
        if (result_dir / f"{matchup_id}.json").exists() and (
            cube_dir / f"{matchup_id}.npz"
        ).exists():
            return result_dir, cube_dir
    return None


def _build_new_diagnostic(
    output: Path,
    receipt: dict[str, Any],
    record: dict[str, Any],
    prior_record: dict[str, Any],
    result_dir: Path,
    cube_dir: Path,
) -> dict[str, Any]:
    matchup_id = receipt["matchup_id"]
    truth = float(receipt["truth"])
    cube_path = cube_dir / f"{matchup_id}.npz"
    diagnostic_record = _load_json(result_dir / f"{matchup_id}.json")
    with np.load(cube_path, allow_pickle=False) as cost:
        x = np.asarray(cost["x"], dtype=float)
        y = np.asarray(cost["y"], dtype=float)
        local_mask, iy, ix = _window_mask(
            x,
            y,
            lon=float(record["lon"]),
            lat=float(record["lat"]),
            crs=str(record["scene_crs"]),
        )
        solve_valid = np.asarray(cost["solve_valid"], dtype=bool)
        pooled_aod, pooled_unc, pooled_cost = _pooled_maps(cost)
        no_backstop_aod, _no_backstop_unc, _no_backstop_cost = _pooled_maps(
            cost,
            prior_unc=np.full(solve_valid.shape, np.inf, dtype=float),
        )
        no_backstop_values = no_backstop_aod[np.isfinite(no_backstop_aod)]
        no_backstop = (
            float(np.mean(no_backstop_values)) if no_backstop_values.size else None
        )
        solver_cams = float(np.nanmedian(np.asarray(cost["aot_prior"])[solve_valid]))
        solver_cams_unc = float(
            np.nanmedian(np.asarray(cost["aot_prior_unc"])[solve_valid])
        )
        candidates = [
            _candidate("Physical SIAC", "retrieval", receipt["baseline"], truth),
            _candidate("ExtraTrees raw", "calibrator", receipt["unadjusted"], truth),
            _candidate("Final global offset", "calibrator", receipt["adjusted"], truth),
            _candidate("Surface-anchor CAMS", "atmosphere", prior_record.get("cams_aot"), truth),
            _candidate("Solver MAIAC backstop", "atmosphere", solver_cams, truth),
            _candidate("No-backstop replay", "diagnostic", no_backstop, truth),
            _candidate(
                "Seed-300 adjusted minimum",
                "robustness",
                receipt["seed300_adjusted_min"],
                truth,
            ),
            _candidate(
                "Seed-300 adjusted maximum",
                "robustness",
                receipt["seed300_adjusted_max"],
                truth,
            ),
        ]
        spatial_path = output / "assets/spatial" / f"{matchup_id}.png"
        diagnostic_path = output / "assets/diagnostic" / f"{matchup_id}.png"
        _save_spatial_figure(
            spatial_path,
            cost=cost,
            record=record,
            pooled_aod=pooled_aod,
            pooled_unc=pooled_unc,
            no_backstop_aod=no_backstop_aod,
            iy=iy,
            ix=ix,
        )
        curves = _save_diagnostic_figure(
            diagnostic_path,
            cost=cost,
            record=record,
            prior_record=prior_record,
            local_mask=local_mask,
            pooled_aod=pooled_aod,
            pooled_unc=pooled_unc,
            candidates=candidates,
        )
        local_cost = pooled_cost[local_mask & solve_valid]
        cube = {
            "shape": list(np.asarray(cost["cube"]).shape),
            "bands": [str(value) for value in np.asarray(cost["band_names"])],
            "aod_nodes": int(np.asarray(cost["aot_axis"]).size),
            "aod_min": float(np.min(cost["aot_axis"])),
            "aod_max": float(np.max(cost["aot_axis"])),
            "pool_window": int(np.asarray(cost["pool_window"]).item()),
            "pool_min_count": int(np.asarray(cost["min_count"]).item()),
            "local_pooled_cost_median": (
                float(np.nanmedian(local_cost)) if np.isfinite(local_cost).any() else None
            ),
        }
    return {
        "source": "final exact diagnostic rerun",
        "spatial_image": f"assets/spatial/{matchup_id}.png",
        "diagnostic_image": f"assets/diagnostic/{matchup_id}.png",
        "canonical_result_url": f"../../{TARGET_RESULTS.name}/{matchup_id}.json",
        "diagnostic_result_url": f"../../{result_dir.name}/{matchup_id}.json",
        "cost_cube_url": f"../../{cube_dir.name}/{matchup_id}.npz",
        "diagnostic_retrieved": diagnostic_record.get("retrieved"),
        "surface_anchor_cams_aod": prior_record.get("cams_aot"),
        "solver_cams_aod": solver_cams,
        "solver_cams_unc": solver_cams_unc,
        "no_backstop_aod": no_backstop,
        "backstop_shift": (
            float(record["retrieved"]) - no_backstop if no_backstop is not None else None
        ),
        "candidates": candidates,
        "curves": curves,
        "cube": cube,
    }


def _copy_legacy_diagnostic(
    output: Path, old_root: Path, old: dict[str, Any]
) -> dict[str, Any]:
    matchup_id = old["matchup_id"]
    for kind, field in (("spatial", "spatial_image"), ("diagnostic", "diagnostic_image")):
        source = old_root / old[field]
        destination = output / "assets" / kind / f"{matchup_id}.png"
        shutil.copy2(source, destination)
    return {
        "source": "2026-07-12 physical cost diagnostic",
        "spatial_image": f"assets/spatial/{matchup_id}.png",
        "diagnostic_image": f"assets/diagnostic/{matchup_id}.png",
        "canonical_result_url": old.get("canonical_result_url"),
        "diagnostic_result_url": old.get("diagnostic_result_url"),
        "cost_cube_url": old.get("cost_cube_url"),
        "diagnostic_retrieved": old.get("diagnostic_retrieved"),
        "surface_anchor_cams_aod": old.get("surface_anchor_cams_aod"),
        "solver_cams_aod": old.get("solver_cams_aod"),
        "solver_cams_unc": old.get("solver_cams_unc"),
        "no_backstop_aod": old.get("no_backstop_aod"),
        "backstop_shift": old.get("backstop_shift"),
        "candidates": old.get("candidates", []),
        "curves": old.get("curves", {}),
        "cube": old.get("cube", {}),
        "legacy_group": {
            "code": old.get("mechanism_code"),
            "label": old.get("mechanism"),
            "evidence": old.get("diagnostic_evidence"),
        },
    }


def _experiment_ledger(root: Path) -> list[dict[str, Any]]:
    ledger = [
        {
            "category": "Physical baseline",
            "method": "Three-band chi-square, native MAIAC adaptive staging",
            "scope": "152 low-cloud cases",
            "hits": 111,
            "count": 152,
            "operational": True,
            "status": "reference",
            "source": "../aod-global-offset-validation-20260713.json",
            "note": "One S2 monthly SWIR/NIR-anchored ExtraTree surface prior; 152/152 terminal OK.",
        },
        {
            "category": "Surface / solver",
            "method": "Historical R2",
            "scope": "152 low-cloud cases",
            "hits": 111,
            "count": 152,
            "operational": True,
            "status": "tested",
            "source": "../aod-low-cloud-20260711/single-prior-scores.csv",
            "note": "Measured historical comparator.",
        },
        {
            "category": "Surface / solver",
            "method": "Fixed one-prior consensus",
            "scope": "152 low-cloud cases",
            "hits": 108,
            "count": 152,
            "operational": True,
            "status": "tested",
            "source": "../aod-low-cloud-20260711/single-prior-scores.csv",
            "note": "Validation-screened fixed posterior.",
        },
        {
            "category": "Band set",
            "method": "Add Sentinel-2 B01",
            "scope": "152 low-cloud cases",
            "hits": 98,
            "count": 152,
            "operational": True,
            "status": "rejected",
            "source": "../aod-b01-band-comparison-20260712/analysis.json",
            "note": "10 gains, 22 losses relative to the 110-hit run used in that report.",
        },
        {
            "category": "Aerosol RT",
            "method": "6S continental",
            "scope": "152 low-cloud cases",
            "hits": 110,
            "count": 152,
            "operational": True,
            "status": "tested",
            "source": "../aod-aerosol-species-20260712/analysis.json",
            "note": "10 gains and 10 losses versus LUT continental.",
        },
        {
            "category": "Aerosol RT",
            "method": "6S CCI-3, one tile-wide species",
            "scope": "152 low-cloud cases",
            "hits": 107,
            "count": 152,
            "operational": True,
            "status": "rejected",
            "source": "../aod-aerosol-species-20260712/analysis.json",
            "note": "Species fixed for the whole tile; selection based on aggregate cost.",
        },
        {
            "category": "Diagnostic ceiling",
            "method": "Expanded one-prior hindsight oracle",
            "scope": "152 low-cloud cases",
            "hits": 140,
            "count": 152,
            "operational": False,
            "status": "oracle",
            "source": "../aod-low-cloud-20260711/single-prior-scores.csv",
            "note": "Truth-selected ceiling, not an operational algorithm.",
        },
    ]

    calibration_sources = (
        (
            "External-selected HistGradientBoosting",
            "aod-generic-logratio-external588-20260713.json",
            "candidate",
            "externally selected",
        ),
        (
            "Operational-domain weighted scoring",
            "aod-generic-operational-domain-score-external588-20260713.json",
            "candidate",
            "rejected",
        ),
        (
            "Operational-domain weighted fit",
            "aod-generic-operational-domain-fit-external588-20260713.json",
            "candidate",
            "rejected",
        ),
        (
            "Externally selected convex blend",
            "aod-generic-convex-blend-external588-20260713.json",
            "candidate",
            "rejected",
        ),
    )
    for method, filename, metric_key, status in calibration_sources:
        record = _load_json(root / "reports" / filename)
        target = record["target"][metric_key]
        holdout = record["external_holdout"][metric_key]
        ledger.append(
            {
                "category": "Generic calibration",
                "method": method,
                "scope": "152 target / 123 external holdout",
                "hits": target["hits"],
                "count": target["count"],
                "holdout_hits": holdout["hits"],
                "holdout_count": holdout["count"],
                "operational": True,
                "status": status,
                "source": f"../{filename}",
                "note": "One uniform prediction rule; no case routing.",
            }
        )
    ledger.extend(
        (
            {
                "category": "Generic calibration",
                "method": "Target-development selected ExtraTrees, 1,500 trees",
                "scope": "152 target / 123 external holdout",
                "hits": 131,
                "count": 152,
                "holdout_hits": 82,
                "holdout_count": 123,
                "operational": True,
                "status": "pre-offset",
                "source": "../aod-selected-et-seed-robustness-20260713.json",
                "note": "Fixed absolute-error log-ratio model; no target label at inference.",
            },
            {
                "category": "Final generic calibration",
                "method": "ExtraTrees plus one global shifted-log offset",
                "scope": "152 target / 123 external holdout",
                "hits": 134,
                "count": 152,
                "holdout_hits": 83,
                "holdout_count": 123,
                "operational": True,
                "status": "final candidate",
                "source": "../aod-global-offset-validation-20260713.json",
                "note": "All tested full-schema seed/tree variants meet 133/152.",
            },
        )
    )
    return ledger


def _save_global_figures(output: Path, rows: Sequence[dict[str, Any]]) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    truth = np.asarray([row["truth"] for row in rows], dtype=float)
    baseline = np.asarray([row["baseline"] for row in rows], dtype=float)
    final = np.asarray([row["adjusted"] for row in rows], dtype=float)
    upper = max(1.65, float(np.max([truth.max(), baseline.max(), final.max()])) * 1.04)
    x = np.linspace(0, upper, 200)
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 6.0), facecolor="white")
    for ax, values, title in zip(
        axes, (baseline, final), ("Physical SIAC retrieval", "Final generic calibration"), strict=True
    ):
        ax.fill_between(x, np.maximum(0, x - (0.05 + 0.15 * x)), x + 0.05 + 0.15 * x, color="#e5e9eb")
        ax.plot(x, x, color="#263238", linewidth=1.2)
        ax.scatter(truth, values, s=22, color="#1769aa", alpha=0.72, edgecolors="white", linewidths=0.3)
        ax.set(xlim=(0, upper), ylim=(0, upper), xlabel="AERONET AOD", ylabel="Retrieved AOD", title=title)
        ax.grid(alpha=0.18)
        ax.set_aspect("equal", adjustable="box")
    fig.suptitle("AOD retrievals against AERONET; grey band is expected error", x=0.02, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(output / "assets/performance_scatter.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    ordered = sorted(rows, key=lambda row: abs(row["baseline_error_over_ee"]), reverse=True)
    baseline_ratio = np.asarray([row["baseline_error_over_ee"] for row in ordered])
    final_ratio = np.asarray([row["adjusted_error_over_ee"] for row in ordered])
    positions = np.arange(len(ordered))
    fig, ax = plt.subplots(figsize=(15.5, 6.2), facecolor="white")
    ax.axhspan(-1, 1, color="#e5e9eb")
    for index in positions:
        ax.plot(
            [index, index],
            [baseline_ratio[index], final_ratio[index]],
            color="#aeb8be",
            linewidth=0.6,
        )
    ax.scatter(positions, baseline_ratio, s=12, color="#8a959b", label="Physical SIAC")
    ax.scatter(positions, final_ratio, s=16, color="#1769aa", label="Final")
    ax.axhline(0, color="#263238", linewidth=0.8)
    ax.set(xlabel="Cases ordered by physical absolute error / EE", ylabel="Signed error / EE", ylim=(-8, 5))
    ax.legend(frameon=False, ncol=2)
    ax.grid(axis="y", alpha=0.18)
    fig.tight_layout()
    fig.savefig(output / "assets/error_transition.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


def build(root: Path, output: Path) -> dict[str, Any]:
    final = _load_json(root / "reports" / FINAL_REPORT.name)
    selection = _load_json(root / "reports" / SELECTION_REPORT.name)
    manifest = _load_json(root / "reports" / MODEL_MANIFEST.name)
    ablation = _load_json(root / "reports" / ABLATION_REPORT.name)
    legacy_data = _load_json(LEGACY_EXPLORER / "data/cases.json")
    legacy = {row["matchup_id"]: row for row in legacy_data["cases"]}
    current = _load_records(TARGET_RESULTS)
    priors = _load_records(PRIOR_RESULTS)
    samples = load_samples(
        TARGET_RESULTS,
        TARGET_CONTEXT,
        root / "matchups/matchups.csv",
        _ids(TARGET_MIDS),
        atmo_context_dir=TARGET_ATMO,
        include_geography=False,
        require_complete=True,
    )
    sample_map = {sample.matchup_id: sample for sample in samples}
    receipts = final["case_receipts"]
    if [row["matchup_id"] for row in receipts] != [sample.matchup_id for sample in samples]:
        raise ValueError("Final receipts and operational sample order do not agree")

    output.mkdir(parents=True, exist_ok=True)
    for directory in ("assets/spatial", "assets/diagnostic", "data"):
        (output / directory).mkdir(parents=True, exist_ok=True)
    shutil.copy2(WEB_ASSETS / "app.css", output / "app.css")
    shutil.copy2(WEB_ASSETS / "app.js", output / "app.js")

    rows: list[dict[str, Any]] = []
    new_diagnostic_count = 0
    legacy_diagnostic_count = 0
    for receipt in receipts:
        matchup_id = receipt["matchup_id"]
        sample: CalibrationSample = sample_map[matchup_id]
        features = {name: float(value) for name, value in sample.features.items()}
        record = current[matchup_id]
        source = _new_diagnostic_source(matchup_id)
        diagnostic = None
        if source is not None:
            diagnostic = _build_new_diagnostic(
                output,
                receipt,
                record,
                priors[matchup_id],
                source[0],
                source[1],
            )
            new_diagnostic_count += 1
        elif matchup_id in legacy:
            diagnostic = _copy_legacy_diagnostic(output, LEGACY_EXPLORER, legacy[matchup_id])
            legacy_diagnostic_count += 1
        rows.append(
            {
                **receipt,
                "cloud_cover": features.get("metadata_scene_cloud_cover"),
                "elevation_m": features.get("metadata_elevation_m"),
                "cams_aod": features.get(
                    "context_cams_total_aerosol_optical_depth_at_550nm_surface"
                ),
                "maiac_aod": features.get("atmo_aot_mean"),
                "band_argmin_spread": features.get("solver_surface_band_argmin_spread"),
                "solver_cost": features.get("solver_surface_cost_per_band"),
                "valid_aot_fraction": features.get("siac_valid_aot_fraction"),
                "surface_b02": features.get("surface_prior_boa_B02"),
                "surface_b03": features.get("surface_prior_boa_B03"),
                "surface_b04": features.get("surface_prior_boa_B04"),
                "surface_unc_b02": features.get("surface_prior_unc_B02"),
                "surface_unc_b03": features.get("surface_prior_unc_B03"),
                "surface_unc_b04": features.get("surface_prior_unc_B04"),
                "operational_features": features,
                "solver": record.get("solver") or {},
                "retrieval_extraction": record.get("retrieval_extraction") or {},
                "prior_boa": record.get("prior_boa") or {},
                "prior_unc": record.get("prior_unc") or {},
                "anchor_iterate": record.get("anchor_iterate") or {},
                "diagnostic": diagnostic,
            }
        )

    changed_or_missed = [row for row in rows if row["transition"] != "retained_hit"]
    missing_diagnostics = [
        row["matchup_id"] for row in changed_or_missed if row["diagnostic"] is None
    ]
    if missing_diagnostics:
        raise ValueError(f"Transition diagnostics are incomplete: {missing_diagnostics}")

    all_metrics = _metric_block(rows)
    cohorts = {
        "development_folds_0_3": _metric_block(
            [row for row in rows if row["site_fold"] != 4]
        ),
        "untouched_confirmation_fold_4": _metric_block(
            [row for row in rows if row["site_fold"] == 4]
        ),
        "external_training_site_seen": _metric_block(
            [row for row in rows if row["external_site_seen"]]
        ),
        "external_training_site_unseen": _metric_block(
            [row for row in rows if not row["external_site_seen"]]
        ),
    }
    folds = [
        {"fold": fold, **_metric_block([row for row in rows if row["site_fold"] == fold])}
        for fold in range(5)
    ]
    segments = {
        "truth_aod": binned_summaries(
            rows,
            "truth",
            (0, 0.1, 0.2, 0.4, 0.6, 1.0, math.inf),
            ("<0.10", "0.10-0.20", "0.20-0.40", "0.40-0.60", "0.60-1.00", ">=1.00"),
        ),
        "cloud_cover": binned_summaries(
            rows,
            "cloud_cover",
            (0, 2, 5, 10, 15, 20),
            ("0-2%", "2-5%", "5-10%", "10-15%", "15-20%"),
        ),
        "baseline_aod": binned_summaries(
            rows,
            "baseline",
            (0, 0.1, 0.2, 0.4, 0.6, 1.0, math.inf),
            ("<0.10", "0.10-0.20", "0.20-0.40", "0.40-0.60", "0.60-1.00", ">=1.00"),
        ),
    }
    offset_curve = _offset_curve(rows, final["selection_policy"]["offset_grid"])
    seed_rows = []
    for name, variant in final["seed_tree_replay"].items():
        seed_rows.append(
            {
                "variant": name,
                "target": variant["target_all"]["candidate"],
                "development": variant["target_development"]["candidate"],
                "confirmation": variant["target_confirmation"]["candidate"],
                "seen_sites": variant["target_seen_sites"]["candidate"],
                "unseen_sites": variant["target_unseen_sites"]["candidate"],
                "external_holdout": variant["external_holdout"]["candidate"],
            }
        )
    data = {
        "schema_version": "siac-final-aod-dashboard-v1",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "title": "Low-cloud AOD performance investigation",
        "subtitle": "One S2 monthly SWIR/NIR-anchored surface prior; one uniform operational calibrator",
        "benchmark": {
            "count": len(rows),
            "cloud_rule": "catalogue scene cloud cover < 20%",
            "expected_error": "abs(retrieved - AERONET) <= 0.05 + 0.15 * AERONET",
            "target_hits": 133,
            "target_rate": 87.0,
            "all_retrieval_jobs_terminal_ok": True,
        },
        "method": {
            "surface_prior": "S2 monthly SWIR/NIR-anchored ExtraTree surface prediction",
            "physical_retrieval": "B02/B03/B04 chi-square with scene-mean pooling",
            "calibrator": selection["selected_recipe"],
            "global_log_offset": final["selected_offset"],
            "formula": manifest["formula"],
            "model_sha256": manifest["model_sha256"],
            "model_url": "../../models/aod_et35_global_offset_20260713.joblib",
            "feature_schema_count": manifest["features"]["schema_count"],
            "selected_feature_count": manifest["features"]["selected_count"],
            "training": manifest["training"],
            "inference_constraints": {
                "one_surface_prior": True,
                "case_routing": False,
                "geography_features": False,
                "aeronet_at_inference": False,
                "target_label_at_inference": False,
            },
        },
        "metrics": {
            "all": all_metrics,
            "development": cohorts["development_folds_0_3"],
            "confirmation": cohorts["untouched_confirmation_fold_4"],
            "seen_sites": cohorts["external_training_site_seen"],
            "unseen_sites": cohorts["external_training_site_unseen"],
            "nested_offset": final["target_nested_offset_selection"]["aggregate"],
            "external_holdout": final["external_holdout"],
        },
        "cohorts": cohorts,
        "folds": folds,
        "segments": segments,
        "robustness": {
            "offset_curve": offset_curve,
            "nested_folds": final["target_nested_offset_selection"]["folds"],
            "seed_tree": seed_rows,
            "all_replayed_variants_meet_target": final[
                "all_replayed_variants_meet_target"
            ],
            "feature_ablation": ablation,
            "model_feature_importance": manifest["features"]["importance"],
            "artifact_reproduction": manifest["validation_reproduction"],
        },
        "transition_counts": {
            transition: sum(row["transition"] == transition for row in rows)
            for transition in ("gain", "loss", "remaining_miss", "retained_hit")
        },
        "diagnostics": {
            "case_count": sum(row["diagnostic"] is not None for row in rows),
            "transition_case_count": len(changed_or_missed),
            "new_exact_count": new_diagnostic_count,
            "legacy_physical_count": legacy_diagnostic_count,
            "missing_transition_diagnostics": missing_diagnostics,
        },
        "experiments": _experiment_ledger(root),
        "cases": _jsonable(rows),
        "sources": {
            "final_validation": "../aod-global-offset-validation-20260713.json",
            "recipe_selection": "../aod-target-domain-recipe-selection-20260713.json",
            "seed_robustness": "../aod-selected-et-seed-robustness-20260713.json",
            "feature_ablation": "../aod-generic-feature-ablation-20260713.json",
            "model_manifest": "../aod-et35-global-offset-model-manifest-20260713.json",
            "case_csv": "data/all-cases.csv",
        },
        "disclosures": [
            "The 588-record external cohort fits the calibrator; target folds 0-3 select the fixed recipe and scalar offset; target fold 4 is untouched confirmation.",
            "The nested 134/152 result validates scalar-offset selection for the already-fixed ExtraTrees recipe; separate recipe-family nested selection was 130/152.",
            "External-site-unseen target cases are 46/53 (86.79%), below the aggregate 87% threshold by one case.",
            "The external held-out fold is 83/123 (67.48%); it improves on its 66/123 physical baseline but is below the target-domain rate.",
            "Feature importance is model dependence, not causal attribution. UTC phase can encode acquisition context and may also proxy longitude.",
            "The no-MAIAC ablation reaches 135/152 on one seed but falls to 131-132 on other seeds; it is shown as sensitivity evidence, not the final method.",
            "Legacy physical diagnostic images use the earlier 2026-07-12 cost-cube snapshot. Ten newly missing gain/loss examples use exact final diagnostic reruns.",
        ],
    }
    (output / "data/dashboard.json").write_text(
        json.dumps(data, separators=(",", ":")), encoding="utf-8"
    )
    with (output / "data/all-cases.csv").open("w", encoding="utf-8", newline="") as stream:
        fields = (
            "matchup_id",
            "site",
            "site_fold",
            "external_site_seen",
            "cloud_cover",
            "truth",
            "baseline",
            "unadjusted",
            "adjusted",
            "baseline_error_over_ee",
            "adjusted_error_over_ee",
            "baseline_within_ee",
            "adjusted_within_ee",
            "transition",
            "cams_aod",
            "maiac_aod",
            "band_argmin_spread",
            "solver_cost",
        )
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows({field: row.get(field) for field in fields} for row in rows)
    _save_global_figures(output, rows)
    css_version = hashlib.sha256((output / "app.css").read_bytes()).hexdigest()[:12]
    js_version = hashlib.sha256((output / "app.js").read_bytes()).hexdigest()[:12]
    index = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <meta name="color-scheme" content="light dark">
  <title>Low-cloud AOD performance investigation</title>
  <link rel="icon" href="data:,">
  <link rel="stylesheet" href="app.css?v={css_version}">
</head>
<body>
  <div id="app" class="app-shell"><div class="loading-state">Loading the 152-case investigation snapshot...</div></div>
  <noscript><main class="noscript"><h1>Low-cloud AOD performance investigation</h1><p>The interactive evidence workspace requires JavaScript.</p><img src="assets/performance_scatter.png" alt="Physical and final AOD scatter plots"><img src="assets/error_transition.png" alt="Per-case error transitions"></main></noscript>
  <script src="app.js?v={js_version}"></script>
</body>
</html>
"""
    (output / "index.html").write_text(index, encoding="utf-8")
    receipt = {
        "output": str(output),
        "cases": len(rows),
        "target_hits": all_metrics["final"]["hits"],
        "diagnostic_cases": data["diagnostics"]["case_count"],
        "transition_diagnostic_cases": data["diagnostics"]["transition_case_count"],
        "spatial_images": len(list((output / "assets/spatial").glob("*.png"))),
        "diagnostic_images": len(list((output / "assets/diagnostic").glob("*.png"))),
        "missing_transition_diagnostics": missing_diagnostics,
    }
    (output / "build-receipt.json").write_text(
        json.dumps(receipt, indent=2) + "\n", encoding="utf-8"
    )
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=ROOT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    print(json.dumps(build(args.root, args.output), indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
