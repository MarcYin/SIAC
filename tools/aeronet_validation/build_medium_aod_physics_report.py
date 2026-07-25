"""Build the locked-development medium-AOD physical investigation webpage."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import shutil
import subprocess
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from tools.aeronet_validation.analyze_medium_aod_physics import (  # noqa: E402
    DEFAULT_EXCLUDED,
    DEFAULT_HOLDOUT_FOLDS,
    DEFAULT_SPLIT_SEED,
    _site_group_folds,
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
MIDS = ROOT / "campaign250_lowcloud20_mids.txt"
ARCHIVED = ROOT / "phaseD_results_lowcloud20_native_maiac_adaptive_b03_chi2_20260713"
FRESH = ROOT / (
    "phaseD_results_lowcloud20_mediumphysics_baseline_costcube_mediumdev_20260713"
)
CUBES = ROOT / (
    "phaseD_cost_cubes_lowcloud20_mediumphysics_baseline_mediumdev_20260713"
)
MATCHUPS = ROOT / "reports/aod-production-reproduction-spec-20260713/downloads/target-matchups.csv"
DEFAULT_OUTPUT = ROOT / "reports/aod-medium-physics-20260713"
WEB_ASSETS = Path(__file__).with_name("medium_aod_physics_report")

IMPLEMENTATION_FILES = (
    REPO_ROOT / "python/siac/adapters/satellite/sentinel2.py",
    REPO_ROOT / "python/siac/adapters/bestpixel_surface_prior.py",
    REPO_ROOT / "python/siac/algorithms/surface/seasonal_predictor.py",
    REPO_ROOT / "python/siac/algorithms/solver/surface_driven.py",
    REPO_ROOT / "python/siac/algorithms/rt/lut/backend.py",
    REPO_ROOT / "python/siac/config/public.py",
    REPO_ROOT / "src/siac_rs/src/surface_driven.rs",
    REPO_ROOT / "tools/aeronet_validation/campaign250_medium_aod_physics_combinations_submit.sbatch",
)

EXCLUDED = set(DEFAULT_EXCLUDED)
TARGET_RATE = 0.87
TARGET_DEV_HITS = 32

VARIANTS: tuple[dict[str, Any], ...] = (
    {
        "id": "archived",
        "label": "Archived current SIAC",
        "path": ARCHIVED,
        "family": "Reference",
        "description": "Frozen 152-case low-cloud snapshot used to lock the split.",
    },
    {
        "id": "fresh",
        "label": "Fresh current-code replay",
        "path": FRESH,
        "family": "Reference",
        "description": "Exact current code with cost-cube capture on the 36 development cases.",
    },
    {
        "id": "conflict",
        "label": "Prior-conflict selector",
        "path": ROOT
        / "analysis/medium_aod_current_end_to_end_prior_conflict_z2576_development_20260713",
        "family": "Atmospheric prior",
        "description": "Uses 50%-relative MAIAC uncertainty only when the surface optimum is >2.576 calibrated sigmas above MAIAC.",
    },
    {
        "id": "footprint_conflict",
        "label": "Matched-footprint conflict",
        "path": ROOT
        / "analysis/medium_aod_footprint0045_prior_conflict_z2576_development_20260713",
        "family": "Atmospheric prior",
        "description": "Same rule after matching the staged MAIAC footprint to the S2 retrieval AOI.",
    },
    {
        "id": "footprint_flat",
        "label": "Matched footprint, loose prior",
        "path": ROOT
        / "phaseD_results_lowcloud20_mediumphysics_maiac_footprint0045_product_unc_mediumdev_20260713",
        "family": "Atmospheric prior",
        "description": "AOI-matched MAIAC centre with 50%-relative uncertainty for every case.",
    },
    {
        "id": "linear_flat",
        "label": "Continuous geometry, loose prior",
        "path": ROOT
        / "phaseD_results_lowcloud20_mediumphysics_lut_linear_geometry_product_unc_mediumdev_20260713",
        "family": "RT geometry",
        "description": "Continuous LUT geometry interpolation with the same surface and loose MAIAC prior.",
    },
    {
        "id": "tau",
        "label": "Tau-dependent surface",
        "path": ROOT
        / "phaseD_results_lowcloud20_mediumphysics_tau_always_mediumdev_20260713",
        "family": "Surface iteration",
        "description": "Recomputes target anchors and visible surface at every candidate AOD.",
    },
    {
        "id": "ridge",
        "label": "Polynomial ridge surface",
        "path": ROOT
        / "phaseD_results_lowcloud20_mediumphysics_ridge_predictor_mediumdev_20260713",
        "family": "Surface model",
        "description": "Coefficient-stable ridge replacement for each historical realization predictor.",
    },
    {
        "id": "oof_global",
        "label": "Global OOF uncertainty",
        "path": ROOT
        / "analysis/medium_aod_surface_oof_unc_adaptive_z2576_development_20260713",
        "family": "Surface uncertainty",
        "description": "Single-tree historical OOF error added by band; retained as the original diagnostic snapshot.",
    },
    {
        "id": "oof_map_et20",
        "label": "Local ET20 OOF map",
        "path": ROOT
        / "analysis/medium_aod_surface_oof_map_et20_bias_unc_adaptive_development_20260713",
        "family": "Surface uncertainty",
        "description": "Local ET20 OOF bias and uncertainty map while retaining the current surface correction.",
    },
    {
        "id": "oof_no_cal",
        "label": "Coefficient-free ET20 OOF",
        "path": ROOT
        / "analysis/medium_aod_surface_oof_map_et20_oof_corrected_unc_adaptive_development_20260713",
        "family": "Surface model",
        "description": "Removes the existing site-derived surface debias and replaces it with historical S2 OOF correction.",
    },
    {
        "id": "pooled",
        "label": "Pooled historical ET20",
        "path": ROOT
        / "analysis/medium_aod_pooled_et20_surface_baseunc_calibrated_development_20260713",
        "family": "Surface model",
        "description": "One coefficient-free forest trained across every usable historical realization pixel.",
    },
    {
        "id": "pooled_blend",
        "label": "Current + pooled surface average",
        "path": ROOT
        / "analysis/medium_aod_pooled_et20_surface_blend50_current_mixunc_adaptive_development_20260713",
        "family": "Surface model",
        "description": "One final equal-weight S2 surface with model disagreement included in uncertainty.",
    },
    {
        "id": "clean_day",
        "label": "Absolute clean-day gate",
        "path": ROOT
        / "phaseD_results_lowcloud20_mediumphysics_clean_day_aod015_mediumdev_20260713",
        "family": "Surface sampling",
        "description": "Historical S2 surface days restricted to staged MAIAC AOD <=0.15.",
    },
    {
        "id": "cloud_bypass",
        "label": "Cloud-mask bypass diagnostic",
        "path": ROOT
        / "phaseD_results_lowcloud20_mediumphysics_allow_cloud_retrieval_mediumdev_20260713",
        "family": "Masking",
        "description": "Targeted six-case diagnostic retaining water/no-data masks while relaxing cloud exclusion.",
    },
    {
        "id": "cci3",
        "label": "Tile-wide CCI-3 species",
        "path": ROOT
        / "phaseD_results_lowcloud20_aerosol_species_sixs_cci3_cci3_tile_mediumdev_adaptive_20260713",
        "family": "Aerosol optics",
        "description": "One aerosol family for the whole tile, selected by aggregate physical cost.",
    },
)

JOB_IDS = (
    "39192707",
    "39200880",
    "39200976",
    "39201188",
    "39202016",
    "39202191",
    "39202192",
    "39202438",
    "39202441",
    "39203399",
    "39203668",
    "39203788",
    "39204054",
)

JOB_REPAIRS = {
    "39200976": ("39202441",),
    "39202016": ("39202438",),
}


def _finite(value: Any) -> float | None:
    try:
        output = float(value)
    except (TypeError, ValueError):
        return None
    return output if math.isfinite(output) else None


def _load_result(path: Path) -> dict[str, Any] | None:
    if not path.exists():
        return None
    try:
        return json.loads(path.read_text())
    except (OSError, json.JSONDecodeError):
        return None


def _load_records(directory: Path, ids: list[str]) -> tuple[dict[str, dict[str, Any]], Counter[str]]:
    records: dict[str, dict[str, Any]] = {}
    statuses: Counter[str] = Counter()
    for matchup_id in ids:
        record = _load_result(directory / f"{matchup_id}.json")
        if record is None:
            statuses["MISSING"] += 1
            continue
        status = str(record.get("status", "UNKNOWN")).upper()
        statuses[status] += 1
        if status == "OK" and _finite(record.get("retrieved")) is not None:
            records[matchup_id] = record
    return records, statuses


def _load_matchups() -> dict[str, dict[str, str]]:
    with MATCHUPS.open(newline="") as stream:
        return {row["matchup_id"]: row for row in csv.DictReader(stream)}


def _metrics(records: dict[str, dict[str, Any]], ids: list[str]) -> dict[str, Any]:
    metric = compute_metrics(records[mid] for mid in ids if mid in records).as_dict()
    if metric["within_ee_rate"] is not None:
        metric["within_ee_percent"] = 100.0 * metric["within_ee_rate"]
    return metric


def _transitions(
    reference: dict[str, dict[str, Any]],
    candidate: dict[str, dict[str, Any]],
    ids: list[str],
) -> dict[str, int]:
    counts = Counter()
    for matchup_id in ids:
        if matchup_id not in reference or matchup_id not in candidate:
            counts["unpaired"] += 1
            continue
        truth = float(reference[matchup_id]["truth"])
        old = within_ee(float(reference[matchup_id]["retrieved"]), truth)
        new = within_ee(float(candidate[matchup_id]["retrieved"]), truth)
        counts[
            "gain" if new and not old else "loss" if old and not new else "stable_hit" if old else "stable_miss"
        ] += 1
    return dict(counts)


def _candidate(
    label: str,
    value: Any,
    truth: float,
    *,
    source: str,
) -> dict[str, Any]:
    number = _finite(value)
    return {
        "label": label,
        "value": number,
        "error": None if number is None else number - truth,
        "error_over_ee": None
        if number is None
        else (number - truth) / (0.05 + 0.15 * truth),
        "within_ee": None if number is None else within_ee(number, truth),
        "source": source,
    }


def _case_diagnostic_figure(
    path: Path,
    *,
    cost: np.lib.npyio.NpzFile,
    record: dict[str, Any],
    local_mask: np.ndarray,
    candidates: list[dict[str, Any]],
    render: bool = True,
) -> dict[str, Any]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    axis = np.asarray(cost["aot_axis"], dtype=np.float64)
    cube = np.asarray(cost["cube"], dtype=np.float64)
    band_cost = np.asarray(cost["band_cost_cube"], dtype=np.float64)
    band_residual = np.asarray(cost["band_residual_cube"], dtype=np.float64)
    band_signed = np.asarray(cost["band_signed_residual_cube"], dtype=np.float64)
    valid = np.asarray(cost["solve_valid"], dtype=bool)
    curve_valid = local_mask & valid
    scope = "1.5 km station window"
    if int(np.count_nonzero(curve_valid)) < 20:
        curve_valid = valid
        scope = "scene support"
    truth = float(record["truth"])
    retrieved = float(record["retrieved"])
    ee = 0.05 + 0.15 * truth
    atmo = float(np.nanmedian(np.asarray(cost["aot_prior"])[curve_valid]))
    anchor = _finite((record.get("anchor_iterate") or {}).get("pass1_scene_mean"))
    total_local = _median_curve(cube, curve_valid)
    total_scene = _median_curve(cube, valid)
    per_band = [_median_curve(band_cost[index], curve_valid) for index in range(3)]
    residual = [_median_curve(band_residual[index], curve_valid) for index in range(3)]
    signed = [_median_curve(band_signed[index], curve_valid) for index in range(3)]
    toa = np.asarray(cost["toa"], dtype=np.float64)
    prior = np.asarray(cost["boa_prior"], dtype=np.float64)
    prior_unc = np.asarray(cost["boa_unc"], dtype=np.float64)
    wavelengths = np.asarray([WAVELENGTHS[band] for band in BANDS])
    toa_median = np.asarray([np.nanmedian(toa[index][curve_valid]) for index in range(3)])
    prior_median = np.asarray([np.nanmedian(prior[index][curve_valid]) for index in range(3)])
    prior_unc_median = np.asarray(
        [np.nanmedian(prior_unc[index][curve_valid]) for index in range(3)]
    )
    summary = {
        "scope": scope,
        "support_count": int(np.count_nonzero(curve_valid)),
        "maiac_aod": atmo,
        "surface_anchor_aod": anchor,
        "toa_median": dict(zip(BANDS, toa_median.tolist(), strict=True)),
        "surface_median": dict(zip(BANDS, prior_median.tolist(), strict=True)),
        "surface_unc_median": dict(zip(BANDS, prior_unc_median.tolist(), strict=True)),
    }
    if not render:
        return summary

    def context(ax: Any) -> None:
        ax.axvspan(truth - ee, truth + ee, color="#e8ecee", alpha=0.9, label="EE")
        ax.axvline(truth, color="#111111", linewidth=1.5, label="AERONET")
        ax.axvline(retrieved, color="#2f66a3", linestyle="--", label="Fresh SIAC")
        ax.axvline(atmo, color="#bb4f2b", linestyle=":", label="MAIAC")
        if anchor is not None:
            ax.axvline(anchor, color="#0b6e69", linestyle="-.", label="Final surface anchor")

    fig, axes = plt.subplots(2, 3, figsize=(16, 9.4), facecolor="white")
    axes = axes.ravel()
    axes[0].plot(axis, _normalise_curve(total_scene), color="#6d7478", label="Scene")
    axes[0].plot(axis, _normalise_curve(total_local), color="#2f66a3", label=scope)
    context(axes[0])
    axes[0].set(title="Total surface likelihood", xlabel="AOD550", ylabel="log10(1 + delta cost)")
    axes[0].legend(fontsize=7, ncol=2)

    for index, band in enumerate(BANDS):
        axes[1].plot(axis, _normalise_curve(per_band[index]), color=BAND_COLORS[band], label=band)
    context(axes[1])
    axes[1].set(title="Per-band cost profiles", xlabel="AOD550", ylabel="log10(1 + delta cost)")
    axes[1].legend(fontsize=7, ncol=3)

    for index, band in enumerate(BANDS):
        axes[2].plot(axis, signed[index], color=BAND_COLORS[band], label=band)
    context(axes[2])
    axes[2].axhline(0.0, color="#777777", linewidth=0.8)
    axes[2].set(title="Signed BOA minus surface residual", xlabel="AOD550", ylabel="Reflectance")
    axes[2].legend(fontsize=7, ncol=3)

    axes[3].plot(wavelengths, toa_median, "o-", color="#222222", label="TOA")
    axes[3].errorbar(
        wavelengths,
        prior_median,
        yerr=prior_unc_median,
        fmt="o-",
        color="#2f66a3",
        capsize=4,
        label="Surface prior",
    )
    axes[3].set_xticks(wavelengths, BANDS)
    axes[3].set(title=f"Reflectance spectrum: {scope}", ylabel="Reflectance")
    axes[3].legend(fontsize=8)

    truth_index = int(np.argmin(np.abs(axis - truth)))
    current_index = int(np.argmin(np.abs(axis - retrieved)))
    positions = np.arange(3)
    axes[4].bar(
        positions - 0.18,
        [residual[index][current_index] for index in range(3)],
        width=0.36,
        color="#2f66a3",
        label="Fresh SIAC node",
    )
    axes[4].bar(
        positions + 0.18,
        [residual[index][truth_index] for index in range(3)],
        width=0.36,
        color="#bb4f2b",
        label="AERONET node",
    )
    axes[4].plot(positions, prior_unc_median, "D", color="#111111", label="Prior sigma")
    axes[4].set_xticks(positions, BANDS)
    axes[4].set(title="Spectral closure", ylabel="Absolute reflectance residual")
    axes[4].legend(fontsize=7)

    plotted = [item for item in candidates if item["value"] is not None][:18]
    ypos = np.arange(len(plotted))
    values = [item["value"] for item in plotted]
    axes[5].axvspan(truth - ee, truth + ee, color="#e8ecee", alpha=0.9)
    axes[5].axvline(truth, color="#111111", linewidth=1.5)
    axes[5].scatter(
        values,
        ypos,
        c=["#0b6e69" if item["within_ee"] else "#bb4f2b" for item in plotted],
        s=30,
        zorder=3,
    )
    for index, value in enumerate(values):
        axes[5].plot([truth, value], [index, index], color="#c8ced1", linewidth=0.8)
    axes[5].set_yticks(ypos, [item["label"] for item in plotted], fontsize=7)
    axes[5].invert_yaxis()
    axes[5].set(title="Scalar estimates", xlabel="AOD550")

    x_max = min(float(axis[-1]), max(1.25, truth * 1.5))
    for index, ax in enumerate(axes):
        if index in {0, 1, 2}:
            ax.set_xlim(float(axis[0]), x_max)
        ax.grid(alpha=0.18)
        ax.tick_params(labelsize=8)
    fig.suptitle(f"{record['matchup_id']} | spectral and cost evidence", x=0.02, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=125, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return summary


def _surface_model_figure(
    path: Path,
    *,
    cost: np.lib.npyio.NpzFile,
    pooled_path: Path,
    oof_path: Path,
    title: str,
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    current = np.asarray(cost["boa_prior"], dtype=np.float64)
    current_unc = np.asarray(cost["boa_unc"], dtype=np.float64)
    with np.load(pooled_path, allow_pickle=False) as pooled_data:
        pooled = np.asarray(pooled_data["surface"], dtype=np.float64)
        pooled_sigma = np.asarray(pooled_data["tree_sigma"], dtype=np.float64)
    with np.load(oof_path, allow_pickle=False) as oof_data:
        oof_bias = np.asarray(oof_data["bias"], dtype=np.float64)
        oof_sigma = np.asarray(oof_data["sigma"], dtype=np.float64)

    fig, axes = plt.subplots(2, 4, figsize=(17.5, 8.2), facecolor="white")
    panels = (
        (current[0], "Current B02 surface", "viridis", 0.0, 0.3),
        (pooled[0], "Pooled B02 surface", "viridis", 0.0, 0.3),
        (pooled[0] - current[0], "Pooled - current B02", "RdBu_r", -0.05, 0.05),
        (oof_bias[0], "Historical OOF B02 bias", "RdBu_r", -0.04, 0.04),
        (current_unc[0], "Current B02 uncertainty", "magma", 0.0, 0.05),
        (pooled_sigma[0], "Pooled tree B02 spread", "magma", 0.0, 0.05),
        (oof_sigma[0], "Historical OOF B02 spread", "magma", 0.0, 0.05),
        (np.asarray(cost["solve_valid"], dtype=float), "Solver support", "Greys", 0.0, 1.0),
    )
    for ax, (values, panel_title, cmap, vmin, vmax) in zip(axes.ravel(), panels, strict=True):
        image = ax.imshow(values, cmap=cmap, vmin=vmin, vmax=vmax, interpolation="nearest")
        ax.set_title(panel_title, fontsize=9)
        ax.set_xticks([])
        ax.set_yticks([])
        fig.colorbar(image, ax=ax, fraction=0.045, pad=0.02)
    fig.suptitle(f"{title} | surface-model evidence", x=0.02, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=125, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def _global_figures(output: Path, cases: list[dict[str, Any]], variants: list[dict[str, Any]]) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    complete = [row for row in variants if row["metrics"]["n"] == 36]
    complete = sorted(complete, key=lambda row: row["metrics"]["within_ee_rate"], reverse=True)
    fig, ax = plt.subplots(figsize=(12.5, 7.2), facecolor="white")
    labels = [row["label"] for row in complete]
    rates = [100.0 * row["metrics"]["within_ee_rate"] for row in complete]
    colors = ["#0b6e69" if rate >= 87.0 else "#2f66a3" for rate in rates]
    ax.barh(np.arange(len(labels)), rates, color=colors)
    ax.axvline(87.0, color="#bb4f2b", linewidth=1.5, label="Acceptance target")
    ax.set_yticks(np.arange(len(labels)), labels, fontsize=8)
    ax.invert_yaxis()
    ax.set(xlabel="Within expected error (%)", xlim=(0, 100), title="Locked medium-AOD development performance")
    ax.grid(axis="x", alpha=0.2)
    ax.legend()
    fig.tight_layout()
    fig.savefig(output / "assets" / "variant-ranking.png", dpi=140, bbox_inches="tight")
    plt.close(fig)

    truth = np.asarray([row["truth"] for row in cases])
    fig, axes = plt.subplots(1, 3, figsize=(15.8, 5.2), facecolor="white", sharex=True, sharey=True)
    for ax, key, label in (
        (axes[0], "archived", "Archived"),
        (axes[1], "fresh", "Fresh current code"),
        (axes[2], "conflict", "Prior conflict"),
    ):
        values = np.asarray([row["candidate_by_id"][key]["value"] for row in cases])
        hit = np.asarray([row["candidate_by_id"][key]["within_ee"] for row in cases])
        grid = np.linspace(0.2, 0.9, 200)
        ax.fill_between(grid, grid - (0.05 + 0.15 * grid), grid + (0.05 + 0.15 * grid), color="#e8ecee")
        ax.plot(grid, grid, color="#111111", linewidth=1)
        ax.scatter(truth[hit], values[hit], color="#0b6e69", s=28, label="Within EE")
        ax.scatter(truth[~hit], values[~hit], color="#bb4f2b", s=28, label="Outside EE")
        ax.set(title=label, xlabel="AERONET AOD550", xlim=(0.2, 0.9), ylim=(0.0, 1.05))
        ax.grid(alpha=0.16)
    axes[0].set_ylabel("Retrieved AOD550")
    axes[0].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(output / "assets" / "retrieval-scatter.png", dpi=140, bbox_inches="tight")
    plt.close(fig)

    columns = ("fresh", "conflict", "maiac", "surface_min", "B02_min", "B03_min", "B04_min")
    matrix = np.asarray(
        [
            [row["candidate_by_id"][key]["error_over_ee"] for key in columns]
            for row in sorted(cases, key=lambda item: item["fresh_error_over_ee"])
        ],
        dtype=np.float64,
    )
    fig, ax = plt.subplots(figsize=(12.8, 10.5), facecolor="white")
    image = ax.imshow(matrix, cmap="RdBu_r", vmin=-3.0, vmax=3.0, aspect="auto")
    ordered = sorted(cases, key=lambda item: item["fresh_error_over_ee"])
    ax.set_xticks(np.arange(len(columns)), ["Fresh", "Conflict", "MAIAC", "Surface min", "B02 min", "B03 min", "B04 min"], rotation=25, ha="right")
    ax.set_yticks(np.arange(len(ordered)), [row["site"] for row in ordered], fontsize=7)
    ax.set_title("Signed offset from AERONET in expected-error units")
    fig.colorbar(image, ax=ax, label="Signed error / EE")
    fig.tight_layout()
    fig.savefig(output / "assets" / "evidence-matrix.png", dpi=140, bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(13.5, 5.8), facecolor="white")
    errors = np.asarray([row["fresh_error"] for row in cases])
    support = np.asarray([row["valid_fraction"] for row in cases])
    sizes = 30.0 + 120.0 * np.clip(support, 0.0, 1.0)
    scatter = ax.scatter(
        [row["lon"] for row in cases],
        [row["lat"] for row in cases],
        c=errors,
        s=sizes,
        cmap="RdBu_r",
        vmin=-0.3,
        vmax=0.3,
        edgecolor="white",
        linewidth=0.5,
    )
    ax.axhline(0, color="#cccccc", linewidth=0.6)
    ax.set(xlabel="Longitude", ylabel="Latitude", title="Spatial distribution of fresh-code medium-AOD error")
    ax.grid(alpha=0.15)
    fig.colorbar(scatter, ax=ax, label="Retrieved - AERONET AOD550")
    fig.tight_layout()
    fig.savefig(output / "assets" / "spatial-error.png", dpi=140, bbox_inches="tight")
    plt.close(fig)


def _job_audit() -> list[dict[str, Any]]:
    rows = []
    for job_id in JOB_IDS:
        result = subprocess.run(
            [
                "sacct",
                "-j",
                job_id,
                "--format=JobIDRaw,State,ExitCode,Elapsed",
                "-n",
                "-X",
            ],
            check=False,
            capture_output=True,
            text=True,
        )
        states: Counter[str] = Counter()
        exit_codes: Counter[str] = Counter()
        for line in result.stdout.splitlines():
            fields = line.split()
            if len(fields) >= 3:
                states[fields[1].split("+")[0]] += 1
                exit_codes[fields[2]] += 1
        repair_ids = list(JOB_REPAIRS.get(job_id, ()))
        rows.append(
            {
                "job_id": job_id,
                "states": dict(states),
                "exit_codes": dict(exit_codes),
                "query_ok": result.returncode == 0,
                "repair_job_ids": repair_ids,
                "final_state": "PENDING_AUDIT",
            }
        )
    by_id = {row["job_id"]: row for row in rows}
    for row in rows:
        failed = int(row["states"].get("FAILED", 0))
        repairs = [by_id.get(job_id) for job_id in row["repair_job_ids"]]
        repair_ok = bool(repairs) and all(
            repair is not None
            and repair["query_ok"]
            and repair["states"].get("COMPLETED", 0) > 0
            and sum(
                count
                for state, count in repair["states"].items()
                if state != "COMPLETED"
            )
            == 0
            for repair in repairs
        )
        row["final_state"] = (
            "REPAIRED" if failed and repair_ok else "UNRESOLVED" if failed else "COMPLETE"
        )
    return rows


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        return
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _copy_implementation_files(output: Path) -> list[dict[str, Any]]:
    rows = []
    for source in IMPLEMENTATION_FILES:
        relative = source.relative_to(REPO_ROOT)
        destination = output / "downloads/implementation" / relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, destination)
        payload = destination.read_bytes()
        rows.append(
            {
                "path": str(relative),
                "url": f"downloads/implementation/{relative}",
                "sha256": hashlib.sha256(payload).hexdigest(),
                "bytes": len(payload),
            }
        )
    manifest = output / "downloads/implementation/manifest.json"
    manifest.write_text(json.dumps(rows, indent=2) + "\n")
    return rows


def _inventory_family(name: str) -> str:
    lowered = name.lower()
    if any(token in lowered for token in ("multisource", "calibrat", "extra_tree", "classifier")):
        return "Learned / multi-source diagnostic"
    if lowered.startswith("lowcloud20_prior_conflict_"):
        return "Mixed-version prior diagnostic"
    if "aerosol_species" in lowered or "cci3" in lowered or "canonical" in lowered:
        return "Aerosol optics"
    if any(token in lowered for token in ("prior_conflict", "robust_cams", "same_cube", "hierarchical_tile")):
        return "Atmospheric prior"
    if any(token in lowered for token in ("oof", "pooled", "temporal_analog", "spectral_covariance", "regularized_offset")):
        return "Surface model / uncertainty"
    if any(token in lowered for token in ("maiac_footprint", "maiac_orbit", "maiac_native")):
        return "Atmospheric sampling"
    if any(token in lowered for token in ("tau", "anchor_converged")):
        return "Surface iteration"
    if any(token in lowered for token in ("nir", "swir", "b01", "additive_offset")):
        return "Bands / cost"
    if any(token in lowered for token in ("dem", "l2awvp", "tcwv")):
        return "Atmospheric state"
    if any(token in lowered for token in ("cloud", "water_buffer")):
        return "Masking"
    if "baseline" in lowered or "current" in lowered:
        return "Reference"
    return "Other physical diagnostic"


def _inventory_description(name: str) -> str:
    lowered = name.lower()
    if "multisource" in lowered:
        return "Learned selector over multiple saved retrieval sources; retained for comparison but outside the requested single-prior generic algorithm constraint."
    if any(token in lowered for token in ("calibrat", "extra_tree", "classifier")):
        return "Learned or fitted diagnostic output; retained for comparison and not treated as a coefficient-free physical recipe."
    if lowered.startswith("lowcloud20_prior_conflict_"):
        return "Diagnostic selector combining the frozen archived baseline with a separate current-code loose-prior rerun; not one end-to-end implementation snapshot."
    return "Saved physical or diagnostic experiment directory; inspect its path and paired transitions before treating it as an end-to-end recipe."


def _experiment_inventory(
    medium_ids: list[str],
    fresh_records: dict[str, dict[str, Any]],
    known_rows: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    known = {Path(row["path"]).resolve(): row for row in known_rows}
    directories = {
        path
        for path in ROOT.iterdir()
        if path.is_dir() and path.name.startswith("phaseD_results")
    }
    directories.update(path for path in (ROOT / "analysis").iterdir() if path.is_dir())
    rows: list[dict[str, Any]] = []
    seen: set[Path] = set()
    for directory in sorted(directories):
        if not directory.is_dir():
            continue
        resolved = directory.resolve()
        if resolved in seen:
            continue
        seen.add(resolved)
        if resolved in known:
            rows.append({**known[resolved], "inventory_kind": "curated"})
            continue
        records, statuses = _load_records(directory, medium_ids)
        if len(records) < 6:
            continue
        metric = _metrics(records, medium_ids)
        digest = hashlib.sha1(str(directory).encode()).hexdigest()[:12]
        rows.append(
            {
                "id": f"inventory_{digest}",
                "label": directory.name.replace("phaseD_results_lowcloud20_", "").replace(
                    "medium_aod_", ""
                ),
                "family": _inventory_family(directory.name),
                "description": _inventory_description(directory.name),
                "path": str(directory),
                "status_counts": dict(statuses),
                "metrics": metric,
                "complete": metric["n"] == len(medium_ids),
                "transitions_vs_fresh": _transitions(fresh_records, records, medium_ids),
                "inventory_kind": "discovered",
            }
        )
    for resolved, row in known.items():
        if resolved not in seen:
            rows.append({**row, "inventory_kind": "curated"})
    return rows


def build(output: Path, *, reuse_figures: bool = False) -> dict[str, Any]:
    mids = [line.strip() for line in MIDS.read_text().splitlines() if line.strip()]
    archived_all, archived_status = _load_records(ARCHIVED, mids)
    if len(archived_all) != 152:
        raise ValueError(f"archived cohort is incomplete: {dict(archived_status)}")
    folds = _site_group_folds(archived_all, mids, seed=DEFAULT_SPLIT_SEED)
    development_ids = [mid for mid in mids if folds[mid] not in DEFAULT_HOLDOUT_FOLDS]
    holdout_ids = [mid for mid in mids if folds[mid] in DEFAULT_HOLDOUT_FOLDS]
    medium_ids = [
        mid
        for mid in development_ids
        if 0.25 <= float(archived_all[mid]["truth"]) <= 0.85
    ]
    if len(medium_ids) != 36:
        raise ValueError(f"expected 36 locked development cases, got {len(medium_ids)}")
    matchups = _load_matchups()

    variant_records: dict[str, dict[str, dict[str, Any]]] = {}
    variant_rows: list[dict[str, Any]] = []
    fresh_records: dict[str, dict[str, Any]] = {}
    for spec in VARIANTS:
        records, statuses = _load_records(Path(spec["path"]), medium_ids)
        variant_records[spec["id"]] = records
        if spec["id"] == "fresh":
            fresh_records = records
        metric = _metrics(records, medium_ids)
        variant_rows.append(
            {
                **{key: spec[key] for key in ("id", "label", "family", "description")},
                "path": str(spec["path"]),
                "status_counts": dict(statuses),
                "metrics": metric,
                "complete": metric["n"] == 36,
            }
        )
    if len(fresh_records) != 36:
        raise ValueError("fresh current-code results must be 36/36 OK")
    for row in variant_rows:
        row["transitions_vs_fresh"] = _transitions(
            fresh_records, variant_records[row["id"]], medium_ids
        )
    inventory_rows = _experiment_inventory(medium_ids, fresh_records, variant_rows)

    output.mkdir(parents=True, exist_ok=True)
    for directory in (
        output / "assets",
        output / "assets/spatial",
        output / "assets/diagnostic",
        output / "assets/surface",
        output / "data",
        output / "downloads",
        output / "downloads/implementation",
    ):
        directory.mkdir(parents=True, exist_ok=True)
    shutil.copy2(WEB_ASSETS / "app.css", output / "app.css")
    shutil.copy2(WEB_ASSETS / "app.js", output / "app.js")
    implementation_files = _copy_implementation_files(output)

    pooled_fields = ROOT / "analysis/medium_aod_pooled_et20_surface_fields_development_20260713"
    oof_fields = ROOT / "analysis/medium_aod_surface_oof_map_et20_fields_development_20260713"
    cases: list[dict[str, Any]] = []
    for index, matchup_id in enumerate(medium_ids, start=1):
        fresh = fresh_records[matchup_id]
        truth = float(fresh["truth"])
        solver = fresh.get("solver") or {}
        atmo = fresh.get("atmo_prior") or {}
        candidate_list = [
            _candidate(
                next(row["label"] for row in variant_rows if row["id"] == variant_id),
                variant_records[variant_id].get(matchup_id, {}).get("retrieved"),
                truth,
                source=variant_id,
            )
            for variant_id in (
                "archived",
                "fresh",
                "conflict",
                "footprint_conflict",
                "linear_flat",
                "tau",
                "ridge",
                "oof_map_et20",
                "oof_no_cal",
                "pooled",
                "pooled_blend",
                "clean_day",
                "cci3",
            )
        ]
        candidate_list.extend(
            (
                _candidate("Staged MAIAC centre", atmo.get("aot_median"), truth, source="maiac"),
                _candidate("Surface-only curve minimum", solver.get("surface_cost_curve_min_aot"), truth, source="surface_min"),
                _candidate("B02 curve minimum", solver.get("surface_band_B02_argmin_aot"), truth, source="B02_min"),
                _candidate("B03 curve minimum", solver.get("surface_band_B03_argmin_aot"), truth, source="B03_min"),
                _candidate("B04 curve minimum", solver.get("surface_band_B04_argmin_aot"), truth, source="B04_min"),
            )
        )
        candidate_by_id = {item["source"]: item for item in candidate_list}

        cube_path = CUBES / f"{matchup_id}.npz"
        if not cube_path.exists():
            raise ValueError(f"missing cost cube: {cube_path}")
        with np.load(cube_path, allow_pickle=False) as cost:
            local_mask, _iy, _ix = _window_mask(
                np.asarray(cost["x"], dtype=np.float64),
                np.asarray(cost["y"], dtype=np.float64),
                lon=float(fresh["lon"]),
                lat=float(fresh["lat"]),
                crs=str(fresh["scene_crs"]),
            )
            pooled_aod, pooled_unc, _pooled_cost = _pooled_maps(cost)
            no_backstop, _no_unc, _no_cost = _pooled_maps(
                cost,
                prior_unc=np.full(np.asarray(cost["solve_valid"]).shape, np.inf),
            )
            spatial_path = output / "assets/spatial" / f"{matchup_id}.png"
            diagnostic_path = output / "assets/diagnostic" / f"{matchup_id}.png"
            surface_path = output / "assets/surface" / f"{matchup_id}.png"
            if not (reuse_figures and spatial_path.exists()):
                _save_spatial_figure(
                    spatial_path,
                    cost=cost,
                    record=fresh,
                    pooled_aod=pooled_aod,
                    pooled_unc=pooled_unc,
                    no_backstop_aod=no_backstop,
                    iy=_iy,
                    ix=_ix,
                )
            diagnostic_summary = _case_diagnostic_figure(
                diagnostic_path,
                cost=cost,
                record=fresh,
                local_mask=local_mask,
                candidates=candidate_list,
                render=not (reuse_figures and diagnostic_path.exists()),
            )
            if not (reuse_figures and surface_path.exists()):
                _surface_model_figure(
                    surface_path,
                    cost=cost,
                    pooled_path=pooled_fields / f"{matchup_id}.npz",
                    oof_path=oof_fields / f"{matchup_id}.npz",
                    title=matchup_id,
                )

        fresh_candidate = candidate_by_id["fresh"]
        conflict_candidate = candidate_by_id["conflict"]
        if not fresh_candidate["within_ee"] and conflict_candidate["within_ee"]:
            transition = "gain"
        elif fresh_candidate["within_ee"] and not conflict_candidate["within_ee"]:
            transition = "loss"
        elif fresh_candidate["within_ee"]:
            transition = "stable_hit"
        else:
            transition = "stable_miss"
        bands = [
            _finite(solver.get(f"surface_band_{band}_argmin_aot")) for band in BANDS
        ]
        finite_bands = [value for value in bands if value is not None]
        case = {
            "case_index": index,
            "matchup_id": matchup_id,
            "site": fresh["site"],
            "truth": truth,
            "ee": 0.05 + 0.15 * truth,
            "ee_low": truth - (0.05 + 0.15 * truth),
            "ee_high": truth + (0.05 + 0.15 * truth),
            "lon": float(fresh["lon"]),
            "lat": float(fresh["lat"]),
            "cloud_fraction": _finite(fresh.get("cloud_frac")),
            "scene_cloud_cover": _finite(matchups.get(matchup_id, {}).get("scene_cloud_cover")),
            "aeronet_count": _finite(matchups.get(matchup_id, {}).get("n_aeronet")),
            "aeronet_std": _finite(matchups.get(matchup_id, {}).get("aeronet_aod550_std")),
            "angstrom": _finite(matchups.get(matchup_id, {}).get("aeronet_angstrom_mean")),
            "elevation_m": _finite(matchups.get(matchup_id, {}).get("elevation_m")),
            "fresh": float(fresh["retrieved"]),
            "fresh_error": float(fresh["retrieved"]) - truth,
            "fresh_error_over_ee": (float(fresh["retrieved"]) - truth) / (0.05 + 0.15 * truth),
            "direction": "under" if float(fresh["retrieved"]) < truth else "over",
            "transition": transition,
            "maiac": _finite(atmo.get("aot_median")),
            "maiac_unc": _finite(atmo.get("aot_unc_median")),
            "surface_anchor_aod": _finite((fresh.get("anchor_iterate") or {}).get("pass1_scene_mean")),
            "surface_min": _finite(solver.get("surface_cost_curve_min_aot")),
            "band_min_spread": (
                max(finite_bands) - min(finite_bands) if finite_bands else None
            ),
            "cost_per_band": _finite(solver.get("surface_static_cost_per_band")),
            "curvature": _finite(solver.get("surface_cost_curve_curvature")),
            "valid_fraction": _finite(solver.get("surface_valid_observation_fraction")),
            "candidate_list": candidate_list,
            "candidate_by_id": candidate_by_id,
            "diagnostic": diagnostic_summary,
            "solver": solver,
            "atmo_prior": atmo,
            "prior_boa": fresh.get("prior_boa") or {},
            "prior_unc": fresh.get("prior_unc") or {},
            "anchor_iterate": fresh.get("anchor_iterate") or {},
            "spatial_image": f"assets/spatial/{matchup_id}.png",
            "diagnostic_image": f"assets/diagnostic/{matchup_id}.png",
            "surface_image": f"assets/surface/{matchup_id}.png",
            "result_url": f"../../{FRESH.name}/{matchup_id}.json",
            "cost_url": f"../../{CUBES.name}/{matchup_id}.npz",
        }
        cases.append(case)

    _global_figures(output, cases, variant_rows)

    subsets = {
        "low": [mid for mid in mids if float(archived_all[mid]["truth"]) < 0.25],
        "medium": [mid for mid in mids if 0.25 <= float(archived_all[mid]["truth"]) <= 0.85],
        "high": [mid for mid in mids if float(archived_all[mid]["truth"]) > 0.85],
        "non_extreme": [mid for mid in mids if mid not in EXCLUDED],
        "all": mids,
    }
    cohort_metrics = {name: _metrics(archived_all, ids) for name, ids in subsets.items()}
    best = max(
        (row for row in variant_rows if row["metrics"]["n"] == 36),
        key=lambda row: row["metrics"]["within_ee_rate"],
    )
    jobs = _job_audit()
    report = {
        "schema_version": 1,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "title": "Medium-AOD physical investigation",
        "subtitle": "Sentinel-2 SIAC against AERONET, cloud fraction below 20%",
        "target": {
            "within_ee_rate": TARGET_RATE,
            "development_required_hits": TARGET_DEV_HITS,
            "development_count": 36,
            "best_variant": best["label"],
            "best_hits": best["metrics"]["within_ee"],
            "best_rate": best["metrics"]["within_ee_rate"],
            "gap_cases": TARGET_DEV_HITS - best["metrics"]["within_ee"],
            "met": best["metrics"]["within_ee"] >= TARGET_DEV_HITS,
        },
        "cohort": {
            "count": 152,
            "definition": "retrieval cloud fraction < 0.20",
            "scene_metadata_cloud_below_20_count": 127,
            "medium_range": [0.25, 0.85],
            "ee_definition": "abs(retrieved - AERONET) <= 0.05 + 0.15*AERONET",
            "excluded_extreme_ids": sorted(EXCLUDED),
            "split_seed": DEFAULT_SPLIT_SEED,
            "holdout_folds": list(DEFAULT_HOLDOUT_FOLDS),
            "development_count": len(development_ids),
            "holdout_count": len(holdout_ids),
            "medium_development_count": len(medium_ids),
            "holdout_status": "locked; never scored because no recipe met development acceptance",
            "archived_metrics": cohort_metrics,
        },
        "variants": variant_rows,
        "experiment_inventory": inventory_rows,
        "cases": cases,
        "jobs": jobs,
        "reproduction": {
            "constraints": [
                "One historical Sentinel-2 L2A surface-prior source type",
                "No case-level source routing",
                "No AERONET value at inference",
                "One tile-wide aerosol family when species selection is enabled",
                "Site-group holdout remains inaccessible until development acceptance",
            ],
            "pipeline": [
                {"stage": "Cohort", "detail": "S2 L1C/AERONET +/-30 minute matchup; retrieval cloud fraction below 0.20; AOD550 EE rule."},
                {"stage": "Radiometry", "detail": "Read L1C DN, RADIO_ADD_OFFSET and QUANTIFICATION_VALUE; reflectance=(DN+offset)/quantification, clipped to [0,1.5]."},
                {"stage": "Masks", "detail": "OmniCloudMask plus water and no-data exclusions; masks resampled to the 60 m solve grid."},
                {"stage": "Surface history", "detail": "One source type: Planetary Computer S2 L2A. Use five complete prior years and target month +/-1 (15 year-month windows), cloud-cover ceiling 90%, top-k 15, and retain the lowest-MAIAC-AOD 60% of days per window."},
                {"stage": "Surface predictor", "detail": "For each realization fit ET20 (20 ExtraTrees, min_samples_leaf=5, random_state=0). Features are scene B8A/B11/B12 BOA plus four historical mean-visible localizer planes; targets are B02/B03/B04. Median predictions and MAD*1.4826 spread across realizations, floor 0.006."},
                {"stage": "Existing surface correction", "detail": "B02=(-0.0003+0.0243*AOD), B03=(-0.0006+0.0235*AOD), B04=(-0.0011+0.0223*AOD). This is a current-method reproduction dependency and was explicitly removed in coefficient-free tests."},
                {"stage": "Anchor iteration", "detail": "Correct target B8A/B11/B12 TOA to BOA in the same LUT space at staged MAIAC AOD, build/solve once, rebuild the same surface type at the pass-1 scene-mean AOD, and solve a second time."},
                {"stage": "Atmospheric prior", "detail": "Staged MCD19A2 MAIAC AOD550 supplies one tile-wide centre. Calibrated sigma is max(0.5m,0.02) below m=0.15; otherwise max(0.07, 0.5m/(1+exp(-(m-0.5)/0.15)))."},
                {"stage": "RT", "detail": "Packaged libRadtran continental-average 1 nm Zarr LUT, Sentinel-2 spectral convolution, scene-mean geometry, fixed TCWV=2.0 cm and sea-level elevation. Sweep the 68-node ACIX-III AOD axis from 0.01 to 4.0."},
                {"stage": "Surface cost", "detail": "At every AOD node, LUT-correct B02/B03/B04 TOA to BOA and sum absolute diagonal chi-square against predicted surface BOA; predicted uncertainty is floored at 0.006."},
                {"stage": "Spatial solve", "detail": "Median-pool each AOD-node cost in a centred 20x20-pixel window (60 m grid, minimum 80 finite samples), add the Gaussian MAIAC penalty, require every node finite, and take the per-pixel argmin."},
                {"stage": "Output", "detail": "Scene output is the arithmetic mean of every finite solved AOD pixel over the retrieval AOI. AOD uncertainty is half the node range within minimum total cost +0.5, floored at 0.02."},
            ],
            "parameters": {
                "surface_source": "Planetary Computer Sentinel-2 L2A only",
                "surface_history": "target year -5 through -1; target month +/-1",
                "surface_windows": 15,
                "surface_top_k": 15,
                "surface_low_aod_fraction": 0.6,
                "surface_max_cloud_cover_percent": 90.0,
                "predictor": "per-realization ExtraTreesRegressor",
                "predictor_trees": 20,
                "predictor_min_samples_leaf": 5,
                "predictor_random_state": 0,
                "predictor_anchor_bands": ["B8A", "B11", "B12"],
                "predictor_localizer_bands": ["B01", "B02", "B03", "B04"],
                "predicted_and_solved_bands": ["B02", "B03", "B04"],
                "surface_uncertainty_floor": 0.006,
                "solver_resolution_m": 60.0,
                "pool_window_pixels": 20,
                "pool_min_finite": 80,
                "aod_axis_nodes": 68,
                "aod_axis_definition": "arange(0.01,0.2,0.01) + arange(0.2,0.5,0.025) + arange(0.5,1.5,0.05) + arange(1.5,2.6,0.1) + arange(2.75,4.01,0.25)",
                "reference_tcwv_cm": 2.0,
                "reference_elevation_km": 0.0,
                "cost_mode": "absolute diagonal chi-square",
                "retrieval_extraction": "scene mean over finite solved AOI pixels",
                "aerosol_species_current": "libRadtran continental average; no per-case family routing",
            },
            "implementation_files": implementation_files,
            "conflict_rule": "Use loose MAIAC uncertainty only when (surface-only AOD minimum - MAIAC AOD) / calibrated MAIAC sigma > 2.576; otherwise retain current calibrated uncertainty.",
            "acceptance": "Freeze only if at least 32/36 medium-development cases are within EE with no unexplained missing outputs; then and only then run folds 2 and 3.",
        },
        "sources": {
            "cases_csv": "downloads/cases.csv",
            "variants_csv": "downloads/variants.csv",
            "report_json": "data/report.json",
            "mids": f"../../{MIDS.name}",
            "archived_results": f"../../{ARCHIVED.name}/",
            "fresh_results": f"../../{FRESH.name}/",
            "cost_cubes": f"../../{CUBES.name}/",
            "split_assignment": "downloads/split-assignment.csv",
            "implementation_manifest": "downloads/implementation/manifest.json",
            "experiment_inventory": "downloads/experiment-inventory.csv",
        },
    }
    (output / "data/report.json").write_text(
        json.dumps(report, separators=(",", ":")) + "\n"
    )
    _write_csv(
        output / "downloads/cases.csv",
        [
            {
                "matchup_id": row["matchup_id"],
                "site": row["site"],
                "truth": row["truth"],
                "fresh": row["fresh"],
                "fresh_error": row["fresh_error"],
                "fresh_within_ee": row["candidate_by_id"]["fresh"]["within_ee"],
                "conflict": row["candidate_by_id"]["conflict"]["value"],
                "conflict_within_ee": row["candidate_by_id"]["conflict"]["within_ee"],
                "maiac": row["maiac"],
                "surface_min": row["surface_min"],
                "B02_min": row["candidate_by_id"]["B02_min"]["value"],
                "B03_min": row["candidate_by_id"]["B03_min"]["value"],
                "B04_min": row["candidate_by_id"]["B04_min"]["value"],
                "transition": row["transition"],
                "cloud_fraction": row["cloud_fraction"],
                "valid_fraction": row["valid_fraction"],
            }
            for row in cases
        ],
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
                "within_ee_percent": row["metrics"].get("within_ee_percent"),
                "bias": row["metrics"]["bias"],
                "mae": row["metrics"]["mae"],
                "rmse": row["metrics"]["rmse"],
                "gains": row["transitions_vs_fresh"].get("gain", 0),
                "losses": row["transitions_vs_fresh"].get("loss", 0),
                "description": row["description"],
                "path": row["path"],
            }
            for row in variant_rows
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
                "within_ee_percent": row["metrics"].get("within_ee_percent"),
                "bias": row["metrics"]["bias"],
                "mae": row["metrics"]["mae"],
                "rmse": row["metrics"]["rmse"],
                "gains": row["transitions_vs_fresh"].get("gain", 0),
                "losses": row["transitions_vs_fresh"].get("loss", 0),
                "status_counts": json.dumps(row["status_counts"], sort_keys=True),
                "description": row["description"],
                "path": row["path"],
            }
            for row in inventory_rows
        ],
    )
    _write_csv(
        output / "downloads/split-assignment.csv",
        [
            {
                "matchup_id": matchup_id,
                "site": archived_all[matchup_id]["site"],
                "truth": archived_all[matchup_id]["truth"],
                "fold": folds[matchup_id],
                "partition": "holdout" if matchup_id in holdout_ids else "development",
                "regime": (
                    "low"
                    if float(archived_all[matchup_id]["truth"]) < 0.25
                    else "medium"
                    if float(archived_all[matchup_id]["truth"]) <= 0.85
                    else "high"
                ),
                "excluded_extreme": matchup_id in EXCLUDED,
            }
            for matchup_id in mids
        ],
    )
    (output / "index.html").write_text(
        """<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width,initial-scale=1">
  <meta name="color-scheme" content="light">
  <title>Medium-AOD physical investigation</title>
  <link rel="icon" href="data:,">
  <link rel="stylesheet" href="app.css">
</head>
<body>
  <div id="app"><div class="loading">Loading the locked development evidence...</div></div>
  <noscript>This report requires JavaScript. CSV and JSON downloads remain available.</noscript>
  <script src="app.js"></script>
</body>
</html>
"""
    )
    receipt = {
        "output": str(output),
        "cases": len(cases),
        "variants": len(variant_rows),
        "spatial_images": len(list((output / "assets/spatial").glob("*.png"))),
        "diagnostic_images": len(list((output / "assets/diagnostic").glob("*.png"))),
        "surface_images": len(list((output / "assets/surface").glob("*.png"))),
        "best_variant": best["label"],
        "best_hits": best["metrics"]["within_ee"],
        "target_met": report["target"]["met"],
        "holdout_scored": False,
    }
    (output / "build-receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
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
