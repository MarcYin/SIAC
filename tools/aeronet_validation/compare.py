"""Stage 4: score retrieved AOT against AERONET AOD at 550 nm per approach.

Joins every ``runs/<approach>/<matchup_id>/result.json`` with the matchup
table, applies quality filters, and reports per-approach statistics plus
scatter plots. The expected-error envelope follows the common satellite-AOD
convention EE = +/-(0.05 + 0.15 * AOD_AERONET).
"""

from __future__ import annotations

import argparse
import logging
import math
from pathlib import Path

import numpy as np
import pandas as pd
from tools.aeronet_validation.common import APPROACHES, ExperimentPaths, read_json

logger = logging.getLogger("aeronet_validation.compare")

EE_OFFSET = 0.05
EE_SLOPE = 0.15


def collect_results(paths: ExperimentPaths) -> pd.DataFrame:
    records = []
    for result_path in sorted(paths.runs_dir.glob("*/*/result.json")):
        record = read_json(result_path)
        record.setdefault("approach", result_path.parent.parent.name)
        record.setdefault("matchup_id", result_path.parent.name)
        records.append(record)
    if not records:
        raise SystemExit(f"No run results under {paths.runs_dir}; run the run stage first.")
    results = pd.DataFrame(records)
    matchups = pd.read_csv(paths.matchups_file)
    merged = results.merge(
        matchups[
            [
                "matchup_id",
                "mgrs_tile",
                "scene_cloud_cover",
                "n_aeronet",
                "mean_abs_time_offset_min",
                "aeronet_aod550_mean",
                "aeronet_aod550_std",
                "aeronet_angstrom_mean",
                "elevation_m",
            ]
        ],
        on="matchup_id",
        how="left",
    )
    return merged


def apply_quality_filters(
    results: pd.DataFrame,
    *,
    min_valid_fraction: float,
    max_cloud_fraction: float,
    aot_field: str,
) -> pd.DataFrame:
    candidates = results[results["status"] == "ok"].copy()
    candidates = candidates.dropna(subset=[aot_field, "aeronet_aod550_mean"])
    n_valid = candidates["aot_window_n_valid"].astype(float)
    n_total = candidates["aot_window_n_total"].astype(float).clip(lower=1.0)
    candidates["aot_valid_fraction"] = n_valid / n_total
    kept = candidates[
        (candidates["aot_valid_fraction"] >= min_valid_fraction)
        & (candidates["aoi_cloud_fraction"].fillna(1.0) <= max_cloud_fraction)
    ]
    logger.info(
        "Quality filter: %d/%d successful runs kept (valid_fraction>=%.2f, cloud<=%.2f)",
        len(kept),
        len(candidates),
        min_valid_fraction,
        max_cloud_fraction,
    )
    return kept


def score_approach(frame: pd.DataFrame, aot_field: str) -> dict[str, float]:
    retrieved = frame[aot_field].astype(float).to_numpy()
    reference = frame["aeronet_aod550_mean"].astype(float).to_numpy()
    error = retrieved - reference
    envelope = EE_OFFSET + EE_SLOPE * reference
    n = len(frame)
    metrics = {
        "n": float(n),
        "bias": float(error.mean()),
        "mae": float(np.abs(error).mean()),
        "rmse": float(math.sqrt((error**2).mean())),
        "median_abs_error": float(np.median(np.abs(error))),
        "within_ee_fraction": float((np.abs(error) <= envelope).mean()),
        "aeronet_mean_aod550": float(reference.mean()),
        "retrieved_mean_aod550": float(retrieved.mean()),
    }
    if n >= 3 and reference.std() > 0 and retrieved.std() > 0:
        metrics["pearson_r"] = float(np.corrcoef(reference, retrieved)[0, 1])
        slope, intercept = np.polyfit(reference, retrieved, 1)
        metrics["slope"] = float(slope)
        metrics["intercept"] = float(intercept)
    return metrics


def plot_scatter(frame: pd.DataFrame, approach: str, aot_field: str, out_path: Path) -> bool:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        logger.warning("matplotlib unavailable; skipping scatter plot for %s", approach)
        return False

    reference = frame["aeronet_aod550_mean"].astype(float).to_numpy()
    retrieved = frame[aot_field].astype(float).to_numpy()
    limit = max(0.5, float(np.nanpercentile(np.concatenate([reference, retrieved]), 99.5)) * 1.1)
    line = np.linspace(0.0, limit, 100)

    figure, axes = plt.subplots(figsize=(5.5, 5.5))
    axes.plot(line, line, "k-", linewidth=1, label="1:1")
    axes.plot(line, line + (EE_OFFSET + EE_SLOPE * line), "k--", linewidth=0.8, label="EE")
    axes.plot(line, line - (EE_OFFSET + EE_SLOPE * line), "k--", linewidth=0.8)
    axes.scatter(reference, retrieved, s=12, alpha=0.5, edgecolors="none")
    metrics = score_approach(frame, aot_field)
    text = (
        f"N = {int(metrics['n'])}\n"
        f"bias = {metrics['bias']:.3f}\n"
        f"RMSE = {metrics['rmse']:.3f}\n"
        f"R = {metrics.get('pearson_r', float('nan')):.3f}\n"
        f"within EE = {100 * metrics['within_ee_fraction']:.1f}%"
    )
    axes.text(
        0.03,
        0.97,
        text,
        transform=axes.transAxes,
        va="top",
        ha="left",
        fontsize=9,
        bbox={"boxstyle": "round", "facecolor": "white", "alpha": 0.8},
    )
    axes.set_xlim(0, limit)
    axes.set_ylim(0, limit)
    axes.set_xlabel("AERONET AOD 550 nm")
    axes.set_ylabel("SIAC retrieved AOT 550 nm")
    axes.set_title(f"Surface prior: {approach}")
    axes.legend(loc="lower right", fontsize=8)
    figure.tight_layout()
    figure.savefig(out_path, dpi=150)
    plt.close(figure)
    return True


def add_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--aot-field",
        default="aot_window_mean",
        choices=("aot_window_mean", "aot_window_median", "aot_nearest"),
        help="Which extracted AOT statistic to validate.",
    )
    parser.add_argument("--min-valid-fraction", type=float, default=0.5)
    parser.add_argument("--max-cloud-fraction", type=float, default=0.4)


def run(args: argparse.Namespace, paths: ExperimentPaths) -> None:
    paths.results_dir.mkdir(parents=True, exist_ok=True)
    results = collect_results(paths)
    results.to_csv(paths.results_dir / "matchup_results_raw.csv", index=False)

    kept = apply_quality_filters(
        results,
        min_valid_fraction=args.min_valid_fraction,
        max_cloud_fraction=args.max_cloud_fraction,
        aot_field=args.aot_field,
    )
    kept.to_csv(paths.results_dir / "matchup_results_filtered.csv", index=False)

    summary_rows = []
    for approach in sorted(APPROACHES):
        frame = kept[kept["approach"] == approach]
        if frame.empty:
            logger.info("No filtered results for approach %s", approach)
            continue
        metrics = score_approach(frame, args.aot_field)
        metrics["approach"] = approach
        summary_rows.append(metrics)
        plot_scatter(frame, approach, args.aot_field, paths.results_dir / f"scatter_{approach}.png")

    if not summary_rows:
        raise SystemExit("No approaches with filtered results to score.")
    summary = pd.DataFrame(summary_rows).set_index("approach")
    ordered = [
        "n",
        "bias",
        "mae",
        "rmse",
        "median_abs_error",
        "within_ee_fraction",
        "pearson_r",
        "slope",
        "intercept",
        "aeronet_mean_aod550",
        "retrieved_mean_aod550",
    ]
    summary = summary[[c for c in ordered if c in summary.columns]]
    summary.to_csv(paths.results_dir / "summary.csv")
    logger.info("Summary written to %s", paths.results_dir / "summary.csv")
    print()
    print(summary.round(4).to_string())
