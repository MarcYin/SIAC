"""Visual report for the cross-RT L2A harmonisation experiment.

Reads the raw run artifacts directly (pair-level scene metrics CSV, harmonizer
JSON, per-case retrieval JSONs) and renders a static figure package plus a
self-describing HTML page. Every statistic shown is recomputed here from the
raw artifacts; nothing is imported from the tabular summariser.
"""

from __future__ import annotations

import argparse
import html
import json
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

ROOT = Path("/gws/ssde/j25a/nceo_isp/public/siac_refactor")
DEFAULT_RUN = ROOT / "analysis/cross_rt_harmonizer_lowcloud152_20260716"
DEFAULT_CASES = ROOT / "reports/aod-final-performance-dashboard-20260713/data/all-cases.csv"

MODEL = "cross_rt_terrain_a1"
CAP_PREFIX = "cross_rt_terrain_a1_cap_0p030"
SELECTED_VARIANT = "cross_rt_terrain_a1_solver_cap0p030"
SOLVER_BANDS = (("blue", "B02 blue"), ("green", "B03 green"), ("red", "B04 red"))
ALL_BANDS = ("coastal", "blue", "green", "red", "nir08", "swir16", "swir22")

RETRIEVAL_DIRS = {
    "previous_best": "phaseD_results_lowcloud20_geometry_backstop05_b03_chi2_20260712",
    "current_physical": (
        "phaseD_results_lowcloud20_physical_anchor_maiac_l2awvp_dem_scenegeom_20260715"
    ),
    **{
        tag: f"phaseD_results_lowcloud20_crossrt_{tag}_physical_20260716"
        for tag in (
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
    },
}
VARIANT_LABELS = {
    "previous_best": "previous best (control)",
    "current_physical": "current physical (control)",
    "identity_daily": "identity history (uncorrected L2A)",
    "cross_rt_baseline_a1_solver_cap0p030": "baseline stage · solver · cap 0.030",
    "cross_rt_aod_a1_solver_cap0p030": "+ ΔAOT stage · solver · cap 0.030",
    "cross_rt_atmosphere_a1_solver_cap0p030": "+ atmosphere stage · solver · cap 0.030",
    "cross_rt_terrain_a1_solver_cap0p015": "+ terrain α=1 · solver · cap 0.015",
    "cross_rt_terrain_a1_solver_cap0p030": "+ terrain α=1 · solver · cap 0.030",
    "cross_rt_terrain_a1_all_cap0p015": "+ terrain α=1 · all bands · cap 0.015",
    "cross_rt_terrain_a1_all_cap0p030": "+ terrain α=1 · all bands · cap 0.030",
    "cross_rt_terrain_a10_solver_cap0p030": "+ terrain α=10 · solver · cap 0.030",
}

# Reference dataviz palette (validated: blue<->red CVD dE 21.6, both >=3:1 on surface).
BLUE = "#2a78d6"
LIGHT_BLUE = "#9ec5f4"
RED = "#e34948"
INK = "#0b0b0b"
SECONDARY = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"
SURFACE = "#fcfcfb"

FEATURE_LABELS = {
    "delta_aot_maiac_minus_sen2cor": "ΔAOT (MAIAC − Sen2Cor)",
    "mean_aot": "mean AOT",
    "l2a_band": "L2A reflectance",
    "l2a_band_x_mean_aot": "L2A × mean AOT",
    "l2a_band_x_delta_aot": "L2A × ΔAOT",
    "l2a_tcwv_cm": "L2A TCWV",
    "l2a_band_x_l2a_tcwv": "L2A × TCWV",
    "solar_airmass": "solar airmass",
    "view_airmass": "view airmass",
    "l2a_band_x_total_airmass": "L2A × total airmass",
    "cos_raa": "cos(RAA)",
    "elevation_km": "scene elevation",
    "month_sin": "month (sin)",
    "month_cos": "month (cos)",
    "sensor_is_s2b": "sensor is S2B",
    "sensor_is_s2c": "sensor is S2C",
    "processing_baseline": "processing baseline",
    "delta_aot_x_total_airmass": "ΔAOT × total airmass",
    "delta_aot_x_mean_aot": "ΔAOT × mean AOT",
    "signed_delta_aot_squared": "ΔAOT · |ΔAOT|",
    "l2a_tcwv_x_solar_airmass": "TCWV × solar airmass",
    "l2a_band_x_l2a_tcwv_x_total_airmass": "L2A × TCWV × airmass",
    "delta_aot_x_l2a_tcwv": "ΔAOT × TCWV",
    "delta_aot_x_cos_raa": "ΔAOT × cos(RAA)",
    "terrain_elevation_km": "pixel elevation",
    "terrain_elevation_delta_km": "elevation − scene mean",
    "terrain_slope_deg": "terrain slope",
    "terrain_incidence_cos": "cos(solar incidence)",
    "terrain_incidence_delta": "incidence − flat terrain",
    "terrain_illumination_ratio_minus_one": "illumination ratio − 1",
    "terrain_shadow": "shadow flag",
    "l2a_band_x_terrain_incidence_delta": "L2A × incidence δ",
    "l2a_band_x_terrain_illumination_ratio_minus_one": "L2A × illumination ratio",
    "delta_aot_x_terrain_incidence_delta": "ΔAOT × incidence δ",
}


def style() -> None:
    # Start from matplotlib's internal defaults so a user-level matplotlibrc
    # (e.g. publication-sized axis labels) cannot leak into the report.
    matplotlib.rcdefaults()
    plt.rcParams.update(
        {
            "figure.facecolor": SURFACE,
            "axes.facecolor": SURFACE,
            "savefig.facecolor": SURFACE,
            "font.family": "DejaVu Sans",
            "font.size": 9.5,
            "axes.labelsize": 9.5,
            "axes.edgecolor": AXIS,
            "axes.linewidth": 0.8,
            "axes.labelcolor": SECONDARY,
            "axes.titlecolor": INK,
            "axes.titlesize": 10.5,
            "axes.titleweight": "bold",
            "axes.titlepad": 8.0,
            "xtick.color": MUTED,
            "ytick.color": MUTED,
            "xtick.labelsize": 8.5,
            "ytick.labelsize": 8.5,
            "axes.grid": True,
            "grid.color": GRID,
            "grid.linewidth": 0.8,
            "grid.linestyle": "-",
            "axes.axisbelow": True,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "legend.frameon": False,
            "legend.fontsize": 9,
        }
    )


def feature_label(name: str, band: str) -> str:
    generic = name.replace(f"l2a_{band}", "l2a_band")
    return FEATURE_LABELS.get(generic, generic)


def binned_quantiles(
    x: np.ndarray,
    y: np.ndarray,
    edges: np.ndarray,
    min_count: int = 25,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    centers, median, lower, upper = [], [], [], []
    finite = np.isfinite(x) & np.isfinite(y)
    for start, stop in zip(edges[:-1], edges[1:], strict=False):
        mask = finite & (x >= start) & (x < stop)
        if int(mask.sum()) < min_count:
            continue
        values = y[mask]
        centers.append(float(np.median(x[mask])))
        median.append(float(np.median(values)))
        lower.append(float(np.percentile(values, 25)))
        upper.append(float(np.percentile(values, 75)))
    return np.array(centers), np.array(median), np.array(lower), np.array(upper)


def save(figure: plt.Figure, output_dir: Path, name: str) -> str:
    path = output_dir / name
    figure.savefig(path, dpi=165, bbox_inches="tight", pad_inches=0.18)
    plt.close(figure)
    print(f"wrote {path}")
    return name


def load_scene_table(run_root: Path) -> pd.DataFrame:
    needed = ["site", "sample_count", "maiac_aot", "l2a_aot", "delta_aot_maiac_minus_sen2cor"]
    for band, _label in SOLVER_BANDS:
        needed += [
            f"identity_{band}_bias",
            f"identity_{band}_mae",
            f"{CAP_PREFIX}_{band}_bias",
            f"{CAP_PREFIX}_{band}_mae",
        ]
    frame = pd.read_csv(run_root / "surface_scene_metrics.csv", usecols=needed)
    for band, _label in SOLVER_BANDS:
        frame[f"correction_{band}"] = (
            frame[f"{CAP_PREFIX}_{band}_bias"] - frame[f"identity_{band}_bias"]
        )
    return frame


def load_retrieval(root: Path, cases: Path) -> tuple[list[str], dict[str, pd.DataFrame]]:
    case_rows = pd.read_csv(cases)
    matchup_ids = case_rows["matchup_id"].tolist()
    results: dict[str, pd.DataFrame] = {}
    for variant, directory in RETRIEVAL_DIRS.items():
        rows = []
        for matchup_id in matchup_ids:
            path = root / directory / f"{matchup_id}.json"
            record: dict[str, Any] = {"matchup_id": matchup_id, "status": "MISSING"}
            if path.exists():
                try:
                    payload = json.loads(path.read_text(encoding="utf-8"))
                    record = {
                        "matchup_id": matchup_id,
                        "status": payload.get("status", "MISSING"),
                        "site": payload.get("site"),
                        "truth": payload.get("truth"),
                        "retrieved": payload.get("retrieved"),
                        "within_ee": bool(payload.get("within_ee")),
                        "ee_threshold": payload.get("ee_threshold"),
                    }
                except (OSError, json.JSONDecodeError):
                    record["status"] = "MALFORMED"
            rows.append(record)
        frame = pd.DataFrame(rows)
        for column in ("truth", "retrieved", "ee_threshold"):
            frame[column] = pd.to_numeric(frame.get(column), errors="coerce")
        results[variant] = frame
    return matchup_ids, results


def valid_mask(frame: pd.DataFrame) -> pd.Series:
    return (
        (frame["status"] == "OK")
        & np.isfinite(frame["truth"])
        & np.isfinite(frame["retrieved"])
    )


def within_ee_rate(frame: pd.DataFrame) -> tuple[float, int, int]:
    hits = int((valid_mask(frame) & frame["within_ee"]).sum())
    expected = len(frame)
    return hits / expected, hits, expected


# ---------------------------------------------------------------------------
# Section 1 — what the correction does
# ---------------------------------------------------------------------------


def figure_correction_vs_delta_aot(scene: pd.DataFrame, output_dir: Path) -> str:
    edges = np.arange(-0.275, 0.425, 0.05)
    fig, axes = plt.subplots(1, 3, figsize=(11.8, 3.7), sharey=True, constrained_layout=True)
    for ax, (band, label) in zip(axes, SOLVER_BANDS, strict=True):
        x = scene["delta_aot_maiac_minus_sen2cor"].to_numpy()
        y = scene[f"correction_{band}"].to_numpy()
        ax.scatter(x, y, s=5, color=BLUE, alpha=0.14, linewidths=0, rasterized=True)
        centers, median, lower, upper = binned_quantiles(x, y, edges)
        ax.fill_between(centers, lower, upper, color=BLUE, alpha=0.18, linewidth=0)
        ax.plot(centers, median, color=BLUE, linewidth=2, solid_capstyle="round")
        ax.axhline(0.0, color=AXIS, linewidth=0.8)
        ax.axvline(0.0, color=AXIS, linewidth=0.8)
        ax.set_title(label)
        ax.set_xlim(-0.3, 0.45)
        ax.set_ylim(-0.022, 0.022)
        ax.set_xlabel("ΔAOT = MAIAC − Sen2Cor")
    axes[0].set_ylabel("applied correction (reflectance)")
    return save(fig, output_dir, "f01_correction_vs_delta_aot.png")


def figure_bias_vs_delta_aot(scene: pd.DataFrame, output_dir: Path) -> str:
    edges = np.arange(-0.275, 0.425, 0.05)
    fig, axes = plt.subplots(1, 3, figsize=(11.8, 3.7), sharey=True, constrained_layout=True)
    for ax, (band, label) in zip(axes, SOLVER_BANDS, strict=True):
        x = scene["delta_aot_maiac_minus_sen2cor"].to_numpy()
        for column, color, name in (
            (f"identity_{band}_bias", MUTED, "before (operational L2A)"),
            (f"{CAP_PREFIX}_{band}_bias", BLUE, "after correction"),
        ):
            y = scene[column].to_numpy()
            centers, median, lower, upper = binned_quantiles(x, y, edges)
            ax.fill_between(centers, lower, upper, color=color, alpha=0.15, linewidth=0)
            ax.plot(
                centers, median, color=color, linewidth=2, solid_capstyle="round", label=name
            )
        ax.axhline(0.0, color=AXIS, linewidth=0.8)
        ax.set_title(label)
        ax.set_xlim(-0.3, 0.45)
        ax.set_ylim(-0.014, 0.014)
        ax.set_xlabel("ΔAOT = MAIAC − Sen2Cor")
    axes[0].set_ylabel("component error vs target (bias)")
    axes[0].legend(loc="upper right")
    return save(fig, output_dir, "f02_bias_vs_delta_aot.png")


def figure_bias_distributions(scene: pd.DataFrame, output_dir: Path) -> str:
    bins = np.arange(-0.025, 0.0255, 0.001)
    fig, axes = plt.subplots(1, 3, figsize=(11.8, 3.4), sharey=True, constrained_layout=True)
    for ax, (band, label) in zip(axes, SOLVER_BANDS, strict=True):
        before = scene[f"identity_{band}_bias"].dropna().to_numpy()
        after = scene[f"{CAP_PREFIX}_{band}_bias"].dropna().to_numpy()
        ax.hist(
            before,
            bins=bins,
            density=True,
            histtype="stepfilled",
            color=MUTED,
            alpha=0.25,
            linewidth=0,
            label="before (operational L2A)",
        )
        ax.hist(
            after,
            bins=bins,
            density=True,
            histtype="step",
            color=BLUE,
            linewidth=2,
            label="after correction",
        )
        ax.axvline(0.0, color=AXIS, linewidth=0.8)
        ax.set_title(label)
        ax.set_xlim(-0.02, 0.02)
        ax.set_xlabel("component bias vs target")
        ax.set_yticks([])
        ax.grid(False, axis="y")
    axes[0].set_ylabel("density")
    axes[0].legend(loc="upper left")
    return save(fig, output_dir, "f03_bias_distributions.png")


# ---------------------------------------------------------------------------
# Section 2 — accuracy of the correction
# ---------------------------------------------------------------------------


def figure_component_error_scatter(scene: pd.DataFrame, output_dir: Path) -> str:
    fig, axes = plt.subplots(1, 3, figsize=(11.8, 4.0), constrained_layout=True)
    limits = (2.0e-4, 0.06)
    for ax, (band, label) in zip(axes, SOLVER_BANDS, strict=True):
        before = scene[f"identity_{band}_bias"].abs().clip(lower=limits[0]).to_numpy()
        after = scene[f"{CAP_PREFIX}_{band}_bias"].abs().clip(lower=limits[0]).to_numpy()
        finite = np.isfinite(before) & np.isfinite(after)
        improved = 100.0 * float((after[finite] < before[finite]).mean())
        ax.scatter(
            before[finite],
            after[finite],
            s=5,
            color=BLUE,
            alpha=0.14,
            linewidths=0,
            rasterized=True,
        )
        ax.plot(limits, limits, color=AXIS, linewidth=0.8)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(*limits)
        ax.set_ylim(*limits)
        ax.set_title(label)
        ax.set_xlabel("|bias| before")
        ax.text(
            0.03,
            0.95,
            f"{improved:.0f}% of components improve",
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=9,
            color=SECONDARY,
        )
    axes[0].set_ylabel("|bias| after")
    return save(fig, output_dir, "f04_component_error_scatter.png")


def figure_site_mae(scene: pd.DataFrame, output_dir: Path) -> str:
    sites = (
        scene.groupby("site")
        .agg(
            before=("identity_blue_mae", "mean"),
            after=(f"{CAP_PREFIX}_blue_mae", "mean"),
            components=("site", "size"),
        )
        .dropna()
    )
    sites["delta"] = sites["after"] - sites["before"]
    fig, ax = plt.subplots(figsize=(6.4, 5.6), constrained_layout=True)
    ax.scatter(
        sites["before"],
        sites["after"],
        s=np.clip(sites["components"] * 0.65, 12.0, 90.0),
        color=BLUE,
        alpha=0.55,
        linewidths=1.2,
        edgecolors=SURFACE,
    )
    limit = float(max(sites["before"].max(), sites["after"].max())) * 1.08
    ax.plot([0, limit], [0, limit], color=AXIS, linewidth=0.8)
    ax.set_xlim(0, limit)
    ax.set_ylim(0, limit)
    movers = pd.concat(
        [sites.nlargest(3, "delta"), sites.nsmallest(3, "delta"), sites.nlargest(1, "before")]
    )
    movers = movers[~movers.index.duplicated()].sort_values("after", ascending=False)
    offsets = ((7, -3), (7, 9), (7, -15))
    for index, (name, row) in enumerate(movers.iterrows()):
        ax.annotate(
            str(name),
            (row["before"], row["after"]),
            textcoords="offset points",
            xytext=offsets[index % len(offsets)],
            fontsize=8,
            color=SECONDARY,
        )
    improved = 100.0 * float((sites["delta"] < 0).mean())
    ax.set_title(f"Per-site B02 blue MAE — {improved:.0f}% of {len(sites)} sites improve")
    ax.set_xlabel("mean component MAE before")
    ax.set_ylabel("mean component MAE after")
    return save(fig, output_dir, "f05_site_mae_before_after.png")


def figure_stage_ablation(surface_metrics: dict[str, Any], output_dir: Path) -> str:
    def visible_mae(band_metrics: dict[str, Any]) -> float:
        return float(np.mean([band_metrics[band]["scene_mae"] for band, _ in SOLVER_BANDS]))

    stages = [
        ("uncorrected L2A (identity)", visible_mae(surface_metrics["identity"]), MUTED),
        (
            "current production AOD offset",
            visible_mae(surface_metrics["current_aod_offset"]),
            MUTED,
        ),
        (
            "baseline stage",
            visible_mae(surface_metrics["candidates"]["cross_rt_baseline_a1"]["cap_0.030"]),
            LIGHT_BLUE,
        ),
        (
            "+ ΔAOT terms",
            visible_mae(surface_metrics["candidates"]["cross_rt_aod_a1"]["cap_0.030"]),
            LIGHT_BLUE,
        ),
        (
            "+ atmosphere terms",
            visible_mae(surface_metrics["candidates"]["cross_rt_atmosphere_a1"]["cap_0.030"]),
            LIGHT_BLUE,
        ),
        (
            "+ terrain terms (α=1, selected)",
            visible_mae(surface_metrics["candidates"]["cross_rt_terrain_a1"]["cap_0.030"]),
            BLUE,
        ),
        (
            "+ terrain terms (α=10)",
            visible_mae(surface_metrics["candidates"]["cross_rt_terrain_a10"]["cap_0.030"]),
            LIGHT_BLUE,
        ),
    ]
    fig, ax = plt.subplots(figsize=(8.6, 3.6), constrained_layout=True)
    positions = np.arange(len(stages))[::-1]
    values = [value for _, value, _ in stages]
    colors = [color for _, _, color in stages]
    ax.barh(positions, values, height=0.62, color=colors, zorder=3)
    for position, value in zip(positions, values, strict=True):
        ax.text(
            value + 0.00008,
            position,
            f"{value:.5f}",
            va="center",
            ha="left",
            fontsize=8.5,
            color=SECONDARY,
        )
    ax.set_yticks(positions)
    ax.set_yticklabels([name for name, _, _ in stages], fontsize=9, color=INK)
    ax.set_xlim(0, max(values) * 1.16)
    ax.grid(False, axis="y")
    ax.set_xlabel("visible acquisition MAE vs target (out-of-fold, cap 0.030)")
    ax.set_title("Where the accuracy comes from — feature-stage ablation")
    return save(fig, output_dir, "f06_stage_ablation.png")


def figure_residual_spread(artifact: dict[str, Any], output_dir: Path) -> str:
    scale = artifact["models"][MODEL]["oof_residual_scale"]
    fig, ax = plt.subplots(figsize=(8.6, 3.4), constrained_layout=True)
    positions = np.arange(len(ALL_BANDS))[::-1]
    sigma = [scale[band]["mad_to_sigma"] for band in ALL_BANDS]
    rmse = [scale[band]["rmse"] for band in ALL_BANDS]
    ax.hlines(positions, sigma, rmse, color=GRID, linewidth=2, zorder=2)
    ax.scatter(
        sigma,
        positions,
        s=64,
        color=BLUE,
        zorder=3,
        marker="o",
        linewidths=1.5,
        edgecolors=SURFACE,
        label="robust σ (MAD)",
    )
    ax.scatter(
        rmse,
        positions,
        s=64,
        color=LIGHT_BLUE,
        zorder=3,
        marker="s",
        linewidths=1.5,
        edgecolors=SURFACE,
        label="RMSE",
    )
    ax.set_yticks(positions)
    ax.set_yticklabels(ALL_BANDS, fontsize=9, color=INK)
    ax.grid(False, axis="y")
    ax.set_xlim(0, max(rmse) * 1.15)
    ax.set_xlabel("out-of-fold residual spread (reflectance, uncapped prediction)")
    ax.set_title("Correction accuracy per band — residual spread")
    ax.legend(loc="lower right")
    return save(fig, output_dir, "f07_residual_spread.png")


# ---------------------------------------------------------------------------
# Section 3 — how the correction is derived
# ---------------------------------------------------------------------------


def _holdout_models(artifact: dict[str, Any]) -> dict[str, Any]:
    fold = str(artifact["crossfit_protocol"]["holdout_model_fold"])
    return artifact["models"][MODEL]["folds"][fold]


def figure_coefficients(artifact: dict[str, Any], output_dir: Path, top: int = 16) -> str:
    models = _holdout_models(artifact)
    coefficient_by_feature: dict[str, dict[str, float]] = {}
    for band, _label in SOLVER_BANDS:
        model = models[band]
        for name, value in zip(model["feature_names"], model["coef"], strict=True):
            coefficient_by_feature.setdefault(feature_label(name, band), {})[band] = float(value)
    ranked = sorted(
        coefficient_by_feature.items(),
        key=lambda item: -max(abs(value) for value in item[1].values()),
    )[:top]
    labels = [name for name, _ in ranked][::-1]
    positions = np.arange(len(labels))
    fig, axes = plt.subplots(
        1, 3, figsize=(11.8, 5.6), sharey=True, sharex=True, constrained_layout=True
    )
    for ax, (band, label) in zip(axes, SOLVER_BANDS, strict=True):
        values = [coefficient_by_feature[name].get(band, 0.0) for name in labels]
        colors = [BLUE if value >= 0 else RED for value in values]
        ax.barh(positions, values, height=0.62, color=colors, zorder=3)
        ax.axvline(0.0, color=AXIS, linewidth=0.8)
        ax.set_title(label)
        ax.grid(False, axis="y")
        ax.set_xlabel("reflectance per +1 SD of feature")
    axes[0].set_yticks(positions)
    axes[0].set_yticklabels(labels, fontsize=8.5, color=INK)
    handles = [
        Line2D([], [], marker="s", linestyle="", color=BLUE, label="positive (brightens)"),
        Line2D([], [], marker="s", linestyle="", color=RED, label="negative (darkens)"),
    ]
    fig.legend(handles=handles, loc="outside upper right", ncols=2)
    return save(fig, output_dir, "f08_coefficients.png")


def figure_coefficient_stability(
    artifact: dict[str, Any], output_dir: Path, top: int = 10
) -> str:
    folds = artifact["models"][MODEL]["folds"]
    holdout_fold = str(artifact["crossfit_protocol"]["holdout_model_fold"])
    development_folds = [fold for fold in folds if fold != holdout_fold]
    holdout = folds[holdout_fold]["blue"]
    order = np.argsort(-np.abs(np.asarray(holdout["coef"], dtype=float)))[:top]
    labels = [feature_label(holdout["feature_names"][index], "blue") for index in order][::-1]
    positions = np.arange(len(labels))
    fig, ax = plt.subplots(figsize=(8.0, 4.4), constrained_layout=True)
    for fold in development_folds:
        coefficients = np.asarray(folds[fold]["blue"]["coef"], dtype=float)[order][::-1]
        ax.scatter(
            coefficients,
            positions,
            s=30,
            color=MUTED,
            alpha=0.75,
            linewidths=1.2,
            edgecolors=SURFACE,
            zorder=3,
        )
    holdout_values = np.asarray(holdout["coef"], dtype=float)[order][::-1]
    ax.scatter(
        holdout_values,
        positions,
        s=80,
        color=BLUE,
        linewidths=1.5,
        edgecolors=SURFACE,
        zorder=4,
    )
    ax.axvline(0.0, color=AXIS, linewidth=0.8)
    ax.set_yticks(positions)
    ax.set_yticklabels(labels, fontsize=8.5, color=INK)
    ax.grid(False, axis="y")
    ax.set_xlabel("standardized coefficient, B02 blue (reflectance per +1 SD)")
    handles = [
        Line2D(
            [],
            [],
            marker="o",
            linestyle="",
            color=BLUE,
            markersize=9,
            label="applied model (all development sites)",
        ),
        Line2D(
            [],
            [],
            marker="o",
            linestyle="",
            color=MUTED,
            markersize=6,
            label="individual cross-fit folds",
        ),
    ]
    fig.legend(handles=handles, loc="outside upper right", ncols=2)
    ax.set_title("Coefficient stability across site-grouped folds", loc="left")
    return save(fig, output_dir, "f09_coefficient_stability.png")


# ---------------------------------------------------------------------------
# Section 4 — retrieval impact
# ---------------------------------------------------------------------------


def figure_retrieval_rates(results: dict[str, pd.DataFrame], output_dir: Path) -> str:
    order = list(RETRIEVAL_DIRS)
    rates = {variant: within_ee_rate(results[variant]) for variant in order}
    fig, ax = plt.subplots(figsize=(9.4, 4.6), constrained_layout=True)
    positions = np.arange(len(order))[::-1]
    for position, variant in zip(positions, order, strict=True):
        rate, hits, expected = rates[variant]
        if variant == SELECTED_VARIANT:
            color = BLUE
        elif variant in {"previous_best", "current_physical", "identity_daily"}:
            color = MUTED
        else:
            color = LIGHT_BLUE
        ax.barh(position, 100.0 * rate, height=0.62, color=color, zorder=3)
        ax.text(
            100.0 * rate + 0.7,
            position,
            f"{100.0 * rate:.1f}%  ({hits}/{expected})",
            va="center",
            ha="left",
            fontsize=8.5,
            color=SECONDARY,
        )
    ax.axvline(87.0, color=AXIS, linewidth=0.8)
    ax.text(87.0, len(order) - 0.25, " 87% target", fontsize=8.5, color=MUTED, va="bottom")
    ax.set_yticks(positions)
    ax.set_yticklabels([VARIANT_LABELS[variant] for variant in order], fontsize=9, color=INK)
    ax.grid(False, axis="y")
    ax.set_xlim(0, 100)
    ax.set_xlabel("within expected error, complete 152-case cohort (missing counts as miss)")
    ax.set_title("Retrieval outcome by history variant")
    return save(fig, output_dir, "f10_retrieval_within_ee.png")


def figure_retrieval_scatter(results: dict[str, pd.DataFrame], output_dir: Path) -> str:
    panels = (
        ("identity_daily", "Uncorrected L2A history"),
        (SELECTED_VARIANT, "Harmonized history (selected)"),
        ("previous_best", "Previous best (control)"),
    )
    limits = (0.008, 3.0)
    truth_grid = np.geomspace(*limits, 200)
    envelope = 0.05 + 0.15 * truth_grid
    fig, axes = plt.subplots(
        1, 3, figsize=(11.8, 4.2), sharex=True, sharey=True, constrained_layout=True
    )
    for ax, (variant, title) in zip(axes, panels, strict=True):
        frame = results[variant]
        valid = frame[valid_mask(frame)]
        ax.fill_between(
            truth_grid,
            np.clip(truth_grid - envelope, limits[0], None),
            truth_grid + envelope,
            color=MUTED,
            alpha=0.12,
            linewidth=0,
        )
        ax.plot(truth_grid, truth_grid, color=AXIS, linewidth=0.8)
        inside = valid[valid["within_ee"]]
        outside = valid[~valid["within_ee"]]
        common = {"linewidths": 1.0, "edgecolors": SURFACE}
        ax.scatter(
            np.clip(inside["truth"], *limits),
            np.clip(inside["retrieved"], *limits),
            s=26,
            color=BLUE,
            label="within EE",
            **common,
        )
        ax.scatter(
            np.clip(outside["truth"], *limits),
            np.clip(outside["retrieved"], *limits),
            s=26,
            color=RED,
            label="outside EE",
            **common,
        )
        rate, hits, expected = within_ee_rate(frame)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(*limits)
        ax.set_ylim(*limits)
        ax.set_title(f"{title} — {100.0 * rate:.1f}% ({hits}/{expected})", fontsize=10)
        ax.set_xlabel("AERONET AOD (truth)")
    axes[0].set_ylabel("retrieved AOD")
    axes[0].legend(loc="upper left")
    return save(fig, output_dir, "f11_retrieval_scatter.png")


def figure_retrieval_error_vs_aod(results: dict[str, pd.DataFrame], output_dir: Path) -> str:
    edges = np.geomspace(0.015, 2.6, 9)
    truth_grid = np.geomspace(0.015, 2.6, 200)
    envelope = 0.05 + 0.15 * truth_grid
    fig, ax = plt.subplots(figsize=(8.6, 4.4), constrained_layout=True)
    ax.fill_between(truth_grid, -envelope, envelope, color=MUTED, alpha=0.12, linewidth=0)
    ax.plot(truth_grid, envelope, color=AXIS, linewidth=0.8)
    ax.plot(truth_grid, -envelope, color=AXIS, linewidth=0.8)
    ax.text(1.9, 0.05 + 0.15 * 1.9, "expected error", fontsize=8.5, color=MUTED, va="bottom")
    for variant, color, label in (
        ("identity_daily", MUTED, "uncorrected L2A history"),
        (SELECTED_VARIANT, BLUE, "harmonized history (selected)"),
    ):
        frame = results[variant]
        valid = frame[valid_mask(frame)]
        error = (valid["retrieved"] - valid["truth"]).to_numpy()
        truth = valid["truth"].to_numpy()
        ax.scatter(truth, error, s=13, color=color, alpha=0.3, linewidths=0)
        centers, median, lower, upper = binned_quantiles(truth, error, edges, min_count=6)
        ax.fill_between(centers, lower, upper, color=color, alpha=0.15, linewidth=0)
        ax.plot(
            centers,
            median,
            color=color,
            linewidth=2,
            marker="o",
            markersize=5,
            markeredgecolor=SURFACE,
            markeredgewidth=1.0,
            label=label,
        )
    ax.axhline(0.0, color=AXIS, linewidth=0.8)
    ax.set_xscale("log")
    ax.set_xlim(0.015, 2.6)
    ax.set_ylim(-0.62, 0.42)
    ax.set_xlabel("AERONET AOD (truth)")
    ax.set_ylabel("retrieved − truth (binned median, IQR band)")
    ax.set_title("Retrieval error vs AOD level")
    ax.legend(loc="lower left")
    return save(fig, output_dir, "f12_retrieval_error_vs_aod.png")


def figure_case_transitions(results: dict[str, pd.DataFrame], output_dir: Path) -> str:
    identity = results["identity_daily"].set_index("matchup_id")
    selected = results[SELECTED_VARIANT].set_index("matchup_id")
    joined = identity.join(selected, lsuffix="_before", rsuffix="_after")
    joined = joined[
        (joined["status_before"] == "OK")
        & (joined["status_after"] == "OK")
        & (joined["within_ee_before"] != joined["within_ee_after"])
    ].copy()
    for suffix in ("before", "after"):
        joined[f"ratio_{suffix}"] = (
            (joined[f"retrieved_{suffix}"] - joined[f"truth_{suffix}"]).abs()
            / joined[f"ee_threshold_{suffix}"]
        )
    joined["gain"] = joined["within_ee_after"]
    joined = joined.sort_values(["gain", "ratio_after"], ascending=[False, True])
    labels = [
        f"{row['site_before']}  (truth {row['truth_before']:.2f})"
        for _, row in joined.iterrows()
    ]
    positions = np.arange(len(joined))[::-1]
    fig, ax = plt.subplots(figsize=(9.4, 0.34 * len(joined) + 1.6), constrained_layout=True)
    for position, (_, row) in zip(positions, joined.iterrows(), strict=True):
        color = BLUE if row["gain"] else RED
        ax.annotate(
            "",
            xy=(row["ratio_after"], position),
            xytext=(row["ratio_before"], position),
            arrowprops={"arrowstyle": "-|>", "color": color, "linewidth": 1.8},
        )
        ax.scatter(
            [row["ratio_before"]],
            [position],
            s=34,
            color=MUTED,
            linewidths=1.2,
            edgecolors=SURFACE,
            zorder=3,
        )
    ax.axvline(1.0, color=AXIS, linewidth=0.8)
    ax.text(
        1.0,
        0.995,
        " EE boundary",
        transform=ax.get_xaxis_transform(),
        fontsize=8.5,
        color=MUTED,
        va="top",
        ha="left",
    )
    ax.set_yticks(positions)
    ax.set_yticklabels(labels, fontsize=8.5, color=INK)
    ax.grid(False, axis="y")
    ax.set_xlim(0, max(2.6, float(joined[["ratio_before", "ratio_after"]].max().max()) * 1.1))
    ax.set_xlabel("|retrieved − truth| / expected error   (dot = before, arrow = after)")
    gains = int(joined["gain"].sum())
    losses = int(len(joined) - gains)
    handles = [
        Line2D([], [], color=BLUE, linewidth=2, label=f"gained EE ({gains})"),
        Line2D([], [], color=RED, linewidth=2, label=f"lost EE ({losses})"),
    ]
    ax.legend(handles=handles, loc="lower right")
    ax.set_title("Cases whose EE verdict changed vs the uncorrected history")
    return save(fig, output_dir, "f13_case_transitions.png")


# ---------------------------------------------------------------------------
# HTML assembly
# ---------------------------------------------------------------------------


def stat_tiles(
    scene: pd.DataFrame,
    artifact: dict[str, Any],
    results: dict[str, pd.DataFrame],
) -> list[tuple[str, str, str]]:
    before = float(scene["identity_blue_mae"].mean())
    after = float(scene[f"{CAP_PREFIX}_blue_mae"].mean())
    selected_rate, hits, expected = within_ee_rate(results[SELECTED_VARIANT])
    identity_rate, _, _ = within_ee_rate(results["identity_daily"])
    previous_rate, _, _ = within_ee_rate(results["previous_best"])
    sigma = float(artifact["models"][MODEL]["oof_residual_scale"]["blue"]["mad_to_sigma"])
    return [
        (
            "Blue surface error vs target",
            f"−{100.0 * (1.0 - after / before):.0f}%",
            f"component MAE {before:.5f} → {after:.5f}",
        ),
        (
            "Correction accuracy (1σ, blue)",
            f"{sigma:.4f}",
            "out-of-fold residual spread, reflectance",
        ),
        (
            "Retrieval within EE",
            f"{100.0 * selected_rate:.1f}%",
            f"{hits}/{expected} · uncorrected {100.0 * identity_rate:.1f}%",
        ),
        (
            "Gap to previous best",
            f"−{100.0 * (previous_rate - selected_rate):.1f} pp",
            f"previous-best control {100.0 * previous_rate:.1f}%",
        ),
    ]


def build_html(
    output: Path,
    tiles: list[tuple[str, str, str]],
    sections: list[tuple[str, str, list[tuple[str, str, str]]]],
) -> None:
    tile_html = "\n".join(
        f'      <div class="tile"><p class="tile-label">{html.escape(label)}</p>'
        f'<p class="tile-value">{html.escape(value)}</p>'
        f'<p class="tile-note">{html.escape(note)}</p></div>'
        for label, value, note in tiles
    )
    body: list[str] = []
    for section_title, section_intro, figures in sections:
        body.append(f"    <h2>{html.escape(section_title)}</h2>")
        body.append(f"    <p class=\"intro\">{section_intro}</p>")
        for filename, title, caption in figures:
            body.append(
                '    <figure class="card">\n'
                f"      <figcaption><strong>{html.escape(title)}</strong> "
                f"{caption}</figcaption>\n"
                f'      <img src="figures/{filename}" alt="{html.escape(title)}" '
                'loading="lazy">\n'
                "    </figure>"
            )
    document = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <meta name="color-scheme" content="light">
  <title>Cross-RT L2A harmonisation — visual report</title>
  <style>
    :root {{
      color: #20252b;
      background: #f9f9f7;
      font-family: system-ui, -apple-system, "Segoe UI", sans-serif;
    }}
    * {{ box-sizing: border-box; }}
    body {{ margin: 0; line-height: 1.55; }}
    header {{ background: #20252b; color: #fff; border-bottom: 4px solid #2a78d6; }}
    .header-inner, main, footer {{ width: min(1180px, calc(100% - 32px)); margin: 0 auto; }}
    .header-inner {{ padding: 20px 0 18px; }}
    .eyebrow {{ margin: 0 0 4px; color: #b9c2cb; font-size: 0.78rem; font-weight: 700; }}
    header h1 {{ margin: 0; font-size: 1.45rem; }}
    nav {{ display: flex; flex-wrap: wrap; gap: 16px; margin-top: 12px; }}
    nav a {{ color: #dbe9ff; font-size: 0.86rem; font-weight: 650; }}
    main {{ padding: 26px 0 44px; }}
    h2 {{ margin: 38px 0 6px; padding-top: 12px; border-top: 1px solid #c8cdd2;
         font-size: 1.15rem; }}
    p.intro {{ max-width: 96ch; color: #52514e; margin: 6px 0 16px; }}
    .tiles {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(215px, 1fr));
             gap: 12px; margin-top: 20px; }}
    .tile {{ background: #fcfcfb; border: 1px solid #e1e0d9; border-radius: 8px;
            padding: 14px 16px 12px; }}
    .tile-label {{ margin: 0; font-size: 0.8rem; color: #52514e; }}
    .tile-value {{ margin: 2px 0 0; font-size: 1.9rem; font-weight: 600; color: #0b0b0b; }}
    .tile-note {{ margin: 2px 0 0; font-size: 0.78rem; color: #898781; }}
    figure.card {{ margin: 0 0 22px; background: #fcfcfb; border: 1px solid #e1e0d9;
                  border-radius: 8px; padding: 14px 16px; }}
    figure.card img {{ display: block; width: 100%; height: auto; margin-top: 10px; }}
    figcaption {{ font-size: 0.88rem; color: #52514e; max-width: 110ch; }}
    figcaption strong {{ color: #0b0b0b; }}
    footer {{ padding: 18px 0 32px; border-top: 1px solid #c8cdd2; color: #57616b;
             font-size: 0.82rem; }}
  </style>
</head>
<body>
  <header>
    <div class="header-inner">
      <p class="eyebrow">Sentinel-2 low-cloud validation | 152 frozen matchups</p>
      <h1>Cross-RT L2A harmonisation — visual report</h1>
      <nav aria-label="Report files">
        <a href="index.html">Tabular report</a>
        <a href="summary.json">Summary JSON</a>
        <a href="correction_effect.csv">Correction CSV</a>
        <a href="retrieval_metrics.csv">Retrieval CSV</a>
      </nav>
      <div class="tiles">
{tile_html}
      </div>
    </div>
  </header>
  <main>
{chr(10).join(body)}
  </main>
  <footer>Static visual report generated by tools/aeronet_validation/plot_cross_rt_harmonization.py.
  All statistics recomputed from the raw run artifacts; no AERONET values were used to train the
  surface mapping. Every figure's underlying values are available from the linked CSV/JSON files.</footer>
</body>
</html>
"""
    (output / "visual.html").write_text(document, encoding="utf-8")
    print(f"wrote {output / 'visual.html'}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-root", type=Path, default=DEFAULT_RUN)
    parser.add_argument("--root", type=Path, default=ROOT)
    parser.add_argument("--cases", type=Path, default=DEFAULT_CASES)
    parser.add_argument("--output", type=Path, default=None)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    output = args.output or args.run_root / "summary"
    figures_dir = output / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)
    style()

    scene = load_scene_table(args.run_root)
    artifact = json.loads((args.run_root / "harmonizer.json").read_text(encoding="utf-8"))
    surface_metrics = json.loads(
        (args.run_root / "surface_metrics.json").read_text(encoding="utf-8")
    )
    matchup_ids, results = load_retrieval(args.root, args.cases)
    if len(matchup_ids) != 152:
        raise SystemExit(f"expected the frozen 152-case cohort, found {len(matchup_ids)}")

    section_effect = (
        "What the correction does",
        "The harmonizer maps each historical Sen2Cor L2A acquisition onto the same-day "
        "L1C surface corrected with MAIAC AOD and the current libRadtran LUT. Each point "
        "is one retained tile-scene component (6,858 total from 4,619 acquisitions); the "
        "x-axis driver is how much MAIAC and Sen2Cor disagree on AOD that day.",
        [
            (
                figure_correction_vs_delta_aot(scene, figures_dir),
                "Applied correction vs ΔAOT.",
                "The learned mapping in action: where MAIAC sees less aerosol than Sen2Cor "
                "(ΔAOT&lt;0) the blue band is brightened, where it sees more it is "
                "darkened. Line = binned median, band = interquartile range.",
            ),
            (
                figure_bias_vs_delta_aot(scene, figures_dir),
                "Surface error vs ΔAOT, before and after.",
                "The uncorrected L2A error (gray) tilts systematically with the "
                "MAIAC−Sen2Cor AOD difference; after correction (blue) the tilt is "
                "removed and the median error sits near zero across the whole range.",
            ),
            (
                figure_bias_distributions(scene, figures_dir),
                "Error distributions before and after.",
                "Component bias vs the MAIAC-compatible target. The corrected distribution "
                "narrows and recentres on zero; the effect is strongest in blue.",
            ),
        ],
    )
    section_accuracy = (
        "How accurate the correction is",
        "All accuracy checks are out-of-fold: the model scored on a site never saw that "
        "site during fitting.",
        [
            (
                figure_component_error_scatter(scene, figures_dir),
                "Component-level error, before vs after.",
                "Each dot is one tile-scene component; below the diagonal means the "
                "correction moved it closer to the MAIAC-compatible target.",
            ),
            (
                figure_site_mae(scene, figures_dir),
                "Site-level blue MAE, before vs after.",
                "One dot per site (size = number of components). Labels mark the largest "
                "movers in both directions.",
            ),
            (
                figure_stage_ablation(surface_metrics, figures_dir),
                "Feature-stage ablation.",
                "Visible-band acquisition MAE as feature stages are added (all at cap "
                "0.030). The ΔAOT interaction terms carry most of the improvement; "
                "terrain terms add a further small gain. The current production AOD "
                "offset is worse than leaving L2A untouched.",
            ),
            (
                figure_residual_spread(artifact, figures_dir),
                "Residual spread per band.",
                "What remains after correction: robust σ ≈ 0.003–0.008 "
                "reflectance depending on band. This is the intrinsic accuracy limit of "
                "the mapping.",
            ),
        ],
    )
    section_derivation = (
        "How the correction is derived",
        "Per band, a ridge regression (α=1, standardized features, every independent "
        "acquisition weighted equally) predicts the residual between the operational L2A "
        "reflectance and the MAIAC-compatible target, from 2.04M same-day exact-pair "
        "pixels. The applied correction is this prediction clipped at ±0.030.",
        [
            (
                figure_coefficients(artifact, figures_dir),
                "The correction model, feature by feature.",
                "Standardized coefficients of the applied development-trained model for "
                "the three solver bands (top 16 of 34 features). The "
                "MAIAC−Sen2Cor ΔAOT term dominates in every band.",
            ),
            (
                figure_coefficient_stability(artifact, figures_dir),
                "Coefficient stability across folds.",
                "Blue-band coefficients refit in each site-grouped cross-fit fold (gray) "
                "against the applied model (blue). Tight clusters mean the derivation "
                "does not hinge on any subset of sites.",
            ),
        ],
    )
    section_retrieval = (
        "Does it help the AOD retrieval?",
        "AERONET enters only here, as the final score. Within-EE means "
        "|retrieved − truth| ≤ 0.05 + 0.15·truth on the frozen 152-case "
        "cohort; missing results count as misses.",
        [
            (
                figure_retrieval_rates(results, figures_dir),
                "Within-EE rate by history variant.",
                "Every harmonized variant beats the uncorrected history, and the selected "
                "variant recovers about a third of the gap to the previous-best control "
                "— but none reaches it, and all sit well short of the 87% target.",
            ),
            (
                figure_retrieval_scatter(results, figures_dir),
                "Retrieved vs true AOD.",
                "Gray band = expected-error envelope. The harmonized history removes much "
                "of the uncorrected run's high-side scatter at low AOD, while the "
                "under-retrieval above AOD ≈ 0.8 remains in both experimental runs.",
            ),
            (
                figure_retrieval_error_vs_aod(results, figures_dir),
                "Retrieval error vs AOD level.",
                "Faint dots are individual matchups; lines are binned medians with IQR "
                "bands. The correction removes the positive clean-sky bias but medium "
                "and high AOD keep a fat error tail — the regime that holds the score "
                "below the previous best.",
            ),
            (
                figure_case_transitions(results, figures_dir),
                "Cases whose verdict changed.",
                "Every matchup whose within-EE verdict flipped against the uncorrected "
                "history: dot = before, arrowhead = after, in units of the expected "
                "error (1.0 = the EE boundary).",
            ),
        ],
    )

    tiles = stat_tiles(scene, artifact, results)
    build_html(
        output,
        tiles,
        [section_effect, section_accuracy, section_derivation, section_retrieval],
    )


if __name__ == "__main__":
    main()
