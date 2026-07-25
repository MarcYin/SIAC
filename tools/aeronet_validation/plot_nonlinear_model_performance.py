"""Extensive performance visuals for the nonlinear L2A harmonizer.

Recomputes out-of-fold per-pixel predictions from the saved fold models for
both the gradient-boosted (target-band) and ridge harmonizers under the same
locked site-grouped folds, joins both runs' per-component tables, and renders:
pixel prediction-vs-truth scatters, error distributions, component-level model
comparison, conditional performance across the sample dimensions (AOD, ΔAOT,
water vapour, geometry, terrain, season), site-level maps, and spatial images
of the applied correction from the harmonized history rasters.
"""

from __future__ import annotations

import argparse
import html
import json
from pathlib import Path
from typing import Any

import joblib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.ticker import NullFormatter
from tools.aeronet_validation.plot_cross_rt_harmonization import (
    AXIS,
    BLUE,
    GRID,
    INK,
    LIGHT_BLUE,
    MUTED,
    RED,
    SECONDARY,
    SURFACE,
    binned_quantiles,
    save,
    style,
)
from tools.aeronet_validation.train_l2a_l1c_harmonizer import (
    BAND_NAMES,
    _feature_kwargs,
    _predict_model,
    build_features,
    load_pairs,
)
from tools.aeronet_validation.train_l2a_l1c_nonlinear_harmonizer import (
    _case_ids,
    build_columns,
)

ROOT = Path("/gws/ssde/j25a/nceo_isp/public/siac_refactor")
RIDGE_RUN = ROOT / "analysis/cross_rt_harmonizer_lowcloud152_20260716"
HGB_RUN = ROOT / "analysis/cross_rt_nonlinear_targetband_lowcloud152_20260716"
DEFAULT_PAIRS = ROOT / "analysis/l2a_l1c_physical_pairs_lowcloud152_20260716"
DEFAULT_CASES = ROOT / "reports/aod-final-performance-dashboard-20260713/data/all-cases.csv"
MATCHUPS = ROOT / "matchups/matchups.csv"

RIDGE_MODEL = "cross_rt_terrain_a1"
RIDGE_PREFIX = "cross_rt_terrain_a1_cap_0p030"
HGB_MODEL = "hgb_target_band"
HGB_PREFIX = "hgb_target_band_cap_0p030"
SOLVER_BANDS = (("blue", "B02 blue"), ("green", "B03 green"), ("red", "B04 red"))

GREEN = "#008300"  # categorical slot 2 — validated adjacent to slot-1 blue
DENSITY = LinearSegmentedColormap.from_list("blues", ["#cde2fb", "#3987e5", "#0d366b"])
DIVERGING = LinearSegmentedColormap.from_list("blue_gray_red", [BLUE, "#f0efec", RED])
MODEL_SERIES = (
    ("identity", "uncorrected L2A", MUTED),
    ("ridge", "ridge (linear)", GREEN),
    ("hgb", "gradient boosting (nonlinear)", BLUE),
)


# ---------------------------------------------------------------------------
# Out-of-fold pixel predictions
# ---------------------------------------------------------------------------


def applied_fold_per_sample(dataset: Any, artifact: dict[str, Any]) -> np.ndarray:
    folds = np.full(dataset.site_code.shape, -1, dtype=np.int16)
    fold_by_matchup = artifact["fold_by_matchup_id"]
    for site_code, matchup_ids in dataset.matchup_by_site_code.items():
        folds[dataset.site_code == site_code] = int(fold_by_matchup[matchup_ids[0]])
    if np.any(folds < 0):
        raise RuntimeError("some samples received no applied fold")
    return folds


def oof_predictions(
    dataset: Any,
    hgb_artifact: dict[str, Any],
    ridge_artifact: dict[str, Any],
) -> dict[str, dict[str, np.ndarray]]:
    """Per-band OOF predicted residuals for both models, plus the true residual."""
    columns = build_columns(dataset)
    kwargs = _feature_kwargs(dataset)
    folds = applied_fold_per_sample(dataset, hgb_artifact)
    ridge_folds = applied_fold_per_sample(dataset, ridge_artifact)
    if not np.array_equal(folds, ridge_folds):
        raise RuntimeError("the two artifacts disagree on the locked fold assignment")
    result: dict[str, dict[str, np.ndarray]] = {}
    for band_index, band in enumerate(BAND_NAMES):
        truth = (dataset.siac[:, band_index] - dataset.l2a[:, band_index]).astype(np.float64)
        hgb_pred = np.full(truth.shape, np.nan)
        ridge_pred = np.full(truth.shape, np.nan)
        names = hgb_artifact["feature_names_by_band"][band]
        features_hgb = np.column_stack([columns[name] for name in names])
        for fold in np.unique(folds):
            mask = folds == fold
            estimator = joblib.load(
                HGB_RUN / "models" / hgb_artifact["model_files"][str(fold)][band]
            )
            hgb_pred[mask] = estimator.predict(features_hgb[mask])
            ridge_band = ridge_artifact["models"][RIDGE_MODEL]["folds"][str(fold)][band]
            ridge_features = build_features(
                dataset.l2a[mask],
                band_index=band_index,
                feature_set="cross_rt_terrain",
                **{key: value[mask] for key, value in kwargs.items()},
            )
            ridge_pred[mask] = _predict_model(ridge_features, ridge_band)
        result[band] = {"truth": truth, "hgb": hgb_pred, "ridge": ridge_pred}
    return result


def figure_pred_vs_truth(predictions: dict[str, dict[str, np.ndarray]], output_dir: Path) -> str:
    rng = np.random.default_rng(0)
    limits = (-0.045, 0.045)
    fig, axes = plt.subplots(2, 3, figsize=(11.8, 7.6), sharex=True, sharey=True)
    for column, (band, label) in enumerate(SOLVER_BANDS):
        truth = predictions[band]["truth"]
        sample = rng.choice(truth.size, size=min(250_000, truth.size), replace=False)
        for row, model_key, model_label in ((0, "hgb", "gradient boosting"), (1, "ridge", "ridge")):
            ax = axes[row, column]
            predicted = predictions[band][model_key]
            residual = predicted - truth
            sigma = float(np.median(np.abs(residual - np.median(residual))) * 1.4826)
            variance = float(np.var(truth))
            r2 = 1.0 - float(np.mean(np.square(residual))) / variance if variance else np.nan
            ax.hexbin(
                truth[sample],
                predicted[sample],
                gridsize=70,
                cmap=DENSITY,
                mincnt=1,
                extent=(*limits, *limits),
                linewidths=0,
            )
            ax.plot(limits, limits, color=AXIS, linewidth=0.8)
            ax.axhline(0.0, color=GRID, linewidth=0.8)
            ax.axvline(0.0, color=GRID, linewidth=0.8)
            ax.set_xlim(*limits)
            ax.set_ylim(*limits)
            ax.set_title(f"{label} — {model_label}", fontsize=10)
            ax.text(
                0.03,
                0.95,
                f"R² {r2:.2f} · robust σ {sigma:.4f}",
                transform=ax.transAxes,
                ha="left",
                va="top",
                fontsize=8.5,
                color=SECONDARY,
            )
            if column == 0:
                ax.set_ylabel("predicted correction")
            if row == 1:
                ax.set_xlabel("true residual (target − L2A)")
    fig.set_constrained_layout(True)
    return save(fig, output_dir, "n01_pred_vs_truth.png")


def figure_pixel_error_distributions(
    predictions: dict[str, dict[str, np.ndarray]], output_dir: Path
) -> str:
    bins = np.arange(-0.03, 0.0305, 0.001)
    fig, axes = plt.subplots(1, 3, figsize=(11.8, 3.5), sharey=True, constrained_layout=True)
    for ax, (band, label) in zip(axes, SOLVER_BANDS, strict=True):
        truth = predictions[band]["truth"]
        series = {
            "identity": -truth,
            "ridge": predictions[band]["ridge"] - truth,
            "hgb": predictions[band]["hgb"] - truth,
        }
        for key, name, color in MODEL_SERIES:
            values = series[key]
            if key == "identity":
                ax.hist(
                    values,
                    bins=bins,
                    density=True,
                    histtype="stepfilled",
                    color=color,
                    alpha=0.25,
                    linewidth=0,
                    label=name,
                )
            else:
                ax.hist(
                    values,
                    bins=bins,
                    density=True,
                    histtype="step",
                    color=color,
                    linewidth=2,
                    label=name,
                )
        ax.axvline(0.0, color=AXIS, linewidth=0.8)
        ax.set_title(label)
        ax.set_xlim(-0.025, 0.025)
        ax.set_yticks([])
        ax.grid(False, axis="y")
        ax.set_xlabel("pixel error after correction")
    axes[0].set_ylabel("density")
    axes[0].legend(loc="upper left", fontsize=8)
    return save(fig, output_dir, "n02_pixel_error_distributions.png")


def _corrected(dataset: Any, predictions: dict[str, Any], band: str, cap: float = 0.03):
    band_index = BAND_NAMES.index(band)
    return dataset.l2a[:, band_index] + np.clip(predictions[band]["hgb"], -cap, cap)


def figure_error_decomposition(
    dataset: Any, predictions: dict[str, dict[str, np.ndarray]], output_dir: Path
) -> str:
    """Split the residual variance into scene-common and within-scene parts.

    The scene-common (between-acquisition) part is what a per-scene error in the
    TARGET itself — one MAIAC AOD value and scene-mean geometry per acquisition —
    would produce; no per-pixel feature can remove it.
    """
    acquisition = np.asarray(dataset.acquisition_code, dtype=np.int64)
    counts = np.bincount(acquisition)
    panels = (
        ("uncorrected L2A", None),
        ("ridge (linear)", "ridge"),
        ("gradient boosting", "hgb"),
    )
    fig, axes = plt.subplots(1, 3, figsize=(11.8, 3.9), sharex=True, constrained_layout=True)
    positions = np.arange(len(BAND_NAMES))[::-1]
    decomposed: list[tuple[Any, list[float], list[float]]] = []
    for ax, (title, model_key) in zip(axes, panels, strict=True):
        between_list, within_list = [], []
        for band in BAND_NAMES:
            truth = predictions[band]["truth"]
            if model_key is None:
                residual = -truth
            else:
                residual = np.clip(predictions[band][model_key], -0.03, 0.03) - truth
            sums = np.bincount(acquisition, weights=residual)
            squares = np.bincount(acquisition, weights=residual**2)
            means = sums / counts
            within = float(np.mean(squares / counts - means**2))
            between = float(np.var(means))
            between_list.append(between * 1.0e6)
            within_list.append(within * 1.0e6)
        decomposed.append((ax, between_list, within_list))
        ax.set_title(title, fontsize=10)
    xmax = 1.3 * max(
        between + within
        for _ax, between_list, within_list in decomposed
        for between, within in zip(between_list, within_list, strict=True)
    )
    for ax, between_list, within_list in decomposed:
        ax.barh(positions, between_list, height=0.62, color=BLUE, zorder=3)
        ax.barh(
            positions,
            within_list,
            left=between_list,
            height=0.62,
            color=LIGHT_BLUE,
            zorder=3,
        )
        for position, between, within in zip(positions, between_list, within_list, strict=True):
            sigma = float(np.sqrt((between + within) * 1.0e-6))
            ax.text(
                between + within + 0.012 * xmax,
                position,
                f"σ {sigma:.4f}",
                va="center",
                ha="left",
                fontsize=8,
                color=SECONDARY,
            )
        ax.set_yticks(positions)
        ax.set_yticklabels(BAND_NAMES, fontsize=8.5, color=INK)
        ax.grid(False, axis="y")
        ax.set_xlabel("residual variance (×10⁻⁶)")
        ax.set_xlim(0, xmax)
    handles = [
        plt.Rectangle((0, 0), 1, 1, color=BLUE, label="scene-common (between acquisitions)"),
        plt.Rectangle((0, 0), 1, 1, color=LIGHT_BLUE, label="within acquisition"),
    ]
    fig.legend(handles=handles, loc="outside upper right", ncols=2, fontsize=8.5)
    return save(fig, output_dir, "n08_error_decomposition.png")


def figure_band_scatter_global(
    dataset: Any, predictions: dict[str, dict[str, np.ndarray]], output_dir: Path
) -> str:
    rng = np.random.default_rng(0)
    total = dataset.l2a.shape[0]
    sample = rng.choice(total, size=min(150_000, total), replace=False)
    limits = (0.0, 0.55)
    fig, axes = plt.subplots(2, len(BAND_NAMES), figsize=(13.4, 4.5), sharex=True, sharey=True)
    for column, band in enumerate(BAND_NAMES):
        band_index = column
        target = dataset.siac[:, band_index].astype(np.float64)
        rows = (
            ("L2A before", dataset.l2a[:, band_index].astype(np.float64)),
            ("corrected after", _corrected(dataset, predictions, band)),
        )
        for row, (name, values) in enumerate(rows):
            ax = axes[row, column]
            ax.hexbin(
                target[sample],
                values[sample],
                gridsize=55,
                cmap=DENSITY,
                mincnt=1,
                extent=(*limits, *limits),
                linewidths=0,
            )
            ax.plot(limits, limits, color=AXIS, linewidth=0.7)
            mae = float(np.nanmean(np.abs(values - target)))
            ax.text(
                0.05,
                0.94,
                f"MAE {mae:.4f}",
                transform=ax.transAxes,
                ha="left",
                va="top",
                fontsize=7.5,
                color=SECONDARY,
            )
            ax.set_xlim(*limits)
            ax.set_ylim(*limits)
            ax.set_xticks([0.0, 0.25, 0.5])
            ax.set_yticks([0.0, 0.25, 0.5])
            ax.tick_params(labelsize=7)
            if row == 0:
                ax.set_title(band, fontsize=9.5)
            if column == 0:
                ax.set_ylabel(name, fontsize=8.5)
            if row == 1:
                ax.set_xlabel("target reflectance", fontsize=8)
    fig.set_constrained_layout(True)
    return save(fig, output_dir, "n09_band_scatter_global.png")


def matchup_of_sample(dataset: Any) -> np.ndarray:
    matchups = np.asarray(
        [str(scene.get("matchup_id")) for scene in dataset.scene_metadata], dtype=object
    )
    return matchups[np.asarray(dataset.scene_code, dtype=np.int64)]


def figure_band_scatter_example(
    dataset: Any,
    predictions: dict[str, dict[str, np.ndarray]],
    sample_matchups: np.ndarray,
    matchup_id: str,
    output_dir: Path,
    index: int,
) -> tuple[str, str] | None:
    mask = sample_matchups == matchup_id
    if int(mask.sum()) < 200:
        return None
    site = matchup_id.split("__", 1)[0]
    fig, axes = plt.subplots(2, 4, figsize=(11.8, 6.0), sharex=True, sharey=True)
    upper = min(0.6, float(np.nanpercentile(dataset.siac[mask], 99.8)) * 1.25 + 0.02)
    for panel, band in enumerate(BAND_NAMES):
        ax = axes.ravel()[panel]
        band_index = panel
        target = dataset.siac[mask, band_index].astype(np.float64)
        before = dataset.l2a[mask, band_index].astype(np.float64)
        after = _corrected(dataset, predictions, band)[mask]
        ax.scatter(target, before, s=4, color=MUTED, alpha=0.3, linewidths=0, rasterized=True)
        ax.scatter(target, after, s=4, color=BLUE, alpha=0.3, linewidths=0, rasterized=True)
        ax.plot([0, upper], [0, upper], color=AXIS, linewidth=0.8)
        ax.set_xlim(0, upper)
        ax.set_ylim(0, upper)
        ax.set_title(band, fontsize=9.5)
        ax.text(
            0.04,
            0.95,
            f"MAE {np.nanmean(np.abs(before - target)):.4f} → "
            f"{np.nanmean(np.abs(after - target)):.4f}",
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=7.5,
            color=SECONDARY,
        )
    legend_ax = axes.ravel()[-1]
    legend_ax.axis("off")
    handles = [
        plt.Line2D([], [], marker="o", linestyle="", color=MUTED, label="L2A before"),
        plt.Line2D([], [], marker="o", linestyle="", color=BLUE, label="corrected (nonlinear)"),
    ]
    legend_ax.legend(handles=handles, loc="center", fontsize=9)
    for ax in axes[1, :3]:
        ax.set_xlabel("target reflectance (L1C + MAIAC AOD)", fontsize=8)
    for ax in axes[:, 0]:
        ax.set_ylabel("L2A reflectance", fontsize=8)
    fig.suptitle(
        f"{site} — {int(mask.sum())} paired pixels, all bands",
        fontsize=11,
        fontweight="bold",
        color=INK,
    )
    fig.set_constrained_layout(True)
    name = save(fig, output_dir, f"n2{index}_bands_{site.lower()[:24]}.png")
    caption = (
        f"{site}: every paired pixel of this matchup, original L2A (gray) and corrected "
        "(blue) against the same-day L1C surface corrected with MAIAC AOD. The "
        "correction tightens the visible bands onto the diagonal; NIR/SWIR sit close "
        "to it already."
    )
    return name, caption


# ---------------------------------------------------------------------------
# Component-level comparison (both runs' scene tables joined)
# ---------------------------------------------------------------------------


def load_joined_components() -> pd.DataFrame:
    keep = [
        "scene_code",
        "matchup_id",
        "site",
        "scene_id",
        "window",
        "sample_count",
        "maiac_aot",
        "l2a_aot",
        "l2a_tcwv_cm",
        "sza_deg",
        "elevation_km",
        "terrain_slope_deg_mean",
        "delta_aot_maiac_minus_sen2cor",
    ]
    hgb_columns = keep + [
        column
        for band, _label in SOLVER_BANDS
        for column in (
            f"identity_{band}_bias",
            f"identity_{band}_mae",
            f"{HGB_PREFIX}_{band}_bias",
            f"{HGB_PREFIX}_{band}_mae",
        )
    ]
    hgb = pd.read_csv(HGB_RUN / "surface_scene_metrics.csv", usecols=hgb_columns)
    ridge_columns = ["scene_code", "scene_id"] + [
        column
        for band, _label in SOLVER_BANDS
        for column in (f"{RIDGE_PREFIX}_{band}_bias", f"{RIDGE_PREFIX}_{band}_mae")
    ]
    ridge = pd.read_csv(RIDGE_RUN / "surface_scene_metrics.csv", usecols=ridge_columns)
    joined = hgb.merge(ridge, on="scene_code", suffixes=("", "_ridge"), validate="1:1")
    if not (joined["scene_id"] == joined["scene_id_ridge"]).all():
        raise RuntimeError("the two runs' component tables are not aligned")
    joined["month"] = pd.to_datetime(
        joined["window"].astype(str), format="%Y-%m", errors="coerce"
    ).dt.month
    return joined


def figure_component_model_scatter(components: pd.DataFrame, output_dir: Path) -> str:
    limits = (2.0e-4, 0.06)
    fig, axes = plt.subplots(1, 3, figsize=(11.8, 4.3), constrained_layout=True)
    norm = Normalize(vmin=-0.2, vmax=0.2)
    for ax, (band, label) in zip(axes, SOLVER_BANDS, strict=True):
        ridge_abs = components[f"{RIDGE_PREFIX}_{band}_bias"].abs().clip(lower=limits[0])
        hgb_abs = components[f"{HGB_PREFIX}_{band}_bias"].abs().clip(lower=limits[0])
        finite = np.isfinite(ridge_abs) & np.isfinite(hgb_abs)
        better = 100.0 * float((hgb_abs[finite] < ridge_abs[finite]).mean())
        scatter = ax.scatter(
            ridge_abs[finite],
            hgb_abs[finite],
            s=6,
            c=components.loc[finite, "delta_aot_maiac_minus_sen2cor"],
            cmap=DIVERGING,
            norm=norm,
            alpha=0.5,
            linewidths=0,
            rasterized=True,
        )
        ax.plot(limits, limits, color=AXIS, linewidth=0.8)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(*limits)
        ax.set_ylim(*limits)
        ax.set_title(label)
        ax.set_xlabel("ridge component |bias|")
        ax.text(
            0.03,
            0.95,
            f"nonlinear better on {better:.0f}%",
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=9,
            color=SECONDARY,
        )
    axes[0].set_ylabel("gradient-boosting component |bias|")
    colorbar = fig.colorbar(scatter, ax=axes, shrink=0.85, pad=0.01)
    colorbar.set_label("ΔAOT = MAIAC − Sen2Cor", fontsize=9, color=SECONDARY)
    colorbar.outline.set_visible(False)
    return save(fig, output_dir, "n03_component_model_scatter.png")


def figure_conditional_performance(components: pd.DataFrame, output_dir: Path) -> str:
    drivers = (
        ("maiac_aot", "MAIAC AOD", np.geomspace(0.02, 1.2, 9), "log"),
        (
            "delta_aot_maiac_minus_sen2cor",
            "ΔAOT (MAIAC − Sen2Cor)",
            np.arange(-0.275, 0.425, 0.05),
            "linear",
        ),
        ("l2a_tcwv_cm", "L2A water vapour (cm)", np.linspace(0.0, 5.5, 10), "linear"),
        ("sza_deg", "solar zenith (deg)", np.linspace(15.0, 70.0, 10), "linear"),
        ("terrain_slope_deg_mean", "terrain slope (deg)", np.linspace(0.0, 18.0, 9), "linear"),
        ("month", "month", np.arange(0.5, 13.0, 1.0), "linear"),
    )
    series = {
        "identity": components["identity_blue_bias"].abs().to_numpy(),
        "ridge": components[f"{RIDGE_PREFIX}_blue_bias"].abs().to_numpy(),
        "hgb": components[f"{HGB_PREFIX}_blue_bias"].abs().to_numpy(),
    }
    fig, axes = plt.subplots(2, 3, figsize=(11.8, 6.6), sharey=True, constrained_layout=True)
    for ax, (column, label, edges, scale) in zip(axes.ravel(), drivers, strict=True):
        x = components[column].to_numpy(dtype=float)
        for key, name, color in MODEL_SERIES:
            centers, median, lower, upper = binned_quantiles(
                x, series[key], np.asarray(edges), min_count=40
            )
            ax.fill_between(centers, lower, upper, color=color, alpha=0.12, linewidth=0)
            ax.plot(centers, median, color=color, linewidth=2, solid_capstyle="round", label=name)
        if scale == "log":
            ax.set_xscale("log")
            ax.xaxis.set_minor_formatter(NullFormatter())
        ax.set_ylim(0, 0.0125)
        ax.set_xlabel(label)
    for ax in axes[:, 0]:
        ax.set_ylabel("component |bias|, B02 blue")
    axes[0, 0].legend(loc="upper left", fontsize=8)
    return save(fig, output_dir, "n04_conditional_performance.png")


def figure_component_ecdf(components: pd.DataFrame, output_dir: Path) -> str:
    fig, axes = plt.subplots(1, 3, figsize=(11.8, 3.6), sharey=True, constrained_layout=True)
    for ax, (band, label) in zip(axes, SOLVER_BANDS, strict=True):
        columns = {
            "identity": f"identity_{band}_bias",
            "ridge": f"{RIDGE_PREFIX}_{band}_bias",
            "hgb": f"{HGB_PREFIX}_{band}_bias",
        }
        for key, name, color in MODEL_SERIES:
            values = np.sort(components[columns[key]].abs().dropna().to_numpy())
            fraction = np.arange(1, values.size + 1) / values.size
            ax.plot(values, 100.0 * fraction, color=color, linewidth=2, label=name)
        ax.set_xlim(0, 0.02)
        ax.set_ylim(0, 100)
        ax.set_title(label)
        ax.set_xlabel("component |bias|")
    axes[0].set_ylabel("% of components below")
    axes[0].legend(loc="lower right", fontsize=8)
    return save(fig, output_dir, "n05_component_ecdf.png")


# ---------------------------------------------------------------------------
# Site-level views
# ---------------------------------------------------------------------------


def site_table(components: pd.DataFrame) -> pd.DataFrame:
    sites = components.groupby("site").agg(
        hgb=(f"{HGB_PREFIX}_blue_mae", "mean"),
        ridge=(f"{RIDGE_PREFIX}_blue_mae", "mean"),
        identity=("identity_blue_mae", "mean"),
        components=("site", "size"),
    )
    coordinates = pd.read_csv(MATCHUPS, usecols=["matchup_id", "latitude", "longitude"])
    coordinates["site"] = coordinates["matchup_id"].str.split("__").str[0]
    coordinates = coordinates.groupby("site")[["latitude", "longitude"]].first()
    return sites.join(coordinates).dropna(subset=["hgb", "ridge"])


def figure_site_map(sites: pd.DataFrame, output_dir: Path) -> str:
    change = 100.0 * (sites["hgb"] / sites["ridge"] - 1.0)
    fig, ax = plt.subplots(figsize=(11.8, 5.6), constrained_layout=True)
    norm = Normalize(vmin=-30.0, vmax=30.0)
    scatter = ax.scatter(
        sites["longitude"],
        sites["latitude"],
        c=np.clip(change, -30.0, 30.0),
        s=np.clip(sites["components"] * 0.8, 18.0, 130.0),
        cmap=DIVERGING,
        norm=norm,
        linewidths=1.2,
        edgecolors=SURFACE,
    )
    for name, row in sites.loc[change.abs().nlargest(8).index].iterrows():
        ax.annotate(
            str(name),
            (row["longitude"], row["latitude"]),
            textcoords="offset points",
            xytext=(6, 4),
            fontsize=7.5,
            color=SECONDARY,
        )
    ax.set_xlim(-180, 180)
    ax.set_ylim(-65, 80)
    ax.set_xlabel("longitude")
    ax.set_ylabel("latitude")
    ax.set_title(
        "Site map — blue-band MAE, gradient boosting vs ridge "
        "(blue = nonlinear better, red = ridge better)"
    )
    colorbar = fig.colorbar(scatter, ax=ax, shrink=0.8, pad=0.01)
    colorbar.set_label("MAE change vs ridge (%)", fontsize=9, color=SECONDARY)
    colorbar.outline.set_visible(False)
    return save(fig, output_dir, "n06_site_map.png")


def figure_site_scatter(sites: pd.DataFrame, output_dir: Path) -> str:
    fig, ax = plt.subplots(figsize=(6.4, 5.8), constrained_layout=True)
    limit = float(max(sites["hgb"].max(), sites["ridge"].max())) * 1.08
    ax.scatter(
        sites["ridge"],
        sites["hgb"],
        s=np.clip(sites["components"] * 0.65, 12.0, 90.0),
        color=BLUE,
        alpha=0.55,
        linewidths=1.2,
        edgecolors=SURFACE,
    )
    ax.plot([0, limit], [0, limit], color=AXIS, linewidth=0.8)
    ax.set_xlim(0, limit)
    ax.set_ylim(0, limit)
    change = sites["hgb"] - sites["ridge"]
    movers = pd.concat([change.nlargest(3), change.nsmallest(3)])
    offsets = ((7, -3), (7, 9), (7, -15))
    for index, name in enumerate(movers.index):
        row = sites.loc[name]
        ax.annotate(
            str(name),
            (row["ridge"], row["hgb"]),
            textcoords="offset points",
            xytext=offsets[index % len(offsets)],
            fontsize=8,
            color=SECONDARY,
        )
    better = 100.0 * float((sites["hgb"] < sites["ridge"]).mean())
    ax.set_title(f"Per-site blue MAE — nonlinear better at {better:.0f}% of {len(sites)} sites")
    ax.set_xlabel("ridge site MAE")
    ax.set_ylabel("gradient-boosting site MAE")
    return save(fig, output_dir, "n07_site_scatter.png")


# ---------------------------------------------------------------------------
# Spatial imagery from the harmonized history rasters
# ---------------------------------------------------------------------------


def available_spatial_cases(components: pd.DataFrame, limit: int = 4) -> list[str]:
    spread = (
        components.groupby("matchup_id")["delta_aot_maiac_minus_sen2cor"]
        .agg(lambda values: float(np.nanpercentile(np.abs(values), 95)))
        .sort_values(ascending=False)
    )
    selected: list[str] = []
    for matchup_id in spread.index:
        paths = (
            HGB_RUN / "daily_histories/identity_daily" / f"{matchup_id}.npz",
            HGB_RUN / "daily_histories/hgb_target_band_blue_cap0p030" / f"{matchup_id}.npz",
            RIDGE_RUN / "daily_histories/cross_rt_terrain_a1_solver_cap0p030" / f"{matchup_id}.npz",
        )
        if all(path.exists() for path in paths):
            selected.append(matchup_id)
        if len(selected) >= limit:
            break
    return selected


def _history(path: Path) -> tuple[np.ndarray, list[str]]:
    with np.load(path, allow_pickle=False) as payload:
        return (
            np.asarray(payload["comp"], dtype=np.float64),
            [str(value) for value in payload["realizations"]],
        )


def figure_spatial_case(
    matchup_id: str, components: pd.DataFrame, output_dir: Path, index: int
) -> tuple[str, str] | None:
    identity, identity_windows = _history(
        HGB_RUN / "daily_histories/identity_daily" / f"{matchup_id}.npz"
    )
    hgb, hgb_windows = _history(
        HGB_RUN / "daily_histories/hgb_target_band_blue_cap0p030" / f"{matchup_id}.npz"
    )
    ridge, ridge_windows = _history(
        RIDGE_RUN / "daily_histories/cross_rt_terrain_a1_solver_cap0p030" / f"{matchup_id}.npz"
    )
    shared = [
        window
        for window in identity_windows
        if window in set(hgb_windows) and window in set(ridge_windows)
    ]
    if not shared:
        return None
    case = components[components["matchup_id"] == matchup_id]
    by_window = case.groupby("window")["delta_aot_maiac_minus_sen2cor"].median()
    window = max(shared, key=lambda name: abs(float(by_window.get(name, 0.0))))
    delta_aot = float(by_window.get(window, float("nan")))
    slices = {
        "identity": identity[identity_windows.index(window)],
        "hgb": hgb[hgb_windows.index(window)],
        "ridge": ridge[ridge_windows.index(window)],
    }
    band_index = {band: position for position, band in enumerate(BAND_NAMES)}
    rgb = np.stack(
        [slices["identity"][band_index[band]] for band in ("red", "green", "blue")], axis=-1
    )
    rgb = np.clip(rgb / 0.25, 0.0, 1.0) ** (1.0 / 1.6)
    rgb = np.where(np.isfinite(rgb), rgb, 1.0)
    fig, axes = plt.subplots(1, 4, figsize=(12.6, 3.6), constrained_layout=True)
    axes[0].imshow(rgb, interpolation="nearest")
    axes[0].set_title("L2A composite (true colour)", fontsize=9.5)
    blue = band_index["blue"]
    axes[1].imshow(
        slices["identity"][blue], cmap="gray", vmin=0.0, vmax=0.18, interpolation="nearest"
    )
    axes[1].set_title("B02 blue reflectance", fontsize=9.5)
    correction_kwargs = {
        "cmap": DIVERGING,
        "vmin": -0.02,
        "vmax": 0.02,
        "interpolation": "nearest",
    }
    axes[2].imshow(slices["hgb"][blue] - slices["identity"][blue], **correction_kwargs)
    axes[2].set_title("nonlinear correction (blue)", fontsize=9.5)
    image = axes[3].imshow(slices["ridge"][blue] - slices["identity"][blue], **correction_kwargs)
    axes[3].set_title("ridge correction (blue)", fontsize=9.5)
    for ax in axes:
        ax.set_xticks([])
        ax.set_yticks([])
        ax.grid(False)
        for spine in ax.spines.values():
            spine.set_visible(False)
    colorbar = fig.colorbar(image, ax=axes[2:], shrink=0.85, pad=0.02)
    colorbar.set_label("applied correction (reflectance)", fontsize=8.5, color=SECONDARY)
    colorbar.outline.set_visible(False)
    site = matchup_id.split("__", 1)[0]
    fig.suptitle(
        f"{site} — window {window}, median ΔAOT {delta_aot:+.3f}",
        fontsize=11,
        fontweight="bold",
        color=INK,
    )
    name = save(fig, output_dir, f"n1{index}_spatial_{site.lower()[:24]}.png")
    caption = (
        f"{site}: monthly composite for {window} (ΔAOT median {delta_aot:+.3f}). The two "
        "right panels show each model's applied blue-band correction on the same "
        "reflectance scale — spatial texture in the correction comes from the per-pixel "
        "reflectance and terrain inputs."
    )
    return name, caption


# ---------------------------------------------------------------------------
# Page assembly
# ---------------------------------------------------------------------------


def build_html(output: Path, sections: list[tuple[str, str, list[tuple[str, str, str]]]]) -> None:
    body: list[str] = []
    for section_title, section_intro, figures in sections:
        body.append(f"    <h2>{html.escape(section_title)}</h2>")
        body.append(f'    <p class="intro">{section_intro}</p>')
        for filename, title, caption in figures:
            body.append(
                '    <figure class="card">\n'
                f"      <figcaption><strong>{html.escape(title)}</strong> {caption}"
                "</figcaption>\n"
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
  <title>Nonlinear L2A harmonizer — model performance</title>
  <style>
    :root {{ color: #20252b; background: #f9f9f7;
            font-family: system-ui, -apple-system, "Segoe UI", sans-serif; }}
    * {{ box-sizing: border-box; }}
    body {{ margin: 0; line-height: 1.55; }}
    header {{ background: #20252b; color: #fff; border-bottom: 4px solid #2a78d6; }}
    .header-inner, main {{ width: min(1180px, calc(100% - 32px)); margin: 0 auto; }}
    .header-inner {{ padding: 20px 0 18px; }}
    .eyebrow {{ margin: 0 0 4px; color: #b9c2cb; font-size: 0.78rem; font-weight: 700; }}
    header h1 {{ margin: 0; font-size: 1.45rem; }}
    nav {{ display: flex; flex-wrap: wrap; gap: 16px; margin-top: 12px; }}
    nav a {{ color: #dbe9ff; font-size: 0.86rem; font-weight: 650; }}
    main {{ padding: 26px 0 44px; }}
    h2 {{ margin: 38px 0 6px; padding-top: 12px; border-top: 1px solid #c8cdd2;
         font-size: 1.15rem; }}
    p.intro {{ max-width: 96ch; color: #52514e; margin: 6px 0 16px; }}
    figure.card {{ margin: 0 0 22px; background: #fcfcfb; border: 1px solid #e1e0d9;
                  border-radius: 8px; padding: 14px 16px; }}
    figure.card img {{ display: block; width: 100%; height: auto; margin-top: 10px; }}
    figcaption {{ font-size: 0.88rem; color: #52514e; max-width: 110ch; }}
    figcaption strong {{ color: #0b0b0b; }}
  </style>
</head>
<body>
  <header>
    <div class="header-inner">
      <p class="eyebrow">Same-day L2A/L1C exact pairs | 2.04M pixels · 6,858 components · 142 sites</p>
      <h1>Nonlinear L2A harmonizer — model performance</h1>
      <nav aria-label="Related">
        <a href="../../cross_rt_harmonizer_lowcloud152_20260716/summary/visual.html">Ridge visual report</a>
        <a href="../ridge_comparison.json">Ridge comparison JSON</a>
        <a href="../surface_metrics.json">Surface metrics JSON</a>
      </nav>
    </div>
  </header>
  <main>
{chr(10).join(body)}
  </main>
</body>
</html>
"""
    (output / "index.html").write_text(document, encoding="utf-8")
    print(f"wrote {output / 'index.html'}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pairs", type=Path, default=DEFAULT_PAIRS)
    parser.add_argument("--cases", type=Path, default=DEFAULT_CASES)
    parser.add_argument("--output", type=Path, default=HGB_RUN / "report")
    parser.add_argument("--spatial-cases", type=int, default=4)
    args = parser.parse_args()
    style()
    figures_dir = args.output / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)

    hgb_artifact = json.loads((HGB_RUN / "harmonizer.json").read_text(encoding="utf-8"))
    ridge_artifact = json.loads((RIDGE_RUN / "harmonizer.json").read_text(encoding="utf-8"))
    dataset = load_pairs(
        args.pairs,
        _case_ids(args.cases),
        scene_day_max=str(hgb_artifact["training_cutoff"]),
        max_samples_per_scene=int(hgb_artifact["max_samples_per_scene"]),
        allow_missing_matchups=True,
    )
    predictions = oof_predictions(dataset, hgb_artifact, ridge_artifact)
    components = load_joined_components()
    sites = site_table(components)

    section_pixel = (
        "Pixel-level model performance",
        "Out-of-fold predictions recomputed from the saved fold models: every pixel is "
        "scored by a model that never saw its site. The gradient boosting maps the raw "
        "state (both AODs, water vapour, geometry, terrain); the ridge maps hand-crafted "
        "ΔAOT interactions.",
        [
            (
                figure_pred_vs_truth(predictions, figures_dir),
                "Predicted correction vs true residual.",
                "Density of 250k sampled pixels per panel; the diagonal is a perfect "
                "correction. The two models trade off differently: the nonlinear model "
                "is tighter in the core (smaller robust σ) but regresses to the mean in "
                "the tails, while the ridge extrapolates the tails linearly (higher R²). "
                "Neither model resolves the red-band scatter.",
            ),
            (
                figure_pixel_error_distributions(predictions, figures_dir),
                "Pixel error after correction.",
                "Distribution of corrected-minus-target error per pixel (uncapped predictions).",
            ),
        ],
    )
    section_anatomy = (
        "Error anatomy — why the accuracy plateaus",
        "The residual after correction splits into a scene-common part (one offset per "
        "acquisition) and a within-scene part. The scene-common part is the signature "
        "of error in the target itself — a single scene-level MAIAC AOD value and "
        "scene-mean geometry per acquisition — which no per-pixel input can predict. "
        "Both models drive the within-scene part down; the scene-common floor is what "
        "neither can pass.",
        [
            (
                figure_error_decomposition(dataset, predictions, figures_dir),
                "Residual variance decomposition per band (corrections capped ±0.030).",
                "Both models remove most of the predictable variance, and what remains "
                "is dominated by the scene-common component in the aerosol-sensitive "
                "bands — consistent with the ~±0.05 uncertainty of the scene MAIAC AOD "
                "used to build the target (≈0.002–0.004 blue reflectance).",
            ),
            (
                figure_band_scatter_global(dataset, predictions, figures_dir),
                "All-band reflectance before and after, against the target.",
                "150k sampled pixels per panel: original L2A (top row) and corrected "
                "(bottom row) against the same-day L1C surface corrected with MAIAC "
                "AOD. The visible bands tighten onto the diagonal; NIR/SWIR are close "
                "to it before any correction.",
            ),
        ],
    )
    section_component = (
        "Performance across the sample dimensions",
        "Each point/curve summarizes one retained tile-scene component (6,858 total). "
        "Bias is the component-mean error vs the MAIAC-compatible target.",
        [
            (
                figure_component_model_scatter(components, figures_dir),
                "Component |bias|: nonlinear vs ridge, coloured by ΔAOT.",
                "Below the diagonal means the nonlinear model is closer to the target "
                "for that component. Large-|ΔAOT| components (saturated colours) sit "
                "mostly below the diagonal in blue; red-band components sit above.",
            ),
            (
                figure_conditional_performance(components, figures_dir),
                "Blue-band |bias| conditioned on each driver.",
                "Binned medians with IQR bands across MAIAC AOD, ΔAOT, water vapour, "
                "solar zenith, terrain slope, and season, for the uncorrected history "
                "and both models.",
            ),
            (
                figure_component_ecdf(components, figures_dir),
                "Cumulative distribution of component |bias|.",
                "Higher curves are better; the models separate in blue and converge "
                "in red where the correction has little atmospheric signal to use.",
            ),
        ],
    )
    section_site = (
        "Site-level performance",
        "One value per AERONET site (mean over its components).",
        [
            (
                figure_site_map(sites, figures_dir),
                "Where the nonlinear model wins and loses.",
                "Blue-band MAE change vs the ridge at each of the 142 sites; marker "
                "size is the number of components behind the estimate.",
            ),
            (
                figure_site_scatter(sites, figures_dir),
                "Per-site blue MAE, nonlinear vs ridge.",
                "Sites below the diagonal prefer the nonlinear model.",
            ),
        ],
    )
    spatial_figures: list[tuple[str, str, str]] = []
    sample_matchups = matchup_of_sample(dataset)
    example_matchups = available_spatial_cases(components, limit=args.spatial_cases)
    for index, matchup_id in enumerate(example_matchups):
        rendered = figure_spatial_case(matchup_id, components, figures_dir, index)
        if rendered is not None:
            name, caption = rendered
            spatial_figures.append((name, "Applied correction in space.", caption))
        band_scatter = figure_band_scatter_example(
            dataset, predictions, sample_matchups, matchup_id, figures_dir, index
        )
        if band_scatter is not None:
            name, caption = band_scatter
            spatial_figures.append((name, "Per-band reflectance before and after.", caption))
    section_spatial = (
        "Spatial imagery of the applied correction",
        "Monthly composite rasters from the harmonized daily histories (60 m grid, "
        "0.12° AOI). Cases are the highest-|ΔAOT| matchups whose history rasters are "
        "built so far; more appear as the history array completes.",
        spatial_figures if spatial_figures else [],
    )

    sections = [section_pixel, section_anatomy, section_component, section_site]
    if spatial_figures:
        sections.append(section_spatial)
    build_html(args.output, sections)
    print(
        json.dumps(
            {
                "output": str(args.output),
                "pixels": int(dataset.l2a.shape[0]),
                "components": int(len(components)),
                "sites": int(len(sites)),
                "spatial_cases": [name for name, _t, _c in spatial_figures],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
