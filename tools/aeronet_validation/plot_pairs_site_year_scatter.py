"""Site-by-year scatters of L2A vs the L1C-with-MAIAC surface target.

One grid per band: rows are the highest-support sites, columns are history
years. Each panel shows the operational L2A reflectance (gray) and the
out-of-fold corrected reflectance (blue) against the same-day L1C surface
corrected with MAIAC AOD, exposing how the Sen2Cor-to-current-RT offset and
the learned correction vary across sites and processing years.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import joblib
import matplotlib.pyplot as plt
import numpy as np
from tools.aeronet_validation.plot_cross_rt_harmonization import (
    AXIS,
    BLUE,
    INK,
    MUTED,
    SECONDARY,
    save,
    style,
)
from tools.aeronet_validation.train_l2a_l1c_harmonizer import (
    BAND_NAMES,
    _case_ids,
    load_pairs,
)
from tools.aeronet_validation.train_l2a_l1c_nonlinear_harmonizer import (
    build_columns,
    sen2cor_ozone_node,
)

ROOT = Path("/gws/ssde/j25a/nceo_isp/public/siac_refactor")
V2_RUN = ROOT / "analysis/cross_rt_nonlinear_v2_lowcloud152_20260717"
V2_PAIRS = ROOT / "analysis/l2a_l1c_physical_pairs_v2_lowcloud152_20260717"
DEFAULT_CASES = ROOT / "reports/aod-final-performance-dashboard-20260713/data/all-cases.csv"
PLOT_BANDS = (("blue", "B02 blue", 0.32), ("swir22", "B12 swir22", 0.5))
YEARS = (2019, 2020, 2021, 2022, 2023)
SITE_COUNT = 6
CAP = 0.03


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=V2_RUN / "report_figures")
    args = parser.parse_args()
    style()
    args.output.mkdir(parents=True, exist_ok=True)

    artifact = json.loads((V2_RUN / "harmonizer.json").read_text(encoding="utf-8"))
    dataset = load_pairs(
        V2_PAIRS,
        _case_ids(DEFAULT_CASES),
        scene_day_max=str(artifact["training_cutoff"]),
        max_samples_per_scene=int(artifact["max_samples_per_scene"]),
        allow_missing_matchups=True,
    )
    columns = build_columns(dataset)
    cams = np.empty(len(dataset.scene_metadata))
    lookup = {
        (row.split(",")[0], row.split(",")[1]): float(row.split(",")[3])
        for row in (V2_PAIRS / "cams_scene_state.csv").read_text().splitlines()[1:]
        if row.strip()
    }
    scene_code = np.asarray(dataset.scene_code, dtype=np.int64)
    for index, scene in enumerate(dataset.scene_metadata):
        cams[index] = lookup[(str(scene["matchup_id"]), str(scene["day"]))]
    columns["ozone_du_cams"] = cams[scene_code]
    columns["ozone_du_sen2cor_node"] = sen2cor_ozone_node(columns["ozone_du_cams"])

    folds = np.full(dataset.site_code.shape, -1, dtype=np.int16)
    for site_code, matchup_ids in dataset.matchup_by_site_code.items():
        folds[dataset.site_code == site_code] = int(artifact["fold_by_matchup_id"][matchup_ids[0]])

    sites = np.asarray(
        [str(scene["matchup_id"]).split("__", 1)[0] for scene in dataset.scene_metadata],
        dtype=object,
    )[scene_code]
    years = np.asarray(
        [int(str(scene["day"])[:4]) for scene in dataset.scene_metadata], dtype=np.int64
    )[scene_code]

    support: dict[str, tuple[int, int]] = {}
    for site in np.unique(sites):
        mask = sites == site
        support[str(site)] = (int(len(set(years[mask]) & set(YEARS))), int(mask.sum()))
    chosen = [site for site, _ in sorted(support.items(), key=lambda kv: (-kv[1][0], -kv[1][1]))][
        :SITE_COUNT
    ]

    for band, label, upper in PLOT_BANDS:
        band_index = BAND_NAMES.index(band)
        names = artifact["feature_names_by_band"][band]
        features = np.column_stack([columns[name] for name in names])
        predicted = np.full(features.shape[0], np.nan)
        for fold in np.unique(folds):
            mask = folds == fold
            estimator = joblib.load(V2_RUN / "models" / artifact["model_files"][str(fold)][band])
            predicted[mask] = estimator.predict(features[mask])
        target = dataset.siac[:, band_index].astype(np.float64)
        before = dataset.l2a[:, band_index].astype(np.float64)
        after = before + np.clip(predicted, -CAP, CAP)

        fig, axes = plt.subplots(
            len(chosen),
            len(YEARS),
            figsize=(12.4, 2.15 * len(chosen) + 0.8),
            sharex=True,
            sharey=True,
        )
        for row, site in enumerate(chosen):
            for column, year in enumerate(YEARS):
                ax = axes[row, column]
                mask = (sites == site) & (years == year)
                count = int(mask.sum())
                ax.plot([0, upper], [0, upper], color=AXIS, linewidth=0.7)
                if count:
                    ax.scatter(
                        target[mask],
                        before[mask],
                        s=3,
                        color=MUTED,
                        alpha=0.35,
                        linewidths=0,
                        rasterized=True,
                    )
                    ax.scatter(
                        target[mask],
                        after[mask],
                        s=3,
                        color=BLUE,
                        alpha=0.35,
                        linewidths=0,
                        rasterized=True,
                    )
                    mae_before = float(np.nanmean(np.abs(before[mask] - target[mask])))
                    mae_after = float(np.nanmean(np.abs(after[mask] - target[mask])))
                    ax.text(
                        0.04,
                        0.95,
                        f"{mae_before:.4f}→{mae_after:.4f}",
                        transform=ax.transAxes,
                        ha="left",
                        va="top",
                        fontsize=6.5,
                        color=SECONDARY,
                    )
                    ax.text(
                        0.96,
                        0.04,
                        f"n={count}",
                        transform=ax.transAxes,
                        ha="right",
                        va="bottom",
                        fontsize=6,
                        color=MUTED,
                    )
                else:
                    ax.text(
                        0.5,
                        0.5,
                        "no pairs",
                        transform=ax.transAxes,
                        ha="center",
                        va="center",
                        fontsize=7,
                        color=MUTED,
                    )
                ax.set_xlim(0, upper)
                ax.set_ylim(0, upper)
                ax.set_xticks([0, upper / 2, upper])
                ax.set_yticks([0, upper / 2, upper])
                ax.tick_params(labelsize=6)
                if row == 0:
                    ax.set_title(str(year), fontsize=9)
                if column == 0:
                    ax.set_ylabel(site, fontsize=7.5, color=INK)
        for ax in axes[-1, :]:
            ax.set_xlabel("target reflectance", fontsize=7)
        handles = [
            plt.Line2D([], [], marker="o", linestyle="", color=MUTED, label="L2A before"),
            plt.Line2D([], [], marker="o", linestyle="", color=BLUE, label="corrected (OOF)"),
        ]
        fig.legend(handles=handles, loc="outside upper right", ncols=2, fontsize=8)
        fig.suptitle(
            f"{label} — L2A vs L1C corrected with MAIAC AOD, by site and year "
            f"(MAE before→after per panel)",
            fontsize=11,
            fontweight="bold",
            color=INK,
        )
        fig.set_constrained_layout(True)
        save(fig, args.output, f"v03_site_year_{band}.png")


if __name__ == "__main__":
    main()
