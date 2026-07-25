"""Retrieved-vs-truth scatters for the ozone-consistency retrieval ladder."""

from __future__ import annotations

import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from tools.aeronet_validation.plot_cross_rt_harmonization import (
    AXIS,
    BLUE,
    INK,
    MUTED,
    RED,
    SECONDARY,
    SURFACE,
    save,
    style,
)

ROOT = Path("/gws/ssde/j25a/nceo_isp/public/siac_refactor")
OUTPUT = ROOT / "analysis/cross_rt_nonlinear_v2_lowcloud152_20260717/report_figures"
PANELS = (
    (
        "identity · solver O3 0.3",
        "phaseD_results_lowcloud20_crossrt_identity_daily_physical_20260716",
    ),
    (
        "identity · solver O3 CAMS",
        "phaseD_results_lowcloud20_crossrt_identity_daily_physical_o3v2_20260717",
    ),
    (
        "blue-only + O3 feature (best)",
        "phaseD_results_lowcloud20_crossrt_hgb_target_band_cams_o3_blue_cap0p030_physical_20260716",
    ),
    (
        "blue-only · fully O3-consistent",
        "phaseD_results_lowcloud20_crossrt_hgb_target_band_cams_o3_blue_cap0p030_physical_o3v2_20260717",
    ),
    (
        "solver-bands · fully O3-consistent",
        "phaseD_results_lowcloud20_crossrt_hgb_target_band_cams_o3_solver_cap0p030_physical_o3v2_20260717",
    ),
    ("previous best (control)", "phaseD_results_lowcloud20_geometry_backstop05_b03_chi2_20260712"),
)


def load(directory: str, cases: list[str]) -> list[dict]:
    rows = []
    for matchup_id in cases:
        path = ROOT / directory / f"{matchup_id}.json"
        payload = json.loads(path.read_text()) if path.exists() else {"status": "MISSING"}
        payload["matchup_id"] = matchup_id
        rows.append(payload)
    return rows


def valid(rows: list[dict]) -> list[dict]:
    return [
        row
        for row in rows
        if row.get("status") == "OK"
        and isinstance(row.get("truth"), float)
        and isinstance(row.get("retrieved"), float)
        and math.isfinite(row["truth"])
        and math.isfinite(row["retrieved"])
    ]


def main() -> None:
    style()
    OUTPUT.mkdir(parents=True, exist_ok=True)
    cases = [
        line.split(",")[0]
        for line in (ROOT / "reports/aod-final-performance-dashboard-20260713/data/all-cases.csv")
        .read_text()
        .splitlines()[1:]
        if line.strip()
    ]
    limits = (0.008, 3.0)
    grid = np.geomspace(*limits, 200)
    envelope = 0.05 + 0.15 * grid

    fig, axes = plt.subplots(2, 3, figsize=(11.8, 8.0), sharex=True, sharey=True)
    for ax, (title, directory) in zip(axes.ravel(), PANELS, strict=True):
        rows = load(directory, cases)
        usable = valid(rows)
        hits = sum(1 for row in usable if row.get("within_ee"))
        ax.fill_between(
            grid,
            np.clip(grid - envelope, limits[0], None),
            grid + envelope,
            color=MUTED,
            alpha=0.12,
            linewidth=0,
        )
        ax.plot(grid, grid, color=AXIS, linewidth=0.8)
        common = {"linewidths": 1.0, "edgecolors": SURFACE, "s": 24}
        for subset, color, label in (
            ([row for row in usable if row.get("within_ee")], BLUE, "within EE"),
            ([row for row in usable if not row.get("within_ee")], RED, "outside EE"),
        ):
            ax.scatter(
                np.clip([row["truth"] for row in subset], *limits),
                np.clip([row["retrieved"] for row in subset], *limits),
                color=color,
                label=label,
                **common,
            )
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(*limits)
        ax.set_ylim(*limits)
        ax.set_title(f"{title}\n{100.0 * hits / len(rows):.1f}% ({hits}/{len(rows)})", fontsize=9.5)
    for ax in axes[1, :]:
        ax.set_xlabel("AERONET AOD (truth)")
    for ax in axes[:, 0]:
        ax.set_ylabel("retrieved AOD")
    axes[0, 0].legend(loc="upper left", fontsize=8)
    fig.set_constrained_layout(True)
    save(fig, OUTPUT, "v01_retrieval_scatter_ladder.png")

    best = {row["matchup_id"]: row for row in valid(load(PANELS[2][1], cases))}
    consistent = {row["matchup_id"]: row for row in valid(load(PANELS[3][1], cases))}
    shared = sorted(set(best) & set(consistent))
    fig, ax = plt.subplots(figsize=(6.6, 6.2), constrained_layout=True)
    ax.plot(limits, limits, color=AXIS, linewidth=0.8)
    changed = [
        m for m in shared if bool(best[m].get("within_ee")) != bool(consistent[m].get("within_ee"))
    ]
    stable = [m for m in shared if m not in set(changed)]
    ax.scatter(
        np.clip([best[m]["retrieved"] for m in stable], *limits),
        np.clip([consistent[m]["retrieved"] for m in stable], *limits),
        s=22,
        color=BLUE,
        alpha=0.55,
        linewidths=1.0,
        edgecolors=SURFACE,
        label="EE verdict unchanged",
    )
    ax.scatter(
        np.clip([best[m]["retrieved"] for m in changed], *limits),
        np.clip([consistent[m]["retrieved"] for m in changed], *limits),
        s=48,
        color=RED,
        linewidths=1.2,
        edgecolors=SURFACE,
        label="EE verdict flipped",
        zorder=3,
    )
    for m in changed:
        ax.annotate(
            m.split("__", 1)[0],
            (best[m]["retrieved"], consistent[m]["retrieved"]),
            textcoords="offset points",
            xytext=(7, 4),
            fontsize=8,
            color=SECONDARY,
        )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(*limits)
    ax.set_ylim(*limits)
    ax.set_xlabel("retrieved AOD — best recipe (solver O3 0.3)")
    ax.set_ylabel("retrieved AOD — fully O3-consistent")
    ax.set_title(
        f"Same recipe, two ozone frames — {len(changed)} of {len(shared)} verdicts flip",
        fontsize=10.5,
        color=INK,
    )
    ax.legend(loc="upper left", fontsize=8.5)
    save(fig, OUTPUT, "v02_frame_comparison.png")


if __name__ == "__main__":
    main()
