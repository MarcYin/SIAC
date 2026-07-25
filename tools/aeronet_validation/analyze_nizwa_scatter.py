"""Spatial anatomy of the University_of_Nizwa pair scatter.

Fetches full same-acquisition L2A/target grids for one Nizwa scene and one
Lahore control, maps the per-pixel disagreement against the local solar
illumination field, and estimates any residual co-registration shift.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from tools.aeronet_validation.build_l2a_l1c_harmonization_pairs import (
    TARGET_LUT,
    _load_pair_grids,
)
from tools.aeronet_validation.plot_cross_rt_harmonization import (
    AXIS,
    BLUE,
    INK,
    RED,
    SECONDARY,
    save,
    style,
)
from tools.aeronet_validation.terrain_features import fetch_glo30_terrain, local_solar_incidence

ROOT = Path("/gws/ssde/j25a/nceo_isp/public/siac_refactor")
PAIRS = ROOT / "analysis/l2a_l1c_physical_pairs_v2_lowcloud152_20260717"
OUTPUT = ROOT / "analysis/cross_rt_nonlinear_v2_lowcloud152_20260717/report_figures"
MATCHUPS_CSV = ROOT / "matchups/matchups.csv"
DIVERGING = None  # built in main from palette poles
CASES = (
    ("University_of_Nizwa__T40QEL_20240105T064251", None),
    ("Lahore__T43RDQ_20240130T055101", None),
)
HALF_SIZE_DEGREES = 0.06


def scene_input(metadata: dict) -> dict:
    return {
        "l1c_id": metadata["scene_id"],
        "day": metadata["day"],
        "tile": metadata["tile"],
        "maiac": metadata["maiac_aot"],
        "wvp": metadata["maiac_tcwv_cm"],
        "sza": metadata["sza_deg"],
        "saa": metadata["saa_deg"],
        "vza": metadata["vza_deg"],
        "vaa": metadata["vaa_deg"],
    }


def main() -> None:
    import csv

    from bestpixel._gee import init_ee, utm_epsg_from_bbox, utm_grid
    from matplotlib.colors import LinearSegmentedColormap
    from skimage.registration import phase_cross_correlation

    from siac.algorithms.rt.lut.backend import ZarrLUTBackend

    style()
    OUTPUT.mkdir(parents=True, exist_ok=True)
    diverging = LinearSegmentedColormap.from_list("bgr", [BLUE, "#f0efec", RED])
    coordinates = {
        row["matchup_id"]: (float(row["latitude"]), float(row["longitude"]))
        for row in csv.DictReader(MATCHUPS_CSV.open())
    }
    ee = init_ee()
    rt_backend = ZarrLUTBackend(TARGET_LUT)
    fig, axes = plt.subplots(2, 4, figsize=(13.0, 7.2), constrained_layout=True)
    for row_index, (matchup_id, _unused) in enumerate(CASES):
        with np.load(PAIRS / f"{matchup_id}.npz", allow_pickle=False) as archive:
            scenes = json.loads(str(archive["scenes_json"].item()))
        # median-AOD scene keeps the atmospheric term ordinary.
        scene_meta = sorted(scenes, key=lambda s: s["maiac_aot"])[len(scenes) // 2]
        lat, lon = coordinates[matchup_id]
        bbox = (
            lon - HALF_SIZE_DEGREES,
            lat - HALF_SIZE_DEGREES,
            lon + HALF_SIZE_DEGREES,
            lat + HALF_SIZE_DEGREES,
        )
        grid = utm_grid(bbox, utm_epsg_from_bbox(bbox), 60.0)
        terrain = fetch_glo30_terrain(ee, grid)
        grids, metadata = _load_pair_grids(
            ee,
            bbox=bbox,
            sidecar_path=Path(scene_meta["sidecar"]),
            scene=scene_input(scene_meta),
            elevation_km=scene_meta["elevation_km"],
            terrain=terrain,
            target_rt_mode="physical_v2",
            rt_backend=rt_backend,
            sensor_cache={},
        )
        band = 6  # swir22
        valid = np.asarray(grids["valid"], dtype=bool) & terrain.valid
        l2a = np.where(valid, grids["l2a_surface"][band], np.nan)
        target = np.where(valid, grids["siac_surface"][band], np.nan)
        difference = l2a - target
        incidence = local_solar_incidence(
            terrain, sza_deg=scene_meta["sza_deg"], saa_deg=scene_meta["saa_deg"]
        )
        incidence = np.where(valid, incidence, np.nan)

        filled_l2a = np.nan_to_num(l2a, nan=float(np.nanmedian(l2a)))
        filled_target = np.nan_to_num(target, nan=float(np.nanmedian(target)))
        shift, _error, _phase = phase_cross_correlation(
            filled_target, filled_l2a, upsample_factor=20
        )
        finite = np.isfinite(difference) & np.isfinite(incidence)
        correlation = float(
            np.corrcoef(incidence[finite].ravel(), difference[finite].ravel())[0, 1]
        )

        site = matchup_id.split("__", 1)[0]
        panels = axes[row_index]
        panels[0].imshow(target, cmap="gray", vmin=0.0, vmax=0.45, interpolation="nearest")
        panels[0].set_title(f"{site} — target B12\n{scene_meta['day']}", fontsize=9)
        image = panels[1].imshow(
            difference, cmap=diverging, vmin=-0.08, vmax=0.08, interpolation="nearest"
        )
        panels[1].set_title(
            f"L2A − target (B12)\nshift ({shift[0]:+.2f}, {shift[1]:+.2f}) px", fontsize=9
        )
        panels[2].imshow(incidence, cmap="gray", vmin=0.0, vmax=1.0, interpolation="nearest")
        panels[2].set_title("cos(local solar incidence)", fontsize=9)
        for ax in panels[:3]:
            ax.set_xticks([])
            ax.set_yticks([])
            ax.grid(False)
            for spine in ax.spines.values():
                spine.set_visible(False)
        scatter_ax = panels[3]
        sample = np.flatnonzero(finite.ravel())
        rng = np.random.default_rng(0)
        sample = rng.choice(sample, size=min(30_000, sample.size), replace=False)
        scatter_ax.scatter(
            incidence.ravel()[sample],
            difference.ravel()[sample],
            s=3,
            color=BLUE,
            alpha=0.15,
            linewidths=0,
            rasterized=True,
        )
        scatter_ax.axhline(0.0, color=AXIS, linewidth=0.8)
        scatter_ax.set_xlim(0.2, 1.05)
        scatter_ax.set_ylim(-0.12, 0.12)
        scatter_ax.set_xlabel("cos(local solar incidence)", fontsize=8)
        scatter_ax.set_ylabel("L2A − target (B12)", fontsize=8)
        scatter_ax.set_title(f"r = {correlation:+.2f}", fontsize=9)
        print(
            json.dumps(
                {
                    "site": site,
                    "day": scene_meta["day"],
                    "maiac_aot": scene_meta["maiac_aot"],
                    "registration_shift_px": [float(shift[0]), float(shift[1])],
                    "incidence_diff_correlation": correlation,
                    "mean_abs_diff_swir22": float(np.nanmean(np.abs(difference))),
                }
            ),
            flush=True,
        )
    colorbar = fig.colorbar(image, ax=axes[:, 1], shrink=0.7, pad=0.01)
    colorbar.set_label("reflectance difference", fontsize=8, color=SECONDARY)
    colorbar.outline.set_visible(False)
    fig.suptitle(
        "Same-acquisition L2A vs target disagreement against the illumination field",
        fontsize=12,
        fontweight="bold",
        color=INK,
    )
    save(fig, OUTPUT, "v04_nizwa_scatter_anatomy.png")


if __name__ == "__main__":
    main()
