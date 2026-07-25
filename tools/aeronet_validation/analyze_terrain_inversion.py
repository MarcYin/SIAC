"""Test analytic inversion of Sen2Cor's rugged-terrain correction.

Sen2Cor normalizes slope reflectance by the local irradiance; our flat target
does not, which paints the illumination field into the L2A-minus-target
difference. This fits a two-parameter first-order inversion — band direct-beam
fraction ``c`` and the empirical shade limiter ``g`` —

    rho_flat = rho_L2A * (c * max(cos beta_i, g * cos theta_s) + 1 - c)
                       / (c * cos theta_s + 1 - c)

per band on one rugged scene, and measures how much of the terrain-driven
disagreement it removes out of sample spatially (fit on one scene half, score
on the other).
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
MATCHUP_ID = "University_of_Nizwa__T40QEL_20240105T064251"
HALF_SIZE_DEGREES = 0.06
BANDS = (("blue", 1), ("red", 3), ("swir22", 6))


def invert(l2a: np.ndarray, incidence: np.ndarray, cos_sza: float, c: float, g: float):
    effective = np.maximum(incidence, g * cos_sza)
    return l2a * (c * effective + 1.0 - c) / (c * cos_sza + 1.0 - c)


def main() -> None:
    import csv

    from bestpixel._gee import init_ee, utm_epsg_from_bbox, utm_grid
    from matplotlib.colors import LinearSegmentedColormap

    from siac.algorithms.rt.lut.backend import ZarrLUTBackend

    style()
    OUTPUT.mkdir(parents=True, exist_ok=True)
    diverging = LinearSegmentedColormap.from_list("bgr", [BLUE, "#f0efec", RED])
    coordinates = {
        row["matchup_id"]: (float(row["latitude"]), float(row["longitude"]))
        for row in csv.DictReader(MATCHUPS_CSV.open())
    }
    with np.load(PAIRS / f"{MATCHUP_ID}.npz", allow_pickle=False) as archive:
        scenes = json.loads(str(archive["scenes_json"].item()))
    scene_meta = sorted(scenes, key=lambda s: s["maiac_aot"])[len(scenes) // 2]
    scene = {
        "l1c_id": scene_meta["scene_id"],
        "day": scene_meta["day"],
        "tile": scene_meta["tile"],
        "maiac": scene_meta["maiac_aot"],
        "wvp": scene_meta["maiac_tcwv_cm"],
        "sza": scene_meta["sza_deg"],
        "saa": scene_meta["saa_deg"],
        "vza": scene_meta["vza_deg"],
        "vaa": scene_meta["vaa_deg"],
    }
    lat, lon = coordinates[MATCHUP_ID]
    bbox = (
        lon - HALF_SIZE_DEGREES,
        lat - HALF_SIZE_DEGREES,
        lon + HALF_SIZE_DEGREES,
        lat + HALF_SIZE_DEGREES,
    )
    ee = init_ee()
    grid = utm_grid(bbox, utm_epsg_from_bbox(bbox), 60.0)
    terrain = fetch_glo30_terrain(ee, grid)
    grids, _metadata = _load_pair_grids(
        ee,
        bbox=bbox,
        sidecar_path=Path(scene_meta["sidecar"]),
        scene=scene,
        elevation_km=scene_meta["elevation_km"],
        terrain=terrain,
        target_rt_mode="physical_v2",
        rt_backend=ZarrLUTBackend(TARGET_LUT),
        sensor_cache={},
    )
    valid = np.asarray(grids["valid"], dtype=bool) & terrain.valid
    incidence = local_solar_incidence(
        terrain, sza_deg=scene_meta["sza_deg"], saa_deg=scene_meta["saa_deg"]
    )
    cos_sza = float(np.cos(np.radians(scene_meta["sza_deg"])))
    rows, columns_count = valid.shape
    west = valid & (np.arange(columns_count)[None, :] < columns_count // 2)
    east = valid & ~west

    figure_band = "swir22"
    report = {}
    fig, axes = plt.subplots(1, 4, figsize=(13.0, 3.7), constrained_layout=True)
    for band_name, band_index in BANDS:
        l2a = grids["l2a_surface"][band_index].astype(np.float64)
        target = grids["siac_surface"][band_index].astype(np.float64)
        best = None
        for c in np.linspace(0.5, 1.0, 26):
            for g in np.linspace(0.0, 0.6, 13):
                fitted = invert(l2a[west], incidence[west], cos_sza, c, g)
                score = float(np.nanmean(np.abs(fitted - target[west])))
                if best is None or score < best[0]:
                    best = (score, float(c), float(g))
        _score, c_fit, g_fit = best
        inverted = invert(l2a, incidence, cos_sza, c_fit, g_fit)
        before = l2a - target
        after = inverted - target
        mask = east
        r_before = float(np.corrcoef(incidence[mask], before[mask])[0, 1])
        r_after = float(np.corrcoef(incidence[mask], after[mask])[0, 1])
        report[band_name] = {
            "c_direct_fraction": c_fit,
            "g_shade_limiter": g_fit,
            "holdout_mae_before": float(np.nanmean(np.abs(before[mask]))),
            "holdout_mae_after": float(np.nanmean(np.abs(after[mask]))),
            "holdout_r_incidence_before": r_before,
            "holdout_r_incidence_after": r_after,
        }
        if band_name == figure_band:
            span = {"cmap": diverging, "vmin": -0.08, "vmax": 0.08, "interpolation": "nearest"}
            axes[0].imshow(np.where(valid, before, np.nan), **span)
            axes[0].set_title("before inversion\nL2A − target (B12)", fontsize=9)
            image = axes[1].imshow(np.where(valid, after, np.nan), **span)
            axes[1].set_title(f"after inversion (c={c_fit:.2f}, g={g_fit:.2f})", fontsize=9)
            for ax in axes[:2]:
                ax.set_xticks([])
                ax.set_yticks([])
                ax.grid(False)
                for spine in ax.spines.values():
                    spine.set_visible(False)
            rng = np.random.default_rng(0)
            sample = np.flatnonzero(mask.ravel())
            sample = rng.choice(sample, size=min(25_000, sample.size), replace=False)
            for ax, values, r_value, name in (
                (axes[2], before, r_before, "before"),
                (axes[3], after, r_after, "after"),
            ):
                ax.scatter(
                    incidence.ravel()[sample],
                    values.ravel()[sample],
                    s=3,
                    color=BLUE,
                    alpha=0.15,
                    linewidths=0,
                    rasterized=True,
                )
                ax.axhline(0.0, color=AXIS, linewidth=0.8)
                ax.set_xlim(0.2, 1.05)
                ax.set_ylim(-0.12, 0.12)
                ax.set_xlabel("cos(local solar incidence)", fontsize=8)
                ax.set_title(f"{name}: r = {r_value:+.2f} (holdout half)", fontsize=9)
            axes[2].set_ylabel("L2A − target (B12)", fontsize=8)
    colorbar = fig.colorbar(image, ax=axes[:2], shrink=0.85, pad=0.02)
    colorbar.set_label("reflectance difference", fontsize=8, color=SECONDARY)
    colorbar.outline.set_visible(False)
    fig.suptitle(
        f"Inverting Sen2Cor's terrain correction — {MATCHUP_ID.split('__')[0]} "
        f"{scene_meta['day']} (fit west half, score east half)",
        fontsize=11,
        fontweight="bold",
        color=INK,
    )
    save(fig, OUTPUT, "v05_terrain_inversion.png")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
