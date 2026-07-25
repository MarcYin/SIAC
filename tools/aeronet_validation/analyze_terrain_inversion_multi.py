"""Stability of the terrain-inversion parameters across scenes and sites.

Fits the two-parameter de-terrain transform per band on many rugged scenes
spanning sites, seasons and AOT, then examines whether the direct-beam
fraction ``c`` varies with the atmospheric state (as physics predicts) and
whether the shade limiter ``g`` behaves as a processor constant.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from tools.aeronet_validation.analyze_terrain_inversion import invert
from tools.aeronet_validation.build_l2a_l1c_harmonization_pairs import (
    TARGET_LUT,
    _load_pair_grids,
)
from tools.aeronet_validation.plot_cross_rt_harmonization import (
    BLUE,
    INK,
    MUTED,
    save,
    style,
)
from tools.aeronet_validation.terrain_features import fetch_glo30_terrain, local_solar_incidence

ROOT = Path("/gws/ssde/j25a/nceo_isp/public/siac_refactor")
PAIRS = ROOT / "analysis/l2a_l1c_physical_pairs_v2_lowcloud152_20260717"
OUTPUT = ROOT / "analysis/cross_rt_nonlinear_v2_lowcloud152_20260717/report_figures"
MATCHUPS_CSV = ROOT / "matchups/matchups.csv"
SITES = (
    ("University_of_Nizwa__T40QEL_20240105T064251", 12),
    ("ICIMOD__T45RUL_20240107T045159", 6),
    ("Huancayo-IGP__T18LVM_20240926T151719", 6),
)
HALF_SIZE_DEGREES = 0.06
BANDS = (("blue", 1), ("red", 3), ("swir22", 6))
BAND_COLORS = {"blue": BLUE, "red": "#e34948", "swir22": "#008300"}
SITE_MARKERS = {"University_of_Nizwa": "o", "ICIMOD": "s", "Huancayo-IGP": "^"}


def main() -> None:
    import csv

    from bestpixel._gee import init_ee, utm_epsg_from_bbox, utm_grid

    from siac.algorithms.rt.lut.backend import ZarrLUTBackend

    style()
    OUTPUT.mkdir(parents=True, exist_ok=True)
    coordinates = {
        row["matchup_id"]: (float(row["latitude"]), float(row["longitude"]))
        for row in csv.DictReader(MATCHUPS_CSV.open())
    }
    ee = init_ee()
    rt_backend = ZarrLUTBackend(TARGET_LUT)
    records: list[dict] = []
    for matchup_id, scene_count in SITES:
        site = matchup_id.split("__", 1)[0]
        with np.load(PAIRS / f"{matchup_id}.npz", allow_pickle=False) as archive:
            scenes = json.loads(str(archive["scenes_json"].item()))
        # spread the selection across the AOT range deterministically
        ordered = sorted(scenes, key=lambda s: s["maiac_aot"])
        picks = [
            ordered[round(i * (len(ordered) - 1) / max(scene_count - 1, 1))]
            for i in range(min(scene_count, len(ordered)))
        ]
        lat, lon = coordinates[matchup_id]
        bbox = (
            lon - HALF_SIZE_DEGREES,
            lat - HALF_SIZE_DEGREES,
            lon + HALF_SIZE_DEGREES,
            lat + HALF_SIZE_DEGREES,
        )
        grid = utm_grid(bbox, utm_epsg_from_bbox(bbox), 60.0)
        terrain = fetch_glo30_terrain(ee, grid)
        sensor_cache: dict = {}
        for scene_meta in picks:
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
            try:
                grids, _meta = _load_pair_grids(
                    ee,
                    bbox=bbox,
                    sidecar_path=Path(scene_meta["sidecar"]),
                    scene=scene,
                    elevation_km=scene_meta["elevation_km"],
                    terrain=terrain,
                    target_rt_mode="physical_v2",
                    rt_backend=rt_backend,
                    sensor_cache=sensor_cache,
                )
            except Exception as exc:  # noqa: BLE001 - skip transient scene failures
                print(f"skip {site} {scene_meta['day']}: {exc}", flush=True)
                continue
            valid = np.asarray(grids["valid"], dtype=bool) & terrain.valid
            if int(valid.sum()) < 3000:
                continue
            incidence = local_solar_incidence(
                terrain, sza_deg=scene_meta["sza_deg"], saa_deg=scene_meta["saa_deg"]
            )
            cos_sza = float(np.cos(np.radians(scene_meta["sza_deg"])))
            record = {
                "site": site,
                "day": scene_meta["day"],
                "maiac_aot": float(scene_meta["maiac_aot"]),
                "sza_deg": float(scene_meta["sza_deg"]),
            }
            for band_name, band_index in BANDS:
                l2a = grids["l2a_surface"][band_index].astype(np.float64)[valid]
                target = grids["siac_surface"][band_index].astype(np.float64)[valid]
                inc = incidence[valid]
                best = None
                for c in np.linspace(0.4, 1.0, 31):
                    for g in np.linspace(0.0, 0.6, 13):
                        score = float(np.nanmean(np.abs(invert(l2a, inc, cos_sza, c, g) - target)))
                        if best is None or score < best[0]:
                            best = (score, float(c), float(g))
                _s, c_fit, g_fit = best
                r_before = float(np.corrcoef(inc, l2a - target)[0, 1])
                record[f"{band_name}_c"] = c_fit
                record[f"{band_name}_g"] = g_fit
                record[f"{band_name}_r_before"] = r_before
            records.append(record)
            print(json.dumps(record), flush=True)

    fig, axes = plt.subplots(1, 3, figsize=(12.6, 4.0), constrained_layout=True)
    axes_map = {"blue": axes[0], "swir22": axes[1]}
    for band_name in ("blue", "swir22"):
        ax = axes_map[band_name]
        for record in records:
            ax.scatter(
                record["maiac_aot"],
                record[f"{band_name}_c"],
                s=46,
                color=BAND_COLORS[band_name],
                marker=SITE_MARKERS[record["site"]],
                linewidths=1.0,
                edgecolors="#fcfcfb",
            )
        ax.set_xlabel("scene MAIAC AOT")
        ax.set_ylabel(f"fitted c — {band_name}")
        ax.set_ylim(0.35, 1.05)
        ax.set_title(f"direct-beam fraction, {band_name}", fontsize=10)
    ax = axes[2]
    for band_name, _index in BANDS:
        values = [record[f"{band_name}_g"] for record in records]
        ax.hist(
            values,
            bins=np.linspace(-0.025, 0.625, 14),
            histtype="step",
            linewidth=2,
            color=BAND_COLORS[band_name],
            label=band_name,
        )
    ax.set_xlabel("fitted g (shade limiter)")
    ax.set_ylabel("scenes")
    ax.set_title("g across scenes", fontsize=10)
    ax.legend(fontsize=8)
    handles = [
        plt.Line2D([], [], marker=marker, linestyle="", color=MUTED, label=site)
        for site, marker in SITE_MARKERS.items()
    ]
    fig.legend(handles=handles, loc="outside upper right", ncols=3, fontsize=8)
    fig.suptitle(
        "Terrain-inversion parameters across scenes and sites",
        fontsize=11,
        fontweight="bold",
        color=INK,
    )
    save(fig, OUTPUT, "v06_terrain_inversion_stability.png")
    (OUTPUT / "terrain_inversion_parameters.json").write_text(
        json.dumps(records, indent=2) + "\n", encoding="utf-8"
    )
    for band_name, _index in BANDS:
        c_values = np.asarray([record[f"{band_name}_c"] for record in records])
        aot = np.asarray([record["maiac_aot"] for record in records])
        print(
            f"{band_name}: c median {np.median(c_values):.2f} "
            f"IQR {np.percentile(c_values, 25):.2f}-{np.percentile(c_values, 75):.2f} "
            f"corr(c, AOT) {np.corrcoef(aot, c_values)[0, 1]:+.2f}"
        )


if __name__ == "__main__":
    main()
