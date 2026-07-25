"""Validate the fully analytic de-terrain factor (no fitted parameters).

For one rugged scene per site, computes the per-band direct-beam fraction from
the LUT at the Sen2Cor scene state, applies ``deterrain_factor`` with the
configured shade limiter, and reports how much of the terrain-driven
L2A-vs-target disagreement is removed — the go/no-go for wiring the step into
the pair and history builders.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from tools.aeronet_validation.build_l2a_l1c_harmonization_pairs import (
    TARGET_LUT,
    _load_pair_grids,
    _sentinel_satellite_id,
)
from tools.aeronet_validation.terrain_deshade import (
    BRDF_EXPONENT_NIR_SWIR,
    BRDF_EXPONENT_VIS,
    G_SHADE_LIMITER,
    NIR_SWIR_BANDS,
    DirectFractionLUT,
    deterrain_factor,
)
from tools.aeronet_validation.terrain_features import fetch_glo30_terrain, local_solar_incidence

ROOT = Path("/gws/ssde/j25a/nceo_isp/public/siac_refactor")
PAIRS = ROOT / "analysis/l2a_l1c_physical_pairs_v2_lowcloud152_20260717"
MATCHUPS_CSV = ROOT / "matchups/matchups.csv"
MATCHUP_IDS = (
    "University_of_Nizwa__T40QEL_20240105T064251",
    "ICIMOD__T45RUL_20240107T045159",
    "Huancayo-IGP__T18LVM_20240926T151719",
)
HALF_SIZE_DEGREES = 0.06
BANDS = (("blue", "B02", 1), ("red", "B04", 3), ("swir22", "B12", 6))


def main() -> None:
    import csv

    from bestpixel._gee import init_ee, utm_epsg_from_bbox, utm_grid

    from siac.adapters.rsrf import load_sensor_config_with_rsrf
    from siac.algorithms.rt.lut.backend import ZarrLUTBackend

    coordinates = {
        row["matchup_id"]: (float(row["latitude"]), float(row["longitude"]))
        for row in csv.DictReader(MATCHUPS_CSV.open())
    }
    ee = init_ee()
    rt_backend = ZarrLUTBackend(TARGET_LUT)
    fraction_lut = DirectFractionLUT(TARGET_LUT)
    for matchup_id in MATCHUP_IDS:
        with np.load(PAIRS / f"{matchup_id}.npz", allow_pickle=False) as archive:
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
            scene=scene,
            elevation_km=scene_meta["elevation_km"],
            terrain=terrain,
            target_rt_mode="physical_v2",
            rt_backend=rt_backend,
            sensor_cache={},
        )
        satellite = _sentinel_satellite_id((scene_meta.get("l2a") or {}).get("spacecraft"))
        sensor = load_sensor_config_with_rsrf("MSI", satellite)
        bands = [sensor.get_band(s2_name) for _key, s2_name, _index in BANDS]
        valid = np.asarray(grids["valid"], dtype=bool) & terrain.valid
        sen2cor_state = {
            "aot": float(np.nanmedian(grids["l2a_aot"][valid])),
            "tcwv_cm": float(np.nanmedian(grids["l2a_tcwv"][valid])),
            "tco3_atm_cm": float(metadata.get("tco3_atm_cm") or 0.30),
        }
        fractions = fraction_lut.band_direct_fractions(
            bands,
            aot=sen2cor_state["aot"],
            sza_deg=scene_meta["sza_deg"],
            tcwv_cm=sen2cor_state["tcwv_cm"],
            tco3_atm_cm=sen2cor_state["tco3_atm_cm"],
            elevation_km=scene_meta["elevation_km"],
        )
        incidence = local_solar_incidence(
            terrain, sza_deg=scene_meta["sza_deg"], saa_deg=scene_meta["saa_deg"]
        )
        cos_sza = float(np.cos(np.radians(scene_meta["sza_deg"])))
        report = {
            "site": matchup_id.split("__", 1)[0],
            "day": scene_meta["day"],
            "sen2cor_state": sen2cor_state,
            "g_limiter": G_SHADE_LIMITER,
        }
        for band_key, s2_name, band_index in BANDS:
            l2a = grids["l2a_surface"][band_index].astype(np.float64)[valid]
            target = grids["siac_surface"][band_index].astype(np.float64)[valid]
            inc = incidence[valid]
            slope = terrain.slope_deg[valid]
            rho_background = float(np.nanmedian(l2a))
            factor = deterrain_factor(
                fractions[s2_name],
                inc,
                slope,
                cos_sza,
                rho_background=rho_background,
                exponent=BRDF_EXPONENT_NIR_SWIR
                if band_key in NIR_SWIR_BANDS
                else BRDF_EXPONENT_VIS,
            )
            before = l2a - target
            after = l2a * factor - target
            report[band_key] = {
                "c_from_lut": round(fractions[s2_name], 3),
                "mae_before": round(float(np.nanmean(np.abs(before))), 5),
                "mae_after": round(float(np.nanmean(np.abs(after))), 5),
                "r_incidence_before": round(float(np.corrcoef(inc, before)[0, 1]), 3),
                "r_incidence_after": round(float(np.corrcoef(inc, after)[0, 1]), 3),
            }
        print(json.dumps(report), flush=True)


if __name__ == "__main__":
    main()
