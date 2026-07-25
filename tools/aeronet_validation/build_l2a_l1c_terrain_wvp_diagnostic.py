"""Diagnose terrain and water-vapour differences in exact L2A/L1C pairs.

This is an attribution diagnostic for the selected spatial pair gallery. It
does not alter the harmonizer or retrieval. The L1C target is evaluated once
with its saved scalar scene water vapour and once with the matching L2A WVP
map, while terrain fields are sampled independently from Copernicus GLO-30.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any

import numpy as np
from tools.aeronet_validation.build_l2a_l1c_harmonization_pairs import (
    CANONICAL_BANDS,
    ELEVATION_CACHE,
    HALF_SIZE_DEGREES,
    MATCHUPS,
    ROOT,
    _load_pair_grids,
    _site_elevation_km,
)
from tools.aeronet_validation.build_l2a_l1c_pair_gallery import (
    EXAMPLES,
    _pair_scene_input,
    _scene_from_archive,
)
from tools.aeronet_validation.terrain_features import (
    fetch_glo30_terrain,
    local_solar_incidence,
)

DEFAULT_OUTPUT = ROOT / "reports/l2a-harmonization-20260714/assets/pair-examples"
VISIBLE_BANDS = (1, 2, 3)


def _correlation(left: np.ndarray, right: np.ndarray, valid: np.ndarray) -> float | None:
    mask = valid & np.isfinite(left) & np.isfinite(right)
    if int(np.count_nonzero(mask)) < 20:
        return None
    x = np.asarray(left[mask], dtype=np.float64)
    y = np.asarray(right[mask], dtype=np.float64)
    x -= np.mean(x)
    y -= np.mean(y)
    denominator = np.sqrt(np.dot(x, x) * np.dot(y, y))
    if not np.isfinite(denominator) or denominator <= 0.0:
        return None
    return float(np.dot(x, y) / denominator)


def _zscore(values: np.ndarray) -> np.ndarray:
    mean = np.mean(values)
    scale = np.std(values)
    if not np.isfinite(scale) or scale <= 1.0e-12:
        return np.zeros(values.shape, dtype=np.float64)
    return (values - mean) / scale


def _scene_model(
    residual: np.ndarray,
    fields: dict[str, np.ndarray],
    valid: np.ndarray,
) -> dict[str, Any]:
    keys = tuple(fields)
    mask = valid & np.isfinite(residual)
    for value in fields.values():
        mask &= np.isfinite(value)
    if int(np.count_nonzero(mask)) < 100:
        return {"n": int(np.count_nonzero(mask)), "r2": None, "coefficients": {}}
    y = _zscore(np.asarray(residual[mask], dtype=np.float64))
    x = np.column_stack(
        [np.ones(y.size, dtype=np.float64)]
        + [_zscore(np.asarray(fields[key][mask], dtype=np.float64)) for key in keys]
    )
    coefficients, _residuals, _rank, _singular = np.linalg.lstsq(x, y, rcond=None)
    fitted = x @ coefficients
    total = float(np.dot(y, y))
    r2 = None if total <= 0.0 else float(1.0 - np.dot(y - fitted, y - fitted) / total)
    return {
        "n": int(y.size),
        "r2": r2,
        "coefficients": {key: float(coefficients[index + 1]) for index, key in enumerate(keys)},
    }


def _pooled_model(
    values: list[tuple[np.ndarray, dict[str, np.ndarray], np.ndarray]]
) -> dict[str, Any]:
    keys = ("elevation", "slope", "incidence", "wvp_departure", "target_visible")
    xs: list[np.ndarray] = []
    ys: list[np.ndarray] = []
    for residual, fields, valid in values:
        mask = valid & np.isfinite(residual)
        for key in keys:
            mask &= np.isfinite(fields[key])
        selected = np.flatnonzero(mask)
        if selected.size > 12_000:
            selected = np.linspace(0, selected.size - 1, 12_000, dtype=np.int64)
            selected = np.flatnonzero(mask)[selected]
        if selected.size < 100:
            continue
        ys.append(_zscore(np.asarray(residual.ravel()[selected], dtype=np.float64)))
        xs.append(
            np.column_stack(
                [
                    _zscore(np.asarray(fields[key].ravel()[selected], dtype=np.float64))
                    for key in keys
                ]
            )
        )
    if not xs:
        return {"n": 0, "r2": None, "coefficients": {}}
    x_values = np.concatenate(xs, axis=0)
    y_values = np.concatenate(ys, axis=0)
    design = np.column_stack((np.ones(y_values.size, dtype=np.float64), x_values))
    coefficients, _residuals, _rank, _singular = np.linalg.lstsq(design, y_values, rcond=None)
    fitted = design @ coefficients
    total = float(np.dot(y_values, y_values))
    return {
        "n": int(y_values.size),
        "r2": float(1.0 - np.dot(y_values - fitted, y_values - fitted) / total),
        "coefficients": {
            key: float(coefficients[index + 1]) for index, key in enumerate(keys)
        },
        "normalization": "Each scene is standardized before equal-capped pooling; coefficients are associations, not causal terrain corrections.",
    }


def _masked(values: np.ndarray, valid: np.ndarray) -> np.ndarray:
    return np.where(valid, np.asarray(values, dtype=np.float32), np.nan)


def _limits(values: list[np.ndarray], *, fallback: float) -> tuple[float, float]:
    finite = np.concatenate([value[np.isfinite(value)] for value in values])
    if finite.size == 0:
        return -fallback, fallback
    limit = max(float(np.percentile(np.abs(finite), 99.0)), fallback)
    return -limit, limit


def _record(
    ee: Any,
    *,
    matchup_id: str,
    scene_id: str,
    setting: str,
    rows: dict[str, dict[str, str]],
    elevations: dict[str, float],
) -> tuple[dict[str, Any], dict[str, np.ndarray]]:
    from bestpixel._gee import utm_epsg_from_bbox, utm_grid
    from bestpixel.atmosphere import AtmoSidecar

    saved_scene = _scene_from_archive(matchup_id, scene_id)
    row = rows[matchup_id]
    lon = float(row["longitude"])
    lat = float(row["latitude"])
    bbox = (
        lon - HALF_SIZE_DEGREES,
        lat - HALF_SIZE_DEGREES,
        lon + HALF_SIZE_DEGREES,
        lat + HALF_SIZE_DEGREES,
    )
    terrain_grid = utm_grid(bbox, utm_epsg_from_bbox(bbox), 60.0)
    grids, metadata = _load_pair_grids(
        ee,
        bbox=bbox,
        sidecar_path=Path(saved_scene["sidecar"]),
        scene=_pair_scene_input(saved_scene),
        elevation_km=_site_elevation_km(row, elevations),
    )
    valid = np.asarray(grids["valid"], dtype=bool)
    l2a_surface = np.asarray(grids["l2a_surface"], dtype=np.float32)
    current_target = np.asarray(grids["siac_surface"], dtype=np.float32)
    l1_toa = np.asarray(grids["l1_toa"], dtype=np.float32)
    l2a_wvp = np.asarray(grids["l2a_tcwv"], dtype=np.float32)

    sidecar = AtmoSidecar.load(str(saved_scene["sidecar"]))
    atmosphere = sidecar.scenes[scene_id]
    l2a_wvp_target = sidecar.correct(
        atmosphere,
        CANONICAL_BANDS,
        l1_toa,
        wvp=l2a_wvp,
    )
    terrain = fetch_glo30_terrain(ee, terrain_grid)
    elevation = terrain.elevation_m
    slope = terrain.slope_deg
    incidence = local_solar_incidence(
        terrain,
        sza_deg=metadata["sza_deg"],
        saa_deg=metadata["saa_deg"],
    )
    l2a_visible = np.mean(l2a_surface[list(VISIBLE_BANDS)], axis=0)
    target_visible = np.mean(current_target[list(VISIBLE_BANDS)], axis=0)
    wvp_target_visible = np.mean(l2a_wvp_target[list(VISIBLE_BANDS)], axis=0)
    residual = l2a_visible - target_visible
    wvp_effect = wvp_target_visible - target_visible
    adjusted_residual = l2a_visible - wvp_target_visible
    wvp_departure = l2a_wvp - np.float32(atmosphere.wvp)
    valid &= terrain.valid & np.isfinite(l2a_wvp_target).all(axis=0)

    fields = {
        "elevation": elevation,
        "slope": slope,
        "incidence": incidence,
        "wvp_departure": wvp_departure,
        "target_visible": target_visible,
    }
    scene_model = _scene_model(residual, fields, valid)
    record = {
        "matchup_id": matchup_id,
        "site": row["site"],
        "setting": setting,
        "scene_id": scene_id,
        "day": metadata["day"],
        "maiac_aot": metadata["maiac_aot"],
        "current_target_wvp_cm": float(atmosphere.wvp),
        "l2a_wvp_median_cm": float(np.nanmedian(l2a_wvp[valid])),
        "l2a_wvp_minus_target_median_cm": float(np.nanmedian(wvp_departure[valid])),
        "valid_pixels": int(np.count_nonzero(valid)),
        "valid_fraction": float(np.mean(valid)),
        "visible_current_target_mae": float(np.mean(np.abs(residual[valid]))),
        "visible_l2a_wvp_target_mae": float(np.mean(np.abs(adjusted_residual[valid]))),
        "visible_wvp_target_effect_mae": float(np.mean(np.abs(wvp_effect[valid]))),
        "visible_wvp_target_effect_median": float(np.median(wvp_effect[valid])),
        "correlations": {
            "elevation": _correlation(residual, elevation, valid),
            "slope": _correlation(residual, slope, valid),
            "solar_incidence": _correlation(residual, incidence, valid),
            "l2a_wvp_minus_target": _correlation(residual, wvp_departure, valid),
        },
        "standardized_scene_model": scene_model,
        "grid": metadata["grid"],
        "terrain_source": "COPERNICUS/DEM/GLO30_2024_1, resampled to the exact 60 m pair grid",
    }
    arrays = {
        "elevation": _masked(elevation, valid),
        "slope": _masked(slope, valid),
        "incidence": _masked(incidence, valid),
        "wvp_departure": _masked(wvp_departure, valid),
        "residual": _masked(residual, valid),
        "wvp_effect": _masked(wvp_effect, valid),
        "model_residual": residual,
        "model_valid": valid,
        "model_elevation": elevation,
        "model_slope": slope,
        "model_incidence": incidence,
        "model_wvp_departure": wvp_departure,
        "model_target_visible": target_visible,
    }
    return record, arrays


def _render(
    path: Path,
    examples: list[tuple[dict[str, Any], dict[str, np.ndarray]]],
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    cmaps = {name: plt.get_cmap(name).copy() for name in ("terrain", "magma", "cividis", "BrBG", "RdBu_r")}
    for cmap in cmaps.values():
        cmap.set_bad("#f3f5f5")
    elevation_values = [arrays["elevation"] for _record, arrays in examples]
    elevation_finite = np.concatenate([value[np.isfinite(value)] for value in elevation_values])
    elevation_limits = tuple(float(value) for value in np.percentile(elevation_finite, (1.0, 99.0)))
    wvp_limits = _limits([arrays["wvp_departure"] for _record, arrays in examples], fallback=0.25)
    effect_limits = _limits([arrays["wvp_effect"] for _record, arrays in examples], fallback=0.001)
    figure, axes = plt.subplots(len(examples), 6, figsize=(23.0, 18.2), facecolor="white")
    columns = (
        ("elevation", "DEM elevation (m)", cmaps["terrain"], elevation_limits),
        ("slope", "DEM slope (degrees)", cmaps["magma"], (0.0, 35.0)),
        ("incidence", "Local solar incidence cos(i)", cmaps["cividis"], (-0.25, 1.0)),
        ("wvp_departure", "L2A WVP - scalar L1C WVP (cm)", cmaps["BrBG"], wvp_limits),
        ("residual", "Visible residual: L2A - current-RT", cmaps["RdBu_r"], (-0.04, 0.04)),
        ("wvp_effect", "Visible target shift from L2A WVP", cmaps["RdBu_r"], effect_limits),
    )
    for row_axes, (record, arrays) in zip(axes, examples, strict=True):
        for axis, (key, _title, cmap, limits) in zip(row_axes, columns, strict=True):
            axis.imshow(arrays[key], cmap=cmap, vmin=limits[0], vmax=limits[1], interpolation="nearest")
            axis.set_xticks([])
            axis.set_yticks([])
        correlations = record["correlations"]
        row_axes[0].set_ylabel(
            (
                f"{record['site']}\n{record['day']} | AOD {record['maiac_aot']:.3f}\n"
                f"r elev {correlations['elevation']:+.2f}; "
                f"slope {correlations['slope']:+.2f}; "
                f"incidence {correlations['solar_incidence']:+.2f}"
            ),
            fontsize=8.5,
            rotation=0,
            ha="right",
            va="center",
        )
    for axis, (_key, title, _cmap, _limits_value) in zip(
        axes[0], columns, strict=True
    ):
        axis.set_title(title, fontsize=9.5, pad=10)
    figure.suptitle(
        "Terrain and water-vapour attribution diagnostic for selected exact pairs",
        x=0.03,
        y=0.993,
        ha="left",
        fontsize=15,
    )
    figure.text(
        0.03,
        0.973,
        "Correlations are within-scene pixel associations with the visible L2A-current-RT residual; they do not establish causality. The final column changes only the L1C correction WVP to the same-pixel L2A WVP.",
        ha="left",
        fontsize=8.8,
    )
    figure.subplots_adjust(left=0.16, right=0.99, top=0.945, bottom=0.035, hspace=0.16, wspace=0.035)
    figure.savefig(path, dpi=155, bbox_inches="tight")
    plt.close(figure)


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    import csv as csv_module

    fields = (
        "site",
        "scene_id",
        "day",
        "setting",
        "maiac_aot",
        "current_target_wvp_cm",
        "l2a_wvp_median_cm",
        "l2a_wvp_minus_target_median_cm",
        "valid_pixels",
        "valid_fraction",
        "visible_current_target_mae",
        "visible_l2a_wvp_target_mae",
        "visible_wvp_target_effect_mae",
        "visible_wvp_target_effect_median",
        "correlation_elevation",
        "correlation_slope",
        "correlation_solar_incidence",
        "correlation_l2a_wvp_minus_target",
        "model_r2",
    )
    flattened = []
    for row in rows:
        flattened.append(
            {
                "site": row["site"],
                "scene_id": row["scene_id"],
                "day": row["day"],
                "setting": row["setting"],
                "maiac_aot": row["maiac_aot"],
                "current_target_wvp_cm": row["current_target_wvp_cm"],
                "l2a_wvp_median_cm": row["l2a_wvp_median_cm"],
                "l2a_wvp_minus_target_median_cm": row["l2a_wvp_minus_target_median_cm"],
                "valid_pixels": row["valid_pixels"],
                "valid_fraction": row["valid_fraction"],
                "visible_current_target_mae": row["visible_current_target_mae"],
                "visible_l2a_wvp_target_mae": row["visible_l2a_wvp_target_mae"],
                "visible_wvp_target_effect_mae": row["visible_wvp_target_effect_mae"],
                "visible_wvp_target_effect_median": row["visible_wvp_target_effect_median"],
                "correlation_elevation": row["correlations"]["elevation"],
                "correlation_slope": row["correlations"]["slope"],
                "correlation_solar_incidence": row["correlations"]["solar_incidence"],
                "correlation_l2a_wvp_minus_target": row["correlations"]["l2a_wvp_minus_target"],
                "model_r2": row["standardized_scene_model"]["r2"],
            }
        )
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv_module.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(flattened)


def build(output: Path) -> dict[str, Any]:
    from bestpixel._gee import init_ee

    output.mkdir(parents=True, exist_ok=True)
    rows = {row["matchup_id"]: row for row in csv.DictReader(MATCHUPS.open())}
    elevations = json.loads(ELEVATION_CACHE.read_text(encoding="utf-8"))
    ee = init_ee()
    examples: list[tuple[dict[str, Any], dict[str, np.ndarray]]] = []
    model_inputs: list[tuple[np.ndarray, dict[str, np.ndarray], np.ndarray]] = []
    for matchup_id, scene_id, setting in EXAMPLES:
        record, arrays = _record(
            ee,
            matchup_id=matchup_id,
            scene_id=scene_id,
            setting=setting,
            rows=rows,
            elevations=elevations,
        )
        examples.append((record, arrays))
        model_inputs.append(
            (
                arrays["model_residual"],
                {
                    "elevation": arrays["model_elevation"],
                    "slope": arrays["model_slope"],
                    "incidence": arrays["model_incidence"],
                    "wvp_departure": arrays["model_wvp_departure"],
                    "target_visible": arrays["model_target_visible"],
                },
                arrays["model_valid"],
            )
        )
        print(
            f"TERRAIN_WVP {record['site']} n={record['valid_pixels']} "
            f"r_slope={record['correlations']['slope']:+.3f} "
            f"wvp_effect={record['visible_wvp_target_effect_mae']:.5f}",
            flush=True,
        )
    _render(output / "terrain-wvp-diagnostic.png", examples)
    records = [record for record, _arrays in examples]
    _write_csv(output / "terrain-wvp-diagnostic.csv", records)
    manifest = {
        "title": "Terrain and water-vapour attribution diagnostic",
        "scope": "Four selected exact same-day pairs used for image-level inspection; this is not a 36-case AOD score.",
        "current_target": "Saved L1C/current-RT correction uses one scene water-vapour value and one site elevation in its sidecar coefficients across each 60 m pair grid.",
        "l2a_product": "Operational COPERNICUS/S2_SR_HARMONIZED L2A includes the per-pixel WVP map used here and operational Sen2Cor terrain correction.",
        "wvp_counterfactual": "The final panel re-corrects the same L1C TOA with unchanged AOD, geometry and coefficient LUT but substitutes L2A WVP per pixel. It is a sensitivity calculation, not a production target because it would reuse an L2A retrieval inside its own harmonization target.",
        "terrain_method": "Copernicus GLO-30 2024.1 DEM is sampled to the common 60 m UTM grid. Slope and local cosine incidence are calculated from DEM gradients and saved Sentinel sun zenith/azimuth.",
        "correlation_method": "Pearson correlations are computed within each scene on the exact clear-land mask. The standardized models also include target visible reflectance as a brightness control; they are observational associations and do not prove terrain causality.",
        "figure": "assets/pair-examples/terrain-wvp-diagnostic.png",
        "csv": "assets/pair-examples/terrain-wvp-diagnostic.csv",
        "pooled_standardized_model": _pooled_model(model_inputs),
        "examples": records,
    }
    (output / "terrain-wvp-diagnostic.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    print(json.dumps(build(args.output), indent=2))


if __name__ == "__main__":
    main()
