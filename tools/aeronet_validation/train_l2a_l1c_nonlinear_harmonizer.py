"""Cross-fit a nonlinear L2A-to-current-RT harmonizer on the raw atmospheric state.

The ridge cross-RT family models the correction as a linear function of
hand-crafted ΔAOT interactions; the observed response is visibly nonlinear.
This trainer learns a per-band gradient-boosted mapping from the full raw
state — all seven L2A bands, MAIAC AOD and Sen2Cor AOD as separate inputs,
water vapour, viewing/solar geometry, month, sensor, and scene plus per-pixel
elevation/terrain — to the same target: the same-day L1C surface corrected
with MAIAC AOD under the current libRadtran LUT. Pair loading, the locked
site-grouped fold protocol, acquisition-balanced weighting, and every metric
definition are imported from the ridge trainer so results are directly
comparable. AERONET is never read.
"""

from __future__ import annotations

import argparse
import csv
import json
import time
from pathlib import Path
from typing import Any

import joblib
import numpy as np
from sklearn.ensemble import HistGradientBoostingRegressor
from tools.aeronet_validation.train_l2a_l1c_harmonizer import (
    BAND_NAMES,
    CURRENT_OFFSETS,
    PairDataset,
    _add_scene_errors,
    _case_ids,
    _distribution,
    _fixed_fold_splits,
    _metrics,
    _scene_balanced_sample_weight,
    _scene_rows,
    load_pairs,
)

ROOT = Path("/gws/ssde/j25a/nceo_isp/public/siac_refactor")
DEFAULT_PAIRS = ROOT / "analysis/l2a_l1c_physical_pairs_lowcloud152_20260716"
DEFAULT_CASES = ROOT / "reports/aod-final-performance-dashboard-20260713/data/all-cases.csv"
DEFAULT_SPLIT = ROOT / "analysis/l2a_harmonization_retrieval_score_20260714/summary.json"
DEFAULT_OUTPUT = ROOT / "analysis/cross_rt_nonlinear_harmonizer_lowcloud152_20260716"
DEFAULT_RIDGE_METRICS = (
    ROOT / "analysis/cross_rt_harmonizer_lowcloud152_20260716/surface_metrics.json"
)

MODEL_NAME = "hgb_full_state"
# Sen2Cor selects the LUT set whose ozone column is NEAREST these nodes (ATBD
# 4.2/4.6, no interpolation); the quantized value mimics the L2A frame's ozone.
SEN2COR_OZONE_NODES = (250.0, 290.0, 331.0, 370.0, 410.0, 450.0)


def sen2cor_ozone_node(values):
    nodes = np.asarray(SEN2COR_OZONE_NODES, dtype=np.float64)
    values = np.asarray(values, dtype=np.float64)
    return nodes[np.argmin(np.abs(values[:, None] - nodes[None, :]), axis=1)]


HGB_PARAMS: dict[str, Any] = {
    "loss": "squared_error",
    "learning_rate": 0.06,
    "max_iter": 400,
    "max_leaf_nodes": 63,
    "min_samples_leaf": 200,
    "l2_regularization": 1.0,
    "max_bins": 255,
    "early_stopping": False,
    "random_state": 0,
}


STATE_FEATURE_NAMES = (
    "maiac_aot",
    "sen2cor_aot",
    "l2a_tcwv_cm",
    "maiac_tcwv_cm",
    "sza_deg",
    "vza_deg",
    "raa_deg",
    "cos_raa",
    "month",
    "elevation_km",
    "terrain_elevation_km",
    "terrain_slope_deg",
    "terrain_incidence_cos",
    "sensor_is_s2b",
    "sensor_is_s2c",
    "processing_baseline",
)


def feature_names(band_input: str = "all", band_name: str | None = None) -> list[str]:
    """Column order for one band model: reflectance inputs first, then the raw state."""
    if band_input == "all":
        return [*[f"l2a_{band}" for band in BAND_NAMES], *STATE_FEATURE_NAMES]
    if band_input == "target":
        if band_name is None:
            raise ValueError("the target-band schema needs the band name")
        return [f"l2a_{band_name}", *STATE_FEATURE_NAMES]
    raise ValueError(f"unsupported band input {band_input!r}")


def build_columns(dataset: PairDataset) -> dict[str, np.ndarray]:
    """Every named per-pixel column; band models assemble their declared subset."""
    if (
        dataset.terrain_elevation_km is None
        or dataset.terrain_slope_deg is None
        or dataset.terrain_incidence_cos is None
    ):
        raise SystemExit("the nonlinear harmonizer requires terrain-enabled pair archives")
    return {
        **{
            f"l2a_{band}": dataset.l2a[:, index].astype(np.float64)
            for index, band in enumerate(BAND_NAMES)
        },
        "maiac_aot": np.asarray(dataset.maiac_aot, dtype=np.float64),
        "sen2cor_aot": np.asarray(dataset.l2a_aot, dtype=np.float64),
        "l2a_tcwv_cm": np.asarray(dataset.l2a_tcwv, dtype=np.float64),
        "maiac_tcwv_cm": np.asarray(dataset.maiac_tcwv, dtype=np.float64),
        "sza_deg": np.asarray(dataset.sza_deg, dtype=np.float64),
        "vza_deg": np.asarray(dataset.vza_deg, dtype=np.float64),
        "raa_deg": np.asarray(dataset.raa_deg, dtype=np.float64),
        "cos_raa": np.cos(np.radians(np.asarray(dataset.raa_deg, dtype=np.float64))),
        "month": np.asarray(dataset.month, dtype=np.float64),
        "elevation_km": np.asarray(dataset.elevation_km, dtype=np.float64),
        "terrain_elevation_km": np.asarray(dataset.terrain_elevation_km, dtype=np.float64),
        "terrain_slope_deg": np.clip(
            np.asarray(dataset.terrain_slope_deg, dtype=np.float64), 0.0, 70.0
        ),
        "terrain_incidence_cos": np.clip(
            np.asarray(dataset.terrain_incidence_cos, dtype=np.float64), -1.0, 1.0
        ),
        "sensor_is_s2b": np.asarray(dataset.sensor_is_s2b, dtype=np.float64),
        "sensor_is_s2c": np.asarray(dataset.sensor_is_s2c, dtype=np.float64),
        "processing_baseline": np.asarray(dataset.processing_baseline, dtype=np.float64),
    }


def van_heuklon_ozone_du(lat_deg: float, lon_deg: float, day_of_year: int) -> float:
    """Van Heuklon (1979) total-ozone climatology in Dobson units.

    A deterministic proxy for the ECMWF ozone column Sen2Cor selects per scene
    (and quantizes to its LUT nodes); it carries the latitudinal/seasonal cycle.
    """
    radians = np.pi / 180.0
    if lat_deg >= 0.0:
        seasonal = 150.0 + 40.0 * np.sin(0.9865 * (day_of_year - 30.0) * radians)
        zonal = 20.0 * np.sin(3.0 * (lon_deg + 20.0) * radians)
        shape = np.sin(1.28 * lat_deg * radians) ** 2
    else:
        seasonal = 100.0 + 30.0 * np.sin(0.9865 * (day_of_year + 152.625) * radians)
        zonal = 20.0 * np.sin(2.0 * (lon_deg - 75.0) * radians)
        shape = np.sin(1.5 * lat_deg * radians) ** 2
    return float(235.0 + (seasonal + zonal) * shape)


def ozone_climatology_column(dataset: PairDataset, matchups_path: Path) -> np.ndarray:
    """Per-sample ozone column from the site coordinates and each scene's day."""
    import datetime

    coordinates = {
        row["matchup_id"]: (float(row["latitude"]), float(row["longitude"]))
        for row in csv.DictReader(matchups_path.open())
    }
    per_scene = np.empty(len(dataset.scene_metadata), dtype=np.float64)
    for index, scene in enumerate(dataset.scene_metadata):
        lat, lon = coordinates[str(scene["matchup_id"])]
        day_of_year = datetime.date.fromisoformat(str(scene["day"])).timetuple().tm_yday
        per_scene[index] = van_heuklon_ozone_du(lat, lon, day_of_year)
    return per_scene[np.asarray(dataset.scene_code, dtype=np.int64)]


def cams_ozone_column(dataset: PairDataset, table_path: Path) -> np.ndarray:
    """Per-sample real CAMS total-ozone column keyed by (matchup, scene day)."""
    lookup = {
        (row["matchup_id"], row["day"]): float(row["cams_tco3_du"])
        for row in csv.DictReader(table_path.open())
    }
    per_scene = np.empty(len(dataset.scene_metadata), dtype=np.float64)
    for index, scene in enumerate(dataset.scene_metadata):
        key = (str(scene["matchup_id"]), str(scene["day"]))
        if key not in lookup:
            raise SystemExit(f"CAMS scene-state table is missing {key}")
        per_scene[index] = lookup[key]
    return per_scene[np.asarray(dataset.scene_code, dtype=np.int64)]


def load_locked_split(path: Path, matchup_ids: list[str]) -> tuple[dict[str, int], set[int]]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    cohort = payload["cohort"]
    fold_by_matchup_id = {key: int(value) for key, value in cohort["fold_by_matchup_id"].items()}
    holdout_folds = {int(value) for value in cohort["holdout_folds"]}
    missing = sorted(set(matchup_ids) - set(fold_by_matchup_id))
    if missing:
        raise SystemExit(f"locked split is missing {len(missing)} matchup(s)")
    return fold_by_matchup_id, holdout_folds


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pairs", type=Path, default=DEFAULT_PAIRS)
    parser.add_argument("--cases", type=Path, default=DEFAULT_CASES)
    parser.add_argument("--evaluation-split-manifest", type=Path, default=DEFAULT_SPLIT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--training-cutoff", default="2023-12-31")
    parser.add_argument("--max-samples-per-scene", type=int, default=300)
    parser.add_argument("--ridge-metrics", type=Path, default=DEFAULT_RIDGE_METRICS)
    parser.add_argument(
        "--band-input",
        choices=("all", "target"),
        default="all",
        help="feed every L2A band to each model, or only the corrected band itself",
    )
    parser.add_argument("--model-name", default=None)
    parser.add_argument("--learning-rate", type=float, default=HGB_PARAMS["learning_rate"])
    parser.add_argument("--max-iter", type=int, default=HGB_PARAMS["max_iter"])
    parser.add_argument("--max-leaf-nodes", type=int, default=HGB_PARAMS["max_leaf_nodes"])
    parser.add_argument("--min-samples-leaf", type=int, default=HGB_PARAMS["min_samples_leaf"])
    parser.add_argument("--l2-regularization", type=float, default=HGB_PARAMS["l2_regularization"])
    parser.add_argument(
        "--ozone-climatology",
        action="store_true",
        help="add a per-scene Van Heuklon total-ozone column (Sen2Cor ozone proxy)",
    )
    parser.add_argument(
        "--cams-state",
        type=Path,
        default=None,
        help="scene-state table from extract_cams_scene_state.py; adds real CAMS ozone",
    )
    parser.add_argument(
        "--training-maiac-aot-max",
        type=float,
        default=None,
        help="fit only on scenes at/below this MAIAC AOD (cleaner targets); "
        "out-of-fold evaluation still covers every scene",
    )
    parser.add_argument("--matchups", type=Path, default=ROOT / "matchups/matchups.csv")
    args = parser.parse_args()
    model_name = args.model_name or (MODEL_NAME if args.band_input == "all" else "hgb_target_band")
    hgb_params = {
        **HGB_PARAMS,
        "learning_rate": args.learning_rate,
        "max_iter": args.max_iter,
        "max_leaf_nodes": args.max_leaf_nodes,
        "min_samples_leaf": args.min_samples_leaf,
        "l2_regularization": args.l2_regularization,
    }

    matchup_ids = _case_ids(args.cases)
    fold_by_matchup_id, holdout_folds = load_locked_split(
        args.evaluation_split_manifest, matchup_ids
    )
    dataset = load_pairs(
        args.pairs,
        matchup_ids,
        scene_day_max=args.training_cutoff,
        max_samples_per_scene=args.max_samples_per_scene,
        allow_missing_matchups=True,
    )

    raw_folds = np.full(dataset.site_code.shape, -1, dtype=np.int16)
    fold_by_site: dict[int, int] = {}
    for site_code, site_matchups in dataset.matchup_by_site_code.items():
        assigned = {fold_by_matchup_id[matchup_id] for matchup_id in site_matchups}
        if len(assigned) != 1:
            raise SystemExit(f"locked split divides one site across folds: {site_matchups}")
        fold_by_site[int(site_code)] = assigned.pop()
    for site_code, fold in fold_by_site.items():
        raw_folds[dataset.site_code == site_code] = fold
    application_folds, split_rows, split_protocol = _fixed_fold_splits(
        raw_folds, holdout_folds=holdout_folds
    )
    if np.any(application_folds < 0):
        raise RuntimeError("cross-fitting did not assign every sample")
    fold_remap = {
        int(raw): int(applied) for raw, applied in zip(raw_folds, application_folds, strict=True)
    }
    applied_fold_by_matchup_id = {
        matchup_id: fold_remap.get(int(fold), int(fold))
        for matchup_id, fold in fold_by_matchup_id.items()
    }

    columns = build_columns(dataset)
    names_by_band = {band: feature_names(args.band_input, band) for band in BAND_NAMES}
    if args.ozone_climatology:
        columns["ozone_du_climatology"] = ozone_climatology_column(dataset, args.matchups)
        names_by_band = {
            band: [*names, "ozone_du_climatology"] for band, names in names_by_band.items()
        }
    if args.cams_state is not None:
        cams_ozone = cams_ozone_column(dataset, args.cams_state)
        columns["ozone_du_cams"] = cams_ozone
        columns["ozone_du_sen2cor_node"] = sen2cor_ozone_node(cams_ozone)
        names_by_band = {
            band: [*names, "ozone_du_cams", "ozone_du_sen2cor_node"]
            for band, names in names_by_band.items()
        }
    sample_weight = _scene_balanced_sample_weight(dataset.acquisition_code)
    scene_rows = _scene_rows(dataset)

    args.output.mkdir(parents=True, exist_ok=True)
    model_dir = args.output / "models"
    model_dir.mkdir(exist_ok=True)

    metrics: dict[str, Any] = {
        "sample_count": int(dataset.l2a.shape[0]),
        "scene_count": len(dataset.scene_metadata),
        "acquisition_count": int(np.unique(dataset.acquisition_code).size),
        "site_count": len(dataset.matchup_by_site_code),
        "training_cutoff": args.training_cutoff,
        "model_family": "nonlinear",
        "crossfit_protocol": split_protocol,
        "scene_metric_unit": "independent Sentinel-2 acquisition",
        "identity": {},
        "current_aod_offset": {},
        "candidates": {},
    }
    scene_delta_aot = np.asarray(
        [row["delta_aot_maiac_minus_sen2cor"] for row in scene_rows], dtype=np.float64
    )
    metrics["scene_distribution"] = {
        "delta_aot_maiac_minus_sen2cor": _distribution(scene_delta_aot),
        "near_equal_aod_abs_delta_le_0p02": int(
            np.count_nonzero(np.isfinite(scene_delta_aot) & (np.abs(scene_delta_aot) <= 0.02))
        ),
    }

    identity_error = dataset.l2a - dataset.siac
    _add_scene_errors(
        scene_rows, prefix="identity", error=identity_error, scene_code=dataset.scene_code
    )
    current_error = np.empty_like(identity_error, dtype=np.float64)
    for band_index, band_name in enumerate(BAND_NAMES):
        metrics["identity"][band_name] = _metrics(
            identity_error[:, band_index], dataset.acquisition_code
        )
        intercept, slope = CURRENT_OFFSETS.get(band_index, (0.0, 0.0))
        current_error[:, band_index] = (
            dataset.l2a[:, band_index]
            + intercept
            + slope * dataset.maiac_aot
            - dataset.siac[:, band_index]
        )
        metrics["current_aod_offset"][band_name] = _metrics(
            current_error[:, band_index], dataset.acquisition_code
        )

    if args.training_maiac_aot_max is not None:
        maiac = np.asarray(dataset.maiac_aot, dtype=np.float64)
        training_gate = np.isfinite(maiac) & (maiac <= float(args.training_maiac_aot_max))
        print(
            json.dumps(
                {
                    "training_maiac_aot_max": args.training_maiac_aot_max,
                    "gated_fraction": float(training_gate.mean()),
                }
            ),
            flush=True,
        )
    else:
        training_gate = None

    model_files: dict[str, dict[str, str]] = {str(fold): {} for fold, _t, _s in split_rows}
    residual_prediction = np.full_like(dataset.siac, np.nan, dtype=np.float64)
    for band_index, band_name in enumerate(BAND_NAMES):
        features = np.column_stack([columns[name] for name in names_by_band[band_name]])
        target = (dataset.siac[:, band_index] - dataset.l2a[:, band_index]).astype(np.float64)
        for fold, train, test in split_rows:
            if training_gate is not None:
                train = train[training_gate[train]]
            started = time.time()
            model = HistGradientBoostingRegressor(**hgb_params)
            model.fit(features[train], target[train], sample_weight=sample_weight[train])
            residual_prediction[test, band_index] = model.predict(features[test])
            model_path = model_dir / f"hgb_fold{fold}_{band_name}.joblib"
            joblib.dump(model, model_path, compress=3)
            model_files[str(fold)][band_name] = model_path.name
            print(
                json.dumps(
                    {
                        "fold": int(fold),
                        "band": band_name,
                        "train_pixels": int(train.size),
                        "test_pixels": int(test.size),
                        "fit_seconds": round(time.time() - started, 1),
                    }
                ),
                flush=True,
            )
    if np.any(~np.isfinite(residual_prediction)):
        raise RuntimeError("out-of-fold prediction did not cover every sample")

    oof_residual_scale: dict[str, Any] = {}
    candidate_metrics: dict[str, Any] = {}
    for cap in (0.015, 0.03, 0.06, None):
        cap_key = "uncapped" if cap is None else f"cap_{cap:.3f}"
        corrected = residual_prediction if cap is None else np.clip(residual_prediction, -cap, cap)
        corrected_error = dataset.l2a + corrected - dataset.siac
        candidate_metrics[cap_key] = {
            band_name: _metrics(corrected_error[:, band_index], dataset.acquisition_code)
            for band_index, band_name in enumerate(BAND_NAMES)
        }
        _add_scene_errors(
            scene_rows,
            prefix=f"{model_name}_{cap_key}".replace(".", "p"),
            error=corrected_error,
            scene_code=dataset.scene_code,
        )
    metrics["candidates"][model_name] = candidate_metrics
    for band_index, band_name in enumerate(BAND_NAMES):
        error = (
            dataset.l2a[:, band_index]
            + residual_prediction[:, band_index]
            - dataset.siac[:, band_index]
        )
        center = float(np.median(error))
        oof_residual_scale[band_name] = {
            "median": center,
            "mad_to_sigma": float(np.median(np.abs(error - center)) * 1.4826),
            "rmse": float(np.sqrt(np.mean(np.square(error)))),
        }

    artifact = {
        "schema_version": 1,
        "model_type": "hist_gradient_boosting",
        "model_name": model_name,
        "uses_aeronet": False,
        "training_cutoff": args.training_cutoff,
        "target": dataset.target,
        "target_rt": dataset.target_rt,
        "bands": list(BAND_NAMES),
        "band_input": args.band_input,
        "ozone_climatology": bool(args.ozone_climatology),
        "training_maiac_aot_max": args.training_maiac_aot_max,
        "feature_names_by_band": names_by_band,
        "hyperparameters": {key: value for key, value in hgb_params.items() if key != "loss"},
        "regression_weighting_unit": "independent Sentinel-2 acquisition",
        "fold_by_matchup_id": applied_fold_by_matchup_id,
        "source_fold_by_matchup_id": {
            matchup_id: int(fold) for matchup_id, fold in fold_by_matchup_id.items()
        },
        "crossfit_protocol": split_protocol,
        "model_files": model_files,
        "oof_residual_scale": oof_residual_scale,
        "pair_archive": str(args.pairs),
        "max_samples_per_scene": args.max_samples_per_scene,
        "evaluation_split": {
            "manifest": str(args.evaluation_split_manifest),
            "holdout_folds": sorted(holdout_folds),
        },
    }
    (args.output / "harmonizer.json").write_text(
        json.dumps(artifact, indent=2) + "\n", encoding="utf-8"
    )
    metrics["pair_archive"] = str(args.pairs)
    metrics["max_samples_per_scene"] = args.max_samples_per_scene
    (args.output / "surface_metrics.json").write_text(
        json.dumps(metrics, indent=2) + "\n", encoding="utf-8"
    )
    with (args.output / "surface_scene_metrics.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(scene_rows[0]))
        writer.writeheader()
        writer.writerows(scene_rows)

    comparison: dict[str, Any] = {}
    if args.ridge_metrics.exists():
        ridge = json.loads(args.ridge_metrics.read_text(encoding="utf-8"))
        for cap_key in ("cap_0.015", "cap_0.030", "uncapped"):
            rows = {}
            for band_name in BAND_NAMES:
                ours = candidate_metrics[cap_key][band_name]
                theirs = ridge["candidates"]["cross_rt_terrain_a1"][cap_key][band_name]
                rows[band_name] = {
                    "hgb_scene_mae": ours["scene_mae"],
                    "ridge_scene_mae": theirs["scene_mae"],
                    "hgb_scene_bias": ours["scene_bias"],
                    "ridge_scene_bias": theirs["scene_bias"],
                    "scene_mae_change_vs_ridge_pct": 100.0
                    * (ours["scene_mae"] / theirs["scene_mae"] - 1.0),
                }
            comparison[cap_key] = rows
        (args.output / "ridge_comparison.json").write_text(
            json.dumps(comparison, indent=2) + "\n", encoding="utf-8"
        )

    def visible_mae(cap_key: str) -> float:
        return float(
            np.mean(
                [candidate_metrics[cap_key][band]["scene_mae"] for band in ("blue", "green", "red")]
            )
        )

    print(
        json.dumps(
            {
                "output": str(args.output),
                "pixels": int(dataset.l2a.shape[0]),
                "acquisitions": metrics["acquisition_count"],
                "sites": metrics["site_count"],
                "visible_scene_mae": {
                    "identity": float(
                        np.mean(
                            [
                                metrics["identity"][band]["scene_mae"]
                                for band in ("blue", "green", "red")
                            ]
                        )
                    ),
                    "hgb_cap_0.015": visible_mae("cap_0.015"),
                    "hgb_cap_0.030": visible_mae("cap_0.030"),
                    "hgb_uncapped": visible_mae("uncapped"),
                },
                "oof_residual_sigma_blue": oof_residual_scale["blue"]["mad_to_sigma"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
