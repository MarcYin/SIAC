"""Train a scene-level L2A-to-current-RT harmonization control.

Unlike the pixel Ridge mapper, this model predicts one per-band correction for
an acquisition.  Its inputs are robust L2A spectral summaries plus the
atmospheric, geometry, elevation, platform, and processing metadata available
at composite-build time.  Site-group cross-fitting prevents a target site's
history from training its own correction, and AERONET is never read.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any

import joblib
import numpy as np
from sklearn.ensemble import ExtraTreesRegressor
from tools.aeronet_validation.train_l2a_l1c_harmonizer import BAND_NAMES

ROOT = Path("/gws/ssde/j25a/nceo_isp/public/siac_refactor")
DEFAULT_PAIRS = ROOT / "analysis/l2a_l1c_exact_pairs_mediumdev_20260713"
DEFAULT_MODEL = ROOT / "analysis/l2a_l1c_harmonizer_mediumdev_20260713/harmonizer.json"
DEFAULT_CASES = ROOT / "reports/aod-medium-physics-20260713/downloads/cases.csv"
DEFAULT_OUTPUT = ROOT / "analysis/l2a_l1c_harmonizer_mediumdev_20260713/scene_model"


def _finite(value: Any, default: float = 0.0) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError):
        return float(default)
    return result if math.isfinite(result) else float(default)


def feature_names() -> list[str]:
    names = [
        f"l2a_{quantile}_{band}" for quantile in ("p10", "median", "p90") for band in BAND_NAMES
    ]
    names.extend(
        [
            "l2a_aot_p10",
            "l2a_aot_median",
            "l2a_aot_p90",
            "l2a_tcwv_p10",
            "l2a_tcwv_median",
            "l2a_tcwv_p90",
            "maiac_aot",
            "maiac_tcwv",
            "maiac_minus_l2a_aot",
            "maiac_minus_l2a_tcwv",
            "sza_deg",
            "vza_deg",
            "cos_raa",
            "airmass",
            "elevation_km",
            "month_sin",
            "month_cos",
            "sensor_is_s2b",
            "sensor_is_s2c",
            "processing_baseline",
            "ndvi_proxy",
            "swir_ratio",
            "visible_slope",
        ]
    )
    return names


def build_scene_features(
    l2a: np.ndarray,
    l2a_aot: np.ndarray,
    l2a_tcwv: np.ndarray,
    scene: dict[str, Any],
) -> np.ndarray:
    surface = np.asarray(l2a, dtype=np.float64)
    if surface.ndim != 2 or surface.shape[1] != len(BAND_NAMES):
        raise ValueError(f"l2a must have shape (sample, {len(BAND_NAMES)})")
    quantiles = np.nanpercentile(surface, (10.0, 50.0, 90.0), axis=0).reshape(-1)
    aot_q = np.nanpercentile(np.asarray(l2a_aot, dtype=np.float64), (10.0, 50.0, 90.0))
    tcwv_q = np.nanpercentile(np.asarray(l2a_tcwv, dtype=np.float64), (10.0, 50.0, 90.0))
    median = quantiles[len(BAND_NAMES) : 2 * len(BAND_NAMES)]
    maiac_aot = _finite(scene.get("maiac_aot"))
    maiac_tcwv = _finite(scene.get("maiac_tcwv_cm"))
    sza = _finite(scene.get("sza_deg"))
    vza = _finite(scene.get("vza_deg"))
    raa = _finite(scene.get("raa_deg"))
    month = max(1.0, min(12.0, _finite(scene.get("month"), 1.0)))
    l2a_meta = scene.get("l2a") or {}
    ancillary = np.asarray(
        [
            *aot_q,
            *tcwv_q,
            maiac_aot,
            maiac_tcwv,
            maiac_aot - float(aot_q[1]),
            maiac_tcwv - float(tcwv_q[1]),
            sza,
            vza,
            math.cos(math.radians(raa)),
            1.0 / max(math.cos(math.radians(sza)), 0.1)
            + 1.0 / max(math.cos(math.radians(vza)), 0.1),
            _finite(scene.get("elevation_km")),
            math.sin(2.0 * math.pi * month / 12.0),
            math.cos(2.0 * math.pi * month / 12.0),
            float("2B" in str(l2a_meta.get("spacecraft", "")).upper()),
            float("2C" in str(l2a_meta.get("spacecraft", "")).upper()),
            _finite(l2a_meta.get("processing_baseline")),
            (median[4] - median[3]) / max(median[4] + median[3], 1.0e-4),
            median[5] / max(median[6], 1.0e-4),
            median[3] - median[1],
        ],
        dtype=np.float64,
    )
    features = np.concatenate([quantiles, ancillary])
    if features.size != len(feature_names()) or not np.all(np.isfinite(features)):
        raise ValueError("scene feature construction produced an invalid vector")
    return features


def load_scene_table(
    pair_dir: Path,
    matchup_ids: list[str],
    *,
    cutoff: str,
) -> tuple[np.ndarray, np.ndarray, list[dict[str, Any]]]:
    features: list[np.ndarray] = []
    targets: list[np.ndarray] = []
    rows: list[dict[str, Any]] = []
    for matchup_id in matchup_ids:
        with np.load(pair_dir / f"{matchup_id}.npz", allow_pickle=False) as data:
            l2a = np.asarray(data["l2a"], dtype=np.float64)
            siac = np.asarray(data["siac"], dtype=np.float64)
            l2a_aot = np.asarray(data["l2a_aot"], dtype=np.float64)
            l2a_tcwv = np.asarray(data["l2a_tcwv"], dtype=np.float64)
            scene_index = np.asarray(data["scene_index"], dtype=np.int32)
            scenes = json.loads(str(data["scenes_json"].item()))
        for local_index, scene in enumerate(scenes):
            if str(scene.get("day", "")) > cutoff:
                continue
            selected = scene_index == local_index
            if int(np.count_nonzero(selected)) < 100:
                continue
            features.append(
                build_scene_features(
                    l2a[selected],
                    l2a_aot[selected],
                    l2a_tcwv[selected],
                    scene,
                )
            )
            residual = np.nanmedian(siac[selected] - l2a[selected], axis=0)
            targets.append(np.asarray(residual, dtype=np.float64))
            rows.append(
                {
                    "matchup_id": matchup_id,
                    "site": matchup_id.split("__")[0],
                    "scene_id": scene.get("scene_id"),
                    "day": scene.get("day"),
                    "window": scene.get("window"),
                }
            )
    return np.stack(features), np.stack(targets), rows


def _model(*, leaf: int, max_features: float, trees: int) -> ExtraTreesRegressor:
    return ExtraTreesRegressor(
        n_estimators=trees,
        min_samples_leaf=leaf,
        max_features=max_features,
        random_state=0,
        n_jobs=-1,
    )


def _metrics(target: np.ndarray, prediction: np.ndarray, cap: float) -> dict[str, Any]:
    error = target - np.clip(prediction, -cap, cap)
    per_band = {
        band: {
            "bias": float(np.mean(error[:, index])),
            "mae": float(np.mean(np.abs(error[:, index]))),
            "rmse": float(np.sqrt(np.mean(np.square(error[:, index])))),
        }
        for index, band in enumerate(BAND_NAMES)
    }
    return {
        "cap": cap,
        "per_band": per_band,
        "visible_mean_mae": float(
            np.mean([per_band[BAND_NAMES[index]]["mae"] for index in (1, 2, 3)])
        ),
        "all_mean_mae": float(np.mean([value["mae"] for value in per_band.values()])),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pairs", type=Path, default=DEFAULT_PAIRS)
    parser.add_argument("--base-model", type=Path, default=DEFAULT_MODEL)
    parser.add_argument("--cases", type=Path, default=DEFAULT_CASES)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--training-cutoff", default="2023-12-31")
    parser.add_argument("--trees", type=int, default=100)
    args = parser.parse_args()

    matchup_ids = [row["matchup_id"] for row in csv.DictReader(args.cases.open())]
    base = json.loads(args.base_model.read_text(encoding="utf-8"))
    folds = {key: int(value) for key, value in base["fold_by_matchup_id"].items()}
    x, y, rows = load_scene_table(args.pairs, matchup_ids, cutoff=args.training_cutoff)
    scene_folds = np.asarray([folds[row["matchup_id"]] for row in rows], dtype=np.int16)

    specs = [(leaf, max_features) for leaf in (3, 8, 20) for max_features in (0.7, 1.0)]
    candidates: dict[str, dict[str, Any]] = {}
    fitted_folds: dict[str, dict[str, ExtraTreesRegressor]] = {}
    best_key = ""
    best_score = float("inf")
    for leaf, max_features in specs:
        key = f"et{args.trees}_leaf{leaf}_mf{max_features:.1f}".replace(".", "p")
        oof = np.full_like(y, np.nan)
        fold_models: dict[str, ExtraTreesRegressor] = {}
        for fold in range(5):
            train = scene_folds != fold
            test = scene_folds == fold
            model = _model(leaf=leaf, max_features=max_features, trees=args.trees).fit(
                x[train], y[train]
            )
            oof[test] = model.predict(x[test])
            fold_models[str(fold)] = model
        metrics = {str(cap): _metrics(y, oof, cap) for cap in (0.015, 0.03)}
        candidates[key] = {"leaf": leaf, "max_features": max_features, "metrics": metrics}
        fitted_folds[key] = fold_models
        score = metrics["0.03"]["visible_mean_mae"]
        if score < best_score:
            best_score = score
            best_key = key

    best_spec = candidates[best_key]
    full_model = _model(
        leaf=int(best_spec["leaf"]),
        max_features=float(best_spec["max_features"]),
        trees=args.trees,
    ).fit(x, y)
    artifact = {
        "schema_version": 1,
        "uses_aeronet": False,
        "target": "scene median current-RT L1C BOA minus operational L2A BOA",
        "training_cutoff": args.training_cutoff,
        "feature_names": feature_names(),
        "band_names": list(BAND_NAMES),
        "fold_by_matchup_id": folds,
        "best_model": best_key,
        "candidates": candidates,
        "fold_models": fitted_folds[best_key],
        "full_model": full_model,
        "training_scene_count": len(rows),
        "training_site_count": len({row["site"] for row in rows}),
    }
    args.output.mkdir(parents=True, exist_ok=True)
    joblib.dump(artifact, args.output / "scene_harmonizer.joblib", compress=3)
    summary = {key: value for key, value in artifact.items() if "model" not in key}
    summary["best_model"] = best_key
    summary["candidates"] = candidates
    (args.output / "metrics.json").write_text(json.dumps(summary, indent=2) + "\n")
    with (args.output / "scenes.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(json.dumps({"best_model": best_key, "best_visible_mae": best_score, "scenes": len(rows)}))


if __name__ == "__main__":
    main()
