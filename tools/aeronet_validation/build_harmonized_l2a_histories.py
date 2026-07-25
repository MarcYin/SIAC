"""Build production-valid L2A histories with per-acquisition harmonization.

Each operational L2A acquisition is corrected before tile mosaicking and
temporal compositing. The correction uses an AERONET-free, site-cross-fitted
model trained on exact same-day L2A and MAIAC-conditioned L1C/current-RT pairs.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import warnings
from collections import defaultdict
from pathlib import Path
from typing import Any

import joblib
import numpy as np
from scipy import ndimage
from tools.aeronet_validation.build_l2a_l1c_harmonization_pairs import (
    CANONICAL_BANDS,
    CANONICAL_TO_S2_BAND,
    GEE_L2A_BANDS,
    L2A_AUX_BANDS,
    LAND_SCL,
    _normalise_aux,
    _normalise_l2a_reflectance,
)
from tools.aeronet_validation.terrain_features import (
    TERRAIN_SOURCE,
    TerrainFields,
    fetch_glo30_terrain,
    local_solar_incidence,
)
from tools.aeronet_validation.train_l2a_l1c_harmonizer import (
    BAND_NAMES,
    VISIBLE_INDICES,
    _predict_model,
    build_features,
    feature_names,
    feature_set_uses_terrain,
)

ROOT = Path("/gws/ssde/j25a/nceo_isp/public/siac_refactor")
DEFAULT_PAIRS = ROOT / "analysis/l2a_l1c_exact_pairs_mediumdev_20260713"
DEFAULT_MODEL = ROOT / "analysis/l2a_l1c_harmonizer_mediumdev_20260713/harmonizer.json"
DEFAULT_OUTPUT = ROOT / "analysis/l2a_l1c_harmonizer_mediumdev_20260713/daily_histories"
DEFAULT_CASES = ROOT / "reports/aod-medium-physics-20260713/downloads/cases.csv"
FALLBACK_HISTORY = ROOT / "l2a_seasonal_pc"
MATCHUPS = ROOT / "matchups/matchups.csv"
HALF_SIZE_DEGREES = 0.06


def _candidate(value: str) -> tuple[str, str, float]:
    parts = value.split(":")
    if len(parts) != 3:
        raise argparse.ArgumentTypeError("candidate must be MODEL:MODE:CAP")
    model_name, mode, cap_text = parts
    if mode not in {"all", "blue", "deepblue", "visible", "solver"}:
        raise argparse.ArgumentTypeError("MODE must be all, blue, deepblue, visible, or solver")
    try:
        cap = float(cap_text)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("CAP must be numeric") from exc
    if not 0.0 < cap <= 0.2:
        raise argparse.ArgumentTypeError("CAP must be in (0, 0.2]")
    return model_name, mode, cap


def _candidate_tag(candidate: tuple[str, str, float]) -> str:
    model_name, mode, cap = candidate
    return f"{model_name}_{mode}_cap{cap:.3f}".replace(".", "p")


def _finite(value: Any, default: float = 0.0) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError):
        return float(default)
    return result if math.isfinite(result) else float(default)


def _load_scenes(pair_path: Path, cutoff: str) -> list[dict[str, Any]]:
    with np.load(pair_path, allow_pickle=False) as pairs:
        scenes = json.loads(str(pairs["scenes_json"].item()))
    eligible = [scene for scene in scenes if str(scene.get("day", "")) <= cutoff]
    if not eligible:
        raise ValueError(f"no scenes on or before {cutoff} in {pair_path}")
    return eligible


def _fetch_scene(ee: Any, grid: dict[str, Any], scene: dict[str, Any]) -> dict[str, np.ndarray]:
    from bestpixel._gee import get_patch

    system_index = str(scene.get("l2a", {}).get("system_index", ""))
    if not system_index:
        raise ValueError(f"scene has no L2A system_index: {scene.get('scene_id')}")
    asset_id = f"COPERNICUS/S2_SR_HARMONIZED/{system_index}"
    raw = get_patch(ee, asset_id, [*GEE_L2A_BANDS, *L2A_AUX_BANDS], grid=grid)
    cloud_score = get_patch(ee, str(scene["cloud_score_asset"]), ["cs"], grid=grid)[0]
    surface = _normalise_l2a_reflectance(raw[: len(CANONICAL_BANDS)])
    l2a_aot = _normalise_aux(raw[len(CANONICAL_BANDS)])
    l2a_tcwv = _normalise_aux(raw[len(CANONICAL_BANDS) + 1])
    scl = np.asarray(raw[len(CANONICAL_BANDS) + 2], dtype=np.int16)
    valid = (
        (np.asarray(cloud_score, dtype=np.float32) > 0.60)
        & np.isin(scl, LAND_SCL)
        & np.all(np.isfinite(surface), axis=0)
        & np.isfinite(l2a_aot)
        & np.isfinite(l2a_tcwv)
        & np.all((surface > 0.001) & (surface < 0.8), axis=0)
    )
    valid = ndimage.binary_erosion(valid, structure=np.ones((3, 3), dtype=bool), iterations=1)
    surface[:, ~valid] = np.nan
    return {"surface": surface, "l2a_aot": l2a_aot, "l2a_tcwv": l2a_tcwv, "valid": valid}


def _band_indices(mode: str) -> tuple[int, ...]:
    if mode == "all":
        return tuple(range(len(BAND_NAMES)))
    if mode == "blue":
        return (1,)
    if mode == "deepblue":
        return (0, 1)
    if mode == "visible":
        return VISIBLE_INDICES
    if mode == "solver":
        return (1, 2, 3)
    raise ValueError(f"unsupported mode {mode!r}")


def _is_sklearn_artifact(model: dict[str, Any]) -> bool:
    return str(model.get("model_type", "")) == "hist_gradient_boosting"


_SKLEARN_MODEL_CACHE: dict[str, Any] = {}


def _sklearn_band_model(
    artifact_dir: Path, model: dict[str, Any], fold: int, band_name: str
) -> Any:
    file_name = model["model_files"][str(fold)][band_name]
    path = artifact_dir / "models" / file_name
    key = str(path)
    if key not in _SKLEARN_MODEL_CACHE:
        _SKLEARN_MODEL_CACHE[key] = joblib.load(path)
    return _SKLEARN_MODEL_CACHE[key]


def _sklearn_features(
    model: dict[str, Any],
    flat_selected: np.ndarray,
    metadata: dict[str, np.ndarray],
    band_name: str,
    scene_state: dict[str, float] | None = None,
) -> np.ndarray:
    """Assemble the nonlinear state matrix in the artifact's declared column order."""
    names = (model.get("feature_names_by_band") or {}).get(band_name) or model["feature_names"]
    count = flat_selected.shape[0]
    extra = {name: np.full(count, float(value)) for name, value in (scene_state or {}).items()}
    columns: dict[str, np.ndarray] = {
        **{f"l2a_{name}": flat_selected[:, index] for index, name in enumerate(BAND_NAMES)},
        "maiac_aot": metadata["maiac_aot"],
        "sen2cor_aot": metadata["l2a_aot"],
        "l2a_tcwv_cm": metadata["l2a_tcwv"],
        "maiac_tcwv_cm": metadata["maiac_tcwv"],
        "sza_deg": metadata["sza_deg"],
        "vza_deg": metadata["vza_deg"],
        "raa_deg": metadata["raa_deg"],
        "cos_raa": np.cos(np.radians(np.asarray(metadata["raa_deg"], dtype=np.float64))),
        "month": metadata["month"],
        "elevation_km": metadata["elevation_km"],
        "terrain_elevation_km": metadata["terrain_elevation_km"],
        "terrain_slope_deg": np.clip(
            np.asarray(metadata["terrain_slope_deg"], dtype=np.float64), 0.0, 70.0
        ),
        "terrain_incidence_cos": np.clip(
            np.asarray(metadata["terrain_incidence_cos"], dtype=np.float64), -1.0, 1.0
        ),
        "sensor_is_s2b": metadata["sensor_is_s2b"],
        "sensor_is_s2c": metadata["sensor_is_s2c"],
        "processing_baseline": metadata["processing_baseline"],
        **extra,
    }
    try:
        return np.column_stack([np.asarray(columns[name], dtype=np.float64) for name in names])
    except KeyError as exc:
        raise ValueError(f"artifact declares an unknown feature column: {exc}") from exc


def _correct_scene(
    fetched: dict[str, np.ndarray],
    scene: dict[str, Any],
    *,
    model: dict[str, Any],
    candidate: tuple[str, str, float],
    model_scope: str,
    fold: int,
    application_maiac_aot_max: float | None = None,
    terrain: TerrainFields | None = None,
    artifact_dir: Path | None = None,
    scene_state: dict[str, float] | None = None,
) -> tuple[np.ndarray, dict[str, Any]]:
    model_name, mode, cap = candidate
    is_sklearn = _is_sklearn_artifact(model)
    if is_sklearn:
        if model_name != str(model.get("model_name")):
            raise ValueError(f"sklearn artifact does not contain model {model_name!r}")
        if artifact_dir is None:
            raise ValueError("applying a sklearn artifact requires its directory")
        candidate_model = None
        uses_terrain = True
    else:
        candidate_model = model["models"][model_name]
        uses_terrain = feature_set_uses_terrain(str(candidate_model["feature_set"]))
    surface = np.asarray(fetched["surface"], dtype=np.float64).copy()
    valid = np.asarray(fetched["valid"], dtype=bool)
    if uses_terrain:
        if terrain is None:
            raise ValueError(f"terrain model {model_name!r} requires a DEM field")
        if terrain.valid.shape != valid.shape:
            raise ValueError("terrain grid does not match fetched L2A surface")
        valid &= terrain.valid
    flat = surface.reshape(len(BAND_NAMES), -1).T
    selected = np.flatnonzero(valid.ravel())
    if selected.size == 0:
        return surface.astype(np.float32), {"valid_pixels": 0, "bands": {}}
    count = selected.size
    scene_maiac_aot = _finite(scene.get("maiac_aot"), default=float("nan"))
    if application_maiac_aot_max is not None and (
        not math.isfinite(scene_maiac_aot) or scene_maiac_aot > float(application_maiac_aot_max)
    ):
        return surface.astype(np.float32), {
            "valid_pixels": int(count),
            "bands": {},
            "mapping_applied": False,
            "skip_reason": "outside_training_maiac_aot_domain",
            "scene_maiac_aot": scene_maiac_aot,
            "application_maiac_aot_max": float(application_maiac_aot_max),
        }
    metadata = {
        "l2a_aot": fetched["l2a_aot"].ravel()[selected],
        "l2a_tcwv": fetched["l2a_tcwv"].ravel()[selected],
        "maiac_aot": np.full(count, scene_maiac_aot),
        "maiac_tcwv": np.full(count, _finite(scene.get("maiac_tcwv_cm"))),
        "sza_deg": np.full(count, _finite(scene.get("sza_deg"))),
        "vza_deg": np.full(count, _finite(scene.get("vza_deg"))),
        "raa_deg": np.full(count, _finite(scene.get("raa_deg"))),
        "elevation_km": np.full(count, _finite(scene.get("elevation_km"))),
        "month": np.full(count, _finite(scene.get("month"))),
        "sensor_is_s2b": np.full(
            count,
            float("2B" in str((scene.get("l2a") or {}).get("spacecraft", "")).upper()),
        ),
        "sensor_is_s2c": np.full(
            count,
            float("2C" in str((scene.get("l2a") or {}).get("spacecraft", "")).upper()),
        ),
        "processing_baseline": np.full(
            count, _finite((scene.get("l2a") or {}).get("processing_baseline"))
        ),
    }
    if uses_terrain:
        if terrain is None:  # pragma: no cover - narrowed above for type checkers
            raise RuntimeError("missing terrain fields")
        incidence = local_solar_incidence(
            terrain,
            sza_deg=float(scene.get("sza_deg", 0.0)),
            saa_deg=float(scene.get("saa_deg", 0.0)),
        )
        metadata.update(
            {
                "terrain_elevation_km": terrain.elevation_m.ravel()[selected] * np.float32(1.0e-3),
                "terrain_slope_deg": terrain.slope_deg.ravel()[selected],
                "terrain_incidence_cos": incidence.ravel()[selected],
            }
        )
    correction_stats: dict[str, Any] = {}
    if is_sklearn:
        if model_scope != "crossfit":
            raise ValueError("the sklearn artifact only ships cross-fit fold models")
        # Feature columns must reflect the ORIGINAL L2A values for every band,
        # so snapshot the selected reflectance before any band is corrected.
        original_selected = flat[selected].copy()
        for band_index in _band_indices(mode):
            band_name = BAND_NAMES[band_index]
            sklearn_features = _sklearn_features(
                model, original_selected, metadata, band_name, scene_state
            )
            estimator = _sklearn_band_model(artifact_dir, model, fold, band_name)
            correction = np.clip(estimator.predict(sklearn_features), -cap, cap)
            flat[selected, band_index] = np.clip(
                flat[selected, band_index] + correction, 0.001, 0.8
            )
            correction_stats[BAND_NAMES[band_index]] = {
                "median": float(np.median(correction)),
                "median_abs": float(np.median(np.abs(correction))),
                "p95_abs": float(np.percentile(np.abs(correction), 95.0)),
                "cap_fraction": float(np.mean(np.abs(correction) >= cap - 1.0e-12)),
            }
    else:
        model_group = (
            candidate_model["full"]
            if model_scope == "full"
            else candidate_model["folds"][str(fold)]
        )
        for band_index in _band_indices(mode):
            features = build_features(
                flat[selected],
                band_index=band_index,
                feature_set=str(candidate_model["feature_set"]),
                **metadata,
            )
            band_model = model_group[BAND_NAMES[band_index]]
            if band_model.get("feature_names") != feature_names(
                str(candidate_model["feature_set"]), band_index
            ):
                raise ValueError(
                    f"feature schema mismatch for {model_name}/{BAND_NAMES[band_index]}"
                )
            correction = np.clip(_predict_model(features, band_model), -cap, cap)
            flat[selected, band_index] = np.clip(
                flat[selected, band_index] + correction, 0.001, 0.8
            )
            correction_stats[BAND_NAMES[band_index]] = {
                "median": float(np.median(correction)),
                "median_abs": float(np.median(np.abs(correction))),
                "p95_abs": float(np.percentile(np.abs(correction), 95.0)),
                "cap_fraction": float(np.mean(np.abs(correction) >= cap - 1.0e-12)),
            }
    return surface.astype(np.float32), {
        "valid_pixels": int(count),
        "bands": correction_stats,
        "mapping_applied": True,
        "scene_maiac_aot": scene_maiac_aot,
        "application_maiac_aot_max": application_maiac_aot_max,
        "terrain": {
            "enabled": uses_terrain,
            "source": TERRAIN_SOURCE if uses_terrain else None,
            "median_slope_deg": float(np.median(metadata["terrain_slope_deg"]))
            if uses_terrain
            else None,
            "median_incidence_cos": float(np.median(metadata["terrain_incidence_cos"]))
            if uses_terrain
            else None,
        },
    }


def _nanmedian(stack: list[np.ndarray]) -> np.ndarray:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        return np.nanmedian(np.stack(stack, axis=0), axis=0).astype(np.float32)


def _write_history(
    destination: Path,
    composites: list[np.ndarray],
    windows: list[str],
    *,
    grid: dict[str, Any],
    scene_year: int,
    scene_month: int,
    provenance: dict[str, Any],
) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_suffix(".tmp.npz")
    np.savez_compressed(
        temporary,
        comp=np.stack(composites, axis=0),
        epsg=np.asarray(int(grid["epsg"])),
        transform=np.asarray(
            [grid["res"], 0.0, grid["x0"], 0.0, -grid["res"], grid["y1"]],
            dtype=float,
        ),
        realizations=np.asarray(windows),
        scene_year=np.asarray(scene_year),
        month=np.asarray(scene_month),
        harmonization_json=np.asarray(json.dumps(provenance)),
    )
    temporary.replace(destination)


def _write_uncorrected_fallback(
    matchup_id: str,
    *,
    outputs: dict[str, Path],
    audit_path: Path,
    model: dict[str, Any],
    cutoff: str,
    model_scope: str,
    skip_reason: str,
    trigger_details: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Retain the same L2A prior when a correction cannot be supported robustly."""
    source_path = FALLBACK_HISTORY / f"{matchup_id}.npz"
    with np.load(source_path, allow_pickle=False) as source:
        source_payload = {key: source[key] for key in source.files}
    realizations = np.asarray(source_payload["realizations"]).astype(str)
    keep = realizations <= cutoff[:7]
    retained_count = int(np.count_nonzero(keep))
    if retained_count < 1:
        status = {
            "status": "FAILED",
            "matchup_id": matchup_id,
            "uses_aeronet": False,
            "reason": "no production-valid fallback history windows",
        }
        audit_path.parent.mkdir(parents=True, exist_ok=True)
        audit_path.write_text(json.dumps(status, indent=2) + "\n", encoding="utf-8")
        return status

    fold = int(model["fold_by_matchup_id"][matchup_id])
    windows = realizations[keep].tolist()
    common_provenance = {
        "uses_aeronet": False,
        "application": "uncorrected single-source L2A fallback",
        "mapping_applied": False,
        "skip_reason": skip_reason,
        "source_history": str(source_path),
        "training_cutoff": model.get("training_cutoff"),
        "history_cutoff": cutoff,
        "model_scope": model_scope,
        "crossfit_fold": fold if model_scope == "crossfit" else None,
        "model_target": model.get("target"),
        "low_temporal_support": retained_count < 2,
    }
    for name, destination in outputs.items():
        payload = dict(source_payload)
        payload["comp"] = np.asarray(source_payload["comp"])[keep]
        payload["realizations"] = realizations[keep]
        payload["harmonization_json"] = np.asarray(
            json.dumps(
                {
                    **common_provenance,
                    "candidate": None if name == "identity" else name,
                }
            )
        )
        destination.parent.mkdir(parents=True, exist_ok=True)
        temporary = destination.with_suffix(".tmp.npz")
        np.savez_compressed(temporary, **payload)
        temporary.replace(destination)
    status = {
        "status": "OK",
        "matchup_id": matchup_id,
        "uses_aeronet": False,
        "windows": windows,
        "window_count": len(windows),
        "scene_count": 0,
        "errors": [],
        "outputs": {name: str(path) for name, path in outputs.items()},
        "fallback_trigger": trigger_details,
        **common_provenance,
    }
    audit_path.parent.mkdir(parents=True, exist_ok=True)
    audit_path.write_text(json.dumps(status, indent=2) + "\n", encoding="utf-8")
    return status


def build_one(
    matchup_id: str,
    row: dict[str, str],
    *,
    pair_dir: Path,
    output_dir: Path,
    model: dict[str, Any],
    candidates: list[tuple[str, str, float]],
    cutoff: str,
    model_scope: str,
    application_maiac_aot_max: float | None,
    force: bool,
    artifact_dir: Path | None = None,
    cams_state: dict[tuple[str, str], dict[str, float]] | None = None,
    deterrain_l2a: bool = False,
) -> dict[str, Any]:
    from bestpixel._gee import init_ee, utm_epsg_from_bbox, utm_grid

    outputs = {
        "identity": output_dir / "identity_daily" / f"{matchup_id}.npz",
        **{
            _candidate_tag(candidate): output_dir / _candidate_tag(candidate) / f"{matchup_id}.npz"
            for candidate in candidates
        },
    }
    audit_path = output_dir / "audits" / f"{matchup_id}.json"
    if not force and audit_path.exists() and all(path.exists() for path in outputs.values()):
        audit = json.loads(audit_path.read_text(encoding="utf-8"))
        if audit.get("status") == "OK":
            return audit

    pair_path = pair_dir / f"{matchup_id}.npz"
    pair_audit_path = pair_dir / f"{matchup_id}.json"
    if not pair_path.exists() and pair_audit_path.exists():
        pair_audit = json.loads(pair_audit_path.read_text(encoding="utf-8"))
        if pair_audit.get("status") == "DATA_UNAVAILABLE":
            return _write_uncorrected_fallback(
                matchup_id,
                outputs=outputs,
                audit_path=audit_path,
                model=model,
                cutoff=cutoff,
                model_scope=model_scope,
                skip_reason="no_exact_pair_atmospheric_support",
            )

    lon = float(row["longitude"])
    lat = float(row["latitude"])
    bbox = (
        lon - HALF_SIZE_DEGREES,
        lat - HALF_SIZE_DEGREES,
        lon + HALF_SIZE_DEGREES,
        lat + HALF_SIZE_DEGREES,
    )
    grid = utm_grid(bbox, utm_epsg_from_bbox(bbox), 60.0)
    scenes = _load_scenes(pair_path, cutoff)
    by_window: dict[str, dict[str, list[dict[str, Any]]]] = defaultdict(lambda: defaultdict(list))
    for scene in scenes:
        by_window[str(scene["window"])][str(scene["day"])].append(scene)

    fold = int(model["fold_by_matchup_id"][matchup_id])
    products: dict[str, list[np.ndarray]] = {name: [] for name in outputs}
    windows: list[str] = []
    scene_audits: list[dict[str, Any]] = []
    errors: list[dict[str, str]] = []
    ee = init_ee()
    needs_terrain = (
        deterrain_l2a
        or _is_sklearn_artifact(model)
        or any(
            feature_set_uses_terrain(str(model["models"][name]["feature_set"]))
            for name, _mode, _cap in candidates
        )
    )
    terrain = fetch_glo30_terrain(ee, grid) if needs_terrain else None
    fraction_lut = None
    deterrain_sensor_cache: dict[str, Any] = {}
    if deterrain_l2a:
        from tools.aeronet_validation.build_l2a_l1c_harmonization_pairs import TARGET_LUT
        from tools.aeronet_validation.terrain_deshade import DirectFractionLUT

        fraction_lut = DirectFractionLUT(TARGET_LUT)
    for window in sorted(by_window):
        window_days: dict[str, list[np.ndarray]] = {name: [] for name in outputs}
        for day in sorted(by_window[window]):
            day_tiles: dict[str, list[np.ndarray]] = {name: [] for name in outputs}
            for scene in by_window[window][day]:
                try:
                    fetched = _fetch_scene(ee, grid, scene)
                except Exception as exc:  # noqa: BLE001 - retain usable acquisitions
                    errors.append(
                        {
                            "window": window,
                            "day": day,
                            "scene_id": str(scene.get("scene_id", "")),
                            "error_type": type(exc).__name__,
                            "reason": str(exc)[:500],
                        }
                    )
                    continue
                if int(np.count_nonzero(fetched["valid"])) < 100:
                    errors.append(
                        {
                            "window": window,
                            "day": day,
                            "scene_id": str(scene.get("scene_id", "")),
                            "error_type": "SparseScene",
                            "reason": "fewer than 100 clear-land pixels",
                        }
                    )
                    continue
                if deterrain_l2a:
                    from tools.aeronet_validation.build_l2a_l1c_harmonization_pairs import (
                        _sentinel_satellite_id,
                    )
                    from tools.aeronet_validation.terrain_deshade import deterrain_l2a_surface

                    scene_tco3 = 0.30
                    if cams_state is not None:
                        state_row = cams_state.get((matchup_id, str(scene.get("day")))) or {}
                        if "ozone_du_cams" in state_row:
                            scene_tco3 = float(state_row["ozone_du_cams"]) / 1000.0
                    deshaded, deterrain_provenance = deterrain_l2a_surface(
                        fetched["surface"],
                        l2a_aot=fetched["l2a_aot"],
                        l2a_tcwv=fetched["l2a_tcwv"],
                        terrain=terrain,
                        sza_deg=_finite(scene.get("sza_deg")),
                        saa_deg=_finite(scene.get("saa_deg")),
                        elevation_km=_finite(scene.get("elevation_km")),
                        satellite_id=_sentinel_satellite_id(
                            (scene.get("l2a") or {}).get("spacecraft")
                        ),
                        fraction_lut=fraction_lut,
                        sensor_cache=deterrain_sensor_cache,
                        band_names=BAND_NAMES,
                        s2_band_map=CANONICAL_TO_S2_BAND,
                        tco3_atm_cm=scene_tco3,
                    )
                    fetched = {**fetched, "surface": deshaded}
                else:
                    deterrain_provenance = {"applied": False}
                day_tiles["identity"].append(fetched["surface"])
                scene_audit: dict[str, Any] = {
                    "window": window,
                    "day": day,
                    "scene_id": scene.get("scene_id"),
                    "l2a_deterrain": deterrain_provenance,
                    "maiac_aot": scene.get("maiac_aot"),
                    "l2a_aot_median": scene.get("l2a_aot_median"),
                    "candidates": {},
                }
                scene_state = None
                if cams_state is not None:
                    key = (matchup_id, str(scene.get("day")))
                    if key not in cams_state:
                        raise KeyError(f"CAMS scene-state table is missing {key}")
                    scene_state = cams_state[key]
                for candidate in candidates:
                    tag = _candidate_tag(candidate)
                    corrected, stats = _correct_scene(
                        fetched,
                        scene,
                        model=model,
                        candidate=candidate,
                        model_scope=model_scope,
                        fold=fold,
                        application_maiac_aot_max=application_maiac_aot_max,
                        terrain=terrain,
                        artifact_dir=artifact_dir,
                        scene_state=scene_state,
                    )
                    day_tiles[tag].append(corrected)
                    scene_audit["candidates"][tag] = stats
                scene_audits.append(scene_audit)
            if not day_tiles["identity"]:
                continue
            for name, tile_surfaces in day_tiles.items():
                window_days[name].append(_nanmedian(tile_surfaces))
        if not window_days["identity"]:
            continue
        windows.append(window)
        for name, day_surfaces in window_days.items():
            products[name].append(_nanmedian(day_surfaces))

    if len(windows) < 2:
        return _write_uncorrected_fallback(
            matchup_id,
            outputs=outputs,
            audit_path=audit_path,
            model=model,
            cutoff=cutoff,
            model_scope=model_scope,
            skip_reason="insufficient_harmonized_history_windows",
            trigger_details={
                "harmonized_windows": windows,
                "harmonized_scene_count": len(scene_audits),
                "errors": errors,
            },
        )

    sensing = row["sensing_time_utc"]
    common_provenance = {
        "uses_aeronet": False,
        "application": "per acquisition before tile mosaic and temporal median",
        "l2a_deterrain": bool(deterrain_l2a),
        "training_cutoff": model.get("training_cutoff"),
        "history_cutoff": cutoff,
        "model_scope": model_scope,
        "crossfit_fold": fold if model_scope == "crossfit" else None,
        "model_target": model.get("target"),
        "application_maiac_aot_max": application_maiac_aot_max,
        "terrain_features": {
            "enabled": needs_terrain,
            "source": TERRAIN_SOURCE if needs_terrain else None,
            "fields": ["elevation_km", "slope_deg", "solar_incidence_cos"] if needs_terrain else [],
        },
    }
    for name, destination in outputs.items():
        provenance = {
            **common_provenance,
            "candidate": None if name == "identity" else name,
        }
        _write_history(
            destination,
            products[name],
            windows,
            grid=grid,
            scene_year=int(sensing[:4]),
            scene_month=int(sensing[5:7]),
            provenance=provenance,
        )
    status = {
        "status": "OK",
        "matchup_id": matchup_id,
        "uses_aeronet": False,
        "windows": windows,
        "window_count": len(windows),
        "scene_count": len(scene_audits),
        "errors": errors,
        "outputs": {name: str(path) for name, path in outputs.items()},
        "scenes": scene_audits,
        **common_provenance,
    }
    audit_path.parent.mkdir(parents=True, exist_ok=True)
    audit_path.write_text(json.dumps(status, indent=2) + "\n", encoding="utf-8")
    return status


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("matchup_id", nargs="+")
    parser.add_argument("--pairs", type=Path, default=DEFAULT_PAIRS)
    parser.add_argument("--model", type=Path, default=DEFAULT_MODEL)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--cases", type=Path, default=DEFAULT_CASES)
    parser.add_argument("--history-cutoff", default="2023-12-31")
    parser.add_argument("--model-scope", choices=("crossfit", "full"), default="crossfit")
    parser.add_argument("--application-maiac-aot-max", type=float)
    parser.add_argument("--candidate", action="append", type=_candidate, default=[])
    parser.add_argument(
        "--cams-state",
        type=Path,
        default=None,
        help="scene-state table supplying extra sklearn feature columns (e.g. real ozone)",
    )
    parser.add_argument(
        "--deterrain-l2a",
        action="store_true",
        help="analytically invert Sen2Cor's rugged-terrain correction on every acquisition",
    )
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    model = json.loads(args.model.read_text(encoding="utf-8"))
    if model.get("uses_aeronet") is not False:
        raise SystemExit("refusing a model without uses_aeronet=false provenance")
    if args.application_maiac_aot_max is not None and not (
        0.0 < args.application_maiac_aot_max <= 1.0
    ):
        raise SystemExit("--application-maiac-aot-max must be in (0, 1]")
    candidates = args.candidate or [
        ("compact_a10", "all", 0.03),
        ("full_a10", "all", 0.03),
    ]
    cams_state: dict[tuple[str, str], dict[str, float]] | None = None
    if _is_sklearn_artifact(model):
        if args.model_scope != "crossfit":
            raise SystemExit("the sklearn artifact only supports --model-scope crossfit")
        known_models = {str(model.get("model_name"))}
        declared = {
            name for names in (model.get("feature_names_by_band") or {}).values() for name in names
        }
        if "ozone_du_cams" in declared:
            if args.cams_state is None:
                raise SystemExit("this artifact needs --cams-state for its ozone column")
            from tools.aeronet_validation.train_l2a_l1c_nonlinear_harmonizer import (
                SEN2COR_OZONE_NODES,
            )

            def _scene_ozone(value: float) -> dict[str, float]:
                node = min(SEN2COR_OZONE_NODES, key=lambda n: abs(n - value))
                return {"ozone_du_cams": value, "ozone_du_sen2cor_node": float(node)}

            cams_state = {
                (row["matchup_id"], row["day"]): _scene_ozone(float(row["cams_tco3_du"]))
                for row in csv.DictReader(args.cams_state.open())
            }
    else:
        known_models = set(model["models"])
    unknown_models = sorted({name for name, _mode, _cap in candidates} - known_models)
    if unknown_models:
        raise SystemExit(f"model artifact does not contain: {', '.join(unknown_models)}")
    requested = set(args.matchup_id)
    case_ids = {row["matchup_id"] for row in csv.DictReader(args.cases.open())}
    if not requested <= case_ids:
        raise SystemExit("requested matchup is outside the locked development case list")
    rows = {row["matchup_id"]: row for row in csv.DictReader(MATCHUPS.open())}
    exit_code = 0
    for matchup_id in args.matchup_id:
        status = build_one(
            matchup_id,
            rows[matchup_id],
            pair_dir=args.pairs,
            output_dir=args.output,
            model=model,
            candidates=candidates,
            cutoff=args.history_cutoff,
            model_scope=args.model_scope,
            application_maiac_aot_max=args.application_maiac_aot_max,
            force=args.force,
            artifact_dir=args.model.parent,
            cams_state=cams_state,
            deterrain_l2a=args.deterrain_l2a,
        )
        print("HISTORY_DONE " + json.dumps(status), flush=True)
        if status["status"] != "OK":
            exit_code = 1
    raise SystemExit(exit_code)


if __name__ == "__main__":
    main()
