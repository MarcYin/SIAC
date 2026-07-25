"""Cross-fit a generic L2A-to-current-RT surface-reflectance harmonizer.

Inputs are exact same-day L2A/L1C pairs produced without AERONET. Models are
grouped by site during cross-fitting and predict a per-band reflectance residual
from L2A reflectance, MAIAC/Sen2Cor AOD state, matching L2A water vapour,
geometry, month, and elevation. The ``cross-rt`` family separates the equal-AOD
RT-frame baseline from AOD-difference, water-vapour/geometry, and terrain terms. Its
regression weights give every acquisition equal influence, rather than treating
spatially correlated pixels as independent scenes. Coefficients are exported as
JSON and can be applied without scikit-learn.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from dataclasses import dataclass
from datetime import date
from pathlib import Path
from typing import Any, Literal

import numpy as np
from sklearn.linear_model import Ridge
from sklearn.model_selection import GroupKFold
from sklearn.preprocessing import StandardScaler
from tools.aeronet_validation.build_l2a_l1c_harmonization_pairs import TARGET_RT

ROOT = Path("/gws/ssde/j25a/nceo_isp/public/siac_refactor")
DEFAULT_PAIRS = ROOT / "analysis/l2a_l1c_exact_pairs_mediumdev_20260713"
DEFAULT_HISTORY = ROOT / "l2a_seasonal_pc"
DEFAULT_OUTPUT = ROOT / "analysis/l2a_l1c_harmonizer_mediumdev_20260713"
DEFAULT_CASES = ROOT / "reports/aod-medium-physics-20260713/downloads/cases.csv"
BAND_NAMES = ("coastal", "blue", "green", "red", "nir08", "swir16", "swir22")
VISIBLE_INDICES = (0, 1, 2, 3)
CURRENT_OFFSETS = {
    1: (-0.0003, 0.0243),
    2: (-0.0006, 0.0235),
    3: (-0.0011, 0.0223),
}


@dataclass(frozen=True)
class PairDataset:
    l2a: np.ndarray
    siac: np.ndarray
    l2a_aot: np.ndarray
    l2a_tcwv: np.ndarray
    maiac_aot: np.ndarray
    maiac_tcwv: np.ndarray
    sza_deg: np.ndarray
    vza_deg: np.ndarray
    raa_deg: np.ndarray
    elevation_km: np.ndarray
    month: np.ndarray
    sensor_is_s2b: np.ndarray
    sensor_is_s2c: np.ndarray
    processing_baseline: np.ndarray
    terrain_elevation_km: np.ndarray | None
    terrain_slope_deg: np.ndarray | None
    terrain_incidence_cos: np.ndarray | None
    terrain_features: bool
    site_code: np.ndarray
    scene_code: np.ndarray
    acquisition_code: np.ndarray
    matchup_by_site_code: dict[int, list[str]]
    scene_metadata: list[dict[str, Any]]
    target: str
    target_rt: dict[str, Any]


@dataclass(frozen=True)
class ModelSpec:
    name: str
    feature_set: Literal[
        "compact",
        "full",
        "terrain",
        "cross_rt_baseline",
        "cross_rt_aod",
        "cross_rt_atmosphere",
        "cross_rt_terrain",
    ]
    alpha: float
    weighting: Literal["pixel", "scene"] = "pixel"


MODEL_SPECS = (
    ModelSpec("compact_a1", "compact", 1.0),
    ModelSpec("compact_a10", "compact", 10.0),
    ModelSpec("full_a10", "full", 10.0),
    ModelSpec("full_a100", "full", 100.0),
)
TERRAIN_MODEL_SPECS = (ModelSpec("terrain_a100", "terrain", 100.0),)
CROSS_RT_MODEL_SPECS = (
    ModelSpec("cross_rt_baseline_a1", "cross_rt_baseline", 1.0, "scene"),
    ModelSpec("cross_rt_aod_a1", "cross_rt_aod", 1.0, "scene"),
    ModelSpec("cross_rt_atmosphere_a1", "cross_rt_atmosphere", 1.0, "scene"),
)
CROSS_RT_TERRAIN_MODEL_SPECS = (
    ModelSpec("cross_rt_terrain_a1", "cross_rt_terrain", 1.0, "scene"),
    ModelSpec("cross_rt_terrain_a10", "cross_rt_terrain", 10.0, "scene"),
)


def _finite(value: Any, default: float = 0.0) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError):
        return float(default)
    return result if math.isfinite(result) else float(default)


def _airmass(sza_deg: np.ndarray, vza_deg: np.ndarray) -> np.ndarray:
    sza = np.radians(np.asarray(sza_deg, dtype=np.float64))
    vza = np.radians(np.asarray(vza_deg, dtype=np.float64))
    return 1.0 / np.maximum(np.cos(sza), 0.1) + 1.0 / np.maximum(np.cos(vza), 0.1)


def feature_set_uses_terrain(feature_set: str) -> bool:
    """Return whether a feature schema requires local DEM-derived fields."""
    return str(feature_set) in {"terrain", "cross_rt_terrain"}


def feature_names(feature_set: str, band_index: int) -> list[str]:
    band = BAND_NAMES[band_index]
    cross_rt_baseline = [
        f"l2a_{band}",
        "mean_aot",
        f"l2a_{band}_x_mean_aot",
        "l2a_tcwv_cm",
        f"l2a_{band}_x_l2a_tcwv",
        "solar_airmass",
        "view_airmass",
        f"l2a_{band}_x_total_airmass",
        "cos_raa",
        "elevation_km",
        "month_sin",
        "month_cos",
        "sensor_is_s2b",
        "sensor_is_s2c",
        "processing_baseline",
    ]
    cross_rt_aod = [
        *cross_rt_baseline,
        "delta_aot_maiac_minus_sen2cor",
        f"l2a_{band}_x_delta_aot",
        "delta_aot_x_total_airmass",
        "delta_aot_x_mean_aot",
        "signed_delta_aot_squared",
    ]
    cross_rt_atmosphere = [
        *cross_rt_aod,
        "l2a_tcwv_x_solar_airmass",
        f"l2a_{band}_x_l2a_tcwv_x_total_airmass",
        "delta_aot_x_l2a_tcwv",
        "delta_aot_x_cos_raa",
    ]
    if feature_set == "cross_rt_baseline":
        return cross_rt_baseline
    if feature_set == "cross_rt_aod":
        return cross_rt_aod
    if feature_set == "cross_rt_atmosphere":
        return cross_rt_atmosphere
    if feature_set == "cross_rt_terrain":
        return [
            *cross_rt_atmosphere,
            "terrain_elevation_km",
            "terrain_elevation_delta_km",
            "terrain_slope_deg",
            "terrain_incidence_cos",
            "terrain_incidence_delta",
            "terrain_illumination_ratio_minus_one",
            "terrain_shadow",
            f"l2a_{band}_x_terrain_incidence_delta",
            f"l2a_{band}_x_terrain_illumination_ratio_minus_one",
            "delta_aot_x_terrain_incidence_delta",
        ]
    compact = [
        f"l2a_{band}",
        "delta_aot_maiac_minus_sen2cor",
        f"l2a_{band}_x_delta_aot",
        "delta_tcwv_maiac_minus_sen2cor",
        f"l2a_{band}_x_delta_tcwv",
        "airmass",
        f"l2a_{band}_x_airmass",
        "cos_raa",
        "elevation_km",
        "month_sin",
        "month_cos",
        "sensor_is_s2b",
        "sensor_is_s2c",
        "processing_baseline",
    ]
    if feature_set == "compact":
        return compact
    full = [
        *[f"l2a_{name}" for name in BAND_NAMES],
        "maiac_aot",
        "sen2cor_aot",
        "delta_aot_maiac_minus_sen2cor",
        "maiac_tcwv_cm",
        "sen2cor_tcwv_cm",
        "delta_tcwv_maiac_minus_sen2cor",
        "airmass",
        "cos_raa",
        "elevation_km",
        "month_sin",
        "month_cos",
        "sensor_is_s2b",
        "sensor_is_s2c",
        "processing_baseline",
        f"l2a_{band}_x_delta_aot",
        f"l2a_{band}_x_delta_tcwv",
        f"l2a_{band}_x_airmass",
        "delta_aot_x_airmass",
    ]
    if feature_set == "full":
        return full
    if feature_set == "terrain":
        return [
            *full,
            "terrain_elevation_km",
            "terrain_slope_deg",
            "terrain_incidence_delta",
            f"l2a_{band}_x_terrain_incidence_delta",
            f"l2a_{band}_x_terrain_slope_deg",
        ]
    raise ValueError(f"unsupported feature set {feature_set!r}")


def build_features(
    l2a: np.ndarray,
    *,
    band_index: int,
    feature_set: str,
    l2a_aot: np.ndarray,
    l2a_tcwv: np.ndarray,
    maiac_aot: np.ndarray,
    maiac_tcwv: np.ndarray,
    sza_deg: np.ndarray,
    vza_deg: np.ndarray,
    raa_deg: np.ndarray,
    elevation_km: np.ndarray,
    month: np.ndarray,
    sensor_is_s2b: np.ndarray,
    sensor_is_s2c: np.ndarray,
    processing_baseline: np.ndarray,
    terrain_elevation_km: np.ndarray | None = None,
    terrain_slope_deg: np.ndarray | None = None,
    terrain_incidence_cos: np.ndarray | None = None,
) -> np.ndarray:
    reflectance = np.asarray(l2a, dtype=np.float64)
    target = reflectance[:, int(band_index)]
    l2a_aot = np.asarray(l2a_aot, dtype=np.float64)
    l2a_tcwv = np.asarray(l2a_tcwv, dtype=np.float64)
    maiac_aot = np.asarray(maiac_aot, dtype=np.float64)
    maiac_tcwv = np.asarray(maiac_tcwv, dtype=np.float64)
    delta_aot = maiac_aot - l2a_aot
    delta_tcwv = maiac_tcwv - l2a_tcwv
    air = _airmass(sza_deg, vza_deg)
    sza_rad = np.radians(np.asarray(sza_deg, dtype=np.float64))
    vza_rad = np.radians(np.asarray(vza_deg, dtype=np.float64))
    solar_air = 1.0 / np.maximum(np.cos(sza_rad), 0.1)
    view_air = 1.0 / np.maximum(np.cos(vza_rad), 0.1)
    cos_raa = np.cos(np.radians(np.asarray(raa_deg, dtype=np.float64)))
    angle = 2.0 * np.pi * np.asarray(month, dtype=np.float64) / 12.0
    mean_aot = 0.5 * (maiac_aot + l2a_aot)
    cross_rt_baseline = np.column_stack(
        [
            target,
            mean_aot,
            target * mean_aot,
            l2a_tcwv,
            target * l2a_tcwv,
            solar_air,
            view_air,
            target * air,
            cos_raa,
            np.asarray(elevation_km, dtype=np.float64),
            np.sin(angle),
            np.cos(angle),
            np.asarray(sensor_is_s2b, dtype=np.float64),
            np.asarray(sensor_is_s2c, dtype=np.float64),
            np.asarray(processing_baseline, dtype=np.float64),
        ]
    )
    if feature_set == "cross_rt_baseline":
        return cross_rt_baseline
    cross_rt_aod = np.column_stack(
        [
            cross_rt_baseline,
            delta_aot,
            target * delta_aot,
            delta_aot * air,
            delta_aot * mean_aot,
            delta_aot * np.abs(delta_aot),
        ]
    )
    if feature_set == "cross_rt_aod":
        return cross_rt_aod
    cross_rt_atmosphere = np.column_stack(
        [
            cross_rt_aod,
            l2a_tcwv * solar_air,
            target * l2a_tcwv * air,
            delta_aot * l2a_tcwv,
            delta_aot * cos_raa,
        ]
    )
    if feature_set == "cross_rt_atmosphere":
        return cross_rt_atmosphere
    if feature_set == "cross_rt_terrain":
        if (
            terrain_elevation_km is None
            or terrain_slope_deg is None
            or terrain_incidence_cos is None
        ):
            raise ValueError("cross-RT terrain features require per-pixel terrain fields")
        local_elevation = np.asarray(terrain_elevation_km, dtype=np.float64)
        local_slope = np.clip(np.asarray(terrain_slope_deg, dtype=np.float64), 0.0, 70.0)
        incidence = np.clip(np.asarray(terrain_incidence_cos, dtype=np.float64), -1.0, 1.0)
        flat_incidence = np.cos(sza_rad)
        incidence_delta = incidence - flat_incidence
        illumination_ratio = np.clip(
            flat_incidence / np.maximum(incidence, 0.15) - 1.0,
            -2.0,
            4.0,
        )
        shadow = (incidence <= 0.05).astype(np.float64)
        return np.column_stack(
            [
                cross_rt_atmosphere,
                local_elevation,
                local_elevation - np.asarray(elevation_km, dtype=np.float64),
                local_slope,
                incidence,
                incidence_delta,
                illumination_ratio,
                shadow,
                target * incidence_delta,
                target * illumination_ratio,
                delta_aot * incidence_delta,
            ]
        )
    common = [
        target,
        delta_aot,
        target * delta_aot,
        delta_tcwv,
        target * delta_tcwv,
        air,
        target * air,
        cos_raa,
        np.asarray(elevation_km, dtype=np.float64),
        np.sin(angle),
        np.cos(angle),
        np.asarray(sensor_is_s2b, dtype=np.float64),
        np.asarray(sensor_is_s2c, dtype=np.float64),
        np.asarray(processing_baseline, dtype=np.float64),
    ]
    if feature_set == "compact":
        return np.column_stack(common)
    full = np.column_stack(
        [
            reflectance,
            maiac_aot,
            l2a_aot,
            delta_aot,
            maiac_tcwv,
            l2a_tcwv,
            delta_tcwv,
            air,
            cos_raa,
            np.asarray(elevation_km, dtype=np.float64),
            np.sin(angle),
            np.cos(angle),
            np.asarray(sensor_is_s2b, dtype=np.float64),
            np.asarray(sensor_is_s2c, dtype=np.float64),
            np.asarray(processing_baseline, dtype=np.float64),
            target * delta_aot,
            target * delta_tcwv,
            target * air,
            delta_aot * air,
        ]
    )
    if feature_set == "full":
        return full
    if feature_set != "terrain":
        raise ValueError(f"unsupported feature set {feature_set!r}")
    if terrain_elevation_km is None or terrain_slope_deg is None or terrain_incidence_cos is None:
        raise ValueError("terrain feature set requires per-pixel terrain fields")
    local_elevation = np.asarray(terrain_elevation_km, dtype=np.float64)
    local_slope = np.clip(np.asarray(terrain_slope_deg, dtype=np.float64), 0.0, 70.0)
    flat_incidence = np.cos(np.radians(np.asarray(sza_deg, dtype=np.float64)))
    incidence_delta = np.clip(
        np.asarray(terrain_incidence_cos, dtype=np.float64) - flat_incidence,
        -1.0,
        1.0,
    )
    return np.column_stack(
        [
            full,
            local_elevation,
            local_slope,
            incidence_delta,
            target * incidence_delta,
            target * local_slope,
        ]
    )


def _case_ids(path: Path) -> list[str]:
    return [row["matchup_id"] for row in csv.DictReader(path.open())]


def _acquisition_identity(scene: dict[str, Any], local_index: int) -> str:
    """Identify one overpass while allowing it to contain adjacent MGRS tiles."""
    system_index = str((scene.get("l2a") or {}).get("system_index", "")).strip()
    if system_index:
        # GEE indices start with the sensing timestamp; the remaining fields
        # include generation time and tile, which differ across one overpass.
        return system_index.split("_", 1)[0]
    scene_id = str(scene.get("scene_id", "")).strip()
    return scene_id or f"local_scene_{local_index}"


def load_pairs(
    pair_dir: Path,
    matchup_ids: list[str],
    *,
    scene_day_max: str | None = None,
    maiac_aot_max: float | None = None,
    max_samples_per_scene: int | None = None,
    allow_missing_matchups: bool = False,
) -> PairDataset:
    l2a_parts: list[np.ndarray] = []
    siac_parts: list[np.ndarray] = []
    l2a_aot_parts: list[np.ndarray] = []
    l2a_tcwv_parts: list[np.ndarray] = []
    terrain_parts: dict[str, list[np.ndarray]] = {
        "terrain_elevation_km": [],
        "terrain_slope_deg": [],
        "terrain_incidence_cos": [],
    }
    terrain_available: bool | None = None
    scalar_parts: dict[str, list[np.ndarray]] = {
        name: []
        for name in (
            "maiac_aot",
            "maiac_tcwv",
            "sza_deg",
            "vza_deg",
            "raa_deg",
            "elevation_km",
            "month",
            "sensor_is_s2b",
            "sensor_is_s2c",
            "processing_baseline",
            "site_code",
            "scene_code",
            "acquisition_code",
        )
    }
    site_by_name: dict[str, int] = {}
    matchup_by_site_code: dict[int, list[str]] = {}
    scene_metadata: list[dict[str, Any]] = []
    loaded_matchups: set[str] = set()
    pair_targets: dict[str, dict[str, Any]] = {}
    next_scene_code = 0
    acquisition_by_site_and_id: dict[tuple[int, str], int] = {}
    for matchup_id in matchup_ids:
        archive_path = pair_dir / f"{matchup_id}.npz"
        audit_path = pair_dir / f"{matchup_id}.json"
        if not audit_path.exists():
            raise FileNotFoundError(f"missing pair audit for {matchup_id}")
        audit = json.loads(audit_path.read_text(encoding="utf-8"))
        if audit.get("status") == "DATA_UNAVAILABLE" and allow_missing_matchups:
            continue
        if not archive_path.exists():
            raise FileNotFoundError(f"missing pair archive for {matchup_id}")
        if audit.get("status") != "OK" or audit.get("uses_aeronet") is not False:
            raise RuntimeError(f"pair archive is not a clean terminal result: {audit_path}")
        target = str(
            audit.get("target") or "same-day L1C corrected at MAIAC AOD with current libRadtran LUT"
        )
        target_rt = audit.get("target_rt") or TARGET_RT
        pair_targets.setdefault(
            json.dumps({"target": target, "target_rt": target_rt}, sort_keys=True),
            {"target": target, "target_rt": target_rt},
        )
        with np.load(archive_path, allow_pickle=False) as pair:
            l2a = np.asarray(pair["l2a"], dtype=np.float32)
            siac = np.asarray(pair["siac"], dtype=np.float32)
            l2a_aot = np.asarray(pair["l2a_aot"], dtype=np.float32)
            l2a_tcwv = np.asarray(pair["l2a_tcwv"], dtype=np.float32)
            local_scene = np.asarray(pair["scene_index"], dtype=np.int32)
            scenes = json.loads(str(pair["scenes_json"].item()))
            bands = tuple(str(value) for value in pair["band_names"])
            has_terrain = all(name in pair.files for name in terrain_parts)
            if any(name in pair.files for name in terrain_parts) and not has_terrain:
                raise ValueError(f"partial terrain fields in {archive_path}")
            terrain_values = (
                {name: np.asarray(pair[name], dtype=np.float32) for name in terrain_parts}
                if has_terrain
                else {}
            )
        if bands != BAND_NAMES or l2a.shape != siac.shape or l2a.shape[1] != len(BAND_NAMES):
            raise ValueError(f"invalid pair shape or bands in {archive_path}")
        if terrain_available is None:
            terrain_available = has_terrain
        elif terrain_available != has_terrain:
            raise ValueError("cannot mix terrain-enabled and legacy pair archives")
        if has_terrain and any(
            values.shape != (l2a.shape[0],) for values in terrain_values.values()
        ):
            raise ValueError(f"invalid terrain feature shape in {archive_path}")
        if np.any(local_scene < 0) or np.any(local_scene >= len(scenes)):
            raise ValueError(f"scene_index is outside scenes_json in {archive_path}")
        allowed_scene = np.ones(len(scenes), dtype=bool)
        if scene_day_max is not None:
            scene_days = [str(scene.get("day", "")) for scene in scenes]
            if any(not day_value for day_value in scene_days):
                raise ValueError(f"scene metadata has no day in {archive_path}")
            allowed_scene &= np.asarray(
                [day_value <= scene_day_max for day_value in scene_days], dtype=bool
            )
        if maiac_aot_max is not None:
            scene_aot = np.asarray(
                [_finite(scene.get("maiac_aot"), default=float("nan")) for scene in scenes]
            )
            allowed_scene &= np.isfinite(scene_aot) & (scene_aot <= float(maiac_aot_max))
        if not np.all(allowed_scene):
            keep = allowed_scene[local_scene]
            l2a = l2a[keep]
            siac = siac[keep]
            l2a_aot = l2a_aot[keep]
            l2a_tcwv = l2a_tcwv[keep]
            if has_terrain:
                terrain_values = {name: values[keep] for name, values in terrain_values.items()}
            local_scene = local_scene[keep]
        if max_samples_per_scene is not None:
            if max_samples_per_scene < 1:
                raise ValueError("max_samples_per_scene must be positive")
            retained: list[np.ndarray] = []
            for local_index in np.unique(local_scene):
                indices = np.flatnonzero(local_scene == local_index)
                if indices.size > int(max_samples_per_scene):
                    positions = np.linspace(
                        0,
                        indices.size - 1,
                        num=int(max_samples_per_scene),
                        dtype=np.int64,
                    )
                    indices = indices[positions]
                retained.append(indices)
            selected = np.sort(np.concatenate(retained)) if retained else np.asarray([], dtype=int)
            l2a = l2a[selected]
            siac = siac[selected]
            l2a_aot = l2a_aot[selected]
            l2a_tcwv = l2a_tcwv[selected]
            if has_terrain:
                terrain_values = {name: values[selected] for name, values in terrain_values.items()}
            local_scene = local_scene[selected]
        if local_scene.size == 0:
            continue
        loaded_matchups.add(matchup_id)
        site_name = matchup_id.split("__", 1)[0]
        site_code = site_by_name.setdefault(site_name, len(site_by_name))
        matchup_by_site_code.setdefault(site_code, []).append(matchup_id)
        count = l2a.shape[0]
        l2a_parts.append(l2a)
        siac_parts.append(siac)
        l2a_aot_parts.append(l2a_aot)
        l2a_tcwv_parts.append(l2a_tcwv)
        if has_terrain:
            for name, values in terrain_values.items():
                terrain_parts[name].append(values)
        scalar_parts["site_code"].append(np.full(count, site_code, dtype=np.int16))
        global_scene = np.full(count, -1, dtype=np.int32)
        global_acquisition = np.full(count, -1, dtype=np.int32)
        archive_scalars = {
            name: np.empty(count, dtype=np.float32)
            for name in (
                "maiac_aot",
                "maiac_tcwv",
                "sza_deg",
                "vza_deg",
                "raa_deg",
                "elevation_km",
                "month",
                "sensor_is_s2b",
                "sensor_is_s2c",
                "processing_baseline",
            )
        }
        for local_index, scene in enumerate(scenes):
            selected = local_scene == local_index
            scene_count = int(np.count_nonzero(selected))
            if scene_count == 0:
                continue
            global_scene[selected] = next_scene_code
            acquisition_id = _acquisition_identity(scene, local_index)
            acquisition_key = (site_code, acquisition_id)
            if acquisition_key not in acquisition_by_site_and_id:
                acquisition_by_site_and_id[acquisition_key] = len(acquisition_by_site_and_id)
            acquisition_code = acquisition_by_site_and_id[acquisition_key]
            global_acquisition[selected] = acquisition_code
            scene_record = {
                "matchup_id": matchup_id,
                "site": site_name,
                "acquisition_id": acquisition_id,
                "acquisition_code": acquisition_code,
                **scene,
            }
            scene_metadata.append(scene_record)
            for name, key in (
                ("maiac_aot", "maiac_aot"),
                ("maiac_tcwv", "maiac_tcwv_cm"),
                ("sza_deg", "sza_deg"),
                ("vza_deg", "vza_deg"),
                ("raa_deg", "raa_deg"),
                ("elevation_km", "elevation_km"),
                ("month", "month"),
            ):
                value = _finite(scene.get(key))
                archive_scalars[name][selected] = value
            l2a_metadata = scene.get("l2a") or {}
            spacecraft = str(l2a_metadata.get("spacecraft", "")).upper()
            archive_scalars["sensor_is_s2b"][selected] = float("2B" in spacecraft)
            archive_scalars["sensor_is_s2c"][selected] = float("2C" in spacecraft)
            archive_scalars["processing_baseline"][selected] = _finite(
                l2a_metadata.get("processing_baseline")
            )
            next_scene_code += 1
        if np.any(global_scene < 0):
            raise ValueError(f"scene_index does not match scenes_json in {archive_path}")
        if np.any(global_acquisition < 0):
            raise ValueError(f"acquisition identity does not match samples in {archive_path}")
        for name, values in archive_scalars.items():
            scalar_parts[name].append(values)
        scalar_parts["scene_code"].append(global_scene)
        scalar_parts["acquisition_code"].append(global_acquisition)
    missing_matchups = sorted(set(matchup_ids) - loaded_matchups)
    if missing_matchups and not allow_missing_matchups:
        raise ValueError(
            f"no eligible exact-pair samples for {len(missing_matchups)} matchup(s): "
            + ", ".join(missing_matchups)
        )
    if len(pair_targets) != 1:
        raise ValueError("cannot mix pair archives built in different target RT frames")
    pair_provenance = next(iter(pair_targets.values()))
    has_terrain = bool(terrain_available)
    return PairDataset(
        l2a=np.concatenate(l2a_parts),
        siac=np.concatenate(siac_parts),
        l2a_aot=np.concatenate(l2a_aot_parts),
        l2a_tcwv=np.concatenate(l2a_tcwv_parts),
        maiac_aot=np.concatenate(scalar_parts["maiac_aot"]),
        maiac_tcwv=np.concatenate(scalar_parts["maiac_tcwv"]),
        sza_deg=np.concatenate(scalar_parts["sza_deg"]),
        vza_deg=np.concatenate(scalar_parts["vza_deg"]),
        raa_deg=np.concatenate(scalar_parts["raa_deg"]),
        elevation_km=np.concatenate(scalar_parts["elevation_km"]),
        month=np.concatenate(scalar_parts["month"]),
        sensor_is_s2b=np.concatenate(scalar_parts["sensor_is_s2b"]),
        sensor_is_s2c=np.concatenate(scalar_parts["sensor_is_s2c"]),
        processing_baseline=np.concatenate(scalar_parts["processing_baseline"]),
        terrain_elevation_km=(
            np.concatenate(terrain_parts["terrain_elevation_km"]) if has_terrain else None
        ),
        terrain_slope_deg=(
            np.concatenate(terrain_parts["terrain_slope_deg"]) if has_terrain else None
        ),
        terrain_incidence_cos=(
            np.concatenate(terrain_parts["terrain_incidence_cos"]) if has_terrain else None
        ),
        terrain_features=has_terrain,
        site_code=np.concatenate(scalar_parts["site_code"]),
        scene_code=np.concatenate(scalar_parts["scene_code"]),
        acquisition_code=np.concatenate(scalar_parts["acquisition_code"]),
        matchup_by_site_code=matchup_by_site_code,
        scene_metadata=scene_metadata,
        target=str(pair_provenance["target"]),
        target_rt=dict(pair_provenance["target_rt"]),
    )


def _feature_kwargs(dataset: PairDataset) -> dict[str, np.ndarray]:
    values: dict[str, np.ndarray] = {
        "l2a_aot": dataset.l2a_aot,
        "l2a_tcwv": dataset.l2a_tcwv,
        "maiac_aot": dataset.maiac_aot,
        "maiac_tcwv": dataset.maiac_tcwv,
        "sza_deg": dataset.sza_deg,
        "vza_deg": dataset.vza_deg,
        "raa_deg": dataset.raa_deg,
        "elevation_km": dataset.elevation_km,
        "month": dataset.month,
        "sensor_is_s2b": dataset.sensor_is_s2b,
        "sensor_is_s2c": dataset.sensor_is_s2c,
        "processing_baseline": dataset.processing_baseline,
    }
    if dataset.terrain_features:
        if (
            dataset.terrain_elevation_km is None
            or dataset.terrain_slope_deg is None
            or dataset.terrain_incidence_cos is None
        ):
            raise RuntimeError("terrain-enabled dataset has missing terrain arrays")
        values.update(
            {
                "terrain_elevation_km": dataset.terrain_elevation_km,
                "terrain_slope_deg": dataset.terrain_slope_deg,
                "terrain_incidence_cos": dataset.terrain_incidence_cos,
            }
        )
    return values


def _scene_balanced_sample_weight(scene_code: np.ndarray) -> np.ndarray:
    """Give every independent acquisition the same total regression weight."""
    scene = np.asarray(scene_code, dtype=np.int64)
    if scene.ndim != 1 or scene.size == 0 or np.any(scene < 0):
        raise ValueError("scene codes must be a non-empty vector of non-negative integers")
    counts = np.bincount(scene)
    if np.any(counts[scene] <= 0):  # pragma: no cover - guarded by bincount semantics
        raise RuntimeError("scene balancing encountered an empty scene")
    return 1.0 / counts[scene].astype(np.float64)


def _fit_model(
    features: np.ndarray,
    residual: np.ndarray,
    alpha: float,
    *,
    sample_weight: np.ndarray | None = None,
) -> dict[str, Any]:
    weights = None if sample_weight is None else np.asarray(sample_weight, dtype=np.float64)
    scaler = StandardScaler().fit(features, sample_weight=weights)
    transformed = scaler.transform(features)
    model = Ridge(alpha=float(alpha)).fit(
        transformed,
        residual,
        sample_weight=weights,
    )
    return {
        "mean": scaler.mean_.tolist(),
        "scale": scaler.scale_.tolist(),
        "coef": np.asarray(model.coef_, dtype=float).tolist(),
        "intercept": float(model.intercept_),
        "effective_weight_sum": float(np.sum(weights))
        if weights is not None
        else float(features.shape[0]),
    }


def _predict_model(features: np.ndarray, model: dict[str, Any]) -> np.ndarray:
    mean = np.asarray(model["mean"], dtype=np.float64)
    scale = np.asarray(model["scale"], dtype=np.float64)
    coef = np.asarray(model["coef"], dtype=np.float64)
    return ((features - mean) / scale) @ coef + float(model["intercept"])


def _metrics(error: np.ndarray, group_code: np.ndarray) -> dict[str, float]:
    """Return pixel metrics plus equal-weight median errors for each supplied group."""
    finite = np.isfinite(error)
    values = np.asarray(error[finite], dtype=np.float64)
    scenes = np.asarray(group_code[finite], dtype=np.int32)
    if values.size and np.any(scenes[1:] < scenes[:-1]):
        order = np.argsort(scenes, kind="stable")
        values = values[order]
        scenes = scenes[order]
    boundaries = np.flatnonzero(scenes[1:] != scenes[:-1]) + 1
    starts = np.concatenate(([0], boundaries))
    stops = np.concatenate((boundaries, [values.size]))
    scene_medians = np.asarray(
        [np.median(values[start:stop]) for start, stop in zip(starts, stops, strict=True)],
        dtype=np.float64,
    )
    return {
        "n": int(values.size),
        "bias": float(np.mean(values)),
        "mae": float(np.mean(np.abs(values))),
        "rmse": float(np.sqrt(np.mean(np.square(values)))),
        "p95_abs": float(np.percentile(np.abs(values), 95.0)),
        "scene_bias": float(np.mean(scene_medians)),
        "scene_mae": float(np.mean(np.abs(scene_medians))),
        "scene_rmse": float(np.sqrt(np.mean(np.square(scene_medians)))),
    }


def _distribution(values: np.ndarray) -> dict[str, float | int]:
    """Summarize acquisition-level support without pixel-count weighting."""
    finite = np.asarray(values, dtype=np.float64)
    finite = finite[np.isfinite(finite)]
    if finite.size == 0:
        return {"n": 0, "p05": float("nan"), "median": float("nan"), "p95": float("nan")}
    return {
        "n": int(finite.size),
        "p05": float(np.percentile(finite, 5.0)),
        "median": float(np.median(finite)),
        "p95": float(np.percentile(finite, 95.0)),
    }


def _scene_rows(dataset: PairDataset) -> list[dict[str, Any]]:
    counts = np.bincount(dataset.scene_code, minlength=len(dataset.scene_metadata))
    terrain_means: dict[str, np.ndarray] = {}
    if dataset.terrain_features:
        terrain_values = {
            "terrain_elevation_km_mean": dataset.terrain_elevation_km,
            "terrain_slope_deg_mean": dataset.terrain_slope_deg,
            "terrain_incidence_cos_mean": dataset.terrain_incidence_cos,
        }
        if any(values is None for values in terrain_values.values()):
            raise RuntimeError("terrain-enabled dataset has missing scene summaries")
        for name, values in terrain_values.items():
            array = np.asarray(values, dtype=np.float64)
            finite = np.isfinite(array)
            terrain_means[name] = np.divide(
                np.bincount(
                    dataset.scene_code[finite],
                    weights=array[finite],
                    minlength=len(dataset.scene_metadata),
                ),
                np.bincount(dataset.scene_code[finite], minlength=len(dataset.scene_metadata)),
                out=np.full(len(dataset.scene_metadata), np.nan),
                where=np.bincount(dataset.scene_code[finite], minlength=len(dataset.scene_metadata))
                > 0,
            )
    rows: list[dict[str, Any]] = []
    for scene_code, metadata in enumerate(dataset.scene_metadata):
        row = {
            "scene_code": scene_code,
            "matchup_id": metadata.get("matchup_id"),
            "site": metadata.get("site"),
            "window": metadata.get("window"),
            "day": metadata.get("day"),
            "scene_id": metadata.get("scene_id"),
            "sample_count": int(counts[scene_code]),
            "spacecraft": (metadata.get("l2a") or {}).get("spacecraft"),
            "processing_baseline": _finite(
                (metadata.get("l2a") or {}).get("processing_baseline"),
                default=float("nan"),
            ),
        }
        for output_name, metadata_name in (
            ("maiac_aot", "maiac_aot"),
            ("l2a_aot", "l2a_aot_median"),
            ("maiac_tcwv_cm", "maiac_tcwv_cm"),
            ("l2a_tcwv_cm", "l2a_tcwv_median_cm"),
            ("sza_deg", "sza_deg"),
            ("vza_deg", "vza_deg"),
            ("raa_deg", "raa_deg"),
            ("elevation_km", "elevation_km"),
        ):
            row[output_name] = _finite(metadata.get(metadata_name), default=float("nan"))
        row["delta_aot_maiac_minus_sen2cor"] = row["maiac_aot"] - row["l2a_aot"]
        row["delta_tcwv_maiac_minus_sen2cor"] = row["maiac_tcwv_cm"] - row["l2a_tcwv_cm"]
        for name, values in terrain_means.items():
            row[name] = float(values[scene_code])
        rows.append(row)
    return rows


def _add_scene_errors(
    rows: list[dict[str, Any]],
    *,
    prefix: str,
    error: np.ndarray,
    scene_code: np.ndarray,
) -> None:
    n_scenes = len(rows)
    for band_index, band_name in enumerate(BAND_NAMES):
        values = np.asarray(error[:, band_index], dtype=np.float64)
        finite = np.isfinite(values)
        counts = np.bincount(scene_code[finite], minlength=n_scenes)
        bias = np.divide(
            np.bincount(scene_code[finite], weights=values[finite], minlength=n_scenes),
            counts,
            out=np.full(n_scenes, np.nan),
            where=counts > 0,
        )
        mae = np.divide(
            np.bincount(scene_code[finite], weights=np.abs(values[finite]), minlength=n_scenes),
            counts,
            out=np.full(n_scenes, np.nan),
            where=counts > 0,
        )
        for code, row in enumerate(rows):
            row[f"{prefix}_{band_name}_bias"] = float(bias[code])
            row[f"{prefix}_{band_name}_mae"] = float(mae[code])


def _fixed_fold_splits(
    raw_folds: np.ndarray,
    *,
    holdout_folds: set[int] | None,
) -> tuple[np.ndarray, list[tuple[int, np.ndarray, np.ndarray]], dict[str, Any]]:
    """Build OOF partitions, optionally isolating a locked holdout from all fitting."""
    raw = np.asarray(raw_folds, dtype=np.int16)
    fold_ids = sorted(int(value) for value in np.unique(raw) if value >= 0)
    if not fold_ids:
        raise ValueError("fixed fold map has no populated folds")
    if not holdout_folds:
        splits = [
            (fold, np.flatnonzero(raw != fold), np.flatnonzero(raw == fold)) for fold in fold_ids
        ]
        return raw.copy(), splits, {"mode": "site_oof", "source_folds": fold_ids}

    holdout = {int(value) for value in holdout_folds}
    unknown = sorted(holdout - set(fold_ids))
    if unknown:
        raise ValueError(f"holdout folds are absent from the fixed fold map: {unknown}")
    development = [fold for fold in fold_ids if fold not in holdout]
    if len(development) < 2:
        raise ValueError("locked evaluation requires at least two development folds")
    holdout_model_fold = max(fold_ids) + 1
    development_mask = np.isin(raw, development)
    holdout_mask = np.isin(raw, sorted(holdout))
    application_folds = raw.copy()
    application_folds[holdout_mask] = holdout_model_fold
    splits = [
        (
            fold,
            np.flatnonzero(development_mask & (raw != fold)),
            np.flatnonzero(raw == fold),
        )
        for fold in development
    ]
    splits.append(
        (
            holdout_model_fold,
            np.flatnonzero(development_mask),
            np.flatnonzero(holdout_mask),
        )
    )
    return (
        application_folds,
        splits,
        {
            "mode": "locked_development_holdout",
            "source_folds": fold_ids,
            "development_folds": development,
            "holdout_folds": sorted(holdout),
            "holdout_model_fold": holdout_model_fold,
        },
    )


def crossfit_models(
    dataset: PairDataset,
    *,
    training_cutoff: str | None,
    training_maiac_aot_max: float | None = None,
    fixed_fold_by_matchup_id: dict[str, int] | None = None,
    holdout_folds: set[int] | None = None,
    model_family: Literal["legacy", "cross-rt", "all"] = "legacy",
) -> tuple[dict[str, Any], dict[str, Any], list[dict[str, Any]]]:
    if model_family not in {"legacy", "cross-rt", "all"}:
        raise ValueError(f"unsupported model family {model_family!r}")
    raw_folds = np.full(dataset.site_code.shape, -1, dtype=np.int16)
    fold_by_site: dict[int, int] = {}
    if fixed_fold_by_matchup_id is None:
        if holdout_folds:
            raise ValueError("holdout_folds requires fixed_fold_by_matchup_id")
        grouped = list(GroupKFold(n_splits=5).split(dataset.l2a, groups=dataset.site_code))
        split_rows = [(fold, train, test) for fold, (train, test) in enumerate(grouped)]
        for fold, _train, test in split_rows:
            raw_folds[test] = fold
            for site_code in np.unique(dataset.site_code[test]):
                fold_by_site[int(site_code)] = fold
        artifact_fold_map = {
            matchup_id: fold_by_site[site_code]
            for site_code, matchup_ids in dataset.matchup_by_site_code.items()
            for matchup_id in matchup_ids
        }
        application_folds = raw_folds.copy()
        split_protocol = {"mode": "site_oof", "source_folds": list(range(5))}
        source_fold_map = artifact_fold_map.copy()
    else:
        source_fold_map = {key: int(value) for key, value in fixed_fold_by_matchup_id.items()}
        for site_code, matchup_ids in dataset.matchup_by_site_code.items():
            assigned = {source_fold_map[matchup_id] for matchup_id in matchup_ids}
            if len(assigned) != 1:
                raise ValueError(f"fixed fold map splits one site across folds: {matchup_ids}")
            fold_by_site[int(site_code)] = assigned.pop()
        for site_code, fold in fold_by_site.items():
            raw_folds[dataset.site_code == site_code] = fold
        application_folds, split_rows, split_protocol = _fixed_fold_splits(
            raw_folds,
            holdout_folds=holdout_folds,
        )
        fold_remap = {
            int(raw): int(applied)
            for raw, applied in zip(raw_folds, application_folds, strict=True)
        }
        artifact_fold_map = {
            matchup_id: fold_remap[source_fold]
            for matchup_id, source_fold in source_fold_map.items()
        }
        if any(train.size == 0 or test.size == 0 for _fold, train, test in split_rows):
            raise ValueError("fixed fold map leaves an empty clean-day train or test partition")
    if np.any(raw_folds < 0) or np.any(application_folds < 0):
        raise RuntimeError("cross-fitting did not assign every sample")

    artifact: dict[str, Any] = {
        "schema_version": 2,
        "feature_schema_version": 5
        if model_family in {"cross-rt", "all"}
        else (3 if dataset.terrain_features else 2),
        "uses_aeronet": False,
        "training_cutoff": training_cutoff,
        "training_maiac_aot_max": training_maiac_aot_max,
        "target": dataset.target,
        "target_rt": dataset.target_rt,
        "model_family": model_family,
        "bands": list(BAND_NAMES),
        "terrain_features": {
            "enabled": dataset.terrain_features,
            "source": "COPERNICUS/DEM/GLO30_2024_1" if dataset.terrain_features else None,
            "fields": ["elevation_km", "slope_deg", "solar_incidence_cos"]
            if dataset.terrain_features
            else [],
        },
        "fold_by_matchup_id": artifact_fold_map,
        "source_fold_by_matchup_id": source_fold_map,
        "crossfit_protocol": split_protocol,
        "models": {},
        "regression_weighting_unit": "independent Sentinel-2 acquisition",
    }
    metrics: dict[str, Any] = {
        "sample_count": int(dataset.l2a.shape[0]),
        "scene_count": len(dataset.scene_metadata),
        "acquisition_count": int(np.unique(dataset.acquisition_code).size),
        "site_count": len(dataset.matchup_by_site_code),
        "training_cutoff": training_cutoff,
        "training_maiac_aot_max": training_maiac_aot_max,
        "model_family": model_family,
        "crossfit_protocol": split_protocol,
        "scene_metric_unit": "independent Sentinel-2 acquisition",
        "identity": {},
        "current_aod_offset": {},
        "candidates": {},
    }
    scene_rows = _scene_rows(dataset)
    scene_delta_aot = np.asarray(
        [row["delta_aot_maiac_minus_sen2cor"] for row in scene_rows], dtype=np.float64
    )
    metrics["scene_distribution"] = {
        "maiac_aot": _distribution(
            np.asarray([row["maiac_aot"] for row in scene_rows], dtype=np.float64)
        ),
        "sen2cor_aot": _distribution(
            np.asarray([row["l2a_aot"] for row in scene_rows], dtype=np.float64)
        ),
        "delta_aot_maiac_minus_sen2cor": _distribution(scene_delta_aot),
        "delta_tcwv_maiac_minus_sen2cor": _distribution(
            np.asarray(
                [row["delta_tcwv_maiac_minus_sen2cor"] for row in scene_rows],
                dtype=np.float64,
            )
        ),
        "near_equal_aod_abs_delta_le_0p02": int(
            np.count_nonzero(np.isfinite(scene_delta_aot) & (np.abs(scene_delta_aot) <= 0.02))
        ),
    }
    identity_error = dataset.l2a - dataset.siac
    _add_scene_errors(
        scene_rows,
        prefix="identity",
        error=identity_error,
        scene_code=dataset.scene_code,
    )
    current_error = np.empty_like(identity_error, dtype=np.float64)
    for band_index, band_name in enumerate(BAND_NAMES):
        raw_error = identity_error[:, band_index]
        metrics["identity"][band_name] = _metrics(raw_error, dataset.acquisition_code)
        intercept, slope = CURRENT_OFFSETS.get(band_index, (0.0, 0.0))
        old_error = (
            dataset.l2a[:, band_index]
            + intercept
            + slope * dataset.maiac_aot
            - dataset.siac[:, band_index]
        )
        current_error[:, band_index] = old_error
        metrics["current_aod_offset"][band_name] = _metrics(
            old_error, dataset.acquisition_code
        )
    _add_scene_errors(
        scene_rows,
        prefix="current_aod_offset",
        error=current_error,
        scene_code=dataset.scene_code,
    )

    kwargs = _feature_kwargs(dataset)
    legacy_specs = (*MODEL_SPECS, *TERRAIN_MODEL_SPECS) if dataset.terrain_features else MODEL_SPECS
    cross_rt_specs = (
        (*CROSS_RT_MODEL_SPECS, *CROSS_RT_TERRAIN_MODEL_SPECS)
        if dataset.terrain_features
        else CROSS_RT_MODEL_SPECS
    )
    if model_family == "legacy":
        model_specs = legacy_specs
    elif model_family == "cross-rt":
        model_specs = cross_rt_specs
    else:
        model_specs = (*legacy_specs, *cross_rt_specs)
    # Border sites can contribute two adjacent tiles from one overpass. Weight
    # their combined pixels once so tile layout cannot change model influence.
    scene_weight = _scene_balanced_sample_weight(dataset.acquisition_code)
    model_fold_ids = [fold for fold, _train, _test in split_rows]
    for spec in model_specs:
        candidate_models: dict[str, Any] = {
            "folds": {str(fold): {} for fold in model_fold_ids},
            "full": {},
        }
        residual_prediction = np.full_like(dataset.siac, np.nan, dtype=np.float64)
        sample_weight = scene_weight if spec.weighting == "scene" else None
        for band_index, band_name in enumerate(BAND_NAMES):
            features = build_features(
                dataset.l2a,
                band_index=band_index,
                feature_set=spec.feature_set,
                **kwargs,
            )
            target = dataset.siac[:, band_index] - dataset.l2a[:, band_index]
            for fold, train, test in split_rows:
                model = _fit_model(
                    features[train],
                    target[train],
                    spec.alpha,
                    sample_weight=None if sample_weight is None else sample_weight[train],
                )
                model["feature_names"] = feature_names(spec.feature_set, band_index)
                residual_prediction[test, band_index] = _predict_model(features[test], model)
                candidate_models["folds"][str(fold)][band_name] = model
            full_model = _fit_model(
                features,
                target,
                spec.alpha,
                sample_weight=sample_weight,
            )
            full_model["feature_names"] = feature_names(spec.feature_set, band_index)
            candidate_models["full"][band_name] = full_model
        candidate_models["feature_set"] = spec.feature_set
        candidate_models["alpha"] = spec.alpha
        candidate_models["weighting"] = spec.weighting
        candidate_models["weighting_group"] = (
            "acquisition" if spec.weighting == "scene" else "pixel"
        )
        candidate_models["oof_residual_scale"] = {}
        for band_index, band_name in enumerate(BAND_NAMES):
            error = (
                dataset.l2a[:, band_index]
                + residual_prediction[:, band_index]
                - dataset.siac[:, band_index]
            )
            center = float(np.median(error))
            candidate_models["oof_residual_scale"][band_name] = {
                "median": center,
                "mad_to_sigma": float(np.median(np.abs(error - center)) * 1.4826),
                "rmse": float(np.sqrt(np.mean(np.square(error)))),
            }
        artifact["models"][spec.name] = candidate_models

        candidate_metrics: dict[str, Any] = {}
        for cap in (0.015, 0.03, 0.06, None):
            cap_key = "uncapped" if cap is None else f"cap_{cap:.3f}"
            corrected = (
                residual_prediction if cap is None else np.clip(residual_prediction, -cap, cap)
            )
            candidate_metrics[cap_key] = {
                band_name: _metrics(
                    dataset.l2a[:, band_index]
                    + corrected[:, band_index]
                    - dataset.siac[:, band_index],
                    dataset.acquisition_code,
                )
                for band_index, band_name in enumerate(BAND_NAMES)
            }
            corrected_error = dataset.l2a + corrected - dataset.siac
            _add_scene_errors(
                scene_rows,
                prefix=f"{spec.name}_{cap_key}".replace(".", "p"),
                error=corrected_error,
                scene_code=dataset.scene_code,
            )
        metrics["candidates"][spec.name] = candidate_metrics
    return artifact, metrics, scene_rows


def _window_values(scenes: list[dict[str, Any]], window: str, count: int) -> dict[str, np.ndarray]:
    selected = [scene for scene in scenes if str(scene.get("window")) == str(window)]
    if not selected:
        raise KeyError(f"no exact-pair metadata for window {window}")

    def values(key: str) -> np.ndarray:
        return np.full(count, np.median([_finite(scene.get(key)) for scene in selected]))

    def l2a_values(key: str) -> np.ndarray:
        return np.full(
            count,
            np.median([_finite((scene.get("l2a") or {}).get(key)) for scene in selected]),
        )

    return {
        "l2a_aot": values("l2a_aot_median"),
        "l2a_tcwv": values("l2a_tcwv_median_cm"),
        "maiac_aot": values("maiac_aot"),
        "maiac_tcwv": values("maiac_tcwv_cm"),
        "sza_deg": values("sza_deg"),
        "vza_deg": values("vza_deg"),
        "raa_deg": values("raa_deg"),
        "elevation_km": values("elevation_km"),
        "month": values("month"),
        "sensor_is_s2b": np.full(
            count,
            np.median(
                [
                    float("2B" in str((scene.get("l2a") or {}).get("spacecraft", "")).upper())
                    for scene in selected
                ]
            ),
        ),
        "sensor_is_s2c": np.full(
            count,
            np.median(
                [
                    float("2C" in str((scene.get("l2a") or {}).get("spacecraft", "")).upper())
                    for scene in selected
                ]
            ),
        ),
        "processing_baseline": l2a_values("processing_baseline"),
    }


def generate_history(
    *,
    matchup_id: str,
    history_dir: Path,
    pair_dir: Path,
    output_dir: Path,
    artifact: dict[str, Any],
    model_name: str,
    mode: str,
    correction_cap: float,
) -> None:
    source_path = history_dir / f"{matchup_id}.npz"
    pair_path = pair_dir / f"{matchup_id}.npz"
    with np.load(source_path, allow_pickle=False) as source:
        payload = {key: source[key] for key in source.files}
    with np.load(pair_path, allow_pickle=False) as pairs:
        scenes = json.loads(str(pairs["scenes_json"].item()))
    composites = np.asarray(payload["comp"], dtype=np.float64)
    realizations = [str(value) for value in payload["realizations"]]
    training_cutoff = artifact.get("training_cutoff")
    if training_cutoff is not None:
        history_keep = np.asarray(
            [window <= str(training_cutoff)[:7] for window in realizations], dtype=bool
        )
        composites = composites[history_keep]
        payload["realizations"] = np.asarray(realizations)[history_keep]
        realizations = [
            window for window, keep in zip(realizations, history_keep, strict=True) if keep
        ]
    if not realizations:
        raise ValueError(f"no production-valid historical periods in {source_path}")
    fold = int(artifact["fold_by_matchup_id"][matchup_id])
    model = artifact["models"][model_name]
    if mode == "all":
        band_indices = range(len(BAND_NAMES))
    elif mode == "visible":
        band_indices = VISIBLE_INDICES
    elif mode == "solver":
        band_indices = (1, 2, 3)
    else:
        raise ValueError(f"unknown harmonization mode {mode!r}")

    correction_abs: list[float] = []
    for realization_index, window in enumerate(realizations):
        flat = composites[realization_index].reshape(len(BAND_NAMES), -1).T
        valid = np.all(np.isfinite(flat) & (flat > 0.001) & (flat < 0.8), axis=1)
        if not np.any(valid):
            continue
        metadata = _window_values(scenes, window, flat.shape[0])
        for band_index in band_indices:
            features = build_features(
                flat,
                band_index=band_index,
                feature_set=str(model["feature_set"]),
                **metadata,
            )
            band_model = model["folds"][str(fold)][BAND_NAMES[band_index]]
            residual = np.clip(
                _predict_model(features[valid], band_model),
                -float(correction_cap),
                float(correction_cap),
            )
            flat[valid, band_index] = np.clip(flat[valid, band_index] + residual, 0.001, 0.8)
            correction_abs.extend(np.abs(residual).tolist())
        composites[realization_index] = flat.T.reshape(composites.shape[1:])
    payload["comp"] = composites.astype(np.float32)
    payload["harmonization_json"] = np.asarray(
        json.dumps(
            {
                "uses_aeronet": False,
                "model": model_name,
                "mode": mode,
                "correction_cap": float(correction_cap),
                "crossfit_fold": fold,
                "target": artifact["target"],
                "pair_source": str(pair_path),
                "median_abs_correction": float(np.median(correction_abs))
                if correction_abs
                else None,
            }
        )
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    destination = output_dir / source_path.name
    temporary = destination.with_suffix(".tmp.npz")
    np.savez_compressed(temporary, **payload)
    temporary.replace(destination)


def _parse_generate(value: str) -> tuple[str, str, float]:
    parts = value.split(":")
    if len(parts) != 3:
        raise argparse.ArgumentTypeError("--generate must be MODEL:MODE:CAP")
    model, mode, cap = parts
    known_specs = (
        *MODEL_SPECS,
        *TERRAIN_MODEL_SPECS,
        *CROSS_RT_MODEL_SPECS,
        *CROSS_RT_TERRAIN_MODEL_SPECS,
    )
    if model not in {spec.name for spec in known_specs}:
        raise argparse.ArgumentTypeError(f"unknown model {model!r}")
    if mode not in {"all", "blue", "visible", "solver"}:
        raise argparse.ArgumentTypeError(f"unknown mode {mode!r}")
    try:
        cap_value = float(cap)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("CAP must be numeric") from exc
    if not 0.0 < cap_value <= 0.2:
        raise argparse.ArgumentTypeError("CAP must be in (0, 0.2]")
    return model, mode, cap_value


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pairs", type=Path, default=DEFAULT_PAIRS)
    parser.add_argument("--history", type=Path, default=DEFAULT_HISTORY)
    parser.add_argument("--cases", type=Path, default=DEFAULT_CASES)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument(
        "--training-cutoff",
        default="2023-12-31",
        help="latest exact-pair day allowed in the backtest model; use 'none' for deployment",
    )
    parser.add_argument(
        "--training-maiac-aot-max",
        type=float,
        default=None,
        help="retain only exact-pair scenes at or below this MAIAC AOD",
    )
    parser.add_argument(
        "--fold-map",
        type=Path,
        default=None,
        help="optional harmonizer artifact whose site-group folds must be preserved",
    )
    parser.add_argument(
        "--evaluation-split-manifest",
        type=Path,
        default=None,
        help="locked retrieval summary providing source folds and holdout-fold IDs",
    )
    parser.add_argument(
        "--model-family",
        choices=("legacy", "cross-rt", "all"),
        default="legacy",
        help="fit legacy models, the physically structured cross-RT ablation, or both",
    )
    parser.add_argument(
        "--max-samples-per-scene",
        type=int,
        default=None,
        help="deterministically cap correlated pixels while retaining every acquisition",
    )
    parser.add_argument("--generate", action="append", type=_parse_generate, default=[])
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    training_cutoff = None if args.training_cutoff.lower() == "none" else args.training_cutoff
    if training_cutoff is not None:
        date.fromisoformat(training_cutoff)
    if args.training_maiac_aot_max is not None and args.training_maiac_aot_max <= 0.0:
        raise SystemExit("--training-maiac-aot-max must be positive")
    if args.max_samples_per_scene is not None and args.max_samples_per_scene < 10:
        raise SystemExit("--max-samples-per-scene must be at least 10")
    matchup_ids = _case_ids(args.cases)
    fixed_folds = None
    holdout_folds = None
    evaluation_provenance = None
    if args.fold_map is not None and args.evaluation_split_manifest is not None:
        raise SystemExit("--fold-map and --evaluation-split-manifest are mutually exclusive")
    if args.fold_map is not None:
        fixed_folds = {
            key: int(value)
            for key, value in json.loads(args.fold_map.read_text(encoding="utf-8"))[
                "fold_by_matchup_id"
            ].items()
        }
        missing_folds = sorted(set(matchup_ids) - set(fixed_folds))
        if missing_folds:
            raise SystemExit(f"fold map is missing {len(missing_folds)} requested matchups")
    elif args.evaluation_split_manifest is not None:
        split_payload = json.loads(args.evaluation_split_manifest.read_text(encoding="utf-8"))
        cohort = split_payload["cohort"]
        fixed_folds = {key: int(value) for key, value in cohort["fold_by_matchup_id"].items()}
        holdout_folds = {int(value) for value in cohort["holdout_folds"]}
        missing_folds = sorted(set(matchup_ids) - set(fixed_folds))
        if missing_folds:
            raise SystemExit(f"evaluation split is missing {len(missing_folds)} requested matchups")
        evaluation_provenance = {
            "manifest": str(args.evaluation_split_manifest),
            "split_seed": cohort.get("split_seed"),
            "holdout_folds": sorted(holdout_folds),
        }
    dataset = load_pairs(
        args.pairs,
        matchup_ids,
        scene_day_max=training_cutoff,
        maiac_aot_max=args.training_maiac_aot_max,
        max_samples_per_scene=args.max_samples_per_scene,
        allow_missing_matchups=fixed_folds is not None,
    )
    artifact, metrics, scene_rows = crossfit_models(
        dataset,
        training_cutoff=training_cutoff,
        training_maiac_aot_max=args.training_maiac_aot_max,
        fixed_fold_by_matchup_id=fixed_folds,
        holdout_folds=holdout_folds,
        model_family=args.model_family,
    )
    artifact["evaluation_split"] = evaluation_provenance
    metrics["evaluation_split"] = evaluation_provenance
    artifact["pair_archive"] = str(args.pairs)
    metrics["pair_archive"] = str(args.pairs)
    artifact["max_samples_per_scene"] = args.max_samples_per_scene
    metrics["max_samples_per_scene"] = args.max_samples_per_scene
    args.output.mkdir(parents=True, exist_ok=True)
    artifact_path = args.output / "harmonizer.json"
    metrics_path = args.output / "surface_metrics.json"
    artifact_path.write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")
    metrics_path.write_text(json.dumps(metrics, indent=2) + "\n", encoding="utf-8")
    scene_metrics_path = args.output / "surface_scene_metrics.csv"
    with scene_metrics_path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(scene_rows[0]))
        writer.writeheader()
        writer.writerows(scene_rows)
    generated: list[str] = []
    for model_name, mode, cap in args.generate:
        tag = f"{model_name}_{mode}_cap{cap:.3f}".replace(".", "p")
        destination = args.output / f"histories_{tag}"
        for matchup_id in matchup_ids:
            generate_history(
                matchup_id=matchup_id,
                history_dir=args.history,
                pair_dir=args.pairs,
                output_dir=destination,
                artifact=artifact,
                model_name=model_name,
                mode=mode,
                correction_cap=cap,
            )
        generated.append(str(destination))
    print(
        json.dumps(
            {
                "pairs": int(dataset.l2a.shape[0]),
                "scenes": len(dataset.scene_metadata),
                "acquisitions": metrics["acquisition_count"],
                "sites": len(dataset.matchup_by_site_code),
                "training_maiac_aot_max": args.training_maiac_aot_max,
                "model_family": args.model_family,
                "max_samples_per_scene": args.max_samples_per_scene,
                "artifact": str(artifact_path),
                "metrics": str(metrics_path),
                "scene_metrics": str(scene_metrics_path),
                "generated": generated,
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
