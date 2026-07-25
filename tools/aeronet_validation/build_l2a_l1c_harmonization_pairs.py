"""Build exact same-day L2A/SIAC-L1C surface-reflectance training pairs.

The target is not AERONET. For every retained historical Sentinel-2 day, this
tool fetches operational L2A and L1C TOA on the same 60 m grid, corrects L1C
with either saved coefficients or the current MAIAC-conditioned physical RT
state, and samples clear land pixels. The output carries atmospheric and acquisition
metadata to learn a generic L2A-to-current-RT harmonization without using the
validation AOD.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import json
import math
import time
from pathlib import Path
from typing import Any

import numpy as np
import xarray as xr
from scipy import ndimage
from tools.aeronet_validation.terrain_features import (
    TERRAIN_SOURCE,
    TerrainFields,
    fetch_glo30_terrain,
    local_solar_incidence,
)

ROOT = Path("/gws/ssde/j25a/nceo_isp/public/siac_refactor")
MATCHUPS = ROOT / "matchups" / "matchups.csv"
L2A_HISTORY = ROOT / "l2a_seasonal_pc"
SIDECAR_ROOT = ROOT / "l1c_test" / "seasonal"
ELEVATION_CACHE = ROOT / "site_elev_cache.json"
DEFAULT_OUTPUT = ROOT / "analysis/l2a_l1c_exact_pairs_mediumdev_20260713"
TARGET_LUT = "/gws/ssde/j25a/nceo_isp/public/libradtran_continental_average_lut_1nm.zarr.zip"
TARGET_RT = {
    "backend": "lut",
    "lut": TARGET_LUT,
    "aerosol_profile": "continental_average",
    "aot_source": "MAIAC",
    "tcwv_source": "scene metadata",
    "tco3_atm_cm": 0.30,
    "elevation_km": 0.05,
    "coefficient_source": "saved per-scene libRadtran LUT sidecar",
}
PHYSICAL_TARGET_RT = {
    "backend": "lut",
    "lut": TARGET_LUT,
    "aerosol_profile": "continental_average",
    "aot_source": "same-day MAIAC scene value",
    "tcwv_source": "matching operational L2A WVP per pixel",
    "tco3_atm_cm": 0.30,
    "elevation_source": "COPERNICUS/DEM/GLO30_2024_1 per pixel",
    "geometry_source": "matching Sentinel-2 acquisition mean angles",
    "coefficient_source": "current libRadtran spectral LUT",
}
PHYSICAL_V2_TARGET_RT = {
    "backend": "lut",
    "lut": TARGET_LUT,
    "aerosol_profile": "continental_average",
    "aot_source": "same-day MAIAC scene value",
    "tcwv_source": "matching operational L2A WVP per pixel",
    "tco3_source": "CAMS analysis column at the site and local overpass hour",
    "elevation_source": "COPERNICUS/DEM/GLO30_2024_1 per pixel",
    "geometry_source": "matching Sentinel-2 acquisition mean angles",
    "coefficient_source": "current libRadtran spectral LUT",
}
TARGET_RT_BY_MODE = {
    "sidecar": TARGET_RT,
    "physical": PHYSICAL_TARGET_RT,
    "physical_v2": PHYSICAL_V2_TARGET_RT,
}
CAMS_DIR = Path("/gws/ssde/j25b/nceo_ard/public/cams")
KG_M2_PER_DU = 2.1415e-5

CANONICAL_BANDS = ("coastal", "blue", "green", "red", "nir08", "swir16", "swir22")
GEE_L2A_BANDS = ("B1", "B2", "B3", "B4", "B8A", "B11", "B12")
CANONICAL_TO_S2_BAND = {
    "coastal": "B01",
    "blue": "B02",
    "green": "B03",
    "red": "B04",
    "nir08": "B8A",
    "swir16": "B11",
    "swir22": "B12",
}
L2A_AUX_BANDS = ("AOT", "WVP", "SCL")
LAND_SCL = (4, 5)
REFLECTANCE_SCALE = 1.0e-4
ATMOSPHERE_SCALE = 1.0e-3
HALF_SIZE_DEGREES = 0.06
DEFAULT_KEEP_FRACTION = 0.6
DEFAULT_MAX_SAMPLES_PER_SCENE = 1200


def _finite(value: Any) -> float | None:
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None


def _stable_seed(*parts: str) -> int:
    digest = hashlib.sha256("/".join(parts).encode("utf-8")).digest()
    return int.from_bytes(digest[:8], "little", signed=False)


def _selected_metadata(rows: list[dict[str, Any]], keep_fraction: float) -> list[dict[str, Any]]:
    usable = [row for row in rows if _finite(row.get("maiac")) is not None]
    usable.sort(key=lambda row: (float(row["maiac"]), str(row.get("l1c_id", ""))))
    count = max(1, int(round(len(usable) * float(keep_fraction)))) if usable else 0
    return usable[:count]


def _has_only_sparse_land_errors(errors: list[dict[str, str]]) -> bool:
    """Separate unavailable land support from infrastructure or RT failures."""
    return all(
        error.get("error_type") == "RuntimeError"
        and str(error.get("reason", "")).startswith("only ")
        and str(error.get("reason", "")).endswith(" paired clear-land pixels")
        for error in errors
    )


def _archive_scene_day_max(archive_path: Path) -> str | None:
    """Return the latest embedded acquisition date, or ``None`` if unreadable."""
    try:
        with np.load(archive_path, allow_pickle=False) as archive:
            scenes = json.loads(str(archive["scenes_json"].item()))
        days = [str(scene["day"]) for scene in scenes]
        for day in days:
            dt.date.fromisoformat(day)
    except (OSError, ValueError, KeyError, TypeError, json.JSONDecodeError):
        return None
    return max(days, default=None)


def _cached_result_is_compatible(
    existing: dict[str, Any],
    archive_path: Path,
    *,
    keep_fraction: float,
    max_samples_per_scene: int,
    include_terrain: bool,
    target_rt_mode: str,
    scene_day_max: str | None,
    deterrain_l2a: bool = False,
) -> bool:
    """Protect experiment changes from silently reusing stale pair archives.

    Archives created before the cutoff field was recorded may still be reused
    when their embedded dates prove that applying the requested cutoff would
    select exactly the same scenes. A cache generated with an earlier cutoff
    cannot satisfy a later or unbounded request because it may omit support.
    """
    if not math.isclose(
        float(existing.get("keep_fraction", DEFAULT_KEEP_FRACTION)),
        keep_fraction,
        rel_tol=0.0,
        abs_tol=1.0e-12,
    ):
        return False
    if int(existing.get("max_samples_per_scene", DEFAULT_MAX_SAMPLES_PER_SCENE)) != int(
        max_samples_per_scene
    ):
        return False
    if bool(existing.get("deterrain_l2a", False)) != bool(deterrain_l2a):
        return False

    status = existing.get("status")
    existing_cutoff = existing.get("scene_day_max")
    if status == "DATA_UNAVAILABLE":
        # Old unavailable audits did not record enough configuration to prove
        # that an RT/terrain change would remain unavailable, so rebuild them.
        if existing.get("target_rt_mode") != target_rt_mode:
            return False
        if bool(existing.get("include_terrain")) != include_terrain:
            return False
        if scene_day_max is None:
            return existing_cutoff is None
        return existing_cutoff is None or str(existing_cutoff) >= scene_day_max
    if status != "OK" or not archive_path.exists():
        return False

    expected_rt = TARGET_RT_BY_MODE[target_rt_mode]
    terrain_enabled = bool((existing.get("terrain_features") or {}).get("enabled", False))
    if existing.get("target_rt") != expected_rt or terrain_enabled != include_terrain:
        return False
    cached_mode = existing.get("target_rt_mode")
    if cached_mode is not None and cached_mode != target_rt_mode:
        return False

    if scene_day_max is None:
        return existing_cutoff is None
    if existing_cutoff is not None and str(existing_cutoff) < scene_day_max:
        return False
    archive_max = _archive_scene_day_max(archive_path)
    return archive_max is not None and archive_max <= scene_day_max


def _normalise_l2a_reflectance(values: np.ndarray) -> np.ndarray:
    result = np.asarray(values, dtype=np.float32) * np.float32(REFLECTANCE_SCALE)
    return np.where((result > 0.0) & (result < 1.2), result, np.nan).astype(np.float32)


def _normalise_aux(values: np.ndarray) -> np.ndarray:
    result = np.asarray(values, dtype=np.float32) * np.float32(ATMOSPHERE_SCALE)
    return np.where(np.isfinite(result) & (result > 0.0), result, np.nan).astype(np.float32)


def _site_elevation_km(row: dict[str, str], cache: dict[str, float]) -> float:
    key = f"{float(row['longitude']):.3f}_{float(row['latitude']):.3f}"
    value = _finite(cache.get(key))
    return float(value) if value is not None else 0.0


def _sentinel_satellite_id(spacecraft: Any) -> str:
    text = str(spacecraft).upper().replace("_", "").replace("-", "")
    for satellite_id in ("S2A", "S2B", "S2C"):
        if satellite_id in text or f"SENTINEL{satellite_id[1:]}" in text:
            return satellite_id
    raise ValueError(f"unsupported Sentinel-2 spacecraft {spacecraft!r}")


def _rt_field(
    values: np.ndarray | float,
    *,
    fallback: float,
    dims: tuple[str, str] = ("y", "x"),
) -> xr.DataArray:
    """Build a finite RT input field while retaining invalid pixels for later masking."""
    array = np.asarray(values, dtype=np.float32)
    finite = np.isfinite(array)
    fill = float(np.nanmedian(array[finite])) if np.any(finite) else float(fallback)
    if not np.isfinite(fill):
        raise ValueError("physical target RT field has no finite fallback")
    return xr.DataArray(np.where(finite, array, fill).astype(np.float32), dims=dims)


_CAMS_TCO3_CACHE: dict[tuple[str, float, float], float] = {}


def _cams_tco3_atm_cm(day: str, lat: float, lon: float) -> float:
    """Real CAMS total-ozone column (atm-cm) at the site and local overpass hour."""
    key = (day, round(lat, 2), round(lon, 2))
    cached = _CAMS_TCO3_CACHE.get(key)
    if cached is not None:
        return cached
    path = CAMS_DIR / f"{day}.nc"
    with xr.open_dataset(path) as data:
        field = data[["gtco3"]]
        if "forecast_reference_time" in field.dims:
            field = field.squeeze("forecast_reference_time", drop=True)
        hour = int(round(10.5 - lon / 15.0)) % 24
        if "forecast_period" in field.dims:
            selector: dict[str, Any] = {"forecast_period": float(hour)}
        else:
            selector = {"time": np.datetime64(f"{day}T{hour:02d}:00")}
        value = float(
            field.sel(latitude=lat, longitude=lon, **selector, method="nearest")["gtco3"].item()
        )
    tco3 = value / KG_M2_PER_DU / 1000.0
    if not 0.15 <= tco3 <= 0.55:
        raise ValueError(f"implausible CAMS ozone column {tco3:.3f} atm-cm for {day}")
    _CAMS_TCO3_CACHE[key] = tco3
    return tco3


def _physical_current_rt_surface(
    *,
    l1_toa: np.ndarray,
    l2a_tcwv: np.ndarray,
    terrain: TerrainFields,
    scene: dict[str, Any],
    spacecraft: Any,
    fallback_tcwv: float,
    rt_backend: Any,
    sensor_cache: dict[str, Any],
    tco3_atm_cm: float | None = None,
) -> np.ndarray:
    """Correct one historical L1C grid in the operational physical RT frame."""
    from siac.adapters.rsrf import load_sensor_config_with_rsrf
    from siac.runtime import AtmosphericState, GeometryAngles

    satellite_id = _sentinel_satellite_id(spacecraft)
    sensor = sensor_cache.get(satellite_id)
    if sensor is None:
        sensor = load_sensor_config_with_rsrf("MSI", satellite_id)
        sensor_cache[satellite_id] = sensor

    template = _rt_field(np.zeros(np.asarray(l1_toa).shape[1:], dtype=np.float32), fallback=0.0)
    aot = _finite(scene.get("maiac"))
    if aot is None or aot < 0.0:
        raise ValueError(f"invalid MAIAC AOD for physical target: {scene.get('maiac')!r}")
    tcwv = _rt_field(l2a_tcwv, fallback=fallback_tcwv)
    elevation = _rt_field(terrain.elevation_m * np.float32(1.0e-3), fallback=0.0)
    geometry = GeometryAngles.from_degrees(
        xr.full_like(template, float(scene["sza"])),
        xr.full_like(template, float(scene["saa"])),
        xr.full_like(template, float(scene["vza"])),
        xr.full_like(template, float(scene["vaa"])),
    )
    atmo = AtmosphericState(
        aot=xr.full_like(template, aot),
        tcwv=tcwv,
        tco3=xr.full_like(template, 0.30 if tco3_atm_cm is None else float(tco3_atm_cm)),
        aot_unc=xr.full_like(template, 0.10),
        tcwv_unc=xr.full_like(template, 0.50),
        tco3_unc=xr.full_like(template, 0.05),
        elevation=elevation,
    )
    sensor_bands = [sensor.get_band(CANONICAL_TO_S2_BAND[name]) for name in CANONICAL_BANDS]
    rt_backend.preload_scene_subset(geometry, atmo, sensor_bands)
    corrected = np.full(np.asarray(l1_toa).shape, np.nan, dtype=np.float32)
    for index, band in enumerate(sensor_bands):
        toa = xr.DataArray(np.asarray(l1_toa[index], dtype=np.float32), dims=template.dims)
        corrected[index] = np.asarray(
            rt_backend.compute_coefficients(geometry, atmo, band).apply_correction(toa).values,
            dtype=np.float32,
        )
    return corrected


def _resolve_l2a_asset(ee: Any, geom: Any, *, day: str, tile: str) -> dict[str, Any]:
    next_day = (dt.date.fromisoformat(day) + dt.timedelta(days=1)).isoformat()
    collection = (
        ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
        .filterBounds(geom)
        .filterDate(day, next_day)
        .filter(ee.Filter.eq("MGRS_TILE", tile))
    )
    count = int(collection.size().getInfo())
    if count == 0:
        raise RuntimeError(f"no L2A image for tile {tile} on {day}")
    image = ee.Image(collection.sort("system:time_start").first())
    info = ee.Dictionary(
        {
            "system_index": image.get("system:index"),
            "product_id": image.get("PRODUCT_ID"),
            "system_time_start": image.get("system:time_start"),
            "processing_baseline": image.get("PROCESSING_BASELINE"),
            "spacecraft": image.get("SPACECRAFT_NAME"),
            "n_images": count,
        }
    ).getInfo()
    if not info.get("system_index"):
        raise RuntimeError(f"L2A image has no system:index for tile {tile} on {day}")
    info["asset_id"] = f"COPERNICUS/S2_SR_HARMONIZED/{info['system_index']}"
    return info


def _load_pair_grids(
    ee: Any,
    *,
    bbox: tuple[float, float, float, float],
    sidecar_path: Path,
    scene: dict[str, Any],
    elevation_km: float,
    terrain: TerrainFields | None = None,
    target_rt_mode: str = "sidecar",
    rt_backend: Any | None = None,
    sensor_cache: dict[str, Any] | None = None,
    deterrain_l2a: bool = False,
    fraction_lut: Any | None = None,
) -> tuple[dict[str, np.ndarray], dict[str, Any]]:
    """Load an exact same-day L2A/current-RT pair on the common 60 m grid.

    The training archive stores a deterministic sample from this grid. The
    gallery uses this helper so it has the same assets, RT correction and
    clear-land test as the archived training pairs.
    """
    import bestpixel as bp
    from bestpixel._gee import get_patch, resolve_assets, utm_epsg_from_bbox, utm_grid
    from bestpixel.atmosphere import AtmoSidecar
    from bestpixel.l1c import SIDECAR_TO_GEE

    sidecar = AtmoSidecar.load(str(sidecar_path))
    scene_id = str(scene["l1c_id"])
    day = str(scene["day"])
    tile = str(scene["tile"])
    if scene_id not in sidecar.scenes:
        raise KeyError(f"{scene_id} is absent from {sidecar_path.name}")

    epsg = utm_epsg_from_bbox(bbox)
    grid = utm_grid(bbox, epsg, 60.0)
    geom = ee.Geometry.Rectangle(list(bbox))
    l1_assets = resolve_assets(ee, [scene_id], geom)
    if scene_id not in l1_assets:
        raise RuntimeError(f"could not resolve GEE L1C asset for {scene_id}")
    l2_info = _resolve_l2a_asset(ee, geom, day=day, tile=tile)

    l1_band_ids = [SIDECAR_TO_GEE[name] for name in CANONICAL_BANDS]
    l1_raw = get_patch(ee, l1_assets[scene_id][0], l1_band_ids, grid=grid)
    cloud_score = get_patch(ee, l1_assets[scene_id][1], ["cs"], grid=grid)[0]
    l2_raw = get_patch(
        ee,
        l2_info["asset_id"],
        [*GEE_L2A_BANDS, *L2A_AUX_BANDS],
        grid=grid,
    )

    atmospheric = sidecar.scenes[scene_id]
    l1_toa = np.asarray(l1_raw, dtype=np.float32) * np.float32(REFLECTANCE_SCALE)
    l2a_surface = _normalise_l2a_reflectance(l2_raw[: len(CANONICAL_BANDS)])
    l2a_aot = _normalise_aux(l2_raw[len(CANONICAL_BANDS)])
    l2a_tcwv = _normalise_aux(l2_raw[len(CANONICAL_BANDS) + 1])
    scl = np.asarray(l2_raw[len(CANONICAL_BANDS) + 2], dtype=np.int16)
    tco3_atm_cm: float | None = None
    if target_rt_mode in {"physical", "physical_v2"}:
        if terrain is None or rt_backend is None or sensor_cache is None:
            raise ValueError("physical target requires terrain, RT backend, and sensor cache")
        if target_rt_mode == "physical_v2":
            tco3_atm_cm = _cams_tco3_atm_cm(
                day, 0.5 * (bbox[1] + bbox[3]), 0.5 * (bbox[0] + bbox[2])
            )
        if deterrain_l2a:
            from tools.aeronet_validation.terrain_deshade import deterrain_l2a_surface

            if fraction_lut is None or sensor_cache is None:
                raise ValueError("deterrain_l2a requires the fraction LUT and sensor cache")
            l2a_surface, deterrain_provenance = deterrain_l2a_surface(
                l2a_surface,
                l2a_aot=l2a_aot,
                l2a_tcwv=l2a_tcwv,
                terrain=terrain,
                sza_deg=float(scene["sza"]),
                saa_deg=float(scene["saa"]),
                elevation_km=float(elevation_km),
                satellite_id=_sentinel_satellite_id(l2_info.get("spacecraft")),
                fraction_lut=fraction_lut,
                sensor_cache=sensor_cache,
                band_names=CANONICAL_BANDS,
                s2_band_map=CANONICAL_TO_S2_BAND,
                tco3_atm_cm=tco3_atm_cm if tco3_atm_cm is not None else 0.30,
            )
        siac_surface = _physical_current_rt_surface(
            l1_toa=l1_toa,
            l2a_tcwv=l2a_tcwv,
            terrain=terrain,
            scene=scene,
            spacecraft=l2_info.get("spacecraft"),
            fallback_tcwv=float(atmospheric.wvp),
            rt_backend=rt_backend,
            sensor_cache=sensor_cache,
            tco3_atm_cm=tco3_atm_cm,
        )
    elif target_rt_mode == "sidecar":
        triples = [
            atmospheric.coeffs(sidecar.band_index(name), atmospheric.wvp)
            for name in CANONICAL_BANDS
        ]
        siac_surface = np.asarray(
            bp.correct_toa(
                np.ascontiguousarray(l1_toa),
                [triple[0] for triple in triples],
                [triple[1] for triple in triples],
                [triple[2] for triple in triples],
            ),
            dtype=np.float32,
        )
    else:
        raise ValueError(f"unsupported target RT mode {target_rt_mode!r}")

    valid = (
        (np.asarray(cloud_score, dtype=np.float32) > 0.60)
        & np.isin(scl, LAND_SCL)
        & np.all(np.isfinite(l2a_surface), axis=0)
        & np.all(np.isfinite(siac_surface), axis=0)
        & np.isfinite(l2a_aot)
        & np.isfinite(l2a_tcwv)
        & np.all((l2a_surface > 0.001) & (l2a_surface < 0.8), axis=0)
        & np.all((siac_surface > 0.001) & (siac_surface < 0.8), axis=0)
    )
    valid = ndimage.binary_erosion(valid, structure=np.ones((3, 3), dtype=bool), iterations=1)

    grids = {
        "l1_toa": l1_toa,
        "l2a_surface": l2a_surface,
        "siac_surface": siac_surface,
        "l2a_aot": l2a_aot,
        "l2a_tcwv": l2a_tcwv,
        "cloud_score": np.asarray(cloud_score, dtype=np.float32),
        "scl": scl,
        "valid": valid,
    }
    raa = abs((float(scene["saa"]) - float(scene["vaa"]) + 180.0) % 360.0 - 180.0)
    metadata = {
        "scene_id": scene_id,
        "day": day,
        "tile": tile,
        "maiac_aot": float(scene["maiac"]),
        "maiac_tcwv_cm": float(scene["wvp"]),
        "sza_deg": float(scene["sza"]),
        "saa_deg": float(scene["saa"]),
        "vza_deg": float(scene["vza"]),
        "vaa_deg": float(scene["vaa"]),
        "raa_deg": raa,
        "elevation_km": float(elevation_km),
        "month": int(day[5:7]),
        "l2a": {key: value for key, value in l2_info.items() if key != "asset_id"},
        "l1c_asset": l1_assets[scene_id][0],
        "cloud_score_asset": l1_assets[scene_id][1],
        "sidecar": str(sidecar_path),
        "target_rt_mode": target_rt_mode,
        "tco3_atm_cm": tco3_atm_cm,
        "l2a_deterrain": deterrain_provenance if deterrain_l2a else {"applied": False},
        "grid": {
            "epsg": int(grid["epsg"]),
            "resolution": float(grid["res"]),
            "width": int(grid["W"]),
            "height": int(grid["H"]),
        },
    }
    return grids, metadata


def _pair_scene(
    ee: Any,
    *,
    bbox: tuple[float, float, float, float],
    sidecar_path: Path,
    scene: dict[str, Any],
    elevation_km: float,
    max_samples: int,
    matchup_id: str,
    terrain: TerrainFields | None = None,
    include_terrain_features: bool = False,
    target_rt_mode: str = "sidecar",
    rt_backend: Any | None = None,
    sensor_cache: dict[str, Any] | None = None,
    deterrain_l2a: bool = False,
    fraction_lut: Any | None = None,
) -> tuple[dict[str, np.ndarray], dict[str, Any]]:
    grids, metadata = _load_pair_grids(
        ee,
        bbox=bbox,
        sidecar_path=sidecar_path,
        scene=scene,
        elevation_km=elevation_km,
        terrain=terrain,
        target_rt_mode=target_rt_mode,
        rt_backend=rt_backend,
        sensor_cache=sensor_cache,
        deterrain_l2a=deterrain_l2a,
        fraction_lut=fraction_lut,
    )
    scene_id = str(scene["l1c_id"])
    day = str(scene["day"])
    valid = np.asarray(grids["valid"], dtype=bool)
    if terrain is not None and (
        include_terrain_features or target_rt_mode in {"physical", "physical_v2"}
    ):
        if terrain.valid.shape != valid.shape:
            raise ValueError("terrain grid does not match the exact pair grid")
        valid &= terrain.valid
    indices = np.flatnonzero(valid)
    if indices.size < 100:
        raise RuntimeError(f"only {indices.size} paired clear-land pixels")
    rng = np.random.default_rng(_stable_seed(matchup_id, day, scene_id))
    if indices.size > max_samples:
        indices = np.sort(rng.choice(indices, size=max_samples, replace=False))

    flat_l2a = np.asarray(grids["l2a_surface"]).reshape(len(CANONICAL_BANDS), -1).T
    flat_siac = np.asarray(grids["siac_surface"]).reshape(len(CANONICAL_BANDS), -1).T
    samples = {
        "l2a": flat_l2a[indices].astype(np.float32),
        "siac": flat_siac[indices].astype(np.float32),
        "l2a_aot": np.asarray(grids["l2a_aot"]).ravel()[indices].astype(np.float32),
        "l2a_tcwv": np.asarray(grids["l2a_tcwv"]).ravel()[indices].astype(np.float32),
    }
    if terrain is not None and include_terrain_features:
        incidence = local_solar_incidence(
            terrain,
            sza_deg=metadata["sza_deg"],
            saa_deg=metadata["saa_deg"],
        )
        samples.update(
            {
                "terrain_elevation_km": (
                    terrain.elevation_m.ravel()[indices] * np.float32(1.0e-3)
                ).astype(np.float32),
                "terrain_slope_deg": terrain.slope_deg.ravel()[indices].astype(np.float32),
                "terrain_incidence_cos": incidence.ravel()[indices].astype(np.float32),
            }
        )
        metadata["terrain_source"] = terrain.source
    metadata["sample_count"] = int(indices.size)
    metadata["l2a_aot_median"] = float(np.nanmedian(samples["l2a_aot"]))
    metadata["l2a_tcwv_median_cm"] = float(np.nanmedian(samples["l2a_tcwv"]))
    return samples, metadata


def build_one(
    matchup_id: str,
    rows: dict[str, dict[str, str]],
    *,
    output_dir: Path,
    keep_fraction: float,
    max_samples_per_scene: int,
    include_terrain: bool,
    target_rt_mode: str,
    scene_day_max: str | None,
    force: bool,
    deterrain_l2a: bool = False,
) -> dict[str, Any]:
    from bestpixel._gee import init_ee, utm_epsg_from_bbox, utm_grid

    if target_rt_mode not in TARGET_RT_BY_MODE:
        raise ValueError(f"unsupported target RT mode {target_rt_mode!r}")
    build_config = {
        "keep_fraction": float(keep_fraction),
        "max_samples_per_scene": int(max_samples_per_scene),
        "include_terrain": bool(include_terrain),
        "target_rt_mode": target_rt_mode,
        "scene_day_max": scene_day_max,
        "deterrain_l2a": bool(deterrain_l2a),
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    archive_path = output_dir / f"{matchup_id}.npz"
    audit_path = output_dir / f"{matchup_id}.json"
    if audit_path.exists() and not force:
        existing = json.loads(audit_path.read_text(encoding="utf-8"))
        if _cached_result_is_compatible(
            existing,
            archive_path,
            keep_fraction=keep_fraction,
            max_samples_per_scene=max_samples_per_scene,
            include_terrain=include_terrain,
            target_rt_mode=target_rt_mode,
            scene_day_max=scene_day_max,
            deterrain_l2a=deterrain_l2a,
        ):
            # Record inferred defaults/cutoff on compatible legacy archives so
            # future runs can validate the cache without inference.
            if any(existing.get(key) != value for key, value in build_config.items()):
                existing.update(build_config)
                audit_path.write_text(json.dumps(existing, indent=2) + "\n", encoding="utf-8")
            return existing
        # Migrate audits produced before sparse land support had its own
        # terminal state. This is deterministic and avoids repeating GEE/RT
        # work merely to reclassify a physically unsupported site.
        if (
            existing.get("status") == "FAILED"
            and int(existing.get("successful_scenes", 0)) < 2
            and _has_only_sparse_land_errors(existing.get("errors") or [])
        ):
            existing.update(
                {
                    "status": "DATA_UNAVAILABLE",
                    "uses_aeronet": False,
                    **build_config,
                    "reason": "insufficient paired clear-land support for harmonizer training",
                }
            )
            audit_path.write_text(json.dumps(existing, indent=2) + "\n", encoding="utf-8")
            return existing

    # An incompatible or forced rebuild must not leave a stale successful
    # archive beside a failed/unavailable replacement audit.
    archive_path.unlink(missing_ok=True)

    row = rows[matchup_id]
    lon = float(row["longitude"])
    lat = float(row["latitude"])
    bbox = (
        lon - HALF_SIZE_DEGREES,
        lat - HALF_SIZE_DEGREES,
        lon + HALF_SIZE_DEGREES,
        lat + HALF_SIZE_DEGREES,
    )
    elevation_cache = json.loads(ELEVATION_CACHE.read_text(encoding="utf-8"))
    elevation_km = _site_elevation_km(row, elevation_cache)
    history_path = L2A_HISTORY / f"{matchup_id}.npz"
    with np.load(history_path, allow_pickle=False) as history:
        windows = [str(value) for value in history["realizations"]]

    rt_backend = None
    sensor_cache: dict[str, Any] | None = None
    if target_rt_mode in {"physical", "physical_v2"}:
        from siac.algorithms.rt.lut.backend import ZarrLUTBackend

        rt_backend = ZarrLUTBackend(TARGET_LUT)
        sensor_cache = {}
    fraction_lut = None
    if deterrain_l2a:
        from tools.aeronet_validation.terrain_deshade import DirectFractionLUT

        if target_rt_mode not in {"physical", "physical_v2"}:
            raise ValueError("--deterrain-l2a requires a physical target mode")
        fraction_lut = DirectFractionLUT(TARGET_LUT)

    ee = init_ee()
    terrain = None
    if include_terrain or target_rt_mode in {"physical", "physical_v2"}:
        terrain = fetch_glo30_terrain(ee, utm_grid(bbox, utm_epsg_from_bbox(bbox), 60.0))
    all_samples: dict[str, list[np.ndarray]] = {
        "l2a": [],
        "siac": [],
        "l2a_aot": [],
        "l2a_tcwv": [],
    }
    if include_terrain:
        all_samples.update(
            {
                "terrain_elevation_km": [],
                "terrain_slope_deg": [],
                "terrain_incidence_cos": [],
            }
        )
    scene_indices: list[np.ndarray] = []
    scenes: list[dict[str, Any]] = []
    errors: list[dict[str, str]] = []
    attempted = 0
    started = time.monotonic()
    for window in windows:
        metadata_path = SIDECAR_ROOT / matchup_id / f"{window}_meta.json"
        sidecar_path = SIDECAR_ROOT / matchup_id / f"{window}_lut_sidecar.json"
        if not metadata_path.exists() or not sidecar_path.exists():
            errors.append({"window": window, "reason": "missing metadata or LUT sidecar"})
            continue
        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
        selected = _selected_metadata(metadata.get("selected", []), keep_fraction)
        if scene_day_max is not None:
            selected = [scene for scene in selected if str(scene.get("day", "")) <= scene_day_max]
        for scene in selected:
            attempted += 1
            try:
                samples, scene_metadata = _pair_scene(
                    ee,
                    bbox=bbox,
                    sidecar_path=sidecar_path,
                    scene=scene,
                    elevation_km=elevation_km,
                    max_samples=max_samples_per_scene,
                    matchup_id=matchup_id,
                    terrain=terrain,
                    include_terrain_features=include_terrain,
                    target_rt_mode=target_rt_mode,
                    rt_backend=rt_backend,
                    sensor_cache=sensor_cache,
                    deterrain_l2a=deterrain_l2a,
                    fraction_lut=fraction_lut,
                )
            except Exception as exc:  # noqa: BLE001 - retain partial acquisition progress
                errors.append(
                    {
                        "window": window,
                        "scene_id": str(scene.get("l1c_id", "")),
                        "day": str(scene.get("day", "")),
                        "error_type": type(exc).__name__,
                        "reason": str(exc)[:500],
                    }
                )
                continue
            scene_index = len(scenes)
            scene_metadata["window"] = window
            scenes.append(scene_metadata)
            count = samples["l2a"].shape[0]
            scene_indices.append(np.full(count, scene_index, dtype=np.int16))
            for name in all_samples:
                all_samples[name].append(samples[name])
            print(
                f"PAIR {matchup_id} {scene_metadata['day']} n={count} "
                f"maiac={scene_metadata['maiac_aot']:.3f} "
                f"l2a={scene_metadata['l2a_aot_median']:.3f}",
                flush=True,
            )

    sample_count = sum(values.shape[0] for values in all_samples["l2a"])
    # Some persistently cloudy sites only expose two valid acquisitions. Two
    # independent days still provide enough data for training when the pixel
    # support is retained and the rejected acquisitions remain in the audit.
    if attempted == 0 and not errors:
        status = {
            "matchup_id": matchup_id,
            "status": "DATA_UNAVAILABLE",
            "uses_aeronet": False,
            **build_config,
            "reason": "no historical L1C sidecar scenes with MAIAC-conditioned metadata",
            "attempted_scenes": 0,
            "successful_scenes": 0,
            "errors": [],
            "runtime_s": time.monotonic() - started,
        }
        audit_path.write_text(json.dumps(status, indent=2) + "\n", encoding="utf-8")
        return status
    if (len(scenes) < 2 or sample_count < 1000) and _has_only_sparse_land_errors(errors):
        status = {
            "matchup_id": matchup_id,
            "status": "DATA_UNAVAILABLE",
            "uses_aeronet": False,
            **build_config,
            "reason": "insufficient paired clear-land support for harmonizer training",
            "attempted_scenes": attempted,
            "successful_scenes": len(scenes),
            "sample_count": sample_count,
            "errors": errors,
            "runtime_s": time.monotonic() - started,
        }
        audit_path.write_text(json.dumps(status, indent=2) + "\n", encoding="utf-8")
        return status
    if len(scenes) < 2 or sample_count < 1000:
        status = {
            "matchup_id": matchup_id,
            "status": "FAILED",
            **build_config,
            "attempted_scenes": attempted,
            "successful_scenes": len(scenes),
            "errors": errors,
            "runtime_s": time.monotonic() - started,
        }
        audit_path.write_text(json.dumps(status, indent=2) + "\n", encoding="utf-8")
        return status

    payload = {name: np.concatenate(values, axis=0) for name, values in all_samples.items()}
    payload["scene_index"] = np.concatenate(scene_indices)
    payload["band_names"] = np.asarray(CANONICAL_BANDS)
    payload["scenes_json"] = np.asarray(json.dumps(scenes))
    payload["schema_version"] = np.asarray(2 if include_terrain else 1, dtype=np.int16)
    temporary = archive_path.with_suffix(".tmp.npz")
    np.savez_compressed(temporary, **payload)
    temporary.replace(archive_path)
    status = {
        "matchup_id": matchup_id,
        "status": "OK",
        "uses_aeronet": False,
        **build_config,
        "target": {
            "sidecar": "same-day L1C corrected with MAIAC AOD and current libRadtran LUT",
            "physical": "same-day L1C corrected with current physical libRadtran LUT state",
            "physical_v2": (
                "same-day L1C corrected with real CAMS ozone and the current "
                "physical libRadtran LUT state"
            ),
        }[target_rt_mode],
        "target_rt": TARGET_RT_BY_MODE[target_rt_mode],
        "terrain_features": {
            "enabled": include_terrain,
            "source": TERRAIN_SOURCE if include_terrain else None,
            "fields": ["elevation_km", "slope_deg", "solar_incidence_cos"]
            if include_terrain
            else [],
        },
        "attempted_scenes": attempted,
        "successful_scenes": len(scenes),
        "sample_count": int(payload["l2a"].shape[0]),
        "errors": errors,
        "runtime_s": time.monotonic() - started,
        "archive": str(archive_path),
    }
    audit_path.write_text(json.dumps(status, indent=2) + "\n", encoding="utf-8")
    return status


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("matchup_id", nargs="+")
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--keep-fraction", type=float, default=DEFAULT_KEEP_FRACTION)
    parser.add_argument("--max-samples-per-scene", type=int, default=DEFAULT_MAX_SAMPLES_PER_SCENE)
    parser.add_argument(
        "--target-rt-mode",
        choices=("sidecar", "physical", "physical_v2"),
        default="sidecar",
        help="sidecar uses the legacy saved coefficients; physical recomputes current-LUT "
        "coefficients; physical_v2 adds the MAIAC 1-km field and real CAMS ozone",
    )
    parser.add_argument(
        "--include-terrain",
        action="store_true",
        help="attach GLO-30 elevation, slope and local solar incidence per sample",
    )
    parser.add_argument(
        "--scene-day-max",
        default=None,
        help="optional latest historical acquisition day (YYYY-MM-DD)",
    )
    parser.add_argument(
        "--deterrain-l2a",
        action="store_true",
        help="analytically invert Sen2Cor's rugged-terrain correction before pairing",
    )
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if not 0.0 < args.keep_fraction <= 1.0:
        raise SystemExit("--keep-fraction must be in (0, 1]")
    if args.max_samples_per_scene < 100:
        raise SystemExit("--max-samples-per-scene must be at least 100")
    if args.scene_day_max is not None:
        dt.date.fromisoformat(args.scene_day_max)
    rows = {row["matchup_id"]: row for row in csv.DictReader(MATCHUPS.open())}
    exit_code = 0
    for matchup_id in args.matchup_id:
        if matchup_id not in rows:
            raise SystemExit(f"unknown matchup_id {matchup_id}")
        result = build_one(
            matchup_id,
            rows,
            output_dir=args.output_dir,
            keep_fraction=args.keep_fraction,
            max_samples_per_scene=args.max_samples_per_scene,
            include_terrain=args.include_terrain,
            target_rt_mode=args.target_rt_mode,
            scene_day_max=args.scene_day_max,
            force=args.force,
            deterrain_l2a=args.deterrain_l2a,
        )
        print("PAIR_DONE " + json.dumps(result), flush=True)
        if result["status"] not in {"OK", "DATA_UNAVAILABLE"}:
            exit_code = 1
    raise SystemExit(exit_code)


if __name__ == "__main__":
    main()
