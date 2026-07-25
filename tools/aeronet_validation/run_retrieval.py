"""Stage 3: SIAC aerosol retrieval per matchup x surface-prior approach.

Each task processes one Sentinel-2 acquisition over a small AOI centered on
the AERONET site with one surface-prior approach, keeping every other config
knob identical so the surface prior is the controlled variable. Results are
written under ``runs/<approach>/<matchup_id>/`` (config snapshot, extracted
site statistics, small AOT/TCWV rasters) and completed runs are skipped, so
the stage is resumable and SLURM-array friendly.
"""

from __future__ import annotations

import logging
import os
import time
import traceback
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Any

import pandas as pd
from tools.aeronet_validation.common import (
    APPROACHES,
    DEFAULT_DEM,
    DEFAULT_WATER_MASK,
    ExperimentPaths,
    map_legacy_earthdata_env,
    read_json,
    require_slurm_execution,
    resolve_cams_data_path,
    resolve_lut_path,
    write_json,
)

if TYPE_CHECKING:
    import argparse

logger = logging.getLogger("aeronet_validation.run")
LOCAL_EDOWN_RUNTIME = Path(__file__).resolve().parents[1] / "edown_runtime"


@dataclass(frozen=True)
class RunSettings:
    aoi_half_width_deg: float
    extract_radius_m: float
    aerosol_resolution_m: float
    s2_backend: str
    max_workers: int
    save_rasters: bool
    overwrite: bool


def _resolve_edown_executable() -> Path:
    configured = os.environ.get("SIAC_EDOWN_EXECUTABLE") or os.environ.get("EDOWN_EXECUTABLE")
    candidates = [Path(configured).expanduser()] if configured else []
    candidates.append(LOCAL_EDOWN_RUNTIME)

    from siac.adapters.brdf.mcd43_gee import DEFAULT_EDOWN_EXECUTABLE

    candidates.append(Path(DEFAULT_EDOWN_EXECUTABLE).expanduser())
    for candidate in candidates:
        if candidate.is_file():
            return candidate
    searched = ", ".join(str(path) for path in candidates)
    raise FileNotFoundError(f"edown executable not found; searched: {searched}")


def build_config_payload(
    paths: ExperimentPaths, approach: str, settings: RunSettings
) -> dict[str, Any]:
    """Nested SIACConfig payload (TOML shape); approach overlay merged last."""
    cache = paths.cache_dir
    payload: dict[str, Any] = {
        "sensor": "s2",
        "paths": {
            "dem": DEFAULT_DEM,
            "water_mask": DEFAULT_WATER_MASK,
            "lut_path": resolve_lut_path(),
            "cache_root": str(cache),
            "caches": {"s2": str(cache / "s2")},
        },
        "providers": {
            "atmo": {
                "kind": "cams",
                "data_path": resolve_cams_data_path(cache),
                "cache_dir": str(cache / "cams-cache"),
                "download_missing": True,
            },
            # BRDF prior sourced from GEE (MODIS/061/MCD43A1) via edown, so the
            # experiment needs no NASA Earthdata credentials. Only the source
            # changes; the three surface-prior approaches compose identically.
            "brdf": {
                "kind": "mcd43_gee",
                "cache_dir": str(cache / "mcd43_gee"),
                "data_path": str(_resolve_edown_executable()),
                "temporal_window": 16,
                "use_cache": True,
            },
            "s2": {
                "backend": settings.s2_backend,
                "cache_dir": str(cache / "s2"),
                "processing_level": "L1C",
                "max_cloud_cover": 100.0,
            },
        },
        "algorithms": {
            "rt": {"backend": "lut"},
            "solver": {"aerosol_resolution": settings.aerosol_resolution_m},
        },
        "runtime": {
            "log_level": "INFO",
            "execution": {"backend": "thread", "max_workers": settings.max_workers},
        },
    }
    overlay = APPROACHES[approach]
    return _deep_merge(payload, overlay)


def _deep_merge(base: dict[str, Any], updates: dict[str, Any]) -> dict[str, Any]:
    merged = dict(base)
    for key, value in updates.items():
        if isinstance(value, dict) and isinstance(merged.get(key), dict):
            merged[key] = _deep_merge(merged[key], value)
        else:
            merged[key] = value
    return merged


def _write_config_snapshot(payload: dict[str, Any], path: Path) -> None:
    import tomli_w

    path.write_text(tomli_w.dumps(_stringify_paths(payload)))


def _stringify_paths(value: Any) -> Any:
    if isinstance(value, dict):
        return {key: _stringify_paths(item) for key, item in value.items() if item is not None}
    if isinstance(value, Path):
        return str(value)
    return value


def extract_site_statistics(
    result: Any, longitude: float, latitude: float, radius_m: float
) -> dict[str, Any]:
    """Site-window statistics from a CorrectionResult on the scene grid."""
    import numpy as np
    from pyproj import Transformer

    statistics: dict[str, Any] = {}
    aot = result.aot
    crs = aot.rio.crs if hasattr(aot, "rio") else None
    if crs is None:
        raise ValueError("Retrieved AOT raster has no CRS; cannot extract site value.")
    transformer = Transformer.from_crs("EPSG:4326", crs, always_xy=True)
    x_site, y_site = transformer.transform(longitude, latitude)
    statistics["site_x"] = float(x_site)
    statistics["site_y"] = float(y_site)

    for name, field in (("aot", aot), ("tcwv", result.tcwv)):
        if field is None:
            continue
        values = field.values
        x_coords = field.coords["x"].values
        y_coords = field.coords["y"].values
        in_window = (np.abs(x_coords[None, :] - x_site) <= radius_m) & (
            np.abs(y_coords[:, None] - y_site) <= radius_m
        )
        window_values = values[in_window]
        finite = window_values[np.isfinite(window_values)]
        nearest = field.sel(x=x_site, y=y_site, method="nearest")
        nearest_value = float(nearest.values)
        statistics[f"{name}_nearest"] = nearest_value if np.isfinite(nearest_value) else None
        statistics[f"{name}_window_mean"] = float(finite.mean()) if finite.size else None
        statistics[f"{name}_window_median"] = float(np.median(finite)) if finite.size else None
        statistics[f"{name}_window_std"] = float(finite.std(ddof=0)) if finite.size else None
        statistics[f"{name}_window_n_valid"] = int(finite.size)
        statistics[f"{name}_window_n_total"] = int(window_values.size)

    # aoi_cloud_fraction must be ACTUAL cloud (raw OmniCloudMask from M1), not
    # result.cloud_mask, which ORs in per-band invalid-BOA pixels — a clear
    # scene with a poorly-matched surface prior would otherwise read as 100%
    # "cloud" and be wrongly dropped by the compare-stage filter. The folded
    # mask is recorded separately as a retrieval-quality metric.
    raw_cloud = result.cloud_mask_m1 if result.cloud_mask_m1 is not None else result.cloud_mask
    if raw_cloud is not None:
        statistics["aoi_cloud_fraction"] = float(raw_cloud.values.astype(bool).mean())
    if result.cloud_mask is not None:
        statistics["aoi_invalid_fraction"] = float(result.cloud_mask.values.astype(bool).mean())
    return statistics


def run_single(
    paths: ExperimentPaths,
    matchup: dict[str, Any],
    approach: str,
    settings: RunSettings,
) -> dict[str, Any]:
    """Run one matchup x approach task; returns the result record."""
    map_legacy_earthdata_env()
    matchup_id = str(matchup["matchup_id"])
    run_dir = paths.run_dir(approach, matchup_id)
    result_path = run_dir / "result.json"
    if result_path.exists() and not settings.overwrite:
        existing_status = read_json(result_path).get("status")
        if existing_status in {"ok", "no_valid_observation"}:
            logger.info("Skipping %s/%s: terminal result exists", approach, matchup_id)
            return {"matchup_id": matchup_id, "approach": approach, "status": "cached"}
        logger.info(
            "Re-running %s/%s: previous result status was %r",
            approach,
            matchup_id,
            existing_status,
        )
    run_dir.mkdir(parents=True, exist_ok=True)

    longitude = float(matchup["longitude"])
    latitude = float(matchup["latitude"])
    half = settings.aoi_half_width_deg
    aoi_bounds = (longitude - half, latitude - half, longitude + half, latitude + half)

    payload = build_config_payload(paths, approach, settings)
    _write_config_snapshot(payload, run_dir / "config.toml")

    record: dict[str, Any] = {
        "matchup_id": matchup_id,
        "approach": approach,
        "site": matchup["site"],
        "product_id": matchup["product_id"],
        "sensing_time_utc": matchup["sensing_time_utc"],
        "longitude": longitude,
        "latitude": latitude,
        "aoi_bounds_wgs84": list(aoi_bounds),
        "aerosol_resolution_m": settings.aerosol_resolution_m,
        "extract_radius_m": settings.extract_radius_m,
    }
    started = time.time()
    try:
        from siac.api import siac_process_s2
        from siac.config import SIACConfig

        config = SIACConfig.model_validate(payload)
        result = siac_process_s2(config, str(matchup["product_id"]), aoi=aoi_bounds)
        record.update(
            extract_site_statistics(result, longitude, latitude, settings.extract_radius_m)
        )
        record["status"] = (
            "ok" if int(record.get("aot_window_n_valid", 0)) > 0 else "no_valid_observation"
        )
        if settings.save_rasters:
            _save_rasters(result, run_dir / "fields.nc")
    except Exception as error:  # noqa: BLE001 - per-task failure is recorded, not fatal
        record["status"] = "failed"
        record["error"] = f"{type(error).__name__}: {error}"
        (run_dir / "traceback.txt").write_text(traceback.format_exc())
        logger.error("Run failed for %s/%s: %s", approach, matchup_id, record["error"])
    record["wall_time_s"] = round(time.time() - started, 1)
    write_json(result_path, record)
    return record


def _save_rasters(result: Any, path: Path) -> None:
    import xarray as xr

    fields = {"aot": result.aot, "tcwv": result.tcwv}
    if result.cloud_mask is not None:
        fields["cloud_mask"] = result.cloud_mask.astype("uint8")
    dataset = xr.Dataset({k: v for k, v in fields.items() if v is not None})
    encoding = {name: {"zlib": True, "complevel": 4} for name in dataset.data_vars}
    dataset.to_netcdf(path, encoding=encoding)


def build_manifest(
    paths: ExperimentPaths,
    approaches: list[str],
    *,
    limit: int | None = None,
    per_site_per_month: int | None = None,
) -> pd.DataFrame:
    matchups = pd.read_csv(paths.matchups_file)
    if per_site_per_month is not None:
        # Seasonally balanced subsample: up to N matchups per site per
        # calendar month, preferring the least cloudy scenes. One S2 SAFE is
        # ~600 MB, so sampling controls download volume as well as compute.
        matchups["month"] = pd.to_datetime(matchups["sensing_time_utc"]).dt.month
        matchups = (
            matchups.sort_values("scene_cloud_cover")
            .groupby(["site", "month"], as_index=False)
            .head(per_site_per_month)
            .sort_values(["site", "sensing_time_utc"])
            .drop(columns="month")
            .reset_index(drop=True)
        )
    if limit is not None:
        matchups = matchups.head(limit)
    tasks = [
        {"task_index": 0, "matchup_id": row.matchup_id, "approach": approach}
        for approach in approaches
        for row in matchups.itertuples(index=False)
    ]
    manifest = pd.DataFrame(tasks)
    manifest["task_index"] = range(len(manifest))
    paths.runs_dir.mkdir(parents=True, exist_ok=True)
    manifest.to_csv(paths.manifest_file, index=False)
    logger.info(
        "Manifest: %d tasks (%d matchups x %d approaches) at %s",
        len(manifest),
        len(matchups),
        len(approaches),
        paths.manifest_file,
    )
    return manifest


def _load_matchups_by_id(paths: ExperimentPaths) -> dict[str, dict[str, Any]]:
    matchups = pd.read_csv(paths.matchups_file)
    return {str(row["matchup_id"]): row.to_dict() for _, row in matchups.iterrows()}


def add_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--approaches",
        nargs="*",
        choices=sorted(APPROACHES),
        default=sorted(APPROACHES),
        help="Surface-prior approaches to run.",
    )
    parser.add_argument(
        "--matchup-id", default=None, help="Run a single matchup id (all selected approaches)."
    )
    parser.add_argument(
        "--task-index",
        type=int,
        default=None,
        help="Run one task from the manifest (for SLURM arrays).",
    )
    parser.add_argument("--limit", type=int, default=None, help="Cap number of matchups.")
    parser.add_argument(
        "--aoi-half-width-deg",
        type=float,
        default=0.05,
        help="Half-width of the retrieval AOI around the site, degrees (~5.5 km).",
    )
    parser.add_argument(
        "--extract-radius-m",
        type=float,
        default=1500.0,
        help="Half-width of the site extraction window on the aerosol grid.",
    )
    parser.add_argument("--aerosol-resolution", type=float, default=120.0)
    parser.add_argument("--s2-backend", choices=("gcs", "cdse"), default="gcs")
    parser.add_argument("--max-workers", type=int, default=4)
    parser.add_argument("--no-rasters", action="store_true", help="Skip writing fields.nc.")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument(
        "--allow-login",
        action="store_true",
        help="Run locally on an interactive node (not recommended).",
    )


def _preflight_gee_source() -> None:
    """Fail fast if the GEE BRDF source (edown + Earth Engine creds) is missing.

    The BRDF prior is downloaded from MODIS/061/MCD43A1 via edown, which needs
    its executable on PATH (or at the configured path) and Earth Engine
    credentials under ~/.config/earthengine. Checking once per job turns an
    obscure mid-pipeline failure into an immediate, actionable error.
    """
    try:
        edown = _resolve_edown_executable()
    except FileNotFoundError as exc:
        raise SystemExit(str(exc)) from exc
    creds = Path.home() / ".config" / "earthengine" / "credentials"
    if not creds.exists():
        raise SystemExit(
            "No Earth Engine credentials at ~/.config/earthengine/credentials. "
            "Run `earthengine authenticate` (or `edown`'s auth) before submitting."
        )
    logger.info("GEE BRDF pre-flight OK: edown=%s, EE creds present", edown)


def run(args: argparse.Namespace, paths: ExperimentPaths) -> None:
    require_slurm_execution(
        "AERONET retrieval stage",
        allow_login=args.allow_login,
        suggestion=(
            "Use `make-slurm` from `tools.aeronet_validation.cli` and submit via "
            "`sbatch` (for example: `sbatch <data-root>/slurm/submit_runs.sbatch`)."
        ),
    )
    paths.ensure()
    _preflight_gee_source()
    settings = RunSettings(
        aoi_half_width_deg=args.aoi_half_width_deg,
        extract_radius_m=args.extract_radius_m,
        aerosol_resolution_m=args.aerosol_resolution,
        s2_backend=args.s2_backend,
        max_workers=args.max_workers,
        save_rasters=not args.no_rasters,
        overwrite=args.overwrite,
    )
    matchups_by_id = _load_matchups_by_id(paths)

    if args.task_index is not None:
        manifest = pd.read_csv(paths.manifest_file)
        task = manifest[manifest["task_index"] == args.task_index]
        if task.empty:
            raise SystemExit(f"task_index {args.task_index} not in {paths.manifest_file}")
        row = task.iloc[0]
        run_single(paths, matchups_by_id[str(row["matchup_id"])], str(row["approach"]), settings)
        return

    if args.matchup_id is not None:
        selected = [args.matchup_id]
    else:
        selected = list(matchups_by_id)
        if args.limit is not None:
            selected = selected[: args.limit]

    total = len(selected) * len(args.approaches)
    done = 0
    for matchup_id in selected:
        for approach in args.approaches:
            record = run_single(paths, matchups_by_id[matchup_id], approach, settings)
            done += 1
            logger.info("[%d/%d] %s/%s -> %s", done, total, approach, matchup_id, record["status"])
