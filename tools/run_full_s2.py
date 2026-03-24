#!/usr/bin/env python3
"""Run an end-to-end SIAC processing workflow with Sentinel-2 input."""

from __future__ import annotations

import argparse
import json
import logging
import re
from datetime import datetime
from pathlib import Path
from typing import Any
from urllib.parse import urlparse

import numpy as np
import xarray as xr

from siac.adapters.satellite.sentinel2 import Sentinel2Preprocessor
from siac.api import resolve_s2_input, siac_process
from siac.config import DEFAULT_LUT_URL, SIACConfig
from siac.domain import AOI, SensorConfig
from siac.geo.reprojection import reproject_match
from siac.runtime import GeometryAngles, ObservationBundle
from siac.storage.product_writers import write_dataset
from siac.storage.raster_writers import write_netcdf
from siac.storage.stac import write_stac_item

logger = logging.getLogger("run_full_s2")
_L2A_FOLDER_PATTERN = re.compile(r"^S2[A-Z]_MSIL2A_.*T.*$")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Fetch/resolve Sentinel-2 input and run full SIAC processing "
            "(M1->M6), writing BOA and atmospheric outputs."
        )
    )
    parser.add_argument(
        "--query",
        required=True,
        help=(
            "Sentinel-2 query: local SAFE path, product ID, or tile+date "
            "string like T31UDQ_20210801."
        ),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="Directory where BOA and atmospheric outputs will be written.",
    )
    parser.add_argument(
        "--backend",
        choices=("gcs", "cdse", "local"),
        default="gcs",
        help="Sentinel-2 data backend used when --query is not a local existing path.",
    )
    parser.add_argument(
        "--cache-dir",
        type=Path,
        default=Path.home() / ".cache" / "siac" / "s2",
        help="Cache directory for downloaded Sentinel-2 SAFE products.",
    )
    parser.add_argument(
        "--max-cloud-cover",
        type=float,
        default=80.0,
        help="Max cloud cover filter used during Sentinel-2 search.",
    )
    parser.add_argument(
        "--processing-level",
        choices=("L1C", "L2A"),
        default="L1C",
        help="Target Sentinel-2 processing level for search/fetch.",
    )
    parser.add_argument(
        "--aoi-bbox",
        nargs=4,
        type=float,
        metavar=("XMIN", "YMIN", "XMAX", "YMAX"),
        help=(
            "Optional AOI bounds. If provided, TOA/geometry/cloud are clipped "
            "before running M2-M6."
        ),
    )
    parser.add_argument(
        "--aoi-crs",
        default="EPSG:4326",
        help="CRS for --aoi-bbox (default: EPSG:4326).",
    )
    parser.add_argument(
        "--cams-dir",
        default=None,
        help=(
            "Directory, file path, HTTP base/file URL, or s3:// URL for CAMS data. "
            "If omitted, uses <cache-dir>/cams and falls back to provider defaults when empty."
        ),
    )
    parser.add_argument(
        "--lut-path",
        default=DEFAULT_LUT_URL,
        help="Path/URL to LUT store used for radiative transfer.",
    )
    parser.add_argument(
        "--aerosol-resolution",
        type=float,
        default=1000.0,
        help="Aerosol retrieval resolution in meters.",
    )
    parser.add_argument(
        "--surface-prior-method",
        choices=("kernel_model", "whittaker", "monthly_database"),
        default="kernel_model",
        help="Surface-prior derivation route to use during M3.",
    )
    parser.add_argument(
        "--whittaker-lambda",
        type=float,
        default=10.0,
        help="Temporal smoothness strength for Route-A Whittaker BRDF priors.",
    )
    parser.add_argument(
        "--log-level",
        choices=("DEBUG", "INFO", "WARNING", "ERROR"),
        default="INFO",
        help="Logging level.",
    )
    return parser.parse_args()


def setup_logging(level: str) -> None:
    logging.basicConfig(
        level=getattr(logging, level),
        format="%(asctime)s %(levelname)s:%(name)s:%(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )


def _ensure_datetime(value: Any) -> datetime:
    if isinstance(value, datetime):
        return value
    if isinstance(value, str):
        # Handle both "...Z" and timezone-aware ISO strings.
        return datetime.fromisoformat(value.replace("Z", "+00:00"))
    raise TypeError(f"Expected observation_time datetime/ISO string, got {type(value).__name__}")


def _bounds_from_xy(da: xr.DataArray) -> tuple[float, float, float, float]:
    if "x" not in da.coords or "y" not in da.coords:
        raise ValueError("DataArray is missing x/y coordinates; cannot infer bounds.")
    if da.sizes.get("x", 0) == 0 or da.sizes.get("y", 0) == 0:
        raise ValueError("DataArray has empty x/y dimensions; cannot infer bounds.")

    x_vals = da.coords["x"].values
    y_vals = da.coords["y"].values
    return (
        float(x_vals.min()),
        float(y_vals.min()),
        float(x_vals.max()),
        float(y_vals.max()),
    )


def _slice_by_bounds(
    da: xr.DataArray,
    bounds: tuple[float, float, float, float],
) -> xr.DataArray:
    if "x" not in da.coords or "y" not in da.coords:
        return da

    xmin, ymin, xmax, ymax = bounds
    x_slice = slice(min(xmin, xmax), max(xmin, xmax))

    y_values = da.coords["y"].values
    descending_y = bool(y_values[0] > y_values[-1])
    y_slice = slice(max(ymin, ymax), min(ymin, ymax)) if descending_y else slice(
        min(ymin, ymax), max(ymin, ymax)
    )
    return da.sel(x=x_slice, y=y_slice)


def _clip_dataarray(
    da: xr.DataArray,
    bounds: tuple[float, float, float, float],
    bounds_crs: str,
) -> xr.DataArray:
    try:
        if hasattr(da, "rio") and da.rio.crs is not None:
            return da.rio.clip_box(*bounds, crs=bounds_crs)
    except Exception as exc:
        logger.debug("rio.clip_box failed; falling back to coord slicing: %s", exc)
    return _slice_by_bounds(da, bounds)


def _clip_dataset(
    ds: xr.Dataset,
    bounds: tuple[float, float, float, float],
    bounds_crs: str,
) -> xr.Dataset:
    clipped = {
        name: _clip_dataarray(var, bounds=bounds, bounds_crs=bounds_crs)
        for name, var in ds.data_vars.items()
    }
    out = xr.Dataset(clipped)
    out.attrs.update(ds.attrs)
    return out


def _coords_match(a: xr.DataArray, b: xr.DataArray) -> bool:
    if "x" not in a.coords or "y" not in a.coords or "x" not in b.coords or "y" not in b.coords:
        return False
    return (
        a.sizes.get("x") == b.sizes.get("x")
        and a.sizes.get("y") == b.sizes.get("y")
        and np.array_equal(a.coords["x"].values, b.coords["x"].values)
        and np.array_equal(a.coords["y"].values, b.coords["y"].values)
    )


def _harmonize_toa_grid(toa: xr.Dataset, template: xr.DataArray) -> xr.Dataset:
    """Resample all TOA bands to a common template grid."""
    out: dict[str, xr.DataArray] = {}

    for name, da in toa.data_vars.items():
        if _coords_match(da, template):
            out[name] = da
            continue
        try:
            out[name] = reproject_match(da, template, resampling="bilinear")
        except Exception:
            out[name] = da.interp(y=template.coords["y"], x=template.coords["x"], method="linear")

    merged = xr.Dataset(out)
    merged.attrs.update(toa.attrs)
    return merged


def _harmonize_mask_grid(mask: xr.DataArray, template: xr.DataArray) -> xr.DataArray:
    """Resample mask to template grid using nearest-neighbour."""
    if _coords_match(mask, template):
        return mask.astype(bool)
    mask_num = mask.astype(np.float32)
    if (
        hasattr(mask_num, "rio")
        and hasattr(template, "rio")
        and mask_num.rio.crs is None
        and template.rio.crs is not None
    ):
        mask_num = mask_num.rio.write_crs(template.rio.crs)
    try:
        out = reproject_match(mask_num, template, resampling="nearest")
        return (out > 0.5).astype(bool)
    except Exception:
        out = mask_num.interp(y=template.coords["y"], x=template.coords["x"], method="nearest")
        return (out > 0.5).astype(bool)


def _clip_geometry(
    geometry: GeometryAngles,
    bounds: tuple[float, float, float, float],
    bounds_crs: str,
) -> GeometryAngles:
    return GeometryAngles(
        sza=_clip_dataarray(geometry.sza, bounds, bounds_crs),
        saa=_clip_dataarray(geometry.saa, bounds, bounds_crs),
        vza=_clip_dataarray(geometry.vza, bounds, bounds_crs),
        vaa=_clip_dataarray(geometry.vaa, bounds, bounds_crs),
    )


def _build_observation_bundle(
    preprocess_output: dict[str, Any],
    sensor_config: SensorConfig,
    aoi: AOI | None,
) -> ObservationBundle:
    toa = preprocess_output["toa"]
    geometry = preprocess_output["geometry"]
    cloud_mask = preprocess_output["cloud_mask"]
    metadata = dict(preprocess_output["metadata"])
    metadata["observation_time"] = _ensure_datetime(metadata["observation_time"])

    first_band_name = next(iter(toa.data_vars))
    first_band = toa[first_band_name]

    if aoi is not None:
        target_crs = str(first_band.rio.crs) if hasattr(first_band, "rio") and first_band.rio.crs else aoi.crs
        bounds = aoi.get_bounds(target_crs=target_crs)
        toa = _clip_dataset(toa, bounds=bounds, bounds_crs=target_crs)
        geometry = _clip_geometry(geometry, bounds=bounds, bounds_crs=target_crs)
        cloud_mask = _clip_dataarray(cloud_mask, bounds=bounds, bounds_crs=target_crs)

    toa = _harmonize_toa_grid(toa, geometry.sza)
    cloud_mask = _harmonize_mask_grid(cloud_mask, geometry.sza)
    first_band_name = next(iter(toa.data_vars))
    first_band = toa[first_band_name]

    if first_band.sizes.get("x", 0) == 0 or first_band.sizes.get("y", 0) == 0:
        raise ValueError("AOI clipping produced an empty raster; widen --aoi-bbox.")

    crs = str(first_band.rio.crs) if hasattr(first_band, "rio") and first_band.rio.crs else "EPSG:4326"
    bounds: tuple[float, float, float, float]
    if hasattr(first_band, "rio") and first_band.rio.crs is not None:
        try:
            bounds = tuple(map(float, first_band.rio.bounds()))
        except Exception:
            bounds = _bounds_from_xy(first_band)
    else:
        bounds = _bounds_from_xy(first_band)

    return ObservationBundle(
        toa=toa,
        geometry=geometry,
        cloud_mask=cloud_mask.astype(bool),
        sensor_config=sensor_config,
        metadata=metadata,
        crs=crs,
        bounds=bounds,
    )


def _build_config(args: argparse.Namespace) -> SIACConfig:
    cams_dir: str | Path = args.cams_dir if args.cams_dir is not None else args.cache_dir / "cams"
    if isinstance(cams_dir, str):
        scheme = urlparse(cams_dir).scheme.lower()
        if scheme not in {"http", "https", "s3"}:
            cams_dir = Path(cams_dir).expanduser()
    if isinstance(cams_dir, Path):
        cams_dir.mkdir(parents=True, exist_ok=True)
    args.cache_dir.mkdir(parents=True, exist_ok=True)
    cams_cache_dir = args.cache_dir / "cams-cache"
    cams_cache_dir.mkdir(parents=True, exist_ok=True)

    return SIACConfig(
        sensor="s2",
        atmo_prior={"provider": "cams", "data_path": cams_dir, "cache_dir": cams_cache_dir},
        brdf={"provider": "mcd43"},
        surface_prior={
            "method": args.surface_prior_method,
            "whittaker_lambda": args.whittaker_lambda,
        },
        rt_model={"backend": "lut", "lut_path": args.lut_path},
        solver={"aerosol_resolution": args.aerosol_resolution},
        s2_data={
            "backend": args.backend,
            "cache_dir": args.cache_dir,
            "max_cloud_cover": args.max_cloud_cover,
            "processing_level": args.processing_level,
        },
    )


def _derive_l2a_folder_name(input_path: str | Path) -> str:
    """Derive an L2A-like scene folder name from SAFE/product identifier."""
    name = Path(str(input_path)).name
    if name.endswith(".SAFE"):
        name = name[:-5]

    if "_MSIL1C_" in name:
        name = name.replace("_MSIL1C_", "_MSIL2A_", 1)
    elif "_MSIL2A_" not in name:
        # Conservative fallback for non-standard IDs.
        parts = name.split("_", 1)
        name = f"{parts[0]}_MSIL2A_{parts[1]}" if len(parts) == 2 else f"{name}_MSIL2A"
    return name


def _resolve_output_dir(base_output_dir: Path, input_path: str | Path) -> Path:
    """Ensure output folder name follows S2?_MSIL2A_*T* scene pattern."""
    if _L2A_FOLDER_PATTERN.match(base_output_dir.name):
        return base_output_dir
    return base_output_dir / _derive_l2a_folder_name(input_path)


def main() -> int:
    args = parse_args()
    setup_logging(args.log_level)
    config = _build_config(args)

    aoi = AOI.from_bounds(tuple(args.aoi_bbox), crs=args.aoi_crs) if args.aoi_bbox else None

    logger.info("Resolving Sentinel-2 input from query: %s", args.query)
    input_path = resolve_s2_input(args.query, config)
    logger.info("Using local SAFE path: %s", input_path)
    output_dir = _resolve_output_dir(args.output_dir, input_path)
    logger.info("Resolved output directory: %s", output_dir)

    cloud_mask_config = config.cloud_mask.model_dump(exclude={"user_callable"})
    pp = Sentinel2Preprocessor(config={"cloud_mask": cloud_mask_config})
    obs_capture: dict[str, ObservationBundle] = {}

    def _preprocess(path: Path, runtime_aoi: AOI | None = None) -> ObservationBundle:
        raw = pp.preprocess(path)
        obs = _build_observation_bundle(raw, pp.sensor_config, runtime_aoi)
        obs_capture["obs"] = obs
        return obs

    logger.info("Starting SIAC full pipeline run")
    result = siac_process(
        config,
        Path(input_path),
        aoi=aoi,
        preprocessor=_preprocess,
    )

    output_dir.mkdir(parents=True, exist_ok=True)
    boa_paths = write_dataset(
        result.boa,
        output_dir / "boa",
        as_cog=True,
        compression=config.output.compression,
        dtype=config.output.boa_dtype,
    )

    atmo_ds = xr.Dataset(
        {
            "aot": result.aot.astype("float32"),
            "tcwv": result.tcwv.astype("float32"),
            "cloud_mask": result.cloud_mask.astype("uint8"),
        }
    )
    atmo_path = write_netcdf(atmo_ds, output_dir / "atmosphere.nc")
    qa_paths = write_dataset(
        xr.Dataset({"cloud_mask": result.cloud_mask.astype("uint8")}),
        output_dir / "qa",
        as_cog=True,
        compression=config.output.compression,
        dtype="uint8",
    )

    summary = {
        "query": args.query,
        "input_safe": str(input_path),
        "output_dir": str(output_dir),
        "n_boa_bands": len(result.boa.data_vars),
        "boa_files": {name: str(path) for name, path in boa_paths.items()},
        "atmosphere_file": str(atmo_path),
        "cloud_mask_file": str(qa_paths["cloud_mask"]) if "cloud_mask" in qa_paths else None,
        "aot_mean": float(result.aot.mean(skipna=True).values),
        "tcwv_mean": float(result.tcwv.mean(skipna=True).values),
    }
    summary_path = output_dir / "run_summary.json"
    summary["stac_item_file"] = str(output_dir / "item.json")
    summary_path.write_text(json.dumps(summary, indent=2))

    obs = obs_capture.get("obs")
    if obs is None:
        raise RuntimeError("Observation bundle was not captured during preprocessing; cannot write STAC item.")
    item_path = write_stac_item(
        obs,
        result,
        output_dir=output_dir,
        boa_assets=boa_paths,
        atmosphere_asset=atmo_path,
        qa_assets=qa_paths,
        summary_asset=summary_path,
        input_href=input_path,
        item_id=output_dir.name,
    )
    summary["stac_item_file"] = str(item_path)
    summary_path.write_text(json.dumps(summary, indent=2))

    print(json.dumps(summary, indent=2))
    logger.info("Run complete. Summary: %s", summary_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
