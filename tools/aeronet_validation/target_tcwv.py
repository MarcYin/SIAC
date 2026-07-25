"""Resolve fixed target-scene TCWV choices for aerosol experiments."""

from __future__ import annotations

import json
import math
import re
from dataclasses import asdict, dataclass
from datetime import datetime, timedelta
from pathlib import Path
from typing import Any

import numpy as np
import xarray as xr

DEFAULT_ROOT = Path("/gws/ssde/j25a/nceo_isp/public/siac_refactor")
DEFAULT_WVP_DIR = "s2_l2a_wvp_campaign250"


def product_acquisition_time(product_id: object) -> datetime | None:
    """Parse the sensing timestamp shared by matching L1C and L2A products."""

    if not isinstance(product_id, str):
        return None
    match = re.search(r"MSIL[12][AC]_(\d{8}T\d{6})", product_id)
    if match is None:
        return None
    return datetime.strptime(match.group(1), "%Y%m%dT%H%M%S")


@dataclass(frozen=True)
class TargetTCWV:
    """One fixed TCWV value and its experiment provenance."""

    value_cm: float
    mode: str
    source: str
    scale: float
    input_value_cm: float
    record_path: str | None = None

    def to_dict(self) -> dict[str, float | str | None]:
        return asdict(self)


def normalise_s2_l2a_wvp(raw: object) -> tuple[float | None, str | None]:
    """Convert a GEE Sentinel-2 L2A WVP value to centimetres.

    ``COPERNICUS/S2_SR_HARMONIZED.WVP`` normally stores centimetres with a
    0.001 scale. Small already-physical values are retained to make the
    collector robust to an upstream API applying the scale automatically.
    """

    if raw is None:
        return None, None
    try:
        value = float(raw)
    except (TypeError, ValueError):
        return None, None
    if not math.isfinite(value) or value <= 0.0:
        return None, None
    if value > 20.0:
        value /= 1000.0
        scale = "dn_div_1000"
    else:
        scale = "physical_cm"
    if not 0.05 <= value <= 15.0:
        return None, None
    return value, scale


def normalise_s2_l2a_wvp_field(
    raw: object,
    *,
    fallback_cm: float,
) -> tuple[np.ndarray, dict[str, float | int | bool | str]]:
    """Normalise an L2A WVP raster and explicitly fill invalid source pixels.

    Sentinel-2 L2A WVP is normally stored as an integer band scaled by 0.001
    cm. Nodata and implausible values are replaced with the validated
    same-scene station-buffer median. The returned statistics make that
    fallback visible in every retrieval receipt.
    """

    fallback, _ = normalise_s2_l2a_wvp(fallback_cm)
    if fallback is None:
        raise ValueError(f"invalid L2A WVP fallback {fallback_cm!r} cm")

    values = np.asarray(raw, dtype=np.float32)
    if values.ndim != 2:
        raise ValueError(f"L2A WVP raster must be two-dimensional, got {values.shape}")
    physical = np.where(values > 20.0, values * np.float32(0.001), values)
    valid = np.isfinite(physical) & (physical >= 0.05) & (physical <= 15.0)
    valid_count = int(np.count_nonzero(valid))
    output = np.where(valid, physical, np.float32(fallback)).astype(np.float32)
    return output, {
        "source_scale": "dn_div_1000_or_physical_cm",
        "valid_pixel_count": valid_count,
        "total_pixel_count": int(values.size),
        "valid_fraction": float(valid_count / values.size) if values.size else 0.0,
        "fallback_pixel_count": int(values.size - valid_count),
        "fallback_cm": float(fallback),
        "all_pixels_fallback": valid_count == 0,
    }


def l2a_wvp_grid_from_template(template: xr.DataArray, *, crs: str) -> dict[str, object]:
    """Build a GEE pixel request grid exactly aligned to an SIAC atmo template."""

    if tuple(template.dims) != ("y", "x"):
        raise ValueError(
            "spatial L2A WVP requires a two-dimensional atmospheric template with y/x dimensions"
        )
    if "x" not in template.coords or "y" not in template.coords:
        raise ValueError("spatial L2A WVP requires x/y coordinates on the atmospheric template")

    x = np.asarray(template.x.values, dtype=np.float64)
    y = np.asarray(template.y.values, dtype=np.float64)
    if x.size < 2 or y.size < 2:
        raise ValueError("spatial L2A WVP requires at least two x/y grid cells")
    dx = np.diff(x)
    dy = np.diff(y)
    resolution = float(np.median(np.abs(dx)))
    tolerance = max(1e-5, resolution * 1e-5)
    if (
        not np.isfinite(resolution)
        or resolution <= 0.0
        or not np.allclose(dx, resolution, rtol=0.0, atol=tolerance)
        or not np.allclose(dy, -resolution, rtol=0.0, atol=tolerance)
    ):
        raise ValueError("spatial L2A WVP requires a regular north-up square atmospheric grid")
    return {
        "crs": str(crs),
        "res": resolution,
        "x0": float(x[0] - 0.5 * resolution),
        "y1": float(y[0] + 0.5 * resolution),
        "W": int(x.size),
        "H": int(y.size),
    }


def _bounds_to_lonlat(
    bounds: tuple[float, float, float, float],
    crs: str,
) -> tuple[float, float, float, float]:
    if str(crs).upper() in {"EPSG:4326", "WGS84", "OGC:CRS84"}:
        return tuple(float(value) for value in bounds)
    from siac.geo.reprojection import transform_bounds

    return tuple(float(value) for value in transform_bounds(bounds, str(crs), "EPSG:4326"))


def fetch_matching_l2a_wvp_field(
    ee: Any,
    *,
    bounds: tuple[float, float, float, float],
    crs: str,
    template: xr.DataArray,
    product_id: str,
    mgrs_tile: str,
    observation_time: datetime,
    fallback_cm: float,
    max_time_delta_seconds: float = 300.0,
) -> tuple[xr.DataArray, dict[str, Any]]:
    """Fetch matching L2A WVP on the exact M2 atmospheric grid.

    The L2A acquisition is selected by the L1C sensing timestamp plus MGRS
    tile, rather than by an arbitrary same-day granule. This field is then
    shared by the surface anchor correction and AOD cost cube.
    """

    from bestpixel._gee import get_patch

    product_sensing_time = product_acquisition_time(product_id)
    target_time = product_sensing_time or observation_time
    day = target_time.date().isoformat()
    next_day = (target_time.date() + timedelta(days=1)).isoformat()
    tile = str(mgrs_tile).removeprefix("T")
    if not tile:
        raise ValueError("matching L2A WVP requires an MGRS tile")
    west, south, east, north = _bounds_to_lonlat(bounds, crs)
    geom = ee.Geometry.Rectangle([west, south, east, north])
    collection = (
        ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
        .filterBounds(geom)
        .filterDate(day, next_day)
        .filter(ee.Filter.eq("MGRS_TILE", tile))
    )
    if product_sensing_time is not None:
        platform = str(product_id).split("_", 1)[0]
        l2a_prefix = f"{platform}_MSIL2A_{product_sensing_time:%Y%m%dT%H%M%S}"
        # GEE's system:time_start can refer to a tile's scan time rather than
        # the product sensing timestamp. PRODUCT_ID carries the invariant L1C
        # to L2A acquisition identity, while its processing timestamp may vary.
        collection = collection.filter(ee.Filter.stringStartsWith("PRODUCT_ID", l2a_prefix))
    count = int(collection.size().getInfo())
    if count == 0:
        raise RuntimeError(f"no matching L2A product for tile {tile} on {day}")

    target_ms = ee.Date(target_time.strftime("%Y-%m-%dT%H:%M:%S")).millis()

    def _with_time_delta(image: Any) -> Any:
        delta = ee.Number(image.get("system:time_start")).subtract(target_ms).abs()
        return image.set("_siac_time_delta_ms", delta)

    image = ee.Image(collection.map(_with_time_delta).sort("_siac_time_delta_ms").first())
    info = ee.Dictionary(
        {
            "system_index": image.get("system:index"),
            "product_id": image.get("PRODUCT_ID"),
            "system_time_start": image.get("system:time_start"),
            "time_delta_ms": image.get("_siac_time_delta_ms"),
        }
    ).getInfo()
    system_index = info.get("system_index")
    if not system_index:
        raise RuntimeError(f"matching L2A image for tile {tile} on {day} has no system:index")
    delta_ms = info.get("time_delta_ms")
    l2a_product_id = info.get("product_id")
    l2a_sensing_time = product_acquisition_time(l2a_product_id)
    if product_sensing_time is not None and l2a_sensing_time != product_sensing_time:
        raise RuntimeError(
            "matching L2A PRODUCT_ID does not share the L1C sensing timestamp: "
            f"{l2a_product_id!r}"
        )
    if (
        product_sensing_time is None
        and (delta_ms is None or float(delta_ms) > float(max_time_delta_seconds) * 1000.0)
    ):
        raise RuntimeError(
            "matching L2A WVP acquisition differs from the observation time by "
            f"{delta_ms!r} ms"
        )

    grid = l2a_wvp_grid_from_template(template, crs=crs)
    raw = np.asarray(
        get_patch(
            ee,
            f"COPERNICUS/S2_SR_HARMONIZED/{system_index}",
            ["WVP"],
            grid=grid,
        )[0],
        dtype=np.float32,
    )
    if raw.shape != template.shape:
        raise RuntimeError(
            f"matching L2A WVP grid has shape {raw.shape}, expected {template.shape}"
        )
    values, field_info = normalise_s2_l2a_wvp_field(raw, fallback_cm=fallback_cm)
    asset_id = f"COPERNICUS/S2_SR_HARMONIZED/{system_index}"
    field = xr.DataArray(
        values,
        dims=template.dims,
        coords=template.coords,
        attrs={
            **template.attrs,
            # M4's reprojection cache must not substitute a scalar CAMS/L2A
            # field for this matching per-pixel WVP source on the same scene.
            "siac_atmo_cache_identity": f"l2a-wvp:{asset_id}",
        },
        name=template.name,
    )
    finite = values[np.isfinite(values)]
    return field, {
        "mode": "l2a_spatial",
        "source": "COPERNICUS/S2_SR_HARMONIZED.WVP",
        "asset_id": asset_id,
        "l2a_product_id": l2a_product_id,
        "l2a_system_time_start": info.get("system_time_start"),
        "time_delta_ms": float(delta_ms),
        "product_sensing_time_match": l2a_sensing_time == product_sensing_time,
        "grid_resolution_m": float(grid["res"]),
        "field_mean_cm": float(np.mean(finite)) if finite.size else None,
        "field_median_cm": float(np.median(finite)) if finite.size else None,
        **field_info,
    }


def resolve_target_tcwv(
    matchup_id: str,
    mode: str,
    *,
    root: Path = DEFAULT_ROOT,
    fixed_cm: float = 2.0,
    wvp_dir: str = DEFAULT_WVP_DIR,
) -> TargetTCWV | None:
    """Resolve a target TCWV mode without silently changing its source.

    ``prior`` preserves the atmospheric provider and returns ``None``.
    ``fixed`` supplies ``fixed_cm``. ``l2a`` and ``l2a_085`` read a collected
    same-scene Sen2Cor WVP record; missing or failed records are hard errors so
    an experiment cannot accidentally mix L2A and CAMS water vapour.
    """

    normalized_mode = mode.strip().lower()
    if normalized_mode in {"", "prior", "cams"}:
        return None
    if normalized_mode == "fixed":
        value = float(fixed_cm)
        if not math.isfinite(value) or not 0.05 <= value <= 15.0:
            raise ValueError(f"invalid fixed target TCWV {fixed_cm!r} cm")
        return TargetTCWV(
            value_cm=value,
            mode="fixed",
            source="experiment_fixed",
            scale=1.0,
            input_value_cm=value,
        )
    if normalized_mode not in {"l2a", "l2a_085"}:
        raise ValueError(f"unsupported target TCWV mode {mode!r}")

    record_path = root / wvp_dir / f"{matchup_id}.json"
    record = json.loads(record_path.read_text(encoding="utf-8"))
    if record.get("status") != "OK":
        raise RuntimeError(f"target L2A WVP record is not OK: {record_path}")
    input_value = float(record["tcwv_cm"])
    if not math.isfinite(input_value) or not 0.05 <= input_value <= 15.0:
        raise ValueError(f"invalid target L2A WVP {input_value!r} cm in {record_path}")
    factor = 0.85 if normalized_mode == "l2a_085" else 1.0
    value = input_value * factor
    return TargetTCWV(
        value_cm=value,
        mode=normalized_mode,
        source=str(record.get("source") or "COPERNICUS/S2_SR_HARMONIZED.WVP"),
        scale=factor,
        input_value_cm=input_value,
        record_path=str(record_path),
    )
