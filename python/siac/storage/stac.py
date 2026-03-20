"""STAC item writer for SIAC outputs."""

from __future__ import annotations

import json
import os
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np
from rasterio.warp import transform_bounds

if TYPE_CHECKING:
    from rasterio.crs import CRS

    from siac.domain import SensorBand
    from siac.runtime import CorrectionResult, ObservationBundle

_STAC_VERSION = "1.0.0"
_EO_EXTENSION = "https://stac-extensions.github.io/eo/v2.0.0/schema.json"
_PROJECTION_EXTENSION = "https://stac-extensions.github.io/projection/v2.0.0/schema.json"
_VIEW_EXTENSION = "https://stac-extensions.github.io/view/v1.1.0/schema.json"
_COG_MEDIA_TYPE = "image/tiff; application=geotiff; profile=cloud-optimized"
_BAND_COMMON_NAMES = {
    "B01": "coastal",
    "B02": "blue",
    "B03": "green",
    "B04": "red",
    "B05": "rededge",
    "B06": "rededge",
    "B07": "rededge",
    "B08": "nir",
    "B8A": "nir08",
    "B10": "cirrus",
    "B11": "swir16",
    "B12": "swir22",
    "B1": "coastal",
    "B2": "blue",
    "B3": "green",
    "B4": "red",
    "B5": "nir",
    "B6": "swir16",
    "B7": "swir22",
}


def _isoformat_utc(value: datetime) -> str:
    if value.tzinfo is None:
        value = value.replace(tzinfo=timezone.utc)
    else:
        value = value.astimezone(timezone.utc)
    return value.isoformat().replace("+00:00", "Z")


def _relative_href(path: str | Path, base_dir: Path) -> str:
    rel = os.path.relpath(Path(path), base_dir)
    return Path(rel).as_posix()


def _parse_satellite_id(input_href: str | Path | None, metadata: dict[str, Any], fallback: str) -> str:
    if input_href is not None:
        match = re.match(r"^(S2[A-Z]|L\d+)", Path(str(input_href)).name)
        if match is not None:
            return match.group(1)
    value = metadata.get("satellite")
    if isinstance(value, str) and value:
        return value
    return fallback


def _platform_name(satellite_id: str) -> str | None:
    if satellite_id.startswith("S2"):
        return f"sentinel-2{satellite_id[-1].lower()}"
    if satellite_id == "L8":
        return "landsat-8"
    return None


def _constellation_name(satellite_id: str) -> str | None:
    if satellite_id.startswith("S2"):
        return "sentinel-2"
    if satellite_id.startswith("L"):
        return "landsat"
    return None


def _safe_float(value: Any) -> float | None:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return None
    if not np.isfinite(out):
        return None
    return out


def _mean_deg(da) -> float | None:
    data = np.asarray(da.values, dtype=np.float64)
    if data.size == 0:
        return None
    mean = np.nanmean(data)
    if not np.isfinite(mean):
        return None
    return float(np.rad2deg(mean))


def _cloud_cover_percent(mask) -> float | None:
    mean = _safe_float(mask.mean(skipna=True).values)
    if mean is None:
        return None
    return float(np.clip(mean * 100.0, 0.0, 100.0))


def _native_bounds(first_band, fallback_bounds: tuple[float, float, float, float]) -> tuple[float, float, float, float]:
    try:
        return tuple(map(float, first_band.rio.bounds()))
    except Exception:
        return tuple(map(float, fallback_bounds))


def _wgs84_bounds_and_geometry(
    native_bounds: tuple[float, float, float, float],
    crs: CRS | str | None,
) -> tuple[list[float], dict[str, Any]]:
    if crs is None:
        xmin, ymin, xmax, ymax = native_bounds
    else:
        transformed = transform_bounds(crs, "EPSG:4326", *native_bounds, densify_pts=21)
        xmin, ymin, xmax, ymax = map(float, transformed)
    bbox = [xmin, ymin, xmax, ymax]
    geometry = {
        "type": "Polygon",
        "coordinates": [[
            [xmin, ymin],
            [xmin, ymax],
            [xmax, ymax],
            [xmax, ymin],
            [xmin, ymin],
        ]],
    }
    return bbox, geometry


def _proj_properties(first_band, native_bounds: tuple[float, float, float, float]) -> dict[str, Any]:
    props: dict[str, Any] = {
        "proj:bbox": [float(v) for v in native_bounds],
        "proj:shape": [int(first_band.sizes["y"]), int(first_band.sizes["x"])],
    }
    try:
        epsg = first_band.rio.crs.to_epsg()
    except Exception:
        epsg = None
    if epsg is not None:
        props["proj:epsg"] = int(epsg)
    try:
        transform = first_band.rio.transform()
        props["proj:transform"] = [float(v) for v in transform[:6]]
    except Exception:
        pass
    return props


def _gsd(first_band) -> float | None:
    try:
        res_x, res_y = first_band.rio.resolution()
    except Exception:
        return None
    if not np.isfinite(res_x) or not np.isfinite(res_y):
        return None
    return float(max(abs(res_x), abs(res_y)))


def _band_metadata(band: SensorBand) -> dict[str, Any]:
    metadata = {
        "name": band.name,
        "center_wavelength": float(band.wavelength_um),
        "full_width_half_max": float(band.bandwidth / 1000.0),
    }
    common_name = _BAND_COMMON_NAMES.get(band.name)
    if common_name is not None:
        metadata["common_name"] = common_name
    return metadata


def _file_size(path: Path) -> int | None:
    try:
        return int(path.stat().st_size)
    except OSError:
        return None


def _asset_dict(
    href: str,
    *,
    title: str,
    media_type: str,
    roles: list[str],
    extra_fields: dict[str, Any] | None = None,
    file_size: int | None = None,
) -> dict[str, Any]:
    asset: dict[str, Any] = {
        "href": href,
        "type": media_type,
        "title": title,
        "roles": roles,
    }
    if file_size is not None:
        asset["file:size"] = file_size
    if extra_fields:
        asset.update(extra_fields)
    return asset


def build_stac_item(
    obs: ObservationBundle,
    result: CorrectionResult,
    *,
    output_dir: str | Path,
    boa_assets: dict[str, Path],
    atmosphere_asset: Path | None = None,
    qa_assets: dict[str, Path] | None = None,
    summary_asset: Path | None = None,
    input_href: str | Path | None = None,
    item_id: str | None = None,
    item_href: str | Path | None = None,
) -> dict[str, Any]:
    """Build a STAC Item for a SIAC output scene."""
    output_dir = Path(output_dir)
    item_id = item_id or output_dir.name
    item_href = Path(item_href) if item_href is not None else output_dir / "item.json"

    first_band_name = next(iter(result.boa.data_vars))
    first_band = result.boa[first_band_name]
    native_bounds = _native_bounds(first_band, obs.bounds)
    crs = None
    try:
        crs = first_band.rio.crs
    except Exception:
        crs = obs.crs

    bbox, geometry = _wgs84_bounds_and_geometry(native_bounds, crs)
    satellite_id = _parse_satellite_id(input_href, obs.metadata, obs.sensor_config.satellite_id)
    platform = _platform_name(satellite_id)
    constellation = _constellation_name(satellite_id)

    observation_time = obs.metadata["observation_time"]
    if not isinstance(observation_time, datetime):
        raise TypeError("ObservationBundle metadata must contain datetime observation_time")

    gsd = _gsd(first_band)
    cloud_cover = _cloud_cover_percent(obs.cloud_mask)
    valid_pixel_percent = _cloud_cover_percent(~result.cloud_mask.astype(bool))
    mean_sza = _mean_deg(obs.geometry.sza)
    mean_saa = _mean_deg(obs.geometry.saa)
    mean_vza = _mean_deg(obs.geometry.vza)
    mean_vaa = _mean_deg(obs.geometry.vaa)
    mean_raa = _mean_deg(obs.geometry.raa)

    properties: dict[str, Any] = {
        "datetime": _isoformat_utc(observation_time),
        "created": _isoformat_utc(datetime.now(timezone.utc)),
        "instruments": [str(obs.sensor_config.sensor_id).lower()],
        "siac:aot_mean": float(result.aot.mean(skipna=True).values),
        "siac:tcwv_mean": float(result.tcwv.mean(skipna=True).values),
    }
    if platform is not None:
        properties["platform"] = platform
    if constellation is not None:
        properties["constellation"] = constellation
    if gsd is not None:
        properties["gsd"] = gsd
    if cloud_cover is not None:
        properties["eo:cloud_cover"] = cloud_cover
    if valid_pixel_percent is not None:
        properties["siac:masked_pixel_percent"] = valid_pixel_percent
    if mean_saa is not None:
        properties["view:sun_azimuth"] = mean_saa
    if mean_sza is not None:
        properties["view:sun_elevation"] = 90.0 - mean_sza
        properties["siac:mean_sun_zenith"] = mean_sza
    if mean_vaa is not None:
        properties["view:azimuth"] = mean_vaa
    if mean_vza is not None:
        properties["view:off_nadir"] = mean_vza
    if mean_raa is not None:
        properties["siac:mean_relative_azimuth"] = mean_raa
    processing_time = _safe_float(result.metadata.get("processing_time_s"))
    if processing_time is not None:
        properties["siac:processing_time_s"] = processing_time

    for source_key, target_key in (
        ("tile_id", "siac:tile_id"),
        ("processing_baseline", "siac:processing_baseline"),
        ("sensor", "siac:sensor"),
    ):
        value = obs.metadata.get(source_key)
        if isinstance(value, (str, int, float)):
            properties[target_key] = value
    properties["siac:satellite"] = satellite_id

    properties.update(_proj_properties(first_band, native_bounds))

    assets: dict[str, Any] = {}
    for name, path in boa_assets.items():
        band = obs.sensor_config.get_band(name)
        assets[name] = _asset_dict(
            _relative_href(path, output_dir),
            title=f"BOA reflectance {name}",
            media_type=_COG_MEDIA_TYPE,
            roles=["data", "reflectance"],
            extra_fields={"eo:bands": [_band_metadata(band)]},
            file_size=_file_size(path),
        )

    if atmosphere_asset is not None:
        assets["atmosphere"] = _asset_dict(
            _relative_href(atmosphere_asset, output_dir),
            title="Solved atmospheric state",
            media_type="application/x-netcdf",
            roles=["data", "metadata"],
            file_size=_file_size(atmosphere_asset),
        )

    if qa_assets:
        for name, path in qa_assets.items():
            if name == "cloud_mask":
                assets["cloud-mask"] = _asset_dict(
                    _relative_href(path, output_dir),
                    title="Cloud mask",
                    media_type=_COG_MEDIA_TYPE,
                    roles=["metadata", "cloud-mask"],
                    file_size=_file_size(path),
                )

    if summary_asset is not None:
        assets["summary"] = _asset_dict(
            _relative_href(summary_asset, output_dir),
            title="Run summary",
            media_type="application/json",
            roles=["metadata"],
            file_size=_file_size(summary_asset),
        )

    links: list[dict[str, Any]] = [
        {
            "rel": "self",
            "href": _relative_href(item_href, output_dir),
            "type": "application/geo+json",
        }
    ]
    if input_href is not None:
        links.append(
            {
                "rel": "derived_from",
                "href": _relative_href(input_href, output_dir),
            }
        )

    return {
        "type": "Feature",
        "stac_version": _STAC_VERSION,
        "stac_extensions": [_EO_EXTENSION, _PROJECTION_EXTENSION, _VIEW_EXTENSION],
        "id": item_id,
        "bbox": bbox,
        "geometry": geometry,
        "properties": properties,
        "links": links,
        "assets": assets,
    }


def write_stac_item(
    obs: ObservationBundle,
    result: CorrectionResult,
    *,
    output_dir: str | Path,
    boa_assets: dict[str, Path],
    atmosphere_asset: Path | None = None,
    qa_assets: dict[str, Path] | None = None,
    summary_asset: Path | None = None,
    input_href: str | Path | None = None,
    item_id: str | None = None,
    item_href: str | Path | None = None,
) -> Path:
    """Write a STAC Item JSON document for a SIAC output scene."""
    output_dir = Path(output_dir)
    item_path = Path(item_href) if item_href is not None else output_dir / "item.json"
    item = build_stac_item(
        obs,
        result,
        output_dir=output_dir,
        boa_assets=boa_assets,
        atmosphere_asset=atmosphere_asset,
        qa_assets=qa_assets,
        summary_asset=summary_asset,
        input_href=input_href,
        item_id=item_id,
        item_href=item_path,
    )
    item_path.write_text(json.dumps(item, indent=2), encoding="utf-8")
    return item_path
