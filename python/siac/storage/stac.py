"""STAC item writer for SIAC outputs."""

from __future__ import annotations

import logging
import os
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np
from rasterio.warp import transform_bounds

logger = logging.getLogger(__name__)

if TYPE_CHECKING:
    from rasterio.crs import CRS

    from siac.domain import SensorBand
    from siac.runtime import CorrectionResult

_STAC_VERSION = "1.1.0"
_EO_EXTENSION = "https://stac-extensions.github.io/eo/v2.0.0/schema.json"
_PROJECTION_EXTENSION = "https://stac-extensions.github.io/projection/v2.0.0/schema.json"
_VIEW_EXTENSION = "https://stac-extensions.github.io/view/v1.1.0/schema.json"
_PROCESSING_EXTENSION = "https://stac-extensions.github.io/processing/v1.2.0/schema.json"
_COG_MEDIA_TYPE = "image/tiff; application=geotiff; profile=cloud-optimized"
_GEOTIFF_MEDIA_TYPE = "image/tiff; application=geotiff"
_NETCDF_MEDIA_TYPE = "application/x-netcdf"
_ZARR_MEDIA_TYPE = "application/vnd+zarr"
_JSON_MEDIA_TYPE = "application/json"

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
    "B09": "nir09",
    "B10": "cirrus",
    "B11": "swir16",
    "B12": "swir22",
    # Landsat band names
    "B1": "coastal",
    "B2": "blue",
    "B3": "green",
    "B4": "red",
    "B5": "nir",
    "B6": "swir16",
    "B7": "swir22",
}


# ---------------------------------------------------------------------------
# Utility helpers
# ---------------------------------------------------------------------------


def _isoformat_utc(value: datetime) -> str:
    if value.tzinfo is None:
        value = value.replace(tzinfo=timezone.utc)
    else:
        value = value.astimezone(timezone.utc)
    return value.isoformat().replace("+00:00", "Z")


def _relative_href(path: str | Path, base_dir: Path) -> str:
    rel = os.path.relpath(Path(path), base_dir)
    return Path(rel).as_posix()


def _parse_satellite_id(
    input_href: str | Path | None, metadata: dict[str, Any], fallback: str
) -> str:
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
    if satellite_id == "L9":
        return "landsat-9"
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


def _mean_deg(da: Any) -> float | None:
    data = np.asarray(da.values, dtype=np.float64)
    if data.size == 0:
        return None
    mean = np.nanmean(data)
    if not np.isfinite(mean):
        return None
    return float(np.rad2deg(mean))


def _cloud_cover_percent(mask: Any) -> float | None:
    mean = _safe_float(mask.mean(skipna=True).values)
    if mean is None:
        return None
    return float(np.clip(mean * 100.0, 0.0, 100.0))


def _native_bounds(
    first_band: Any, fallback_bounds: tuple[float, float, float, float]
) -> tuple[float, float, float, float]:
    try:
        return tuple(map(float, first_band.rio.bounds()))
    except (AttributeError, ValueError, RuntimeError, OSError) as exc:
        logger.debug("rio.bounds() unavailable, using fallback bounds: %s", exc)
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
        "coordinates": [
            [
                [xmin, ymin],
                [xmin, ymax],
                [xmax, ymax],
                [xmax, ymin],
                [xmin, ymin],
            ]
        ],
    }
    return bbox, geometry


def _proj_properties(
    first_band: Any, native_bounds: tuple[float, float, float, float]
) -> dict[str, Any]:
    props: dict[str, Any] = {
        "proj:bbox": [float(v) for v in native_bounds],
        "proj:shape": [int(first_band.sizes["y"]), int(first_band.sizes["x"])],
    }
    try:
        epsg = first_band.rio.crs.to_epsg()
    except (AttributeError, ValueError, RuntimeError) as exc:
        logger.debug("rio.crs.to_epsg() unavailable: %s", exc)
        epsg = None
    if epsg is not None:
        props["proj:epsg"] = int(epsg)
    try:
        transform = first_band.rio.transform()
        props["proj:transform"] = [float(v) for v in transform[:6]]
    except (AttributeError, ValueError, RuntimeError) as exc:
        logger.debug("rio.transform() unavailable: %s", exc)
    return props


def _gsd(first_band: Any) -> float | None:
    try:
        res_x, res_y = first_band.rio.resolution()
    except (AttributeError, ValueError, RuntimeError) as exc:
        logger.debug("rio.resolution() unavailable: %s", exc)
        return None
    if not np.isfinite(res_x) or not np.isfinite(res_y):
        return None
    return float(max(abs(res_x), abs(res_y)))


def _band_metadata(band: SensorBand) -> dict[str, Any]:
    metadata: dict[str, Any] = {
        "name": band.name,
        "center_wavelength": float(band.wavelength_um),
        "full_width_half_max": float(band.bandwidth / 1000.0),
    }
    common_name = _BAND_COMMON_NAMES.get(band.name)
    if common_name is not None:
        metadata["common_name"] = common_name
    return metadata


def _band_metadata_from_name(band_name: str) -> dict[str, Any]:
    """Minimal band metadata when only the band name is available."""
    metadata: dict[str, Any] = {"name": band_name}
    common_name = _BAND_COMMON_NAMES.get(band_name)
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


def _media_type_for_artifact(path: Path) -> str:
    suffix = path.suffix.lower()
    if suffix == ".tif":
        return _COG_MEDIA_TYPE
    if suffix == ".nc":
        return _NETCDF_MEDIA_TYPE
    if suffix == ".zarr" or str(path).endswith(".zarr"):
        return _ZARR_MEDIA_TYPE
    if suffix == ".json":
        return _JSON_MEDIA_TYPE
    if suffix == ".png":
        return "image/png"
    return "application/octet-stream"


# ---------------------------------------------------------------------------
# Build STAC Item from CorrectionResult
# ---------------------------------------------------------------------------


def build_stac_item_from_result(
    result: CorrectionResult,
    *,
    output_dir: str | Path,
    artifacts: dict[str, Path],
    item_id: str | None = None,
) -> dict[str, Any]:
    """Build a STAC Item from a ``CorrectionResult`` with enriched metadata.

    The ``result.metadata`` dict is expected to contain observation info
    (sensor_config, geometry, crs, bounds, observation_time, satellite, etc.)
    as populated by the pipeline.
    """
    output_dir = Path(output_dir)
    metadata = result.metadata

    # Spatial reference from first BOA band (or AOT if BOA is empty)
    if result.boa.data_vars:
        first_band_name = next(iter(result.boa.data_vars))
        first_band = result.boa[first_band_name]
    else:
        first_band = result.aot
    fallback_bounds = metadata.get("bounds", (0.0, 0.0, 1.0, 1.0))
    native_bounds = _native_bounds(first_band, fallback_bounds)
    crs = None
    try:
        crs = first_band.rio.crs
    except Exception:
        crs = metadata.get("crs")

    bbox, geometry_geojson = _wgs84_bounds_and_geometry(native_bounds, crs)

    # Identifiers
    satellite_id = _parse_satellite_id(
        metadata.get("input_path"), metadata, metadata.get("satellite", "SAT")
    )
    platform = _platform_name(satellite_id)
    constellation = _constellation_name(satellite_id)

    # Observation time
    observation_time = metadata.get("observation_time")
    dt_str = _isoformat_utc(observation_time) if isinstance(observation_time, datetime) else None

    # Item ID
    if item_id is None:
        item_id = output_dir.name

    # View geometry (from metadata if available)
    obs_geometry = metadata.get("geometry")
    mean_sza = _mean_deg(obs_geometry.sza) if obs_geometry is not None else None
    mean_saa = _mean_deg(obs_geometry.saa) if obs_geometry is not None else None
    mean_vza = _mean_deg(obs_geometry.vza) if obs_geometry is not None else None
    mean_vaa = _mean_deg(obs_geometry.vaa) if obs_geometry is not None else None

    # Cloud cover
    cloud_cover = _cloud_cover_percent(result.cloud_mask)
    gsd_val = _gsd(first_band)

    # Sensor config
    sensor_config = metadata.get("sensor_config")

    # ---- Properties ----
    properties: dict[str, Any] = {
        "created": _isoformat_utc(datetime.now(timezone.utc)),
    }
    if dt_str is not None:
        properties["datetime"] = dt_str
    else:
        properties["datetime"] = None
    if platform is not None:
        properties["platform"] = platform
    if constellation is not None:
        properties["constellation"] = constellation
    if sensor_config is not None:
        properties["instruments"] = [str(sensor_config.sensor_id).lower()]
    if gsd_val is not None:
        properties["gsd"] = gsd_val

    # eo extension
    if cloud_cover is not None:
        properties["eo:cloud_cover"] = cloud_cover

    # eo:bands summary (all BOA bands)
    eo_bands: list[dict[str, Any]] = []
    for band_name in result.boa.data_vars:
        if sensor_config is not None:
            try:
                band_obj = sensor_config.get_band(band_name)
                eo_bands.append(_band_metadata(band_obj))
                continue
            except (KeyError, AttributeError):
                pass
        eo_bands.append(_band_metadata_from_name(band_name))
    if eo_bands:
        properties["eo:bands"] = eo_bands

    # view extension
    if mean_saa is not None:
        properties["view:sun_azimuth"] = mean_saa
    if mean_sza is not None:
        properties["view:sun_elevation"] = 90.0 - mean_sza
    if mean_vaa is not None:
        properties["view:azimuth"] = mean_vaa
    if mean_vza is not None:
        properties["view:off_nadir"] = mean_vza

    # projection extension
    properties.update(_proj_properties(first_band, native_bounds))

    # processing extension
    # Read version from package metadata rather than a hard-coded literal so
    # we don't drift on every release (REVIEW.md §3.7 stac.py:383).
    try:
        from siac import __version__ as _siac_version
    except ImportError:  # pragma: no cover - defensive
        _siac_version = "unknown"
    properties["processing:software"] = {"SIAC": _siac_version}
    processing_time = _safe_float(result.diagnostics.processing_time_s)
    if processing_time is not None:
        properties["siac:processing_time_s"] = processing_time

    # SIAC-specific properties. Drop NaN values (which are not valid JSON via
    # the stdlib ``json`` module) so the STAC item is serialisable
    # (REVIEW.md §3.7 stac.py:389-390).
    aot_mean = float(result.aot.mean(skipna=True).values)
    tcwv_mean = float(result.tcwv.mean(skipna=True).values)
    if np.isfinite(aot_mean):
        properties["siac:aot_mean"] = aot_mean
    if np.isfinite(tcwv_mean):
        properties["siac:tcwv_mean"] = tcwv_mean
    properties["siac:satellite"] = satellite_id
    for source_key, target_key in (
        ("tile_id", "siac:tile_id"),
        ("processing_baseline", "siac:processing_baseline"),
        ("sensor", "siac:sensor"),
    ):
        value = metadata.get(source_key)
        if isinstance(value, (str, int, float)):
            properties[target_key] = value

    # ---- Assets (partition artifacts once) ----
    assets: dict[str, Any] = {}
    grouped: dict[str, list[tuple[str, Path]]] = {}
    for key, path in sorted(artifacts.items()):
        prefix_key = key.split(".")[0]
        grouped.setdefault(prefix_key, []).append((key, path))

    # BOA reflectance bands
    for key, path in grouped.get("boa", []):
        band_name = key.removeprefix("boa.")
        band_fields: dict[str, Any] = {}
        if sensor_config is not None:
            try:
                band_obj = sensor_config.get_band(band_name)
                band_fields["eo:bands"] = [_band_metadata(band_obj)]
            except (KeyError, AttributeError):
                band_fields["eo:bands"] = [_band_metadata_from_name(band_name)]
        else:
            band_fields["eo:bands"] = [_band_metadata_from_name(band_name)]
        assets[band_name] = _asset_dict(
            _relative_href(path, output_dir),
            title=f"BOA reflectance {band_name}",
            media_type=_media_type_for_artifact(path),
            roles=["data", "reflectance"],
            extra_fields=band_fields,
            file_size=_file_size(path),
        )

    # BOA uncertainty bands
    for key, path in grouped.get("boa_unc", []):
        band_name = key.removeprefix("boa_unc.")
        assets[f"{band_name}_unc"] = _asset_dict(
            _relative_href(path, output_dir),
            title=f"BOA uncertainty {band_name}",
            media_type=_media_type_for_artifact(path),
            roles=["data", "uncertainty"],
            file_size=_file_size(path),
        )

    # Auxiliary: AOT
    if "auxiliary.aot" in artifacts:
        assets["aot"] = _asset_dict(
            _relative_href(artifacts["auxiliary.aot"], output_dir),
            title="Aerosol optical thickness at 550 nm",
            media_type=_media_type_for_artifact(artifacts["auxiliary.aot"]),
            roles=["data", "metadata"],
            file_size=_file_size(artifacts["auxiliary.aot"]),
        )

    # Auxiliary: TCWV
    if "auxiliary.tcwv" in artifacts:
        assets["tcwv"] = _asset_dict(
            _relative_href(artifacts["auxiliary.tcwv"], output_dir),
            title="Total column water vapour (cm)",
            media_type=_media_type_for_artifact(artifacts["auxiliary.tcwv"]),
            roles=["data", "metadata"],
            file_size=_file_size(artifacts["auxiliary.tcwv"]),
        )

    # Auxiliary: cloud mask
    if "auxiliary.cloud_mask" in artifacts:
        assets["cloud-mask"] = _asset_dict(
            _relative_href(artifacts["auxiliary.cloud_mask"], output_dir),
            title="Cloud mask",
            media_type=_media_type_for_artifact(artifacts["auxiliary.cloud_mask"]),
            roles=["data", "cloud-mask"],
            file_size=_file_size(artifacts["auxiliary.cloud_mask"]),
        )

    # QA assets
    for key, path in grouped.get("auxiliary", []):
        if not key.startswith("auxiliary.qa."):
            continue
        qa_name = key.removeprefix("auxiliary.qa.")
        assets[f"qa_{qa_name}"] = _asset_dict(
            _relative_href(path, output_dir),
            title=f"QA: {qa_name.replace('_', ' ')}",
            media_type=_media_type_for_artifact(path),
            roles=["metadata", "quality"],
            file_size=_file_size(path),
        )

    # Surface prior
    for key, path in grouped.get("surface_prior", []):
        band_name = key.removeprefix("surface_prior.")
        assets[f"surf_{band_name}"] = _asset_dict(
            _relative_href(path, output_dir),
            title=f"Surface prior {band_name}",
            media_type=_media_type_for_artifact(path),
            roles=["data"],
            file_size=_file_size(path),
        )

    # RGB quicklook
    if "quicklook.rgb" in artifacts:
        assets["rendered"] = _asset_dict(
            _relative_href(artifacts["quicklook.rgb"], output_dir),
            title="RGB quicklook",
            media_type=_media_type_for_artifact(artifacts["quicklook.rgb"]),
            roles=["thumbnail", "visual"],
            file_size=_file_size(artifacts["quicklook.rgb"]),
        )

    # Preview PNGs
    _preview_specs: list[tuple[str, str, str, list[str]]] = [
        (
            "preview.false_colour",
            "preview_false_colour",
            "False colour composite (NIR-R-G)",
            ["overview", "visual"],
        ),
        ("preview.aot", "preview_aot", "AOT colour map", ["overview"]),
        ("preview.tcwv", "preview_tcwv", "TCWV colour map", ["overview"]),
        ("preview.cloud_mask", "preview_cloud_mask", "Cloud mask overlay", ["overview"]),
    ]
    for artifact_key, asset_id, title, roles in _preview_specs:
        if artifact_key in artifacts:
            assets[asset_id] = _asset_dict(
                _relative_href(artifacts[artifact_key], output_dir),
                title=title,
                media_type="image/png",
                roles=roles,
                file_size=_file_size(artifacts[artifact_key]),
            )

    # Scatter plot previews
    for key, path in grouped.get("preview", []):
        if key.startswith("preview.scatter."):
            band = key.removeprefix("preview.scatter.")
            assets[f"preview_scatter_{band}"] = _asset_dict(
                _relative_href(path, output_dir),
                title=f"AOT scatter plot {band}",
                media_type="image/png",
                roles=["overview"],
                file_size=_file_size(path),
            )

    # ---- Links ----
    links: list[dict[str, Any]] = [
        {
            "rel": "self",
            "href": "./",
            "type": "application/geo+json",
        }
    ]
    input_href = metadata.get("input_path")
    if input_href is not None:
        links.append({"rel": "derived_from", "href": str(input_href)})

    # ---- Extensions ----
    extensions = [_EO_EXTENSION, _PROJECTION_EXTENSION]
    if any(k.startswith("view:") for k in properties):
        extensions.append(_VIEW_EXTENSION)
    extensions.append(_PROCESSING_EXTENSION)

    return {
        "type": "Feature",
        "stac_version": _STAC_VERSION,
        "stac_extensions": extensions,
        "id": item_id,
        "bbox": bbox,
        "geometry": geometry_geojson,
        "properties": properties,
        "links": links,
        "assets": assets,
    }
