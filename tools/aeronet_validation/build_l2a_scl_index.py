"""Build a GEE-free Sentinel-2 winner index from Planetary Computer L2A SCL.

The production index was made from L1C/SR metadata in Earth Engine and a
Cloud Score+ confidence plane.  This experiment removes that dependency:

* candidate L2A items come from the Planetary Computer STAC API;
* the per-pixel clear decision comes only from the L2A SCL COG;
* the existing per-day MAIAC weight is retained when a reference index is
  available, so the experiment isolates the replacement of Cloud Score+;
* the output uses the normal ``months``/``winners``/``image_table``/
  ``day_scalars`` schema and can be consumed by either live L1C correction or
  the live L2A harmonizer;
* expanded L1C runs use measured MAIAC AOD where available and a spatially
  sampled daily CAMS value only for MAIAC retrieval gaps; there is no constant
  assumed-AOD fallback.

No Earth Engine module or credential is imported by this file.

The SCL mask has two modes:

``exact``
    Matches the existing clean-mask SCL exclusions ({1, 3, 8, 10}) while
    removing the Cloud Score+ requirement.
``standard``
    Keeps clear land/water/unclassified and snow/ice classes {4, 5, 6, 7,
    11}; snow is a valid bright surface target, while dark/shadow, cloud,
    cirrus and no-data classes remain excluded.

The default ``stac`` candidate mode is fully independent of the old winner
index.  ``reference`` is useful for a causal, same-candidate ablation: it
uses the old index's candidate acquisitions but reads their SCL from
Planetary Computer instead of reading Cloud Score+ from Earth Engine.

SCL replacement defaults to no erosion.  An optional erosion radius is
available for a conservative sensitivity run, but applying the historical
Cloud Score+ erosion to a categorical SCL mask is not the direct replacement.
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import math
import os
import re
import time
import zipfile
from collections import defaultdict
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path("/gws/ssde/j25a/nceo_isp/public/siac_refactor")
DEFAULT_CASES = ROOT / "reports/aod-final-performance-dashboard-20260713/data/all-cases.csv"
DEFAULT_DUMPS = ROOT / "calib_dumps_crossrt_dt_20260719"
DEFAULT_REFERENCE_INDEX = ROOT / "acix3_mosaic_index_proto"
DEFAULT_L2A_DICTIONARY = ROOT / "acix3_l2a_pc_dict"
DEFAULT_OUTPUT = ROOT / "analysis/l2a_scl_indices_pc_20260805"
DEFAULT_MAIAC_CACHE = ROOT / "cache/maiac_day_aod"
DEFAULT_CAMS = Path("/gws/ssde/j25b/nceo_ard/public/cams")

STAC_URL = "https://planetarycomputer.microsoft.com/api/stac/v1"
SAS_URL = "https://planetarycomputer.microsoft.com/api/sas/v1/token/sentinel-2-l2a"
LIBRARY_YEARS = range(2018, 2024)
LIBRARY_YEAR_MODES = ("fixed_2018_2023", "previous_5")
WINDOW_HALF_DAYS = 45
DEFAULT_WORKERS = 6

_SENSING_RE = re.compile(r"MSIL2A_(\d{8}T\d{6})")
_DATASTRIP_RE = re.compile(r"_S(\d{8}T\d{6})(?:_|$)")
_TOKEN_RE = re.compile(r"(\d{8}T\d{6})")

SCL_EXCLUDED_EXACT = {1, 3, 8, 10}
SCL_CLEAR = {
    # Match the historical SCL exclusions, but unlike the GEE-linked mask
    # also reject 0 because a standalone SCL raster has no Cloud Score+ mask
    # to suppress pixels outside the source footprint.
    "exact": set(range(1, 12)) - SCL_EXCLUDED_EXACT,
    "standard": {4, 5, 6, 7, 11},
}


def _valid_existing_index(
    path: Path, *, expected_policy: dict[str, Any] | None = None
) -> bool:
    """Return true only for a complete, readable winner-index archive."""

    if not path.is_file():
        return False
    try:
        with np.load(path, allow_pickle=False) as payload:
            required = {"months", "winners", "image_table", "day_scalars"}
            if not required.issubset(payload.files):
                return False
            months = np.asarray(payload["months"])
            winners = np.asarray(payload["winners"])
            if months.ndim != 1 or winners.ndim != 3:
                return False
            if winners.shape[0] != months.size or months.size < 1:
                return False
            # Force reads of the object-free string arrays so CRC/truncation
            # errors are detected before a previous output is reused.
            np.asarray(payload["image_table"])
            np.asarray(payload["day_scalars"])
            if expected_policy is not None:
                if "index_policy_json" not in payload.files:
                    return False
                observed = json.loads(str(np.asarray(payload["index_policy_json"]).item()))
                if observed != expected_policy:
                    return False
    except (OSError, ValueError, KeyError, EOFError, zipfile.BadZipFile):
        return False
    return True


def _imports() -> dict[str, Any]:
    """Import optional geospatial dependencies only in the worker process."""

    import requests
    from affine import Affine
    from osgeo import ogr, osr
    from pyproj import Transformer
    from rasterio.enums import Resampling
    from rasterio.vrt import WarpedVRT

    return {
        "Affine": Affine,
        "ogr": ogr,
        "osr": osr,
        "Resampling": Resampling,
        "Transformer": Transformer,
        "WarpedVRT": WarpedVRT,
        "requests": requests,
    }


def _session(requests: Any) -> Any:
    from requests.adapters import HTTPAdapter
    from urllib3.util.retry import Retry

    session = requests.Session()
    retry = Retry(
        total=7,
        connect=7,
        read=7,
        status=7,
        backoff_factor=1.5,
        backoff_max=45.0,
        backoff_jitter=2.0,
        # Planetary Computer can transiently answer anonymous STAC searches
        # with 403 when a large array fans out.  A persistent permission
        # error still fails after the bounded retry budget, while temporary
        # throttling no longer creates holes in an otherwise valid replay.
        status_forcelist=(403, 429, 500, 502, 503, 504),
        allowed_methods=frozenset({"GET", "POST"}),
        respect_retry_after_header=True,
    )
    session.mount("https://", HTTPAdapter(max_retries=retry, pool_connections=32, pool_maxsize=32))
    return session


def _mid_date(matchup_id: str) -> dt.date:
    match = _TOKEN_RE.search(matchup_id.rsplit("_", 1)[-1])
    if match is None:
        raise ValueError(f"Cannot parse observation date from {matchup_id!r}")
    return dt.datetime.strptime(match.group(1), "%Y%m%dT%H%M%S").date()


def _library_years(obs_date: dt.date, mode: str) -> tuple[int, ...]:
    """Resolve historical years without silently changing the frozen archive.

    ``fixed_2018_2023`` preserves the expanded-training archive.  ``previous_5``
    reproduces the production bestpixel contract used by the locked 2024
    benchmark: the five complete calendar years before the observation.
    """

    if mode == "fixed_2018_2023":
        return tuple(LIBRARY_YEARS)
    if mode == "previous_5":
        return tuple(range(obs_date.year - 5, obs_date.year))
    raise ValueError(f"Unsupported library-year mode {mode!r}")


def _windows(
    obs_date: dt.date,
    library_year_mode: str = "fixed_2018_2023",
) -> list[tuple[dt.date, dt.date]]:
    result = []
    for year in _library_years(obs_date, library_year_mode):
        middle = dt.date(year, obs_date.month, 15)
        result.append(
            (
                middle - dt.timedelta(days=WINDOW_HALF_DAYS),
                middle + dt.timedelta(days=WINDOW_HALF_DAYS),
            )
        )
    return result


def _load_dump(
    path: Path,
    imports: dict[str, Any],
    *,
    grid_name: str = "local60",
) -> dict[str, Any]:
    with np.load(path, allow_pickle=False) as payload:
        if "template_shape" in payload.files:
            if grid_name != "local60":
                raise ValueError(f"{path}: calibration dump has no {grid_name} grid")
            shape_key, transform_key, crs_key = (
                "template_shape",
                "template_transform",
                "template_crs",
            )
        elif f"{grid_name}_shape" in payload.files:
            shape_key, transform_key, crs_key = (
                f"{grid_name}_shape",
                f"{grid_name}_transform",
                "crs",
            )
        else:
            raise ValueError(
                f"{path} has neither a calibration template nor a {grid_name} mixed-resolution grid"
            )
        shape = tuple(int(value) for value in np.asarray(payload[shape_key]).ravel())
        transform_values = tuple(
            float(value) for value in np.asarray(payload[transform_key]).ravel()
        )
        crs = str(np.asarray(payload[crs_key]).item())
    # ``Affine`` serialises as either its six GDAL coefficients or a 3x3
    # homogeneous matrix.  The final row of the latter is always (0, 0, 1).
    transform_values = transform_values[:6]
    if len(shape) != 2 or len(transform_values) != 6:
        raise ValueError(f"Malformed calibration dump {path}")
    affine = imports["Affine"](*transform_values)
    height, width = shape
    # The calibration grid is north-up in this campaign.  The polygon is kept
    # in the native scene CRS for accurate overlap ratios and transformed to
    # WGS84 only for STAC's intersects query.
    corners = [
        (0.0, 0.0),
        (float(width), 0.0),
        (float(width), float(height)),
        (0.0, float(height)),
    ]
    xy = [(affine * point) for point in corners]
    scene_bounds = (
        min(point[0] for point in xy),
        min(point[1] for point in xy),
        max(point[0] for point in xy),
        max(point[1] for point in xy),
    )
    to_wgs84 = imports["Transformer"].from_crs(crs, "EPSG:4326", always_xy=True).transform
    wgs84_xy = [to_wgs84(*point) for point in xy]
    wgs84_xy.append(wgs84_xy[0])
    wgs84_polygon = {
        "type": "Polygon",
        "coordinates": [[list(point) for point in wgs84_xy]],
    }
    return {
        "height": height,
        "width": width,
        "transform": affine,
        "crs": crs,
        "scene_bounds": scene_bounds,
        "wgs84_polygon": wgs84_polygon,
    }


def _reference_payload(
    path: Path,
) -> tuple[list[str], np.ndarray, list[dict[str, Any]], list[dict[str, Any]]]:
    with np.load(path, allow_pickle=False) as payload:
        months = [str(value) for value in np.asarray(payload["months"]).ravel()]
        winners = np.asarray(payload["winners"], dtype=np.int16)
        table = [json.loads(str(value)) for value in np.asarray(payload["image_table"]).ravel()]
        scalars = [json.loads(str(value)) for value in np.asarray(payload["day_scalars"]).ravel()]
    return months, winners, table, scalars


def _reference_maps(path: Path) -> tuple[dict[str, dict[str, Any]], dict[str, dict[str, float]]]:
    """Return local angle and day-state maps from the old index artifact."""

    if not path.is_file():
        return {}, {}
    _, _, table, scalars = _reference_payload(path)
    images = {str(record.get("idx")): dict(record) for record in table if record.get("idx")}
    days: dict[str, dict[str, float]] = {}
    for record in scalars:
        day = str(record.get("day", ""))
        if not day:
            continue
        values: dict[str, float] = {}
        for key, value in record.items():
            try:
                number = float(value)
            except (TypeError, ValueError):
                continue
            if np.isfinite(number):
                values[str(key)] = number
        days[day] = values
    return images, days


def _load_wvp_by_month(path: Path, months: list[str]) -> dict[str, float]:
    if not path.is_file():
        return {}
    with np.load(path, allow_pickle=False) as payload:
        realizations = [str(value) for value in np.asarray(payload["realizations"]).ravel()]
        wvp = np.asarray(payload["wvp"], dtype=np.float64)
    values: dict[str, float] = {}
    for month, plane in zip(realizations, wvp, strict=True):
        finite = plane[np.isfinite(plane)]
        if finite.size:
            values[month] = float(np.nanmedian(finite))
    by_calendar: dict[str, list[float]] = defaultdict(list)
    for month, value in values.items():
        by_calendar[month[5:7]].append(value)
    fallback = float(np.median(list(values.values()))) if values else float("nan")
    for month in months:
        values.setdefault(month, float(np.median(by_calendar.get(month[5:7], [fallback]))))
    return values


def _numeric_baseline(value: Any) -> float:
    text = str(value or "").strip().upper().lstrip("N")
    if not text:
        return -1.0
    try:
        return float(int(text)) / 100.0 if len(text) == 4 and "." not in text else float(text)
    except ValueError:
        return -1.0


def _sensing_token(item: dict[str, Any]) -> str:
    match = _SENSING_RE.search(str(item.get("id", "")))
    if match:
        return match.group(1)
    value = str(item.get("properties", {}).get("datetime", ""))
    return dt.datetime.fromisoformat(value.replace("Z", "+00:00")).strftime("%Y%m%dT%H%M%S")


def _item_idx(item: dict[str, Any]) -> str:
    properties = item.get("properties") or {}
    sensing = _sensing_token(item)
    datastrip = str(properties.get("s2:datastrip_id") or "")
    match = _DATASTRIP_RE.search(datastrip)
    datastrip_token = match.group(1) if match else sensing
    tile = str(properties.get("s2:mgrs_tile") or "").upper()
    if tile and not tile.startswith("T"):
        tile = f"T{tile}"
    return f"{sensing}_{datastrip_token}_{tile}"


def _item_day(item: dict[str, Any]) -> str:
    token = _sensing_token(item)
    return f"{token[:4]}-{token[4:6]}-{token[6:8]}"


def _query_items(
    session: Any,
    wgs84_polygon: Any,
    obs_date: dt.date,
    library_year_mode: str = "fixed_2018_2023",
) -> list[dict[str, Any]]:
    items: dict[str, dict[str, Any]] = {}
    intersects = wgs84_polygon
    for start, end in _windows(obs_date, library_year_mode):
        payload = {
            "collections": ["sentinel-2-l2a"],
            "datetime": f"{start.isoformat()}T00:00:00Z/{end.isoformat()}T23:59:59Z",
            "intersects": intersects,
            "limit": 1000,
        }
        response = session.post(f"{STAC_URL}/search", json=payload, timeout=120)
        response.raise_for_status()
        features = response.json().get("features", [])
        if not isinstance(features, list):
            continue
        for feature in features:
            if not isinstance(feature, dict) or not feature.get("id"):
                continue
            assets = feature.get("assets") or {}
            if "SCL" not in assets:
                continue
            item_id = str(feature["id"])
            previous = items.get(item_id)
            if previous is None:
                items[item_id] = feature
    # Keep one reprocessing of each sensing/tile pair.  The highest baseline,
    # then newest generation, is the deterministic PC equivalent of selecting
    # the harmonized archive's current product.
    selected: dict[tuple[str, str], dict[str, Any]] = {}
    for item in items.values():
        p = item.get("properties") or {}
        key = (_sensing_token(item), str(p.get("s2:mgrs_tile") or ""))
        old = selected.get(key)
        rank = (
            _numeric_baseline(p.get("s2:processing_baseline")),
            str(p.get("s2:generation_time") or ""),
        )
        if old is None:
            selected[key] = item
            continue
        old_p = old.get("properties") or {}
        old_rank = (
            _numeric_baseline(old_p.get("s2:processing_baseline")),
            str(old_p.get("s2:generation_time") or ""),
        )
        if rank > old_rank:
            selected[key] = item
    return sorted(selected.values(), key=lambda item: (_item_day(item), _item_idx(item)))


def _reference_candidates(
    path: Path,
) -> tuple[list[dict[str, Any]], list[str], dict[str, dict[str, float]]]:
    if not path.is_file():
        return [], [], {}
    months, _, table, scalars = _reference_payload(path)
    return (
        table,
        months,
        {
            str(record["day"]): {
                str(key): float(value)
                for key, value in record.items()
                if key != "day" and isinstance(value, (int, float)) and np.isfinite(float(value))
            }
            for record in scalars
            if record.get("day")
        },
    )


def _reference_item_for_idx(items: list[dict[str, Any]]) -> dict[str, dict[str, Any]]:
    return {_item_idx(item): item for item in items}


def _angle_from_xml(session: Any, item: dict[str, Any], sas: str) -> tuple[float, float]:
    import xml.etree.ElementTree as ET

    href = str((item.get("assets") or {}).get("granule-metadata", {}).get("href") or "")
    if not href:
        raise RuntimeError(f"No granule metadata asset for {item.get('id')}")
    url = f"{href}&{sas}" if "?" in href else f"{href}?{sas}"
    response = session.get(url, timeout=120)
    response.raise_for_status()
    root = ET.fromstring(response.content)

    def local(tag: str) -> str:
        return tag.rsplit("}", 1)[-1]

    for element in root.iter():
        if local(element.tag) != "Mean_Viewing_Incidence_Angle":
            continue
        if str(element.attrib.get("bandId")) != "1":  # Sentinel-2 B02
            continue
        zenith = azimuth = None
        for child in element:
            name = local(child.tag)
            if name == "ZENITH_ANGLE":
                zenith = float(child.text)
            elif name == "AZIMUTH_ANGLE":
                azimuth = float(child.text)
        if zenith is not None and azimuth is not None:
            return zenith, azimuth
    raise RuntimeError(f"B02 mean viewing angle missing from {item.get('id')}")


def _record_from_item(
    item: dict[str, Any],
    *,
    ratio: float,
    reference_angles: dict[str, dict[str, Any]],
    session: Any,
    sas: str,
) -> dict[str, Any]:
    properties = item.get("properties") or {}
    idx = _item_idx(item)
    old = reference_angles.get(idx, {})
    if all(key in old for key in ("sza", "saa", "vza", "vaa")):
        sza, saa, vza, vaa = (float(old[key]) for key in ("sza", "saa", "vza", "vaa"))
    else:
        sza = float(properties.get("s2:mean_solar_zenith"))
        saa = float(properties.get("s2:mean_solar_azimuth"))
        vza, vaa = _angle_from_xml(session, item, sas)
    product = str(properties.get("s2:product_uri") or "") or None
    platform = str(properties.get("platform") or properties.get("s2:platform") or "") or None
    baseline = str(properties.get("s2:processing_baseline") or "") or None
    return {
        "idx": idx,
        "day": _item_day(item),
        "ratio": float(max(0.0, min(1.0, ratio))),
        "sza": sza,
        "saa": saa,
        "vza": vza,
        "vaa": vaa,
        "l2a_asset_id": str(item["id"]),
        "l2a_product_id": product,
        "l2a_processing_baseline": baseline,
        "l2a_spacecraft": platform,
    }


def _disk(radius: int) -> np.ndarray:
    yy, xx = np.ogrid[-radius : radius + 1, -radius : radius + 1]
    return (xx * xx + yy * yy) <= radius * radius


def _clear_distance_quality(
    mask: np.ndarray, *, pixel_size_m: float, radius_m: float
) -> np.ndarray:
    """Continuous [0,1] confidence from distance to the nearest SCL exclusion."""

    if pixel_size_m <= 0.0 or radius_m <= 0.0:
        raise ValueError("pixel_size_m and radius_m must be positive")
    from scipy.ndimage import distance_transform_edt

    distance_m = distance_transform_edt(np.asarray(mask, dtype=bool)) * float(pixel_size_m)
    return np.clip(distance_m / float(radius_m), 0.0, 1.0).astype(np.float32)


def _ocm_mask_and_quality(
    confidence: np.ndarray,
    source_valid: np.ndarray,
    *,
    thin_quality_weight: float,
    erosion: np.ndarray | None,
) -> tuple[np.ndarray, np.ndarray]:
    """Retain OCM clear/thin support and rank clear pixels above thin ones."""

    probabilities = np.asarray(confidence, dtype=np.float32)
    if probabilities.ndim != 3 or probabilities.shape[0] != 4:
        raise ValueError(f"Expected four OCM confidence planes, got {probabilities.shape}")
    valid = np.asarray(source_valid, dtype=bool)
    if probabilities.shape[1:] != valid.shape:
        raise ValueError(
            f"OCM/source grid mismatch: {probabilities.shape[1:]} vs {valid.shape}"
        )
    if not 0.0 <= thin_quality_weight <= 1.0:
        raise ValueError("thin_quality_weight must be in [0, 1]")

    # Native OCM classes: 0 clear, 1 thick cloud, 2 thin cloud, 3 shadow.
    classes = np.argmax(probabilities, axis=0)
    valid &= (classes != 1) & (classes != 3)
    if erosion is not None:
        from scipy.ndimage import binary_erosion

        valid = binary_erosion(valid, structure=erosion, border_value=0)
    quality = probabilities[0] + float(thin_quality_weight) * probabilities[2]
    quality = np.clip(quality, 0.0, 1.0)
    quality[~valid] = 0.0
    return valid, quality.astype(np.float32, copy=False)


def _read_ocm_input(
    item: dict[str, Any], *, sas: str, grid: dict[str, Any]
) -> np.ndarray:
    """Read L2A red/green/NIR COG windows on the winner-index grid."""

    import rasterio

    bands: list[np.ndarray] = []
    for band_name in ("B04", "B03", "B08"):
        href = str(item["assets"][band_name]["href"])
        url = f"{href}&{sas}" if "?" in href else f"{href}?{sas}"
        with (
            rasterio.Env(
                GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR",
                GDAL_HTTP_MAX_RETRY="3",
                GDAL_HTTP_RETRY_DELAY="2",
            ),
            rasterio.open(url) as source,
            grid["WarpedVRT"](
                source,
                crs=grid["crs"],
                transform=grid["transform"],
                width=grid["width"],
                height=grid["height"],
                resampling=grid["Resampling"].bilinear,
                nodata=0,
            ) as vrt,
        ):
            bands.append(vrt.read(1).astype(np.float32, copy=False))
    return np.stack(bands)


def _load_ocm_models() -> list[Any]:
    """Load the deterministic CPU ensemble once per matchup process."""

    import torch
    from omnicloudmask.cloud_mask import collect_models

    return collect_models(
        custom_models=None,
        inference_device=torch.device("cpu"),
        inference_dtype=torch.float32,
        source="hugging_face",
    )


def _predict_ocm_confidence(values: np.ndarray, models: list[Any]) -> np.ndarray:
    import omnicloudmask

    return np.asarray(
        omnicloudmask.predict_from_array(
            np.asarray(values, dtype=np.float32),
            inference_device="cpu",
            mosaic_device="cpu",
            export_confidence=True,
            softmax_output=True,
            custom_models=models,
        ),
        dtype=np.float32,
    )


def _read_scl(
    item: dict[str, Any],
    *,
    sas: str,
    grid: dict[str, Any],
    clear_classes: set[int],
    erosion: np.ndarray | None,
) -> np.ndarray:
    import rasterio
    from scipy.ndimage import binary_erosion

    href = str(item["assets"]["SCL"]["href"])
    url = f"{href}&{sas}" if "?" in href else f"{href}?{sas}"
    with (
        rasterio.Env(
            GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR",
            GDAL_HTTP_MAX_RETRY="3",
            GDAL_HTTP_RETRY_DELAY="2",
        ),
        rasterio.open(url) as source,
        grid["WarpedVRT"](
            source,
            crs=grid["crs"],
            transform=grid["transform"],
            width=grid["width"],
            height=grid["height"],
            resampling=grid["Resampling"].nearest,
            nodata=0,
        ) as vrt,
    ):
        scl = vrt.read(1)
    clean = np.isin(scl, list(clear_classes))
    if erosion is None:
        return clean
    # An optional focal-min is applied after warping to the 60 m target grid.
    # The default direct SCL replacement does not use it.
    return binary_erosion(clean, structure=erosion, border_value=0)


def _coverage_ratio(
    item: dict[str, Any],
    scene_bounds: tuple[float, float, float, float],
    crs: str,
    imports: dict[str, Any],
) -> float:
    geometry = item.get("geometry")
    if not geometry:
        return 0.0
    try:
        footprint = imports["ogr"].CreateGeometryFromJson(json.dumps(geometry))
        source_srs = imports["osr"].SpatialReference()
        source_srs.ImportFromEPSG(4326)
        target_srs = imports["osr"].SpatialReference()
        target_srs.SetFromUserInput(crs)
        source_srs.SetAxisMappingStrategy(imports["osr"].OAMS_TRADITIONAL_GIS_ORDER)
        target_srs.SetAxisMappingStrategy(imports["osr"].OAMS_TRADITIONAL_GIS_ORDER)
        # CoordinateTransformation lives in osr, not ogr. Calling it on ogr raises
        # AttributeError, which the guard below swallows into a 0.0 overlap for
        # every acquisition, silently zeroing the aoi_overlap_ratio ordering term.
        footprint.Transform(imports["osr"].CoordinateTransformation(source_srs, target_srs))
        west, south, east, north = scene_bounds
        ring = imports["ogr"].Geometry(imports["ogr"].wkbLinearRing)
        ring.AddPoint(west, south, 0.0)
        ring.AddPoint(east, south, 0.0)
        ring.AddPoint(east, north, 0.0)
        ring.AddPoint(west, north, 0.0)
        ring.AddPoint(west, south, 0.0)
        scene = imports["ogr"].Geometry(imports["ogr"].wkbPolygon)
        scene.AddGeometry(ring)
        area = scene.GetArea()
        intersection = footprint.Intersection(scene)
        return float(intersection.GetArea() / area) if area > 0 else 0.0
    except Exception:  # noqa: BLE001 - an invalid footprint cannot win pixels
        return 0.0


def _day_weight(day_state: dict[str, float]) -> float:
    value = day_state.get("weight")
    return float(value) if value is not None and np.isfinite(value) else 0.0


def _aod_selection_weights(days: list[str], aod_by_day: dict[str, float]) -> np.ndarray:
    """Return month-local low-AOD scores in [0, 1].

    MAIAC values are already in physical AOD units.  Ranking within each
    calendar window is invariant to product scale and gives the cleanest day
    the largest score.  Unknown days remain NaN and must not win a pixel when
    the measured-AOD path is requested.
    """

    weights = np.full(len(days), np.nan, dtype=np.float64)
    by_month: dict[str, list[int]] = defaultdict(list)
    for index, day in enumerate(days):
        value = aod_by_day.get(day)
        if value is not None and np.isfinite(value):
            by_month[day[:7]].append(index)
    for indices in by_month.values():
        ordered = sorted(indices, key=lambda index: (float(aod_by_day[days[index]]), days[index]))
        if len(ordered) == 1:
            weights[ordered[0]] = 1.0
            continue
        denominator = float(len(ordered) - 1)
        for rank, index in enumerate(ordered):
            weights[index] = 1.0 - rank / denominator
    return weights


def _locked_aod_selection_weights(
    days: list[str],
    aod_by_day: dict[str, float],
    *,
    threshold_percentile: float = 60.0,
) -> np.ndarray:
    """Reproduce the locked Cloud-Score+ library's raw-MAIAC preference.

    The committed reference calculated one sigmoid over every historical day,
    using the unscaled ``Optical_Depth_055`` integers and clipping at the 60th
    percentile.  Public MAIAC adapters return physical AOD, so multiply by
    1000 before applying the frozen expression.  Missing-MAIAC days remain
    NaN here and receive zero AOD weight when the complete day score is built,
    exactly as in the locked implementation.
    """

    physical = np.asarray([aod_by_day.get(day, np.nan) for day in days], dtype=np.float64)
    finite = np.isfinite(physical)
    weights = np.full(physical.shape, np.nan, dtype=np.float64)
    if not np.any(finite):
        return weights

    raw = physical * 1000.0
    aod_min = float(np.nanmin(raw))
    aod_max = float(np.nanpercentile(raw, threshold_percentile))
    if np.isclose(aod_max, aod_min):
        weights[finite] = 1.0
        return weights

    exponent = -0.2 * (np.minimum(raw[finite], aod_max) - (aod_max - aod_min) / 2.0)
    exponent = np.clip(exponent, -700.0, 700.0)
    weights[finite] = 1.0 - 1.0 / (1.0 + np.exp(exponent))
    return weights


def _candidate_mask(
    coverage: np.ndarray,
    months_by_day: np.ndarray,
    *,
    eligible: np.ndarray | None = None,
) -> np.ndarray:
    """Choose above-mean positive-coverage days within the eligible pool."""

    values = np.asarray(coverage, dtype=np.float64)
    months = np.asarray(months_by_day)
    allowed = (
        np.ones(values.shape, dtype=bool) if eligible is None else np.asarray(eligible, dtype=bool)
    )
    if values.shape != months.shape or values.shape != allowed.shape:
        raise ValueError("coverage, month and eligibility arrays must align")
    candidate = np.zeros(values.shape, dtype=bool)
    for month in np.unique(months[allowed]):
        pool = (months == month) & allowed
        mean = float(np.nanmean(values[pool]))
        candidate[pool] = (values[pool] >= mean) & (values[pool] > 0)
    return candidate


def _resolve_day_aod(
    *,
    source: str,
    days: list[str],
    reference_scalars: dict[str, dict[str, float]],
    scene_bounds: tuple[float, float, float, float],
    crs: str,
    cache_dir: Path,
) -> dict[str, float]:
    if source == "none":
        return {}
    if source == "reference":
        return {
            day: float(reference_scalars[day]["aod"])
            for day in days
            if day in reference_scalars
            and "aod" in reference_scalars[day]
            and np.isfinite(reference_scalars[day]["aod"])
        }
    if source != "maiac":
        raise ValueError(f"Unsupported day-AOD source {source!r}")
    from siac.adapters.atmo.maiac_day_aod import MAIACDayAODProvider

    periods = sorted({(int(day[:4]), int(day[5:7])) for day in days})
    measured = MAIACDayAODProvider(cache_dir=cache_dir).day_aod_map(
        scene_bounds,
        crs,
        periods,
    )
    return {
        day: float(measured[day]) for day in days if day in measured and np.isfinite(measured[day])
    }


def _cams_day_aod(
    days: list[str],
    *,
    scene_bounds: tuple[float, float, float, float],
    crs: str,
    cams_dir: Path | None,
) -> dict[str, float]:
    """Read daily median CAMS AOD only where the MAIAC day gate is unavailable."""

    if cams_dir is None:
        return {}
    import rasterio
    from pyproj import Transformer
    from rasterio.windows import Window

    west, south, east, north = scene_bounds
    lon, lat = Transformer.from_crs(crs, "EPSG:4326", always_xy=True).transform(
        0.5 * (west + east), 0.5 * (south + north)
    )
    values: dict[str, float] = {}
    for day in dict.fromkeys(days):
        token = day.replace("-", "_")
        path = cams_dir / token / f"{token}_aod550.tif"
        if not path.is_file():
            continue
        try:
            with rasterio.open(path) as source:
                row, column = source.index(lon, lat)
                row0 = max(0, min(source.height - 1, row) - 1)
                column0 = max(0, min(source.width - 1, column) - 1)
                window = Window(
                    column0,
                    row0,
                    min(3, source.width - column0),
                    min(3, source.height - row0),
                )
                raw = source.read(window=window, masked=True).astype(np.float64)
                scales = np.asarray(source.scales, dtype=np.float64)[:, None, None]
                offsets = np.asarray(source.offsets, dtype=np.float64)[:, None, None]
                physical = raw * scales + offsets
                finite = np.asarray(physical.compressed(), dtype=np.float64)
        except (OSError, ValueError, IndexError):
            continue
        finite = finite[np.isfinite(finite) & (finite >= 0.0) & (finite <= 5.0)]
        if finite.size:
            values[day] = float(np.median(finite))
    return values


def _constant_aod_for_unresolved_days(
    days: list[str],
    *resolved_sources: dict[str, float],
    value: float = 0.1,
) -> dict[str, float]:
    """Guarantee an RT state for candidate days unresolved by measured sources."""

    resolved = {day for source in resolved_sources for day in source}
    return {day: float(value) for day in dict.fromkeys(days) if day not in resolved}


def run_one(
    matchup_id: str,
    *,
    dumps_dir: Path,
    reference_dir: Path,
    l2a_dir: Path,
    output_dir: Path,
    candidate_source: str,
    scl_mode: str,
    erosion_radius_m: float,
    workers: int,
    runtime_source: str,
    day_aod_source: str,
    maiac_cache: Path,
    cams_dir: Path | None,
    force: bool,
    dump_grid: str = "local60",
    library_year_mode: str = "fixed_2018_2023",
    distance_quality_radius_m: float = 0.0,
    aod_quality_mode: str = "month_rank",
    ocm_confidence_quality: bool = False,
    ocm_thin_quality_weight: float = 0.25,
    maiac_gap_policy: str = "cams_fallback",
) -> dict[str, Any]:
    imports = _imports()
    target = output_dir / f"{matchup_id}.npz"
    summary_path = output_dir / f"{matchup_id}.json"
    index_policy = {
        "aod_quality_mode": str(aod_quality_mode),
        "candidate_source": str(candidate_source),
        "day_aod_source": str(day_aod_source),
        "distance_quality_radius_m": float(distance_quality_radius_m),
        "dump_grid": str(dump_grid),
        "erosion_radius_m": float(erosion_radius_m),
        "library_year_mode": str(library_year_mode),
        "maiac_gap_policy": str(maiac_gap_policy),
        "ocm_confidence_quality": bool(ocm_confidence_quality),
        "ocm_thin_quality_weight": float(ocm_thin_quality_weight),
        "runtime_source": str(runtime_source),
        "scl_mode": str(scl_mode),
    }
    if not force and _valid_existing_index(target, expected_policy=index_policy):
        return {"matchup_id": matchup_id, "status": "exists", "output": str(target)}

    dump = _load_dump(
        dumps_dir / f"{matchup_id}.npz",
        imports,
        grid_name=dump_grid,
    )
    if ocm_confidence_quality and min(int(dump["height"]), int(dump["width"])) < 200:
        raise ValueError(
            "OmniCloudMask requires a grid at least 200 pixels on each side; "
            "use --dump-grid detail20"
        )
    reference_path = reference_dir / f"{matchup_id}.npz"
    reference_images, reference_days = _reference_maps(reference_path)
    reference_table, reference_months, reference_scalar_days = _reference_candidates(reference_path)
    if not reference_days:
        reference_days = reference_images
    session = _session(imports["requests"])
    sas_response = session.get(SAS_URL, timeout=60)
    sas_response.raise_for_status()
    sas = str(sas_response.json()["token"])

    t0 = time.perf_counter()
    observation_date = _mid_date(matchup_id)
    stac_items = _query_items(
        session,
        dump["wgs84_polygon"],
        observation_date,
        library_year_mode,
    )
    by_idx = _reference_item_for_idx(stac_items)
    if candidate_source == "reference":
        wanted = {str(record.get("idx")) for record in reference_table if record.get("idx")}
        items = [by_idx[idx] for idx in sorted(wanted) if idx in by_idx]
    else:
        items = stac_items
    if not items:
        raise RuntimeError(f"No Planetary Computer L2A SCL items for {matchup_id}")

    # The old index is only an angle/day-weight cache in fully-STAC mode.  A
    # matching local record avoids an XML request; otherwise B02 angles are
    # obtained directly from L2A tile metadata, still without GEE.
    records_by_id: dict[str, dict[str, Any]] = {}
    item_by_id: dict[str, dict[str, Any]] = {}
    for item in items:
        item_id = str(item["id"])
        item_by_id[item_id] = item
        idx = _item_idx(item)
        records_by_id[item_id] = {
            "idx": idx,
            "day": _item_day(item),
            "ratio": _coverage_ratio(item, dump["scene_bounds"], dump["crs"], imports),
        }
    t_stac = time.perf_counter() - t0

    clear_classes = set(SCL_CLEAR[scl_mode])
    if ocm_confidence_quality:
        # In this mode SCL only supplies source validity/saturation. OCM is the
        # cloud classifier, so SCL cloud labels must not remove thin support.
        clear_classes = set(range(2, 12))
    if erosion_radius_m > 0.0:
        radius = max(1, int(math.ceil(erosion_radius_m / abs(float(dump["transform"].a)))))
        erosion: np.ndarray | None = _disk(radius)
    else:
        erosion = None
    grid = {
        **dump,
        "WarpedVRT": imports["WarpedVRT"],
        "Resampling": imports["Resampling"],
    }
    t0 = time.perf_counter()
    masks: dict[str, np.ndarray] = {}
    ocm_inputs: dict[str, np.ndarray] = {}
    failures: dict[str, str] = {}

    def fetch(
        item: dict[str, Any],
    ) -> tuple[str, np.ndarray | None, np.ndarray | None, str | None]:
        item_id = str(item["id"])
        try:
            return (
                item_id,
                _read_scl(
                    item,
                    sas=sas,
                    grid=grid,
                    clear_classes=clear_classes,
                    erosion=None if ocm_confidence_quality else erosion,
                ),
                _read_ocm_input(item, sas=sas, grid=grid)
                if ocm_confidence_quality
                else None,
                None,
            )
        except Exception as exc:  # noqa: BLE001 - retain other acquisitions
            return item_id, None, None, f"{type(exc).__name__}: {exc}"

    from concurrent.futures import ThreadPoolExecutor

    with ThreadPoolExecutor(max_workers=max(1, int(workers))) as pool:
        for item_id, mask, ocm_input, error in pool.map(fetch, items):
            if mask is not None:
                masks[item_id] = mask
                if ocm_input is not None:
                    ocm_inputs[item_id] = ocm_input
            elif error:
                failures[item_id] = error
    t_scl = time.perf_counter() - t0

    ocm_quality: dict[str, np.ndarray] = {}
    t_ocm = 0.0
    if ocm_confidence_quality and masks:
        t0 = time.perf_counter()
        models = _load_ocm_models()
        for item_id in list(masks):
            try:
                confidence = _predict_ocm_confidence(ocm_inputs[item_id], models)
                masks[item_id], ocm_quality[item_id] = _ocm_mask_and_quality(
                    confidence,
                    masks[item_id],
                    thin_quality_weight=float(ocm_thin_quality_weight),
                    erosion=erosion,
                )
            except Exception as exc:  # noqa: BLE001 - retain other acquisitions
                failures[item_id] = f"OCM {type(exc).__name__}: {exc}"
                masks.pop(item_id, None)
                ocm_quality.pop(item_id, None)
        t_ocm = time.perf_counter() - t0

    usable_items = [item for item in items if str(item["id"]) in masks]
    if not usable_items:
        raise RuntimeError(
            f"All SCL reads failed for {matchup_id}: {next(iter(failures.values()), 'unknown')}"
        )

    # Cloud Score+ supplied a continuous per-pixel confidence to qualityMosaic.
    # A binary SCL mask cannot reproduce that ordering: every valid pixel ties.
    # This optional, GEE-free proxy rewards distance from the nearest SCL
    # exclusion while retaining valid near-edge pixels when no better date is
    # available.  It is deliberately separate from hard binary erosion.
    distance_quality: dict[str, np.ndarray] = {}
    if distance_quality_radius_m > 0.0:
        pixel_size_m = max(abs(float(dump["transform"].a)), abs(float(dump["transform"].e)))
        for item_id, mask in masks.items():
            distance_quality[item_id] = _clear_distance_quality(
                mask,
                pixel_size_m=pixel_size_m,
                radius_m=float(distance_quality_radius_m),
            )

    by_day: dict[str, list[str]] = defaultdict(list)
    for item in usable_items:
        by_day[_item_day(item)].append(str(item["id"]))
    days = sorted(by_day)
    coverage = np.zeros(len(days), dtype=np.float32)
    for i, day in enumerate(days):
        stack = np.stack([masks[item_id] for item_id in by_day[day]], axis=0)
        # This is the local equivalent of the original same-day median mask.
        coverage[i] = float((np.mean(stack, axis=0) >= 0.5).mean())
    months_by_day = np.asarray([day[:7] for day in days])
    day_index = {day: i for i, day in enumerate(days)}
    maiac_aod_by_day = _resolve_day_aod(
        source=day_aod_source,
        days=days,
        reference_scalars=reference_scalar_days,
        scene_bounds=dump["scene_bounds"],
        crs=dump["crs"],
        cache_dir=maiac_cache,
    )
    cams_aod_by_day: dict[str, float] = {}
    constant_aod_by_day: dict[str, float] = {}
    if day_aod_source == "maiac":
        missing_days = [day for day in days if day not in maiac_aod_by_day]
        if maiac_gap_policy == "cams_fallback":
            cams_aod_by_day = _cams_day_aod(
                missing_days,
                scene_bounds=dump["scene_bounds"],
                crs=dump["crs"],
                cams_dir=cams_dir,
            )
            # The locked GEE recipe permits missing-MAIAC days to enter the
            # quality mosaic and corrects those winners at AOD=0.1. CAMS is
            # the preferred GEE-free fill, but its local archive can also have
            # a gap; never publish a winning day without a physical RT state.
            if aod_quality_mode == "locked_raw_sigmoid":
                constant_aod_by_day = _constant_aod_for_unresolved_days(
                    missing_days,
                    cams_aod_by_day,
                )
        elif maiac_gap_policy == "locked_constant_0p1":
            constant_aod_by_day = dict.fromkeys(missing_days, 0.1)
        else:
            raise ValueError(f"Unsupported MAIAC gap policy {maiac_gap_policy!r}")
    aod_by_day = {**constant_aod_by_day, **cams_aod_by_day, **maiac_aod_by_day}
    if day_aod_source == "maiac":
        if aod_quality_mode == "locked_raw_sigmoid":
            # Match the locked GEE recipe: cloud coverage defines candidate
            # days, missing MAIAC contributes zero to the final score, and the
            # AOD preference is global rather than month-local.  CAMS remains
            # available as the historical RT fallback but must not alter the
            # MAIAC-only image-quality ranking.
            candidate = _candidate_mask(coverage, months_by_day)
            day_weights = _locked_aod_selection_weights(days, maiac_aod_by_day)
            day_weights[~np.isfinite(day_weights)] = 0.0
        elif aod_quality_mode == "month_rank":
            # Prefer measured MAIAC; daily CAMS fills only MAIAC retrieval gaps
            # (notably snow/ice), so every eligible winner still has an RT state.
            measured = np.asarray([day in aod_by_day for day in days], dtype=bool)
            candidate = _candidate_mask(coverage, months_by_day, eligible=measured)
            day_weights = _aod_selection_weights(days, aod_by_day)
        else:
            raise ValueError(f"Unsupported AOD quality mode {aod_quality_mode!r}")
    else:
        candidate = _candidate_mask(coverage, months_by_day)
        # Preserve the reference ranking for controlled legacy/SCL ablations.
        day_weights = np.asarray(
            [_day_weight(reference_scalar_days.get(day, {})) for day in days],
            dtype=np.float64,
        )
        day_weights[~np.isfinite(day_weights)] = 0.0
    coverage_weights = 1.0 / (1.0 + np.exp(-20.0 * (coverage.astype(np.float64) - 0.5)))

    t0 = time.perf_counter()
    candidate_records: dict[str, dict[str, Any]] = {}
    winners: list[np.ndarray] = []
    output_months: list[str] = []
    for month in sorted(set(months_by_day[candidate])):
        month_days = [day for day in days if day[:7] == month and candidate[day_index[day]]]
        month_items = [item for day in month_days for item in by_day[day]]
        month_items.sort(
            key=lambda item_id: (
                str(item_by_id[item_id].get("properties", {}).get("datetime", "")),
                item_id,
            )
        )
        score = np.full(
            (len(month_items), dump["height"], dump["width"]), -np.inf, dtype=np.float32
        )
        for row, item_id in enumerate(month_items):
            item = item_by_id[item_id]
            day = _item_day(item)
            score[row] = np.where(
                masks[item_id],
                day_weights[day_index[day]]
                + coverage_weights[day_index[day]]
                + float(records_by_id[item_id]["ratio"])
                + ocm_quality.get(item_id, distance_quality.get(item_id, 0.0)),
                -np.inf,
            )
        win = np.argmax(score, axis=0).astype(np.int16)
        win[~np.isfinite(np.max(score, axis=0))] = -1
        winners.append(win)
        output_months.append(month)
        for item_id in month_items:
            candidate_records[item_id] = records_by_id[item_id]

    if not winners:
        raise RuntimeError(f"No SCL candidate month for {matchup_id}; coverage={coverage.tolist()}")

    # Resolve angles only for acquisitions that can actually win pixels.  The
    # local reference cache covers the normal case; XML is the GEE-free
    # fallback for newly discovered STAC acquisitions.
    output_records: list[dict[str, Any]] = []
    for month in output_months:
        month_ids = [
            item_id for item_id, record in candidate_records.items() if record["day"][:7] == month
        ]
        month_ids.sort(key=lambda item_id: (candidate_records[item_id]["day"], item_id))
        for item_id in month_ids:
            record = _record_from_item(
                item_by_id[item_id],
                ratio=float(candidate_records[item_id]["ratio"]),
                reference_angles=reference_images,
                session=session,
                sas=sas,
            )
            output_records.append(record)

    # L1C runtime reads the exact winning L2A WVP raster and needs no prepared
    # L2A dictionary.  The legacy L2A runtime keeps its monthly WVP fallback.
    wvp_by_month = (
        _load_wvp_by_month(l2a_dir / f"{matchup_id}.npz", output_months)
        if runtime_source == "l2a"
        else {}
    )
    day_scalars: list[dict[str, Any]] = []
    for day in days:
        state: dict[str, Any] = (
            dict(reference_scalar_days.get(day, {})) if day_aod_source == "reference" else {}
        )
        state["day"] = day
        if day in aod_by_day:
            state["aod"] = float(aod_by_day[day])
            if day_aod_source == "maiac":
                if day in maiac_aod_by_day:
                    state["aod_source"] = "maiac"
                elif day in cams_aod_by_day:
                    state["aod_source"] = "cams_fallback"
                else:
                    state["aod_source"] = "locked_constant_0p1"
        weight = day_weights[day_index[day]]
        if np.isfinite(weight):
            state["weight"] = float(weight)
        if runtime_source == "l2a":
            state["maiac_tcwv_cm"] = float(wvp_by_month.get(day[:7], np.nan))
        day_scalars.append(state)

    output_dir.mkdir(parents=True, exist_ok=True)
    temporary = target.with_name(f".{target.stem}.{os.getpid()}.partial.npz")
    try:
        np.savez_compressed(
            temporary,
            months=np.asarray(output_months),
            winners=np.stack(winners),
            image_table=np.asarray(
                [json.dumps(record, sort_keys=True) for record in output_records]
            ),
            day_scalars=np.asarray([json.dumps(state, sort_keys=True) for state in day_scalars]),
            index_policy_json=np.asarray(json.dumps(index_policy, sort_keys=True)),
        )
        if not _valid_existing_index(temporary, expected_policy=index_policy):
            raise RuntimeError(f"Generated invalid winner-index archive {temporary}")
        os.replace(temporary, target)
    finally:
        temporary.unlink(missing_ok=True)
    summary = {
        "matchup_id": matchup_id,
        "status": "ok",
        "candidate_source": candidate_source,
        "runtime_source": runtime_source,
        "dump_grid": dump_grid,
        "day_aod_source": day_aod_source,
        "day_aod_policy": (
            (
                (
                    "maiac_with_daily_cams_then_constant_0p1_gap_fill"
                    if aod_quality_mode == "locked_raw_sigmoid"
                    else "maiac_with_daily_cams_gap_fill"
                )
                if maiac_gap_policy == "cams_fallback"
                else "maiac_with_locked_constant_0p1_gap_fill"
            )
            if day_aod_source == "maiac"
            else day_aod_source
        ),
        "aod_quality_mode": aod_quality_mode,
        "maiac_gap_policy": maiac_gap_policy,
        "library_year_mode": library_year_mode,
        "library_years": list(_library_years(observation_date, library_year_mode)),
        "days_with_aod": len(aod_by_day),
        "days_with_maiac_aod": len(maiac_aod_by_day),
        "days_with_cams_fallback_aod": len(cams_aod_by_day),
        "days_with_constant_0p1_aod": len(constant_aod_by_day),
        "scl_mode": scl_mode,
        "erosion_radius_m": float(erosion_radius_m),
        "distance_quality_radius_m": float(distance_quality_radius_m),
        "ocm_confidence_quality": bool(ocm_confidence_quality),
        "ocm_thin_quality_weight": float(ocm_thin_quality_weight),
        "stac_items": len(stac_items),
        "items_used_for_scl": len(items),
        "scl_reads_ok": len(masks),
        "scl_reads_failed": len(failures),
        "days": len(days),
        "candidate_days": int(candidate.sum()),
        "months": len(output_months),
        "winning_images": len(output_records),
        "clear_coverage_median": float(np.median(coverage)),
        "clear_coverage_mean": float(np.mean(coverage)),
        "t_stac_s": round(t_stac, 2),
        "t_scl_s": round(t_scl, 2),
        "t_ocm_s": round(t_ocm, 2),
        "t_index_s": round(time.perf_counter() - t0, 2),
        "output": str(target),
    }
    summary_temporary = summary_path.with_name(f".{summary_path.stem}.{os.getpid()}.partial.json")
    summary_temporary.write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    os.replace(summary_temporary, summary_path)
    print(json.dumps(summary), flush=True)
    return summary


def _matchup_ids(cases_path: Path) -> list[str]:
    lines = cases_path.read_text(encoding="utf-8").splitlines()
    return [line.split(",", 1)[0].strip() for line in lines[1:] if line.strip()]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("matchup_ids", nargs="*")
    parser.add_argument("--cases", type=Path, default=DEFAULT_CASES)
    parser.add_argument("--dumps", type=Path, default=DEFAULT_DUMPS)
    parser.add_argument("--reference-index", type=Path, default=DEFAULT_REFERENCE_INDEX)
    parser.add_argument("--l2a-dictionary", type=Path, default=DEFAULT_L2A_DICTIONARY)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--candidate-source", choices=("stac", "reference"), default="stac")
    parser.add_argument("--runtime-source", choices=("l1c", "l2a"), default="l2a")
    parser.add_argument(
        "--day-aod-source",
        choices=("maiac", "reference", "none"),
        default="reference",
    )
    parser.add_argument(
        "--aod-quality-mode",
        choices=("month_rank", "locked_raw_sigmoid"),
        default="month_rank",
        help=(
            "Historical image AOD preference: the expanded month-local rank or "
            "the locked Cloud-Score+ raw-MAIAC sigmoid contract."
        ),
    )
    parser.add_argument(
        "--maiac-gap-policy",
        choices=("cams_fallback", "locked_constant_0p1"),
        default="cams_fallback",
        help=(
            "AOD used to correct a historical winner when same-day MAIAC is "
            "missing: measured daily CAMS or the frozen reference's 0.1 constant."
        ),
    )
    parser.add_argument("--maiac-cache", type=Path, default=DEFAULT_MAIAC_CACHE)
    parser.add_argument("--cams", type=Path, default=DEFAULT_CAMS)
    parser.add_argument("--scl-mode", choices=("exact", "standard"), default="standard")
    parser.add_argument("--erosion-radius-m", type=float, default=0.0)
    parser.add_argument(
        "--distance-quality-radius-m",
        type=float,
        default=0.0,
        help=(
            "Optional continuous SCL clear-distance score, capped at this radius; 0 disables it."
        ),
    )
    parser.add_argument(
        "--ocm-confidence-quality",
        action="store_true",
        help=(
            "Use deterministic CPU OmniCloudMask confidence for monthly winner "
            "quality; thick cloud and shadow are excluded while thin cloud remains eligible."
        ),
    )
    parser.add_argument(
        "--ocm-thin-quality-weight",
        type=float,
        default=0.25,
        help="Relative OCM thin-cloud confidence in winner ranking (0-1; default 0.25).",
    )
    parser.add_argument(
        "--library-year-mode",
        choices=LIBRARY_YEAR_MODES,
        default="fixed_2018_2023",
        help=(
            "Historical calendar years: preserve the expanded fixed archive or "
            "use the production five complete years before each observation."
        ),
    )
    parser.add_argument("--workers", type=int, default=DEFAULT_WORKERS)
    parser.add_argument(
        "--dump-grid",
        choices=("local60", "detail20"),
        default="local60",
        help="Grid to read from mixed-resolution dump archives.",
    )
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    if args.distance_quality_radius_m < 0.0:
        parser.error("--distance-quality-radius-m must be >= 0")
    if args.ocm_confidence_quality and args.distance_quality_radius_m > 0.0:
        parser.error(
            "--ocm-confidence-quality and --distance-quality-radius-m are alternative quality scores"
        )
    if not 0.0 <= args.ocm_thin_quality_weight <= 1.0:
        parser.error("--ocm-thin-quality-weight must be in [0, 1]")
    matchup_ids = args.matchup_ids or _matchup_ids(args.cases)
    for matchup_id in matchup_ids:
        try:
            run_one(
                matchup_id,
                dumps_dir=args.dumps,
                reference_dir=args.reference_index,
                l2a_dir=args.l2a_dictionary,
                output_dir=args.output,
                candidate_source=args.candidate_source,
                scl_mode=args.scl_mode,
                erosion_radius_m=args.erosion_radius_m,
                workers=args.workers,
                runtime_source=args.runtime_source,
                day_aod_source=args.day_aod_source,
                maiac_cache=args.maiac_cache,
                cams_dir=args.cams,
                force=args.force,
                dump_grid=args.dump_grid,
                library_year_mode=args.library_year_mode,
                distance_quality_radius_m=args.distance_quality_radius_m,
                aod_quality_mode=args.aod_quality_mode,
                ocm_confidence_quality=args.ocm_confidence_quality,
                ocm_thin_quality_weight=args.ocm_thin_quality_weight,
                maiac_gap_policy=args.maiac_gap_policy,
            )
        except Exception as exc:  # noqa: BLE001 - one array task must report its own failure
            print(
                json.dumps(
                    {
                        "matchup_id": matchup_id,
                        "status": "failed",
                        "error": f"{type(exc).__name__}: {exc}",
                    }
                ),
                flush=True,
            )
            raise


if __name__ == "__main__":
    main()
