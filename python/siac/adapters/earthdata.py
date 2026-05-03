"""Earthdata adapter helpers."""

from __future__ import annotations

import logging
from contextlib import ExitStack
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Any, cast
from uuid import uuid4

import numpy as np
import xarray as xr

from siac.adapters.data.earthaccess_catalog import EarthAccessCatalog
from siac.adapters.data.earthaccess_source import EarthAccessSource
from siac.adapters.earthdata_common import (
    build_target_template,
    granule_intersects_bounds,
    parse_granule_date,
    parse_tile_indices,
)
from siac.geo.reprojection import transform_bounds

if TYPE_CHECKING:
    from collections.abc import Sequence
    from datetime import datetime

    from rasterio.enums import Resampling

    from siac.adapters.auth import CredentialManager

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class EarthAccessRuntime:
    """Normalized Earthaccess runtime dependencies for provider wrappers."""

    cache_dir: Path | None
    source: EarthAccessSource
    catalog: EarthAccessCatalog


@dataclass(frozen=True)
class _VirtualStackSpec:
    """Temporary GDAL VRT stack plus cleanup metadata."""

    path: str
    cleanup_paths: tuple[str, ...]
    expected_layers: int


def earthaccess_source_from_auth(
    auth: CredentialManager | None,
    *,
    provider: str | None = None,
) -> EarthAccessSource:
    """Build an ``EarthAccessSource`` using configured credentials when available."""
    if auth is None:
        return EarthAccessSource(provider=provider)
    return cast("EarthAccessSource", auth.earthdata().build_earthaccess_source(provider=provider))


def build_earthaccess_runtime(
    *,
    cache_dir: str | Path | None,
    source: EarthAccessSource | None,
    catalog: EarthAccessCatalog | None,
    provider: str | None,
) -> EarthAccessRuntime:
    """Build normalized source/catalog/cache dependencies for Earthaccess providers."""
    resolved_cache_dir = Path(cache_dir).expanduser() if cache_dir is not None else None
    resolved_source = source or EarthAccessSource(provider=provider)
    resolved_catalog = catalog or EarthAccessCatalog(source=resolved_source)
    return EarthAccessRuntime(
        cache_dir=resolved_cache_dir,
        source=resolved_source,
        catalog=resolved_catalog,
    )


def earthaccess_cache_dir(cache_dir: Path | None, short_name: str) -> Path:
    """Return the on-disk cache destination for a downloaded Earthaccess product."""
    return cache_dir or Path.home() / ".cache" / "siac" / "earthdata" / short_name


def select_candidate_paths(
    paths: list[Path],
    *,
    obs_time: datetime,
    bounds: tuple[float, float, float, float],
    crs: str,
    sample_dates: np.ndarray | None = None,
) -> list[Path]:
    """Rank downloaded granules by tile, temporal proximity, and filename."""
    if not paths:
        return []

    sample_day_set: set[np.datetime64] | None = None
    if sample_dates is not None:
        sample_day_set = {np.datetime64(day, "D") for day in sample_dates.tolist()}

    selected: list[tuple[tuple[int, int], float, str, Path]] = []
    for path in paths:
        try:
            granule_dt = parse_granule_date(path)
            granule_day = np.datetime64(granule_dt.date(), "D")
            if sample_day_set is not None and granule_day not in sample_day_set:
                continue
            delta = abs((granule_dt - obs_time).total_seconds())
            intersects = granule_intersects_bounds(path, bounds=bounds, crs=crs)
        except Exception:
            return paths

        if not intersects:
            continue

        try:
            tile = parse_tile_indices(path)
        except Exception:
            tile = (999, 999)
        selected.append((tile, delta, path.name, path))

    if not selected:
        return []
    return [item[3] for item in sorted(selected, key=lambda value: (value[0], value[1], value[2]))]


def target_grid_coords(
    bounds: tuple[float, float, float, float],
    resolution: float,
    *,
    resolution_name: str = "resolution",
) -> tuple[np.ndarray, np.ndarray]:
    """Build center-point y/x coordinates for a target raster grid."""
    xmin, ymin, xmax, ymax = bounds
    if resolution <= 0:
        raise ValueError(f"{resolution_name} must be > 0, got {resolution}")

    nx = max(1, int(np.ceil((xmax - xmin) / resolution)))
    ny = max(1, int(np.ceil((ymax - ymin) / resolution)))
    x = xmin + (np.arange(nx, dtype=np.float32) + 0.5) * resolution
    y = ymax - (np.arange(ny, dtype=np.float32) + 0.5) * resolution
    return y, x


def constant_target_array(
    bounds: tuple[float, float, float, float],
    resolution: float,
    value: float,
    *,
    resolution_name: str = "resolution",
) -> xr.DataArray:
    """Return a constant 2-D target-grid raster."""
    y, x = target_grid_coords(bounds, resolution, resolution_name=resolution_name)
    arr = cast("np.ndarray", np.full((y.size, x.size), value, dtype=np.float32))
    return xr.DataArray(arr, dims=["y", "x"], coords={"y": y, "x": x})


def constant_target_band_array(
    bands: Sequence[str],
    bounds: tuple[float, float, float, float],
    resolution: float,
    value: float,
    *,
    resolution_name: str = "resolution",
) -> xr.DataArray:
    """Return a constant 3-D ``(band, y, x)`` target-grid raster."""
    y, x = target_grid_coords(bounds, resolution, resolution_name=resolution_name)
    arr = cast("np.ndarray", np.full((len(bands), y.size, x.size), value, dtype=np.float32))
    return xr.DataArray(
        arr,
        dims=["band", "y", "x"],
        coords={"band": np.array(bands, dtype=object), "y": y, "x": x},
    )


def merge_reprojected_tiles(
    arrays: list[xr.DataArray],
    *,
    bounds: tuple[float, float, float, float],
    crs: str,
    resolution: float | None,
    resampling: Resampling,
    nodata: float | None,
    target_template: xr.DataArray | None = None,
) -> xr.DataArray:
    """Merge native tiles by reading a virtual mosaic directly onto the target AOI."""
    if not arrays:
        raise ValueError("Expected at least one array to merge")

    reference = arrays[0]
    if target_template is not None:
        target = target_template
        resolved_resolution = _resolve_target_template_resolution(
            target_template,
            resolution,
            resolution_name="resolution",
        )
    else:
        normalized_resolution = _normalize_resolution_value(
            resolution,
            resolution_name="resolution",
        )
        if normalized_resolution is None:
            target, resolved_resolution = _build_source_aligned_target_template_from_array(
                reference,
                bounds=bounds,
                crs=crs,
                resolution=resolution,
                resolution_name="resolution",
            )
        else:
            resolved_resolution = normalized_resolution
            target = build_target_template(bounds, crs, resolved_resolution)
    if reproject_native_to_target is not _DEFAULT_REPROJECT_NATIVE_TO_TARGET:
        reprojected = [
            reproject_native_to_target(
                arr,
                bounds=bounds,
                crs=crs,
                resolution=resolved_resolution,
                resampling=resampling,
                nodata=nodata,
                target_template=target,
            )
            for arr in arrays
        ]
        merged = reprojected[0].astype(np.float32)
        for arr in reprojected[1:]:
            merged = merged.where(np.isfinite(merged), arr.astype(np.float32))
        return merged.astype(np.float32)
    return _merge_tiles_via_vrt(
        arrays,
        bounds=bounds,
        crs=crs,
        resolution=resolved_resolution,
        resampling=resampling,
        nodata=nodata,
        target=target,
    )


def reproject_native_to_target(
    data: xr.DataArray,
    *,
    bounds: tuple[float, float, float, float],
    crs: str,
    resolution: float | None,
    resampling: Resampling,
    nodata: float | None,
    target_template: xr.DataArray | None = None,
) -> xr.DataArray:
    """Compatibility wrapper for callers that still monkeypatch per-array reprojection."""
    if target_template is not None:
        target = target_template
        resolved_resolution = _resolve_target_template_resolution(
            target_template,
            resolution,
            resolution_name="resolution",
        )
    else:
        normalized_resolution = _normalize_resolution_value(
            resolution,
            resolution_name="resolution",
        )
        if normalized_resolution is None:
            target, resolved_resolution = _build_source_aligned_target_template_from_array(
                data,
                bounds=bounds,
                crs=crs,
                resolution=resolution,
                resolution_name="resolution",
            )
        else:
            resolved_resolution = normalized_resolution
            target = build_target_template(bounds, crs, resolved_resolution)
    return _merge_tiles_via_vrt(
        [data],
        bounds=bounds,
        crs=crs,
        resolution=resolved_resolution,
        resampling=resampling,
        nodata=nodata,
        target=target,
    )


_DEFAULT_REPROJECT_NATIVE_TO_TARGET = reproject_native_to_target


def read_virtual_stack_to_target(
    source_groups: Sequence[Sequence[str]],
    *,
    group_band_counts: Sequence[int],
    bounds: tuple[float, float, float, float],
    crs: str,
    resolution: float | None,
    resampling: Resampling,
    nodata: float | None,
    target_template: xr.DataArray | None = None,
) -> xr.DataArray:
    """Read a stacked GDAL VRT directly onto the target AOI."""
    if not source_groups:
        raise ValueError("Expected at least one source group")
    if len(source_groups) != len(group_band_counts):
        raise ValueError("source_groups and group_band_counts must have the same length")

    from osgeo import gdal  # type: ignore[import-untyped]

    gdal.UseExceptions()
    if target_template is not None:
        target = target_template
        resolved_resolution = _resolve_target_template_resolution(
            target_template,
            resolution,
            resolution_name="resolution",
        )
    else:
        normalized_resolution = _normalize_resolution_value(
            resolution,
            resolution_name="resolution",
        )
        if normalized_resolution is None:
            first_source = next((group[0] for group in source_groups if group), None)
            if first_source is None:
                raise ValueError("Expected at least one source path")
            target, resolved_resolution = build_source_aligned_target_template(
                first_source,
                bounds=bounds,
                crs=crs,
                resolution=resolution,
                resolution_name="resolution",
            )
        else:
            resolved_resolution = normalized_resolution
            target = build_target_template(bounds, crs, resolved_resolution)
    target_transform = _target_transform(target, resolution=resolved_resolution)
    target_bounds = _target_bounds_from_template(target, resolution=resolved_resolution)
    width = int(target.sizes["x"])
    height = int(target.sizes["y"])
    warped = None

    logger.info(
        "GDAL virtual stack start: groups=%d layers=%d target=%dx%d crs=%s resolution=%s bounds=%s",
        len(source_groups),
        int(sum(group_band_counts)),
        width,
        height,
        str(target.rio.crs),
        resolved_resolution,
        target_bounds,
    )

    stack_spec = _build_virtual_stack_vrt(
        source_groups,
        group_band_counts=group_band_counts,
        target_bounds=target_bounds,
        target_crs=str(target.rio.crs),
    )
    try:
        warp_kwargs: dict[str, object] = {
            "format": "MEM",
            "outputBounds": target_bounds,
            "width": width,
            "height": height,
            "dstSRS": str(target.rio.crs),
            "resampleAlg": _normalize_gdal_resampling(resampling),
            "outputType": gdal.GDT_Float32,
            "multithread": True,
            "warpOptions": ["NUM_THREADS=ALL_CPUS"],
        }
        if nodata is not None:
            warp_kwargs["dstNodata"] = float(nodata)

        logger.info(
            "GDAL virtual stack: warping master VRT to AOI with resampling=%s dst_nodata=%s",
            warp_kwargs["resampleAlg"],
            warp_kwargs.get("dstNodata"),
        )
        warped = gdal.Warp("", stack_spec.path, **warp_kwargs)
        if warped is None:
            raise RuntimeError("GDAL failed to warp stacked VRT")

        values = np.asarray(warped.ReadAsArray(), dtype=np.float32)
        logger.info("GDAL virtual stack complete: warped shape=%s", tuple(values.shape))
    finally:
        warped = None
        for path in stack_spec.cleanup_paths:
            gdal.Unlink(path)

    if values.ndim == 2:
        values = values[np.newaxis, :, :]

    if values.shape[0] != stack_spec.expected_layers:
        raise ValueError(
            f"Expected {stack_spec.expected_layers} stacked layer(s), got {values.shape[0]}"
        )

    stacked = xr.DataArray(
        values,
        dims=("layer", "y", "x"),
        coords={
            "layer": np.arange(values.shape[0], dtype=np.int32),
            "y": target.coords["y"],
            "x": target.coords["x"],
        },
    )
    stacked = stacked.rio.set_spatial_dims(x_dim="x", y_dim="y")
    stacked = stacked.rio.write_crs(target.rio.crs)
    return stacked.rio.write_transform(target_transform)


def _build_virtual_stack_vrt(
    source_groups: Sequence[Sequence[str]],
    *,
    group_band_counts: Sequence[int],
    target_bounds: tuple[float, float, float, float],
    target_crs: str,
) -> _VirtualStackSpec:
    """Build a cropped, stacked VRT in ``/vsimem`` without warping it."""
    if not source_groups:
        raise ValueError("Expected at least one source group")
    if len(source_groups) != len(group_band_counts):
        raise ValueError("source_groups and group_band_counts must have the same length")

    from osgeo import gdal

    gdal.UseExceptions()
    expected_layers = int(sum(group_band_counts))
    vrt_prefix = f"/vsimem/siac_{uuid4().hex}"
    cleanup_paths = [f"{vrt_prefix}_stack.vrt"]
    group_vrt_paths: list[str] = []
    master_vrt_path = cleanup_paths[0]

    for index, sources in enumerate(source_groups):
        if not sources:
            raise ValueError(f"source group {index} is empty")
        group_vrt_path = f"{vrt_prefix}_group_{index}.vrt"
        cleanup_paths.append(group_vrt_path)
        logger.debug(
            "GDAL virtual stack group %d/%d: sources=%d bands=%d first=%s last=%s",
            index + 1,
            len(source_groups),
            len(sources),
            group_band_counts[index],
            sources[0],
            sources[-1],
        )
        group_vrt = gdal.BuildVRT(group_vrt_path, list(sources))
        if group_vrt is None:
            raise RuntimeError(f"GDAL failed to build VRT for source group {index}")

        cropped_group_vrt_path = _crop_virtual_group_to_target_bounds(
            gdal,
            source_dataset=group_vrt,
            source_path=group_vrt_path,
            target_bounds=target_bounds,
            target_crs=target_crs,
            output_path=f"{vrt_prefix}_group_{index}_crop.vrt",
        )
        if cropped_group_vrt_path != group_vrt_path:
            cleanup_paths.append(cropped_group_vrt_path)
        group_vrt = None
        group_vrt_paths.append(cropped_group_vrt_path)

    logger.info(
        "GDAL virtual stack: building master VRT from %d group VRT(s)", len(group_vrt_paths)
    )
    master_vrt = gdal.BuildVRT(
        master_vrt_path,
        group_vrt_paths,
        options=gdal.BuildVRTOptions(separate=True),
    )
    if master_vrt is None:
        raise RuntimeError("GDAL failed to build stacked VRT")
    master_vrt = None
    return _VirtualStackSpec(
        path=master_vrt_path,
        cleanup_paths=tuple(dict.fromkeys(cleanup_paths)),
        expected_layers=expected_layers,
    )


def _crop_virtual_group_to_target_bounds(
    gdal_module: Any,
    *,
    source_dataset: Any,
    source_path: str,
    target_bounds: tuple[float, float, float, float],
    target_crs: str,
    output_path: str,
) -> str:
    source_crs = _dataset_projection(source_dataset)
    if not source_crs:
        return source_path

    dataset_bounds = _dataset_bounds(source_dataset)
    if dataset_bounds is None:
        return source_path

    try:
        projected_bounds = transform_bounds(target_bounds, target_crs, source_crs)
    except Exception:
        logger.debug(
            "GDAL virtual stack: native bounds transform failed for %s", source_path, exc_info=True
        )
        return source_path

    clipped_bounds = _intersect_bounds(projected_bounds, dataset_bounds)
    if clipped_bounds is None or _bounds_close(clipped_bounds, dataset_bounds):
        return source_path

    logger.debug(
        "GDAL virtual stack: cropping %s to native bounds=%s in crs=%s",
        source_path,
        clipped_bounds,
        source_crs,
    )
    cropped = gdal_module.Translate(
        output_path,
        source_dataset,
        format="VRT",
        projWin=[clipped_bounds[0], clipped_bounds[3], clipped_bounds[2], clipped_bounds[1]],
        projWinSRS=source_crs,
    )
    if cropped is None:
        raise RuntimeError(f"GDAL failed to crop VRT for {source_path}")
    cropped = None
    return output_path


def _dataset_projection(dataset: Any) -> str | None:
    get_projection = getattr(dataset, "GetProjection", None)
    if callable(get_projection):
        projection = get_projection()
        if projection:
            return str(projection)

    get_projection_ref = getattr(dataset, "GetProjectionRef", None)
    if callable(get_projection_ref):
        projection_ref = get_projection_ref()
        if projection_ref:
            return str(projection_ref)
    return None


def build_source_aligned_target_template(
    source_path: str | Path,
    *,
    bounds: tuple[float, float, float, float],
    crs: str,
    resolution: float | None,
    resolution_name: str = "resolution",
) -> tuple[xr.DataArray, float]:
    from osgeo import gdal

    gdal.UseExceptions()
    dataset = gdal.Open(str(source_path))
    if dataset is None:
        raise RuntimeError(f"GDAL failed to open source dataset: {source_path}")
    try:
        source_crs = _dataset_projection(dataset) or str(crs)
        resolved_resolution = _resolve_dataset_resolution(
            dataset,
            resolution,
            resolution_name=resolution_name,
        )
        source_bounds = _transform_bounds_to_crs(bounds, crs, source_crs)
        return build_target_template(
            source_bounds, source_crs, resolved_resolution
        ), resolved_resolution
    finally:
        dataset = None


def _dataset_bounds(dataset: Any) -> tuple[float, float, float, float] | None:
    get_geotransform = getattr(dataset, "GetGeoTransform", None)
    if not callable(get_geotransform):
        return None

    try:
        geotransform = get_geotransform(can_return_null=True)
    except TypeError:
        geotransform = get_geotransform()
    except Exception:
        return None

    if geotransform is None or len(geotransform) != 6:
        return None
    if abs(float(geotransform[2])) > 1e-9 or abs(float(geotransform[4])) > 1e-9:
        return None

    width = int(getattr(dataset, "RasterXSize", 0))
    height = int(getattr(dataset, "RasterYSize", 0))
    if width <= 0 or height <= 0:
        return None

    x0 = float(geotransform[0])
    y0 = float(geotransform[3])
    x1 = x0 + float(geotransform[1]) * width
    y1 = y0 + float(geotransform[5]) * height
    return (
        min(x0, x1),
        min(y0, y1),
        max(x0, x1),
        max(y0, y1),
    )


def _dataset_resolution(dataset: Any) -> float | None:
    get_geotransform = getattr(dataset, "GetGeoTransform", None)
    if not callable(get_geotransform):
        return None

    try:
        geotransform = get_geotransform(can_return_null=True)
    except TypeError:
        geotransform = get_geotransform()
    except Exception:
        return None

    if geotransform is None or len(geotransform) != 6:
        return None
    if abs(float(geotransform[2])) > 1e-9 or abs(float(geotransform[4])) > 1e-9:
        return None
    return _coerce_square_resolution(float(geotransform[1]), float(geotransform[5]))


def _intersect_bounds(
    left: tuple[float, float, float, float],
    right: tuple[float, float, float, float],
) -> tuple[float, float, float, float] | None:
    xmin = max(float(left[0]), float(right[0]))
    ymin = max(float(left[1]), float(right[1]))
    xmax = min(float(left[2]), float(right[2]))
    ymax = min(float(left[3]), float(right[3]))
    if xmin >= xmax or ymin >= ymax:
        return None
    return (xmin, ymin, xmax, ymax)


def _bounds_close(
    left: tuple[float, float, float, float],
    right: tuple[float, float, float, float],
    *,
    atol: float = 1e-6,
) -> bool:
    return bool(
        np.allclose(
            np.asarray(left, dtype=np.float64),
            np.asarray(right, dtype=np.float64),
            atol=atol,
            rtol=0.0,
        )
    )


def _resolve_target_template_resolution(
    target: xr.DataArray,
    resolution: float | None,
    *,
    resolution_name: str,
) -> float:
    resolved = _normalize_resolution_value(resolution, resolution_name=resolution_name)
    if resolved is not None:
        return resolved
    inferred = _data_array_resolution(target)
    if inferred is None:
        raise ValueError(f"Could not infer {resolution_name} from target template")
    return inferred


def _resolve_dataset_resolution(
    dataset: Any,
    resolution: float | None,
    *,
    resolution_name: str,
) -> float:
    resolved = _normalize_resolution_value(resolution, resolution_name=resolution_name)
    if resolved is not None:
        return resolved
    inferred = _dataset_resolution(dataset)
    if inferred is None:
        raise ValueError(f"Could not infer {resolution_name} from source dataset")
    return inferred


def _build_source_aligned_target_template_from_array(
    data: xr.DataArray,
    *,
    bounds: tuple[float, float, float, float],
    crs: str,
    resolution: float | None,
    resolution_name: str,
) -> tuple[xr.DataArray, float]:
    source_crs = str(data.rio.crs) if data.rio.crs is not None else str(crs)
    resolved_resolution = _normalize_resolution_value(resolution, resolution_name=resolution_name)
    if resolved_resolution is None:
        resolved_resolution = _data_array_resolution(data)
    if resolved_resolution is None:
        raise ValueError(f"Could not infer {resolution_name} from source array")
    source_bounds = _transform_bounds_to_crs(bounds, crs, source_crs)
    return build_target_template(
        source_bounds, source_crs, resolved_resolution
    ), resolved_resolution


def _normalize_resolution_value(
    resolution: float | None,
    *,
    resolution_name: str,
) -> float | None:
    if resolution is None:
        return None
    value = float(resolution)
    if value <= 0:
        raise ValueError(f"{resolution_name} must be > 0, got {resolution}")
    return value


def _data_array_resolution(data: xr.DataArray) -> float | None:
    x_values = np.asarray(data.coords.get("x", xr.DataArray([])).values, dtype=np.float64)
    y_values = np.asarray(data.coords.get("y", xr.DataArray([])).values, dtype=np.float64)
    x_resolution = _axis_resolution(x_values)
    y_resolution = _axis_resolution(y_values)
    if x_resolution is not None and y_resolution is not None:
        return _coerce_square_resolution(x_resolution, y_resolution)
    if x_resolution is not None:
        return float(abs(x_resolution))
    if y_resolution is not None:
        return float(abs(y_resolution))
    return None


def _axis_resolution(values: np.ndarray) -> float | None:
    if values.size <= 1:
        return None
    diffs = np.diff(values)
    magnitudes = np.abs(diffs[np.nonzero(diffs)])
    if magnitudes.size == 0:
        return None
    resolution = float(np.median(magnitudes))
    if not np.allclose(magnitudes, resolution, rtol=0.0, atol=1e-9):
        return None
    return resolution


def _coerce_square_resolution(x_resolution: float, y_resolution: float) -> float | None:
    x_magnitude = abs(float(x_resolution))
    y_magnitude = abs(float(y_resolution))
    if x_magnitude <= 0 or y_magnitude <= 0:
        return None
    if not np.isclose(x_magnitude, y_magnitude, rtol=0.0, atol=1e-9):
        raise ValueError(
            f"Source dataset uses non-square pixels ({x_magnitude}, {y_magnitude}); scalar resolution is ambiguous"
        )
    return x_magnitude


def _transform_bounds_to_crs(
    bounds: tuple[float, float, float, float],
    source_crs: str,
    target_crs: str,
) -> tuple[float, float, float, float]:
    if str(source_crs) == str(target_crs):
        return (
            float(bounds[0]),
            float(bounds[1]),
            float(bounds[2]),
            float(bounds[3]),
        )
    return transform_bounds(bounds, source_crs, target_crs)


def _merge_tiles_via_vrt(
    arrays: list[xr.DataArray],
    *,
    bounds: tuple[float, float, float, float],
    crs: str,
    resolution: float,
    resampling: Resampling,
    nodata: float | None,
    target: xr.DataArray,
) -> xr.DataArray:
    from rasterio.enums import Resampling as RasterioResampling
    from rasterio.io import MemoryFile
    from rasterio.merge import merge
    from rasterio.transform import array_bounds
    from rasterio.vrt import WarpedVRT

    reference = arrays[0]
    if reference.dims[-2:] != ("y", "x"):
        raise ValueError("Expected spatial dimensions to be the trailing ('y', 'x') axes")
    for arr in arrays[1:]:
        if not _merge_layout_compatible(reference, arr):
            raise ValueError(
                "All arrays must share the same band layout and non-spatial coordinates"
            )

    resampling_enum = _normalize_resampling(resampling, RasterioResampling)
    target_transform = _target_transform(target, resolution=resolution)
    source_tiles = [arr for arr in arrays if _array_intersects_target(arr, bounds=bounds, crs=crs)]
    if not source_tiles:
        empty = _empty_target_like(target, reference=reference, nodata=nodata)
        return empty.rio.write_transform(target_transform)
    target_height = int(target.sizes["y"])
    target_width = int(target.sizes["x"])
    target_bounds = array_bounds(target_height, target_width, target_transform)
    band_count = _band_count(reference)
    dst_nodata = np.float32(np.nan if nodata is None else nodata)

    with ExitStack() as stack:
        sources = []
        for tile in source_tiles:
            memfile = stack.enter_context(MemoryFile())
            _write_dataarray_to_memory_dataset(memfile, tile, nodata=dst_nodata)
            dataset = stack.enter_context(memfile.open())
            vrt = stack.enter_context(
                WarpedVRT(
                    dataset,
                    crs=str(target.rio.crs),
                    transform=target_transform,
                    width=target_width,
                    height=target_height,
                    resampling=resampling_enum,
                    src_nodata=dst_nodata,
                    nodata=dst_nodata,
                )
            )
            sources.append(vrt)

        merged_values, _ = merge(
            sources,
            bounds=target_bounds,
            res=(abs(float(target_transform.a)), abs(float(target_transform.e))),
            nodata=dst_nodata,
            dtype="float32",
            output_count=band_count,
            resampling=resampling_enum,
            method="first",
        )

    merged = _restore_merged_dataarray(
        merged_values,
        reference=reference,
        target=target,
    )
    merged = merged.rio.write_transform(target_transform)
    if nodata is not None:
        merged = merged.rio.write_nodata(float(nodata))
    return merged


def _array_intersects_target(
    data: xr.DataArray,
    *,
    bounds: tuple[float, float, float, float],
    crs: str,
) -> bool:
    from rasterio.transform import array_bounds

    if data.dims[-2:] != ("y", "x"):
        raise ValueError("Expected spatial dimensions to be the trailing ('y', 'x') axes")

    source_bounds = transform_bounds(bounds, crs, str(data.rio.crs))
    data_bounds = array_bounds(
        int(data.sizes["y"]),
        int(data.sizes["x"]),
        data.rio.transform(recalc=True),
    )
    return not (
        float(data_bounds[2]) < float(source_bounds[0])
        or float(data_bounds[0]) > float(source_bounds[2])
        or float(data_bounds[3]) < float(source_bounds[1])
        or float(data_bounds[1]) > float(source_bounds[3])
    )


def _merge_layout_compatible(reference: xr.DataArray, candidate: xr.DataArray) -> bool:
    if tuple(candidate.dims) != tuple(reference.dims):
        return False
    for dim in reference.dims[:-2]:
        if not _coord_layout_compatible(reference.coords[dim], candidate.coords[dim]):
            return False
    grid_mappings = {
        name for name in (reference.rio.grid_mapping, candidate.rio.grid_mapping) if name
    }
    for name, coord in reference.coords.items():
        if (
            name in grid_mappings
            or "x" in coord.dims
            or "y" in coord.dims
            or name in reference.dims
        ):
            continue
        candidate_coord = candidate.coords.get(name)
        if candidate_coord is None or not _coord_layout_compatible(coord, candidate_coord):
            return False
    return True


def _coord_layout_compatible(reference: xr.DataArray, candidate: xr.DataArray) -> bool:
    if tuple(reference.dims) != tuple(candidate.dims):
        return False
    if tuple(reference.shape) != tuple(candidate.shape):
        return False
    left = np.asarray(reference.values)
    right = np.asarray(candidate.values)
    if left.shape != right.shape:
        return False
    try:
        return bool(np.array_equal(left, right, equal_nan=True))
    except TypeError:
        return bool(np.array_equal(left, right))


def _normalize_resampling(resampling: Resampling, enum_type: type) -> object:
    if isinstance(resampling, str):
        try:
            return getattr(enum_type, resampling)
        except AttributeError as exc:
            raise ValueError(f"Unknown resampling mode: {resampling!r}") from exc
    return resampling


def _normalize_gdal_resampling(resampling: Resampling) -> str:
    name = (
        resampling if isinstance(resampling, str) else getattr(resampling, "name", str(resampling))
    )
    normalized = str(name).lower()
    mapping = {
        "nearest": "near",
        "near": "near",
        "bilinear": "bilinear",
        "cubic": "cubic",
        "cubicspline": "cubicspline",
        "lanczos": "lanczos",
        "average": "average",
        "mode": "mode",
    }
    if normalized not in mapping:
        raise ValueError(f"Unknown GDAL resampling mode: {resampling!r}")
    return mapping[normalized]


def _band_count(data: xr.DataArray) -> int:
    return int(np.prod(data.shape[:-2], dtype=np.int64)) if data.ndim > 2 else 1


def _affine_from_coords(
    x_values: np.ndarray,
    y_values: np.ndarray,
    *,
    x_step: float,
    y_step: float,
) -> Any:
    from rasterio.transform import Affine

    return Affine.translation(
        x_values[0] - x_step / 2.0, y_values[0] - y_step / 2.0
    ) * Affine.scale(x_step, y_step)


def _target_bounds_from_template(
    target: xr.DataArray,
    *,
    resolution: float,
) -> tuple[float, float, float, float]:
    transform = _target_transform(target, resolution=resolution)
    x_values = np.asarray(target.coords["x"].values, dtype=np.float64)
    y_values = np.asarray(target.coords["y"].values, dtype=np.float64)
    width = int(x_values.size)
    height = int(y_values.size)
    xmin = float(transform.c)
    ymax = float(transform.f)
    xmax = xmin + abs(float(transform.a)) * width
    ymin = ymax + float(transform.e) * height
    return (xmin, ymin, xmax, ymax)


def _target_transform(data: xr.DataArray, *, resolution: float) -> Any:
    x_values = np.asarray(data.coords["x"].values, dtype=np.float64)
    y_values = np.asarray(data.coords["y"].values, dtype=np.float64)
    x_step = float(x_values[1] - x_values[0]) if x_values.size > 1 else float(resolution)
    y_step = float(y_values[1] - y_values[0]) if y_values.size > 1 else -float(resolution)
    return _affine_from_coords(x_values, y_values, x_step=x_step, y_step=y_step)


def _write_dataarray_to_memory_dataset(
    memfile: Any,
    data: xr.DataArray,
    *,
    nodata: np.float32,
) -> None:
    values = np.asarray(data.values, dtype=np.float32)
    band_count = _band_count(data)
    reshaped = values.reshape((band_count, values.shape[-2], values.shape[-1]))
    writer = memfile.open(
        driver="GTiff",
        height=reshaped.shape[1],
        width=reshaped.shape[2],
        count=band_count,
        dtype="float32",
        crs=str(data.rio.crs),
        transform=data.rio.transform(),
        nodata=float(nodata),
    )
    try:
        writer.write(reshaped)
    finally:
        writer.close()


def _restore_merged_dataarray(
    values: np.ndarray,
    *,
    reference: xr.DataArray,
    target: xr.DataArray,
) -> xr.DataArray:
    restored = np.asarray(values, dtype=np.float32).reshape(
        reference.shape[:-2] + values.shape[-2:]
    )
    coords: dict[str, object] = {dim: reference.coords[dim] for dim in reference.dims[:-2]}
    coords["y"] = target.coords["y"]
    coords["x"] = target.coords["x"]
    for name, coord in reference.coords.items():
        if name in coords or "x" in coord.dims or "y" in coord.dims:
            continue
        coords[name] = coord
    merged = xr.DataArray(restored, dims=reference.dims, coords=coords)
    merged = merged.rio.set_spatial_dims(x_dim="x", y_dim="y")
    return merged.rio.write_crs(target.rio.crs)


def _empty_target_like(
    target: xr.DataArray,
    *,
    reference: xr.DataArray,
    nodata: float | None,
) -> xr.DataArray:
    fill_value = np.nan if nodata is None else float(nodata)
    values = np.full(reference.shape[:-2] + target.shape, fill_value, dtype=np.float32)
    coords: dict[str, object] = {dim: reference.coords[dim] for dim in reference.dims[:-2]}
    coords["y"] = target.coords["y"]
    coords["x"] = target.coords["x"]
    for name, coord in reference.coords.items():
        if name in coords or "x" in coord.dims or "y" in coord.dims:
            continue
        coords[name] = coord
    empty = xr.DataArray(values, dims=reference.dims, coords=coords)
    empty = empty.rio.set_spatial_dims(x_dim="x", y_dim="y")
    return empty.rio.write_crs(target.rio.crs)
