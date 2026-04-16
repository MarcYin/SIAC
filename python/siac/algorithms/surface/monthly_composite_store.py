"""Persistence helpers for prepared monthly composite collections."""

from __future__ import annotations

import json
import logging
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Any, cast

import numpy as np
import rasterio
import rioxarray  # noqa: F401
import xarray as xr
from affine import Affine
from rasterio.transform import from_bounds

from siac.algorithms.surface.brdf_monthly_composite import (
    MonthlyBestPixelComposite,
    MonthlyCompositeCollection,
    MonthlyKernelWeightComposite,
)
from siac.domain import SensorBand
from siac.geo.reprojection import reproject_match
from siac.runtime import BRDFKernelWeights
from siac.runtime.models import copy_spatial_metadata_like
from siac.storage.writers import write_raster

logger = logging.getLogger(__name__)

_MANIFEST_NAME = "manifest.json"
_STORE_VERSION = 3
_SUPPORTED_STORE_VERSIONS = frozenset({1, 2, 3})
_GEOTIFF_FORMAT = "geotiff"
_GEOTIFF_EXTENSION = ".tif"


@dataclass(frozen=True)
class MonthlyCompositeStoreGridSpec:
    bounds: tuple[float, float, float, float]
    crs: str
    resolution: float
    width: int
    height: int

    @classmethod
    def from_bounds(
        cls,
        bounds: tuple[float, float, float, float],
        *,
        crs: str,
        resolution: float,
    ) -> MonthlyCompositeStoreGridSpec:
        xmin, ymin, xmax, ymax = (float(value) for value in bounds)
        res = float(resolution)
        width = max(1, int(np.ceil((xmax - xmin) / res)))
        height = max(1, int(np.ceil((ymax - ymin) / res)))
        aligned_xmax = xmin + width * res
        aligned_ymin = ymax - height * res
        return cls(
            bounds=(xmin, aligned_ymin, aligned_xmax, ymax),
            crs=str(crs),
            resolution=res,
            width=width,
            height=height,
        )


@dataclass(frozen=True)
class MonthlyCompositeStoreEntry:
    label: str
    path: str
    year: int
    month: int
    kind: str
    format: str = _GEOTIFF_FORMAT
    assets: dict[str, str] | None = None


@dataclass(frozen=True)
class MonthlyCompositeStoreManifest:
    version: int
    source_name: str | None
    source_bands: tuple[SensorBand, ...]
    entries: tuple[MonthlyCompositeStoreEntry, ...]
    grid: MonthlyCompositeStoreGridSpec | None = None


def filter_materialized_monthly_composite_collection(
    collection: MonthlyCompositeCollection,
) -> MonthlyCompositeCollection:
    """Drop composites whose scientific payload is entirely NaN."""
    return MonthlyCompositeCollection(
        composites=tuple(
            composite
            for composite in collection.composites
            if _composite_has_any_finite_payload(composite)
        ),
        source_bands=collection.source_bands,
        source_name=collection.source_name,
    )


def write_monthly_composite_collection(
    collection: MonthlyCompositeCollection,
    output_path: str | Path,
    *,
    grid: MonthlyCompositeStoreGridSpec | None = None,
) -> Path:
    """Write a prepared monthly composite collection to a directory store."""
    materialized = filter_materialized_monthly_composite_collection(collection)
    root = Path(output_path)
    root.mkdir(parents=True, exist_ok=True)
    stale_paths = _load_existing_store_entry_paths(root)

    entries: list[MonthlyCompositeStoreEntry] = []
    for composite in materialized.composites:
        label = _composite_label(composite.year, composite.month)
        relative_path = label
        assets = _write_monthly_geotiff_period(root / relative_path, composite, grid=grid)
        entries.append(
            MonthlyCompositeStoreEntry(
                label=label,
                path=relative_path,
                year=int(composite.year),
                month=int(composite.month),
                kind=_composite_kind(composite),
                format=_GEOTIFF_FORMAT,
                assets=assets,
            )
        )

    manifest = MonthlyCompositeStoreManifest(
        version=_STORE_VERSION,
        source_name=materialized.source_name,
        source_bands=tuple(materialized.source_bands),
        entries=tuple(entries),
        grid=grid or _infer_collection_grid(materialized),
    )
    _write_store_manifest(root, manifest)
    _remove_stale_store_paths(root, stale_paths - {entry.path for entry in entries})
    return root


def read_monthly_composite_store_manifest(
    input_path: str | Path,
) -> MonthlyCompositeStoreManifest:
    """Read a prepared monthly composite manifest without loading the datasets."""
    root = Path(input_path)
    manifest_path = root / _MANIFEST_NAME
    if not manifest_path.exists():
        raise FileNotFoundError(f"Monthly composite manifest not found: {manifest_path}")

    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    version = int(manifest.get("version", 0))
    if version not in _SUPPORTED_STORE_VERSIONS:
        raise ValueError(
            "Unsupported monthly composite store version "
            f"{manifest.get('version')!r}; expected one of {sorted(_SUPPORTED_STORE_VERSIONS)!r}"
        )
    return MonthlyCompositeStoreManifest(
        version=version,
        source_name=manifest.get("source_name"),
        source_bands=tuple(_deserialize_sensor_band(item) for item in manifest.get("source_bands", ())),
        entries=tuple(
            MonthlyCompositeStoreEntry(
                label=str(entry["label"]),
                path=str(entry["path"]),
                year=int(entry["year"]),
                month=int(entry["month"]),
                kind=str(entry["kind"]),
                format=str(entry.get("format") or _infer_entry_format(str(entry["path"]))),
                assets=(
                    {str(name): str(path) for name, path in cast("dict[str, Any]", entry["assets"]).items()}
                    if isinstance(entry.get("assets"), dict)
                    else None
                ),
            )
            for entry in manifest.get("entries", ())
        ),
        grid=(
            _deserialize_grid_spec(manifest["grid"])
            if manifest.get("grid") is not None
            else None
        ),
    )


def read_monthly_composite_collection(
    input_path: str | Path,
) -> MonthlyCompositeCollection:
    """Read a prepared monthly composite collection from a directory store."""
    root = Path(input_path)
    manifest = read_monthly_composite_store_manifest(root)
    composites = tuple(
        _read_store_entry(
            root,
            entry=entry,
            source_bands=manifest.source_bands,
        )
        for entry in manifest.entries
    )
    return MonthlyCompositeCollection(
        composites=composites,
        source_bands=manifest.source_bands,
        source_name=manifest.source_name,
    )


def _read_store_entry(
    root: Path,
    *,
    entry: MonthlyCompositeStoreEntry,
    source_bands: tuple[SensorBand, ...],
) -> MonthlyBestPixelComposite | MonthlyKernelWeightComposite:
    if entry.format == _GEOTIFF_FORMAT:
        return _read_geotiff_composite(
            root / entry.path,
            kind=entry.kind,
            year=entry.year,
            month=entry.month,
            source_bands=source_bands,
            assets=entry.assets,
        )
    dataset = _open_monthly_dataset(root / entry.path)
    return _dataset_to_composite(
        dataset,
        kind=entry.kind,
        year=entry.year,
        month=entry.month,
    )


def _open_monthly_dataset(path: Path) -> xr.Dataset:
    if path.suffix == ".nc":
        dataset = xr.open_dataset(path)
    elif path.suffix == ".zarr":
        dataset = xr.open_zarr(path, consolidated=False)
    else:
        raise ValueError(f"Unsupported monthly composite dataset format: {path.suffix}")
    return dataset.load()


def _infer_entry_format(path_text: str) -> str:
    suffix = Path(path_text).suffix.lower()
    if suffix in {".tif", ".tiff"}:
        return _GEOTIFF_FORMAT
    if suffix == ".zarr":
        return "zarr"
    if suffix == ".nc":
        return "netcdf"
    return _GEOTIFF_FORMAT


def _write_monthly_geotiff_period(
    period_root: Path,
    composite: MonthlyBestPixelComposite | MonthlyKernelWeightComposite,
    *,
    grid: MonthlyCompositeStoreGridSpec | None = None,
) -> dict[str, str]:
    if period_root.exists():
        if period_root.is_dir():
            shutil.rmtree(period_root)
        else:
            period_root.unlink()
    period_root.mkdir(parents=True, exist_ok=True)

    assets: dict[str, str] = {}
    if isinstance(composite, MonthlyKernelWeightComposite):
        asset_sources: list[tuple[str, xr.DataArray, str]] = [
            ("f0", composite.kernels.f0, "float32"),
            ("f1", composite.kernels.f1, "float32"),
            ("f2", composite.kernels.f2, "float32"),
            ("f0_unc", composite.kernels.f0_unc, "float32"),
            ("f1_unc", composite.kernels.f1_unc, "float32"),
            ("f2_unc", composite.kernels.f2_unc, "float32"),
        ]
        if composite.kernels.reflectance_unc is not None:
            asset_sources.append(("reflectance_unc", composite.kernels.reflectance_unc, "float32"))
    else:
        asset_sources = [("reflectance", composite.reflectance, "float32")]

    asset_sources.extend(
        [
            ("quality", composite.quality, "float32"),
            ("sample_index", composite.sample_index, "int16"),
        ]
    )
    for asset_name, data, dtype in asset_sources:
        relative_asset = f"{asset_name}{_GEOTIFF_EXTENSION}"
        _write_geotiff_asset(
            period_root / relative_asset,
            data,
            asset_name=asset_name,
            dtype=dtype,
            grid=grid,
        )
        assets[asset_name] = relative_asset
    return assets


def _write_geotiff_asset(
    path: Path,
    data: xr.DataArray,
    *,
    asset_name: str,
    dtype: str,
    grid: MonthlyCompositeStoreGridSpec | None = None,
) -> None:
    nodata = _asset_nodata_value(asset_name=asset_name, dtype=dtype)
    prepared = _prepare_raster_asset(data, grid=grid, nodata=nodata)
    write_raster(
        prepared,
        path,
        compression="deflate",
        dtype=dtype,
        nodata=nodata,
        blockxsize=_tiff_block_size(int(prepared.sizes["x"])),
        blockysize=_tiff_block_size(int(prepared.sizes["y"])),
    )
    descriptions = _asset_band_descriptions(prepared, asset_name=asset_name)
    with rasterio.open(path, "r+") as dst:
        dst.descriptions = descriptions
        dst.update_tags(siac_asset=asset_name)


def _prepare_raster_asset(
    data: xr.DataArray,
    *,
    grid: MonthlyCompositeStoreGridSpec | None = None,
    nodata: float | int | None = None,
) -> xr.DataArray:
    source_crs = _data_array_crs(data)
    extra_coords = [name for name in data.coords if name not in data.dims and name != "spatial_ref"]
    if extra_coords:
        data = data.drop_vars(extra_coords, errors="ignore")
    if "x" not in data.coords or "y" not in data.coords:
        raise ValueError("Monthly composite assets require x/y coordinates for GeoTIFF persistence")
    prepared = data.transpose("band", "y", "x") if "band" in data.dims else data.transpose("y", "x")
    transform = _affine_from_coords(prepared.coords["x"].values, prepared.coords["y"].values)
    prepared = prepared.rio.set_spatial_dims(x_dim="x", y_dim="y")
    if source_crs is not None:
        prepared = prepared.rio.write_crs(source_crs)
    prepared = prepared.rio.write_transform(transform)
    if nodata is not None:
        prepared = prepared.rio.write_nodata(nodata)
    if grid is None:
        return prepared
    return _align_raster_asset_to_grid(prepared, grid, nodata=nodata)


def _align_raster_asset_to_grid(
    data: xr.DataArray,
    grid: MonthlyCompositeStoreGridSpec,
    *,
    nodata: float | int | None,
) -> xr.DataArray:
    template = _grid_template(grid)
    if _data_array_crs(data) is None:
        if _matches_template_grid(data, template):
            aligned = copy_spatial_metadata_like(data, template)
            return aligned.astype(data.dtype, copy=False)
        raise ValueError(
            "Monthly composite assets require source CRS metadata for explicit-grid reprojection"
        )
    aligned = reproject_match(
        data,
        template,
        resampling="nearest",
        nodata=nodata,
    )
    return aligned.astype(data.dtype, copy=False)


def _grid_template(grid: MonthlyCompositeStoreGridSpec) -> xr.DataArray:
    xmin, ymin, xmax, ymax = (float(value) for value in grid.bounds)
    width = int(grid.width)
    height = int(grid.height)
    transform = from_bounds(xmin, ymin, xmax, ymax, width, height)
    y, x = _coords_from_transform(transform, height=height, width=width)
    template = xr.DataArray(
        np.full((height, width), np.nan, dtype=np.float32),
        dims=("y", "x"),
        coords={"x": x, "y": y},
    )
    template = template.rio.set_spatial_dims(x_dim="x", y_dim="y")
    template = template.rio.write_crs(grid.crs)
    return template.rio.write_transform(transform)


def _matches_template_grid(data: xr.DataArray, template: xr.DataArray) -> bool:
    if data.sizes.get("x") != template.sizes.get("x") or data.sizes.get("y") != template.sizes.get("y"):
        return False
    return np.array_equal(np.asarray(data.coords["x"].values), np.asarray(template.coords["x"].values)) and np.array_equal(
        np.asarray(data.coords["y"].values),
        np.asarray(template.coords["y"].values),
    )


def _asset_nodata_value(
    *,
    asset_name: str,
    dtype: str,
) -> float | int | None:
    if np.dtype(dtype).kind == "f":
        return np.nan
    if asset_name == "sample_index":
        return -1
    if np.dtype(dtype).kind in {"i", "u"}:
        return 0 if np.dtype(dtype).kind == "u" else -9999
    return None


def _asset_band_descriptions(
    data: xr.DataArray,
    *,
    asset_name: str,
) -> tuple[str, ...]:
    if "band" not in data.dims:
        return (asset_name,)
    return tuple(str(value) for value in data.coords["band"].values)


def _read_geotiff_composite(
    period_root: Path,
    *,
    kind: str,
    year: int,
    month: int,
    source_bands: tuple[SensorBand, ...],
    assets: dict[str, str] | None,
) -> MonthlyBestPixelComposite | MonthlyKernelWeightComposite:
    resolved_assets = _resolve_geotiff_assets(period_root, kind=kind, assets=assets)
    default_band_names = tuple(band.name for band in source_bands)

    if kind == "kernel_weights":
        first_band_data = _read_geotiff_asset(
            period_root / resolved_assets["f0"],
            band_names=default_band_names or None,
        ).astype(np.float32)
        band_names = tuple(str(value) for value in first_band_data.coords["band"].values)
        kernel_arrays = {
            "f0": first_band_data,
            "f1": _read_geotiff_asset(period_root / resolved_assets["f1"], band_names=band_names).astype(np.float32),
            "f2": _read_geotiff_asset(period_root / resolved_assets["f2"], band_names=band_names).astype(np.float32),
            "f0_unc": _read_geotiff_asset(period_root / resolved_assets["f0_unc"], band_names=band_names).astype(np.float32),
            "f1_unc": _read_geotiff_asset(period_root / resolved_assets["f1_unc"], band_names=band_names).astype(np.float32),
            "f2_unc": _read_geotiff_asset(period_root / resolved_assets["f2_unc"], band_names=band_names).astype(np.float32),
        }
        reflectance_unc = (
            _read_geotiff_asset(
                period_root / resolved_assets["reflectance_unc"],
                band_names=band_names,
            ).astype(np.float32)
            if "reflectance_unc" in resolved_assets
            else None
        )
        quality = _read_geotiff_asset(period_root / resolved_assets["quality"]).astype(np.float32)
        sample_index = _read_sample_index_asset(period_root / resolved_assets["sample_index"])
        return MonthlyKernelWeightComposite(
            kernels=BRDFKernelWeights(
                f0=kernel_arrays["f0"],
                f1=kernel_arrays["f1"],
                f2=kernel_arrays["f2"],
                f0_unc=kernel_arrays["f0_unc"],
                f1_unc=kernel_arrays["f1_unc"],
                f2_unc=kernel_arrays["f2_unc"],
                reflectance_unc=reflectance_unc,
            ),
            quality=quality,
            sample_index=sample_index,
            year=year,
            month=month,
        )
    if kind != "reflectance":
        raise ValueError(f"Unsupported monthly composite kind: {kind!r}")

    reflectance = _read_geotiff_asset(
        period_root / resolved_assets["reflectance"],
        band_names=default_band_names or None,
    ).astype(np.float32)
    quality = _read_geotiff_asset(period_root / resolved_assets["quality"]).astype(np.float32)
    sample_index = _read_sample_index_asset(period_root / resolved_assets["sample_index"])
    return MonthlyBestPixelComposite(
        reflectance=reflectance,
        quality=quality,
        sample_index=sample_index,
        year=year,
        month=month,
    )


def _resolve_geotiff_assets(
    period_root: Path,
    *,
    kind: str,
    assets: dict[str, str] | None,
) -> dict[str, str]:
    if assets is not None:
        return {str(name): str(path) for name, path in assets.items()}
    resolved = {
        "quality": f"quality{_GEOTIFF_EXTENSION}",
        "sample_index": f"sample_index{_GEOTIFF_EXTENSION}",
    }
    if kind == "kernel_weights":
        for name in ("f0", "f1", "f2", "f0_unc", "f1_unc", "f2_unc"):
            resolved[name] = f"{name}{_GEOTIFF_EXTENSION}"
        reflectance_unc_path = period_root / f"reflectance_unc{_GEOTIFF_EXTENSION}"
        if reflectance_unc_path.exists():
            resolved["reflectance_unc"] = reflectance_unc_path.name
        return resolved
    if kind == "reflectance":
        resolved["reflectance"] = f"reflectance{_GEOTIFF_EXTENSION}"
        return resolved
    raise ValueError(f"Unsupported monthly composite kind: {kind!r}")


def _read_geotiff_asset(
    path: Path,
    *,
    band_names: tuple[str, ...] | None = None,
) -> xr.DataArray:
    with rasterio.open(path) as src:
        values = src.read()
        y_coords, x_coords = _coords_from_transform(src.transform, height=src.height, width=src.width)
        if src.count == 1:
            data = xr.DataArray(
                values[0],
                dims=("y", "x"),
                coords={"y": y_coords, "x": x_coords},
            )
        else:
            resolved_band_names = _resolve_asset_band_names(src, band_names=band_names)
            data = xr.DataArray(
                values,
                dims=("band", "y", "x"),
                coords={
                    "band": np.asarray(resolved_band_names, dtype=object),
                    "y": y_coords,
                    "x": x_coords,
                },
            )
        data = data.rio.set_spatial_dims(x_dim="x", y_dim="y")
        if src.crs is not None:
            data = data.rio.write_crs(src.crs.to_string())
        data = data.rio.write_transform(src.transform)
    return data


def _read_sample_index_asset(path: Path) -> xr.DataArray:
    data = _read_geotiff_asset(path)
    rounded = np.rint(np.asarray(data.values, dtype=np.float32)).astype(np.int16)
    return copy_spatial_metadata_like(
        xr.DataArray(rounded, dims=data.dims, coords=data.coords),
        data,
    )


def _resolve_asset_band_names(
    src: rasterio.io.DatasetReader,
    *,
    band_names: tuple[str, ...] | None,
) -> tuple[str, ...]:
    if band_names is not None and len(band_names) == src.count:
        return band_names
    descriptions = tuple(str(description) for description in src.descriptions if description)
    if len(descriptions) == src.count:
        return descriptions
    raise ValueError(f"Could not determine multiband names for monthly composite asset {src.name}")


def _coords_from_transform(
    transform: Affine,
    *,
    height: int,
    width: int,
) -> tuple[np.ndarray, np.ndarray]:
    if not np.isclose(transform.b, 0.0, rtol=0.0, atol=1e-12) or not np.isclose(transform.d, 0.0, rtol=0.0, atol=1e-12):
        raise ValueError("Monthly composite GeoTIFF assets must use axis-aligned transforms")
    x = transform.c + (np.arange(width, dtype=np.float32) + 0.5) * np.float32(transform.a)
    y = transform.f + (np.arange(height, dtype=np.float32) + 0.5) * np.float32(transform.e)
    return y.astype(np.float32), x.astype(np.float32)


def _affine_from_coords(
    x_values: Any,
    y_values: Any,
) -> Affine:
    x = np.asarray(x_values, dtype=np.float64)
    y = np.asarray(y_values, dtype=np.float64)
    if x.size == 0 or y.size == 0:
        raise ValueError("Monthly composite GeoTIFF assets require non-empty x/y coordinates")

    x_resolution = _axis_resolution(x)
    y_resolution = _axis_resolution(y)
    if x.size > 1 and x_resolution is None:
        raise ValueError("Monthly composite GeoTIFF assets require regularly spaced x coordinates")
    if y.size > 1 and y_resolution is None:
        raise ValueError("Monthly composite GeoTIFF assets require regularly spaced y coordinates")
    x_sign = _axis_sign(x) or 1.0
    y_sign = _axis_sign(y) or -1.0

    x_step = (x_sign * x_resolution) if x_resolution is not None else (abs(y_resolution) if y_resolution is not None else 1.0)
    y_step = (y_sign * y_resolution) if y_resolution is not None else (-abs(x_step) if x_step != 0 else -1.0)
    return Affine(
        float(x_step),
        0.0,
        float(x[0] - (x_step / 2.0)),
        0.0,
        float(y_step),
        float(y[0] - (y_step / 2.0)),
    )


def _axis_sign(values: np.ndarray) -> float | None:
    if values.size <= 1:
        return None
    diffs = np.diff(values)
    nonzero = diffs[np.abs(diffs) > 1e-9]
    if nonzero.size == 0:
        return None
    return 1.0 if float(np.median(nonzero)) >= 0.0 else -1.0


def _tiff_block_size(length: int) -> int:
    if length <= 16:
        return 16
    return int(min(512, ((length + 15) // 16) * 16))


def _data_array_crs(data: xr.DataArray) -> str | None:
    try:
        return str(data.rio.crs) if data.rio.crs is not None else None
    except (AttributeError, ValueError, RuntimeError) as exc:
        logger.debug("rio.crs accessor failed on DataArray: %s", exc)
        return None


def _composite_label(year: int, month: int) -> str:
    return f"{int(year):04d}_{int(month):02d}"


def _composite_kind(
    composite: MonthlyBestPixelComposite | MonthlyKernelWeightComposite,
) -> str:
    if isinstance(composite, MonthlyKernelWeightComposite):
        return "kernel_weights"
    return "reflectance"


def _composite_to_dataset(
    composite: MonthlyBestPixelComposite | MonthlyKernelWeightComposite,
) -> xr.Dataset:
    if isinstance(composite, MonthlyKernelWeightComposite):
        dataset = xr.Dataset(
            {
                "f0": composite.kernels.f0.astype(np.float32),
                "f1": composite.kernels.f1.astype(np.float32),
                "f2": composite.kernels.f2.astype(np.float32),
                "f0_unc": composite.kernels.f0_unc.astype(np.float32),
                "f1_unc": composite.kernels.f1_unc.astype(np.float32),
                "f2_unc": composite.kernels.f2_unc.astype(np.float32),
                "quality": composite.quality.astype(np.float32),
                "sample_index": composite.sample_index.astype(np.int16),
            }
        )
        if composite.kernels.reflectance_unc is not None:
            dataset["reflectance_unc"] = composite.kernels.reflectance_unc.astype(np.float32)
    else:
        dataset = xr.Dataset(
            {
                "reflectance": composite.reflectance.astype(np.float32),
                "quality": composite.quality.astype(np.float32),
                "sample_index": composite.sample_index.astype(np.int16),
            }
        )

    dataset.attrs.update(
        {
            "year": int(composite.year),
            "month": int(composite.month),
            "composite_kind": _composite_kind(composite),
        }
    )
    return dataset


def _dataset_to_composite(
    dataset: xr.Dataset,
    *,
    kind: str,
    year: int,
    month: int,
) -> MonthlyBestPixelComposite | MonthlyKernelWeightComposite:
    quality = cast("xr.DataArray", dataset["quality"]).astype(np.float32)
    sample_index = cast("xr.DataArray", dataset["sample_index"]).astype(np.int16)
    if kind == "kernel_weights":
        reflectance_unc = (
            cast("xr.DataArray", dataset["reflectance_unc"]).astype(np.float32)
            if "reflectance_unc" in dataset.data_vars
            else None
        )
        return MonthlyKernelWeightComposite(
            kernels=BRDFKernelWeights(
                f0=cast("xr.DataArray", dataset["f0"]).astype(np.float32),
                f1=cast("xr.DataArray", dataset["f1"]).astype(np.float32),
                f2=cast("xr.DataArray", dataset["f2"]).astype(np.float32),
                f0_unc=cast("xr.DataArray", dataset["f0_unc"]).astype(np.float32),
                f1_unc=cast("xr.DataArray", dataset["f1_unc"]).astype(np.float32),
                f2_unc=cast("xr.DataArray", dataset["f2_unc"]).astype(np.float32),
                reflectance_unc=reflectance_unc,
            ),
            quality=quality,
            sample_index=sample_index,
            year=year,
            month=month,
        )
    if kind != "reflectance":
        raise ValueError(f"Unsupported monthly composite kind: {kind!r}")
    return MonthlyBestPixelComposite(
        reflectance=cast("xr.DataArray", dataset["reflectance"]).astype(np.float32),
        quality=quality,
        sample_index=sample_index,
        year=year,
        month=month,
    )


def _composite_has_any_finite_payload(
    composite: MonthlyBestPixelComposite | MonthlyKernelWeightComposite,
) -> bool:
    sample_index = np.asarray(composite.sample_index.values)
    if sample_index.size > 0 and np.all(sample_index < 0):
        return False
    if isinstance(composite, MonthlyKernelWeightComposite):
        return any(
            _data_array_has_any_finite_value(data)
            for data in (composite.kernels.f0, composite.kernels.f1, composite.kernels.f2)
        )
    return _data_array_has_any_finite_value(composite.reflectance)


def _data_array_has_any_finite_value(data: xr.DataArray) -> bool:
    values = np.asarray(data.values)
    if not np.issubdtype(values.dtype, np.floating):
        return False
    return bool(np.isfinite(values).any())


def _serialize_sensor_band(band: SensorBand) -> dict[str, Any]:
    return {
        "name": band.name,
        "center_wavelength": float(band.center_wavelength),
        "bandwidth": float(band.bandwidth),
        "resolution": float(band.resolution),
        "band_index": int(band.band_index),
        "rsrf_wavelengths_nm": (
            np.asarray(band.rsrf_wavelengths_nm, dtype=np.float64).tolist()
            if band.rsrf_wavelengths_nm is not None
            else None
        ),
        "rsrf_response": (
            np.asarray(band.rsrf_response, dtype=np.float64).tolist()
            if band.rsrf_response is not None
            else None
        ),
        "rsrf_sensor_unit_id": band.rsrf_sensor_unit_id,
        "rsrf_representation_variant": band.rsrf_representation_variant,
        "rsrf_band_id": band.rsrf_band_id,
    }


def _deserialize_sensor_band(payload: dict[str, Any]) -> SensorBand:
    rsrf_wavelengths_nm = payload.get("rsrf_wavelengths_nm")
    rsrf_response = payload.get("rsrf_response")
    return SensorBand(
        str(payload["name"]),
        float(payload["center_wavelength"]),
        float(payload["bandwidth"]),
        float(payload["resolution"]),
        int(payload["band_index"]),
        rsrf_wavelengths_nm=(
            np.asarray(rsrf_wavelengths_nm, dtype=np.float64)
            if rsrf_wavelengths_nm is not None
            else None
        ),
        rsrf_response=(
            np.asarray(rsrf_response, dtype=np.float64)
            if rsrf_response is not None
            else None
        ),
        rsrf_sensor_unit_id=cast("str | None", payload.get("rsrf_sensor_unit_id")),
        rsrf_representation_variant=cast("str | None", payload.get("rsrf_representation_variant")),
        rsrf_band_id=cast("str | None", payload.get("rsrf_band_id")),
    )


def _write_store_manifest(root: Path, manifest: MonthlyCompositeStoreManifest) -> None:
    payload = {
        "version": manifest.version,
        "source_name": manifest.source_name,
        "source_bands": [_serialize_sensor_band(band) for band in manifest.source_bands],
        "entries": [
            {
                "label": entry.label,
                "path": entry.path,
                "year": entry.year,
                "month": entry.month,
                "kind": entry.kind,
                "format": entry.format,
                "assets": entry.assets,
            }
            for entry in manifest.entries
        ],
        "grid": _serialize_grid_spec(manifest.grid),
    }
    (root / _MANIFEST_NAME).write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")


def _load_existing_store_entry_paths(root: Path) -> set[str]:
    manifest_path = root / _MANIFEST_NAME
    if not manifest_path.exists():
        return set()
    try:
        manifest = read_monthly_composite_store_manifest(root)
    except (OSError, json.JSONDecodeError, ValueError, KeyError, TypeError) as exc:
        logger.warning(
            "Could not read existing monthly composite manifest at %s (%s: %s); "
            "treating store as empty, stale entries will not be cleaned up.",
            manifest_path, type(exc).__name__, exc,
        )
        return set()
    return {entry.path for entry in manifest.entries}


def _remove_stale_store_paths(root: Path, relative_paths: set[str]) -> None:
    for relative_path in relative_paths:
        target = root / relative_path
        if target.is_dir():
            shutil.rmtree(target)
        elif target.exists():
            target.unlink()


def _serialize_grid_spec(grid: MonthlyCompositeStoreGridSpec | None) -> dict[str, Any] | None:
    if grid is None:
        return None
    return {
        "bounds": [float(value) for value in grid.bounds],
        "crs": grid.crs,
        "resolution": float(grid.resolution),
        "width": int(grid.width),
        "height": int(grid.height),
    }


def _deserialize_grid_spec(payload: dict[str, Any]) -> MonthlyCompositeStoreGridSpec:
    bounds = cast("tuple[float, float, float, float]", tuple(float(value) for value in payload["bounds"]))
    return MonthlyCompositeStoreGridSpec(
        bounds=bounds,
        crs=str(payload["crs"]),
        resolution=float(payload["resolution"]),
        width=int(payload["width"]),
        height=int(payload["height"]),
    )


def _infer_collection_grid(
    collection: MonthlyCompositeCollection,
) -> MonthlyCompositeStoreGridSpec | None:
    if not collection.composites:
        return None
    return _infer_grid_from_data(_representative_spatial_data(collection.composites[0]))


def _representative_spatial_data(
    composite: MonthlyBestPixelComposite | MonthlyKernelWeightComposite,
) -> xr.DataArray:
    if isinstance(composite, MonthlyKernelWeightComposite):
        return composite.kernels.f0
    return composite.reflectance


def _infer_grid_from_data(data: xr.DataArray) -> MonthlyCompositeStoreGridSpec | None:
    if "x" not in data.coords or "y" not in data.coords:
        return None
    if data.sizes.get("x", 0) == 0 or data.sizes.get("y", 0) == 0:
        return None

    x = np.asarray(data.coords["x"].values, dtype=np.float64)
    y = np.asarray(data.coords["y"].values, dtype=np.float64)

    try:
        crs = str(data.rio.crs) if data.rio.crs is not None else None
    except (AttributeError, ValueError, RuntimeError) as exc:
        logger.debug("rio.crs accessor failed while inferring grid: %s", exc)
        crs = None
    if crs is None:
        return None

    x_res = _axis_resolution(x)
    y_res = _axis_resolution(y)
    if x_res is None or y_res is None or not np.isclose(x_res, y_res, rtol=0.0, atol=1e-9):
        return None

    bounds = (
        float(x[0] - x_res / 2.0),
        float(y[-1] - y_res / 2.0),
        float(x[-1] + x_res / 2.0),
        float(y[0] + y_res / 2.0),
    )
    return MonthlyCompositeStoreGridSpec(
        bounds=bounds,
        crs=crs,
        resolution=float(x_res),
        width=int(data.sizes["x"]),
        height=int(data.sizes["y"]),
    )


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


__all__ = [
    "MonthlyCompositeStoreEntry",
    "MonthlyCompositeStoreGridSpec",
    "MonthlyCompositeStoreManifest",
    "filter_materialized_monthly_composite_collection",
    "read_monthly_composite_collection",
    "read_monthly_composite_store_manifest",
    "write_monthly_composite_collection",
]
