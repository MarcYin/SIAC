"""Persistence helpers for prepared monthly composite collections."""

from __future__ import annotations

import json
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Any, cast

import numpy as np
import xarray as xr

from siac.algorithms.surface.brdf_monthly_composite import (
    MonthlyBestPixelComposite,
    MonthlyCompositeCollection,
    MonthlyKernelWeightComposite,
)
from siac.domain import SensorBand
from siac.runtime import BRDFKernelWeights

_MANIFEST_NAME = "manifest.json"
_STORE_VERSION = 2
_SUPPORTED_STORE_VERSIONS = frozenset({1, 2})


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
        return cls(
            bounds=(xmin, ymin, xmax, ymax),
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
        relative_path = f"{label}.zarr"
        dataset = _composite_to_dataset(composite)
        dataset.to_zarr(
            root / relative_path,
            mode="w",
            consolidated=False,
            zarr_format=2,
            write_empty_chunks=True,
        )
        entries.append(
            MonthlyCompositeStoreEntry(
                label=label,
                path=relative_path,
                year=int(composite.year),
                month=int(composite.month),
                kind=_composite_kind(composite),
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
        _dataset_to_composite(
            _open_monthly_dataset(root / entry.path),
            kind=entry.kind,
            year=entry.year,
            month=entry.month,
        )
        for entry in manifest.entries
    )
    return MonthlyCompositeCollection(
        composites=composites,
        source_bands=manifest.source_bands,
        source_name=manifest.source_name,
    )


def _open_monthly_dataset(path: Path) -> xr.Dataset:
    if path.suffix == ".nc":
        dataset = xr.open_dataset(path)
    elif path.suffix == ".zarr":
        dataset = xr.open_zarr(path, consolidated=False)
    else:
        raise ValueError(f"Unsupported monthly composite dataset format: {path.suffix}")
    return dataset.load()


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
    except Exception:
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
    except Exception:
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
