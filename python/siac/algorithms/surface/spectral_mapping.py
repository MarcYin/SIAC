"""Spectral mapping between differing source and target sensor band sets."""

from __future__ import annotations

import csv
import hashlib
import importlib
import json
import logging
import os
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from time import perf_counter
from typing import TYPE_CHECKING, Protocol, TypeAlias, cast

import numpy as np
import xarray as xr
from numpy import typing as npt

from siac.algorithms.surface._spectral_curve_utils import (
    canonicalize_curve as _canonicalize_curve,
)
from siac.algorithms.surface._spectral_curve_utils import (
    gaussian_curve_from_band as _gaussian_curve_from_band,
)
from siac.algorithms.surface._spectral_curve_utils import (
    normalized_band_response as _normalized_band_response,
)
from siac.algorithms.surface._spectral_curve_utils import (
    primary_nir_band_index as _primary_nir_band_index,
)
from siac.algorithms.surface._spectral_curve_utils import (
    segment_for_band as _segment_for_band,
)
from siac.algorithms.surface._spectral_curve_utils import (
    segmentize_curve as _segmentize_curve,
)
from siac.runtime.models import copy_spatial_metadata_like
from siac.storage.writers import write_raster

if TYPE_CHECKING:
    from collections.abc import Callable, Sequence

    from siac.domain import SensorBand


logger = logging.getLogger(__name__)


Float32Array: TypeAlias = npt.NDArray[np.float32]
Float64Array: TypeAlias = npt.NDArray[np.float64]
BoolArray: TypeAlias = npt.NDArray[np.bool_]
Int64Array: TypeAlias = npt.NDArray[np.int64]


class _RSRFCurveProtocol(Protocol):
    wavelength_nm: object
    response: object


class _RSRFModuleProtocol(Protocol):
    def load_response_definition(
        self,
        sensor_unit_id: str,
        band_id: str,
        representation_variant: str,
        *,
        root: Path | None = None,
    ) -> object: ...

    def realize_curve(self, response_definition: object) -> _RSRFCurveProtocol: ...


class _BatchArrayResultProtocol(Protocol):
    """Protocol matching ``spectral_library.BatchMappingArrayResult`` (v0.3.0+)."""

    reflectance: np.ndarray
    source_fit_rmse: np.ndarray
    output_columns: tuple[str, ...]


class _LazyModuleProxy:
    """Import-on-first-use proxy used to avoid eager heavy imports."""

    def __init__(self, module_name: str) -> None:
        self._module_name = module_name
        self._module_cache: object | None = None

    def _module(self) -> object:
        if self._module_cache is None:
            self._module_cache = importlib.import_module(self._module_name)
        return self._module_cache

    def __getattr__(self, name: str) -> object:
        return getattr(self._module(), name)


@lru_cache(maxsize=1)
def _spectral_library_module() -> object:
    return importlib.import_module("spectral_library")


def PackageSpectralMapper(*args: object, **kwargs: object) -> object:
    factory = cast("Callable[..., object]", _spectral_library_module().SpectralMapper)  # type: ignore[union-attr]
    return factory(*args, **kwargs)


def prepare_package_mapping_library(*args: object, **kwargs: object) -> object:
    factory = cast("Callable[..., object]", _spectral_library_module().prepare_mapping_library)  # type: ignore[union-attr]
    return factory(*args, **kwargs)


rsrf = cast("_RSRFModuleProtocol", _LazyModuleProxy("rsrf"))

_DEFAULT_K_NEIGHBORS = 5
_UNCERTAINTY_FLOOR = 0.005
_CANONICAL_WAVELENGTHS_NM: Float32Array = np.arange(400.0, 2501.0, 1.0, dtype=np.float32)
_RSRF_ROOT_ENV = "RSRF_ROOT"
_SIAC_SPECTRAL_LIBRARY_ROOT_ENV = "SIAC_SPECTRAL_LIBRARY_ROOT"
_SIAC_SPECTRAL_MAPPING_CACHE_DIR_ENV = "SIAC_SPECTRAL_MAPPING_CACHE_DIR"
_DEFAULT_NEIGHBOR_ESTIMATOR = "distance_weighted_mean"
_DEFAULT_KNN_BACKEND = "scipy_ckdtree"

@dataclass(frozen=True)
class HyperspectralLibrary:
    """Hyperspectral library sampled on a common wavelength axis."""

    wavelengths_nm: np.ndarray
    spectra: np.ndarray
    sample_ids: tuple[str, ...]
    source_id: str = "siac-default-spectral-library"
    source_version: str = "1"

    def __post_init__(self) -> None:
        wavelengths = np.asarray(self.wavelengths_nm, dtype=np.float32)
        spectra = np.asarray(self.spectra, dtype=np.float32)
        if wavelengths.ndim != 1 or wavelengths.size < 2:
            raise ValueError("wavelengths_nm must be a 1-D array with at least two samples")
        if spectra.ndim != 2 or spectra.shape[1] != wavelengths.size:
            raise ValueError("spectra must have shape (n_samples, n_wavelengths)")
        if spectra.shape[0] < 1:
            raise ValueError("spectra must contain at least one sample")
        if len(self.sample_ids) != spectra.shape[0]:
            raise ValueError("sample_ids must match the number of spectra")
        if not np.all(np.isfinite(wavelengths)) or np.any(np.diff(wavelengths) <= 0.0):
            raise ValueError("wavelengths_nm must be finite and strictly increasing")
        if not np.all(np.isfinite(spectra)):
            raise ValueError("spectra must be finite")
        if np.any((spectra < 0.0) | (spectra > 1.5)):
            raise ValueError("spectra must be bounded to a physically plausible reflectance range")
        object.__setattr__(self, "wavelengths_nm", wavelengths)
        object.__setattr__(self, "spectra", spectra)


@dataclass(frozen=True)
class SpectralMappingConfig:
    """Configuration for the package-backed spectral-mapping adapter."""

    siac_library_root: Path | str | None = None
    rsrf_root: Path | str | None = None
    cache_dir: Path | str | None = None
    neighbor_estimator: str = _DEFAULT_NEIGHBOR_ESTIMATOR
    knn_backend: str = _DEFAULT_KNN_BACKEND
    knn_eps: float = 0.0
    min_valid_bands: int = 1

    def normalized(self) -> SpectralMappingConfig:
        def _path(value: Path | str | None) -> Path | None:
            if value is None:
                return None
            return Path(value).expanduser().resolve()

        return SpectralMappingConfig(
            siac_library_root=_path(self.siac_library_root),
            rsrf_root=_path(self.rsrf_root),
            cache_dir=_path(self.cache_dir),
            neighbor_estimator=str(self.neighbor_estimator).strip() or _DEFAULT_NEIGHBOR_ESTIMATOR,
            knn_backend=str(self.knn_backend).strip() or _DEFAULT_KNN_BACKEND,
            knn_eps=float(self.knn_eps),
            min_valid_bands=max(1, int(self.min_valid_bands)),
        )


@dataclass(frozen=True)
class _SchemaBand:
    original_name: str
    schema_band_id: str
    segment: str
    original_index: int


@dataclass(frozen=True)
class _PreparedRuntime:
    prepared_root: Path
    source_sensor_id: str
    target_sensor_id: str
    source_bands: tuple[_SchemaBand, ...]
    target_bands: tuple[_SchemaBand, ...]


@dataclass(frozen=True)
class _FlattenedSourceCube:
    original_dims: tuple[str, ...]
    transpose_dims: tuple[str, ...]
    spatial_dims: tuple[str, ...]
    source_data: xr.DataArray
    flat_values: Float32Array
    valid_rows: BoolArray


@dataclass(frozen=True)
class _DeduplicatedQueries:
    queries: Float64Array
    valid_masks: BoolArray
    inverse_indices: Int64Array


def needs_spectral_mapping(
    source_bands: Sequence[SensorBand],
    target_bands: Sequence[SensorBand],
    *,
    response_tolerance: float = 5.0e-3,
) -> bool:
    """Return True when source and target band sets differ materially."""
    if len(source_bands) != len(target_bands):
        return True

    for source_band, target_band in zip(source_bands, target_bands, strict=True):
        if source_band.name != target_band.name:
            return True
        common_min = max(
            source_band.center_wavelength - 4.0 * source_band.bandwidth,
            target_band.center_wavelength - 4.0 * target_band.bandwidth,
            350.0,
        )
        common_max = min(
            source_band.center_wavelength + 4.0 * source_band.bandwidth,
            target_band.center_wavelength + 4.0 * target_band.bandwidth,
            2500.0,
        )
        if common_max <= common_min:
            return True
        grid = np.arange(common_min, common_max + 1.0, 1.0, dtype=np.float32)
        source_resp = _normalized_band_response(source_band, grid)
        target_resp = _normalized_band_response(target_band, grid)
        if float(np.max(np.abs(source_resp - target_resp))) > response_tolerance:
            return True
    return False


def convolve_hyperspectral_reflectance(
    reflectance: xr.DataArray,
    wavelengths_nm: np.ndarray,
    target_bands: Sequence[SensorBand],
) -> xr.DataArray:
    """Project a hyperspectral reflectance cube onto multispectral target bands."""
    if "wavelength" not in reflectance.dims:
        raise ValueError("reflectance must have a 'wavelength' dimension")
    source = reflectance.transpose(*[dim for dim in reflectance.dims if dim != "wavelength"], "wavelength")
    if source.sizes["wavelength"] != len(wavelengths_nm):
        raise ValueError("wavelength dimension does not match wavelengths_nm")

    matrix = np.stack(
        [_normalized_band_response(band, np.asarray(wavelengths_nm, dtype=np.float32)) for band in target_bands],
        axis=0,
    ).astype(np.float32)
    values = np.asarray(source.values, dtype=np.float32)
    flat = values.reshape(-1, values.shape[-1])
    projected = flat @ matrix.T
    out_shape = values.shape[:-1] + (len(target_bands),)
    projected = projected.reshape(out_shape)

    other_dims = [dim for dim in source.dims if dim != "wavelength"]
    coords = {dim: source.coords[dim] for dim in other_dims if dim in source.coords}
    coords["band"] = [band.name for band in target_bands]
    return xr.DataArray(
        np.moveaxis(projected, -1, 0),
        dims=["band", *other_dims],
        coords=coords,
    )


def _mapping_config_from_env() -> SpectralMappingConfig:
    return SpectralMappingConfig(
        siac_library_root=os.getenv(_SIAC_SPECTRAL_LIBRARY_ROOT_ENV),
        rsrf_root=os.getenv(_RSRF_ROOT_ENV),
        cache_dir=os.getenv(_SIAC_SPECTRAL_MAPPING_CACHE_DIR_ENV),
    ).normalized()


def _split_mapping_inputs(
    spectral_library: HyperspectralLibrary | SpectralMappingConfig | None,
) -> tuple[HyperspectralLibrary | None, SpectralMappingConfig]:
    if spectral_library is None:
        return None, _mapping_config_from_env()
    if isinstance(spectral_library, HyperspectralLibrary):
        return spectral_library, _mapping_config_from_env()
    if isinstance(spectral_library, SpectralMappingConfig):
        return None, spectral_library.normalized()
    raise TypeError(
        "spectral_library must be a HyperspectralLibrary, SpectralMappingConfig, or None"
    )


def _cache_root(config: SpectralMappingConfig) -> Path:
    if config.cache_dir is not None:
        return Path(config.cache_dir)
    return Path(os.getenv(_SIAC_SPECTRAL_MAPPING_CACHE_DIR_ENV, Path.home() / ".cache" / "siac" / "spectral_mapping"))


def _hash_bytes(*chunks: bytes) -> str:
    digest = hashlib.sha256()
    for chunk in chunks:
        digest.update(chunk)
    return digest.hexdigest()


def _json_hash(payload: dict[str, object]) -> str:
    return _hash_bytes(json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8"))


def _library_signature(library: HyperspectralLibrary | None, config: SpectralMappingConfig) -> str:
    if config.siac_library_root is not None:
        return _hash_bytes(str(config.siac_library_root).encode("utf-8"))
    if library is None:
        raise ValueError(
            "Spectral mapping requires an explicit SIAC spectral library. "
            "Provide spectral_library=HyperspectralLibrary(...) or set "
            "SIAC_SPECTRAL_LIBRARY_ROOT / SpectralMappingConfig.siac_library_root."
        )
    return _hash_bytes(
        np.asarray(library.wavelengths_nm, dtype=np.float32).tobytes(),
        np.asarray(library.spectra, dtype=np.float32).tobytes(),
        "\n".join(library.sample_ids).encode("utf-8"),
        str(library.source_id).encode("utf-8"),
        str(library.source_version).encode("utf-8"),
    )


def _rsrf_curve_for_band(
    band: SensorBand,
    *,
    rsrf_root: Path | None,
) -> tuple[Float32Array, Float32Array] | None:
    if band.rsrf_sensor_unit_id is None:
        return None

    try:
        response_definition = rsrf.load_response_definition(
            band.rsrf_sensor_unit_id,
            band.rsrf_band_id or band.name,
            band.rsrf_representation_variant or "band_average",
            root=rsrf_root,
        )
        if hasattr(response_definition, "wavelength_nm") and hasattr(response_definition, "response"):
            curve = response_definition
        else:
            curve = rsrf.realize_curve(response_definition)
    except Exception as exc:  # pragma: no cover - exercised against real repos
        logger.warning(
            "Falling back to a Gaussian response approximation for band %s because RSRF lookup failed (%s)",
            band.name,
            exc,
        )
        return None

    return cast(
        "tuple[Float32Array, Float32Array]",
        _canonicalize_curve(
            np.asarray(curve.wavelength_nm, dtype=np.float32),
            np.asarray(curve.response, dtype=np.float32),
        ),
    )


def _curve_for_band(
    band: SensorBand,
    *,
    rsrf_root: Path | None,
) -> tuple[Float32Array, Float32Array]:
    if band.has_rsrf:
        return cast(
            "tuple[Float32Array, Float32Array]",
            _canonicalize_curve(
                np.asarray(band.rsrf_wavelengths_nm, dtype=np.float32),
                np.asarray(band.rsrf_response, dtype=np.float32),
            ),
        )
    rsrf_curve = _rsrf_curve_for_band(band, rsrf_root=rsrf_root)
    if rsrf_curve is not None:
        return rsrf_curve
    return cast("tuple[Float32Array, Float32Array]", _gaussian_curve_from_band(band))


def _sensor_id_for_bands(bands: Sequence[SensorBand], *, prefix: str) -> str:
    sensor_unit_ids = {band.rsrf_sensor_unit_id for band in bands if band.rsrf_sensor_unit_id}
    if len(sensor_unit_ids) == 1:
        sensor_id = next(iter(sensor_unit_ids))
        return str(sensor_id)
    payload = [
        {
            "name": band.name,
            "center_wavelength": float(band.center_wavelength),
            "bandwidth": float(band.bandwidth),
            "resolution": float(band.resolution),
            "rsrf_sensor_unit_id": band.rsrf_sensor_unit_id,
            "rsrf_variant": band.rsrf_representation_variant,
            "rsrf_band_id": band.rsrf_band_id,
        }
        for band in bands
    ]
    return f"siac_{prefix}_{_json_hash({'bands': payload})[:12]}"


def _schema_payload_for_bands(
    sensor_id: str,
    bands: Sequence[SensorBand],
    *,
    rsrf_root: Path | None,
) -> tuple[dict[str, object], tuple[_SchemaBand, ...]]:
    primary_nir_index = _primary_nir_band_index(bands)
    used_band_ids: set[str] = set()
    payload_bands: list[dict[str, object]] = []
    schema_bands: list[_SchemaBand] = []

    for index, band in enumerate(bands):
        schema_band_id = "nir" if primary_nir_index == index else band.name
        if schema_band_id in used_band_ids:
            schema_band_id = f"{schema_band_id}_{index}"
        used_band_ids.add(schema_band_id)

        segment = _segment_for_band(band)
        wavelengths_nm, response = _curve_for_band(band, rsrf_root=rsrf_root)
        wavelengths_nm, response = _segmentize_curve(
            wavelengths_nm,
            response,
            segment=segment,
        )
        positive = np.flatnonzero(response > 0.0)
        support_min_nm = float(wavelengths_nm[positive[0]])
        support_max_nm = float(wavelengths_nm[positive[-1]])
        payload_bands.append(
            {
                "band_id": schema_band_id,
                "segment": segment,
                "wavelength_nm": wavelengths_nm.astype(float).tolist(),
                "rsr": response.astype(float).tolist(),
                "center_nm": float(band.center_wavelength),
                "fwhm_nm": float(band.bandwidth),
                "support_min_nm": support_min_nm,
                "support_max_nm": support_max_nm,
            }
        )
        schema_bands.append(
            _SchemaBand(
                original_name=band.name,
                schema_band_id=schema_band_id,
                segment=segment,
                original_index=index,
            )
        )

    return {"sensor_id": sensor_id, "bands": payload_bands}, tuple(schema_bands)


def _export_hyperspectral_library_root(root: Path, library: HyperspectralLibrary) -> None:
    tabular_root = root / "tabular"
    tabular_root.mkdir(parents=True, exist_ok=True)
    metadata_path = tabular_root / "siac_spectra_metadata.csv"
    spectra_path = tabular_root / "siac_normalized_spectra.csv"

    if np.array_equal(library.wavelengths_nm.astype(np.float32), _CANONICAL_WAVELENGTHS_NM):
        canonical = np.asarray(library.spectra, dtype=np.float32)
    else:
        wavelengths = np.asarray(library.wavelengths_nm, dtype=np.float32)
        if wavelengths[0] > float(_CANONICAL_WAVELENGTHS_NM[0]) or wavelengths[-1] < float(_CANONICAL_WAVELENGTHS_NM[-1]):
            raise ValueError(
                "HyperspectralLibrary must cover the canonical 400-2500 nm range to drive spectral-library mapping"
            )
        canonical = np.vstack(
            [
                np.interp(
                    _CANONICAL_WAVELENGTHS_NM,
                    wavelengths,
                    np.asarray(spectrum, dtype=np.float32),
                ).astype(np.float32)
                for spectrum in np.asarray(library.spectra, dtype=np.float32)
            ]
        )

    with metadata_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["source_id", "spectrum_id", "sample_name"])
        writer.writeheader()
        for sample_id in library.sample_ids:
            writer.writerow(
                {
                    "source_id": library.source_id,
                    "spectrum_id": sample_id,
                    "sample_name": sample_id,
                }
            )

    fieldnames = ["source_id", "spectrum_id", "sample_name", *[f"nm_{int(wl)}" for wl in _CANONICAL_WAVELENGTHS_NM]]
    with spectra_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for sample_id, spectrum in zip(library.sample_ids, canonical, strict=True):
            row: dict[str, object] = {
                "source_id": str(library.source_id),
                "spectrum_id": str(sample_id),
                "sample_name": str(sample_id),
            }
            row.update({f"nm_{int(wl)}": float(value) for wl, value in zip(_CANONICAL_WAVELENGTHS_NM, spectrum, strict=True)})
            writer.writerow(row)


def _ensure_siac_library_root(
    cache_root: Path,
    signature: str,
    library: HyperspectralLibrary | None,
    config: SpectralMappingConfig,
) -> Path:
    if config.siac_library_root is not None:
        return Path(config.siac_library_root)
    if library is None:
        raise ValueError(
            "Spectral mapping requires an explicit SIAC spectral library. "
            "Provide spectral_library=HyperspectralLibrary(...) or set "
            "SIAC_SPECTRAL_LIBRARY_ROOT / SpectralMappingConfig.siac_library_root."
        )
    export_root = cache_root / signature / "siac_library"
    metadata_path = export_root / "tabular" / "siac_spectra_metadata.csv"
    spectra_path = export_root / "tabular" / "siac_normalized_spectra.csv"
    if metadata_path.exists() and spectra_path.exists():
        return export_root
    _export_hyperspectral_library_root(export_root, library)
    return export_root


def _write_sensor_schema(path: Path, payload: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _atomic_write_json(path: Path, payload: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    tmp_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    tmp_path.replace(path)


def _atomic_write_npz(path: Path, arrays: dict[str, Float32Array]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    with tmp_path.open("wb") as handle:
        np.savez_compressed(handle, **arrays)
    tmp_path.replace(path)


def _atomic_write_geotiff(path: Path, data: xr.DataArray) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = path.with_name(f".{path.stem}.{os.getpid()}.tmp{path.suffix}")
    write_raster(data, tmp_path, compression="deflate", dtype="float32", tiled=False)
    tmp_path.replace(path)


def _diagnostic_signature(
    metrics: dict[str, xr.DataArray],
    metadata: dict[str, object],
) -> str:
    payload = [
        _json_hash(metadata),
        *[
            _hash_bytes(
                name.encode("utf-8"),
                str(array.dtype).encode("utf-8"),
                "|".join(array.dims).encode("utf-8"),
                np.asarray(array.shape, dtype=np.int64).tobytes(),
                np.asarray(array.values, dtype=np.float32).tobytes(),
                *[
                    _hash_bytes(
                        dim.encode("utf-8"),
                        str(np.asarray(array.coords[dim].values).dtype).encode("utf-8"),
                        np.asarray(np.asarray(array.coords[dim].values).shape, dtype=np.int64).tobytes(),
                        np.asarray(array.coords[dim].values).tobytes(),
                    ).encode("utf-8")
                    for dim in array.dims
                    if dim in array.coords
                ],
                str(array.rio.crs).encode("utf-8"),
                str(array.rio.transform(recalc=True)).encode("utf-8"),
            )
            for name, array in sorted(metrics.items())
        ],
    ]
    return _hash_bytes(*[item.encode("utf-8") for item in payload])


def _write_distance_metric_diagnostics(
    cache_root: Path | str | None,
    *,
    prefix: str,
    metrics: dict[str, xr.DataArray],
    metadata: dict[str, object],
) -> Path | None:
    if cache_root is None or not metrics:
        return None

    normalized_metrics: dict[str, xr.DataArray] = {}
    for name, values in metrics.items():
        metric = values.astype(np.float32)
        if metric.size == 0:
            continue
        try:
            crs = metric.rio.crs
            transform = metric.rio.transform(recalc=True)
        except Exception:
            crs = None
            transform = None
        if crs is None or transform is None:
            logger.warning("Skipping diagnostic GeoTIFF for %s: spatial metadata is incomplete", name)
            continue
        normalized_metrics[str(name)] = metric
    if not normalized_metrics:
        return None

    diagnostics_root = Path(cache_root).expanduser().resolve() / "diagnostics"
    signature = _diagnostic_signature(normalized_metrics, metadata)
    metadata_path = diagnostics_root / f"{prefix}_{signature[:16]}.json"
    metric_paths = {
        name: diagnostics_root / f"{prefix}_{signature[:16]}_{name}.tif"
        for name in sorted(normalized_metrics)
    }
    for name, data_path in metric_paths.items():
        if data_path.exists():
            continue
        _atomic_write_geotiff(data_path, normalized_metrics[name])
    if not metadata_path.exists():
        payload = dict(metadata)
        payload["metrics"] = {
            name: {
                "dtype": str(metric.dtype),
                "shape": [int(size) for size in metric.shape],
                "path": str(metric_paths[name]),
            }
            for name, metric in sorted(normalized_metrics.items())
        }
        _atomic_write_json(metadata_path, payload)
    return metadata_path


def _prepare_runtime(
    source_bands: Sequence[SensorBand],
    target_bands: Sequence[SensorBand],
    *,
    library: HyperspectralLibrary | None,
    config: SpectralMappingConfig,
) -> _PreparedRuntime:
    normalized_config = config.normalized()
    cache_root = _cache_root(normalized_config).expanduser().resolve()
    cache_root.mkdir(parents=True, exist_ok=True)

    source_sensor_id = _sensor_id_for_bands(source_bands, prefix="source")
    target_sensor_id = _sensor_id_for_bands(target_bands, prefix="target")
    source_payload, source_schema_bands = _schema_payload_for_bands(
        source_sensor_id,
        source_bands,
        rsrf_root=cast("Path | None", normalized_config.rsrf_root),
    )
    target_payload, target_schema_bands = _schema_payload_for_bands(
        target_sensor_id,
        target_bands,
        rsrf_root=cast("Path | None", normalized_config.rsrf_root),
    )

    signature = _hash_bytes(
        _library_signature(library, normalized_config).encode("utf-8"),
        _json_hash(source_payload).encode("utf-8"),
        _json_hash(target_payload).encode("utf-8"),
        json.dumps(
            {
                "neighbor_estimator": normalized_config.neighbor_estimator,
                "knn_backend": normalized_config.knn_backend,
                "knn_eps": normalized_config.knn_eps,
                "min_valid_bands": normalized_config.min_valid_bands,
            },
            sort_keys=True,
        ).encode("utf-8"),
    )

    runtime_root = cache_root / signature
    prepared_root = runtime_root / "prepared"
    manifest_path = prepared_root / "manifest.json"
    if manifest_path.exists():
        return _PreparedRuntime(
            prepared_root=prepared_root,
            source_sensor_id=source_sensor_id,
            target_sensor_id=target_sensor_id,
            source_bands=source_schema_bands,
            target_bands=target_schema_bands,
        )

    siac_library_root = _ensure_siac_library_root(cache_root, signature, library, normalized_config)
    srf_root = runtime_root / "srfs"
    _write_sensor_schema(srf_root / f"{source_sensor_id}.json", source_payload)
    _write_sensor_schema(srf_root / f"{target_sensor_id}.json", target_payload)
    prepare_package_mapping_library(
        siac_root=siac_library_root,
        srf_root=srf_root,
        output_root=prepared_root,
        source_sensors=[source_sensor_id],
    )
    return _PreparedRuntime(
        prepared_root=prepared_root,
        source_sensor_id=source_sensor_id,
        target_sensor_id=target_sensor_id,
        source_bands=source_schema_bands,
        target_bands=target_schema_bands,
    )


class SpectralMapper:
    """Package-backed mapper from one multispectral basis to another.

    Uses ``spectral_library >= 0.3.0`` batch array API for efficient
    KNN-based spectral reconstruction.
    """

    def __init__(
        self,
        source_bands: Sequence[SensorBand],
        target_bands: Sequence[SensorBand],
        *,
        spectral_library: HyperspectralLibrary | SpectralMappingConfig | None = None,
        k_neighbors: int = _DEFAULT_K_NEIGHBORS,
    ) -> None:
        if k_neighbors < 1:
            raise ValueError("k_neighbors must be >= 1")

        self.source_bands = tuple(source_bands)
        self.target_bands = tuple(target_bands)
        self.k_neighbors = int(k_neighbors)
        self._identity = not needs_spectral_mapping(self.source_bands, self.target_bands)
        self._target_band_names = [band.name for band in self.target_bands]

        library, mapping_config = _split_mapping_inputs(spectral_library)
        self._mapping_config = mapping_config.normalized()
        self._spectral_library = library
        self._runtime: _PreparedRuntime | None = None
        self._package_mapper: object | None = None
        self._target_internal_to_output_index: dict[str, int] = {}

        if not self._identity:
            self._runtime = _prepare_runtime(
                self.source_bands,
                self.target_bands,
                library=self._spectral_library,
                config=self._mapping_config,
            )
            self._package_mapper = PackageSpectralMapper(self._runtime.prepared_root, verify_checksums=False)
            self._target_internal_to_output_index = {
                band.schema_band_id: band.original_index for band in self._runtime.target_bands
            }

    def map(
        self,
        source_reflectance: xr.DataArray,
        *,
        source_uncertainty: xr.DataArray | None = None,
    ) -> tuple[xr.DataArray, xr.DataArray, xr.DataArray]:
        """Map source-basis multispectral reflectance to target bands."""
        source_data = self._align_source_data(source_reflectance)
        original_dims = tuple(source_data.dims)
        if source_uncertainty is not None:
            source_unc = self._align_source_data(source_uncertainty)
        else:
            source_unc = None

        if self._identity:
            identity = source_data.assign_coords(band=self._target_band_names)
            if source_unc is not None:
                unc = source_unc.assign_coords(band=self._target_band_names)
            else:
                unc = xr.zeros_like(identity, dtype=np.float32) + _UNCERTAINTY_FLOOR
            zero_fit = xr.zeros_like(cast("xr.DataArray", source_data.isel(band=0, drop=True)), dtype=np.float32)
            return (
                identity.transpose(*original_dims).astype(np.float32),
                unc.transpose(*original_dims).astype(np.float32),
                zero_fit.astype(np.float32),
            )

        if self._package_mapper is None or self._runtime is None:
            raise RuntimeError("spectral-library runtime was not initialized")

        flattened = self._flatten_source_cube(source_data)
        unc_flat = self._flatten_optional_uncertainty(source_unc, flattened.transpose_dims)
        n_pixels = flattened.flat_values.shape[0]
        n_target = len(self.target_bands)
        target_flat = np.full((n_pixels, n_target), np.nan, dtype=np.float32)
        target_unc_flat = np.full_like(target_flat, np.nan)
        source_fit_flat = np.full(n_pixels, np.nan, dtype=np.float32)

        if np.any(flattened.valid_rows):
            valid_count = int(np.count_nonzero(flattened.valid_rows))
            logger.info(
                "Spectral mapping start: valid_queries=%d total_queries=%d source_bands=%d target_bands=%d k=%d",
                valid_count,
                n_pixels,
                len(self.source_bands),
                n_target,
                self.k_neighbors,
            )
            batch_started = perf_counter()
            valid_queries = np.asarray(flattened.flat_values[flattened.valid_rows], dtype=np.float64)
            valid_masks = np.isfinite(valid_queries)
            deduplicated = self._deduplicate_queries(valid_queries, valid_masks)
            unique_query_count = int(deduplicated.queries.shape[0])
            if unique_query_count != valid_count:
                logger.info(
                    "Spectral mapping deduplicated queries: %d -> %d unique rows",
                    valid_count,
                    unique_query_count,
                )

            # spectral_library v0.3.0 batch array API — returns arrays directly.
            batch_result = cast("_BatchArrayResultProtocol", self._package_mapper.map_reflectance_batch_arrays_ndarray(  # type: ignore[union-attr]
                source_sensor=self._runtime.source_sensor_id,
                reflectance_rows=deduplicated.queries,
                valid_mask_rows=deduplicated.valid_masks,
                output_mode="target_sensor",
                target_sensor=self._runtime.target_sensor_id,
                k=self.k_neighbors,
                min_valid_bands=self._mapping_config.min_valid_bands,
                neighbor_estimator=self._mapping_config.neighbor_estimator,
                knn_backend=self._mapping_config.knn_backend,
                knn_eps=self._mapping_config.knn_eps,
            ))
            logger.info(
                "Spectral mapping package batch complete: queried_rows=%d elapsed=%.2fs",
                unique_query_count,
                perf_counter() - batch_started,
            )

            valid_indices = np.flatnonzero(flattened.valid_rows)

            # Map batch result columns to target band positions.
            band_id_to_col = self._target_internal_to_output_index
            result_cols: list[int] = []
            target_cols: list[int] = []
            for j, col_id in enumerate(batch_result.output_columns):
                target_j = band_id_to_col.get(col_id)
                if target_j is not None:
                    result_cols.append(j)
                    target_cols.append(target_j)
            result_col_arr = np.array(result_cols, dtype=np.intp)
            target_col_arr = np.array(target_cols, dtype=np.intp)

            # Expand deduplicated results to per-pixel arrays.
            dedup_indices = deduplicated.inverse_indices
            reflectance = np.asarray(batch_result.reflectance, dtype=np.float32)
            fit_rmse = np.asarray(batch_result.source_fit_rmse, dtype=np.float32)

            fill_started = perf_counter()
            target_flat[np.ix_(valid_indices, target_col_arr)] = reflectance[dedup_indices][:, result_col_arr]
            source_fit_flat[valid_indices] = fit_rmse[dedup_indices]
            logger.info(
                "Spectral mapping reflectance fill complete: %d pixels elapsed=%.2fs",
                int(valid_indices.size),
                perf_counter() - fill_started,
            )

            # Vectorized uncertainty: sqrt(floor² + fit_rmse² + input_unc²).
            unc_started = perf_counter()
            per_pixel_fit = fit_rmse[dedup_indices].astype(np.float64)
            input_unc = np.zeros(valid_count, dtype=np.float64)
            if unc_flat is not None:
                unc_sq = np.square(unc_flat[valid_indices].astype(np.float64))
                input_unc = np.sqrt(np.nanmean(unc_sq, axis=1))
            unc_per_pixel = np.sqrt(
                _UNCERTAINTY_FLOOR**2 + per_pixel_fit**2 + input_unc**2
            ).astype(np.float32)
            mapped_mask = np.isfinite(target_flat[valid_indices])
            target_unc_flat[valid_indices] = np.where(
                mapped_mask,
                unc_per_pixel[:, np.newaxis],
                np.nan,
            )
            logger.info(
                "Spectral mapping uncertainty complete: %d pixels elapsed=%.2fs",
                int(valid_indices.size),
                perf_counter() - unc_started,
            )

        reflectance_da = self._restore_target_cube(target_flat, flattened)
        uncertainty_da = self._restore_target_cube(target_unc_flat, flattened)
        source_fit_da = self._restore_spatial_field(source_fit_flat, flattened)
        spatial_reference = cast("xr.DataArray", flattened.source_data.isel(band=0, drop=True))
        source_fit_da = copy_spatial_metadata_like(source_fit_da.astype(np.float32), spatial_reference)
        self._cache_distance_metrics({"source_fit_rmse": source_fit_da})
        logger.info("Spectral mapping complete: output_shape=%s", tuple(reflectance_da.shape))
        return (
            copy_spatial_metadata_like(
                reflectance_da.transpose(*original_dims).astype(np.float32),
                flattened.source_data,
            ),
            copy_spatial_metadata_like(
                uncertainty_da.transpose(*original_dims).astype(np.float32),
                flattened.source_data,
            ),
            source_fit_da,
        )

    def _cache_distance_metrics(self, metrics: dict[str, xr.DataArray]) -> None:
        if self._identity or self._runtime is None:
            return
        prepared_root = getattr(self._runtime, "prepared_root", None)
        if prepared_root is None:
            return
        _write_distance_metric_diagnostics(
            Path(prepared_root).parent,
            prefix="spectral_mapping_distances",
            metrics=metrics,
            metadata={
                "source_sensor_id": self._runtime.source_sensor_id,
                "target_sensor_id": self._runtime.target_sensor_id,
                "source_band_names": [band.name for band in self.source_bands],
                "target_band_names": [band.name for band in self.target_bands],
                "k_neighbors": int(self.k_neighbors),
                "neighbor_estimator": self._mapping_config.neighbor_estimator,
                "knn_backend": self._mapping_config.knn_backend,
            },
        )

    def _align_source_data(self, data: xr.DataArray) -> xr.DataArray:
        if "band" not in data.dims:
            raise ValueError("source_reflectance must have a 'band' dimension")
        missing = [band.name for band in self.source_bands if band.name not in set(data.coords["band"].values.tolist())]
        if missing:
            raise KeyError(f"source_reflectance is missing source bands: {missing}")
        return data.sel(band=[band.name for band in self.source_bands])

    def _flatten_source_cube(self, source_data: xr.DataArray) -> _FlattenedSourceCube:
        transpose_dims = ("band", *[dim for dim in source_data.dims if dim != "band"])
        source_values = np.asarray(source_data.transpose(*transpose_dims).values, dtype=np.float32)
        flat_values = np.asarray(source_values.reshape(source_values.shape[0], -1).T, dtype=np.float32)
        valid_rows = np.asarray(np.any(np.isfinite(flat_values), axis=1), dtype=np.bool_)
        return _FlattenedSourceCube(
            original_dims=tuple(source_data.dims),
            transpose_dims=transpose_dims,
            spatial_dims=transpose_dims[1:],
            source_data=source_data,
            flat_values=flat_values,
            valid_rows=valid_rows,
        )

    @staticmethod
    def _flatten_optional_uncertainty(
        source_uncertainty: xr.DataArray | None,
        transpose_dims: tuple[str, ...],
    ) -> Float32Array | None:
        if source_uncertainty is None:
            return None
        unc_values = np.asarray(source_uncertainty.transpose(*transpose_dims).values, dtype=np.float32)
        return np.asarray(unc_values.reshape(unc_values.shape[0], -1).T, dtype=np.float32)

    @staticmethod
    def _deduplicate_queries(
        valid_queries: Float64Array,
        valid_masks: BoolArray,
    ) -> _DeduplicatedQueries:
        if valid_queries.shape[0] <= 1:
            return _DeduplicatedQueries(
                queries=valid_queries,
                valid_masks=valid_masks,
                inverse_indices=np.arange(valid_queries.shape[0], dtype=np.int64),
            )
        sentinel = np.float32(-9999.0)
        dedup_basis = np.where(valid_masks, valid_queries, sentinel).astype(np.float32, copy=False)
        unique_basis, inverse_indices = np.unique(dedup_basis, axis=0, return_inverse=True)
        unique_masks = unique_basis != sentinel
        unique_queries = unique_basis.astype(np.float64, copy=False)
        unique_queries[~unique_masks] = np.nan
        return _DeduplicatedQueries(
            queries=unique_queries,
            valid_masks=np.asarray(unique_masks, dtype=np.bool_),
            inverse_indices=np.asarray(inverse_indices, dtype=np.int64),
        )

    def _restore_target_cube(
        self,
        flat_values: Float32Array,
        flattened: _FlattenedSourceCube,
    ) -> xr.DataArray:
        shape = (len(self.target_bands), *[flattened.source_data.sizes[dim] for dim in flattened.spatial_dims])
        values = np.asarray(flat_values.T.reshape(shape), dtype=np.float32)
        coords: dict[str, object] = {"band": self._target_band_names}
        coords.update(
            {
                dim: flattened.source_data.coords[dim]
                for dim in flattened.spatial_dims
                if dim in flattened.source_data.coords
            }
        )
        return copy_spatial_metadata_like(
            xr.DataArray(values, dims=flattened.transpose_dims, coords=coords),
            flattened.source_data,
        )

    @staticmethod
    def _restore_spatial_field(
        flat_values: Float32Array,
        flattened: _FlattenedSourceCube,
    ) -> xr.DataArray:
        shape = tuple(flattened.source_data.sizes[dim] for dim in flattened.spatial_dims)
        values = np.asarray(flat_values.reshape(shape), dtype=np.float32)
        coords = {
            dim: flattened.source_data.coords[dim]
            for dim in flattened.spatial_dims
            if dim in flattened.source_data.coords
        }
        return xr.DataArray(values, dims=flattened.spatial_dims, coords=coords)


def map_multispectral_reflectance(
    source_reflectance: xr.DataArray,
    *,
    source_bands: Sequence[SensorBand],
    target_bands: Sequence[SensorBand],
    source_uncertainty: xr.DataArray | None = None,
    spectral_library: HyperspectralLibrary | SpectralMappingConfig | None = None,
    k_neighbors: int = _DEFAULT_K_NEIGHBORS,
) -> tuple[xr.DataArray, xr.DataArray]:
    """Convenience wrapper around :class:`SpectralMapper`."""
    mapper = SpectralMapper(
        source_bands,
        target_bands,
        spectral_library=spectral_library,
        k_neighbors=k_neighbors,
    )
    mapped_reflectance, mapped_uncertainty, _source_fit_rmse = mapper.map(
        source_reflectance,
        source_uncertainty=source_uncertainty,
    )
    return mapped_reflectance, mapped_uncertainty


__all__ = [
    "HyperspectralLibrary",
    "SpectralMappingConfig",
    "SpectralMapper",
    "convolve_hyperspectral_reflectance",
    "map_multispectral_reflectance",
    "needs_spectral_mapping",
]
