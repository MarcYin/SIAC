"""Spectral mapping between differing source and target sensor band sets."""

from __future__ import annotations

import csv
import hashlib
import json
import logging
import os
from dataclasses import dataclass
from pathlib import Path
from time import perf_counter
from typing import TYPE_CHECKING, Protocol, TypeAlias, cast

import numpy as np
import xarray as xr
from numpy import typing as npt
from spectral_library import (
    BandInput,
    SensorInput,
    build_mapping_runtime,
)
from spectral_library import (
    HyperspectralLibraryInput as PackageHyperspectralLibraryInput,
)
from spectral_library import (
    PreparedRuntime as PackagePreparedRuntime,
)

from siac.algorithms.surface import _spectral_curve_utils as curve_utils
from siac.algorithms.surface._spectral_curve_utils import (
    normalized_band_response as _normalized_band_response,
)
from siac.algorithms.surface._spectral_curve_utils import (
    primary_nir_band_index as _primary_nir_band_index,
)
from siac.runtime.models import copy_spatial_metadata_like
from siac.storage.writers import write_raster

if TYPE_CHECKING:
    from collections.abc import Sequence

    from siac.domain import SensorBand


logger = logging.getLogger(__name__)


Float32Array: TypeAlias = npt.NDArray[np.float32]
Float64Array: TypeAlias = npt.NDArray[np.float64]
BoolArray: TypeAlias = npt.NDArray[np.bool_]
Int64Array: TypeAlias = npt.NDArray[np.int64]


def _canonicalize_curve(
    wavelengths_nm: Float32Array,
    response: Float32Array,
) -> tuple[Float32Array, Float32Array]:
    return cast(
        "tuple[Float32Array, Float32Array]",
        curve_utils.canonicalize_curve(wavelengths_nm, response),
    )


def _segmentize_curve(
    wavelengths_nm: Float32Array,
    response: Float32Array,
    *,
    segment: str,
) -> tuple[Float32Array, Float32Array]:
    return cast(
        "tuple[Float32Array, Float32Array]",
        curve_utils.segmentize_curve(wavelengths_nm, response, segment=segment),
    )


class _BatchArrayResultProtocol(Protocol):
    """Protocol matching ``spectral_library.BatchMappingArrayResult`` (v0.5.0+)."""

    reflectance: np.ndarray
    source_fit_rmse: np.ndarray
    output_columns: tuple[str, ...]


class _PackageRuntimeProtocol(Protocol):
    def map_reflectance_batch_arrays_ndarray(
        self,
        *,
        source_sensor: str,
        reflectance_rows: Float64Array,
        valid_mask_rows: BoolArray,
        output_mode: str,
        target_sensor: str,
        k: int,
        min_valid_bands: int,
        neighbor_estimator: str,
        knn_backend: str,
        knn_eps: float,
    ) -> _BatchArrayResultProtocol: ...

_DEFAULT_K_NEIGHBORS = 5
_UNCERTAINTY_FLOOR = 0.005
_CANONICAL_WAVELENGTHS_NM: Float32Array = np.arange(400.0, 2501.0, 1.0, dtype=np.float32)
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
            cache_dir=_path(self.cache_dir),
            neighbor_estimator=str(self.neighbor_estimator).strip() or _DEFAULT_NEIGHBOR_ESTIMATOR,
            knn_backend=str(self.knn_backend).strip() or _DEFAULT_KNN_BACKEND,
            knn_eps=float(self.knn_eps),
            min_valid_bands=max(1, int(self.min_valid_bands)),
        )


@dataclass(frozen=True)
class _PreparedRuntime:
    runtime: PackagePreparedRuntime
    source_sensor_id: str
    target_sensor_id: str
    target_band_ids: tuple[str, ...]

    @property
    def prepared_root(self) -> Path:
        return Path(self.runtime.prepared_root)


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
    source = reflectance.transpose(
        *[dim for dim in reflectance.dims if dim != "wavelength"], "wavelength"
    )
    if source.sizes["wavelength"] != len(wavelengths_nm):
        raise ValueError("wavelength dimension does not match wavelengths_nm")

    matrix = np.stack(
        [
            _normalized_band_response(band, np.asarray(wavelengths_nm, dtype=np.float32))
            for band in target_bands
        ],
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
    return Path(
        os.getenv(
            _SIAC_SPECTRAL_MAPPING_CACHE_DIR_ENV,
            Path.home() / ".cache" / "siac" / "spectral_mapping",
        )
    )


def _hash_bytes(*chunks: bytes) -> str:
    digest = hashlib.sha256()
    for chunk in chunks:
        digest.update(chunk)
    return digest.hexdigest()


def _json_hash(payload: dict[str, object]) -> str:
    return _hash_bytes(json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8"))


def _runtime_band_id(
    band: SensorBand,
    *,
    band_index: int,
    primary_nir_index: int | None,
    used_band_ids: set[str],
) -> str:
    band_id = "nir" if primary_nir_index == band_index else str(band.name)
    if band_id in used_band_ids:
        band_id = f"{band_id}_{band_index}"
    used_band_ids.add(band_id)
    return band_id


def _canonical_hyperspectral_spectra(library: HyperspectralLibrary) -> Float32Array:
    if np.array_equal(library.wavelengths_nm.astype(np.float32), _CANONICAL_WAVELENGTHS_NM):
        return np.asarray(library.spectra, dtype=np.float32)

    wavelengths = np.asarray(library.wavelengths_nm, dtype=np.float32)
    if wavelengths[0] > float(_CANONICAL_WAVELENGTHS_NM[0]) or wavelengths[-1] < float(
        _CANONICAL_WAVELENGTHS_NM[-1]
    ):
        raise ValueError(
            "HyperspectralLibrary must cover the canonical 400-2500 nm range to drive spectral-library mapping"
        )
    return np.vstack(
        [
            np.interp(
                _CANONICAL_WAVELENGTHS_NM,
                wavelengths,
                np.asarray(spectrum, dtype=np.float32),
            ).astype(np.float32)
            for spectrum in np.asarray(library.spectra, dtype=np.float32)
        ]
    )


def _package_library_input_from_hyperspectral_library(
    library: HyperspectralLibrary,
) -> PackageHyperspectralLibraryInput:
    source_version = str(library.source_version).strip()
    metadata_rows = (
        [{"source_version": source_version} for _ in library.sample_ids]
        if source_version
        else None
    )
    return PackageHyperspectralLibraryInput(
        wavelengths_nm=np.asarray(_CANONICAL_WAVELENGTHS_NM, dtype=np.float64),
        spectra=_canonical_hyperspectral_spectra(library),
        sample_ids=tuple(str(sample_id) for sample_id in library.sample_ids),
        metadata_rows=metadata_rows,
        source_id=str(library.source_id),
    )


def _read_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        return [dict(row) for row in reader]


def _package_library_input_from_root(root: Path) -> PackageHyperspectralLibraryInput:
    tabular_root = Path(root).expanduser().resolve() / "tabular"
    metadata_path = tabular_root / "siac_spectra_metadata.csv"
    spectra_path = tabular_root / "siac_normalized_spectra.csv"
    if not metadata_path.exists() or not spectra_path.exists():
        raise ValueError(
            "SIAC spectral library root must contain tabular/siac_spectra_metadata.csv "
            "and tabular/siac_normalized_spectra.csv."
        )

    metadata_rows = _read_csv_rows(metadata_path)
    spectra_rows = _read_csv_rows(spectra_path)
    if not spectra_rows:
        raise ValueError("SIAC spectral library root must contain at least one spectrum row.")
    if metadata_rows and len(metadata_rows) != len(spectra_rows):
        raise ValueError("SIAC spectral library metadata and spectra row counts must match.")

    spectra_fieldnames = list(spectra_rows[0])
    wavelength_fields = [field for field in spectra_fieldnames if field.startswith("nm_")]
    if not wavelength_fields:
        raise ValueError("SIAC spectral library spectra CSV must include nm_<wavelength> columns.")

    wavelengths_nm = np.asarray(
        [float(field.removeprefix("nm_")) for field in wavelength_fields],
        dtype=np.float64,
    )
    spectra = np.asarray(
        [
            [float(row[field]) for field in wavelength_fields]
            for row in spectra_rows
        ],
        dtype=np.float32,
    )

    resolved_metadata_rows = metadata_rows or [
        {
            "source_id": str(row.get("source_id") or "siac-root"),
            "spectrum_id": str(row.get("spectrum_id") or row.get("sample_name") or f"row_{index}"),
            "sample_name": str(row.get("sample_name") or row.get("spectrum_id") or f"row_{index}"),
        }
        for index, row in enumerate(spectra_rows)
    ]
    sample_ids = tuple(
        str(row.get("spectrum_id") or row.get("sample_name") or f"row_{index}")
        for index, row in enumerate(resolved_metadata_rows)
    )
    source_id = str(resolved_metadata_rows[0].get("source_id") or "siac-root")
    return PackageHyperspectralLibraryInput(
        wavelengths_nm=wavelengths_nm,
        spectra=spectra,
        sample_ids=sample_ids,
        metadata_rows=resolved_metadata_rows,
        source_id=source_id,
    )


def _resolved_package_library_input(
    library: HyperspectralLibrary | None,
    config: SpectralMappingConfig,
) -> PackageHyperspectralLibraryInput:
    if library is not None:
        return _package_library_input_from_hyperspectral_library(library)
    if config.siac_library_root is not None:
        return _package_library_input_from_root(Path(config.siac_library_root))
    raise ValueError(
        "Spectral mapping requires an explicit SIAC spectral library. "
        "Provide spectral_library=HyperspectralLibrary(...) or set "
        "SIAC_SPECTRAL_LIBRARY_ROOT / SpectralMappingConfig.siac_library_root."
    )


def _band_input_for_band(
    band: SensorBand,
    *,
    band_id: str,
) -> BandInput:
    if band.has_rsrf:
        return BandInput(
            band_id=band_id,
            center_wavelength_nm=float(band.center_wavelength),
            fwhm_nm=float(band.bandwidth),
            response_definition={
                "kind": "sampled",
                "wavelength_nm": np.asarray(band.rsrf_wavelengths_nm, dtype=np.float64).tolist(),
                "response": np.asarray(band.rsrf_response, dtype=np.float64).tolist(),
            },
        )
    if band.rsrf_sensor_unit_id is not None:
        return BandInput(
            band_id=band_id,
            center_wavelength_nm=float(band.center_wavelength),
            fwhm_nm=float(band.bandwidth),
            rsrf_sensor_id=str(band.rsrf_sensor_unit_id),
            rsrf_band_id=str(band.rsrf_band_id or band.name),
            rsrf_representation_variant=band.rsrf_representation_variant,
        )
    return BandInput(
        band_id=band_id,
        center_wavelength_nm=float(band.center_wavelength),
        fwhm_nm=float(band.bandwidth),
    )


def _sensor_input_for_bands(
    bands: Sequence[SensorBand],
) -> tuple[SensorInput, tuple[str, ...]]:
    primary_nir_index = _primary_nir_band_index(bands)
    used_band_ids: set[str] = set()
    band_ids: list[str] = []
    band_inputs: list[BandInput] = []

    for index, band in enumerate(bands):
        band_id = _runtime_band_id(
            band,
            band_index=index,
            primary_nir_index=primary_nir_index,
            used_band_ids=used_band_ids,
        )
        band_ids.append(band_id)
        band_inputs.append(_band_input_for_band(band, band_id=band_id))

    return SensorInput(bands=tuple(band_inputs)), tuple(band_ids)


def _atomic_write_json(path: Path, payload: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    tmp_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
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
                        np.asarray(
                            np.asarray(array.coords[dim].values).shape, dtype=np.int64
                        ).tobytes(),
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
            logger.warning(
                "Skipping diagnostic GeoTIFF for %s: spatial metadata is incomplete", name
            )
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
    source_sensor_input, _source_band_ids = _sensor_input_for_bands(source_bands)
    target_sensor_input, _target_band_ids = _sensor_input_for_bands(target_bands)
    runtime = build_mapping_runtime(
        library=_resolved_package_library_input(library, normalized_config),
        source_sensors=[source_sensor_input],
        target_sensors=[target_sensor_input],
        cache_root=cache_root,
        verify_checksums=False,
    )
    source_sensor_id = runtime.source_sensor_ids[0]
    target_sensor_id = runtime.target_sensor_ids[0]
    return _PreparedRuntime(
        runtime=runtime,
        source_sensor_id=source_sensor_id,
        target_sensor_id=target_sensor_id,
        target_band_ids=runtime.target_band_ids[target_sensor_id],
    )


class SpectralMapper:
    """Package-backed mapper from one multispectral basis to another.

    Uses ``spectral_library >= 0.6.0`` direct runtime construction plus the
    batch array API for efficient
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
        self._package_mapper: _PackageRuntimeProtocol | None = None
        self._target_internal_to_output_index: dict[str, int] = {}

        if not self._identity:
            self._runtime = _prepare_runtime(
                self.source_bands,
                self.target_bands,
                library=self._spectral_library,
                config=self._mapping_config,
            )
            self._package_mapper = cast("_PackageRuntimeProtocol", self._runtime.runtime)
            self._target_internal_to_output_index = {
                band_id: index for index, band_id in enumerate(self._runtime.target_band_ids)
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
            zero_fit = xr.zeros_like(
                cast("xr.DataArray", source_data.isel(band=0, drop=True)), dtype=np.float32
            )
            return (
                identity.transpose(*original_dims).astype(np.float32),
                unc.transpose(*original_dims).astype(np.float32),
                zero_fit.astype(np.float32),
            )

        if self._package_mapper is None or self._runtime is None:
            raise RuntimeError("spectral-library runtime was not initialized")

        flattened = self._flatten_source_cube(source_data)
        package_queries, package_valid_masks, package_valid_rows = self._package_query_rows(flattened)
        unc_flat = self._flatten_optional_uncertainty(source_unc, flattened.transpose_dims)
        n_pixels = flattened.flat_values.shape[0]
        n_target = len(self.target_bands)
        target_flat = np.full((n_pixels, n_target), np.nan, dtype=np.float32)
        target_unc_flat = np.full_like(target_flat, np.nan)
        source_fit_flat = np.full(n_pixels, np.nan, dtype=np.float32)

        if np.any(package_valid_rows):
            valid_count = int(np.count_nonzero(package_valid_rows))
            logger.info(
                "Spectral mapping start: valid_queries=%d total_queries=%d source_bands=%d target_bands=%d k=%d",
                valid_count,
                n_pixels,
                len(self.source_bands),
                n_target,
                self.k_neighbors,
            )
            batch_started = perf_counter()
            valid_queries = np.asarray(package_queries[package_valid_rows], dtype=np.float64)
            valid_masks = np.asarray(package_valid_masks[package_valid_rows], dtype=np.bool_)
            deduplicated = self._deduplicate_queries(valid_queries, valid_masks)
            unique_query_count = int(deduplicated.queries.shape[0])
            if unique_query_count != valid_count:
                logger.info(
                    "Spectral mapping deduplicated queries: %d -> %d unique rows",
                    valid_count,
                    unique_query_count,
                )

            # spectral_library batch array API returns arrays directly.
            batch_result = self._package_mapper.map_reflectance_batch_arrays_ndarray(
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
                )
            logger.info(
                "Spectral mapping package batch complete: queried_rows=%d elapsed=%.2fs",
                unique_query_count,
                perf_counter() - batch_started,
            )

            valid_indices = np.flatnonzero(package_valid_rows)

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
            target_flat[np.ix_(valid_indices, target_col_arr)] = reflectance[dedup_indices][
                :, result_col_arr
            ]
            source_fit_flat[valid_indices] = fit_rmse[dedup_indices]
            logger.info(
                "Spectral mapping reflectance fill complete: %d pixels elapsed=%.2fs",
                int(valid_indices.size),
                perf_counter() - fill_started,
            )

            # Vectorized uncertainty: sqrt(floor² + fit_rmse² + input_unc²).
            unc_started = perf_counter()
            per_pixel_fit = fit_rmse[dedup_indices].astype(np.float64)
            input_unc: Float64Array = np.zeros(valid_count, dtype=np.float64)
            if unc_flat is not None:
                unc_sq = np.square(unc_flat[valid_indices].astype(np.float64))
                input_unc = np.asarray(np.sqrt(np.nanmean(unc_sq, axis=1)), dtype=np.float64)
            unc_per_pixel = np.sqrt(_UNCERTAINTY_FLOOR**2 + per_pixel_fit**2 + input_unc**2).astype(
                np.float32
            )
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
        source_fit_da = copy_spatial_metadata_like(
            source_fit_da.astype(np.float32), spatial_reference
        )
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
        missing = [
            band.name
            for band in self.source_bands
            if band.name not in set(data.coords["band"].values.tolist())
        ]
        if missing:
            raise KeyError(f"source_reflectance is missing source bands: {missing}")
        return data.sel(band=[band.name for band in self.source_bands])

    def _flatten_source_cube(self, source_data: xr.DataArray) -> _FlattenedSourceCube:
        transpose_dims = ("band", *[dim for dim in source_data.dims if dim != "band"])
        source_values = np.asarray(source_data.transpose(*transpose_dims).values, dtype=np.float32)
        flat_values = np.asarray(
            source_values.reshape(source_values.shape[0], -1).T, dtype=np.float32
        )
        valid_rows = np.asarray(np.any(np.isfinite(flat_values), axis=1), dtype=np.bool_)
        return _FlattenedSourceCube(
            original_dims=tuple(source_data.dims),
            transpose_dims=transpose_dims,
            spatial_dims=transpose_dims[1:],
            source_data=source_data,
            flat_values=flat_values,
            valid_rows=valid_rows,
        )

    def _package_query_rows(
        self,
        flattened: _FlattenedSourceCube,
    ) -> tuple[Float64Array, BoolArray, BoolArray]:
        query_rows = np.asarray(flattened.flat_values, dtype=np.float64)
        valid_masks = np.isfinite(query_rows)
        valid_rows = np.asarray(np.any(valid_masks, axis=1), dtype=np.bool_)
        return query_rows, valid_masks, valid_rows

    @staticmethod
    def _flatten_optional_uncertainty(
        source_uncertainty: xr.DataArray | None,
        transpose_dims: tuple[str, ...],
    ) -> Float32Array | None:
        if source_uncertainty is None:
            return None
        unc_values = np.asarray(
            source_uncertainty.transpose(*transpose_dims).values, dtype=np.float32
        )
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
        shape = (
            len(self.target_bands),
            *[flattened.source_data.sizes[dim] for dim in flattened.spatial_dims],
        )
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
