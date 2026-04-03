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
    from collections.abc import Callable, Mapping, Sequence

    from siac.domain import SensorBand


logger = logging.getLogger(__name__)


Float32Array: TypeAlias = npt.NDArray[np.float32]
Float64Array: TypeAlias = npt.NDArray[np.float64]
BoolArray: TypeAlias = npt.NDArray[np.bool_]
Int32Array: TypeAlias = npt.NDArray[np.int32]
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


class _PackageSchemaBandProtocol(Protocol):
    band_id: str
    segment: str


class _PackageSensorSchemaProtocol(Protocol):
    bands: Sequence[_PackageSchemaBandProtocol]


class _PackageResultProtocol(Protocol):
    target_reflectance: Sequence[float] | None
    target_band_ids: Sequence[str]
    diagnostics: Mapping[str, object]


class _PackageBatchProtocol(Protocol):
    results: Sequence[_PackageResultProtocol]


class _PackageRetrievalProtocol(Protocol):
    success: bool
    reconstructed: object | None
    neighbor_ids: Sequence[object]
    neighbor_weights: object
    query_valid_mask: object
    source_fit_rmse: float


class _PackageMapperProtocol(Protocol):
    _row_index_by_id: Mapping[str, int]

    def get_sensor_schema(self, sensor_id: str) -> _PackageSensorSchemaProtocol: ...

    def map_reflectance_batch(
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
    ) -> _PackageBatchProtocol: ...

    def _load_hyperspectral(self, segment: str) -> object: ...

    def _band_response(
        self,
        sensor_id: str,
        schema_band: _PackageSchemaBandProtocol,
        *,
        segment_only: bool = False,
    ) -> object: ...


class _PackageMapperMinimalProtocol(_PackageMapperProtocol, Protocol):
    def _candidate_rows(self, candidate_rows: object | None) -> object: ...

    def _load_source_matrix(self, source_sensor: str, segment: str) -> np.ndarray: ...

    def _load_hyperspectral(self, segment: str) -> np.ndarray: ...

    def _retrieve_segment_batch(
        self,
        *,
        source_sensor: str,
        segment: str,
        query_values: Float64Array,
        valid_mask: BoolArray,
        k: int,
        min_valid_bands: int,
        neighbor_estimator: str,
        knn_backend: str,
        knn_eps: float,
        candidate_row_indices: Int64Array,
    ) -> Sequence[_PackageRetrievalProtocol]: ...


class _PackageSpectralMapperFactory(Protocol):
    def __call__(self, prepared_root: Path, *, verify_checksums: bool = False) -> _PackageMapperProtocol: ...


class _PackagePrepareMappingLibrary(Protocol):
    def __call__(
        self,
        *,
        siac_root: Path,
        srf_root: Path,
        output_root: Path,
        source_sensors: Sequence[str],
    ) -> object: ...


class _SpectralLibraryModuleProtocol(Protocol):
    SpectralMapper: _PackageSpectralMapperFactory
    prepare_mapping_library: _PackagePrepareMappingLibrary


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
def _spectral_library_module() -> _SpectralLibraryModuleProtocol:
    return cast(
        "_SpectralLibraryModuleProtocol",
        importlib.import_module("spectral_library"),
    )


def PackageSpectralMapper(*args: object, **kwargs: object) -> _PackageMapperProtocol:
    factory = cast("Callable[..., object]", _spectral_library_module().SpectralMapper)
    return cast("_PackageMapperProtocol", factory(*args, **kwargs))


def prepare_package_mapping_library(*args: object, **kwargs: object) -> object:
    factory = cast("Callable[..., object]", _spectral_library_module().prepare_mapping_library)
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
class _PreparedSegmentDiagnostics:
    neighbor_ids: tuple[str, ...]
    neighbor_weights: Float64Array
    query_valid_mask: BoolArray
    source_fit_rmse: float


@dataclass(frozen=True)
class _DeduplicatedQueries:
    queries: Float64Array
    valid_masks: BoolArray
    inverse_indices: Int64Array


@dataclass(frozen=True)
class _MinimalMappingResult:
    target_reflectance: Float64Array | None
    target_band_ids: tuple[str, ...]
    prepared_segments: dict[str, _PreparedSegmentDiagnostics]


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
    """Package-backed mapper from one multispectral basis to another."""

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
        self._package_mapper: _PackageMapperProtocol | None = None
        self._target_internal_to_output_index: dict[str, int] = {}
        self._source_retrieval_indices_by_segment: dict[str, Int32Array] = {}
        self._target_schema_by_band_id: dict[str, _PackageSchemaBandProtocol] = {}
        self._segment_hyperspectral_cache: dict[str, Float64Array] = {}
        self._target_projection_response_cache: dict[str, Float64Array | None] = {}
        self._segment_target_projection_cache: dict[str, tuple[tuple[str, ...], Int32Array, Float64Array] | None] = {}

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
            runtime_target_schema = self._package_mapper.get_sensor_schema(self._runtime.target_sensor_id)
            self._target_schema_by_band_id = {
                band.band_id: band for band in runtime_target_schema.bands
            }
            runtime_source_schema = self._package_mapper.get_sensor_schema(self._runtime.source_sensor_id)
            for segment in ("vnir", "swir"):
                self._source_retrieval_indices_by_segment[segment] = np.array(
                    [
                        index
                        for index, schema_band in enumerate(runtime_source_schema.bands)
                        if schema_band.segment == segment or (segment == "swir" and schema_band.band_id == "nir")
                    ],
                    dtype=np.int32,
                )

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
        target_flat = np.full((flattened.flat_values.shape[0], len(self.target_bands)), np.nan, dtype=np.float32)
        target_unc_flat = np.full_like(target_flat, np.nan)
        source_fit_flat = np.full(flattened.flat_values.shape[0], np.nan, dtype=np.float32)
        segment_fit_flat = {
            segment: np.full(flattened.flat_values.shape[0], np.nan, dtype=np.float32)
            for segment in ("vnir", "swir")
        }
        if np.any(flattened.valid_rows):
            valid_count = int(np.count_nonzero(flattened.valid_rows))
            logger.info(
                "Spectral mapping start: valid_queries=%d total_queries=%d source_bands=%d target_bands=%d k=%d",
                valid_count,
                int(flattened.flat_values.shape[0]),
                len(self.source_bands),
                len(self.target_bands),
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
            batch_results = self._map_reflectance_batch(
                deduplicated.queries,
                deduplicated.valid_masks,
            )
            logger.info(
                "Spectral mapping package batch complete: queried_rows=%d elapsed=%.2fs",
                unique_query_count,
                perf_counter() - batch_started,
            )

            valid_indices = np.flatnonzero(flattened.valid_rows)
            uncertainty_started = perf_counter()

            # Build band-id → output-column mapping once.
            band_id_to_col: dict[str, int] = self._target_internal_to_output_index

            # Vectorized reflectance fill: gather target reflectance from
            # deduplicated batch results via inverse indices.
            for result_index, row_index in enumerate(valid_indices):
                result = batch_results[int(deduplicated.inverse_indices[result_index])]
                if result.target_reflectance is not None:
                    for band_id, value in zip(result.target_band_ids, result.target_reflectance, strict=True):
                        target_flat[row_index, band_id_to_col[band_id]] = np.float32(value)
                source_fit_flat[row_index] = np.float32(self._estimate_source_fit_rmse(result))
                for segment in ("vnir", "swir"):
                    segment_rmse = self._segment_source_fit_rmse(result, segment)
                    if segment_rmse is not None:
                        segment_fit_flat[segment][row_index] = np.float32(segment_rmse)

            logger.info(
                "Spectral mapping reflectance fill complete: %d pixels elapsed=%.2fs",
                int(valid_indices.size),
                perf_counter() - uncertainty_started,
            )

            # Uncertainty estimation — use vectorized path when possible.
            unc_started = perf_counter()
            for result_index, row_index in enumerate(valid_indices):
                result = batch_results[int(deduplicated.inverse_indices[result_index])]
                target_unc_flat[row_index] = self._estimate_uncertainty(
                    result,
                    source_uncertainty=None if unc_flat is None else unc_flat[row_index],
                )
            logger.info(
                "Spectral mapping uncertainty complete: %d pixels elapsed=%.2fs",
                int(valid_indices.size),
                perf_counter() - unc_started,
            )

        reflectance_da = self._restore_target_cube(target_flat, flattened)
        uncertainty_da = self._restore_target_cube(target_unc_flat, flattened)
        source_fit_da = self._restore_spatial_field(source_fit_flat, flattened)
        segment_fit_da = {
            segment: self._restore_spatial_field(flat_values, flattened)
            for segment, flat_values in segment_fit_flat.items()
        }
        spatial_reference = cast("xr.DataArray", flattened.source_data.isel(band=0, drop=True))
        source_fit_da = copy_spatial_metadata_like(source_fit_da.astype(np.float32), spatial_reference)
        segment_fit_da = {
            segment: copy_spatial_metadata_like(metric.astype(np.float32), spatial_reference)
            for segment, metric in segment_fit_da.items()
        }
        self._cache_distance_metrics(
            {
                "source_fit_rmse": source_fit_da,
                "vnir_source_fit_rmse": segment_fit_da["vnir"],
                "swir_source_fit_rmse": segment_fit_da["swir"],
            }
        )
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

    def _map_reflectance_batch(
        self,
        queries: Float64Array,
        valid_masks: BoolArray,
    ) -> tuple[_PackageResultProtocol | _MinimalMappingResult, ...]:
        if self._package_mapper is None or self._runtime is None:
            raise RuntimeError("spectral-library runtime was not initialized")
        # Fast vectorized path: build cKDTree once, query all rows in
        # parallel, do distance-weighted reconstruction as batch matmuls.
        if self._supports_minimal_batch_path():
            try:
                logger.info(
                    "Spectral mapping batch path: fast vectorized KNN for %d unique row(s)",
                    int(queries.shape[0]),
                )
                return self._map_reflectance_batch_fast(queries, valid_masks)
            except Exception as exc:  # pragma: no cover - compatibility fallback
                logger.warning("Spectral mapping fast batch path failed; falling back to package API (%s)", exc)
        batch = self._package_mapper.map_reflectance_batch(
            source_sensor=self._runtime.source_sensor_id,
            reflectance_rows=queries,
            valid_mask_rows=valid_masks,
            output_mode="target_sensor",
            target_sensor=self._runtime.target_sensor_id,
            k=self.k_neighbors,
            min_valid_bands=self._mapping_config.min_valid_bands,
            neighbor_estimator=self._mapping_config.neighbor_estimator,
            knn_backend=self._mapping_config.knn_backend,
            knn_eps=self._mapping_config.knn_eps,
        )
        return tuple(batch.results)

    def _supports_minimal_batch_path(self) -> bool:
        if self._package_mapper is None:
            return False
        required = (
            "_candidate_rows",
            "_retrieve_segment_batch",
            "_load_source_matrix",
            "_load_hyperspectral",
        )
        return all(callable(getattr(self._package_mapper, name, None)) for name in required)

    def _map_reflectance_batch_fast(
        self,
        queries: Float64Array,
        valid_masks: BoolArray,
    ) -> tuple[_MinimalMappingResult, ...]:
        """Fully vectorized KNN + distance-weighted reconstruction.

        Builds one ``cKDTree`` per segment (instead of per valid-band pattern),
        queries all rows in parallel, and computes the distance-weighted
        reconstruction as batch matrix multiplies.
        """
        from scipy.spatial import cKDTree

        if self._package_mapper is None or self._runtime is None:
            raise RuntimeError("spectral-library runtime was not initialized")

        package_mapper = cast("_PackageMapperMinimalProtocol", self._package_mapper)
        candidate_rows = np.asarray(package_mapper._candidate_rows(None), dtype=np.int64)  # noqa: SLF001
        k = min(self.k_neighbors, int(candidate_rows.size))
        n_queries = int(queries.shape[0])
        estimator = self._mapping_config.neighbor_estimator

        # Per-segment: build tree once → batch query → vectorized reconstruction.
        segment_reconstructed: dict[str, Float64Array] = {}
        segment_source_fit: dict[str, Float32Array] = {}

        for segment in ("vnir", "swir"):
            t0 = perf_counter()
            segment_indices = self._source_retrieval_indices_by_segment[segment]
            if segment_indices.size == 0:
                segment_reconstructed[segment] = np.full((n_queries, 1), np.nan, dtype=np.float64)
                segment_source_fit[segment] = np.full(n_queries, np.nan, dtype=np.float32)
                continue
            query_segment = np.asarray(queries[:, segment_indices], dtype=np.float64)
            mask_segment = np.asarray(valid_masks[:, segment_indices], dtype=np.bool_)

            # Load candidate source matrix and hyperspectral for this segment.
            source_matrix = np.asarray(
                package_mapper._load_source_matrix(self._runtime.source_sensor_id, segment),  # noqa: SLF001
                dtype=np.float64,
            )
            candidate_matrix = source_matrix[candidate_rows]
            hyperspectral = np.asarray(
                package_mapper._load_hyperspectral(segment),  # noqa: SLF001
                dtype=np.float64,
            )

            # Per-row validity: which bands are finite for each query.
            row_valid = np.all(mask_segment, axis=1)
            min_valid = self._mapping_config.min_valid_bands
            enough_bands = mask_segment.sum(axis=1) >= min_valid

            # Build cKDTree on full candidate matrix (all bands).
            # For rows where some bands are NaN, we query on valid-band
            # subsets only if all rows share the same pattern; otherwise
            # we use the full-band tree and mask NaN columns to zero.
            all_valid = np.all(row_valid)

            reconstructed = np.full((n_queries, hyperspectral.shape[1]), np.nan, dtype=np.float64)
            source_fit = np.full(n_queries, np.nan, dtype=np.float32)

            if all_valid:
                # Fast path: all queries have all bands valid → single tree.
                tree = cKDTree(candidate_matrix)
                distances, local_indices = tree.query(query_segment, k=k, workers=-1)
                if distances.ndim == 1:
                    distances = distances[:, np.newaxis]
                    local_indices = local_indices[:, np.newaxis]
                global_indices = candidate_rows[local_indices]

                # Vectorized distance-weighted reconstruction.
                valid_query_mask = enough_bands
                if np.any(valid_query_mask):
                    recon, fit = self._vectorized_dwm_reconstruction(
                        hyperspectral, candidate_matrix, query_segment,
                        global_indices, distances, valid_query_mask, estimator,
                    )
                    reconstructed[valid_query_mask] = recon
                    source_fit[valid_query_mask] = fit
            else:
                # Group by valid-band pattern, build one tree per pattern.
                patterns, inverse = np.unique(mask_segment, axis=0, return_inverse=True)
                for pat_idx, pattern in enumerate(patterns):
                    group_mask = (inverse == pat_idx) & enough_bands
                    if not np.any(group_mask):
                        continue
                    group_indices = np.flatnonzero(group_mask)
                    valid_cols = np.flatnonzero(pattern)
                    cand_sub = candidate_matrix[:, valid_cols]
                    query_sub = query_segment[group_indices][:, valid_cols]

                    tree = cKDTree(cand_sub)
                    distances, local_indices = tree.query(query_sub, k=k, workers=-1)
                    if distances.ndim == 1:
                        distances = distances[:, np.newaxis]
                        local_indices = local_indices[:, np.newaxis]
                    global_indices = candidate_rows[local_indices]

                    recon, fit = self._vectorized_dwm_reconstruction(
                        hyperspectral, candidate_matrix, query_segment[group_indices],
                        global_indices, distances, np.ones(len(group_indices), dtype=bool), estimator,
                    )
                    reconstructed[group_indices] = recon
                    source_fit[group_indices] = fit

            segment_reconstructed[segment] = reconstructed
            segment_source_fit[segment] = source_fit
            logger.info(
                "Spectral mapping segment %s KNN+reconstruction: %.1fs (%d queries, k=%d)",
                segment, perf_counter() - t0, n_queries, k,
            )

        # Project reconstructed hyperspectral → target sensor bands.
        all_band_ids: list[str] = []
        all_target_columns: list[np.ndarray] = []
        for segment in ("vnir", "swir"):
            projection = self._segment_target_projection(segment)
            if projection is None:
                continue
            segment_band_ids, _target_positions, response_matrix = projection
            recon = segment_reconstructed[segment]
            valid_rows = np.all(np.isfinite(recon), axis=1)
            projected = np.full((n_queries, response_matrix.shape[1]), np.nan, dtype=np.float64)
            if np.any(valid_rows):
                projected[valid_rows] = recon[valid_rows] @ response_matrix
            all_band_ids.extend(str(b) for b in segment_band_ids)
            all_target_columns.append(projected)

        target_matrix = np.concatenate(all_target_columns, axis=1) if all_target_columns else None
        target_band_ids_tuple = tuple(all_band_ids)

        # Build lightweight result objects (no per-row _prepare_segment_retrieval).
        results: list[_MinimalMappingResult] = []
        for i in range(n_queries):
            target_ref = target_matrix[i] if target_matrix is not None and np.all(np.isfinite(target_matrix[i])) else None
            results.append(_MinimalMappingResult(
                target_reflectance=target_ref,
                target_band_ids=target_band_ids_tuple,
                prepared_segments={},
            ))
        return tuple(results)

    @staticmethod
    def _vectorized_dwm_reconstruction(
        hyperspectral: Float64Array,
        candidate_matrix: Float64Array,
        query_values: Float64Array,
        global_indices: np.ndarray,
        distances: np.ndarray,
        valid_mask: np.ndarray,
        estimator: str,
    ) -> tuple[Float64Array, Float32Array]:
        """Vectorized distance-weighted-mean reconstruction for a batch of queries.

        Returns (reconstructed_spectra, source_fit_rmse) arrays.
        """
        n = int(np.count_nonzero(valid_mask))
        n_wl = hyperspectral.shape[1]
        gi = global_indices[valid_mask]  # (n, k)
        dist = distances[valid_mask]  # (n, k)
        k = gi.shape[1]

        # Compute weights.
        if estimator == "mean":
            weights = np.full((n, k), 1.0 / k, dtype=np.float64)
        else:
            # distance_weighted_mean (default).
            exact = dist <= 1e-12
            has_exact = np.any(exact, axis=1)
            inv_dist = np.where(dist > 1e-12, 1.0 / dist, 0.0)
            weights = np.where(
                has_exact[:, np.newaxis],
                exact.astype(np.float64),
                inv_dist,
            )
            weight_sums = weights.sum(axis=1, keepdims=True)
            weight_sums = np.where(weight_sums > 0, weight_sums, 1.0)
            weights /= weight_sums  # (n, k)

        # Gather neighbor hyperspectral: (n, k, n_wl)
        neighbor_spectra = hyperspectral[gi.ravel()].reshape(n, k, n_wl)

        # Weighted reconstruction: (n, k) · (n, k, n_wl) → (n, n_wl)
        reconstructed = np.einsum("nk,nkw->nw", weights, neighbor_spectra)

        # Source fit RMSE: reconstruct in source-band space.
        n_src = candidate_matrix.shape[1]
        neighbor_source = candidate_matrix[gi.ravel() % candidate_matrix.shape[0]].reshape(n, k, n_src)
        source_pred = np.einsum("nk,nks->ns", weights, neighbor_source)
        query_valid = query_values[valid_mask]
        source_fit_rmse = np.sqrt(np.nanmean((source_pred - query_valid) ** 2, axis=1)).astype(np.float32)

        return reconstructed, source_fit_rmse

    def _map_reflectance_batch_minimal(
        self,
        queries: Float64Array,
        valid_masks: BoolArray,
    ) -> tuple[_MinimalMappingResult, ...]:
        if self._package_mapper is None or self._runtime is None:
            raise RuntimeError("spectral-library runtime was not initialized")

        package_mapper = cast("_PackageMapperMinimalProtocol", self._package_mapper)
        candidate_rows = np.asarray(package_mapper._candidate_rows(None), dtype=np.int64)  # noqa: SLF001
        retrievals_by_segment: dict[str, tuple[_PackageRetrievalProtocol, ...]] = {}
        for segment in ("vnir", "swir"):
            segment_indices = self._source_retrieval_indices_by_segment[segment]
            retrievals_by_segment[segment] = tuple(
                package_mapper._retrieve_segment_batch(  # noqa: SLF001
                    source_sensor=self._runtime.source_sensor_id,
                    segment=segment,
                    query_values=np.asarray(queries[:, segment_indices], dtype=np.float64),
                    valid_mask=np.asarray(valid_masks[:, segment_indices], dtype=np.bool_),
                    k=self.k_neighbors,
                    min_valid_bands=self._mapping_config.min_valid_bands,
                    neighbor_estimator=self._mapping_config.neighbor_estimator,
                    knn_backend=self._mapping_config.knn_backend,
                    knn_eps=self._mapping_config.knn_eps,
                    candidate_row_indices=candidate_rows,
                )
            )

        # Vectorized post-processing: extract reconstructed spectra into
        # contiguous arrays and apply response-matrix projection as a single
        # batch matrix multiply instead of a per-row Python loop.
        n_queries = int(queries.shape[0])
        all_band_ids: list[str] = []
        all_target_columns: list[np.ndarray] = []

        for segment in ("vnir", "swir"):
            projection = self._segment_target_projection(segment)
            if projection is None:
                continue
            segment_band_ids, _target_positions, response_matrix = projection
            retrievals = retrievals_by_segment[segment]

            # Gather reconstructed spectra into (n_queries, n_wavelengths) matrix.
            n_wl = int(response_matrix.shape[0])
            reconstructed_matrix = np.full((n_queries, n_wl), np.nan, dtype=np.float64)
            for i, retrieval in enumerate(retrievals):
                if retrieval.success and retrieval.reconstructed is not None:
                    reconstructed_matrix[i] = np.asarray(retrieval.reconstructed, dtype=np.float64)

            # Batch projection: (n_queries, n_wl) @ (n_wl, n_bands) → (n_queries, n_bands)
            valid_rows = np.all(np.isfinite(reconstructed_matrix), axis=1)
            projected = np.full((n_queries, response_matrix.shape[1]), np.nan, dtype=np.float64)
            if np.any(valid_rows):
                projected[valid_rows] = reconstructed_matrix[valid_rows] @ response_matrix

            all_band_ids.extend(str(b) for b in segment_band_ids)
            all_target_columns.append(projected)

        # Build results with pre-computed target reflectance.
        if all_target_columns:
            target_matrix = np.concatenate(all_target_columns, axis=1)
            target_band_ids_tuple = tuple(all_band_ids)
        else:
            target_matrix = None
            target_band_ids_tuple = ()

        results: list[_MinimalMappingResult] = []
        for sample_index in range(n_queries):
            prepared_segments: dict[str, _PreparedSegmentDiagnostics] = {}
            for segment in ("vnir", "swir"):
                retrieval = retrievals_by_segment[segment][sample_index]
                prepared = self._prepare_segment_retrieval(retrieval)
                if prepared is not None:
                    prepared_segments[segment] = prepared

            if target_matrix is not None and np.all(np.isfinite(target_matrix[sample_index])):
                target_reflectance = target_matrix[sample_index]
            else:
                target_reflectance = None

            results.append(
                _MinimalMappingResult(
                    target_reflectance=target_reflectance,
                    target_band_ids=target_band_ids_tuple,
                    prepared_segments=prepared_segments,
                )
            )
        return tuple(results)

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

    def _segment_hyperspectral_cache_store(self) -> dict[str, Float64Array]:
        cache = getattr(self, "_segment_hyperspectral_cache", None)
        if cache is None:
            cache = {}
            self._segment_hyperspectral_cache = cache
        return cast("dict[str, Float64Array]", cache)

    def _target_projection_response_cache_store(self) -> dict[str, Float64Array | None]:
        cache = getattr(self, "_target_projection_response_cache", None)
        if cache is None:
            cache = {}
            self._target_projection_response_cache = cache
        return cast("dict[str, Float64Array | None]", cache)

    def _segment_target_projection_cache_store(
        self,
    ) -> dict[str, tuple[tuple[str, ...], Int32Array, Float64Array] | None]:
        cache = getattr(self, "_segment_target_projection_cache", None)
        if cache is None:
            cache = {}
            self._segment_target_projection_cache = cache
        return cast("dict[str, tuple[tuple[str, ...], Int32Array, Float64Array] | None]", cache)

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

    @staticmethod
    def _prepared_segment_from_result(
        result: _PackageResultProtocol | _MinimalMappingResult,
        segment: str,
    ) -> _PreparedSegmentDiagnostics | None:
        if isinstance(result, _MinimalMappingResult):
            return result.prepared_segments.get(segment)
        diagnostics = result.diagnostics.get("segments", {})
        if not isinstance(diagnostics, dict):
            diagnostics = {}
        return SpectralMapper._prepare_segment_diagnostics(diagnostics.get(segment))

    def _estimate_source_fit_rmse(
        self,
        result: _PackageResultProtocol | _MinimalMappingResult,
    ) -> float:
        segment_rmse: list[tuple[int, float]] = []
        for segment in ("vnir", "swir"):
            prepared = self._prepared_segment_from_result(result, segment)
            if prepared is None:
                continue
            valid_count = int(np.count_nonzero(prepared.query_valid_mask))
            if valid_count < 1:
                continue
            segment_rmse.append((valid_count, max(0.0, float(prepared.source_fit_rmse))))
        if not segment_rmse:
            return 0.0
        weights = np.asarray([item[0] for item in segment_rmse], dtype=np.float64)
        rmse = np.asarray([item[1] for item in segment_rmse], dtype=np.float64)
        return float(np.sqrt(np.average(np.square(rmse, dtype=np.float64), weights=weights)))

    def _segment_source_fit_rmse(
        self,
        result: _PackageResultProtocol | _MinimalMappingResult,
        segment: str,
    ) -> float | None:
        prepared = self._prepared_segment_from_result(result, segment)
        if prepared is None:
            return None
        if int(np.count_nonzero(prepared.query_valid_mask)) < 1:
            return None
        return max(0.0, float(prepared.source_fit_rmse))

    def _estimate_uncertainty(
        self,
        result: _PackageResultProtocol | _MinimalMappingResult,
        *,
        source_uncertainty: Float32Array | None,
    ) -> Float32Array:
        if self._package_mapper is None or self._runtime is None:
            raise RuntimeError("spectral-library runtime was not initialized")

        output: Float32Array = np.full(len(self.target_bands), np.nan, dtype=np.float32)
        target_values_by_band_id = self._target_values_by_band_id(result)
        for segment in ("vnir", "swir"):
            prepared = self._prepared_segment_from_result(result, segment)
            if prepared is None:
                continue
            neighbor_spectra = self._load_neighbor_spectra(
                prepared.neighbor_ids,
                segment,
            )
            neighbor_weights = prepared.neighbor_weights
            input_unc = self._segment_input_uncertainty(
                segment,
                prepared.query_valid_mask,
                source_uncertainty,
            )
            segment_projection = self._segment_target_projection(segment)
            if segment_projection is None:
                continue
            band_ids, target_positions, response_matrix = segment_projection
            neighbor_target_matrix = np.asarray(neighbor_spectra @ response_matrix, dtype=np.float64)

            for column_index, (band_id, target_index) in enumerate(zip(band_ids, target_positions, strict=True)):
                neighbor_target = neighbor_target_matrix[:, column_index]
                estimate = target_values_by_band_id.get(band_id)
                if estimate is None:
                    estimate = float(np.dot(neighbor_weights, neighbor_target))
                spread = float(
                    np.sqrt(np.dot(neighbor_weights, np.square(neighbor_target - estimate, dtype=np.float64)))
                )
                output[int(target_index)] = np.float32(
                    np.sqrt(max(spread, _UNCERTAINTY_FLOOR) ** 2 + prepared.source_fit_rmse**2 + input_unc**2)
                )

        mapped_positions = [
            self._target_internal_to_output_index[band_id]
            for band_id in target_values_by_band_id
            if band_id in self._target_internal_to_output_index
        ]
        for position in mapped_positions:
            if not np.isfinite(output[position]):
                output[position] = _UNCERTAINTY_FLOOR
        return output

    @staticmethod
    def _target_values_by_band_id(result: _PackageResultProtocol | _MinimalMappingResult) -> dict[str, float]:
        if result.target_reflectance is None:
            return {}
        return {
            str(band_id): float(value)
            for band_id, value in zip(result.target_band_ids, result.target_reflectance, strict=True)
        }

    @staticmethod
    def _prepare_segment_diagnostics(segment_diag: object) -> _PreparedSegmentDiagnostics | None:
        if not isinstance(segment_diag, dict) or segment_diag.get("status") != "ok":
            return None
        neighbor_ids = tuple(str(value) for value in segment_diag.get("neighbor_ids", ()))
        if not neighbor_ids:
            return None
        neighbor_weights = np.asarray(segment_diag.get("neighbor_weights", ()), dtype=np.float64)
        if neighbor_weights.size == 0:
            return None
        weight_sum = float(np.sum(neighbor_weights))
        if weight_sum <= 0.0:
            return None
        normalized_weights = np.asarray(neighbor_weights / weight_sum, dtype=np.float64)
        query_valid_mask = np.asarray(segment_diag.get("query_valid_mask", ()), dtype=np.bool_)
        source_fit_rmse = float(segment_diag.get("source_fit_rmse") or 0.0)
        return _PreparedSegmentDiagnostics(
            neighbor_ids=neighbor_ids,
            neighbor_weights=normalized_weights,
            query_valid_mask=query_valid_mask,
            source_fit_rmse=source_fit_rmse,
        )

    @staticmethod
    def _prepare_segment_retrieval(retrieval: _PackageRetrievalProtocol) -> _PreparedSegmentDiagnostics | None:
        if not retrieval.success:
            return None
        neighbor_ids = tuple(str(value) for value in retrieval.neighbor_ids)
        if not neighbor_ids:
            return None
        neighbor_weights = np.asarray(retrieval.neighbor_weights, dtype=np.float64)
        if neighbor_weights.size == 0:
            return None
        weight_sum = float(np.sum(neighbor_weights))
        if weight_sum <= 0.0:
            return None
        normalized_weights = np.asarray(neighbor_weights / weight_sum, dtype=np.float64)
        query_valid_mask = np.asarray(retrieval.query_valid_mask, dtype=np.bool_)
        source_fit_rmse = float(retrieval.source_fit_rmse or 0.0)
        return _PreparedSegmentDiagnostics(
            neighbor_ids=neighbor_ids,
            neighbor_weights=normalized_weights,
            query_valid_mask=query_valid_mask,
            source_fit_rmse=source_fit_rmse,
        )

    def _simulate_target_sensor_outputs(
        self,
        segment_outputs: dict[str, Float64Array],
    ) -> tuple[Float64Array | None, tuple[str, ...]]:
        band_ids: list[str] = []
        values: list[float] = []
        for segment in ("vnir", "swir"):
            reconstructed = segment_outputs.get(segment)
            if reconstructed is None:
                continue
            projection = self._segment_target_projection(segment)
            if projection is None:
                continue
            segment_band_ids, _target_positions, response_matrix = projection
            projected = np.asarray(reconstructed @ response_matrix, dtype=np.float64)
            for band_id, value in zip(segment_band_ids, projected, strict=True):
                band_ids.append(str(band_id))
                values.append(float(value))
        if not values:
            return None, ()
        return np.asarray(values, dtype=np.float64), tuple(band_ids)

    def _load_neighbor_spectra(
        self,
        neighbor_ids: tuple[str, ...],
        segment: str,
    ) -> Float64Array:
        if self._package_mapper is None:
            raise RuntimeError("spectral-library runtime was not initialized")
        row_indices: Int64Array = np.asarray(
            [self._package_mapper._row_index_by_id[row_id] for row_id in neighbor_ids],  # noqa: SLF001
            dtype=np.int64,
        )
        cache = self._segment_hyperspectral_cache_store()
        neighbor_spectra = cache.get(segment)
        if neighbor_spectra is None:
            neighbor_spectra = np.asarray(
                self._package_mapper._load_hyperspectral(segment),  # noqa: SLF001
                dtype=np.float64,
            )
            cache[segment] = neighbor_spectra
        return np.asarray(neighbor_spectra[row_indices], dtype=np.float64)

    def _segment_input_uncertainty(
        self,
        segment: str,
        query_valid_mask: BoolArray,
        source_uncertainty: Float32Array | None,
    ) -> float:
        if source_uncertainty is None:
            return 0.0
        retrieval_indices = self._source_retrieval_indices_by_segment[segment]
        if query_valid_mask.size == 0:
            return 0.0
        valid_indices = retrieval_indices[query_valid_mask]
        if valid_indices.size == 0:
            return 0.0
        return float(
            np.sqrt(np.mean(np.square(source_uncertainty[valid_indices], dtype=np.float64)))
        )

    def _neighbor_target_projection(
        self,
        neighbor_spectra: Float64Array,
        target_schema_band: _PackageSchemaBandProtocol,
    ) -> Float64Array | None:
        if self._package_mapper is None or self._runtime is None:
            raise RuntimeError("spectral-library runtime was not initialized")
        response = self._target_projection_response(getattr(target_schema_band, "band_id", ""), target_schema_band)
        if response is None:
            return None
        return np.asarray(neighbor_spectra @ response, dtype=np.float64)

    def _target_projection_response(
        self,
        band_id: str,
        target_schema_band: _PackageSchemaBandProtocol,
    ) -> Float64Array | None:
        if self._package_mapper is None or self._runtime is None:
            raise RuntimeError("spectral-library runtime was not initialized")
        if not band_id:
            response = np.asarray(
                self._package_mapper._band_response(  # noqa: SLF001
                    self._runtime.target_sensor_id,
                    target_schema_band,
                    segment_only=True,
                ),
                dtype=np.float64,
            )
            denominator = float(np.sum(response))
            return None if denominator <= 0.0 else np.asarray(response / denominator, dtype=np.float64)
        cache = self._target_projection_response_cache_store()
        if band_id in cache:
            return cache[band_id]
        response = np.asarray(
            self._package_mapper._band_response(  # noqa: SLF001
                self._runtime.target_sensor_id,
                target_schema_band,
                segment_only=True,
            ),
            dtype=np.float64,
        )
        denominator = float(np.sum(response))
        normalized = None if denominator <= 0.0 else np.asarray(response / denominator, dtype=np.float64)
        cache[band_id] = normalized
        return normalized

    def _segment_target_projection(
        self,
        segment: str,
    ) -> tuple[tuple[str, ...], Int32Array, Float64Array] | None:
        cache = self._segment_target_projection_cache_store()
        if segment in cache:
            return cache[segment]

        entries: list[tuple[str, int, Float64Array]] = []
        for band_id, target_schema_band in self._target_schema_by_band_id.items():
            if target_schema_band.segment != segment:
                continue
            target_index = self._target_internal_to_output_index.get(band_id)
            if target_index is None:
                continue
            response = self._target_projection_response(band_id, target_schema_band)
            if response is None:
                continue
            entries.append((band_id, int(target_index), response))

        if not entries:
            cache[segment] = None
            return None

        band_ids = tuple(entry[0] for entry in entries)
        target_positions = np.asarray([entry[1] for entry in entries], dtype=np.int32)
        response_matrix = np.stack([entry[2] for entry in entries], axis=1).astype(np.float64, copy=False)
        cached = (band_ids, target_positions, response_matrix)
        cache[segment] = cached
        return cached


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
