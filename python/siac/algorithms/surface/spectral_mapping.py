"""Spectral mapping between differing source and target sensor band sets."""

from __future__ import annotations

import csv
import hashlib
import json
import logging
import os
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np
import rsrf
import xarray as xr
from spectral_library import SpectralMapper as PackageSpectralMapper
from spectral_library import prepare_mapping_library as prepare_package_mapping_library

if TYPE_CHECKING:
    from collections.abc import Sequence

    from siac.domain import SensorBand


logger = logging.getLogger(__name__)

_VISIBLE_UPPER_NM = 700.0
_DEFAULT_K_NEIGHBORS = 5
_UNCERTAINTY_FLOOR = 0.005
_CANONICAL_WAVELENGTHS_NM = np.arange(400.0, 2501.0, 1.0, dtype=np.float32)
_RSRF_ROOT_ENV = "RSRF_ROOT"
_SIAC_SPECTRAL_LIBRARY_ROOT_ENV = "SIAC_SPECTRAL_LIBRARY_ROOT"
_SIAC_SPECTRAL_MAPPING_CACHE_DIR_ENV = "SIAC_SPECTRAL_MAPPING_CACHE_DIR"
_DEFAULT_NEIGHBOR_ESTIMATOR = "distance_weighted_mean"
_DEFAULT_KNN_BACKEND = "numpy"
_SEGMENT_RANGES_NM: dict[str, tuple[float, float]] = {
    "vnir": (400.0, 1000.0),
    "swir": (800.0, 2500.0),
}


def _trapezoid(y: np.ndarray, x: np.ndarray) -> float:
    return float(np.trapezoid(y, x))


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

def _normalized_band_response(band: SensorBand, wavelengths_nm: np.ndarray) -> np.ndarray:
    response = np.asarray(band.effective_response(wavelengths_nm), dtype=np.float32)
    area = _trapezoid(response, wavelengths_nm)
    if not np.isfinite(area) or area <= 0.0:
        raise ValueError(f"Band {band.name!r} has zero support on the requested wavelength grid")
    return response / area


def _classify_band_region(band: SensorBand) -> str:
    return "visible" if band.center_wavelength < _VISIBLE_UPPER_NM else "infrared"


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


def _segment_for_band(band: SensorBand) -> str:
    return "swir" if float(band.center_wavelength) >= 1000.0 else "vnir"


def _segmentize_curve(
    wavelengths_nm: np.ndarray,
    response: np.ndarray,
    *,
    segment: str,
) -> tuple[np.ndarray, np.ndarray]:
    segment_min_nm, segment_max_nm = _SEGMENT_RANGES_NM[segment]
    segment_wavelengths = np.arange(
        segment_min_nm,
        segment_max_nm + 1.0,
        1.0,
        dtype=np.float32,
    )
    segment_response = np.interp(
        segment_wavelengths,
        np.asarray(wavelengths_nm, dtype=np.float32),
        np.asarray(response, dtype=np.float32),
        left=0.0,
        right=0.0,
    ).astype(np.float32, copy=False)
    segment_response = np.clip(segment_response, 0.0, None)

    positive = np.flatnonzero(segment_response > 0.0)
    if positive.size == 0:
        raise ValueError(
            f"Band response does not overlap the supported {segment!r} segment range "
            f"{segment_min_nm:.0f}-{segment_max_nm:.0f} nm"
        )

    start = max(int(positive[0]) - 1, 0)
    stop = min(int(positive[-1]) + 2, segment_response.size)
    return segment_wavelengths[start:stop], segment_response[start:stop]


def _primary_nir_band_index(bands: Sequence[SensorBand]) -> int | None:
    candidates = [
        (idx, abs(float(band.center_wavelength) - 865.0))
        for idx, band in enumerate(bands)
        if 750.0 <= float(band.center_wavelength) <= 1000.0
    ]
    if not candidates:
        return None
    candidates.sort(key=lambda item: (item[1], item[0]))
    return int(candidates[0][0])


def _canonicalize_curve(
    wavelengths_nm: np.ndarray,
    response: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    wavelengths = np.asarray(wavelengths_nm, dtype=np.float32).reshape(-1)
    weights = np.asarray(response, dtype=np.float32).reshape(-1)
    if wavelengths.size != weights.size or wavelengths.size < 2:
        raise ValueError("Spectral response curves require at least two wavelength samples")
    order = np.argsort(wavelengths, kind="stable")
    wavelengths = wavelengths[order]
    weights = np.clip(weights[order], 0.0, None)
    unique_wavelengths, unique_idx = np.unique(wavelengths, return_index=True)
    wavelengths = unique_wavelengths.astype(np.float32, copy=False)
    weights = weights[unique_idx].astype(np.float32, copy=False)
    positive = np.flatnonzero(weights > 0.0)
    if positive.size == 0:
        raise ValueError("Spectral response curves must contain at least one positive response sample")
    start = max(int(positive[0]) - 1, 0)
    stop = min(int(positive[-1]) + 2, weights.size)
    return wavelengths[start:stop], weights[start:stop]


def _gaussian_curve_from_band(band: SensorBand) -> tuple[np.ndarray, np.ndarray]:
    sigma_nm = float(band.bandwidth) / (2.0 * np.sqrt(2.0 * np.log(2.0)))
    half_window = max(1.0, 4.0 * sigma_nm)
    start = max(350.0, float(band.center_wavelength) - half_window)
    stop = min(2500.0, float(band.center_wavelength) + half_window)
    wavelengths = np.arange(start, stop + 1.0, 1.0, dtype=np.float32)
    response = np.asarray(band.gaussian_response(wavelengths), dtype=np.float32)
    if response.size >= 2:
        response[0] = 0.0
        response[-1] = 0.0
    return _canonicalize_curve(wavelengths, response)


def _rsrf_curve_for_band(
    band: SensorBand,
    *,
    rsrf_root: Path | None,
) -> tuple[np.ndarray, np.ndarray] | None:
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

    return _canonicalize_curve(
        np.asarray(curve.wavelength_nm, dtype=np.float32),
        np.asarray(curve.response, dtype=np.float32),
    )


def _curve_for_band(
    band: SensorBand,
    *,
    rsrf_root: Path | None,
) -> tuple[np.ndarray, np.ndarray]:
    if band.has_rsrf:
        return _canonicalize_curve(
            np.asarray(band.rsrf_wavelengths_nm, dtype=np.float32),
            np.asarray(band.rsrf_response, dtype=np.float32),
        )
    rsrf_curve = _rsrf_curve_for_band(band, rsrf_root=rsrf_root)
    if rsrf_curve is not None:
        return rsrf_curve
    return _gaussian_curve_from_band(band)


def _sensor_id_for_bands(bands: Sequence[SensorBand], *, prefix: str) -> str:
    sensor_unit_ids = {band.rsrf_sensor_unit_id for band in bands if band.rsrf_sensor_unit_id}
    if len(sensor_unit_ids) == 1:
        return next(iter(sensor_unit_ids))
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
            row = {
                "source_id": library.source_id,
                "spectrum_id": sample_id,
                "sample_name": sample_id,
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
        rsrf_root=normalized_config.rsrf_root,
    )
    target_payload, target_schema_bands = _schema_payload_for_bands(
        target_sensor_id,
        target_bands,
        rsrf_root=normalized_config.rsrf_root,
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
        self._package_mapper: Any | None = None
        self._target_internal_to_output_index: dict[str, int] = {}
        self._source_retrieval_indices_by_segment: dict[str, np.ndarray] = {}
        self._target_schema_by_band_id: dict[str, Any] = {}

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
    ) -> tuple[xr.DataArray, xr.DataArray]:
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
            return (
                identity.transpose(*original_dims).astype(np.float32),
                unc.transpose(*original_dims).astype(np.float32),
            )

        if self._package_mapper is None or self._runtime is None:
            raise RuntimeError("spectral-library runtime was not initialized")

        transpose_dims = ["band", *[dim for dim in source_data.dims if dim != "band"]]
        source_values = np.asarray(source_data.transpose(*transpose_dims).values, dtype=np.float32)
        spatial_dims = transpose_dims[1:]
        flat = source_values.reshape(source_values.shape[0], -1).T
        valid_rows = np.any(np.isfinite(flat), axis=1)

        if source_unc is not None:
            unc_values = np.asarray(source_unc.transpose(*transpose_dims).values, dtype=np.float32)
            unc_flat = unc_values.reshape(unc_values.shape[0], -1).T
        else:
            unc_flat = None

        target_flat = np.full((flat.shape[0], len(self.target_bands)), np.nan, dtype=np.float32)
        target_unc_flat = np.full_like(target_flat, np.nan)
        if np.any(valid_rows):
            valid_queries = np.asarray(flat[valid_rows], dtype=np.float64)
            valid_masks = np.isfinite(valid_queries)
            batch = self._package_mapper.map_reflectance_batch(
                source_sensor=self._runtime.source_sensor_id,
                reflectance_rows=valid_queries,
                valid_mask_rows=valid_masks,
                output_mode="target_sensor",
                target_sensor=self._runtime.target_sensor_id,
                k=self.k_neighbors,
                min_valid_bands=self._mapping_config.min_valid_bands,
                neighbor_estimator=self._mapping_config.neighbor_estimator,
                knn_backend=self._mapping_config.knn_backend,
                knn_eps=self._mapping_config.knn_eps,
            )

            valid_indices = np.flatnonzero(valid_rows)
            for row_index, result in zip(valid_indices, batch.results, strict=True):
                if result.target_reflectance is not None:
                    for band_id, value in zip(result.target_band_ids, result.target_reflectance, strict=True):
                        target_index = self._target_internal_to_output_index[band_id]
                        target_flat[row_index, target_index] = np.float32(value)
                target_unc_flat[row_index] = self._estimate_uncertainty(
                    result,
                    source_uncertainty=None if unc_flat is None else unc_flat[row_index],
                )

        shape = (len(self.target_bands), *[source_data.sizes[dim] for dim in spatial_dims])
        reflectance = target_flat.T.reshape(shape)
        uncertainty = target_unc_flat.T.reshape(shape)
        coords = {"band": self._target_band_names}
        coords.update({dim: source_data.coords[dim] for dim in spatial_dims if dim in source_data.coords})
        reflectance_da = xr.DataArray(reflectance, dims=transpose_dims, coords=coords).astype(np.float32)
        uncertainty_da = xr.DataArray(uncertainty, dims=transpose_dims, coords=coords).astype(np.float32)
        return reflectance_da.transpose(*original_dims), uncertainty_da.transpose(*original_dims)

    def _align_source_data(self, data: xr.DataArray) -> xr.DataArray:
        if "band" not in data.dims:
            raise ValueError("source_reflectance must have a 'band' dimension")
        missing = [band.name for band in self.source_bands if band.name not in set(data.coords["band"].values.tolist())]
        if missing:
            raise KeyError(f"source_reflectance is missing source bands: {missing}")
        return data.sel(band=[band.name for band in self.source_bands])

    def _estimate_uncertainty(
        self,
        result: Any,
        *,
        source_uncertainty: np.ndarray | None,
    ) -> np.ndarray:
        if self._package_mapper is None or self._runtime is None:
            raise RuntimeError("spectral-library runtime was not initialized")

        output = np.full(len(self.target_bands), np.nan, dtype=np.float32)
        target_values_by_band_id = {}
        if result.target_reflectance is not None:
            target_values_by_band_id = {
                band_id: float(value)
                for band_id, value in zip(result.target_band_ids, result.target_reflectance, strict=True)
            }

        diagnostics = result.diagnostics.get("segments", {})
        for segment in ("vnir", "swir"):
            segment_diag = diagnostics.get(segment)
            if not isinstance(segment_diag, dict) or segment_diag.get("status") != "ok":
                continue

            neighbor_ids = tuple(str(value) for value in segment_diag.get("neighbor_ids", ()))
            if not neighbor_ids:
                continue
            neighbor_weights = np.asarray(segment_diag.get("neighbor_weights", ()), dtype=np.float64)
            if neighbor_weights.size == 0:
                continue
            weight_sum = float(np.sum(neighbor_weights))
            if weight_sum <= 0.0:
                continue
            neighbor_weights = neighbor_weights / weight_sum

            row_indices = np.asarray(
                [self._package_mapper._row_index_by_id[row_id] for row_id in neighbor_ids],  # noqa: SLF001
                dtype=np.int64,
            )
            neighbor_spectra = np.asarray(self._package_mapper._load_hyperspectral(segment)[row_indices], dtype=np.float64)  # noqa: SLF001

            if source_uncertainty is not None:
                retrieval_indices = self._source_retrieval_indices_by_segment[segment]
                query_valid_mask = np.asarray(segment_diag.get("query_valid_mask", ()), dtype=bool)
                valid_indices = retrieval_indices[query_valid_mask]
                if valid_indices.size > 0:
                    input_unc = float(np.sqrt(np.mean(np.square(source_uncertainty[valid_indices]), dtype=np.float64)))
                else:
                    input_unc = 0.0
            else:
                input_unc = 0.0
            source_fit_rmse = float(segment_diag.get("source_fit_rmse") or 0.0)

            for band_id, target_schema_band in self._target_schema_by_band_id.items():
                if getattr(target_schema_band, "segment", None) != segment:
                    continue
                target_index = self._target_internal_to_output_index.get(band_id)
                if target_index is None:
                    continue
                response = np.asarray(
                    self._package_mapper._band_response(  # noqa: SLF001
                        self._runtime.target_sensor_id,
                        target_schema_band,
                        segment_only=True,
                    ),
                    dtype=np.float64,
                )
                denominator = float(np.sum(response))
                if denominator <= 0.0:
                    continue
                neighbor_target = neighbor_spectra @ response / denominator
                estimate = target_values_by_band_id.get(band_id)
                if estimate is None:
                    estimate = float(np.dot(neighbor_weights, neighbor_target))
                spread = float(
                    np.sqrt(np.dot(neighbor_weights, np.square(neighbor_target - estimate, dtype=np.float64)))
                )
                output[target_index] = np.float32(
                    np.sqrt(max(spread, _UNCERTAINTY_FLOOR) ** 2 + source_fit_rmse**2 + input_unc**2)
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
    return mapper.map(source_reflectance, source_uncertainty=source_uncertainty)


__all__ = [
    "HyperspectralLibrary",
    "SpectralMappingConfig",
    "SpectralMapper",
    "convolve_hyperspectral_reflectance",
    "map_multispectral_reflectance",
    "needs_spectral_mapping",
]
