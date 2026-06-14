"""Spectral mapping between differing source and target sensor band sets."""

from __future__ import annotations

import hashlib
import json
import logging
import os
from dataclasses import dataclass
from pathlib import Path
from time import perf_counter
from typing import TYPE_CHECKING, Any, Protocol, TypeAlias, cast

import numpy as np
import xarray as xr
from numpy import typing as npt
from spectral_library import SpectralMapper as PackageSpectralMapper
from spectral_library.distribution import resolve_prepared_library_root
from spectral_library.mapping.engine.core import RSRF_SENSOR_BAND_SELECTIONS

from siac.algorithms.surface import _spectral_curve_utils as curve_utils
from siac.algorithms.surface._spectral_curve_utils import (
    normalized_band_response as _normalized_band_response,
)
from siac.geo._spatial import copy_spatial_metadata_like
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
    """Protocol matching ``spectral_library.BatchMappingArrayResult``."""

    reflectance: np.ndarray
    source_fit_rmse: np.ndarray
    output_columns: tuple[str, ...]


class _PackageMapperProtocol(Protocol):
    prepared_root: Path
    manifest: object

    def get_sensor_schema(self, sensor_id: str) -> object: ...

    def map_reflectance_batch_arrays_ndarray(
        self,
        *,
        source_sensor: str,
        reflectance_rows: Float64Array,
        valid_mask_rows: BoolArray,
        output_mode: str,
        target_sensor: str | None = None,
        k: int,
        min_valid_bands: int,
        neighbor_estimator: str,
        knn_backend: str,
        knn_eps: float,
    ) -> _BatchArrayResultProtocol: ...


_DEFAULT_K_NEIGHBORS = 5
_UNCERTAINTY_FLOOR = 0.005
_DEFAULT_NEIGHBOR_ESTIMATOR = "distance_weighted_mean"
_DEFAULT_KNN_BACKEND = "scipy_ckdtree"


@dataclass(frozen=True)
class SpectralMappingConfig:
    """Configuration for the published-runtime spectral-mapping adapter."""

    neighbor_estimator: str = _DEFAULT_NEIGHBOR_ESTIMATOR
    knn_backend: str = _DEFAULT_KNN_BACKEND
    knn_eps: float = 0.0
    min_valid_bands: int = 1

    def normalized(self) -> SpectralMappingConfig:
        return SpectralMappingConfig(
            neighbor_estimator=str(self.neighbor_estimator).strip() or _DEFAULT_NEIGHBOR_ESTIMATOR,
            knn_backend=str(self.knn_backend).strip() or _DEFAULT_KNN_BACKEND,
            knn_eps=float(self.knn_eps),
            min_valid_bands=max(1, int(self.min_valid_bands)),
        )


@dataclass(frozen=True)
class _PublishedRuntime:
    mapper: _PackageMapperProtocol
    prepared_root: Path
    source_sensor_id: str
    source_schema_band_ids: tuple[str, ...]
    source_schema_indices: tuple[int | None, ...]
    target_sensor_id: str
    target_output_band_ids: tuple[str, ...]
    ignored_source_band_names: tuple[str, ...]


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


def _split_mapping_inputs(
    spectral_library: object | None,
) -> SpectralMappingConfig:
    if spectral_library is None:
        return SpectralMappingConfig()
    if isinstance(spectral_library, SpectralMappingConfig):
        return spectral_library.normalized()
    raise TypeError(
        "spectral_library must be a SpectralMappingConfig or None. "
        "The published spectral-library runtime now owns the hyperspectral prior."
    )


def _hash_bytes(*chunks: bytes) -> str:
    digest = hashlib.sha256()
    for chunk in chunks:
        digest.update(chunk)
    return digest.hexdigest()


def _json_hash(payload: dict[str, object]) -> str:
    return _hash_bytes(json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8"))


def _native_rsrf_band_id(band: SensorBand) -> str:
    return str(band.rsrf_band_id or band.name)


def _source_sensor_id(source_bands: Sequence[SensorBand]) -> str:
    sensor_ids = {
        str(band.rsrf_sensor_unit_id)
        for band in source_bands
        if band.rsrf_sensor_unit_id is not None
    }
    if len(sensor_ids) != 1:
        raise ValueError(
            "Published spectral-library mapping requires all source bands to share one rsrf_sensor_unit_id."
        )
    return next(iter(sensor_ids))


def _schema_band_id_for_band(sensor_id: str, band: SensorBand) -> str | None:
    native_band_id = _native_rsrf_band_id(band)
    selections = RSRF_SENSOR_BAND_SELECTIONS.get(sensor_id)
    if selections is not None:
        for selection in selections:
            if selection.rsrf_band_id == native_band_id:
                return str(selection.band_id)
        return None
    return native_band_id


def _source_schema_indices_for_bands(
    source_bands: Sequence[SensorBand],
    *,
    sensor_id: str,
    source_schema_band_ids: tuple[str, ...],
) -> tuple[tuple[int | None, ...], tuple[str, ...]]:
    schema_index = {band_id: index for index, band_id in enumerate(source_schema_band_ids)}
    indices: list[int | None] = []
    ignored: list[str] = []
    seen_schema_ids: set[str] = set()

    for band in source_bands:
        canonical_band_id = _schema_band_id_for_band(sensor_id, band)
        if canonical_band_id is None or canonical_band_id not in schema_index:
            indices.append(None)
            ignored.append(str(band.name))
            continue
        if canonical_band_id in seen_schema_ids:
            raise ValueError(
                "Published spectral-library mapping resolved duplicate source bands onto the same canonical band.",
            )
        seen_schema_ids.add(canonical_band_id)
        indices.append(schema_index[canonical_band_id])

    return tuple(indices), tuple(ignored)


def _target_sensor_id(target_bands: Sequence[SensorBand]) -> str:
    sensor_ids = {
        str(band.rsrf_sensor_unit_id)
        for band in target_bands
        if band.rsrf_sensor_unit_id is not None
    }
    if len(sensor_ids) != 1:
        raise ValueError(
            "Published spectral-library target mapping requires all target bands to share one rsrf_sensor_unit_id."
        )
    return next(iter(sensor_ids))


def _target_output_band_ids_for_bands(
    target_bands: Sequence[SensorBand],
    *,
    sensor_id: str,
    target_schema_band_ids: tuple[str, ...],
) -> tuple[str, ...]:
    schema_index = {band_id: index for index, band_id in enumerate(target_schema_band_ids)}
    output_band_ids: list[str] = []
    seen_schema_ids: set[str] = set()

    for band in target_bands:
        schema_band_id = _schema_band_id_for_band(sensor_id, band)
        if schema_band_id is None or schema_band_id not in schema_index:
            raise ValueError(
                "Published spectral-library target mapping cannot resolve the requested target bands "
                "onto the upstream target sensor schema."
            )
        if schema_band_id in seen_schema_ids:
            raise ValueError(
                "Published spectral-library target mapping resolved duplicate target bands onto the same target band."
            )
        seen_schema_ids.add(schema_band_id)
        output_band_ids.append(schema_band_id)
    return tuple(output_band_ids)


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
        except (AttributeError, ValueError, RuntimeError):
            # rioxarray raises one of these when spatial metadata is
            # missing or coords aren't strictly monotonic. Other
            # exception classes propagate. (REVIEW.md §2.1)
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
    config: SpectralMappingConfig,
) -> _PublishedRuntime:
    del config
    source_sensor_id = _source_sensor_id(source_bands)
    prepared_root = resolve_prepared_library_root()
    mapper = cast(
        "_PackageMapperProtocol", PackageSpectralMapper(prepared_root, verify_checksums=False)
    )
    manifest_source_sensors = tuple(getattr(mapper.manifest, "source_sensors", ()))
    if source_sensor_id not in manifest_source_sensors:
        raise ValueError(
            "Published spectral-library runtime does not support the requested source sensor. "
            f"source_sensor={source_sensor_id!r} supported={manifest_source_sensors!r}"
        )
    schema: Any = mapper.get_sensor_schema(source_sensor_id)
    source_schema_band_ids = tuple(str(bid) for bid in schema.band_ids())
    source_schema_indices, ignored_source_band_names = _source_schema_indices_for_bands(
        source_bands,
        sensor_id=source_sensor_id,
        source_schema_band_ids=source_schema_band_ids,
    )
    if not any(index is not None for index in source_schema_indices):
        raise ValueError(
            "None of the requested source bands map onto the published spectral-library source basis."
        )
    target_sensor_id = _target_sensor_id(target_bands)
    target_schema: Any = mapper.get_sensor_schema(target_sensor_id)
    target_schema_band_ids = tuple(str(bid) for bid in target_schema.band_ids())
    target_output_band_ids = _target_output_band_ids_for_bands(
        target_bands,
        sensor_id=target_sensor_id,
        target_schema_band_ids=target_schema_band_ids,
    )
    return _PublishedRuntime(
        mapper=mapper,
        prepared_root=Path(prepared_root),
        source_sensor_id=source_sensor_id,
        source_schema_band_ids=source_schema_band_ids,
        source_schema_indices=source_schema_indices,
        target_sensor_id=target_sensor_id,
        target_output_band_ids=target_output_band_ids,
        ignored_source_band_names=ignored_source_band_names,
    )


class SpectralMapper:
    """Map multispectral reflectance using the published spectral-library runtime."""

    def __init__(
        self,
        source_bands: Sequence[SensorBand],
        target_bands: Sequence[SensorBand],
        *,
        spectral_library: SpectralMappingConfig | None = None,
        k_neighbors: int = _DEFAULT_K_NEIGHBORS,
    ) -> None:
        if k_neighbors < 1:
            raise ValueError("k_neighbors must be >= 1")

        self.source_bands = tuple(source_bands)
        self.target_bands = tuple(target_bands)
        self.k_neighbors = int(k_neighbors)
        self._identity = not needs_spectral_mapping(self.source_bands, self.target_bands)
        self._target_band_names = [band.name for band in self.target_bands]
        self._mapping_config = _split_mapping_inputs(spectral_library).normalized()
        self._runtime: _PublishedRuntime | None = None
        self._package_mapper: _PackageMapperProtocol | None = None
        self._supported_source_input_indices: tuple[int, ...] = ()

        if not self._identity:
            self._runtime = _prepare_runtime(
                self.source_bands,
                self.target_bands,
                config=self._mapping_config,
            )
            self._package_mapper = self._runtime.mapper
            self._supported_source_input_indices = tuple(
                index
                for index, schema_index in enumerate(self._runtime.source_schema_indices)
                if schema_index is not None
            )

    def _nearest_supported_source_per_target(self) -> npt.NDArray[np.intp]:
        """Column (into the supported-source uncertainty) of the spectrally
        nearest source band for each target band.

        Lets each target band inherit the uncertainty of the source band it
        chiefly derives from (e.g. S2 B02/B04 from the MODIS visible bands)
        instead of an average over all source bands, which dilutes the tight
        visible uncertainty with bright NIR/SWIR source bands.
        """
        supported_wl = np.array(
            [self.source_bands[i].center_wavelength for i in self._supported_source_input_indices],
            dtype=np.float64,
        )
        return np.array(
            [
                int(np.argmin(np.abs(supported_wl - band.center_wavelength)))
                for band in self.target_bands
            ],
            dtype=np.intp,
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
        source_unc = (
            self._align_source_data(source_uncertainty) if source_uncertainty is not None else None
        )

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
        package_queries, package_valid_masks, package_valid_rows = self._package_query_rows(
            flattened
        )
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

            output_index = {
                str(band_id): index for index, band_id in enumerate(batch_result.output_columns)
            }
            missing_output_band_ids = [
                band_id
                for band_id in self._runtime.target_output_band_ids
                if band_id not in output_index
            ]
            if missing_output_band_ids:
                raise RuntimeError(
                    "spectral-library target_sensor output did not include all requested target bands: "
                    f"{missing_output_band_ids!r}"
                )
            mapped_unique = np.asarray(
                batch_result.reflectance[
                    :,
                    [output_index[band_id] for band_id in self._runtime.target_output_band_ids],
                ],
                dtype=np.float32,
            )
            fit_rmse = np.asarray(batch_result.source_fit_rmse, dtype=np.float32)
            dedup_indices = deduplicated.inverse_indices
            valid_indices = np.flatnonzero(package_valid_rows)

            fill_started = perf_counter()
            target_flat[valid_indices] = mapped_unique[dedup_indices]
            source_fit_flat[valid_indices] = fit_rmse[dedup_indices]
            logger.info(
                "Spectral mapping reflectance fill complete: %d pixels elapsed=%.2fs",
                int(valid_indices.size),
                perf_counter() - fill_started,
            )

            unc_started = perf_counter()
            per_pixel_fit = fit_rmse[dedup_indices].astype(np.float64)
            n_target = target_flat.shape[1]
            # Per-target-band BRDF input uncertainty: each target band inherits the
            # uncertainty of the spectrally-nearest *supported* source band, rather
            # than the mean across all source bands. The mean let bright NIR/SWIR
            # source uncertainty dilute the tight visible uncertainty the
            # dark-target aerosol solve depends on — a dark visible target (B02)
            # was assigned ~0.013 instead of the source visible ~0.003, which kept
            # the surface prior too loose and let AOT collapse onto the CAMS prior.
            # target_flat has exactly one column per target band (n_target ==
            # len(self.target_bands)), so each target column maps to the nearest
            # source band.
            input_unc: npt.NDArray[np.float64] = np.zeros((valid_count, n_target), dtype=np.float64)
            if unc_flat is not None and self._supported_source_input_indices:
                supported_unc = np.nan_to_num(
                    np.asarray(
                        unc_flat[valid_indices][:, self._supported_source_input_indices],
                        dtype=np.float64,
                    ),
                    nan=0.0,
                    posinf=0.0,
                    neginf=0.0,
                )
                input_unc = supported_unc[:, self._nearest_supported_source_per_target()]
            unc_per_band = np.sqrt(
                _UNCERTAINTY_FLOOR**2 + per_pixel_fit[:, np.newaxis] ** 2 + input_unc**2
            ).astype(np.float32)
            mapped_mask = np.isfinite(target_flat[valid_indices])
            target_unc_flat[valid_indices] = np.where(mapped_mask, unc_per_band, np.nan)
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
        # Intentionally a no-op: the package mapper's distance metrics are
        # discarded under the upstream runtime path. ``test_spectral_mapping_extra
        # .test_map_does_not_cache_distance_metrics_under_upstream_runtime`` pins
        # this contract. Override in subclasses if persistence is needed.
        del metrics

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
        if self._runtime is None:
            raise RuntimeError("spectral-library runtime was not initialized")
        query_rows = np.full(
            (flattened.flat_values.shape[0], len(self._runtime.source_schema_band_ids)),
            np.nan,
            dtype=np.float64,
        )
        for source_index, schema_index in enumerate(self._runtime.source_schema_indices):
            if schema_index is None:
                continue
            query_rows[:, schema_index] = np.asarray(
                flattened.flat_values[:, source_index], dtype=np.float64
            )
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
    spectral_library: SpectralMappingConfig | None = None,
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
