"""Earthaccess-backed BRDF providers for MCD43, VNP43, and MCD19."""

from __future__ import annotations

import logging
from datetime import timedelta
from pathlib import Path
from typing import TYPE_CHECKING, Any, TypeAlias, cast

import numpy as np
import xarray as xr
from rasterio.enums import Resampling

from siac.adapters.data.earthaccess_source import EarthAccessSource
from siac.adapters.earthdata import (
    _build_virtual_stack_vrt,
    _target_bounds_from_template,
    build_earthaccess_runtime,
    build_source_aligned_target_template,
    constant_target_band_array,
    earthaccess_cache_dir,
    merge_reprojected_tiles,
    read_virtual_stack_to_target,
    select_candidate_paths,
    target_grid_coords,
)
from siac.adapters.earthdata_common import (
    ProductBandDefinition,
    apply_scale_and_mask,
    build_target_template,
    gdal_available,
    make_native_grid_dataarray,
    parse_granule_date,
    read_hdf4_dataset,
    read_hdf4_datasets,
    read_hdf4_datasets_attrs,
    read_hdf5_dataset,
    read_hdf5_datasets,
    read_hdf5_datasets_attrs,
    resolve_gdal_subdataset_path,
)
from siac.domain import SensorBand
from siac.runtime import BRDFKernelWeights

if TYPE_CHECKING:
    from collections.abc import Sequence
    from datetime import date, datetime

    from siac.adapters.data.earthaccess_catalog import EarthAccessCatalog

logger = logging.getLogger(__name__)

# HDF4/HDF5 libraries raise their own exception types on I/O failures.
# Import them defensively so they can be included in except clauses.
_HDF4Error: type[BaseException]
try:
    from pyhdf.error import HDF4Error as _HDF4Error  # type: ignore[import-untyped,no-redef]
except ImportError:  # pragma: no cover
    _HDF4Error = OSError
_HDF5Error: type[BaseException]
try:
    from h5py import HDF5ExtError as _HDF5Error  # type: ignore[attr-defined,no-redef]
except (ImportError, AttributeError):  # pragma: no cover
    _HDF5Error = OSError

# Combined tuple used in except clauses throughout this module.
_DATA_READ_ERRORS: tuple[type[BaseException], ...] = (
    OSError, KeyError, ValueError, TypeError, RuntimeError, _HDF4Error, _HDF5Error,
)

_BEST_QA_REFLECTANCE_UNCERTAINTY = 0.015
_QA_UNCERTAINTY_POWER = 1.6
_ORIGINAL_HDF4_READER = read_hdf4_dataset
_ORIGINAL_HDF5_READER = read_hdf5_dataset
_RequestedBand: TypeAlias = int | SensorBand
_RequestedBandCoord: TypeAlias = int | str
_RequestedBandSpec: TypeAlias = tuple[_RequestedBandCoord, ProductBandDefinition]
_PAYLOAD_FIELDS = ("f0", "f1", "f2", "unc")
_FLOAT32_NBYTES = np.dtype(np.float32).itemsize
_MAX_ONE_SHOT_TEMPORAL_VRT_OUTPUT_BYTES = 1024**3


def _granule_day(granule: object) -> np.datetime64 | None:
    """Best-effort extraction of the product day from an Earthaccess granule object."""
    render_dict = getattr(granule, "render_dict", None)
    if render_dict is not None:
        meta = render_dict.get("meta", {})
        umm = render_dict.get("umm", {})
    elif isinstance(granule, dict):
        meta = granule.get("meta", {})
        umm = granule.get("umm", {})
    else:
        meta = {}
        umm = {}

    filename_candidates: list[str] = []
    if isinstance(umm, dict):
        for key in ("GranuleUR", "ProducerGranuleId"):
            value = umm.get(key)
            if isinstance(value, str):
                filename_candidates.append(value)
    if isinstance(meta, dict):
        for key in ("native-id", "producer-granule-id"):
            value = meta.get(key)
            if isinstance(value, str):
                filename_candidates.append(value)

    for value in filename_candidates:
        try:
            return np.datetime64(parse_granule_date(Path(value)).date(), "D")
        except (ValueError, TypeError, IndexError):
            continue

    temporal_candidates: list[str] = []
    temporal = umm.get("TemporalExtent", {}) if isinstance(umm, dict) else {}
    if isinstance(temporal, dict):
        range_dt = temporal.get("RangeDateTime", {})
        if isinstance(range_dt, dict):
            for key in ("BeginningDateTime", "EndingDateTime"):
                value = range_dt.get(key)
                if isinstance(value, str):
                    temporal_candidates.append(value)
        single_dt = temporal.get("SingleDateTime")
        if isinstance(single_dt, str):
            temporal_candidates.append(single_dt)

    for value in temporal_candidates:
        try:
            return np.datetime64(value[:10], "D")
        except (ValueError, TypeError, IndexError):
            continue
    return None


def _granule_key(granule: object) -> str:
    """Best-effort stable key for deduplicating Earthaccess granules."""
    render_dict = getattr(granule, "render_dict", None)
    if render_dict is not None:
        meta = render_dict.get("meta", {})
        umm = render_dict.get("umm", {})
    elif isinstance(granule, dict):
        meta = granule.get("meta", {})
        umm = granule.get("umm", {})
    else:
        meta = {}
        umm = {}

    candidates: list[str] = []
    if isinstance(umm, dict):
        for key in ("GranuleUR", "ProducerGranuleId"):
            value = umm.get(key)
            if isinstance(value, str):
                candidates.append(value)
    if isinstance(meta, dict):
        for key in ("native-id", "producer-granule-id", "concept-id"):
            value = meta.get(key)
            if isinstance(value, str):
                candidates.append(value)

    for value in candidates:
        if value:
            return value
    return repr(granule)


class _EarthAccessBRDFProvider:
    """Shared logic for Earthaccess BRDF product wrappers."""

    product_key: str = ""
    _source_name: str = ""
    _product_bands: tuple[ProductBandDefinition, ...] = ()
    _legacy_band_map: dict[int, str] = {}
    _rsrf_sensor_unit_id: str | None = None
    _rsrf_representation_variant: str | None = "band_average"

    def __init__(
        self,
        cache_dir: str | Path | None = None,
        *,
        source: EarthAccessSource | None = None,
        catalog: EarthAccessCatalog | None = None,
        short_name: str | None = None,
        provider: str | None = None,
        probe_earthdata: bool = True,
        max_granules: int = 64,
    ) -> None:
        runtime = build_earthaccess_runtime(
            cache_dir=cache_dir,
            source=source,
            catalog=catalog,
            provider=provider,
        )
        self.cache_dir = runtime.cache_dir
        self.source = runtime.source
        self.catalog = runtime.catalog
        self.short_name = short_name
        self.provider = provider
        self.probe_earthdata = probe_earthdata
        resolved_max_granules = int(max_granules)
        self.max_granules = -1 if resolved_max_granules < 0 else max(1, resolved_max_granules)

    @property
    def source_name(self) -> str:
        return self._source_name

    @property
    def _bands_by_label(self) -> dict[str, ProductBandDefinition]:
        return {band.label: band for band in self._product_bands}

    @property
    def source_bands(self) -> tuple[SensorBand, ...]:
        return tuple(
            SensorBand(
                name=band.label,
                center_wavelength=band.wavelength_nm,
                bandwidth=band.bandwidth_nm,
                resolution=500.0,
                band_index=index,
                rsrf_sensor_unit_id=self._rsrf_sensor_unit_id,
                rsrf_representation_variant=self._rsrf_representation_variant,
                rsrf_band_id=band.rsrf_band_id or band.label,
            )
            for index, band in enumerate(self._product_bands)
        )

    def get_brdf_parameters(
        self,
        bounds: tuple[float, float, float, float],
        crs: str,
        obs_time: datetime,
        target_resolution: float,
        bands: Sequence[int] | Sequence[SensorBand],
        temporal_window: int = 16,
    ) -> BRDFKernelWeights:
        """Return BRDF kernel weights on the requested grid."""
        if target_resolution <= 0:
            raise ValueError(f"target_resolution must be > 0, got {target_resolution}")
        if not bands:
            raise ValueError("bands must be a non-empty sequence")

        requested = self._resolve_requested_bands(list(bands))
        if self.probe_earthdata:
            paths = self._download_granules(bounds, crs, obs_time, temporal_window)
            if paths:
                try:
                    return self._load_from_granules(
                        paths,
                        requested=requested,
                        bounds=bounds,
                        crs=crs,
                        target_resolution=target_resolution,
                    )
                except _DATA_READ_ERRORS as exc:  # pragma: no cover - external/system dependent
                    logger.warning(
                        "%s BRDF granule parsing failed; using defaults (%s)",
                        self._source_name,
                        exc,
                    )

        return self._default_weights(
            bounds,
            target_resolution,
            [coord for coord, _band in requested],
        )

    def get_temporal_brdf_parameters(
        self,
        bounds: tuple[float, float, float, float],
        crs: str,
        obs_time: datetime,
        target_resolution: float,
        bands: Sequence[int] | Sequence[SensorBand],
        temporal_window: int = 16,
        sample_dates: Sequence[date | datetime | np.datetime64] | None = None,
    ) -> BRDFKernelWeights:
        """Return temporal BRDF kernel weights with dims ``(time, band, y, x)``."""
        if target_resolution <= 0:
            raise ValueError(f"target_resolution must be > 0, got {target_resolution}")
        if not bands:
            raise ValueError("bands must be a non-empty sequence")

        requested = self._resolve_requested_bands(list(bands))
        time_axis = self._coerce_sample_time_axis(sample_dates) if sample_dates is not None else self._time_axis(obs_time, temporal_window)
        logger.info(
            "%s temporal BRDF request: obs_time=%s sample_days=%d resolved_bands=%d",
            self._source_name,
            obs_time.isoformat(),
            len(time_axis),
            len(requested),
        )

        if self.probe_earthdata:
            paths = self._download_granules(
                bounds,
                crs,
                obs_time,
                temporal_window,
                sample_dates=time_axis,
            )
            if paths:
                try:
                    logger.info(
                        "%s temporal BRDF request: selected %d candidate granule(s) for %s",
                        self._source_name,
                        len(paths),
                        obs_time.date().isoformat(),
                    )
                    return self._load_temporal_from_granules(
                        paths,
                        requested=requested,
                        bounds=bounds,
                        crs=crs,
                        target_resolution=target_resolution,
                        time_axis=time_axis,
                    )
                except _DATA_READ_ERRORS as exc:  # pragma: no cover - external/system dependent
                    logger.warning(
                        "%s temporal BRDF granule parsing failed; using defaults (%s)",
                        self._source_name,
                        exc,
                    )

        return self._default_temporal_weights(
            bounds,
            target_resolution,
            [coord for coord, _band in requested],
            time_axis,
            crs=crs,
        )

    def get_temporal_brdf_parameters_batch(
        self,
        *,
        bounds: tuple[float, float, float, float],
        crs: str,
        obs_times: Sequence[datetime],
        target_resolution: float,
        bands: Sequence[int] | Sequence[SensorBand],
        temporal_windows: Sequence[int],
        sample_date_sets: Sequence[Sequence[date | datetime | np.datetime64] | None],
    ) -> list[BRDFKernelWeights]:
        """Return multiple temporal BRDF stacks after one batched download."""
        if target_resolution <= 0:
            raise ValueError(f"target_resolution must be > 0, got {target_resolution}")
        if not bands:
            raise ValueError("bands must be a non-empty sequence")
        if not (len(obs_times) == len(temporal_windows) == len(sample_date_sets)):
            raise ValueError("obs_times, temporal_windows, and sample_date_sets must have the same length")

        requested = self._resolve_requested_bands(list(bands))
        requested_coords = [coord for coord, _band in requested]
        request_specs: list[tuple[datetime, np.ndarray]] = []
        for obs_time, temporal_window, sample_dates in zip(obs_times, temporal_windows, sample_date_sets, strict=True):
            time_axis = (
                self._coerce_sample_time_axis(sample_dates)
                if sample_dates is not None
                else self._time_axis(obs_time, temporal_window)
            )
            request_specs.append((obs_time, time_axis))

        logger.info(
            "%s batched temporal BRDF request: requests=%d resolved_bands=%d target_resolution=%s",
            self._source_name,
            len(request_specs),
            len(requested),
            target_resolution,
        )

        if self.probe_earthdata:
            downloaded = self._download_granules_batch(
                bounds=bounds,
                crs=crs,
                request_specs=request_specs,
                temporal_windows=list(temporal_windows),
            )
            if downloaded:
                outputs: list[BRDFKernelWeights] = []
                try:
                    selected_path_groups: list[list[Path]] = []
                    for request_index, (obs_time, time_axis) in enumerate(request_specs, start=1):
                        logger.info(
                            "%s temporal batch request %d/%d: obs_time=%s sample_days=%d first_day=%s last_day=%s",
                            self._source_name,
                            request_index,
                            len(request_specs),
                            obs_time.isoformat(),
                            len(time_axis),
                            str(time_axis[0]),
                            str(time_axis[-1]),
                        )
                        paths = self._select_candidate_paths(
                            downloaded,
                            obs_time,
                            bounds,
                            crs,
                            sample_dates=time_axis,
                        )
                        if paths:
                            logger.info(
                                "%s temporal batch request %d/%d: selected %d candidate granule(s)",
                                self._source_name,
                                request_index,
                                len(request_specs),
                                len(paths),
                            )
                        else:
                            logger.info(
                                "%s temporal batch request %d/%d: no candidate granules matched sampled days",
                                self._source_name,
                                request_index,
                                len(request_specs),
                            )
                        selected_path_groups.append(paths)

                    combined_time_axis = np.unique(
                        np.concatenate([time_axis for _obs_time, time_axis in request_specs])
                    ).astype("datetime64[D]")
                    combined_paths = list(
                        dict.fromkeys(path for paths in selected_path_groups for path in paths)
                    )
                    combined_temporal: BRDFKernelWeights | None = None
                    if combined_paths:
                        logger.info(
                            "%s temporal batch request: loading combined temporal stack for %d day(s) across %d request(s)",
                            self._source_name,
                            len(combined_time_axis),
                            len(request_specs),
                        )
                        combined_temporal = self._load_temporal_from_granules(
                            combined_paths,
                            requested=requested,
                            bounds=bounds,
                            crs=crs,
                            target_resolution=target_resolution,
                            time_axis=combined_time_axis,
                        )

                    for time_axis, paths in zip(
                        [time_axis for _obs_time, time_axis in request_specs],
                        selected_path_groups,
                        strict=True,
                    ):
                        if combined_temporal is not None and paths:
                            output = self._slice_temporal_weights(
                                combined_temporal,
                                source_time_axis=combined_time_axis,
                                request_time_axis=time_axis,
                            )
                            outputs.append(output)
                            continue
                        outputs.append(
                            self._nan_temporal_weights(
                                bounds,
                                target_resolution,
                                requested_coords,
                                time_axis,
                                crs=crs,
                            )
                        )
                    return outputs
                except _DATA_READ_ERRORS as exc:  # pragma: no cover - external/system dependent
                    logger.warning(
                        "%s batched temporal BRDF granule parsing failed; preserving NaN weights (%s)",
                        self._source_name,
                        exc,
                    )
                    return [
                        self._nan_temporal_weights(
                            bounds,
                            target_resolution,
                            requested_coords,
                            time_axis,
                            crs=crs,
                        )
                        for _obs_time, time_axis in request_specs
                    ]

        return [
            self._default_temporal_weights(
                bounds,
                target_resolution,
                [coord for coord, _band in requested],
                time_axis,
                crs=crs,
            )
            for _obs_time, time_axis in request_specs
        ]

    def _resolve_requested_bands(
        self,
        bands: Sequence[_RequestedBand],
    ) -> list[_RequestedBandSpec]:
        if not bands:
            raise ValueError("bands must be a non-empty sequence")

        product_bands = self._product_bands
        resolved: list[_RequestedBandSpec] = []
        for band in bands:
            if isinstance(band, SensorBand):
                match = min(
                    product_bands,
                    key=lambda candidate: abs(candidate.wavelength_nm - band.center_wavelength),
                )
                resolved.append((band.name, match))
                continue

            label = self._legacy_band_map.get(int(band))
            if label is None:
                raise KeyError(f"Band {band!r} is not available in {self._source_name}")
            resolved.append((int(band), self._bands_by_label[label]))
        return resolved

    @staticmethod
    def _slice_temporal_weights(
        weights: BRDFKernelWeights,
        *,
        source_time_axis: np.ndarray,
        request_time_axis: np.ndarray,
    ) -> BRDFKernelWeights:
        resolved_request_time_axis = np.asarray(request_time_axis, dtype="datetime64[D]")
        source_lookup = {
            np.datetime64(day, "D"): index
            for index, day in enumerate(np.asarray(source_time_axis, dtype="datetime64[D]"))
        }
        time_index = np.array(
            [source_lookup.get(np.datetime64(day, "D"), -1) for day in resolved_request_time_axis],
            dtype=np.intp,
        )
        coords = {"time": xr.IndexVariable("time", resolved_request_time_axis)}

        def _slice(data: xr.DataArray | None) -> xr.DataArray | None:
            if data is None:
                return None
            if time_index.size == 0:
                return data.isel(time=slice(0, 0)).assign_coords(coords)
            missing = time_index < 0
            filled_index = np.where(missing, 0, time_index)
            sliced = data.isel(time=filled_index).astype(np.float32, copy=False)
            if np.any(missing):
                sliced = sliced.copy(deep=True)
                sliced.values[missing] = np.nan
            return sliced.assign_coords(coords)

        return BRDFKernelWeights(
            f0=cast("xr.DataArray", _slice(weights.f0)),
            f1=cast("xr.DataArray", _slice(weights.f1)),
            f2=cast("xr.DataArray", _slice(weights.f2)),
            f0_unc=cast("xr.DataArray", _slice(weights.f0_unc)),
            f1_unc=cast("xr.DataArray", _slice(weights.f1_unc)),
            f2_unc=cast("xr.DataArray", _slice(weights.f2_unc)),
            reflectance_unc=_slice(weights.reflectance_unc),
        )

    def _download_granules(
        self,
        bounds: tuple[float, float, float, float],
        crs: str,
        obs_time: datetime,
        temporal_window: int,
        sample_dates: np.ndarray | None = None,
    ) -> list[Path]:
        short_name = self._resolved_short_name()
        temporal = self._temporal_search_window(obs_time, temporal_window, sample_dates)
        granules = self._search_granules(
            short_name=short_name,
            bounds=bounds,
            crs=crs,
            temporal=temporal,
        )
        if not granules:
            logger.warning(
                "%s granule probe returned no results for AOI/time window",
                self._source_name,
            )
            return []

        granules = self._filter_sampled_granules(
            granules,
            sample_dates,
            empty_message="%s granule probe returned no sampled results for AOI/time window",
        )
        if not granules:
            return []

        downloaded = self._download_granules_to_cache(granules, short_name=short_name)
        return self._select_candidate_paths(
            downloaded,
            obs_time,
            bounds,
            crs,
            sample_dates=sample_dates,
        )

    def _download_granules_batch(
        self,
        *,
        bounds: tuple[float, float, float, float],
        crs: str,
        request_specs: Sequence[tuple[datetime, np.ndarray]],
        temporal_windows: Sequence[int],
    ) -> list[Path]:
        short_name = self._resolved_short_name()
        unique_granules: dict[str, object] = {}
        for start_day, end_day, sample_dates in self._merge_search_batches(request_specs, temporal_windows):
            temporal = (
                f"{str(start_day)}T00:00:00Z",
                f"{str(end_day)}T23:59:59Z",
            )
            granules = self._filter_sampled_granules(
                self._search_granules(
                    short_name=short_name,
                    bounds=bounds,
                    crs=crs,
                    temporal=temporal,
                    count=-1,
                ),
                sample_dates,
            )
            for granule in granules:
                unique_granules[_granule_key(granule)] = granule

        if not unique_granules:
            logger.warning(
                "%s batched granule probe returned no sampled results for AOI/time windows",
                self._source_name,
            )
            return []

        return self._download_granules_to_cache(list(unique_granules.values()), short_name=short_name)

    def _resolved_short_name(self) -> str:
        return self.short_name or self.catalog.resolve_short_name(self.product_key)

    def _search_granules(
        self,
        *,
        short_name: str,
        bounds: tuple[float, float, float, float],
        crs: str,
        temporal: tuple[str, str],
        count: int | None = None,
    ) -> list[object]:
        from siac.adapters._retry import retry_transient

        return cast(
            "list[object]",
            retry_transient(
                self.source.search_granules,
                short_name=short_name,
                bounds=bounds,
                crs=crs,
                temporal=temporal,
                provider=self.provider,
                count=self.max_granules if count is None else count,
                max_attempts=3,
                base_delay_s=2.0,
                label=f"{self._source_name}.search_granules",
            ),
        )

    def _filter_sampled_granules(
        self,
        granules: list[object],
        sample_dates: np.ndarray | None,
        *,
        empty_message: str | None = None,
    ) -> list[object]:
        if sample_dates is None:
            return granules
        filtered = self._filter_granules_to_sample_dates(granules, sample_dates)
        if not filtered and empty_message is not None:
            logger.warning(empty_message, self._source_name)
        return filtered

    def _download_granules_to_cache(
        self,
        granules: list[object],
        *,
        short_name: str,
    ) -> list[Path]:
        from siac.adapters._retry import retry_transient

        dest = earthaccess_cache_dir(self.cache_dir, short_name)
        downloaded = cast(
            "list[Path]",
            retry_transient(
                self.source.download_granules,
                granules,
                dest,
                max_attempts=3,
                base_delay_s=3.0,
                label=f"{self._source_name}.download_granules",
            ),
        )
        # Warn if any downloaded path falls outside the expected cache directory.
        cache_root = Path(dest).resolve()
        for p in downloaded:
            resolved = Path(p).resolve()
            if not resolved.is_relative_to(cache_root):
                logger.warning(
                    "Downloaded file %s is outside cache directory %s.",
                    resolved, cache_root,
                )
        return downloaded

    @staticmethod
    def _merge_search_batches(
        request_specs: Sequence[tuple[datetime, np.ndarray]],
        temporal_windows: Sequence[int],
    ) -> list[tuple[Any, Any, np.ndarray]]:
        windows: list[tuple[Any, Any, np.ndarray]] = []
        for (obs_time, sample_dates), temporal_window in zip(request_specs, temporal_windows, strict=True):
            days = np.asarray(sample_dates, dtype="datetime64[D]")
            start_day: Any
            end_day: Any
            if days.size == 0:
                temporal = EarthAccessSource.temporal_window(obs_time, temporal_window)
                start_day = np.datetime64(temporal[0][:10], "D")
                end_day = np.datetime64(temporal[1][:10], "D")
            else:
                start_month = cast("np.datetime64", days.min().astype("datetime64[M]"))
                end_month = cast("np.datetime64", days.max().astype("datetime64[M]"))
                start_day = start_month.astype("datetime64[D]")
                end_month_day: np.datetime64 = (end_month + np.timedelta64(1, "M")).astype("datetime64[D]")
                end_day = end_month_day - np.timedelta64(1, "D")
                days = np.unique(days).astype("datetime64[D]")

            windows.append((start_day, end_day, days))

        windows.sort(key=lambda item: item[0])
        batches: list[tuple[Any, Any, np.ndarray]] = []
        for window in windows:
            start_day = window[0]
            end_day = window[1]
            sample_dates = window[2]
            if not batches:
                batches.append((start_day, end_day, sample_dates))
                continue

            prev_start, prev_end, prev_sample_dates = batches[-1]
            if start_day <= (prev_end + np.timedelta64(1, "D")):
                merged_start = min(prev_start, start_day)
                merged_end = max(prev_end, end_day)
                merged_sample_dates = np.unique(np.concatenate([prev_sample_dates, sample_dates])).astype("datetime64[D]")
                batches[-1] = (merged_start, merged_end, merged_sample_dates)
            else:
                batches.append((start_day, end_day, sample_dates))
        return batches

    @staticmethod
    def _select_candidate_paths(
        paths: list[Path],
        obs_time: datetime,
        bounds: tuple[float, float, float, float],
        crs: str,
        sample_dates: np.ndarray | None = None,
    ) -> list[Path]:
        return select_candidate_paths(
            paths,
            obs_time=obs_time,
            bounds=bounds,
            crs=crs,
            sample_dates=sample_dates,
        )

    def _load_from_granules(
        self,
        paths: list[Path],
        *,
        requested: list[tuple[int | str, ProductBandDefinition]],
        bounds: tuple[float, float, float, float],
        crs: str,
        target_resolution: float,
    ) -> BRDFKernelWeights:
        target_template, resolved_target_resolution = self._build_source_target_template(
            paths,
            requested=requested,
            bounds=bounds,
            crs=crs,
            target_resolution=target_resolution,
        )
        payload = self._merge_requested_payload(
            paths,
            requested=requested,
            bounds=bounds,
            crs=crs,
            target_resolution=resolved_target_resolution,
            target_template=target_template,
        )
        params, unc = self._unpack_payload_stack(payload, requested=requested)
        if not np.isfinite(params.values).any():
            logger.warning(
                "%s BRDF payload contained no finite values for bounds=%s crs=%s; preserving NaN weights",
                self._source_name,
                bounds,
                crs,
            )
            return self._weights_from_layers(
                params.astype(np.float32),
                unc.astype(np.float32),
            )
        return self._weights_from_layers(
            self._fill_parameter_defaults(params),
            unc.fillna(0.08),
        )

    @staticmethod
    def _merge_reprojected_tiles(
        arrays: list[xr.DataArray],
        *,
        bounds: tuple[float, float, float, float],
        crs: str,
        target_resolution: float,
        resampling: Resampling,
        nodata: float | None,
        target_template: xr.DataArray | None = None,
    ) -> xr.DataArray:
        return merge_reprojected_tiles(
            arrays,
            bounds=bounds,
            crs=crs,
            resolution=target_resolution,
            resampling=resampling,
            nodata=nodata,
            target_template=target_template,
        )

    @staticmethod
    def _build_source_target_template(
        paths: Sequence[str | Path],
        *,
        requested: Sequence[_RequestedBandSpec],
        bounds: tuple[float, float, float, float],
        crs: str,
        target_resolution: float,
    ) -> tuple[xr.DataArray, float]:
        try:
            first_dataset = resolve_gdal_subdataset_path(paths[0], requested[0][1].parameter_dataset)
            return build_source_aligned_target_template(
                first_dataset,
                bounds=bounds,
                crs=crs,
                resolution=target_resolution,
                resolution_name="target_resolution",
            )
        except _DATA_READ_ERRORS as exc:
            logger.debug(
                "%s could not build source-aligned BRDF target template from %s; using AOI grid (%s)",
                getattr(paths[0], "name", paths[0]),
                requested[0][1].parameter_dataset,
                exc,
            )
            return build_target_template(bounds, crs, target_resolution), float(target_resolution)

    @staticmethod
    def _virtual_stack_output_bytes(
        target_template: xr.DataArray,
        *,
        expected_layers: int,
    ) -> int:
        return (
            int(target_template.sizes["y"])
            * int(target_template.sizes["x"])
            * int(expected_layers)
            * _FLOAT32_NBYTES
        )

    def _allow_one_shot_temporal_vrt(
        self,
        *,
        target_template: xr.DataArray,
        expected_layers: int,
        label: str,
    ) -> bool:
        estimated_bytes = self._virtual_stack_output_bytes(
            target_template,
            expected_layers=expected_layers,
        )
        if estimated_bytes <= _MAX_ONE_SHOT_TEMPORAL_VRT_OUTPUT_BYTES:
            return True

        logger.info(
            "%s %s: estimated warped payload %.2f GiB across %d layer(s) exceeds the %.2f GiB one-shot budget; "
            "falling back to smaller merges",
            self._source_name,
            label,
            float(estimated_bytes) / (1024.0**3),
            int(expected_layers),
            float(_MAX_ONE_SHOT_TEMPORAL_VRT_OUTPUT_BYTES) / (1024.0**3),
        )
        return False

    @staticmethod
    def _coerce_sample_day(value: date | datetime | np.datetime64) -> np.datetime64:
        if isinstance(value, np.datetime64):
            return cast("np.datetime64", value.astype("datetime64[D]"))
        if hasattr(value, "date"):
            return np.datetime64(value.date(), "D")
        return np.datetime64(value, "D")

    @staticmethod
    def _time_axis(obs_time: datetime, temporal_window: int) -> np.ndarray:
        start = obs_time.date() - timedelta(days=int(temporal_window))
        axis: np.ndarray = np.array(
            [
                np.datetime64(start + timedelta(days=offset))
                for offset in range(2 * int(temporal_window) + 1)
            ],
            dtype="datetime64[D]",
        )
        return axis

    @staticmethod
    def _coerce_sample_time_axis(
        sample_dates: Sequence[date | datetime | np.datetime64],
    ) -> np.ndarray:
        if not sample_dates:
            raise ValueError("sample_dates must not be empty")
        axis: np.ndarray = np.unique(
            np.array(
                [_EarthAccessBRDFProvider._coerce_sample_day(value) for value in sample_dates],
                dtype="datetime64[D]",
            )
        )
        return axis

    @staticmethod
    def _temporal_search_window(
        obs_time: datetime,
        temporal_window: int,
        sample_dates: np.ndarray | None,
    ) -> tuple[str, str]:
        if sample_dates is None or sample_dates.size == 0:
            temporal: tuple[str, str] = EarthAccessSource.temporal_window(obs_time, temporal_window)
            return temporal
        min_day = cast("Any", sample_dates.min())
        max_day = cast("Any", sample_dates.max())
        start = f"{str(min_day)}T00:00:00Z"
        end = f"{str(max_day)}T23:59:59Z"
        return (start, end)

    @staticmethod
    def _filter_granules_to_sample_dates(
        granules: list[object],
        sample_dates: np.ndarray,
    ) -> list[object]:
        sample_day_set = {np.datetime64(day, "D") for day in sample_dates.tolist()}
        filtered: list[object] = []
        for granule in granules:
            granule_day = _granule_day(granule)
            if granule_day is None or granule_day in sample_day_set:
                filtered.append(granule)
        return filtered

    @staticmethod
    def _qa_to_uncertainty(qa: xr.DataArray) -> xr.DataArray:
        qa_values = _EarthAccessBRDFProvider._qa_values_to_uncertainty(qa.values)
        return xr.DataArray(qa_values, dims=qa.dims, coords=qa.coords)

    @staticmethod
    def _qa_values_to_uncertainty(qa_values: np.ndarray) -> np.ndarray:
        qa_values = np.asarray(qa_values, dtype=np.float32)
        unc = np.full(qa_values.shape, np.nan, dtype=np.float32)
        valid = np.isfinite(qa_values) & (qa_values >= 0.0)
        unc = np.where(
            valid,
            _BEST_QA_REFLECTANCE_UNCERTAINTY * np.power(qa_values + 1.0, _QA_UNCERTAINTY_POWER),
            unc,
        )
        unc_array: np.ndarray = np.asarray(unc, dtype=np.float32)
        return unc_array

    def _empty_spatial_array(
        self,
        bounds: tuple[float, float, float, float],
        resolution: float,
        *,
        fill_value: float = np.nan,
    ) -> xr.DataArray:
        y, x = self._grid(bounds, resolution)
        return xr.DataArray(
            np.full((y.size, x.size), fill_value, dtype=np.float32),
            dims=["y", "x"],
            coords={"y": y, "x": x},
        )

    def _coerce_to_target_grid(
        self,
        data: xr.DataArray,
        bounds: tuple[float, float, float, float],
        resolution: float,
    ) -> xr.DataArray:
        y, x = self._grid(bounds, resolution)
        extra_coords = [name for name in data.coords if name not in data.dims]
        if extra_coords:
            data = data.drop_vars(extra_coords, errors="ignore")
        return xr.DataArray(
            np.asarray(data.values, dtype=np.float32),
            dims=["y", "x"],
            coords={"y": y, "x": x},
        )

    def _default_temporal_weights(
        self,
        bounds: tuple[float, float, float, float],
        resolution: float,
        bands: list[int | str],
        time_axis: np.ndarray,
        *,
        crs: str | None = None,
        target_template: xr.DataArray | None = None,
    ) -> BRDFKernelWeights:
        if target_template is None:
            if crs is None:
                raise ValueError("crs is required when target_template is not provided")
            target_template = build_target_template(bounds, crs, resolution)
        time_coords = xr.IndexVariable("time", time_axis)
        y_coords = target_template.coords["y"]
        x_coords = target_template.coords["x"]
        shape = (
            len(time_axis),
            len(bands),
            int(target_template.sizes["y"]),
            int(target_template.sizes["x"]),
        )

        def _constant(fill_value: float) -> xr.DataArray:
            return self._target_array_like(
                target_template,
                np.full(shape, np.float32(fill_value), dtype=np.float32),
                dims=("time", "band", "y", "x"),
                coords={
                    "time": time_coords,
                    "band": xr.IndexVariable("band", bands),
                    "y": y_coords,
                    "x": x_coords,
                },
            )

        return BRDFKernelWeights(
            f0=_constant(0.20),
            f1=_constant(0.05),
            f2=_constant(0.02),
            f0_unc=_constant(0.03),
            f1_unc=_constant(0.02),
            f2_unc=_constant(0.02),
            reflectance_unc=_constant(0.08),
        )

    def _nan_temporal_weights(
        self,
        bounds: tuple[float, float, float, float],
        resolution: float,
        bands: list[int | str],
        time_axis: np.ndarray,
        *,
        crs: str | None = None,
        target_template: xr.DataArray | None = None,
    ) -> BRDFKernelWeights:
        base = self._default_temporal_weights(
            bounds,
            resolution,
            bands,
            time_axis,
            crs=crs,
            target_template=target_template,
        )

        def _nan_like(data: xr.DataArray | None) -> xr.DataArray | None:
            if data is None:
                return None
            return xr.full_like(data, np.nan, dtype=np.float32)

        return BRDFKernelWeights(
            f0=cast("xr.DataArray", _nan_like(base.f0)),
            f1=cast("xr.DataArray", _nan_like(base.f1)),
            f2=cast("xr.DataArray", _nan_like(base.f2)),
            f0_unc=cast("xr.DataArray", _nan_like(base.f0_unc)),
            f1_unc=cast("xr.DataArray", _nan_like(base.f1_unc)),
            f2_unc=cast("xr.DataArray", _nan_like(base.f2_unc)),
            reflectance_unc=_nan_like(base.reflectance_unc),
        )

    @staticmethod
    def _allocate_temporal_payload_arrays(
        time_axis: np.ndarray,
        requested: Sequence[_RequestedBandSpec],
        *,
        y_size: int,
        x_size: int,
    ) -> tuple[np.ndarray, np.ndarray]:
        params_values: np.ndarray = np.full(
            (len(time_axis), len(requested), 3, y_size, x_size),
            np.nan,
            dtype=np.float32,
        )
        unc_values: np.ndarray = np.full(
            (len(time_axis), len(requested), y_size, x_size),
            np.nan,
            dtype=np.float32,
        )
        return params_values, unc_values

    @staticmethod
    def _allocate_spatial_payload_arrays(
        requested: Sequence[_RequestedBandSpec],
        target_template: xr.DataArray,
    ) -> tuple[np.ndarray, np.ndarray]:
        y_size = int(target_template.sizes["y"])
        x_size = int(target_template.sizes["x"])
        params_values: np.ndarray = np.full(
            (len(requested), 3, y_size, x_size),
            np.nan,
            dtype=np.float32,
        )
        unc_values: np.ndarray = np.full(
            (len(requested), y_size, x_size),
            np.nan,
            dtype=np.float32,
        )
        return params_values, unc_values

    @classmethod
    def _temporal_weights_from_arrays(
        cls,
        params_values: np.ndarray,
        unc_values: np.ndarray,
        *,
        requested: Sequence[_RequestedBandSpec],
        time_axis: np.ndarray,
        target_template: xr.DataArray,
    ) -> BRDFKernelWeights:
        coords = {
            "time": xr.IndexVariable("time", time_axis),
            "band": xr.IndexVariable("band", [band_coord for band_coord, _band in requested]),
            "y": target_template.coords["y"],
            "x": target_template.coords["x"],
        }

        def _wrap(values: np.ndarray) -> xr.DataArray:
            return cls._target_array_like(
                target_template,
                values,
                dims=("time", "band", "y", "x"),
                coords=coords,
            )

        scaled_unc = unc_values * np.float32(1.1)
        return BRDFKernelWeights(
            f0=_wrap(params_values[:, :, 0, :, :]),
            f1=_wrap(params_values[:, :, 1, :, :]),
            f2=_wrap(params_values[:, :, 2, :, :]),
            f0_unc=_wrap(unc_values),
            f1_unc=_wrap(scaled_unc),
            f2_unc=_wrap(scaled_unc),
            reflectance_unc=_wrap(unc_values),
        )

    def _load_temporal_payload_vrt(
        self,
        grouped_paths: dict[np.datetime64, list[Path]],
        *,
        requested: Sequence[_RequestedBandSpec],
        bounds: tuple[float, float, float, float],
        crs: str,
        target_resolution: float,
        time_axis: np.ndarray,
        target_template: xr.DataArray,
    ) -> tuple[np.ndarray, np.ndarray] | None:
        del grouped_paths, requested, bounds, crs, target_resolution, time_axis, target_template
        return None

    def _load_temporal_payload_vrt_from_daily_payload_vrts(
        self,
        grouped_paths: dict[np.datetime64, list[Path]],
        *,
        requested: Sequence[_RequestedBandSpec],
        bounds: tuple[float, float, float, float],
        crs: str,
        target_resolution: float,
        time_axis: np.ndarray,
        target_template: xr.DataArray,
    ) -> tuple[np.ndarray, np.ndarray] | None:
        del grouped_paths, requested, bounds, crs, target_resolution, time_axis, target_template
        return None

    def _load_temporal_from_granules(
        self,
        paths: list[Path],
        *,
        requested: list[tuple[int | str, ProductBandDefinition]],
        bounds: tuple[float, float, float, float],
        crs: str,
        target_resolution: float,
        time_axis: np.ndarray,
    ) -> BRDFKernelWeights:
        grouped_paths: dict[np.datetime64, list[Path]] = {}
        for path in paths:
            granule_day = np.datetime64(parse_granule_date(path).date())
            grouped_paths.setdefault(granule_day, []).append(path)
        logger.info(
            "%s temporal BRDF stack: %d/%d day%s available in search window.",
            self._source_name,
            len(grouped_paths),
            len(time_axis),
            "" if len(time_axis) == 1 else "s",
        )
        logger.debug(
            "%s temporal BRDF stack days: requested=%s available=%s",
            self._source_name,
            [str(day) for day in time_axis.tolist()],
            sorted(str(day) for day in grouped_paths),
        )

        target_template, resolved_target_resolution = self._build_source_target_template(
            paths,
            requested=requested,
            bounds=bounds,
            crs=crs,
            target_resolution=target_resolution,
        )
        logger.info(
            "%s temporal BRDF stack: attempting direct temporal VRT read for %d day(s) x %d band(s)",
            self._source_name,
            len(grouped_paths),
            len(requested),
        )
        temporal_payload = self._load_temporal_payload_vrt(
            grouped_paths,
            requested=requested,
            bounds=bounds,
            crs=crs,
            target_resolution=resolved_target_resolution,
            time_axis=time_axis,
            target_template=target_template,
        )
        if temporal_payload is None:
            logger.info(
                "%s temporal BRDF stack: direct temporal VRT path unavailable, attempting one-warp daily-payload fallback",
                self._source_name,
            )
            temporal_payload = self._load_temporal_payload_vrt_from_daily_payload_vrts(
                grouped_paths,
                requested=requested,
                bounds=bounds,
                crs=crs,
                target_resolution=resolved_target_resolution,
                time_axis=time_axis,
                target_template=target_template,
            )
            if temporal_payload is None:
                logger.info(
                    "%s temporal BRDF stack: daily-payload fallback unavailable, falling back to per-day merges",
                    self._source_name,
                )
                params_values, unc_values = self._allocate_temporal_payload_arrays(
                    time_axis,
                    requested,
                    y_size=int(target_template.sizes["y"]),
                    x_size=int(target_template.sizes["x"]),
                )
            else:
                logger.info("%s temporal BRDF stack: daily-payload fallback completed", self._source_name)
                params_values, unc_values = temporal_payload
        else:
            logger.info("%s temporal BRDF stack: direct temporal VRT read completed", self._source_name)
            params_values, unc_values = temporal_payload

        def _coerce_daily(data: xr.DataArray) -> np.ndarray:
            extra_coords = [name for name in data.coords if name not in data.dims]
            if extra_coords:
                data = data.drop_vars(extra_coords, errors="ignore")
            daily_values: np.ndarray = np.asarray(data.transpose("band", "y", "x").values, dtype=np.float32)
            return daily_values

        if temporal_payload is None:
            available_day_total = int(sum(1 for day in time_axis if grouped_paths.get(day)))
            completed_days = 0
            for time_index, day in enumerate(time_axis):
                day_paths = grouped_paths.get(day, [])
                if not day_paths:
                    continue
                completed_days += 1
                logger.info(
                    "%s temporal BRDF stack: per-day fallback %d/%d day=%s tiles=%d",
                    self._source_name,
                    completed_days,
                    available_day_total,
                    str(day),
                    len(day_paths),
                )

                payload = self._merge_requested_payload(
                    day_paths,
                    requested=requested,
                    bounds=bounds,
                    crs=crs,
                    target_resolution=resolved_target_resolution,
                    target_template=target_template,
                )
                params, unc = self._unpack_payload_stack(payload, requested=requested)
                params_values[time_index] = np.asarray(
                    params.transpose("band", "parameter", "y", "x").values,
                    dtype=np.float32,
                )
                unc_values[time_index] = _coerce_daily(unc)

        temporal = self._temporal_weights_from_arrays(
            params_values,
            unc_values,
            requested=requested,
            time_axis=time_axis,
            target_template=target_template,
        )
        if np.isfinite(temporal.f0.values).any():
            return temporal
        logger.warning(
            "%s temporal BRDF payload contained no finite values for bounds=%s crs=%s days=%d; preserving NaN weights",
            self._source_name,
            bounds,
            crs,
            len(time_axis),
        )
        return temporal

    def _load_native_band_stack(
        self,
        path: str | Path,
        product_band: ProductBandDefinition,
    ) -> tuple[tuple[xr.DataArray, xr.DataArray, xr.DataArray], xr.DataArray]:
        raise NotImplementedError

    def _load_native_requested_payload(
        self,
        path: str | Path,
        requested: Sequence[_RequestedBandSpec],
    ) -> xr.DataArray:
        params, qa = self._load_native_requested_stack(path, requested)
        # Carry uncertainty alongside the BRDF parameters so each granule only
        # goes through one spatial merge/reprojection path.
        return self._pack_payload_stack(params, self._qa_to_uncertainty(qa))

    def _load_native_requested_stack(
        self,
        path: str | Path,
        requested: Sequence[_RequestedBandSpec],
    ) -> tuple[xr.DataArray, xr.DataArray]:
        band_coords = [band_coord for band_coord, _product_band in requested]
        param_tiles = []
        qa_tiles = []
        for _band_coord, product_band in requested:
            params, qa = self._load_native_band_stack(path, product_band)
            param_tiles.append(self._stack_parameter_cube(params))
            qa_tiles.append(qa)
        return (
            xr.concat(param_tiles, dim=xr.IndexVariable("band", band_coords)).transpose("band", "parameter", "y", "x"),
            xr.concat(qa_tiles, dim=xr.IndexVariable("band", band_coords)).transpose("band", "y", "x"),
        )

    @staticmethod
    def _grid(bounds: tuple[float, float, float, float], resolution: float) -> tuple[np.ndarray, np.ndarray]:
        return target_grid_coords(bounds, resolution, resolution_name="target_resolution")

    def _constant_band_array(
        self,
        bands: list[int | str],
        bounds: tuple[float, float, float, float],
        resolution: float,
        value: float,
    ) -> xr.DataArray:
        return constant_target_band_array(
            bands,
            bounds,
            resolution,
            value,
            resolution_name="target_resolution",
        )

    def _merge_requested_payload(
        self,
        paths: Sequence[str | Path],
        *,
        requested: Sequence[_RequestedBandSpec],
        bounds: tuple[float, float, float, float],
        crs: str,
        target_resolution: float,
        target_template: xr.DataArray,
    ) -> xr.DataArray:
        payload = self._load_requested_payload_vrt(
            paths,
            requested=requested,
            bounds=bounds,
            crs=crs,
            target_resolution=target_resolution,
            target_template=target_template,
        )
        if payload is not None:
            return payload

        payload_tiles = [self._load_native_requested_payload(path, requested) for path in paths]
        return self._merge_reprojected_tiles(
            payload_tiles,
            bounds=bounds,
            crs=crs,
            target_resolution=target_resolution,
            resampling=Resampling.nearest,
            nodata=np.nan,
            target_template=target_template,
        )

    def _load_requested_payload_vrt(
        self,
        paths: Sequence[str | Path],
        *,
        requested: Sequence[_RequestedBandSpec],
        bounds: tuple[float, float, float, float],
        crs: str,
        target_resolution: float,
        target_template: xr.DataArray,
    ) -> xr.DataArray | None:
        del paths, requested, bounds, crs, target_resolution, target_template
        return None

    @staticmethod
    def _native_array_like(
        reference: xr.DataArray,
        values: np.ndarray,
        *,
        dims: tuple[str, ...],
        coords: dict[str, object],
    ) -> xr.DataArray:
        out = xr.DataArray(np.asarray(values, dtype=np.float32), dims=dims, coords=coords)
        out = out.rio.set_spatial_dims(x_dim="x", y_dim="y")
        reference_crs = reference.rio.crs
        if reference_crs is not None:
            out = out.rio.write_crs(reference_crs)
        try:
            return out.rio.write_transform(reference.rio.transform(recalc=True))
        except _DATA_READ_ERRORS:
            return out

    @classmethod
    def _target_array_like(
        cls,
        target_template: xr.DataArray,
        values: np.ndarray,
        *,
        dims: tuple[str, ...],
        coords: dict[str, object],
    ) -> xr.DataArray:
        out = xr.DataArray(np.asarray(values, dtype=np.float32), dims=dims, coords=coords)
        out = out.rio.set_spatial_dims(x_dim="x", y_dim="y")
        out = out.rio.write_crs(target_template.rio.crs)
        return out.rio.write_transform(target_template.rio.transform(recalc=True))

    @classmethod
    def _pack_payload_stack(
        cls,
        params: xr.DataArray,
        unc: xr.DataArray,
    ) -> xr.DataArray:
        params = params.transpose("band", "parameter", "y", "x")
        unc = unc.transpose("band", "y", "x")
        band_coords = np.asarray(params.coords["band"].values, dtype=object)
        payload_values = np.concatenate(
            [
                np.asarray(params.values, dtype=np.float32),
                np.asarray(unc.values, dtype=np.float32)[:, np.newaxis, :, :],
            ],
            axis=1,
        )
        layer_count = band_coords.size * len(_PAYLOAD_FIELDS)
        layer_values = payload_values.reshape(layer_count, *payload_values.shape[-2:])
        return cls._native_array_like(
            params,
            layer_values,
            dims=("layer", "y", "x"),
            coords={
                "layer": xr.IndexVariable("layer", np.arange(layer_count, dtype=np.int32)),
                "y": params.coords["y"],
                "x": params.coords["x"],
            },
        )

    @classmethod
    def _unpack_payload_stack(
        cls,
        payload: xr.DataArray,
        *,
        requested: Sequence[_RequestedBandSpec],
    ) -> tuple[xr.DataArray, xr.DataArray]:
        payload = payload.transpose("layer", "y", "x")
        band_coords = [band_coord for band_coord, _band in requested]
        expected_layers = len(band_coords) * len(_PAYLOAD_FIELDS)
        if payload.sizes["layer"] != expected_layers:
            raise ValueError(
                f"Expected {expected_layers} payload layers for {len(band_coords)} band(s), "
                f"got {payload.sizes['layer']}"
            )
        values = np.asarray(payload.values, dtype=np.float32).reshape(
            len(band_coords),
            len(_PAYLOAD_FIELDS),
            payload.sizes["y"],
            payload.sizes["x"],
        )
        coords = {
            "band": xr.IndexVariable("band", band_coords),
            "y": payload.coords["y"],
            "x": payload.coords["x"],
        }
        params = cls._native_array_like(
            payload,
            values[:, :3, :, :],
            dims=("band", "parameter", "y", "x"),
            coords={
                **coords,
                "parameter": xr.IndexVariable("parameter", ["f0", "f1", "f2"]),
            },
        )
        unc = cls._native_array_like(
            payload,
            values[:, 3, :, :],
            dims=("band", "y", "x"),
            coords=coords,
        )
        return params, unc

    @staticmethod
    def _stack_parameter_cube(
        params: tuple[xr.DataArray, xr.DataArray, xr.DataArray],
    ) -> xr.DataArray:
        return xr.concat(
            list(params),
            dim=xr.IndexVariable("parameter", ["f0", "f1", "f2"]),
        )

    @staticmethod
    def _fill_parameter_defaults(params: xr.DataArray) -> xr.DataArray:
        defaults = xr.DataArray(
            np.array([0.20, 0.05, 0.02], dtype=np.float32),
            dims=["parameter"],
            coords={"parameter": ["f0", "f1", "f2"]},
        )
        return params.fillna(defaults)

    @staticmethod
    def _weights_from_layers(
        params: xr.DataArray,
        unc: xr.DataArray,
    ) -> BRDFKernelWeights:
        unc = unc.transpose("band", "y", "x")
        return BRDFKernelWeights(
            f0=params.sel(parameter="f0", drop=True).transpose("band", "y", "x"),
            f1=params.sel(parameter="f1", drop=True).transpose("band", "y", "x"),
            f2=params.sel(parameter="f2", drop=True).transpose("band", "y", "x"),
            f0_unc=unc,
            f1_unc=(unc * np.float32(1.1)).transpose("band", "y", "x"),
            f2_unc=(unc * np.float32(1.1)).transpose("band", "y", "x"),
            reflectance_unc=unc,
        )

    def _default_weights(
        self,
        bounds: tuple[float, float, float, float],
        resolution: float,
        bands: list[int | str],
    ) -> BRDFKernelWeights:
        if not bands:
            raise ValueError("bands must be a non-empty sequence")

        f0 = self._constant_band_array(bands, bounds, resolution, 0.20)
        f1 = self._constant_band_array(bands, bounds, resolution, 0.05)
        f2 = self._constant_band_array(bands, bounds, resolution, 0.02)

        return BRDFKernelWeights(
            f0=f0,
            f1=f1,
            f2=f2,
            f0_unc=xr.full_like(f0, 0.03),
            f1_unc=xr.full_like(f1, 0.02),
            f2_unc=xr.full_like(f2, 0.02),
            reflectance_unc=xr.full_like(f0, 0.08),
        )


class _StackParameterProvider(_EarthAccessBRDFProvider):
    """BRDF provider where one dataset stores all three kernel parameters."""

    _read_dataset = staticmethod(read_hdf4_dataset)

    def _read_named_datasets(
        self,
        path: str | Path,
        dataset_names: tuple[str, ...],
    ) -> dict[str, tuple[np.ndarray, dict[str, Any]]]:
        if self._read_dataset is _ORIGINAL_HDF4_READER:
            return cast("dict[str, tuple[np.ndarray, dict[str, Any]]]", read_hdf4_datasets(path, dataset_names))
        if self._read_dataset is _ORIGINAL_HDF5_READER:
            return cast("dict[str, tuple[np.ndarray, dict[str, Any]]]", read_hdf5_datasets(path, dataset_names))
        return {
            dataset_name: self._read_dataset(path, dataset_name)
            for dataset_name in dataset_names
        }

    def _read_named_dataset_attrs(
        self,
        path: str | Path,
        dataset_names: tuple[str, ...],
    ) -> dict[str, dict[str, Any]]:
        if self._read_dataset is _ORIGINAL_HDF4_READER:
            return cast("dict[str, dict[str, Any]]", read_hdf4_datasets_attrs(path, dataset_names))
        if self._read_dataset is _ORIGINAL_HDF5_READER:
            return cast("dict[str, dict[str, Any]]", read_hdf5_datasets_attrs(path, dataset_names))
        return {
            dataset_name: self._read_dataset(path, dataset_name)[1]
            for dataset_name in dataset_names
        }

    def _load_native_band_stack(
        self,
        path: str | Path,
        product_band: ProductBandDefinition,
    ) -> tuple[tuple[xr.DataArray, xr.DataArray, xr.DataArray], xr.DataArray]:
        dataset_names = (product_band.parameter_dataset, product_band.qa_dataset or "")
        datasets = self._read_named_datasets(path, dataset_names)
        params_raw, params_attrs = datasets[dataset_names[0]]
        qa_raw, qa_attrs = datasets[dataset_names[1]]

        params = apply_scale_and_mask(params_raw, params_attrs)
        if params.ndim != 3 or params.shape[-1] != 3:
            raise ValueError(
                f"Expected a 3-parameter BRDF stack for {product_band.parameter_dataset}, "
                f"got shape={params.shape}"
            )

        qa = apply_scale_and_mask(qa_raw, qa_attrs)
        param_da = make_native_grid_dataarray(
            np.moveaxis(params, -1, 0),
            granule_path=path,
            dims=("parameter", "y", "x"),
            coords={"parameter": ["f0", "f1", "f2"]},
        )
        qa_da = make_native_grid_dataarray(qa, granule_path=path)
        return (
            param_da.sel(parameter="f0", drop=True),
            param_da.sel(parameter="f1", drop=True),
            param_da.sel(parameter="f2", drop=True),
        ), qa_da

    def _load_native_requested_stack(
        self,
        path: str | Path,
        requested: Sequence[_RequestedBandSpec],
    ) -> tuple[xr.DataArray, xr.DataArray]:
        band_coords = [band_coord for band_coord, _product_band in requested]
        dataset_names = tuple(
            dict.fromkeys(
                [
                    dataset_name
                    for _band_coord, product_band in requested
                    for dataset_name in (product_band.parameter_dataset, product_band.qa_dataset or "")
                ]
            )
        )
        datasets = self._read_named_datasets(path, dataset_names)

        param_tiles = []
        qa_tiles = []
        for _band_coord, product_band in requested:
            params_raw, params_attrs = datasets[product_band.parameter_dataset]
            qa_raw, qa_attrs = datasets[product_band.qa_dataset or ""]
            params = apply_scale_and_mask(params_raw, params_attrs)
            if params.ndim != 3 or params.shape[-1] != 3:
                raise ValueError(
                    f"Expected a 3-parameter BRDF stack for {product_band.parameter_dataset}, "
                    f"got shape={params.shape}"
                )
            qa = apply_scale_and_mask(qa_raw, qa_attrs)
            param_tiles.append(
                make_native_grid_dataarray(
                    np.moveaxis(params, -1, 0),
                    granule_path=path,
                    dims=("parameter", "y", "x"),
                    coords={"parameter": ["f0", "f1", "f2"]},
                )
            )
            qa_tiles.append(make_native_grid_dataarray(qa, granule_path=path))

        return (
            xr.concat(param_tiles, dim=xr.IndexVariable("band", band_coords)).transpose("band", "parameter", "y", "x"),
            xr.concat(qa_tiles, dim=xr.IndexVariable("band", band_coords)).transpose("band", "y", "x"),
        )

    def _load_requested_payload_vrt(
        self,
        paths: Sequence[str | Path],
        *,
        requested: Sequence[_RequestedBandSpec],
        bounds: tuple[float, float, float, float],
        crs: str,
        target_resolution: float,
        target_template: xr.DataArray,
    ) -> xr.DataArray | None:
        if not gdal_available():
            return None

        dataset_names = tuple(
            dict.fromkeys(
                [
                    dataset_name
                    for _band_coord, product_band in requested
                    for dataset_name in (product_band.parameter_dataset, product_band.qa_dataset or "")
                ]
            )
        )
        try:
            dataset_attrs = self._read_named_dataset_attrs(paths[0], dataset_names)
            source_groups: list[list[str]] = []
            group_band_counts: list[int] = []
            param_attrs: list[dict[str, Any]] = []
            qa_attrs: list[dict[str, Any]] = []
            band_coords = [band_coord for band_coord, _product_band in requested]

            for _band_coord, product_band in requested:
                qa_name = product_band.qa_dataset or ""
                source_groups.append(
                    [resolve_gdal_subdataset_path(path, product_band.parameter_dataset) for path in paths]
                )
                group_band_counts.append(3)
                param_attrs.append(dataset_attrs[product_band.parameter_dataset])

                source_groups.append([resolve_gdal_subdataset_path(path, qa_name) for path in paths])
                group_band_counts.append(1)
                qa_attrs.append(dataset_attrs[qa_name])

            stacked = read_virtual_stack_to_target(
                source_groups,
                group_band_counts=group_band_counts,
                bounds=bounds,
                crs=crs,
                resolution=target_resolution,
                resampling=Resampling.nearest,
                nodata=np.nan,
                target_template=target_template,
            )
        except _DATA_READ_ERRORS:
            logger.debug("%s direct VRT payload path unavailable; falling back to array merge.", self._source_name, exc_info=True)
            return None

        layer_values = np.asarray(stacked.values, dtype=np.float32)
        params_values, unc_values = self._allocate_spatial_payload_arrays(requested, target_template)

        offset = 0
        for band_index in range(len(requested)):
            params_values[band_index] = apply_scale_and_mask(layer_values[offset : offset + 3], param_attrs[band_index])
            offset += 3
            qa_values = apply_scale_and_mask(layer_values[offset], qa_attrs[band_index])
            offset += 1
            unc_values[band_index] = self._qa_values_to_uncertainty(qa_values)

        params = self._target_array_like(
            target_template,
            params_values,
            dims=("band", "parameter", "y", "x"),
            coords={
                "band": xr.IndexVariable("band", band_coords),
                "parameter": xr.IndexVariable("parameter", ["f0", "f1", "f2"]),
                "y": target_template.coords["y"],
                "x": target_template.coords["x"],
            },
        )
        unc = self._target_array_like(
            target_template,
            unc_values,
            dims=("band", "y", "x"),
            coords={
                "band": xr.IndexVariable("band", band_coords),
                "y": target_template.coords["y"],
                "x": target_template.coords["x"],
            },
        )
        return self._pack_payload_stack(params, unc)

    def _load_temporal_payload_vrt(
        self,
        grouped_paths: dict[np.datetime64, list[Path]],
        *,
        requested: Sequence[_RequestedBandSpec],
        bounds: tuple[float, float, float, float],
        crs: str,
        target_resolution: float,
        time_axis: np.ndarray,
        target_template: xr.DataArray,
    ) -> tuple[np.ndarray, np.ndarray] | None:
        if not gdal_available():
            return None

        params_values, unc_values = self._allocate_temporal_payload_arrays(
            time_axis,
            requested,
            y_size=int(target_template.sizes["y"]),
            x_size=int(target_template.sizes["x"]),
        )
        available_days = [(time_index, grouped_paths[day]) for time_index, day in enumerate(time_axis) if grouped_paths.get(day)]
        if not available_days:
            return params_values, unc_values

        dataset_names = tuple(
            dict.fromkeys(
                [
                    dataset_name
                    for _band_coord, product_band in requested
                    for dataset_name in (product_band.parameter_dataset, product_band.qa_dataset or "")
                ]
            )
        )
        try:
            dataset_attrs = self._read_named_dataset_attrs(available_days[0][1][0], dataset_names)
            source_groups: list[list[str]] = []
            group_band_counts: list[int] = []
            entries: list[tuple[str, int, int, str]] = []

            logger.info(
                "%s direct temporal VRT payload: available_days=%d requested_bands=%d",
                self._source_name,
                len(available_days),
                len(requested),
            )

            for time_index, day_paths in available_days:
                logger.debug(
                    "%s direct temporal VRT payload day=%s tiles=%d",
                    self._source_name,
                    str(time_axis[time_index]),
                    len(day_paths),
                )
                for band_index, (_band_coord, product_band) in enumerate(requested):
                    qa_name = product_band.qa_dataset or ""
                    source_groups.append(
                        [resolve_gdal_subdataset_path(path, product_band.parameter_dataset) for path in day_paths]
                    )
                    group_band_counts.append(3)
                    entries.append(("params", time_index, band_index, product_band.parameter_dataset))

                    source_groups.append([resolve_gdal_subdataset_path(path, qa_name) for path in day_paths])
                    group_band_counts.append(1)
                    entries.append(("qa", time_index, band_index, qa_name))

            expected_layers = int(sum(group_band_counts))
            logger.info(
                "%s direct temporal VRT payload: source_groups=%d expected_layers=%d",
                self._source_name,
                len(source_groups),
                expected_layers,
            )
            if not self._allow_one_shot_temporal_vrt(
                target_template=target_template,
                expected_layers=expected_layers,
                label="direct temporal VRT payload",
            ):
                return None
            stacked = read_virtual_stack_to_target(
                source_groups,
                group_band_counts=group_band_counts,
                bounds=bounds,
                crs=crs,
                resolution=target_resolution,
                resampling=Resampling.nearest,
                nodata=np.nan,
                target_template=target_template,
            )
        except _DATA_READ_ERRORS:
            logger.debug(
                "%s direct temporal VRT payload path unavailable; falling back to per-day array merge.",
                self._source_name,
                exc_info=True,
            )
            return None

        layer_values = np.asarray(stacked.values, dtype=np.float32)
        logger.info("%s direct temporal VRT payload: stacked shape=%s", self._source_name, tuple(layer_values.shape))
        offset = 0
        for kind, time_index, band_index, dataset_name in entries:
            if kind == "params":
                params_values[time_index, band_index] = apply_scale_and_mask(
                    layer_values[offset : offset + 3],
                    dataset_attrs[dataset_name],
                )
                offset += 3
                continue

            qa_values = apply_scale_and_mask(layer_values[offset], dataset_attrs[dataset_name])
            offset += 1
            unc_values[time_index, band_index] = self._qa_values_to_uncertainty(qa_values)

        return params_values, unc_values

    def _load_temporal_payload_vrt_from_daily_payload_vrts(
        self,
        grouped_paths: dict[np.datetime64, list[Path]],
        *,
        requested: Sequence[_RequestedBandSpec],
        bounds: tuple[float, float, float, float],
        crs: str,
        target_resolution: float,
        time_axis: np.ndarray,
        target_template: xr.DataArray,
    ) -> tuple[np.ndarray, np.ndarray] | None:
        if not gdal_available():
            return None

        params_values, unc_values = self._allocate_temporal_payload_arrays(
            time_axis,
            requested,
            y_size=int(target_template.sizes["y"]),
            x_size=int(target_template.sizes["x"]),
        )
        available_days = [(time_index, grouped_paths[day]) for time_index, day in enumerate(time_axis) if grouped_paths.get(day)]
        if not available_days:
            return params_values, unc_values

        dataset_names = tuple(
            dict.fromkeys(
                [
                    dataset_name
                    for _band_coord, product_band in requested
                    for dataset_name in (product_band.parameter_dataset, product_band.qa_dataset or "")
                ]
            )
        )
        day_specs: list[Any] = []
        try:
            dataset_attrs = self._read_named_dataset_attrs(available_days[0][1][0], dataset_names)
            param_attrs = [dataset_attrs[product_band.parameter_dataset] for _band_coord, product_band in requested]
            qa_attrs = [dataset_attrs[product_band.qa_dataset or ""] for _band_coord, product_band in requested]
            target_bounds = _target_bounds_from_template(target_template, resolution=target_resolution)
            target_crs = str(target_template.rio.crs)
            source_groups: list[list[str]] = []
            group_band_counts: list[int] = []

            for time_index, day_paths in available_days:
                logger.debug(
                    "%s daily-payload temporal VRT day=%s tiles=%d",
                    self._source_name,
                    str(time_axis[time_index]),
                    len(day_paths),
                )
                day_source_groups: list[list[str]] = []
                day_group_band_counts: list[int] = []
                for _band_coord, product_band in requested:
                    qa_name = product_band.qa_dataset or ""
                    day_source_groups.append(
                        [resolve_gdal_subdataset_path(path, product_band.parameter_dataset) for path in day_paths]
                    )
                    day_group_band_counts.append(3)
                    day_source_groups.append([resolve_gdal_subdataset_path(path, qa_name) for path in day_paths])
                    day_group_band_counts.append(1)

                day_spec = _build_virtual_stack_vrt(
                    day_source_groups,
                    group_band_counts=day_group_band_counts,
                    target_bounds=target_bounds,
                    target_crs=target_crs,
                )
                day_specs.append(day_spec)
                source_groups.append([day_spec.path])
                group_band_counts.append(day_spec.expected_layers)

            expected_layers = int(sum(group_band_counts))
            logger.info(
                "%s daily-payload temporal VRT: available_days=%d requested_bands=%d expected_layers=%d",
                self._source_name,
                len(available_days),
                len(requested),
                expected_layers,
            )
            if not self._allow_one_shot_temporal_vrt(
                target_template=target_template,
                expected_layers=expected_layers,
                label="daily-payload temporal VRT",
            ):
                return None
            stacked = read_virtual_stack_to_target(
                source_groups,
                group_band_counts=group_band_counts,
                bounds=bounds,
                crs=crs,
                resolution=target_resolution,
                resampling=Resampling.nearest,
                nodata=np.nan,
                target_template=target_template,
            )
        except _DATA_READ_ERRORS:
            logger.debug(
                "%s daily-payload temporal VRT path unavailable; falling back to per-day array merge.",
                self._source_name,
                exc_info=True,
            )
            return None
        finally:
            if day_specs:
                from osgeo import gdal  # type: ignore[import-untyped]

                for day_spec in day_specs:
                    for path in day_spec.cleanup_paths:
                        gdal.Unlink(path)

        layer_values = np.asarray(stacked.values, dtype=np.float32)
        offset = 0
        for (time_index, _day_paths), day_spec in zip(available_days, day_specs, strict=True):
            day_layers = layer_values[offset : offset + day_spec.expected_layers]
            day_offset = 0
            for band_index in range(len(requested)):
                params_values[time_index, band_index] = apply_scale_and_mask(
                    day_layers[day_offset : day_offset + 3],
                    param_attrs[band_index],
                )
                day_offset += 3
                qa_values = apply_scale_and_mask(day_layers[day_offset], qa_attrs[band_index])
                day_offset += 1
                unc_values[time_index, band_index] = self._qa_values_to_uncertainty(qa_values)
            offset += day_spec.expected_layers

        return params_values, unc_values


class MCD43EarthAccessProvider(_StackParameterProvider):
    """MODIS MCD43 BRDF provider using real MCD43A1 kernel parameters."""

    product_key = "mcd43_brdf"
    _source_name = "MCD43"
    _rsrf_sensor_unit_id = "terra_modis"
    _product_bands = (
        ProductBandDefinition("Band1", 645.0, 50.0, "BRDF_Albedo_Parameters_Band1", "BRDF_Albedo_Band_Mandatory_Quality_Band1", rsrf_band_id="B1"),
        ProductBandDefinition("Band2", 858.5, 35.0, "BRDF_Albedo_Parameters_Band2", "BRDF_Albedo_Band_Mandatory_Quality_Band2", rsrf_band_id="B2"),
        ProductBandDefinition("Band3", 469.0, 20.0, "BRDF_Albedo_Parameters_Band3", "BRDF_Albedo_Band_Mandatory_Quality_Band3", rsrf_band_id="B3"),
        ProductBandDefinition("Band4", 555.0, 20.0, "BRDF_Albedo_Parameters_Band4", "BRDF_Albedo_Band_Mandatory_Quality_Band4", rsrf_band_id="B4"),
        ProductBandDefinition("Band5", 1240.0, 20.0, "BRDF_Albedo_Parameters_Band5", "BRDF_Albedo_Band_Mandatory_Quality_Band5", rsrf_band_id="B5"),
        ProductBandDefinition("Band6", 1640.0, 24.0, "BRDF_Albedo_Parameters_Band6", "BRDF_Albedo_Band_Mandatory_Quality_Band6", rsrf_band_id="B6"),
        ProductBandDefinition("Band7", 2130.0, 50.0, "BRDF_Albedo_Parameters_Band7", "BRDF_Albedo_Band_Mandatory_Quality_Band7", rsrf_band_id="B7"),
    )
    _legacy_band_map = {index + 1: band.label for index, band in enumerate(_product_bands)}


class VNP43EarthAccessProvider(_StackParameterProvider):
    """VIIRS VNP43 BRDF provider using real VNP43MA1 kernel parameters."""

    product_key = "vnp43_brdf"
    _source_name = "VNP43"
    _rsrf_sensor_unit_id = "snpp_viirs"
    _read_dataset = staticmethod(read_hdf5_dataset)
    _product_bands = (
        ProductBandDefinition(
            "M1",
            412.0,
            20.0,
            "HDFEOS/GRIDS/VIIRS_Grid_BRDF/Data Fields/BRDF_Albedo_Parameters_M1",
            "HDFEOS/GRIDS/VIIRS_Grid_BRDF/Data Fields/BRDF_Albedo_Band_Mandatory_Quality_M1",
        ),
        ProductBandDefinition(
            "M2",
            445.0,
            18.0,
            "HDFEOS/GRIDS/VIIRS_Grid_BRDF/Data Fields/BRDF_Albedo_Parameters_M2",
            "HDFEOS/GRIDS/VIIRS_Grid_BRDF/Data Fields/BRDF_Albedo_Band_Mandatory_Quality_M2",
        ),
        ProductBandDefinition(
            "M3",
            488.0,
            20.0,
            "HDFEOS/GRIDS/VIIRS_Grid_BRDF/Data Fields/BRDF_Albedo_Parameters_M3",
            "HDFEOS/GRIDS/VIIRS_Grid_BRDF/Data Fields/BRDF_Albedo_Band_Mandatory_Quality_M3",
        ),
        ProductBandDefinition(
            "M4",
            555.0,
            20.0,
            "HDFEOS/GRIDS/VIIRS_Grid_BRDF/Data Fields/BRDF_Albedo_Parameters_M4",
            "HDFEOS/GRIDS/VIIRS_Grid_BRDF/Data Fields/BRDF_Albedo_Band_Mandatory_Quality_M4",
        ),
        ProductBandDefinition(
            "M5",
            672.0,
            20.0,
            "HDFEOS/GRIDS/VIIRS_Grid_BRDF/Data Fields/BRDF_Albedo_Parameters_M5",
            "HDFEOS/GRIDS/VIIRS_Grid_BRDF/Data Fields/BRDF_Albedo_Band_Mandatory_Quality_M5",
        ),
        ProductBandDefinition(
            "M7",
            865.0,
            39.0,
            "HDFEOS/GRIDS/VIIRS_Grid_BRDF/Data Fields/BRDF_Albedo_Parameters_M7",
            "HDFEOS/GRIDS/VIIRS_Grid_BRDF/Data Fields/BRDF_Albedo_Band_Mandatory_Quality_M7",
        ),
        ProductBandDefinition(
            "M8",
            1240.0,
            20.0,
            "HDFEOS/GRIDS/VIIRS_Grid_BRDF/Data Fields/BRDF_Albedo_Parameters_M8",
            "HDFEOS/GRIDS/VIIRS_Grid_BRDF/Data Fields/BRDF_Albedo_Band_Mandatory_Quality_M8",
        ),
        ProductBandDefinition(
            "M10",
            1610.0,
            60.0,
            "HDFEOS/GRIDS/VIIRS_Grid_BRDF/Data Fields/BRDF_Albedo_Parameters_M10",
            "HDFEOS/GRIDS/VIIRS_Grid_BRDF/Data Fields/BRDF_Albedo_Band_Mandatory_Quality_M10",
        ),
        ProductBandDefinition(
            "M11",
            2250.0,
            50.0,
            "HDFEOS/GRIDS/VIIRS_Grid_BRDF/Data Fields/BRDF_Albedo_Parameters_M11",
            "HDFEOS/GRIDS/VIIRS_Grid_BRDF/Data Fields/BRDF_Albedo_Band_Mandatory_Quality_M11",
        ),
    )
    _legacy_band_map = {
        1: "M5",
        2: "M7",
        3: "M3",
        4: "M4",
        5: "M8",
        6: "M10",
        7: "M11",
    }


class MCD19EarthAccessProvider(_EarthAccessBRDFProvider):
    """MODIS MAIAC BRDF provider using real MCD19A3 kernel parameters."""

    product_key = "mcd19_brdf"
    _source_name = "MCD19"
    _rsrf_sensor_unit_id = "terra_modis"
    _product_bands = (
        ProductBandDefinition("Band1", 645.0, 50.0, "Kiso_Band1", "Status_QA", rsrf_band_id="B1"),
        ProductBandDefinition("Band2", 858.5, 35.0, "Kiso_Band2", "Status_QA", rsrf_band_id="B2"),
        ProductBandDefinition("Band3", 469.0, 20.0, "Kiso_Band3", "Status_QA", rsrf_band_id="B3"),
        ProductBandDefinition("Band4", 555.0, 20.0, "Kiso_Band4", "Status_QA", rsrf_band_id="B4"),
        ProductBandDefinition("Band5", 1240.0, 20.0, "Kiso_Band5", "Status_QA", rsrf_band_id="B5"),
        ProductBandDefinition("Band6", 1640.0, 24.0, "Kiso_Band6", "Status_QA", rsrf_band_id="B6"),
        ProductBandDefinition("Band7", 2130.0, 50.0, "Kiso_Band7", "Status_QA", rsrf_band_id="B7"),
        ProductBandDefinition("Band8", 412.0, 20.0, "Kiso_Band8", "Status_QA", rsrf_band_id="B8"),
    )
    _legacy_band_map = {index + 1: band.label for index, band in enumerate(_product_bands)}

    def _load_native_band_stack(
        self,
        path: str | Path,
        product_band: ProductBandDefinition,
    ) -> tuple[tuple[xr.DataArray, xr.DataArray, xr.DataArray], xr.DataArray]:
        label = product_band.label.replace("Band", "")
        dataset_names = (
            f"Kiso_Band{label}",
            f"Kvol_Band{label}",
            f"Kgeo_Band{label}",
            "Status_QA",
        )
        if read_hdf4_dataset is _ORIGINAL_HDF4_READER:
            datasets = read_hdf4_datasets(path, dataset_names)
        else:
            datasets = {
                dataset_name: read_hdf4_dataset(path, dataset_name)
                for dataset_name in dataset_names
            }
        f0_raw, f0_attrs = datasets[dataset_names[0]]
        f1_raw, f1_attrs = datasets[dataset_names[1]]
        f2_raw, f2_attrs = datasets[dataset_names[2]]
        qa_raw, qa_attrs = datasets[dataset_names[3]]

        f0 = make_native_grid_dataarray(apply_scale_and_mask(f0_raw, f0_attrs), granule_path=path)
        f1 = make_native_grid_dataarray(apply_scale_and_mask(f1_raw, f1_attrs), granule_path=path)
        f2 = make_native_grid_dataarray(apply_scale_and_mask(f2_raw, f2_attrs), granule_path=path)
        qa = make_native_grid_dataarray(apply_scale_and_mask(qa_raw, qa_attrs), granule_path=path)
        return (f0, f1, f2), qa

    def _load_native_requested_stack(
        self,
        path: str | Path,
        requested: Sequence[_RequestedBandSpec],
    ) -> tuple[xr.DataArray, xr.DataArray]:
        band_coords = [band_coord for band_coord, _product_band in requested]
        dataset_names = ["Status_QA"]
        labels = [product_band.label.replace("Band", "") for _band_coord, product_band in requested]
        for label in labels:
            dataset_names.extend(
                [
                    f"Kiso_Band{label}",
                    f"Kvol_Band{label}",
                    f"Kgeo_Band{label}",
                ]
            )
        unique_dataset_names = tuple(dict.fromkeys(dataset_names))
        if read_hdf4_dataset is _ORIGINAL_HDF4_READER:
            datasets = read_hdf4_datasets(path, unique_dataset_names)
        else:
            datasets = {
                dataset_name: read_hdf4_dataset(path, dataset_name)
                for dataset_name in unique_dataset_names
            }

        qa_raw, qa_attrs = datasets["Status_QA"]
        qa_da = make_native_grid_dataarray(apply_scale_and_mask(qa_raw, qa_attrs), granule_path=path)
        param_tiles = []
        qa_tiles = []
        for label in labels:
            f0_raw, f0_attrs = datasets[f"Kiso_Band{label}"]
            f1_raw, f1_attrs = datasets[f"Kvol_Band{label}"]
            f2_raw, f2_attrs = datasets[f"Kgeo_Band{label}"]
            param_tiles.append(
                self._stack_parameter_cube(
                    (
                        make_native_grid_dataarray(apply_scale_and_mask(f0_raw, f0_attrs), granule_path=path),
                        make_native_grid_dataarray(apply_scale_and_mask(f1_raw, f1_attrs), granule_path=path),
                        make_native_grid_dataarray(apply_scale_and_mask(f2_raw, f2_attrs), granule_path=path),
                    )
                )
            )
            qa_tiles.append(qa_da)

        return (
            xr.concat(param_tiles, dim=xr.IndexVariable("band", band_coords)).transpose("band", "parameter", "y", "x"),
            xr.concat(qa_tiles, dim=xr.IndexVariable("band", band_coords)).transpose("band", "y", "x"),
        )

    def _load_requested_payload_vrt(
        self,
        paths: Sequence[str | Path],
        *,
        requested: Sequence[_RequestedBandSpec],
        bounds: tuple[float, float, float, float],
        crs: str,
        target_resolution: float,
        target_template: xr.DataArray,
    ) -> xr.DataArray | None:
        if not gdal_available():
            return None

        labels = [product_band.label.replace("Band", "") for _band_coord, product_band in requested]
        dataset_names = ["Status_QA"]
        for label in labels:
            dataset_names.extend(
                [
                    f"Kiso_Band{label}",
                    f"Kvol_Band{label}",
                    f"Kgeo_Band{label}",
                ]
            )
        unique_dataset_names = tuple(dict.fromkeys(dataset_names))

        try:
            dataset_attrs = read_hdf4_datasets_attrs(paths[0], unique_dataset_names)
            source_groups: list[list[str]] = []
            group_band_counts: list[int] = []
            for label in labels:
                for dataset_name in (
                    f"Kiso_Band{label}",
                    f"Kvol_Band{label}",
                    f"Kgeo_Band{label}",
                ):
                    source_groups.append([resolve_gdal_subdataset_path(path, dataset_name) for path in paths])
                    group_band_counts.append(1)

            source_groups.append([resolve_gdal_subdataset_path(path, "Status_QA") for path in paths])
            group_band_counts.append(1)

            stacked = read_virtual_stack_to_target(
                source_groups,
                group_band_counts=group_band_counts,
                bounds=bounds,
                crs=crs,
                resolution=target_resolution,
                resampling=Resampling.nearest,
                nodata=np.nan,
                target_template=target_template,
            )
        except _DATA_READ_ERRORS:
            logger.debug("%s direct VRT payload path unavailable; falling back to array merge.", self._source_name, exc_info=True)
            return None

        layer_values = np.asarray(stacked.values, dtype=np.float32)
        params_values, _unused_unc_values = self._allocate_spatial_payload_arrays(requested, target_template)
        offset = 0
        for band_index, label in enumerate(labels):
            params_values[band_index, 0] = apply_scale_and_mask(layer_values[offset], dataset_attrs[f"Kiso_Band{label}"])
            offset += 1
            params_values[band_index, 1] = apply_scale_and_mask(layer_values[offset], dataset_attrs[f"Kvol_Band{label}"])
            offset += 1
            params_values[band_index, 2] = apply_scale_and_mask(layer_values[offset], dataset_attrs[f"Kgeo_Band{label}"])
            offset += 1

        qa_values = apply_scale_and_mask(layer_values[offset], dataset_attrs["Status_QA"])
        unc_base = self._qa_values_to_uncertainty(qa_values)
        unc_values = np.repeat(unc_base[np.newaxis, :, :], len(requested), axis=0)
        band_coords = [band_coord for band_coord, _product_band in requested]

        params = self._target_array_like(
            target_template,
            params_values,
            dims=("band", "parameter", "y", "x"),
            coords={
                "band": xr.IndexVariable("band", band_coords),
                "parameter": xr.IndexVariable("parameter", ["f0", "f1", "f2"]),
                "y": target_template.coords["y"],
                "x": target_template.coords["x"],
            },
        )
        unc = self._target_array_like(
            target_template,
            unc_values,
            dims=("band", "y", "x"),
            coords={
                "band": xr.IndexVariable("band", band_coords),
                "y": target_template.coords["y"],
                "x": target_template.coords["x"],
            },
        )
        return self._pack_payload_stack(params, unc)

    def _load_temporal_payload_vrt(
        self,
        grouped_paths: dict[np.datetime64, list[Path]],
        *,
        requested: Sequence[_RequestedBandSpec],
        bounds: tuple[float, float, float, float],
        crs: str,
        target_resolution: float,
        time_axis: np.ndarray,
        target_template: xr.DataArray,
    ) -> tuple[np.ndarray, np.ndarray] | None:
        if not gdal_available():
            return None

        params_values, unc_values = self._allocate_temporal_payload_arrays(
            time_axis,
            requested,
            y_size=int(target_template.sizes["y"]),
            x_size=int(target_template.sizes["x"]),
        )
        available_days = [(time_index, grouped_paths[day]) for time_index, day in enumerate(time_axis) if grouped_paths.get(day)]
        if not available_days:
            return params_values, unc_values

        labels = [product_band.label.replace("Band", "") for _band_coord, product_band in requested]
        dataset_names = ["Status_QA"]
        for label in labels:
            dataset_names.extend([f"Kiso_Band{label}", f"Kvol_Band{label}", f"Kgeo_Band{label}"])
        unique_dataset_names = tuple(dict.fromkeys(dataset_names))

        try:
            dataset_attrs = read_hdf4_datasets_attrs(available_days[0][1][0], unique_dataset_names)
            source_groups: list[list[str]] = []
            group_band_counts: list[int] = []
            entries: list[tuple[str, int, int, int, str]] = []

            logger.info(
                "%s direct temporal VRT payload: available_days=%d requested_bands=%d shared_qa=yes",
                self._source_name,
                len(available_days),
                len(requested),
            )

            for time_index, day_paths in available_days:
                logger.debug(
                    "%s direct temporal VRT payload day=%s tiles=%d",
                    self._source_name,
                    str(time_axis[time_index]),
                    len(day_paths),
                )
                for band_index, label in enumerate(labels):
                    for parameter_index, dataset_name in enumerate(
                        (
                            f"Kiso_Band{label}",
                            f"Kvol_Band{label}",
                            f"Kgeo_Band{label}",
                        )
                    ):
                        source_groups.append([resolve_gdal_subdataset_path(path, dataset_name) for path in day_paths])
                        group_band_counts.append(1)
                        entries.append(("param", time_index, band_index, parameter_index, dataset_name))

                source_groups.append([resolve_gdal_subdataset_path(path, "Status_QA") for path in day_paths])
                group_band_counts.append(1)
                entries.append(("qa", time_index, -1, -1, "Status_QA"))

            expected_layers = int(sum(group_band_counts))
            logger.info(
                "%s direct temporal VRT payload: source_groups=%d expected_layers=%d",
                self._source_name,
                len(source_groups),
                expected_layers,
            )
            if not self._allow_one_shot_temporal_vrt(
                target_template=target_template,
                expected_layers=expected_layers,
                label="direct temporal VRT payload",
            ):
                return None
            stacked = read_virtual_stack_to_target(
                source_groups,
                group_band_counts=group_band_counts,
                bounds=bounds,
                crs=crs,
                resolution=target_resolution,
                resampling=Resampling.nearest,
                nodata=np.nan,
                target_template=target_template,
            )
        except _DATA_READ_ERRORS:
            logger.debug(
                "%s direct temporal VRT payload path unavailable; falling back to per-day array merge.",
                self._source_name,
                exc_info=True,
            )
            return None

        layer_values = np.asarray(stacked.values, dtype=np.float32)
        logger.info("%s direct temporal VRT payload: stacked shape=%s", self._source_name, tuple(layer_values.shape))
        offset = 0
        for kind, time_index, band_index, parameter_index, dataset_name in entries:
            if kind == "param":
                params_values[time_index, band_index, parameter_index] = apply_scale_and_mask(
                    layer_values[offset],
                    dataset_attrs[dataset_name],
                )
                offset += 1
                continue

            qa_values = apply_scale_and_mask(layer_values[offset], dataset_attrs[dataset_name])
            offset += 1
            unc_base = self._qa_values_to_uncertainty(qa_values)
            unc_values[time_index] = np.repeat(unc_base[np.newaxis, :, :], len(requested), axis=0)

        return params_values, unc_values

    def _load_temporal_payload_vrt_from_daily_payload_vrts(
        self,
        grouped_paths: dict[np.datetime64, list[Path]],
        *,
        requested: Sequence[_RequestedBandSpec],
        bounds: tuple[float, float, float, float],
        crs: str,
        target_resolution: float,
        time_axis: np.ndarray,
        target_template: xr.DataArray,
    ) -> tuple[np.ndarray, np.ndarray] | None:
        if not gdal_available():
            return None

        params_values, unc_values = self._allocate_temporal_payload_arrays(
            time_axis,
            requested,
            y_size=int(target_template.sizes["y"]),
            x_size=int(target_template.sizes["x"]),
        )
        available_days = [(time_index, grouped_paths[day]) for time_index, day in enumerate(time_axis) if grouped_paths.get(day)]
        if not available_days:
            return params_values, unc_values

        labels = [product_band.label.replace("Band", "") for _band_coord, product_band in requested]
        dataset_names = ["Status_QA"]
        for label in labels:
            dataset_names.extend([f"Kiso_Band{label}", f"Kvol_Band{label}", f"Kgeo_Band{label}"])
        unique_dataset_names = tuple(dict.fromkeys(dataset_names))

        day_specs: list[Any] = []
        try:
            dataset_attrs = read_hdf4_datasets_attrs(available_days[0][1][0], unique_dataset_names)
            target_bounds = _target_bounds_from_template(target_template, resolution=target_resolution)
            target_crs = str(target_template.rio.crs)
            source_groups: list[list[str]] = []
            group_band_counts: list[int] = []

            for time_index, day_paths in available_days:
                logger.debug(
                    "%s daily-payload temporal VRT day=%s tiles=%d",
                    self._source_name,
                    str(time_axis[time_index]),
                    len(day_paths),
                )
                day_source_groups: list[list[str]] = []
                day_group_band_counts: list[int] = []
                for label in labels:
                    for dataset_name in (
                        f"Kiso_Band{label}",
                        f"Kvol_Band{label}",
                        f"Kgeo_Band{label}",
                    ):
                        day_source_groups.append([resolve_gdal_subdataset_path(path, dataset_name) for path in day_paths])
                        day_group_band_counts.append(1)

                day_source_groups.append([resolve_gdal_subdataset_path(path, "Status_QA") for path in day_paths])
                day_group_band_counts.append(1)

                day_spec = _build_virtual_stack_vrt(
                    day_source_groups,
                    group_band_counts=day_group_band_counts,
                    target_bounds=target_bounds,
                    target_crs=target_crs,
                )
                day_specs.append(day_spec)
                source_groups.append([day_spec.path])
                group_band_counts.append(day_spec.expected_layers)

            expected_layers = int(sum(group_band_counts))
            logger.info(
                "%s daily-payload temporal VRT: available_days=%d requested_bands=%d shared_qa=yes expected_layers=%d",
                self._source_name,
                len(available_days),
                len(requested),
                expected_layers,
            )
            if not self._allow_one_shot_temporal_vrt(
                target_template=target_template,
                expected_layers=expected_layers,
                label="daily-payload temporal VRT",
            ):
                return None
            stacked = read_virtual_stack_to_target(
                source_groups,
                group_band_counts=group_band_counts,
                bounds=bounds,
                crs=crs,
                resolution=target_resolution,
                resampling=Resampling.nearest,
                nodata=np.nan,
                target_template=target_template,
            )
        except _DATA_READ_ERRORS:
            logger.debug(
                "%s daily-payload temporal VRT path unavailable; falling back to per-day array merge.",
                self._source_name,
                exc_info=True,
            )
            return None
        finally:
            if day_specs:
                from osgeo import gdal

                for day_spec in day_specs:
                    for path in day_spec.cleanup_paths:
                        gdal.Unlink(path)

        layer_values = np.asarray(stacked.values, dtype=np.float32)
        offset = 0
        for (time_index, _day_paths), day_spec in zip(available_days, day_specs, strict=True):
            day_layers = layer_values[offset : offset + day_spec.expected_layers]
            day_offset = 0
            for band_index, label in enumerate(labels):
                params_values[time_index, band_index, 0] = apply_scale_and_mask(
                    day_layers[day_offset],
                    dataset_attrs[f"Kiso_Band{label}"],
                )
                day_offset += 1
                params_values[time_index, band_index, 1] = apply_scale_and_mask(
                    day_layers[day_offset],
                    dataset_attrs[f"Kvol_Band{label}"],
                )
                day_offset += 1
                params_values[time_index, band_index, 2] = apply_scale_and_mask(
                    day_layers[day_offset],
                    dataset_attrs[f"Kgeo_Band{label}"],
                )
                day_offset += 1

            qa_values = apply_scale_and_mask(day_layers[day_offset], dataset_attrs["Status_QA"])
            unc_base = self._qa_values_to_uncertainty(qa_values)
            unc_values[time_index] = np.repeat(unc_base[np.newaxis, :, :], len(requested), axis=0)
            offset += day_spec.expected_layers

        return params_values, unc_values
