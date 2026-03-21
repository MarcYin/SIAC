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
    build_earthaccess_runtime,
    constant_target_band_array,
    earthaccess_cache_dir,
    merge_reprojected_tiles,
    select_candidate_paths,
    target_grid_coords,
)
from siac.adapters.earthdata_common import (
    ProductBandDefinition,
    apply_scale_and_mask,
    make_native_grid_dataarray,
    parse_granule_date,
    read_hdf4_dataset,
    read_hdf5_dataset,
)
from siac.domain import SensorBand
from siac.runtime import BRDFKernelWeights

if TYPE_CHECKING:
    from collections.abc import Sequence
    from datetime import date, datetime

    from siac.adapters.data.earthaccess_catalog import EarthAccessCatalog

logger = logging.getLogger(__name__)

_BEST_QA_REFLECTANCE_UNCERTAINTY = 0.015
_QA_UNCERTAINTY_POWER = 1.6
_RequestedBand: TypeAlias = int | SensorBand
_RequestedBandCoord: TypeAlias = int | str
_RequestedBandSpec: TypeAlias = tuple[_RequestedBandCoord, ProductBandDefinition]


def _granule_day(granule: object) -> np.datetime64 | None:
    """Best-effort extraction of a daily timestamp from an Earthaccess granule object."""
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
    temporal = umm.get("TemporalExtent", {}) if isinstance(umm, dict) else {}
    if isinstance(temporal, dict):
        range_dt = temporal.get("RangeDateTime", {})
        if isinstance(range_dt, dict):
            for key in ("BeginningDateTime", "EndingDateTime"):
                value = range_dt.get(key)
                if isinstance(value, str):
                    candidates.append(value)
        single_dt = temporal.get("SingleDateTime")
        if isinstance(single_dt, str):
            candidates.append(single_dt)

    if isinstance(umm, dict):
        for key in ("GranuleUR", "ProducerGranuleId"):
            value = umm.get(key)
            if isinstance(value, str):
                candidates.append(value)
    if isinstance(meta, dict):
        for key in ("native-id", "producer-granule-id"):
            value = meta.get(key)
            if isinstance(value, str):
                candidates.append(value)

    for value in candidates:
        try:
            return np.datetime64(value[:10], "D")
        except Exception:
            pass
        try:
            return np.datetime64(parse_granule_date(Path(value)).date(), "D")
        except Exception:
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
        self.max_granules = max(1, int(max_granules))

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
                rsrf_band_id=band.label,
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
                except Exception as exc:  # pragma: no cover - external/system dependent
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
                    return self._load_temporal_from_granules(
                        paths,
                        requested=requested,
                        bounds=bounds,
                        crs=crs,
                        target_resolution=target_resolution,
                        time_axis=time_axis,
                    )
                except Exception as exc:  # pragma: no cover - external/system dependent
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
        request_specs: list[tuple[datetime, np.ndarray]] = []
        for obs_time, temporal_window, sample_dates in zip(obs_times, temporal_windows, sample_date_sets, strict=True):
            time_axis = (
                self._coerce_sample_time_axis(sample_dates)
                if sample_dates is not None
                else self._time_axis(obs_time, temporal_window)
            )
            request_specs.append((obs_time, time_axis))

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
                    for obs_time, time_axis in request_specs:
                        paths = self._select_candidate_paths(
                            downloaded,
                            obs_time,
                            bounds,
                            crs,
                            sample_dates=time_axis,
                        )
                        if paths:
                            outputs.append(
                                self._load_temporal_from_granules(
                                    paths,
                                    requested=requested,
                                    bounds=bounds,
                                    crs=crs,
                                    target_resolution=target_resolution,
                                    time_axis=time_axis,
                                )
                            )
                        else:
                            outputs.append(
                                self._default_temporal_weights(
                                    bounds,
                                    target_resolution,
                                    [coord for coord, _band in requested],
                                    time_axis,
                                )
                            )
                    return outputs
                except Exception as exc:  # pragma: no cover - external/system dependent
                    logger.warning(
                        "%s batched temporal BRDF granule parsing failed; using defaults (%s)",
                        self._source_name,
                        exc,
                    )

        return [
            self._default_temporal_weights(
                bounds,
                target_resolution,
                [coord for coord, _band in requested],
                time_axis,
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
    ) -> list[object]:
        return cast(
            "list[object]",
            self.source.search_granules(
                short_name=short_name,
                bounds=bounds,
                crs=crs,
                temporal=temporal,
                provider=self.provider,
                count=self.max_granules,
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
        dest = earthaccess_cache_dir(self.cache_dir, short_name)
        return cast("list[Path]", self.source.download_granules(granules, dest))

    @staticmethod
    def _merge_search_batches(
        request_specs: Sequence[tuple[datetime, np.ndarray]],
        temporal_windows: Sequence[int],
    ) -> list[tuple[Any, Any, np.ndarray]]:
        windows: list[tuple[Any, Any, np.ndarray]] = []
        for (obs_time, sample_dates), temporal_window in zip(request_specs, temporal_windows, strict=True):
            temporal = EarthAccessSource.temporal_window(obs_time, temporal_window)
            start_day = np.datetime64(temporal[0][:10], "D")
            end_day = np.datetime64(temporal[1][:10], "D")
            windows.append((start_day, end_day, sample_dates))

        windows.sort(key=lambda item: item[0])
        batches: list[tuple[Any, Any, np.ndarray]] = []
        for start_day, end_day, sample_dates in windows:
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
        f0_list: list[xr.DataArray] = []
        f1_list: list[xr.DataArray] = []
        f2_list: list[xr.DataArray] = []
        f0_unc_list: list[xr.DataArray] = []
        f1_unc_list: list[xr.DataArray] = []
        f2_unc_list: list[xr.DataArray] = []
        reflectance_unc_list: list[xr.DataArray] = []

        for band_coord, product_band in requested:
            param_tiles: list[tuple[xr.DataArray, xr.DataArray, xr.DataArray]] = []
            qa_tiles: list[xr.DataArray] = []
            for path in paths:
                params, qa = self._load_native_band_stack(path, product_band)
                qa_tiles.append(qa)
                param_tiles.append(params)

            f0 = self._merge_reprojected_tiles(
                [params[0] for params in param_tiles],
                bounds=bounds,
                crs=crs,
                target_resolution=target_resolution,
                resampling=Resampling.bilinear,
                nodata=np.nan,
            )
            f1 = self._merge_reprojected_tiles(
                [params[1] for params in param_tiles],
                bounds=bounds,
                crs=crs,
                target_resolution=target_resolution,
                resampling=Resampling.bilinear,
                nodata=np.nan,
            )
            f2 = self._merge_reprojected_tiles(
                [params[2] for params in param_tiles],
                bounds=bounds,
                crs=crs,
                target_resolution=target_resolution,
                resampling=Resampling.bilinear,
                nodata=np.nan,
            )
            qa = self._merge_reprojected_tiles(
                qa_tiles,
                bounds=bounds,
                crs=crs,
                target_resolution=target_resolution,
                resampling=Resampling.nearest,
                nodata=np.nan,
            )
            f0 = f0.fillna(0.20)
            f1 = f1.fillna(0.05)
            f2 = f2.fillna(0.02)
            unc = self._qa_to_uncertainty(qa).fillna(0.08)

            f0_list.append(f0.expand_dims(band=[band_coord]))
            f1_list.append(f1.expand_dims(band=[band_coord]))
            f2_list.append(f2.expand_dims(band=[band_coord]))
            f0_unc_list.append(unc.expand_dims(band=[band_coord]))
            f1_unc_list.append((unc * 1.1).expand_dims(band=[band_coord]))
            f2_unc_list.append((unc * 1.1).expand_dims(band=[band_coord]))
            reflectance_unc_list.append(unc.expand_dims(band=[band_coord]))

        return BRDFKernelWeights(
            f0=xr.concat(f0_list, dim="band").transpose("band", "y", "x"),
            f1=xr.concat(f1_list, dim="band").transpose("band", "y", "x"),
            f2=xr.concat(f2_list, dim="band").transpose("band", "y", "x"),
            f0_unc=xr.concat(f0_unc_list, dim="band").transpose("band", "y", "x"),
            f1_unc=xr.concat(f1_unc_list, dim="band").transpose("band", "y", "x"),
            f2_unc=xr.concat(f2_unc_list, dim="band").transpose("band", "y", "x"),
            reflectance_unc=xr.concat(reflectance_unc_list, dim="band").transpose("band", "y", "x"),
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
    ) -> xr.DataArray:
        return merge_reprojected_tiles(
            arrays,
            bounds=bounds,
            crs=crs,
            resolution=target_resolution,
            resampling=resampling,
            nodata=nodata,
        )

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
        qa_values = qa.values.astype(np.float32)
        unc = np.full(qa_values.shape, np.nan, dtype=np.float32)
        valid = np.isfinite(qa_values) & (qa_values >= 0.0)
        unc = np.where(
            valid,
            _BEST_QA_REFLECTANCE_UNCERTAINTY * np.power(qa_values + 1.0, _QA_UNCERTAINTY_POWER),
            unc,
        )
        return xr.DataArray(unc, dims=qa.dims, coords=qa.coords)

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
    ) -> BRDFKernelWeights:
        base = self._default_weights(bounds, resolution, bands)
        time_coords = xr.IndexVariable("time", time_axis)

        def _repeat(data: xr.DataArray) -> xr.DataArray:
            repeated = np.repeat(data.values[np.newaxis, ...], len(time_axis), axis=0)
            return xr.DataArray(
                repeated,
                dims=["time", "band", "y", "x"],
                coords={
                    "time": time_coords,
                    "band": data.coords["band"],
                    "y": data.coords["y"],
                    "x": data.coords["x"],
                },
            )

        return BRDFKernelWeights(
            f0=_repeat(base.f0),
            f1=_repeat(base.f1),
            f2=_repeat(base.f2),
            f0_unc=_repeat(base.f0_unc),
            f1_unc=_repeat(base.f1_unc),
            f2_unc=_repeat(base.f2_unc),
            reflectance_unc=_repeat(base.reflectance_unc) if base.reflectance_unc is not None else None,
        )

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

        f0_list: list[xr.DataArray] = []
        f1_list: list[xr.DataArray] = []
        f2_list: list[xr.DataArray] = []
        f0_unc_list: list[xr.DataArray] = []
        f1_unc_list: list[xr.DataArray] = []
        f2_unc_list: list[xr.DataArray] = []
        reflectance_unc_list: list[xr.DataArray] = []

        for band_coord, product_band in requested:
            band_f0_days: list[xr.DataArray] = []
            band_f1_days: list[xr.DataArray] = []
            band_f2_days: list[xr.DataArray] = []
            band_f0_unc_days: list[xr.DataArray] = []
            band_f1_unc_days: list[xr.DataArray] = []
            band_f2_unc_days: list[xr.DataArray] = []
            band_reflectance_unc_days: list[xr.DataArray] = []

            for day in time_axis:
                day_paths = grouped_paths.get(day, [])
                if not day_paths:
                    empty = self._empty_spatial_array(bounds, target_resolution)
                    unc = self._empty_spatial_array(bounds, target_resolution)
                    band_f0_days.append(empty.expand_dims(time=[day]))
                    band_f1_days.append(empty.expand_dims(time=[day]))
                    band_f2_days.append(empty.expand_dims(time=[day]))
                    band_f0_unc_days.append(unc.expand_dims(time=[day]))
                    band_f1_unc_days.append(unc.expand_dims(time=[day]))
                    band_f2_unc_days.append(unc.expand_dims(time=[day]))
                    band_reflectance_unc_days.append(unc.expand_dims(time=[day]))
                    continue

                param_tiles: list[tuple[xr.DataArray, xr.DataArray, xr.DataArray]] = []
                qa_tiles: list[xr.DataArray] = []
                for path in day_paths:
                    params, qa = self._load_native_band_stack(path, product_band)
                    qa_tiles.append(qa)
                    param_tiles.append(params)

                f0 = self._merge_reprojected_tiles(
                    [params[0] for params in param_tiles],
                    bounds=bounds,
                    crs=crs,
                    target_resolution=target_resolution,
                    resampling=Resampling.bilinear,
                    nodata=np.nan,
                )
                f0 = self._coerce_to_target_grid(f0, bounds, target_resolution)
                f1 = self._merge_reprojected_tiles(
                    [params[1] for params in param_tiles],
                    bounds=bounds,
                    crs=crs,
                    target_resolution=target_resolution,
                    resampling=Resampling.bilinear,
                    nodata=np.nan,
                )
                f1 = self._coerce_to_target_grid(f1, bounds, target_resolution)
                f2 = self._merge_reprojected_tiles(
                    [params[2] for params in param_tiles],
                    bounds=bounds,
                    crs=crs,
                    target_resolution=target_resolution,
                    resampling=Resampling.bilinear,
                    nodata=np.nan,
                )
                f2 = self._coerce_to_target_grid(f2, bounds, target_resolution)
                qa = self._merge_reprojected_tiles(
                    qa_tiles,
                    bounds=bounds,
                    crs=crs,
                    target_resolution=target_resolution,
                    resampling=Resampling.nearest,
                    nodata=np.nan,
                )
                qa = self._coerce_to_target_grid(qa, bounds, target_resolution)
                unc = self._qa_to_uncertainty(qa).fillna(np.nan)
                band_f0_days.append(f0.expand_dims(time=[day]))
                band_f1_days.append(f1.expand_dims(time=[day]))
                band_f2_days.append(f2.expand_dims(time=[day]))
                band_f0_unc_days.append(unc.expand_dims(time=[day]))
                band_f1_unc_days.append((unc * 1.1).expand_dims(time=[day]))
                band_f2_unc_days.append((unc * 1.1).expand_dims(time=[day]))
                band_reflectance_unc_days.append(unc.expand_dims(time=[day]))

            f0_list.append(xr.concat(band_f0_days, dim="time").expand_dims(band=[band_coord]))
            f1_list.append(xr.concat(band_f1_days, dim="time").expand_dims(band=[band_coord]))
            f2_list.append(xr.concat(band_f2_days, dim="time").expand_dims(band=[band_coord]))
            f0_unc_list.append(xr.concat(band_f0_unc_days, dim="time").expand_dims(band=[band_coord]))
            f1_unc_list.append(xr.concat(band_f1_unc_days, dim="time").expand_dims(band=[band_coord]))
            f2_unc_list.append(xr.concat(band_f2_unc_days, dim="time").expand_dims(band=[band_coord]))
            reflectance_unc_list.append(xr.concat(band_reflectance_unc_days, dim="time").expand_dims(band=[band_coord]))

        temporal = BRDFKernelWeights(
            f0=xr.concat(f0_list, dim="band").transpose("time", "band", "y", "x"),
            f1=xr.concat(f1_list, dim="band").transpose("time", "band", "y", "x"),
            f2=xr.concat(f2_list, dim="band").transpose("time", "band", "y", "x"),
            f0_unc=xr.concat(f0_unc_list, dim="band").transpose("time", "band", "y", "x"),
            f1_unc=xr.concat(f1_unc_list, dim="band").transpose("time", "band", "y", "x"),
            f2_unc=xr.concat(f2_unc_list, dim="band").transpose("time", "band", "y", "x"),
            reflectance_unc=xr.concat(reflectance_unc_list, dim="band").transpose("time", "band", "y", "x"),
        )
        if np.isfinite(temporal.f0.values).any():
            return temporal
        return self._default_temporal_weights(
            bounds,
            target_resolution,
            [coord for coord, _band in requested],
            time_axis,
        )

    def _load_native_band_stack(
        self,
        path: str | Path,
        product_band: ProductBandDefinition,
    ) -> tuple[tuple[xr.DataArray, xr.DataArray, xr.DataArray], xr.DataArray]:
        raise NotImplementedError

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

    def _load_native_band_stack(
        self,
        path: str | Path,
        product_band: ProductBandDefinition,
    ) -> tuple[tuple[xr.DataArray, xr.DataArray, xr.DataArray], xr.DataArray]:
        params_raw, params_attrs = self._read_dataset(path, product_band.parameter_dataset)
        qa_raw, qa_attrs = self._read_dataset(path, product_band.qa_dataset or "")

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


class MCD43EarthAccessProvider(_StackParameterProvider):
    """MODIS MCD43 BRDF provider using real MCD43A1 kernel parameters."""

    product_key = "mcd43_brdf"
    _source_name = "MCD43"
    _rsrf_sensor_unit_id = "terra_modis"
    _product_bands = (
        ProductBandDefinition("Band1", 645.0, 50.0, "BRDF_Albedo_Parameters_Band1", "BRDF_Albedo_Band_Mandatory_Quality_Band1"),
        ProductBandDefinition("Band2", 858.5, 35.0, "BRDF_Albedo_Parameters_Band2", "BRDF_Albedo_Band_Mandatory_Quality_Band2"),
        ProductBandDefinition("Band3", 469.0, 20.0, "BRDF_Albedo_Parameters_Band3", "BRDF_Albedo_Band_Mandatory_Quality_Band3"),
        ProductBandDefinition("Band4", 555.0, 20.0, "BRDF_Albedo_Parameters_Band4", "BRDF_Albedo_Band_Mandatory_Quality_Band4"),
        ProductBandDefinition("Band5", 1240.0, 20.0, "BRDF_Albedo_Parameters_Band5", "BRDF_Albedo_Band_Mandatory_Quality_Band5"),
        ProductBandDefinition("Band6", 1640.0, 24.0, "BRDF_Albedo_Parameters_Band6", "BRDF_Albedo_Band_Mandatory_Quality_Band6"),
        ProductBandDefinition("Band7", 2130.0, 50.0, "BRDF_Albedo_Parameters_Band7", "BRDF_Albedo_Band_Mandatory_Quality_Band7"),
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
        ProductBandDefinition("Band1", 645.0, 50.0, "Kiso_Band1", "Status_QA"),
        ProductBandDefinition("Band2", 858.5, 35.0, "Kiso_Band2", "Status_QA"),
        ProductBandDefinition("Band3", 469.0, 20.0, "Kiso_Band3", "Status_QA"),
        ProductBandDefinition("Band4", 555.0, 20.0, "Kiso_Band4", "Status_QA"),
        ProductBandDefinition("Band5", 1240.0, 20.0, "Kiso_Band5", "Status_QA"),
        ProductBandDefinition("Band6", 1640.0, 24.0, "Kiso_Band6", "Status_QA"),
        ProductBandDefinition("Band7", 2130.0, 50.0, "Kiso_Band7", "Status_QA"),
        ProductBandDefinition("Band8", 412.0, 20.0, "Kiso_Band8", "Status_QA"),
    )
    _legacy_band_map = {index + 1: band.label for index, band in enumerate(_product_bands)}

    def _load_native_band_stack(
        self,
        path: str | Path,
        product_band: ProductBandDefinition,
    ) -> tuple[tuple[xr.DataArray, xr.DataArray, xr.DataArray], xr.DataArray]:
        label = product_band.label.replace("Band", "")
        f0_raw, f0_attrs = read_hdf4_dataset(path, f"Kiso_Band{label}")
        f1_raw, f1_attrs = read_hdf4_dataset(path, f"Kvol_Band{label}")
        f2_raw, f2_attrs = read_hdf4_dataset(path, f"Kgeo_Band{label}")
        qa_raw, qa_attrs = read_hdf4_dataset(path, "Status_QA")

        f0 = make_native_grid_dataarray(apply_scale_and_mask(f0_raw, f0_attrs), granule_path=path)
        f1 = make_native_grid_dataarray(apply_scale_and_mask(f1_raw, f1_attrs), granule_path=path)
        f2 = make_native_grid_dataarray(apply_scale_and_mask(f2_raw, f2_attrs), granule_path=path)
        qa = make_native_grid_dataarray(apply_scale_and_mask(qa_raw, qa_attrs), granule_path=path)
        return (f0, f1, f2), qa
