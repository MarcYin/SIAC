"""MAIAC atmospheric prior providers using Earthaccess-backed granule parsing."""

from __future__ import annotations

import logging
import re
from datetime import datetime, timedelta, timezone
from functools import lru_cache
from typing import TYPE_CHECKING, TypedDict

import numpy as np
import xarray as xr
from rasterio.enums import Resampling

from siac.adapters.data.earthaccess_source import EarthAccessSource
from siac.adapters.earthdata import (
    build_earthaccess_runtime,
    constant_target_array,
    earthaccess_cache_dir,
    merge_reprojected_tiles,
    select_candidate_paths,
    target_grid_coords,
)
from siac.adapters.earthdata_common import (
    apply_scale_and_mask,
    make_native_grid_dataarray,
    read_hdf4_dataset,
    read_hdf5_dataset,
    reduce_orbit_stack,
)
from siac.runtime import AtmosphericState

if TYPE_CHECKING:
    from pathlib import Path

    from siac.adapters.data.earthaccess_catalog import EarthAccessCatalog

logger = logging.getLogger(__name__)

_MCD19_ORBIT_TIMESTAMP_RE = re.compile(
    r"(?P<year>\d{4})(?P<day_of_year>\d{3})(?P<hour>\d{2})(?P<minute>\d{2})[TA]"
)


@lru_cache(maxsize=256)
def _read_mcd19_orbit_times(path: str) -> tuple[datetime, ...]:
    """Read per-layer UTC acquisition times from MCD19A2 global metadata."""
    from osgeo import gdal

    gdal.UseExceptions()
    dataset = gdal.OpenEx(path, gdal.OF_RASTER | gdal.OF_READONLY)
    if dataset is None:
        return ()
    raw = str(dataset.GetMetadata().get("Orbit_time_stamp", ""))
    times: list[datetime] = []
    for match in _MCD19_ORBIT_TIMESTAMP_RE.finditer(raw):
        stamp = (
            f"{match.group('year')}{match.group('day_of_year')}"
            f"{match.group('hour')}{match.group('minute')}"
        )
        times.append(datetime.strptime(stamp, "%Y%j%H%M").replace(tzinfo=timezone.utc))
    return tuple(times)


def _nearest_valid_orbit_indices(
    valid: np.ndarray,
    orbit_times: tuple[datetime, ...],
    obs_time: datetime,
) -> np.ndarray | None:
    """Choose the nearest-in-time valid orbit independently at each pixel."""
    mask = np.asarray(valid, dtype=bool)
    if mask.ndim != 3 or len(orbit_times) != mask.shape[0]:
        return None
    target = obs_time
    if target.tzinfo is None:
        target = target.replace(tzinfo=timezone.utc)
    else:
        target = target.astimezone(timezone.utc)

    def distance_seconds(value: datetime) -> float:
        stamp = value
        if stamp.tzinfo is None:
            stamp = stamp.replace(tzinfo=timezone.utc)
        else:
            stamp = stamp.astimezone(timezone.utc)
        return abs((stamp - target).total_seconds())

    order = sorted(
        range(len(orbit_times)),
        key=lambda index: distance_seconds(orbit_times[index]),
    )
    selected = np.full(mask.shape[1:], -1, dtype=np.int16)
    for index in order:
        take = (selected < 0) & mask[index]
        selected[take] = index
    return selected


def _select_orbit_values(values: np.ndarray, indices: np.ndarray) -> np.ndarray:
    """Select one layer from an orbit stack using per-pixel layer indices."""
    arr = np.asarray(values, dtype=np.float32)
    if arr.ndim <= 2:
        return arr
    if arr.ndim != 3 or arr.shape[1:] != indices.shape:
        return reduce_orbit_stack(arr)
    out = np.full(indices.shape, np.nan, dtype=np.float32)
    valid = indices >= 0
    iy, ix = np.nonzero(valid)
    out[iy, ix] = arr[indices[iy, ix], iy, ix]
    return out


class _NativeTileState(TypedDict):
    aot: xr.DataArray
    aot_unc: xr.DataArray
    tcwv: xr.DataArray | None


def _maiac_best_quality_mask(qa_raw: np.ndarray) -> np.ndarray:
    """Best-quality MCD19A2/VNP19 ``AOD_QA`` mask from the RAW 16-bit field.

    The bit field must be decoded from the raw integer — scaling / valid-range
    masking (``apply_scale_and_mask``) would corrupt it. Keeps clear +
    best-quality + no-adjacency pixels, matching the harness ``maiac_qa.py``::

        bits 0-2  cloud mask     (001 = clear)
        bits 5-7  adjacency mask (000 = normal)
        bits 8-11 AOD QA level   (0000 = best quality)
    """
    qa_int = np.asarray(qa_raw, dtype=np.int64)
    mask = (
        ((qa_int & 0b111) == 1) & (((qa_int >> 5) & 0b111) == 0) & (((qa_int >> 8) & 0b1111) == 0)
    )
    mask_array: np.ndarray = np.asarray(mask, dtype=bool)
    return mask_array


class _EarthAccessMAIACAODProvider:
    """Shared MAIAC AOD prior logic for MODIS and VIIRS products."""

    product_keys: tuple[str, ...] = ()
    _source_name: str = ""
    _read_dataset = staticmethod(read_hdf4_dataset)
    aod_dataset: str = "Optical_Depth_055"
    aod_unc_dataset: str = "AOD_Uncertainty"
    qa_dataset: str = "AOD_QA"
    tcwv_dataset: str | None = "Column_WV"

    def __init__(
        self,
        cache_dir: str | Path | None = None,
        *,
        source: EarthAccessSource | None = None,
        catalog: EarthAccessCatalog | None = None,
        short_name: str | None = None,
        provider: str | None = None,
        probe_earthdata: bool = True,
        temporal_window_days: int = 2,
        max_granules: int = 8,
        best_quality_qa: bool = False,
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
        self.temporal_window_days = max(0, int(temporal_window_days))
        self.max_granules = max(1, int(max_granules))
        self.best_quality_qa = bool(best_quality_qa)

    @property
    def source_name(self) -> str:
        return self._source_name

    def get_prior(
        self,
        bounds: tuple[float, float, float, float],
        crs: str,
        obs_time: datetime,
        resolution: float,
    ) -> AtmosphericState:
        if resolution <= 0:
            raise ValueError(f"resolution must be > 0, got {resolution}")

        if self.probe_earthdata:
            paths, used_short_name = self._download_granules(bounds, crs, obs_time)
            if paths:
                try:
                    return self._load_from_granules(
                        paths,
                        bounds=bounds,
                        crs=crs,
                        resolution=resolution,
                        short_name=used_short_name,
                        obs_time=obs_time,
                    )
                except (OSError, ValueError, KeyError, RuntimeError):
                    # Narrowed (REVIEW.md §2.1, §3.3 mcd19_earthaccess.py:108).
                    # OSError covers HDF/NetCDF I/O; ValueError/KeyError cover
                    # malformed metadata; RuntimeError covers gdal/rasterio
                    # wrapping. exc_info=True so the cause is visible.
                    logger.warning(
                        "%s granule parsing failed; using defaults",
                        self._source_name,
                        exc_info=True,
                    )

        return self._default_prior(bounds, resolution)

    def _download_granules(
        self,
        bounds: tuple[float, float, float, float],
        crs: str,
        obs_time: datetime,
    ) -> tuple[list[Path], str | None]:
        temporal = EarthAccessSource.temporal_window(obs_time, self.temporal_window_days)
        candidate_short_names = (
            (self.short_name,)
            if self.short_name is not None
            else tuple(self.catalog.resolve_short_name(key) for key in self.product_keys)
        )

        for short_name in candidate_short_names:
            granules = self.source.search_granules(
                short_name=short_name,
                bounds=bounds,
                crs=crs,
                temporal=temporal,
                provider=self.provider,
                count=self.max_granules,
            )
            if not granules:
                continue

            dest = earthaccess_cache_dir(self.cache_dir, short_name)
            downloaded = self.source.download_granules(granules, dest)
            selected = self._select_candidate_paths(downloaded, obs_time, bounds, crs)
            if selected:
                return selected, short_name

        logger.warning(
            "No %s granules found via Earthaccess for requested AOI/time", self._source_name
        )
        return [], None

    @staticmethod
    def _select_candidate_paths(
        paths: list[Path],
        obs_time: datetime,
        bounds: tuple[float, float, float, float],
        crs: str,
    ) -> list[Path]:
        return select_candidate_paths(
            paths,
            obs_time=obs_time,
            bounds=bounds,
            crs=crs,
        )

    def _load_from_granules(
        self,
        paths: list[Path],
        *,
        bounds: tuple[float, float, float, float],
        crs: str,
        resolution: float,
        short_name: str | None,
        obs_time: datetime | None = None,
    ) -> AtmosphericState:
        aot_tiles: list[xr.DataArray] = []
        aot_unc_tiles: list[xr.DataArray] = []
        tcwv_tiles: list[xr.DataArray] = []

        for path in paths:
            tile_state = (
                self._load_native_tile(path, obs_time=obs_time)
                if obs_time is not None
                else self._load_native_tile(path)
            )
            aot_tiles.append(tile_state["aot"])
            aot_unc_tiles.append(tile_state["aot_unc"])
            if tile_state["tcwv"] is not None:
                tcwv_tiles.append(tile_state["tcwv"])

        aot = self._merge_tiles(aot_tiles, bounds=bounds, crs=crs, resolution=resolution)
        aot_unc = self._merge_tiles(
            aot_unc_tiles,
            bounds=bounds,
            crs=crs,
            resolution=resolution,
        )
        finite_aot = np.asarray(aot.values, dtype=np.float64)
        finite_aot = finite_aot[np.isfinite(finite_aot)]
        finite_unc = np.asarray(aot_unc.values, dtype=np.float64)
        finite_unc = finite_unc[np.isfinite(finite_unc)]
        if finite_aot.size == 0:
            raise ValueError(
                f"{self._source_name} has no QA-valid AOD after reprojection to the requested AOI"
            )
        aot_fill = float(np.median(finite_aot))
        unc_fill = max(float(np.median(finite_unc)) if finite_unc.size else 0.10, 0.05)
        aot = aot.fillna(aot_fill)
        aot_unc = aot_unc.fillna(unc_fill)

        if tcwv_tiles:
            tcwv = self._merge_tiles(tcwv_tiles, bounds=bounds, crs=crs, resolution=resolution)
            tcwv = tcwv.fillna(1.5)
            tcwv_unc = xr.full_like(tcwv, 0.3)
        else:
            tcwv = xr.full_like(aot, 1.5)
            tcwv_unc = xr.full_like(aot, 0.3)
            if short_name is not None:
                logger.info(
                    "%s %s does not provide TCWV; using default prior values.",
                    self._source_name,
                    short_name,
                )

        tco3 = xr.full_like(aot, 0.30)
        tco3_unc = xr.full_like(aot, 0.03)
        elevation = xr.zeros_like(aot)

        return AtmosphericState(
            aot=aot.astype(np.float32),
            tcwv=tcwv.astype(np.float32),
            tco3=tco3.astype(np.float32),
            aot_unc=aot_unc.astype(np.float32),
            tcwv_unc=tcwv_unc.astype(np.float32),
            tco3_unc=tco3_unc.astype(np.float32),
            elevation=elevation.astype(np.float32),
        )

    def _load_native_tile(
        self,
        path: str | Path,
        *,
        obs_time: datetime | None = None,
    ) -> _NativeTileState:
        aod_raw, aod_attrs = self._read_dataset(path, self.aod_dataset)
        aod_unc_raw, aod_unc_attrs = self._read_dataset(path, self.aod_unc_dataset)
        qa_raw, qa_attrs = self._read_dataset(path, self.qa_dataset)

        aod = apply_scale_and_mask(aod_raw, aod_attrs)
        aod_unc = apply_scale_and_mask(aod_unc_raw, aod_unc_attrs)

        finite = np.isfinite(aod) & np.isfinite(aod_unc)
        qa = apply_scale_and_mask(qa_raw, qa_attrs)
        loose = finite & np.isfinite(qa) & (qa > 0)

        if self.best_quality_qa:
            # Best-quality AOD_QA bit decode (clear + best + no-adjacency), matching
            # the harness maiac_qa.py. The loose ``qa > 0`` below keeps lower-quality
            # retrievals and reads systematically higher at clean-coastal/polar sites.
            valid = finite & _maiac_best_quality_mask(qa_raw)
            # Over a small AOI on a hazy day the strict mask can reject every
            # pixel. Falling through to no AOD at all is worse than a
            # lower-quality retrieval: the caller then substitutes a default,
            # which reads far too low exactly where aerosol is thick and the
            # prior matters most. Degrade to the loose mask instead.
            if not valid.any() and loose.any():
                logger.info(
                    "%s: no best-quality AOD pixels in this granule; "
                    "falling back to the loose QA mask (%d pixels).",
                    self._source_name,
                    int(loose.sum()),
                )
                valid = loose
        else:
            valid = loose
        aod = np.where(valid, aod, np.nan)
        aod_unc = np.where(valid, aod_unc, np.nan)

        orbit_indices: np.ndarray | None = None
        if obs_time is not None and self._source_name == "MCD19":
            try:
                orbit_indices = _nearest_valid_orbit_indices(
                    valid,
                    _read_mcd19_orbit_times(str(path)),
                    obs_time,
                )
            except (OSError, RuntimeError, ValueError):
                logger.warning(
                    "Could not read MCD19 orbit timestamps from %s; using orbit mean",
                    path,
                    exc_info=True,
                )

        if orbit_indices is None:
            aod_2d = reduce_orbit_stack(aod)
            aod_unc_2d = reduce_orbit_stack(aod_unc)
        else:
            aod_2d = _select_orbit_values(aod, orbit_indices)
            aod_unc_2d = _select_orbit_values(aod_unc, orbit_indices)

        tcwv_da: xr.DataArray | None = None
        if self.tcwv_dataset is not None:
            try:
                tcwv_raw, tcwv_attrs = self._read_dataset(path, self.tcwv_dataset)
                tcwv = apply_scale_and_mask(tcwv_raw, tcwv_attrs)
                tcwv = np.where(np.isfinite(tcwv), tcwv, np.nan)
                tcwv_2d = (
                    reduce_orbit_stack(tcwv)
                    if orbit_indices is None
                    else _select_orbit_values(tcwv, orbit_indices)
                )
                tcwv_da = make_native_grid_dataarray(
                    tcwv_2d,
                    granule_path=path,
                )
            except (OSError, KeyError, ValueError, RuntimeError):
                # Narrowed (REVIEW.md §2.1): OSError covers HDF I/O,
                # KeyError missing TCWV dataset, ValueError malformed
                # scale/offset, RuntimeError gdal wrappers. Log so an
                # operator notices missing TCWV.
                logger.warning(
                    "%s TCWV dataset %r unavailable in granule %s; using AOD-only prior",
                    self._source_name,
                    self.tcwv_dataset,
                    path,
                    exc_info=True,
                )
                tcwv_da = None

        return {
            "aot": make_native_grid_dataarray(aod_2d, granule_path=path),
            "aot_unc": make_native_grid_dataarray(aod_unc_2d, granule_path=path),
            "tcwv": tcwv_da,
        }

    @staticmethod
    def _merge_tiles(
        arrays: list[xr.DataArray],
        *,
        bounds: tuple[float, float, float, float],
        crs: str,
        resolution: float,
    ) -> xr.DataArray:
        return merge_reprojected_tiles(
            arrays,
            bounds=bounds,
            crs=crs,
            resolution=resolution,
            resampling=Resampling.nearest,
            nodata=np.nan,
        )

    @staticmethod
    def _grid(
        bounds: tuple[float, float, float, float], resolution: float
    ) -> tuple[np.ndarray, np.ndarray]:
        return target_grid_coords(bounds, resolution, resolution_name="resolution")

    def _constant_array(
        self,
        bounds: tuple[float, float, float, float],
        resolution: float,
        value: float,
    ) -> xr.DataArray:
        return constant_target_array(bounds, resolution, value, resolution_name="resolution")

    def _default_prior(
        self,
        bounds: tuple[float, float, float, float],
        resolution: float,
    ) -> AtmosphericState:
        aot = self._constant_array(bounds, resolution, 0.12)
        tcwv = self._constant_array(bounds, resolution, 1.5)
        tco3 = self._constant_array(bounds, resolution, 0.30)

        return AtmosphericState(
            aot=aot,
            tcwv=tcwv,
            tco3=tco3,
            aot_unc=xr.full_like(aot, 0.05),
            tcwv_unc=xr.full_like(tcwv, 0.3),
            tco3_unc=xr.full_like(tco3, 0.03),
            elevation=xr.zeros_like(aot),
        )


class MCD19AODProvider(_EarthAccessMAIACAODProvider):
    """Atmospheric prior provider using real MCD19A2 MAIAC fields."""

    product_keys = ("mcd19_aod",)
    _source_name = "MCD19"


class CachedMCD19AODProvider(MCD19AODProvider):
    """Read same-day MCD19A2 fields from an already staged local cache.

    This is deliberately network-free.  The normal Earthaccess provider can
    silently fall back to a climatological field when a remote query or a
    reprojection fails; callers that need to distinguish a real MAIAC field
    from that fallback can use :meth:`get_cached_prior` and provide their own
    fallback policy.
    """

    def __init__(
        self,
        cache_dir: str | Path,
        *,
        temporal_window_days: int = 0,
        max_granules: int = 40,
        best_quality_qa: bool = True,
    ) -> None:
        super().__init__(
            cache_dir=cache_dir,
            probe_earthdata=False,
            temporal_window_days=temporal_window_days,
            max_granules=max_granules,
            best_quality_qa=best_quality_qa,
        )
        if self.cache_dir is None:
            raise ValueError("CachedMCD19AODProvider requires a cache_dir.")

    def _cached_paths(self, obs_time: datetime) -> list[Path]:
        """Return local MCD19A2 granules in the configured temporal window."""
        timestamp = (
            obs_time.replace(tzinfo=timezone.utc)
            if obs_time.tzinfo is None
            else obs_time.astimezone(timezone.utc)
        )
        paths: list[Path] = []
        for offset in range(-self.temporal_window_days, self.temporal_window_days + 1):
            day = timestamp.date() + timedelta(days=offset)
            stamp = day.strftime("%Y%j")
            paths.extend(self.cache_dir.glob(f"MCD19A2.A{stamp}.*"))
        return sorted(set(paths), key=lambda path: path.name)

    def get_cached_prior(
        self,
        bounds: tuple[float, float, float, float],
        crs: str,
        obs_time: datetime,
        resolution: float,
    ) -> AtmosphericState | None:
        """Return a QA-valid cached MAIAC field, or ``None`` when unavailable."""
        paths = self._cached_paths(obs_time)
        selected = self._select_candidate_paths(paths, obs_time, bounds, crs)
        if not selected:
            logger.info(
                "No cached %s granules cover the requested AOI/time.", self._source_name
            )
            return None
        try:
            return self._load_from_granules(
                selected[: self.max_granules],
                bounds=bounds,
                crs=crs,
                resolution=resolution,
                short_name="MCD19A2",
                obs_time=obs_time,
            )
        except ValueError as exc:
            # A same-day granule can legitimately have no QA-valid AOD at a
            # small AOI.  That is an expected coverage condition for callers
            # which supply a staged MAIAC fallback, not a parsing failure.
            if "no QA-valid AOD" in str(exc):
                logger.info(
                    "Cached %s has no QA-valid AOD over the requested AOI.",
                    self._source_name,
                )
                return None
            logger.warning(
                "Cached %s granule parsing failed; no spatial MAIAC field available",
                self._source_name,
                exc_info=True,
            )
            return None
        except (OSError, KeyError, RuntimeError):
            logger.warning(
                "Cached %s granule parsing failed; no spatial MAIAC field available",
                self._source_name,
                exc_info=True,
            )
            return None

    def get_prior(
        self,
        bounds: tuple[float, float, float, float],
        crs: str,
        obs_time: datetime,
        resolution: float,
    ) -> AtmosphericState:
        state = self.get_cached_prior(bounds, crs, obs_time, resolution)
        return state if state is not None else self._default_prior(bounds, resolution)


class VNP19AODProvider(_EarthAccessMAIACAODProvider):
    """Atmospheric prior provider using real VIIRS MAIAC aerosol fields."""

    product_keys = ("vnp19_aod", "vj119_aod", "vj219_aod")
    _source_name = "VNP19"
    _read_dataset = staticmethod(read_hdf5_dataset)
    aod_dataset = "HDFEOS/GRIDS/grid750m/Data Fields/Optical_Depth_055"
    aod_unc_dataset = "HDFEOS/GRIDS/grid750m/Data Fields/AOD_Uncertainty"
    qa_dataset = "HDFEOS/GRIDS/grid750m/Data Fields/AOD_QA"
    tcwv_dataset = None
