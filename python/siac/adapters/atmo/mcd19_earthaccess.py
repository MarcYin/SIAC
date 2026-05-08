"""MAIAC atmospheric prior providers using Earthaccess-backed granule parsing."""

from __future__ import annotations

import logging
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
    from datetime import datetime
    from pathlib import Path

    from siac.adapters.data.earthaccess_catalog import EarthAccessCatalog

logger = logging.getLogger(__name__)


class _NativeTileState(TypedDict):
    aot: xr.DataArray
    aot_unc: xr.DataArray
    tcwv: xr.DataArray | None


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
                    )
                except (OSError, ValueError, KeyError, RuntimeError) as exc:
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
    ) -> AtmosphericState:
        aot_tiles: list[xr.DataArray] = []
        aot_unc_tiles: list[xr.DataArray] = []
        tcwv_tiles: list[xr.DataArray] = []

        for path in paths:
            tile_state = self._load_native_tile(path)
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
        aot = aot.fillna(0.12)
        aot_unc = aot_unc.fillna(0.05)

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

    def _load_native_tile(self, path: str | Path) -> _NativeTileState:
        aod_raw, aod_attrs = self._read_dataset(path, self.aod_dataset)
        aod_unc_raw, aod_unc_attrs = self._read_dataset(path, self.aod_unc_dataset)
        qa_raw, qa_attrs = self._read_dataset(path, self.qa_dataset)

        aod = apply_scale_and_mask(aod_raw, aod_attrs)
        aod_unc = apply_scale_and_mask(aod_unc_raw, aod_unc_attrs)
        qa = apply_scale_and_mask(qa_raw, qa_attrs)

        valid = np.isfinite(aod) & np.isfinite(aod_unc) & np.isfinite(qa) & (qa > 0)
        aod = np.where(valid, aod, np.nan)
        aod_unc = np.where(valid, aod_unc, np.nan)

        tcwv_da: xr.DataArray | None = None
        if self.tcwv_dataset is not None:
            try:
                tcwv_raw, tcwv_attrs = self._read_dataset(path, self.tcwv_dataset)
                tcwv = apply_scale_and_mask(tcwv_raw, tcwv_attrs)
                tcwv = np.where(np.isfinite(tcwv), tcwv, np.nan)
                tcwv_da = make_native_grid_dataarray(
                    reduce_orbit_stack(tcwv),
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
            "aot": make_native_grid_dataarray(reduce_orbit_stack(aod), granule_path=path),
            "aot_unc": make_native_grid_dataarray(reduce_orbit_stack(aod_unc), granule_path=path),
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


class VNP19AODProvider(_EarthAccessMAIACAODProvider):
    """Atmospheric prior provider using real VIIRS MAIAC aerosol fields."""

    product_keys = ("vnp19_aod", "vj119_aod", "vj219_aod")
    _source_name = "VNP19"
    _read_dataset = staticmethod(read_hdf5_dataset)
    aod_dataset = "HDFEOS/GRIDS/grid750m/Data Fields/Optical_Depth_055"
    aod_unc_dataset = "HDFEOS/GRIDS/grid750m/Data Fields/AOD_Uncertainty"
    qa_dataset = "HDFEOS/GRIDS/grid750m/Data Fields/AOD_QA"
    tcwv_dataset = None
