"""MAIAC atmospheric prior providers using Earthaccess-backed granule parsing."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
import xarray as xr
from rasterio.enums import Resampling

from siac.core.types import AtmosphericState
from siac.io.earthaccess_catalog import EarthAccessCatalog
from siac.io.earthaccess_source import EarthAccessSource
from siac.io.reprojection import transform_bounds
from siac.priors.earthdata_common import (
    MODLAND_SINUSOIDAL_CRS,
    apply_scale_and_mask,
    make_native_grid_dataarray,
    modland_tile_bounds,
    parse_granule_date,
    parse_tile_indices,
    read_hdf4_dataset,
    read_hdf5_dataset,
    reduce_orbit_stack,
    reproject_native_to_target,
)

if TYPE_CHECKING:
    from datetime import datetime

logger = logging.getLogger(__name__)


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
        self.cache_dir = Path(cache_dir).expanduser() if cache_dir is not None else None
        self.source = source or EarthAccessSource(provider=provider)
        self.catalog = catalog or EarthAccessCatalog(source=self.source)
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
                except Exception as exc:  # pragma: no cover - external/system dependent
                    logger.warning(
                        "%s granule parsing failed; using defaults (%s)",
                        self._source_name,
                        exc,
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
            (self.short_name,) if self.short_name is not None else tuple(
                self.catalog.resolve_short_name(key) for key in self.product_keys
            )
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

            dest = self.cache_dir or Path.home() / ".cache" / "siac" / "earthdata" / short_name
            downloaded = self.source.download_granules(granules, dest)
            selected = self._select_candidate_paths(downloaded, obs_time, bounds, crs)
            if selected:
                return selected, short_name

        logger.warning("No %s granules found via Earthaccess for requested AOI/time", self._source_name)
        return [], None

    @staticmethod
    def _select_candidate_paths(
        paths: list[Path],
        obs_time: datetime,
        bounds: tuple[float, float, float, float],
        crs: str,
    ) -> list[Path]:
        if not paths:
            return []

        target_bounds_native = transform_bounds(bounds, crs, MODLAND_SINUSOIDAL_CRS)
        selected: list[tuple[tuple[int, int], float, Path]] = []
        for path in paths:
            try:
                tile = parse_tile_indices(path)
                delta = abs((parse_granule_date(path) - obs_time).total_seconds())
                tile_bounds = modland_tile_bounds(*tile)
                intersects = not (
                    tile_bounds[2] <= target_bounds_native[0]
                    or tile_bounds[0] >= target_bounds_native[2]
                    or tile_bounds[3] <= target_bounds_native[1]
                    or tile_bounds[1] >= target_bounds_native[3]
                )
            except Exception:
                return paths

            if intersects:
                selected.append((tile, delta, path))

        if not selected:
            return []
        return [item[2] for item in sorted(selected, key=lambda value: (value[0], value[1]))]

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
                logger.info("%s %s does not provide TCWV; using default prior values.", self._source_name, short_name)

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

    def _load_native_tile(self, path: str | Path) -> dict[str, xr.DataArray | None]:
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
            except Exception:
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
        reprojected = [
            reproject_native_to_target(
                arr,
                target_bounds=bounds,
                target_crs=crs,
                target_resolution=resolution,
                resampling=Resampling.bilinear,
                nodata=np.nan,
            )
            for arr in arrays
        ]
        merged = reprojected[0]
        for arr in reprojected[1:]:
            merged = merged.combine_first(arr)
        return merged

    @staticmethod
    def _grid(bounds: tuple[float, float, float, float], resolution: float) -> tuple[np.ndarray, np.ndarray]:
        xmin, ymin, xmax, ymax = bounds
        if resolution <= 0:
            raise ValueError(f"resolution must be > 0, got {resolution}")

        nx = max(1, int(np.ceil((xmax - xmin) / resolution)))
        ny = max(1, int(np.ceil((ymax - ymin) / resolution)))
        x = xmin + (np.arange(nx, dtype=np.float32) + 0.5) * resolution
        y = ymax - (np.arange(ny, dtype=np.float32) + 0.5) * resolution
        return y, x

    def _constant_array(
        self,
        bounds: tuple[float, float, float, float],
        resolution: float,
        value: float,
    ) -> xr.DataArray:
        y, x = self._grid(bounds, resolution)
        arr = np.full((y.size, x.size), value, dtype=np.float32)
        return xr.DataArray(arr, dims=["y", "x"], coords={"y": y, "x": x})

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
