"""Per-day MAIAC AOD extraction for the bestpixel surface-prior gate.

The surface-driven retrieval's bestpixel prior wants to composite only the
clearest (lowest-aerosol) acquisition days per (year, month) window. The
existing :class:`~siac.adapters.atmo.mcd19_earthaccess._EarthAccessMAIACAODProvider`
reduces a whole temporal window into a single gridded prior; here we instead
need one scalar AOD value *per acquisition day* so bestpixel can drop the
high-AOD days before compositing.

:class:`MAIACDayAODProvider` searches + downloads MCD19A2 (MAIAC) granules for
each requested month, parses each granule's acquisition day from its filename,
reprojects that day's tiles onto the AOI, and reduces to one area-median
``Optical_Depth_055`` value per day. earthaccess is reached lazily (only when
:meth:`day_aod_map` is called); construction does no network I/O.
"""

from __future__ import annotations

import logging
from collections import defaultdict
from typing import TYPE_CHECKING

import numpy as np
from rasterio.enums import Resampling

from siac.adapters.earthdata import (
    build_earthaccess_runtime,
    earthaccess_cache_dir,
    merge_reprojected_tiles,
)
from siac.adapters.earthdata_common import (
    apply_scale_and_mask,
    make_native_grid_dataarray,
    parse_granule_date,
    read_hdf4_dataset,
    reduce_orbit_stack,
)

if TYPE_CHECKING:
    from collections.abc import Iterable, Sequence
    from pathlib import Path

    import xarray as xr

    from siac.adapters.data.earthaccess_catalog import EarthAccessCatalog
    from siac.adapters.data.earthaccess_source import EarthAccessSource

logger = logging.getLogger(__name__)

#: MODIS/MAIAC product the day-AOD gate sources ``Optical_Depth_055`` from.
_MCD19_PRODUCT_KEY = "mcd19_aod"
#: AOI-clip resolution (metres) for the per-day area-median reduction. Coarse
#: on purpose — the output is a single scalar per day, so the grid only sets
#: the AOI clip granularity, not retrieval resolution.
_DEFAULT_GATE_RESOLUTION_M = 1000.0


class MAIACDayAODProvider:
    """Resolve per-acquisition-day MAIAC ``Optical_Depth_055`` over an AOI.

    Mirrors the granule search/download path of the MAIAC atmospheric-prior
    provider but exposes per-day scalars instead of a window-reduced grid.
    """

    aod_dataset: str = "Optical_Depth_055"

    def __init__(
        self,
        cache_dir: str | Path | None = None,
        *,
        source: EarthAccessSource | None = None,
        catalog: EarthAccessCatalog | None = None,
        provider: str | None = None,
        gate_resolution_m: float = _DEFAULT_GATE_RESOLUTION_M,
        max_granules_per_month: int = 64,
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
        self.provider = provider
        self.gate_resolution_m = float(gate_resolution_m)
        self.max_granules_per_month = max(1, int(max_granules_per_month))

    def day_aod_map(
        self,
        bounds: tuple[float, float, float, float],
        crs: str,
        periods: Iterable[tuple[int, int]],
    ) -> dict[str, float]:
        """Return ``{"YYYY-MM-DD": aod}`` for all candidate days in ``periods``.

        ``periods`` is an iterable of ``(year, month)`` windows. Each window is
        searched independently; a window that yields no granules (or whose
        granule parsing fails) is simply omitted so the gate keeps every scene
        day for it. The result is resilient: a per-window failure is logged and
        skipped rather than aborting the whole map.
        """
        short_name = self.catalog.resolve_short_name(_MCD19_PRODUCT_KEY)
        day_aod: dict[str, float] = {}
        for year, month in dict.fromkeys((int(y), int(m)) for y, m in periods):
            try:
                day_aod.update(
                    self._window_day_aod(
                        bounds=bounds,
                        crs=crs,
                        year=year,
                        month=month,
                        short_name=short_name,
                    )
                )
            except Exception:  # noqa: BLE001 - per-window resilience (network/HDF)
                logger.warning(
                    "MAIAC day-AOD gate: window %04d-%02d failed; keeping all its days",
                    year,
                    month,
                    exc_info=True,
                )
        return day_aod

    def _window_day_aod(
        self,
        *,
        bounds: tuple[float, float, float, float],
        crs: str,
        year: int,
        month: int,
        short_name: str,
    ) -> dict[str, float]:
        temporal = _month_temporal_range(year, month)
        granules = self.source.search_granules(
            short_name=short_name,
            bounds=bounds,
            crs=crs,
            temporal=temporal,
            provider=self.provider,
            count=self.max_granules_per_month,
        )
        if not granules:
            logger.info("MAIAC day-AOD gate: no granules for %04d-%02d", year, month)
            return {}

        dest = earthaccess_cache_dir(self.cache_dir, short_name)
        paths = self.source.download_granules(granules, dest)
        return self._reduce_paths_to_day_aod(paths, bounds=bounds, crs=crs)

    def _reduce_paths_to_day_aod(
        self,
        paths: Sequence[Path],
        *,
        bounds: tuple[float, float, float, float],
        crs: str,
    ) -> dict[str, float]:
        by_day: dict[str, list[xr.DataArray]] = defaultdict(list)
        for path in paths:
            try:
                day = parse_granule_date(path).date().isoformat()
                tile = self._load_aod_tile(path)
            except Exception:  # noqa: BLE001 - skip an unreadable granule
                logger.warning("MAIAC day-AOD gate: skipping granule %s", path, exc_info=True)
                continue
            by_day[day].append(tile)

        day_aod: dict[str, float] = {}
        for day, tiles in by_day.items():
            merged = merge_reprojected_tiles(
                tiles,
                bounds=bounds,
                crs=crs,
                resolution=self.gate_resolution_m,
                resampling=Resampling.nearest,
                nodata=np.nan,
            )
            values = np.asarray(merged.values, dtype=np.float64)
            finite = values[np.isfinite(values)]
            if finite.size:
                day_aod[day] = float(np.median(finite))
        return day_aod

    def _load_aod_tile(self, path: str | Path) -> xr.DataArray:
        aod_raw, aod_attrs = read_hdf4_dataset(path, self.aod_dataset)
        aod = apply_scale_and_mask(aod_raw, aod_attrs)
        aod = np.where(np.isfinite(aod), aod, np.nan)
        return make_native_grid_dataarray(reduce_orbit_stack(aod), granule_path=path)


def _month_temporal_range(year: int, month: int) -> tuple[str, str]:
    """Build a ``(start, end)`` ISO range spanning the whole calendar month."""
    import calendar

    last_day = calendar.monthrange(year, month)[1]
    start = f"{year:04d}-{month:02d}-01T00:00:00Z"
    end = f"{year:04d}-{month:02d}-{last_day:02d}T23:59:59Z"
    return start, end
