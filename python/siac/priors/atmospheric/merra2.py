"""MERRA-2 atmospheric prior provider with Earthaccess discovery hooks.

Phase M0 implementation: this provider validates queryability through
Earthaccess and currently falls back to climatological prior fields until full
MERRA-2 granule parsing is implemented.
"""

from __future__ import annotations

import logging
from datetime import datetime
from pathlib import Path

import numpy as np
import xarray as xr

from siac.core.types import AtmosphericState
from siac.io.earthaccess_catalog import EarthAccessCatalog
from siac.io.earthaccess_source import EarthAccessSource

logger = logging.getLogger(__name__)


class MERRA2Provider:
    """Atmospheric prior provider using Earthaccess product discovery."""

    def __init__(
        self,
        cache_dir: str | Path | None = None,
        *,
        source: EarthAccessSource | None = None,
        catalog: EarthAccessCatalog | None = None,
        short_name: str | None = None,
        provider: str | None = None,
        probe_earthdata: bool = True,
        temporal_window_days: int = 1,
    ) -> None:
        self.cache_dir = Path(cache_dir).expanduser() if cache_dir is not None else None
        self.source = source or EarthAccessSource(provider=provider)
        self.catalog = catalog or EarthAccessCatalog(source=self.source)
        self.short_name = short_name
        self.provider = provider
        self.probe_earthdata = probe_earthdata
        self.temporal_window_days = max(0, int(temporal_window_days))

    @property
    def source_name(self) -> str:
        return "MERRA-2"

    def get_prior(
        self,
        bounds: tuple[float, float, float, float],
        crs: str,
        obs_time: datetime,
        resolution: float,
    ) -> AtmosphericState:
        """Return atmospheric prior on the requested spatial grid."""
        if self.probe_earthdata:
            self._probe_granules(bounds, crs, obs_time)

        # Phase M0 fallback: return stable climatological priors while
        # granule variable parsing is being wired in later milestones.
        return self._default_prior(bounds, resolution)

    def _probe_granules(
        self,
        bounds: tuple[float, float, float, float],
        crs: str,
        obs_time: datetime,
    ) -> None:
        """Best-effort Earthdata availability probe for monitoring/debugging."""
        try:
            short_name = self.short_name or self.catalog.resolve_short_name("merra2_atmo")
            temporal = EarthAccessSource.temporal_window(obs_time, self.temporal_window_days)
            granules = self.source.search_granules(
                short_name=short_name,
                bounds=bounds,
                crs=crs,
                temporal=temporal,
                provider=self.provider,
                count=1,
            )
            if not granules:
                logger.warning("No MERRA-2 granules found via Earthaccess for requested AOI/time")
        except Exception as exc:  # pragma: no cover - external/system dependent
            logger.warning("MERRA-2 Earthaccess probe failed; using defaults (%s)", exc)

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
        # Conservative defaults similar to CAMS fallback branch.
        aot = self._constant_array(bounds, resolution, 0.15)
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
