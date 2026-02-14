"""Earthaccess-backed BRDF providers (MCD43 + VNP43).

Phase M0 implementation: providers probe Earthdata availability through
Earthaccess and return stable default BRDF kernel weights until full granule
parsing is added.
"""

from __future__ import annotations

import logging
from datetime import datetime
from pathlib import Path
from typing import Sequence

import numpy as np
import xarray as xr

from siac.core.types import BRDFKernelWeights
from siac.io.earthaccess_catalog import EarthAccessCatalog
from siac.io.earthaccess_source import EarthAccessSource

logger = logging.getLogger(__name__)


class _EarthAccessBRDFProvider:
    """Shared logic for Earthaccess BRDF product wrappers."""

    product_key: str = ""
    _source_name: str = ""

    def __init__(
        self,
        cache_dir: str | Path | None = None,
        *,
        source: EarthAccessSource | None = None,
        catalog: EarthAccessCatalog | None = None,
        short_name: str | None = None,
        provider: str | None = None,
        probe_earthdata: bool = True,
    ) -> None:
        self.cache_dir = Path(cache_dir).expanduser() if cache_dir is not None else None
        self.source = source or EarthAccessSource(provider=provider)
        self.catalog = catalog or EarthAccessCatalog(source=self.source)
        self.short_name = short_name
        self.provider = provider
        self.probe_earthdata = probe_earthdata

    @property
    def source_name(self) -> str:
        return self._source_name

    def get_brdf_parameters(
        self,
        bounds: tuple[float, float, float, float],
        crs: str,
        obs_time: datetime,
        target_resolution: float,
        bands: Sequence[int],
        temporal_window: int = 16,
    ) -> BRDFKernelWeights:
        """Return BRDF kernel weights on the requested grid."""
        if self.probe_earthdata:
            self._probe_granules(bounds, crs, obs_time, temporal_window)

        # Phase M0 fallback values that are physically plausible and stable.
        return self._default_weights(bounds, target_resolution, list(bands))

    def _probe_granules(
        self,
        bounds: tuple[float, float, float, float],
        crs: str,
        obs_time: datetime,
        temporal_window: int,
    ) -> None:
        try:
            short_name = self.short_name or self.catalog.resolve_short_name(self.product_key)
            temporal = EarthAccessSource.temporal_window(obs_time, temporal_window)
            granules = self.source.search_granules(
                short_name=short_name,
                bounds=bounds,
                crs=crs,
                temporal=temporal,
                provider=self.provider,
                count=1,
            )
            if not granules:
                logger.warning(
                    "%s granule probe returned no results for AOI/time window", self._source_name
                )
        except Exception as exc:  # pragma: no cover - external/system dependent
            logger.warning("%s Earthaccess probe failed; using defaults (%s)", self._source_name, exc)

    @staticmethod
    def _grid(bounds: tuple[float, float, float, float], resolution: float) -> tuple[np.ndarray, np.ndarray]:
        xmin, ymin, xmax, ymax = bounds
        if resolution <= 0:
            raise ValueError(f"target_resolution must be > 0, got {resolution}")

        nx = max(1, int(np.ceil((xmax - xmin) / resolution)))
        ny = max(1, int(np.ceil((ymax - ymin) / resolution)))
        x = xmin + (np.arange(nx, dtype=np.float32) + 0.5) * resolution
        y = ymax - (np.arange(ny, dtype=np.float32) + 0.5) * resolution
        return y, x

    def _constant_band_array(
        self,
        bands: list[int],
        bounds: tuple[float, float, float, float],
        resolution: float,
        value: float,
    ) -> xr.DataArray:
        y, x = self._grid(bounds, resolution)
        arr = np.full((len(bands), y.size, x.size), value, dtype=np.float32)
        return xr.DataArray(
            arr,
            dims=["band", "y", "x"],
            coords={"band": np.array(bands, dtype=np.int16), "y": y, "x": x},
        )

    def _default_weights(
        self,
        bounds: tuple[float, float, float, float],
        resolution: float,
        bands: list[int],
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
        )


class MCD43EarthAccessProvider(_EarthAccessBRDFProvider):
    """MODIS MCD43 BRDF provider."""

    product_key = "mcd43_brdf"
    _source_name = "MCD43"


class VNP43EarthAccessProvider(_EarthAccessBRDFProvider):
    """VIIRS VNP43 BRDF provider."""

    product_key = "vnp43_brdf"
    _source_name = "VNP43"
