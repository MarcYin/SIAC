"""Landsat SRF source provider."""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from pathlib import Path

    from siac.domain import SensorConfig

USGS_LANDSAT_SCV_URL = "https://landsat.usgs.gov/spectral-characteristics-viewer"


def load_landsat_sensor_config(
    sensor_id: str,
    satellite_id: str,
    *,
    cache_dir: str | Path | None = None,
    refresh: bool = False,
    local_path: str | Path | None = None,
) -> SensorConfig:
    """Load Landsat sensor SRFs from the official USGS source or a local file."""
    del cache_dir, refresh, local_path
    raise NotImplementedError(
        f"Landsat SRF loading is not implemented yet for {sensor_id}/{satellite_id}. "
        f"Official source: {USGS_LANDSAT_SCV_URL}"
    )
