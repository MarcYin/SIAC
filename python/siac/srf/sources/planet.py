"""PlanetScope / SuperDove SRF source provider."""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from pathlib import Path

    from siac.core.types import SensorConfig

PLANET_RSR_URL = "https://support.planet.com/hc/en-us/articles/4411132050451-How-Can-I-Access-Relative-Spectral-Responses-RSRs-"


def load_planet_sensor_config(
    sensor_id: str,
    satellite_id: str,
    *,
    cache_dir: str | Path | None = None,
    refresh: bool = False,
    local_path: str | Path | None = None,
) -> SensorConfig:
    """Load PlanetScope / SuperDove SRFs from the official Planet source or a local file."""
    del cache_dir, refresh, local_path
    raise NotImplementedError(
        f"Planet SRF loading is not implemented yet for {sensor_id}/{satellite_id}. "
        f"Official source: {PLANET_RSR_URL}"
    )
