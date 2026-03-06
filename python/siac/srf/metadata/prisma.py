"""PRISMA metadata-driven spectral characterization loader."""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from pathlib import Path

    from siac.core.types import SensorConfig

PRISMA_MISSION_URL = "https://www.asi.it/en/earth-science/prisma/"


def load_prisma_sensor_config_from_metadata(metadata_path: str | Path) -> SensorConfig:
    """Load PRISMA spectral characterization from product metadata."""
    del metadata_path
    raise NotImplementedError(
        f"PRISMA metadata-driven SRF loading is not implemented yet. Official source: {PRISMA_MISSION_URL}"
    )
