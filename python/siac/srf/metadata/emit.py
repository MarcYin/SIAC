"""EMIT metadata-driven spectral characterization loader."""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from pathlib import Path

    from siac.core.types import SensorConfig

EMIT_ATBD_URL = "https://lpdaac.usgs.gov/documents/2147/EMIT_L2A-RFL_ATBD_V1.pdf"


def load_emit_sensor_config_from_metadata(metadata_path: str | Path) -> SensorConfig:
    """Load EMIT spectral characterization from granule metadata."""
    del metadata_path
    raise NotImplementedError(
        f"EMIT metadata-driven SRF loading is not implemented yet. Official source: {EMIT_ATBD_URL}"
    )
