"""EnMAP metadata-driven spectral characterization loader."""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from pathlib import Path

    from siac.core.types import SensorConfig

ENMAP_SPEC_URL = "https://www.enmap.org/data/doc/EN-PCV-ICD-2009-2_HSI_Product_Specification_Level1_Level2.pdf"


def load_enmap_sensor_config_from_metadata(metadata_path: str | Path) -> SensorConfig:
    """Load EnMAP spectral characterization from product metadata."""
    del metadata_path
    raise NotImplementedError(
        f"EnMAP metadata-driven SRF loading is not implemented yet. Official source: {ENMAP_SPEC_URL}"
    )
