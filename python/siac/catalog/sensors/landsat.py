"""Landsat OLI sensor catalog entries."""

from siac.catalog.sensors._common import rsrf_band
from siac.domain.sensors import SensorConfig

# Landsat OLI/OLI-2 aerosol retrieval picks the visible coastal-aerosol +
# blue anchors. Setting this on the catalog entry (rather than relying on
# the 400-900 nm generic fallback) preserves the prior 400-520 nm behaviour
# specifically for Landsat. See REVIEW.md §3.2.
_OLI_AEROSOL_BANDS = ("B1", "B2")

LANDSAT8_OLI_CONFIG = SensorConfig(
    sensor_id="OLI",
    satellite_id="L8",
    bands=(
        rsrf_band("B1", 443.0, 16.0, 30.0, 0, sensor_unit_id="landsat-8_oli"),
        rsrf_band("B2", 482.0, 60.0, 30.0, 1, sensor_unit_id="landsat-8_oli"),
        rsrf_band("B3", 561.5, 57.0, 30.0, 2, sensor_unit_id="landsat-8_oli"),
        rsrf_band("B4", 654.5, 37.0, 30.0, 3, sensor_unit_id="landsat-8_oli"),
        rsrf_band("B5", 865.0, 28.0, 30.0, 4, sensor_unit_id="landsat-8_oli"),
        rsrf_band("B6", 1608.5, 85.0, 30.0, 5, sensor_unit_id="landsat-8_oli"),
        rsrf_band("B7", 2200.5, 187.0, 30.0, 6, sensor_unit_id="landsat-8_oli"),
    ),
    default_ref_scale=2.75e-5,
    default_ref_offset=-0.2,
    aerosol_solver_band_names=_OLI_AEROSOL_BANDS,
)

LANDSAT9_OLI2_CONFIG = SensorConfig(
    sensor_id="OLI",
    satellite_id="L9",
    bands=(
        rsrf_band("B1", 443.0, 16.0, 30.0, 0, sensor_unit_id="landsat-9_oli2"),
        rsrf_band("B2", 482.0, 60.0, 30.0, 1, sensor_unit_id="landsat-9_oli2"),
        rsrf_band("B3", 561.5, 57.0, 30.0, 2, sensor_unit_id="landsat-9_oli2"),
        rsrf_band("B4", 654.5, 37.0, 30.0, 3, sensor_unit_id="landsat-9_oli2"),
        rsrf_band("B5", 865.0, 28.0, 30.0, 4, sensor_unit_id="landsat-9_oli2"),
        rsrf_band("B6", 1608.5, 85.0, 30.0, 5, sensor_unit_id="landsat-9_oli2"),
        rsrf_band("B7", 2200.5, 187.0, 30.0, 6, sensor_unit_id="landsat-9_oli2"),
    ),
    default_ref_scale=2.75e-5,
    default_ref_offset=-0.2,
    aerosol_solver_band_names=_OLI_AEROSOL_BANDS,
)
