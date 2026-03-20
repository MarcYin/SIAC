"""Sentinel-2 MSI sensor catalog entries."""

from siac.catalog.sensors._common import rsrf_band
from siac.domain.sensors import SensorConfig

SENTINEL2A_CONFIG = SensorConfig(
    sensor_id="MSI",
    satellite_id="S2A",
    bands=(
        rsrf_band("B01", 443.0, 20.0, 60.0, 0, sensor_unit_id="sentinel-2a_msi"),
        rsrf_band("B02", 490.0, 65.0, 10.0, 1, sensor_unit_id="sentinel-2a_msi"),
        rsrf_band("B03", 560.0, 35.0, 10.0, 2, sensor_unit_id="sentinel-2a_msi"),
        rsrf_band("B04", 665.0, 30.0, 10.0, 3, sensor_unit_id="sentinel-2a_msi"),
        rsrf_band("B05", 705.0, 15.0, 20.0, 4, sensor_unit_id="sentinel-2a_msi"),
        rsrf_band("B06", 740.0, 15.0, 20.0, 5, sensor_unit_id="sentinel-2a_msi"),
        rsrf_band("B07", 783.0, 20.0, 20.0, 6, sensor_unit_id="sentinel-2a_msi"),
        rsrf_band("B08", 842.0, 115.0, 10.0, 7, sensor_unit_id="sentinel-2a_msi"),
        rsrf_band("B8A", 865.0, 20.0, 20.0, 8, sensor_unit_id="sentinel-2a_msi"),
        rsrf_band("B09", 945.0, 20.0, 60.0, 9, sensor_unit_id="sentinel-2a_msi"),
        rsrf_band("B10", 1375.0, 30.0, 60.0, 10, sensor_unit_id="sentinel-2a_msi"),
        rsrf_band("B11", 1610.0, 90.0, 20.0, 11, sensor_unit_id="sentinel-2a_msi"),
        rsrf_band("B12", 2190.0, 180.0, 20.0, 12, sensor_unit_id="sentinel-2a_msi"),
    ),
    default_ref_scale=1.0 / 10000.0,
    default_ref_offset=0.0,
)

SENTINEL2B_CONFIG = SensorConfig(
    sensor_id="MSI",
    satellite_id="S2B",
    bands=(
        rsrf_band("B01", 442.0, 20.0, 60.0, 0, sensor_unit_id="sentinel-2b_msi"),
        rsrf_band("B02", 492.0, 65.0, 10.0, 1, sensor_unit_id="sentinel-2b_msi"),
        rsrf_band("B03", 559.0, 35.0, 10.0, 2, sensor_unit_id="sentinel-2b_msi"),
        rsrf_band("B04", 665.0, 30.0, 10.0, 3, sensor_unit_id="sentinel-2b_msi"),
        rsrf_band("B05", 704.0, 15.0, 20.0, 4, sensor_unit_id="sentinel-2b_msi"),
        rsrf_band("B06", 739.0, 15.0, 20.0, 5, sensor_unit_id="sentinel-2b_msi"),
        rsrf_band("B07", 780.0, 20.0, 20.0, 6, sensor_unit_id="sentinel-2b_msi"),
        rsrf_band("B08", 833.0, 115.0, 10.0, 7, sensor_unit_id="sentinel-2b_msi"),
        rsrf_band("B8A", 864.0, 20.0, 20.0, 8, sensor_unit_id="sentinel-2b_msi"),
        rsrf_band("B09", 943.0, 20.0, 60.0, 9, sensor_unit_id="sentinel-2b_msi"),
        rsrf_band("B10", 1377.0, 30.0, 60.0, 10, sensor_unit_id="sentinel-2b_msi"),
        rsrf_band("B11", 1610.0, 90.0, 20.0, 11, sensor_unit_id="sentinel-2b_msi"),
        rsrf_band("B12", 2186.0, 180.0, 20.0, 12, sensor_unit_id="sentinel-2b_msi"),
    ),
    default_ref_scale=1.0 / 10000.0,
    default_ref_offset=0.0,
)

SENTINEL2C_CONFIG = SensorConfig(
    sensor_id="MSI",
    satellite_id="S2C",
    bands=(
        rsrf_band("B01", 443.0, 20.0, 60.0, 0, sensor_unit_id="sentinel-2c_msi"),
        rsrf_band("B02", 490.0, 65.0, 10.0, 1, sensor_unit_id="sentinel-2c_msi"),
        rsrf_band("B03", 560.0, 35.0, 10.0, 2, sensor_unit_id="sentinel-2c_msi"),
        rsrf_band("B04", 665.0, 30.0, 10.0, 3, sensor_unit_id="sentinel-2c_msi"),
        rsrf_band("B05", 705.0, 15.0, 20.0, 4, sensor_unit_id="sentinel-2c_msi"),
        rsrf_band("B06", 740.0, 15.0, 20.0, 5, sensor_unit_id="sentinel-2c_msi"),
        rsrf_band("B07", 783.0, 20.0, 20.0, 6, sensor_unit_id="sentinel-2c_msi"),
        rsrf_band("B08", 842.0, 115.0, 10.0, 7, sensor_unit_id="sentinel-2c_msi"),
        rsrf_band("B8A", 865.0, 20.0, 20.0, 8, sensor_unit_id="sentinel-2c_msi"),
        rsrf_band("B09", 945.0, 20.0, 60.0, 9, sensor_unit_id="sentinel-2c_msi"),
        rsrf_band("B10", 1375.0, 30.0, 60.0, 10, sensor_unit_id="sentinel-2c_msi"),
        rsrf_band("B11", 1610.0, 90.0, 20.0, 11, sensor_unit_id="sentinel-2c_msi"),
        rsrf_band("B12", 2190.0, 180.0, 20.0, 12, sensor_unit_id="sentinel-2c_msi"),
    ),
    default_ref_scale=1.0 / 10000.0,
    default_ref_offset=0.0,
)


__all__ = ["SENTINEL2A_CONFIG", "SENTINEL2B_CONFIG", "SENTINEL2C_CONFIG"]
