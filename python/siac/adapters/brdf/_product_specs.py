"""Static Earthaccess BRDF product band definitions."""

from __future__ import annotations

from siac.adapters.earthdata_common import ProductBandDefinition

MCD43_PRODUCT_BANDS = (
    ProductBandDefinition(
        "Band1",
        645.0,
        50.0,
        "BRDF_Albedo_Parameters_Band1",
        "BRDF_Albedo_Band_Mandatory_Quality_Band1",
        rsrf_band_id="B1",
    ),
    ProductBandDefinition(
        "Band2",
        858.5,
        35.0,
        "BRDF_Albedo_Parameters_Band2",
        "BRDF_Albedo_Band_Mandatory_Quality_Band2",
        rsrf_band_id="B2",
    ),
    ProductBandDefinition(
        "Band3",
        469.0,
        20.0,
        "BRDF_Albedo_Parameters_Band3",
        "BRDF_Albedo_Band_Mandatory_Quality_Band3",
        rsrf_band_id="B3",
    ),
    ProductBandDefinition(
        "Band4",
        555.0,
        20.0,
        "BRDF_Albedo_Parameters_Band4",
        "BRDF_Albedo_Band_Mandatory_Quality_Band4",
        rsrf_band_id="B4",
    ),
    ProductBandDefinition(
        "Band5",
        1240.0,
        20.0,
        "BRDF_Albedo_Parameters_Band5",
        "BRDF_Albedo_Band_Mandatory_Quality_Band5",
        rsrf_band_id="B5",
    ),
    ProductBandDefinition(
        "Band6",
        1640.0,
        24.0,
        "BRDF_Albedo_Parameters_Band6",
        "BRDF_Albedo_Band_Mandatory_Quality_Band6",
        rsrf_band_id="B6",
    ),
    ProductBandDefinition(
        "Band7",
        2130.0,
        50.0,
        "BRDF_Albedo_Parameters_Band7",
        "BRDF_Albedo_Band_Mandatory_Quality_Band7",
        rsrf_band_id="B7",
    ),
)

VNP43_PRODUCT_BANDS = (
    ProductBandDefinition(
        "M1",
        412.0,
        20.0,
        "HDFEOS/GRIDS/VIIRS_Grid_BRDF/Data Fields/BRDF_Albedo_Parameters_M1",
        "HDFEOS/GRIDS/VIIRS_Grid_BRDF/Data Fields/BRDF_Albedo_Band_Mandatory_Quality_M1",
    ),
    ProductBandDefinition(
        "M2",
        445.0,
        18.0,
        "HDFEOS/GRIDS/VIIRS_Grid_BRDF/Data Fields/BRDF_Albedo_Parameters_M2",
        "HDFEOS/GRIDS/VIIRS_Grid_BRDF/Data Fields/BRDF_Albedo_Band_Mandatory_Quality_M2",
    ),
    ProductBandDefinition(
        "M3",
        488.0,
        20.0,
        "HDFEOS/GRIDS/VIIRS_Grid_BRDF/Data Fields/BRDF_Albedo_Parameters_M3",
        "HDFEOS/GRIDS/VIIRS_Grid_BRDF/Data Fields/BRDF_Albedo_Band_Mandatory_Quality_M3",
    ),
    ProductBandDefinition(
        "M4",
        555.0,
        20.0,
        "HDFEOS/GRIDS/VIIRS_Grid_BRDF/Data Fields/BRDF_Albedo_Parameters_M4",
        "HDFEOS/GRIDS/VIIRS_Grid_BRDF/Data Fields/BRDF_Albedo_Band_Mandatory_Quality_M4",
    ),
    ProductBandDefinition(
        "M5",
        672.0,
        20.0,
        "HDFEOS/GRIDS/VIIRS_Grid_BRDF/Data Fields/BRDF_Albedo_Parameters_M5",
        "HDFEOS/GRIDS/VIIRS_Grid_BRDF/Data Fields/BRDF_Albedo_Band_Mandatory_Quality_M5",
    ),
    ProductBandDefinition(
        "M7",
        865.0,
        39.0,
        "HDFEOS/GRIDS/VIIRS_Grid_BRDF/Data Fields/BRDF_Albedo_Parameters_M7",
        "HDFEOS/GRIDS/VIIRS_Grid_BRDF/Data Fields/BRDF_Albedo_Band_Mandatory_Quality_M7",
    ),
    ProductBandDefinition(
        "M8",
        1240.0,
        20.0,
        "HDFEOS/GRIDS/VIIRS_Grid_BRDF/Data Fields/BRDF_Albedo_Parameters_M8",
        "HDFEOS/GRIDS/VIIRS_Grid_BRDF/Data Fields/BRDF_Albedo_Band_Mandatory_Quality_M8",
    ),
    ProductBandDefinition(
        "M10",
        1610.0,
        60.0,
        "HDFEOS/GRIDS/VIIRS_Grid_BRDF/Data Fields/BRDF_Albedo_Parameters_M10",
        "HDFEOS/GRIDS/VIIRS_Grid_BRDF/Data Fields/BRDF_Albedo_Band_Mandatory_Quality_M10",
    ),
    ProductBandDefinition(
        "M11",
        2250.0,
        50.0,
        "HDFEOS/GRIDS/VIIRS_Grid_BRDF/Data Fields/BRDF_Albedo_Parameters_M11",
        "HDFEOS/GRIDS/VIIRS_Grid_BRDF/Data Fields/BRDF_Albedo_Band_Mandatory_Quality_M11",
    ),
)

MCD19_PRODUCT_BANDS = (
    ProductBandDefinition("Band1", 645.0, 50.0, "Kiso_Band1", "Status_QA", rsrf_band_id="B1"),
    ProductBandDefinition("Band2", 858.5, 35.0, "Kiso_Band2", "Status_QA", rsrf_band_id="B2"),
    ProductBandDefinition("Band3", 469.0, 20.0, "Kiso_Band3", "Status_QA", rsrf_band_id="B3"),
    ProductBandDefinition("Band4", 555.0, 20.0, "Kiso_Band4", "Status_QA", rsrf_band_id="B4"),
    ProductBandDefinition("Band5", 1240.0, 20.0, "Kiso_Band5", "Status_QA", rsrf_band_id="B5"),
    ProductBandDefinition("Band6", 1640.0, 24.0, "Kiso_Band6", "Status_QA", rsrf_band_id="B6"),
    ProductBandDefinition("Band7", 2130.0, 50.0, "Kiso_Band7", "Status_QA", rsrf_band_id="B7"),
    ProductBandDefinition("Band8", 412.0, 20.0, "Kiso_Band8", "Status_QA", rsrf_band_id="B8"),
)
