/*************************************************************************

NAME: landsat8.c

PURPOSE: Defines the attributes for Landsat 8

**************************************************************************/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "ias_logging.h"
#include "ias_satellite_attributes.h"  /* Function prototypes and additional
                                          required headers */
#include "local_defines.h"

#define OLI_DETECTORS_PER_SCA_NORMAL 494
#define OLI_DETECTORS_PER_SCA_PAN 988
#define OLI_DETECTORS_PER_SCA_VRP 12
#define OLI_DETECTORS_PER_SCA_SWIR_BLIND 104 
#define OLI_DETECTORS_PER_SCA_CIRRUS_BLIND 103
#define OLI_DETECTORS_PER_SCA_BLIND_VRP 65
#define OLI_DETECTORS_PER_SCA_PAN_VRP 24
#define TIRS_DETECTORS_PER_SCA_NORMAL 640
#define TIRS_DETECTORS_PER_SCA_BLIND 640
#define OLI_PIXEL_RESOLUTION_NORMAL 30
#define OLI_PIXEL_RESOLUTION_PAN 15
#define TIRS_PIXEL_RESOLUTION 100
#define OLI_NUMBER_OF_SCA 14
#define TIRS_NUMBER_OF_SCA 3
#define L8_TOTAL_NUMBER_OF_NORMAL_BANDS 11
#define L8_TOTAL_NUMBER_OF_BANDS 30
/* The QUANTIZE_CAL_MIN and QUANTIZE_CAL_MAX defines are used by 
   the radiance rescaling to adjustthe output range of the data. */
#define QUANTIZE_CAL_MIN 1
#define QUANTIZE_CAL_MAX 65535

static IAS_SATELLITE_ATTRIBUTES l8_sat_attribs;
static IAS_BAND_ATTRIBUTES l8_band_attribs[L8_TOTAL_NUMBER_OF_BANDS];
static IAS_SENSOR_ID sensor_id[2];

/*************************************************************************

NAME: ias_sat_attr_initialize_Landsat8

PURPOSE: Initializes an IAS_SATELLITE_ATTRIBUTES structure for
         Landsat 8.

RETURNS: Constant pointer to the created IAS_SATELLITE_ATTRIBUTES
         structure or NULL if an error occurs

NOTES:
    The band numbers for Landsat 8 are defined as:
        1: OLI_Coastal_Aerosol
        2: OLI_Blue
        3: OLI_Green
        4: OLI_Red
        5: OLI_NIR
        6: OLI_SWIR1
        7: OLI_SWIR2
        8: OLI_PAN
        9: OLI_CIRRUS
        10: TIRS_THERMAL1
        11: TIRS_THERMAL2
        12: OLI_BLIND_SWIR1
        13: OLI_BLIND_SWIR2
        14: OLI_BLIND_CIRRUS
        15: TIRS_BLIND
        16: TIRS_THERMAL1_SECONDARY
        17: TIRS_THERMAL2_SECONDARY
        18: TIRS_BLIND_SECONDARY
        19: OLI_VRP_Coastal_Aerosol
        20: OLI_VRP_Blue
        21: OLI_VRP_Green
        22: OLI_VRP_Red
        23: OLI_VRP_NIR
        24: OLI_VRP_SWIR1
        25: OLI_VRP_SWIR2
        26: OLI_VRP_PAN
        27: OLI_VRP_CIRRUS
        28: OLI_VRP_BLIND_SWIR1
        29: OLI_VRP_BLIND_SWIR2
        30: OLI_VRP_BLIND_CIRRUS

**************************************************************************/
IAS_SATELLITE_ATTRIBUTES *ias_sat_attr_initialize_landsat8()
{
    int band_index;

    /* Allocate sensor_ids pointer to static structure object */
    l8_sat_attribs.sensor_ids = &sensor_id[0]; 

    /* Allocate l8_band_attributes pointer to static structure object */
    l8_sat_attribs.band_attributes = &l8_band_attribs[0];

    /* Initialization the Satellite Attributes structure */
    l8_sat_attribs.satellite_id = IAS_L8;
    l8_sat_attribs.sensors = 2;
    l8_sat_attribs.bands = L8_TOTAL_NUMBER_OF_NORMAL_BANDS;
    l8_sat_attribs.total_bands = L8_TOTAL_NUMBER_OF_BANDS;

    strncpy(l8_sat_attribs.satellite_name, IAS_SATELLITE_NAME_L8,
        sizeof(l8_sat_attribs.satellite_name));
    l8_sat_attribs.satellite_name[sizeof(l8_sat_attribs.satellite_name) - 1]
        = '\0';
    strncpy(l8_sat_attribs.sensor_name, IAS_SENSOR_NAME_L8,
        sizeof(l8_sat_attribs.sensor_name));
    l8_sat_attribs.sensor_name[sizeof(l8_sat_attribs.sensor_name) - 1] = '\0';

    sensor_id[0] = IAS_OLI;
    sensor_id[1] = IAS_TIRS;

    for ( band_index = 0; band_index < L8_TOTAL_NUMBER_OF_BANDS; band_index++)
    {
        l8_band_attribs[band_index].band_number = band_index + 1;
        l8_band_attribs[band_index].band_index = band_index;  
        l8_band_attribs[band_index].qcal_min = QUANTIZE_CAL_MIN;
        l8_band_attribs[band_index].qcal_max = QUANTIZE_CAL_MAX;
        l8_band_attribs[band_index].can_saturate = FALSE;
    }

    /* Set up the OLI Coastal Aerosol band (band number 1) */
    band_index = 0;
    strcpy(l8_band_attribs[band_index].band_name, "OLI_Coastal_Aerosol");
    l8_band_attribs[band_index].sensor_id = IAS_OLI;
    l8_band_attribs[band_index].sensor_type = IAS_PUSHBROOM_SENSOR;
    l8_band_attribs[band_index].band_type = IAS_COASTAL_AEROSOL_BAND;
    l8_band_attribs[band_index].band_classification = IAS_NORMAL_BAND;
    l8_band_attribs[band_index].spectral_type = IAS_SPECTRAL_VNIR;
    l8_band_attribs[band_index].normal_band_number = 1; 
    l8_band_attribs[band_index].vrp_band_number = 19; 
    l8_band_attribs[band_index].blind_band_number = 0; 
    l8_band_attribs[band_index].secondary_band_number = 0;
    l8_band_attribs[band_index].scas = OLI_NUMBER_OF_SCA;
    l8_band_attribs[band_index].detectors_per_sca = 
                       OLI_DETECTORS_PER_SCA_NORMAL;
    l8_band_attribs[band_index].lines_per_frame = 1;
    l8_band_attribs[band_index].pixel_resolution = OLI_PIXEL_RESOLUTION_NORMAL;
    l8_band_attribs[band_index].wavelength_nm_range[0]=433;
    l8_band_attribs[band_index].wavelength_nm_range[1]=453;
    l8_band_attribs[band_index].can_saturate = TRUE;

    /* Set up the OLI blue band (band number 2) */
    band_index++;
    strcpy(l8_band_attribs[band_index].band_name, "OLI_Blue");
    l8_band_attribs[band_index].sensor_id = IAS_OLI;
    l8_band_attribs[band_index].sensor_type = IAS_PUSHBROOM_SENSOR;
    l8_band_attribs[band_index].band_type = IAS_BLUE_BAND;
    l8_band_attribs[band_index].band_classification = IAS_NORMAL_BAND;
    l8_band_attribs[band_index].spectral_type = IAS_SPECTRAL_VNIR;
    l8_band_attribs[band_index].normal_band_number = 2; 
    l8_band_attribs[band_index].vrp_band_number = 20; 
    l8_band_attribs[band_index].blind_band_number = 0; 
    l8_band_attribs[band_index].secondary_band_number = 0;
    l8_band_attribs[band_index].scas = OLI_NUMBER_OF_SCA;
    l8_band_attribs[band_index].detectors_per_sca = 
                       OLI_DETECTORS_PER_SCA_NORMAL;
    l8_band_attribs[band_index].lines_per_frame = 1;
    l8_band_attribs[band_index].pixel_resolution = OLI_PIXEL_RESOLUTION_NORMAL;
    l8_band_attribs[band_index].wavelength_nm_range[0]=450;
    l8_band_attribs[band_index].wavelength_nm_range[1]=515;
    l8_band_attribs[band_index].can_saturate = TRUE;

    /* Set up the OLI Green band (band number 3) */
    band_index++;
    strcpy(l8_band_attribs[band_index].band_name, "OLI_Green");
    l8_band_attribs[band_index].sensor_id = IAS_OLI;
    l8_band_attribs[band_index].sensor_type = IAS_PUSHBROOM_SENSOR;
    l8_band_attribs[band_index].band_type = IAS_GREEN_BAND;
    l8_band_attribs[band_index].band_classification = IAS_NORMAL_BAND;
    l8_band_attribs[band_index].spectral_type = IAS_SPECTRAL_VNIR;
    l8_band_attribs[band_index].normal_band_number = 3; 
    l8_band_attribs[band_index].vrp_band_number = 21; 
    l8_band_attribs[band_index].blind_band_number = 0; 
    l8_band_attribs[band_index].secondary_band_number = 0;
    l8_band_attribs[band_index].scas = OLI_NUMBER_OF_SCA;
    l8_band_attribs[band_index].detectors_per_sca = 
                            OLI_DETECTORS_PER_SCA_NORMAL;
    l8_band_attribs[band_index].lines_per_frame = 1;
    l8_band_attribs[band_index].pixel_resolution = OLI_PIXEL_RESOLUTION_NORMAL;
    l8_band_attribs[band_index].wavelength_nm_range[0]=525;
    l8_band_attribs[band_index].wavelength_nm_range[1]=600;
    l8_band_attribs[band_index].can_saturate = TRUE;

    /* Set up the OLI Red band (band number 4) */
    band_index++;
    strcpy(l8_band_attribs[band_index].band_name, "OLI_Red");
    l8_band_attribs[band_index].sensor_id = IAS_OLI;
    l8_band_attribs[band_index].sensor_type = IAS_PUSHBROOM_SENSOR;
    l8_band_attribs[band_index].band_type = IAS_RED_BAND;
    l8_band_attribs[band_index].band_classification = IAS_NORMAL_BAND;
    l8_band_attribs[band_index].spectral_type = IAS_SPECTRAL_VNIR;
    l8_band_attribs[band_index].normal_band_number = 4; 
    l8_band_attribs[band_index].vrp_band_number = 22; 
    l8_band_attribs[band_index].blind_band_number = 0; 
    l8_band_attribs[band_index].secondary_band_number = 0;
    l8_band_attribs[band_index].scas = OLI_NUMBER_OF_SCA;
    l8_band_attribs[band_index].detectors_per_sca = 
                            OLI_DETECTORS_PER_SCA_NORMAL;
    l8_band_attribs[band_index].lines_per_frame = 1;
    l8_band_attribs[band_index].pixel_resolution = OLI_PIXEL_RESOLUTION_NORMAL;
    l8_band_attribs[band_index].wavelength_nm_range[0]=630;
    l8_band_attribs[band_index].wavelength_nm_range[1]=680;
    l8_band_attribs[band_index].can_saturate = TRUE;

    /* Set up the OLI Near Infrared (NIR) band (band number 5) */
    band_index++;
    strcpy(l8_band_attribs[band_index].band_name, "OLI_NIR");
    l8_band_attribs[band_index].sensor_id = IAS_OLI;
    l8_band_attribs[band_index].sensor_type = IAS_PUSHBROOM_SENSOR;
    l8_band_attribs[band_index].band_type = IAS_NIR_BAND;
    l8_band_attribs[band_index].band_classification = IAS_NORMAL_BAND;
    l8_band_attribs[band_index].spectral_type = IAS_SPECTRAL_VNIR;
    l8_band_attribs[band_index].normal_band_number = 5; 
    l8_band_attribs[band_index].vrp_band_number = 23; 
    l8_band_attribs[band_index].blind_band_number = 0; 
    l8_band_attribs[band_index].secondary_band_number = 0;
    l8_band_attribs[band_index].scas = OLI_NUMBER_OF_SCA;
    l8_band_attribs[band_index].detectors_per_sca = 
                            OLI_DETECTORS_PER_SCA_NORMAL;
    l8_band_attribs[band_index].lines_per_frame = 1;
    l8_band_attribs[band_index].pixel_resolution = OLI_PIXEL_RESOLUTION_NORMAL;
    l8_band_attribs[band_index].wavelength_nm_range[0]=845;
    l8_band_attribs[band_index].wavelength_nm_range[1]=885;
    l8_band_attribs[band_index].can_saturate = TRUE;

    /* Set up the first OLI Shortwave Infrared (SWIR) band (band number 6) */
    band_index++;
    strcpy(l8_band_attribs[band_index].band_name, "OLI_SWIR1");
    l8_band_attribs[band_index].sensor_id = IAS_OLI;
    l8_band_attribs[band_index].sensor_type = IAS_PUSHBROOM_SENSOR;
    l8_band_attribs[band_index].band_type = IAS_SWIR1_BAND;
    l8_band_attribs[band_index].band_classification = IAS_NORMAL_BAND;
    l8_band_attribs[band_index].spectral_type = IAS_SPECTRAL_SWIR;
    l8_band_attribs[band_index].normal_band_number = 6; 
    l8_band_attribs[band_index].vrp_band_number = 24; 
    l8_band_attribs[band_index].blind_band_number = 12; 
    l8_band_attribs[band_index].secondary_band_number = 0;
    l8_band_attribs[band_index].scas = OLI_NUMBER_OF_SCA;
    l8_band_attribs[band_index].detectors_per_sca = 
                            OLI_DETECTORS_PER_SCA_NORMAL;
    l8_band_attribs[band_index].lines_per_frame = 1;
    l8_band_attribs[band_index].pixel_resolution = OLI_PIXEL_RESOLUTION_NORMAL;
    l8_band_attribs[band_index].wavelength_nm_range[0]=1560;
    l8_band_attribs[band_index].wavelength_nm_range[1]=1660;
    l8_band_attribs[band_index].can_saturate = TRUE;

    /* Set up the second OLI Shortwave Infrared (SWIR) band (band number 7) */
    band_index++;
    strcpy(l8_band_attribs[band_index].band_name, "OLI_SWIR2");
    l8_band_attribs[band_index].sensor_id = IAS_OLI;
    l8_band_attribs[band_index].sensor_type = IAS_PUSHBROOM_SENSOR;
    l8_band_attribs[band_index].band_type = IAS_SWIR2_BAND;
    l8_band_attribs[band_index].band_classification = IAS_NORMAL_BAND;
    l8_band_attribs[band_index].spectral_type = IAS_SPECTRAL_SWIR;
    l8_band_attribs[band_index].normal_band_number = 7; 
    l8_band_attribs[band_index].vrp_band_number = 25; 
    l8_band_attribs[band_index].blind_band_number = 13; 
    l8_band_attribs[band_index].secondary_band_number = 0;
    l8_band_attribs[band_index].scas = OLI_NUMBER_OF_SCA;
    l8_band_attribs[band_index].detectors_per_sca = 
                            OLI_DETECTORS_PER_SCA_NORMAL;
    l8_band_attribs[band_index].lines_per_frame = 1;
    l8_band_attribs[band_index].pixel_resolution = OLI_PIXEL_RESOLUTION_NORMAL;
    l8_band_attribs[band_index].wavelength_nm_range[0]=2100;
    l8_band_attribs[band_index].wavelength_nm_range[1]=2300;
    l8_band_attribs[band_index].can_saturate = TRUE;

    /* Set up the OLI Panchromatic band (band number 8).
       Note: Band 8 can saturate, but it will never be reported
             so its 'can_saturate' attribute it set to false */
    band_index++;
    strcpy(l8_band_attribs[band_index].band_name, "OLI_PAN");
    l8_band_attribs[band_index].sensor_id = IAS_OLI;
    l8_band_attribs[band_index].sensor_type = IAS_PUSHBROOM_SENSOR;
    l8_band_attribs[band_index].band_type = IAS_PAN_BAND;
    l8_band_attribs[band_index].band_classification = IAS_NORMAL_BAND;
    l8_band_attribs[band_index].spectral_type = IAS_SPECTRAL_PAN;
    l8_band_attribs[band_index].normal_band_number = 8; 
    l8_band_attribs[band_index].vrp_band_number = 26; 
    l8_band_attribs[band_index].blind_band_number = 0; 
    l8_band_attribs[band_index].secondary_band_number = 0;
    l8_band_attribs[band_index].scas = OLI_NUMBER_OF_SCA;
    l8_band_attribs[band_index].detectors_per_sca = OLI_DETECTORS_PER_SCA_PAN;
    l8_band_attribs[band_index].lines_per_frame = 2;
    l8_band_attribs[band_index].pixel_resolution = OLI_PIXEL_RESOLUTION_PAN;
    l8_band_attribs[band_index].wavelength_nm_range[0]=500;
    l8_band_attribs[band_index].wavelength_nm_range[1]=680;
    l8_band_attribs[band_index].can_saturate = FALSE;

    /* Set up the OLI Cirrus Cloud band (band number 9) */
    band_index++;
    strcpy(l8_band_attribs[band_index].band_name, "OLI_CIRRUS");
    l8_band_attribs[band_index].sensor_id = IAS_OLI;
    l8_band_attribs[band_index].sensor_type = IAS_PUSHBROOM_SENSOR;
    l8_band_attribs[band_index].band_type = IAS_CIRRUS_BAND;
    l8_band_attribs[band_index].band_classification = IAS_NORMAL_BAND;
    l8_band_attribs[band_index].spectral_type = IAS_SPECTRAL_SWIR;
    l8_band_attribs[band_index].normal_band_number = 9; 
    l8_band_attribs[band_index].vrp_band_number = 27; 
    l8_band_attribs[band_index].blind_band_number = 14; 
    l8_band_attribs[band_index].secondary_band_number = 0;
    l8_band_attribs[band_index].scas = OLI_NUMBER_OF_SCA;
    l8_band_attribs[band_index].detectors_per_sca = 
                            OLI_DETECTORS_PER_SCA_NORMAL;
    l8_band_attribs[band_index].lines_per_frame = 1;
    l8_band_attribs[band_index].pixel_resolution = OLI_PIXEL_RESOLUTION_NORMAL;
    l8_band_attribs[band_index].wavelength_nm_range[0]=1360;
    l8_band_attribs[band_index].wavelength_nm_range[1]=1390;
    l8_band_attribs[band_index].can_saturate = TRUE;

    /* Set up the first TIRS thermal band (band number 10) */
    band_index++;
    strcpy(l8_band_attribs[band_index].band_name, "TIRS_THERMAL1");
    l8_band_attribs[band_index].sensor_id = IAS_TIRS;
    l8_band_attribs[band_index].sensor_type = IAS_PUSHBROOM_SENSOR;
    l8_band_attribs[band_index].band_type = IAS_THERMAL1_BAND;
    l8_band_attribs[band_index].band_classification = IAS_NORMAL_BAND;
    l8_band_attribs[band_index].spectral_type = IAS_SPECTRAL_THERMAL;
    l8_band_attribs[band_index].normal_band_number = 10; 
    l8_band_attribs[band_index].vrp_band_number = 0; 
    l8_band_attribs[band_index].blind_band_number = 15; 
    l8_band_attribs[band_index].secondary_band_number = 16;
    l8_band_attribs[band_index].scas = TIRS_NUMBER_OF_SCA;
    l8_band_attribs[band_index].detectors_per_sca = 
                            TIRS_DETECTORS_PER_SCA_NORMAL;
    l8_band_attribs[band_index].lines_per_frame = 1;
    l8_band_attribs[band_index].pixel_resolution = TIRS_PIXEL_RESOLUTION;
    l8_band_attribs[band_index].wavelength_nm_range[0]=1030;
    l8_band_attribs[band_index].wavelength_nm_range[1]=1130;

    /* Set up the second TIRS thermal band (band number 11) */
    band_index++;
    strcpy(l8_band_attribs[band_index].band_name, "TIRS_THERMAL2");
    l8_band_attribs[band_index].sensor_id = IAS_TIRS;
    l8_band_attribs[band_index].sensor_type = IAS_PUSHBROOM_SENSOR;
    l8_band_attribs[band_index].band_type = IAS_THERMAL2_BAND;
    l8_band_attribs[band_index].band_classification = IAS_NORMAL_BAND;
    l8_band_attribs[band_index].spectral_type = IAS_SPECTRAL_THERMAL;
    l8_band_attribs[band_index].normal_band_number = 11; 
    l8_band_attribs[band_index].vrp_band_number = 0; 
    l8_band_attribs[band_index].blind_band_number = 15;
    l8_band_attribs[band_index].secondary_band_number = 17; 
    l8_band_attribs[band_index].scas = TIRS_NUMBER_OF_SCA;
    l8_band_attribs[band_index].detectors_per_sca = 
                            TIRS_DETECTORS_PER_SCA_NORMAL;
    l8_band_attribs[band_index].lines_per_frame = 1;
    l8_band_attribs[band_index].pixel_resolution = TIRS_PIXEL_RESOLUTION;
    l8_band_attribs[band_index].wavelength_nm_range[0]=1150;
    l8_band_attribs[band_index].wavelength_nm_range[1]=1250;

    /* Set up the OLI SWIR1 blind band (band number 12) */
    band_index++;
    strcpy(l8_band_attribs[band_index].band_name, "OLI_BLIND_SWIR1");
    l8_band_attribs[band_index].sensor_id = IAS_OLI;
    l8_band_attribs[band_index].sensor_type = IAS_PUSHBROOM_SENSOR;
    l8_band_attribs[band_index].band_type = IAS_SWIR1_BAND;
    l8_band_attribs[band_index].band_classification = IAS_BLIND_BAND;
    l8_band_attribs[band_index].spectral_type = IAS_SPECTRAL_SWIR;
    l8_band_attribs[band_index].normal_band_number = 6; 
    l8_band_attribs[band_index].vrp_band_number = 28; 
    l8_band_attribs[band_index].blind_band_number = 12; 
    l8_band_attribs[band_index].secondary_band_number = 0;
    l8_band_attribs[band_index].scas = OLI_NUMBER_OF_SCA;
    l8_band_attribs[band_index].detectors_per_sca = 
                            OLI_DETECTORS_PER_SCA_SWIR_BLIND;
    l8_band_attribs[band_index].lines_per_frame = 1;
    l8_band_attribs[band_index].pixel_resolution = OLI_PIXEL_RESOLUTION_NORMAL;
    l8_band_attribs[band_index].wavelength_nm_range[0]=0;
    l8_band_attribs[band_index].wavelength_nm_range[1]=0;

    /* Set up the OLI SWIR2 blind band (band number 13) */
    band_index++;
    strcpy(l8_band_attribs[band_index].band_name, "OLI_BLIND_SWIR2");
    l8_band_attribs[band_index].sensor_id = IAS_OLI;
    l8_band_attribs[band_index].sensor_type = IAS_PUSHBROOM_SENSOR;
    l8_band_attribs[band_index].band_type = IAS_SWIR2_BAND;
    l8_band_attribs[band_index].band_classification = IAS_BLIND_BAND;
    l8_band_attribs[band_index].spectral_type = IAS_SPECTRAL_SWIR;
    l8_band_attribs[band_index].normal_band_number = 7; 
    l8_band_attribs[band_index].vrp_band_number = 29; 
    l8_band_attribs[band_index].blind_band_number = 13; 
    l8_band_attribs[band_index].secondary_band_number = 0;
    l8_band_attribs[band_index].scas = OLI_NUMBER_OF_SCA;
    l8_band_attribs[band_index].detectors_per_sca = 
                            OLI_DETECTORS_PER_SCA_SWIR_BLIND;
    l8_band_attribs[band_index].lines_per_frame = 1;
    l8_band_attribs[band_index].pixel_resolution = OLI_PIXEL_RESOLUTION_NORMAL;
    l8_band_attribs[band_index].wavelength_nm_range[0]=0;
    l8_band_attribs[band_index].wavelength_nm_range[1]=0;

    /* Set up the OLI Cirrus blind band (band number 14) */
    band_index++;
    strcpy(l8_band_attribs[band_index].band_name, "OLI_BLIND_CIRRUS");
    l8_band_attribs[band_index].sensor_id = IAS_OLI;
    l8_band_attribs[band_index].sensor_type = IAS_PUSHBROOM_SENSOR;
    l8_band_attribs[band_index].band_type = IAS_CIRRUS_BAND;
    l8_band_attribs[band_index].band_classification = IAS_BLIND_BAND;
    l8_band_attribs[band_index].spectral_type = IAS_SPECTRAL_SWIR;
    l8_band_attribs[band_index].normal_band_number = 9; 
    l8_band_attribs[band_index].vrp_band_number = 30; 
    l8_band_attribs[band_index].blind_band_number = 14; 
    l8_band_attribs[band_index].secondary_band_number = 0;
    l8_band_attribs[band_index].scas = OLI_NUMBER_OF_SCA;
    l8_band_attribs[band_index].detectors_per_sca = 
                           OLI_DETECTORS_PER_SCA_CIRRUS_BLIND;
    l8_band_attribs[band_index].lines_per_frame = 1;
    l8_band_attribs[band_index].pixel_resolution = OLI_PIXEL_RESOLUTION_NORMAL;
    l8_band_attribs[band_index].wavelength_nm_range[0]=0;
    l8_band_attribs[band_index].wavelength_nm_range[1]=0;

    /* Set up the TIRS blind band (band number 15) */
    band_index++;
    strcpy(l8_band_attribs[band_index].band_name, "TIRS_BLIND");
    l8_band_attribs[band_index].sensor_id = IAS_TIRS;
    l8_band_attribs[band_index].sensor_type = IAS_PUSHBROOM_SENSOR;
    l8_band_attribs[band_index].band_type = IAS_THERMAL1_BAND;
    l8_band_attribs[band_index].band_classification = IAS_BLIND_BAND;
    l8_band_attribs[band_index].spectral_type = IAS_SPECTRAL_THERMAL;
    /* the blind band is associated with both bands 10 and 11, but pick the
       first one for this normal band number */
    l8_band_attribs[band_index].normal_band_number = 10; 
    l8_band_attribs[band_index].vrp_band_number = 0; 
    l8_band_attribs[band_index].blind_band_number = 15; 
    l8_band_attribs[band_index].secondary_band_number = 18;
    l8_band_attribs[band_index].scas = TIRS_NUMBER_OF_SCA;
    l8_band_attribs[band_index].detectors_per_sca = 
                           TIRS_DETECTORS_PER_SCA_BLIND;
    l8_band_attribs[band_index].lines_per_frame = 1;
    l8_band_attribs[band_index].pixel_resolution = TIRS_PIXEL_RESOLUTION;
    l8_band_attribs[band_index].wavelength_nm_range[0]=0;
    l8_band_attribs[band_index].wavelength_nm_range[1]=0;

    /* Set up the first TIRS secondary thermal band (band number 16) */
    band_index++;
    strcpy(l8_band_attribs[band_index].band_name, "TIRS_THERMAL1_SECONDARY");
    l8_band_attribs[band_index].sensor_id = IAS_TIRS;
    l8_band_attribs[band_index].sensor_type = IAS_PUSHBROOM_SENSOR;
    l8_band_attribs[band_index].band_type = IAS_THERMAL1_BAND;
    l8_band_attribs[band_index].band_classification = IAS_SECONDARY_BAND;
    l8_band_attribs[band_index].spectral_type = IAS_SPECTRAL_THERMAL;
    l8_band_attribs[band_index].normal_band_number = 16; 
    l8_band_attribs[band_index].vrp_band_number = 0; 
    l8_band_attribs[band_index].blind_band_number = 18; 
    l8_band_attribs[band_index].secondary_band_number = 0;
    l8_band_attribs[band_index].scas = TIRS_NUMBER_OF_SCA;
    l8_band_attribs[band_index].detectors_per_sca = 
                            TIRS_DETECTORS_PER_SCA_NORMAL;
    l8_band_attribs[band_index].lines_per_frame = 1;
    l8_band_attribs[band_index].pixel_resolution = TIRS_PIXEL_RESOLUTION;
    l8_band_attribs[band_index].wavelength_nm_range[0]=1030;
    l8_band_attribs[band_index].wavelength_nm_range[1]=1130;

    /* Set up the second TIRS secondary thermal band (band number 17) */
    band_index++;
    strcpy(l8_band_attribs[band_index].band_name, "TIRS_THERMAL2_SECONDARY");
    l8_band_attribs[band_index].sensor_id = IAS_TIRS;
    l8_band_attribs[band_index].sensor_type = IAS_PUSHBROOM_SENSOR;
    l8_band_attribs[band_index].band_type = IAS_THERMAL2_BAND;
    l8_band_attribs[band_index].band_classification = IAS_SECONDARY_BAND;
    l8_band_attribs[band_index].spectral_type = IAS_SPECTRAL_THERMAL;
    l8_band_attribs[band_index].normal_band_number = 17; 
    l8_band_attribs[band_index].vrp_band_number = 0; 
    l8_band_attribs[band_index].blind_band_number = 18; 
    l8_band_attribs[band_index].secondary_band_number = 0;
    l8_band_attribs[band_index].scas = TIRS_NUMBER_OF_SCA;
    l8_band_attribs[band_index].detectors_per_sca = 
                            TIRS_DETECTORS_PER_SCA_NORMAL;
    l8_band_attribs[band_index].lines_per_frame = 1;
    l8_band_attribs[band_index].pixel_resolution = TIRS_PIXEL_RESOLUTION;
    l8_band_attribs[band_index].wavelength_nm_range[0]=1150;
    l8_band_attribs[band_index].wavelength_nm_range[1]=1250;

    /* Set up the TIRS secondary blind band (band number 18) */
    band_index++;
    strcpy(l8_band_attribs[band_index].band_name, "TIRS_BLIND_SECONDARY");
    l8_band_attribs[band_index].sensor_id = IAS_TIRS;
    l8_band_attribs[band_index].sensor_type = IAS_PUSHBROOM_SENSOR;
    l8_band_attribs[band_index].band_type = IAS_THERMAL1_BAND;
    l8_band_attribs[band_index].band_classification
                = IAS_BLIND_BAND | IAS_SECONDARY_BAND;
    l8_band_attribs[band_index].spectral_type = IAS_SPECTRAL_THERMAL;
    /* the blind band is associated with both bands 16 and 17, but pick the
       first one for this normal band number */
    l8_band_attribs[band_index].normal_band_number = 16; 
    l8_band_attribs[band_index].vrp_band_number = 0; 
    l8_band_attribs[band_index].blind_band_number = 18; 
    l8_band_attribs[band_index].secondary_band_number = 0;
    l8_band_attribs[band_index].scas = TIRS_NUMBER_OF_SCA;
    l8_band_attribs[band_index].detectors_per_sca = 
                           TIRS_DETECTORS_PER_SCA_BLIND;
    l8_band_attribs[band_index].lines_per_frame = 1;
    l8_band_attribs[band_index].pixel_resolution = TIRS_PIXEL_RESOLUTION;
    l8_band_attribs[band_index].wavelength_nm_range[0]=0;
    l8_band_attribs[band_index].wavelength_nm_range[1]=0;

    /* Set up the OLI Coastal Aerosol Video Reference Pixel (VRP) band
       (band number 19) */
    band_index++;
    strcpy(l8_band_attribs[band_index].band_name, "OLI_VRP_Coastal_Aerosol");
    l8_band_attribs[band_index].sensor_id = IAS_OLI;
    l8_band_attribs[band_index].sensor_type = IAS_PUSHBROOM_SENSOR;
    l8_band_attribs[band_index].band_type = IAS_COASTAL_AEROSOL_BAND;
    l8_band_attribs[band_index].band_classification = IAS_VRP_BAND;
    l8_band_attribs[band_index].spectral_type = IAS_SPECTRAL_VNIR;
    l8_band_attribs[band_index].normal_band_number = 1; 
    l8_band_attribs[band_index].vrp_band_number = 19; 
    l8_band_attribs[band_index].blind_band_number = 0; 
    l8_band_attribs[band_index].secondary_band_number = 0;
    l8_band_attribs[band_index].scas = OLI_NUMBER_OF_SCA;
    l8_band_attribs[band_index].detectors_per_sca = OLI_DETECTORS_PER_SCA_VRP;
    l8_band_attribs[band_index].lines_per_frame = 1;
    l8_band_attribs[band_index].pixel_resolution = OLI_PIXEL_RESOLUTION_NORMAL;
    l8_band_attribs[band_index].wavelength_nm_range[0]=433;
    l8_band_attribs[band_index].wavelength_nm_range[1]=453;

    /* Set up the OLI Blue Video Reference Pixel (VRP) band (band number 20) */
    band_index++;
    strcpy(l8_band_attribs[band_index].band_name, "OLI_VRP_Blue");
    l8_band_attribs[band_index].sensor_id = IAS_OLI;
    l8_band_attribs[band_index].sensor_type = IAS_PUSHBROOM_SENSOR;
    l8_band_attribs[band_index].band_type = IAS_BLUE_BAND;
    l8_band_attribs[band_index].band_classification = IAS_VRP_BAND;
    l8_band_attribs[band_index].spectral_type = IAS_SPECTRAL_VNIR;
    l8_band_attribs[band_index].normal_band_number = 2; 
    l8_band_attribs[band_index].vrp_band_number = 20; 
    l8_band_attribs[band_index].blind_band_number = 0; 
    l8_band_attribs[band_index].secondary_band_number = 0;
    l8_band_attribs[band_index].scas = OLI_NUMBER_OF_SCA;
    l8_band_attribs[band_index].detectors_per_sca = OLI_DETECTORS_PER_SCA_VRP;
    l8_band_attribs[band_index].lines_per_frame = 1;
    l8_band_attribs[band_index].pixel_resolution = OLI_PIXEL_RESOLUTION_NORMAL;
    l8_band_attribs[band_index].wavelength_nm_range[0]=450;
    l8_band_attribs[band_index].wavelength_nm_range[1]=515;

    /* Set up the OLI Green Video Reference Pixel (VRP) band (band number 21) */
    band_index++;
    strcpy(l8_band_attribs[band_index].band_name, "OLI_VRP_Green");
    l8_band_attribs[band_index].sensor_id = IAS_OLI;
    l8_band_attribs[band_index].sensor_type = IAS_PUSHBROOM_SENSOR;
    l8_band_attribs[band_index].band_type = IAS_GREEN_BAND;
    l8_band_attribs[band_index].band_classification = IAS_VRP_BAND;
    l8_band_attribs[band_index].spectral_type = IAS_SPECTRAL_VNIR;
    l8_band_attribs[band_index].normal_band_number = 3; 
    l8_band_attribs[band_index].vrp_band_number = 21; 
    l8_band_attribs[band_index].blind_band_number = 0; 
    l8_band_attribs[band_index].secondary_band_number = 0;
    l8_band_attribs[band_index].scas = OLI_NUMBER_OF_SCA;
    l8_band_attribs[band_index].detectors_per_sca = OLI_DETECTORS_PER_SCA_VRP;
    l8_band_attribs[band_index].lines_per_frame = 1;
    l8_band_attribs[band_index].pixel_resolution = OLI_PIXEL_RESOLUTION_NORMAL;
    l8_band_attribs[band_index].wavelength_nm_range[0]=525;
    l8_band_attribs[band_index].wavelength_nm_range[1]=600;

    /* Set up the OLI Red Video Reference Pixel (VRP) band (band number 22) */
    band_index++;
    strcpy(l8_band_attribs[band_index].band_name, "OLI_VRP_Red");
    l8_band_attribs[band_index].sensor_id = IAS_OLI;
    l8_band_attribs[band_index].sensor_type = IAS_PUSHBROOM_SENSOR;
    l8_band_attribs[band_index].band_type = IAS_RED_BAND;
    l8_band_attribs[band_index].band_classification = IAS_VRP_BAND;
    l8_band_attribs[band_index].spectral_type = IAS_SPECTRAL_VNIR;
    l8_band_attribs[band_index].normal_band_number = 4; 
    l8_band_attribs[band_index].vrp_band_number = 22; 
    l8_band_attribs[band_index].blind_band_number = 0; 
    l8_band_attribs[band_index].secondary_band_number = 0;
    l8_band_attribs[band_index].scas = OLI_NUMBER_OF_SCA;
    l8_band_attribs[band_index].detectors_per_sca = OLI_DETECTORS_PER_SCA_VRP;
    l8_band_attribs[band_index].lines_per_frame = 1;
    l8_band_attribs[band_index].pixel_resolution = OLI_PIXEL_RESOLUTION_NORMAL;
    l8_band_attribs[band_index].wavelength_nm_range[0]=630;
    l8_band_attribs[band_index].wavelength_nm_range[1]=680;

    /* Set up the OLI NIR Video Reference Pixel (VRP) band (band number 23) */
    band_index++;
    strcpy(l8_band_attribs[band_index].band_name, "OLI_VRP_NIR");
    l8_band_attribs[band_index].sensor_id = IAS_OLI;
    l8_band_attribs[band_index].sensor_type = IAS_PUSHBROOM_SENSOR;
    l8_band_attribs[band_index].band_type = IAS_NIR_BAND;
    l8_band_attribs[band_index].band_classification = IAS_VRP_BAND;
    l8_band_attribs[band_index].spectral_type = IAS_SPECTRAL_VNIR;
    l8_band_attribs[band_index].normal_band_number = 5; 
    l8_band_attribs[band_index].vrp_band_number = 23; 
    l8_band_attribs[band_index].blind_band_number = 0; 
    l8_band_attribs[band_index].secondary_band_number = 0;
    l8_band_attribs[band_index].scas = OLI_NUMBER_OF_SCA;
    l8_band_attribs[band_index].detectors_per_sca = OLI_DETECTORS_PER_SCA_VRP;
    l8_band_attribs[band_index].lines_per_frame = 1;
    l8_band_attribs[band_index].pixel_resolution = OLI_PIXEL_RESOLUTION_NORMAL;
    l8_band_attribs[band_index].wavelength_nm_range[0]=845;
    l8_band_attribs[band_index].wavelength_nm_range[1]=885;

    /* Set up the OLI SWIR1 Video Reference Pixel (VRP) band (band number 24) */
    band_index++;
    strcpy(l8_band_attribs[band_index].band_name, "OLI_VRP_SWIR1");
    l8_band_attribs[band_index].sensor_id = IAS_OLI;
    l8_band_attribs[band_index].sensor_type = IAS_PUSHBROOM_SENSOR;
    l8_band_attribs[band_index].band_type = IAS_SWIR1_BAND;
    l8_band_attribs[band_index].band_classification = IAS_VRP_BAND;
    l8_band_attribs[band_index].spectral_type = IAS_SPECTRAL_SWIR;
    l8_band_attribs[band_index].normal_band_number = 6; 
    l8_band_attribs[band_index].vrp_band_number = 24; 
    l8_band_attribs[band_index].blind_band_number = 28; 
    l8_band_attribs[band_index].secondary_band_number = 0;
    l8_band_attribs[band_index].scas = OLI_NUMBER_OF_SCA;
    l8_band_attribs[band_index].detectors_per_sca = OLI_DETECTORS_PER_SCA_VRP;
    l8_band_attribs[band_index].lines_per_frame = 1;
    l8_band_attribs[band_index].pixel_resolution = OLI_PIXEL_RESOLUTION_NORMAL;
    l8_band_attribs[band_index].wavelength_nm_range[0]=1560;
    l8_band_attribs[band_index].wavelength_nm_range[1]=1660;

    /* Set up the OLI SWIR2 Video Reference Pixel (VRP) band (band number 25) */
    band_index++;
    strcpy(l8_band_attribs[band_index].band_name, "OLI_VRP_SWIR2");
    l8_band_attribs[band_index].sensor_id = IAS_OLI;
    l8_band_attribs[band_index].sensor_type = IAS_PUSHBROOM_SENSOR;
    l8_band_attribs[band_index].band_type = IAS_SWIR2_BAND;
    l8_band_attribs[band_index].band_classification = IAS_VRP_BAND;
    l8_band_attribs[band_index].spectral_type = IAS_SPECTRAL_SWIR;
    l8_band_attribs[band_index].normal_band_number = 7; 
    l8_band_attribs[band_index].vrp_band_number = 25; 
    l8_band_attribs[band_index].blind_band_number = 29; 
    l8_band_attribs[band_index].secondary_band_number = 0;
    l8_band_attribs[band_index].scas = OLI_NUMBER_OF_SCA;
    l8_band_attribs[band_index].detectors_per_sca = OLI_DETECTORS_PER_SCA_VRP;
    l8_band_attribs[band_index].lines_per_frame = 1;
    l8_band_attribs[band_index].pixel_resolution = OLI_PIXEL_RESOLUTION_NORMAL;
    l8_band_attribs[band_index].wavelength_nm_range[0]=2100;
    l8_band_attribs[band_index].wavelength_nm_range[1]=2300;

    /* Set up the OLI Pan Video Reference Pixel (VRP) band (band number 26) */
    band_index++;
    strcpy(l8_band_attribs[band_index].band_name, "OLI_VRP_PAN");
    l8_band_attribs[band_index].sensor_id = IAS_OLI;
    l8_band_attribs[band_index].sensor_type = IAS_PUSHBROOM_SENSOR;
    l8_band_attribs[band_index].band_type = IAS_PAN_BAND;
    l8_band_attribs[band_index].band_classification = IAS_VRP_BAND;
    l8_band_attribs[band_index].spectral_type = IAS_SPECTRAL_PAN;
    l8_band_attribs[band_index].normal_band_number = 8; 
    l8_band_attribs[band_index].vrp_band_number = 26; 
    l8_band_attribs[band_index].blind_band_number = 0;
    l8_band_attribs[band_index].secondary_band_number = 0; 
    l8_band_attribs[band_index].scas = OLI_NUMBER_OF_SCA;
    l8_band_attribs[band_index].detectors_per_sca = 
                            OLI_DETECTORS_PER_SCA_PAN_VRP;
    l8_band_attribs[band_index].lines_per_frame = 2;
    l8_band_attribs[band_index].pixel_resolution = OLI_PIXEL_RESOLUTION_PAN;
    l8_band_attribs[band_index].wavelength_nm_range[0]=500;
    l8_band_attribs[band_index].wavelength_nm_range[1]=680;

    /* Set up the OLI Cirrus Video Reference Pixel (VRP) band
       (band number 27) */
    band_index++;
    strcpy(l8_band_attribs[band_index].band_name, "OLI_VRP_CIRRUS");
    l8_band_attribs[band_index].sensor_id = IAS_OLI;
    l8_band_attribs[band_index].sensor_type = IAS_PUSHBROOM_SENSOR;
    l8_band_attribs[band_index].band_type = IAS_CIRRUS_BAND;
    l8_band_attribs[band_index].band_classification = IAS_VRP_BAND;
    l8_band_attribs[band_index].spectral_type = IAS_SPECTRAL_SWIR;
    l8_band_attribs[band_index].normal_band_number = 9; 
    l8_band_attribs[band_index].vrp_band_number = 27; 
    l8_band_attribs[band_index].blind_band_number = 30; 
    l8_band_attribs[band_index].secondary_band_number = 0;
    l8_band_attribs[band_index].scas = OLI_NUMBER_OF_SCA;
    l8_band_attribs[band_index].detectors_per_sca = OLI_DETECTORS_PER_SCA_VRP;
    l8_band_attribs[band_index].lines_per_frame = 1;
    l8_band_attribs[band_index].pixel_resolution = OLI_PIXEL_RESOLUTION_NORMAL;
    l8_band_attribs[band_index].wavelength_nm_range[0]=1360;
    l8_band_attribs[band_index].wavelength_nm_range[1]=1390;

    /* Set up the OLI SWIR1 Blind Video Reference Pixel (VRP) band
       (band number 28) */
    /* The blind band listed for OLI_VRP_BLIND_SWIR1 is the band number 
       associated with the masked detectors the VRPs were obtained from as 
       described by the pixel output sequence.*/
    band_index++;
    strcpy(l8_band_attribs[band_index].band_name, "OLI_VRP_BLIND_SWIR1");
    l8_band_attribs[band_index].sensor_id = IAS_OLI;
    l8_band_attribs[band_index].sensor_type = IAS_PUSHBROOM_SENSOR;
    l8_band_attribs[band_index].band_type = IAS_SWIR1_BAND;
    l8_band_attribs[band_index].band_classification 
                = IAS_VRP_BAND | IAS_BLIND_BAND;
    l8_band_attribs[band_index].spectral_type = IAS_SPECTRAL_SWIR;
    l8_band_attribs[band_index].normal_band_number = 12; 
    l8_band_attribs[band_index].vrp_band_number = 28; 
    l8_band_attribs[band_index].blind_band_number = 12; 
    l8_band_attribs[band_index].secondary_band_number = 0;
    l8_band_attribs[band_index].scas = OLI_NUMBER_OF_SCA;
    l8_band_attribs[band_index].detectors_per_sca = 
                            OLI_DETECTORS_PER_SCA_BLIND_VRP;
    l8_band_attribs[band_index].lines_per_frame = 1;
    l8_band_attribs[band_index].pixel_resolution = OLI_PIXEL_RESOLUTION_NORMAL;
    l8_band_attribs[band_index].wavelength_nm_range[0]=1560;
    l8_band_attribs[band_index].wavelength_nm_range[1]=1660;

    /* Set up the OLI SWIR2 Blind Video Reference Pixel (VRP) band
       (band number 29) */
    /* The blind band listed for OLI_VRP_BLIND_SWIR2 is the band number 
       associated with the masked detectors the VRPs were obtained from as 
       described by the pixel output sequence.*/
    band_index++;
    strcpy(l8_band_attribs[band_index].band_name, "OLI_VRP_BLIND_SWIR2");
    l8_band_attribs[band_index].sensor_id = IAS_OLI;
    l8_band_attribs[band_index].sensor_type = IAS_PUSHBROOM_SENSOR;
    l8_band_attribs[band_index].band_type = IAS_SWIR2_BAND;
    l8_band_attribs[band_index].band_classification
                = IAS_VRP_BAND | IAS_BLIND_BAND;
    l8_band_attribs[band_index].spectral_type = IAS_SPECTRAL_SWIR;
    l8_band_attribs[band_index].normal_band_number = 13; 
    l8_band_attribs[band_index].vrp_band_number = 29; 
    l8_band_attribs[band_index].blind_band_number = 13; 
    l8_band_attribs[band_index].secondary_band_number = 0;
    l8_band_attribs[band_index].scas = OLI_NUMBER_OF_SCA;
    l8_band_attribs[band_index].detectors_per_sca = 
                            OLI_DETECTORS_PER_SCA_BLIND_VRP;
    l8_band_attribs[band_index].lines_per_frame = 1;
    l8_band_attribs[band_index].pixel_resolution = OLI_PIXEL_RESOLUTION_NORMAL;
    l8_band_attribs[band_index].wavelength_nm_range[0]=2100;
    l8_band_attribs[band_index].wavelength_nm_range[1]=2300;

    /* Set up the OLI Cirrus Blind Video Reference Pixel (VRP) band
       (band number 30) */
    /* The blind band listed for OLI_VRP_BLIND_CIRRUS is the band number 
       associated with the masked detectors the VRPs were obtained from as 
       described by the pixel output sequence.*/
    band_index++;
    strcpy(l8_band_attribs[band_index].band_name, "OLI_VRP_BLIND_CIRRUS");
    l8_band_attribs[band_index].sensor_id = IAS_OLI;
    l8_band_attribs[band_index].sensor_type = IAS_PUSHBROOM_SENSOR;
    l8_band_attribs[band_index].band_type = IAS_CIRRUS_BAND;
    l8_band_attribs[band_index].band_classification
                = IAS_VRP_BAND | IAS_BLIND_BAND;
    l8_band_attribs[band_index].spectral_type = IAS_SPECTRAL_SWIR;
    l8_band_attribs[band_index].normal_band_number = 14; 
    l8_band_attribs[band_index].vrp_band_number = 30; 
    l8_band_attribs[band_index].blind_band_number = 14; 
    l8_band_attribs[band_index].secondary_band_number = 0;
    l8_band_attribs[band_index].scas = OLI_NUMBER_OF_SCA;
    l8_band_attribs[band_index].detectors_per_sca = 
                            OLI_DETECTORS_PER_SCA_BLIND_VRP;
    l8_band_attribs[band_index].lines_per_frame = 1;
    l8_band_attribs[band_index].pixel_resolution = OLI_PIXEL_RESOLUTION_NORMAL;
    l8_band_attribs[band_index].wavelength_nm_range[0]=1360;
    l8_band_attribs[band_index].wavelength_nm_range[1]=1390;

    return &l8_sat_attribs;
}

