/*************************************************************************

NAME: ias_satellite_attributes.c

PURPOSE: Implements the Satellite Attributes library for applications

**************************************************************************/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "ias_logging.h"
#include "ias_satellite_attributes.h"  /* Function prototypes and additional
                                          required headers */
#include "local_defines.h"

static IAS_SATELLITE_ATTRIBUTES *curr_satellite;

/*************************************************************************

NAME: ias_sat_attr_initialize

PURPOSE: Initializes the attributes for the satellite requested and sets
         that satellite as the default to use for the other routines in the 
         library.

RETURNS: SUCCESS or ERROR

**************************************************************************/
int ias_sat_attr_initialize
(
    int satellite_id        /* I: Satellite ID */
)
{
    /* Do satellite attribute initialization based on satellite ID */
    switch (satellite_id)
    {
        case IAS_L8 :
            curr_satellite = ias_sat_attr_initialize_landsat8();
            return SUCCESS;
        default:
            IAS_LOG_ERROR("Unrecognized satellite type");
            return ERROR;
    }
}

/*************************************************************************

NAME: ias_sat_attr_get_attributes

PURPOSE: Returns the attributes structure for the current satellite.
    This should not be used very often.

RETURNS: The satellite attributes pointer, or NULL if not initialized.

**************************************************************************/
const IAS_SATELLITE_ATTRIBUTES *ias_sat_attr_get_attributes()
{
    if (!curr_satellite)
        IAS_LOG_ERROR("Satellite Attributes haven't been initialized");

    return curr_satellite;
}

/*************************************************************************

NAME: ias_sat_attr_get_satellite_id

PURPOSE: Returns the satellite ID of the currently configured satellite.

RETURNS: The satellite ID or ERROR if it is not initialized.

**************************************************************************/
IAS_SATELLITE_ID ias_sat_attr_get_satellite_id()
{
    if (!curr_satellite)
    {
        IAS_LOG_ERROR("Satellite Attributes haven't been initialized");
        return ERROR;
    }

    return curr_satellite->satellite_id;
}

/*************************************************************************

NAME: ias_sat_attr_get_satellite_id_from_satellite_number

PURPOSE: Returns the satellite ID for a given satellite number.

RETURNS: The satellite ID or ERROR if the satellite number is not supported

**************************************************************************/
IAS_SATELLITE_ID ias_sat_attr_get_satellite_id_from_satellite_number
(
    int satellite_number /* I: The satellite number */
)
{
    if (satellite_number == 8)
        return IAS_L8;

    IAS_LOG_ERROR("Unsupported satellite number: %d", satellite_number);

    return ERROR;
}

/*************************************************************************

NAME: ias_sat_attr_get_satellite_number

PURPOSE: Returns the satellite number of the currently configured satellite.

RETURNS: The satellite ID or ERROR if it is not initialized or the satellite
    ID is not recognized.

**************************************************************************/
int ias_sat_attr_get_satellite_number()
{
    if (!curr_satellite)
    {
        IAS_LOG_ERROR("Satellite Attributes haven't been initialized");
        return ERROR;
    }

    switch (curr_satellite->satellite_id)
    {
        case IAS_L8:
            return 8;
        default:
            return ERROR;
    }
}

/*************************************************************************

NAME: ias_sat_attr_get_sensor_count

PURPOSE: Returns the number of sensors for the currently configured satellite.

RETURNS: The number of sensors or ERROR if it is not initialized.

**************************************************************************/
int ias_sat_attr_get_sensor_count()
{
    if (!curr_satellite)
    {
        IAS_LOG_ERROR("Satellite Attributes haven't been initialized");
        return ERROR;
    }

    return curr_satellite->sensors;
}

/*************************************************************************

NAME: ias_sat_attr_get_normal_band_count

PURPOSE: Returns the number of normal (non-blind/VRP) bands available on
         the satellite

RETURNS: The number of normal bands on the satellite

**************************************************************************/
int ias_sat_attr_get_normal_band_count()
{
    /* Check whether the Satellite Attributes haven't been initialized */
    if (!curr_satellite)
    {
        IAS_LOG_ERROR("Satellite Attributes haven't been initialized");
        return ERROR;
    }         

    return curr_satellite->bands;
}

/*************************************************************************

NAME: ias_sat_attr_get_scas_per_sensor

PURPOSE: Returns the number of SCAs available per sensor

RETURNS: The number of SCAs per sensor
*************************************************************************/
int ias_sat_attr_get_sensor_sca_count
(
    int sensor_id           /* I: sensor ID */
)
{
    int band_index;

    /* Has the library been initialized? */
    if (!curr_satellite)
    {
        IAS_LOG_ERROR("Satellite Attributes have not been initialized");
        return ERROR;
    }

    /* Loop through the bands to find a band for that sensor and return the
       number of SCAs it has */
    for (band_index = 0; band_index < curr_satellite->total_bands; band_index++)
    {
        const IAS_BAND_ATTRIBUTES *band_attr 
                = &curr_satellite->band_attributes[band_index];

        if (band_attr->sensor_id == sensor_id)
        {
            return band_attr->scas;
        }
    }

    IAS_LOG_ERROR("Unknown sensor ID: %d", sensor_id);
    return ERROR;

}

/*************************************************************************

NAME: ias_sat_attr_get_band_attributes

PURPOSE: Get the attributes for a band

RETURNS: Constant pointer to the correct IAS_BAND_ATTRIBUTES entry or
         NULL if it wasn't found

**************************************************************************/
const IAS_BAND_ATTRIBUTES *ias_sat_attr_get_band_attributes
(
    int band_number             /* I: Band number */
)
{
    int band_index = ias_sat_attr_convert_band_number_to_index(band_number);

    if (band_index == ERROR)
    {
        IAS_LOG_ERROR("Failed to convert the band number %d to an index",
            band_number);
        return NULL;
    }

    return &curr_satellite->band_attributes[band_index];
}

/*************************************************************************

NAME: ias_sat_attr_get_sensor_band_numbers

PURPOSE: Populate an array of bands with the band numbers that belong to a
         particular sensor.

RETURNS: SUCCESS or ERROR.

**************************************************************************/
int ias_sat_attr_get_sensor_band_numbers
(
    int sensor_id,              /* I: Sensor ID, or IAS_MAX_SENSORS if any
                                      sensor */
    int band_class,             /* I: band classification that must be
                                      included */
    int band_class_exclusion,   /* I: Band classification exclusion(s) an OR
                                      of classifications to exclude or the
                                      value 0 for no exclusions */
    int *band_number,           /* O: Pointer to an integer band number array 
                                      to populate */
    int size,                   /* I: Size of the band number array */
    int *number_of_bands        /* O: Pointer to the number of bands returned 
                                      in the band number array */
)
{
    int band_index;
    int band_count;

    /* Check whether the Satellite Attributes haven't been initialized */
    if (!curr_satellite)
    {
        IAS_LOG_ERROR("Satellite Attributes haven't been initialized");
        return ERROR;
    }         

    /* Initializing number_of_bands counter */
    band_count = 0;

    /* Loop to get all the bands */
    for (band_index = 0; band_index < curr_satellite->total_bands; band_index++)
    {
        const IAS_BAND_ATTRIBUTES *band_attr 
                = &curr_satellite->band_attributes[band_index];

        if (((sensor_id == IAS_MAX_SENSORS) 
             || (band_attr->sensor_id == sensor_id))
            && ((band_attr->band_classification & band_class) == band_class)
            && ((band_attr->band_classification & band_class_exclusion) == 0 ))
        {
            if (band_count < size)
            {
                band_number[band_count] = band_attr->band_number;
                band_count++;
            }
            else
                return ERROR;
        }
    }

    *number_of_bands = band_count;

    return SUCCESS;
}

/*************************************************************************

NAME: ias_sat_attr_get_any_sensor_band_numbers

PURPOSE: Populate an array of bands with the band numbers that belong to a
         particular sensor that match the any of the bands included in the 
         band_class.  This is different from
         ias_sat_get_attr_get_sensor_band_numbers since it returns bands that
         exactly match the band_class value requested.

         
RETURNS: SUCCESS or ERROR.

**************************************************************************/
int ias_sat_attr_get_any_sensor_band_numbers
(
    int sensor_id,              /* I: Sensor ID, or IAS_MAX_SENSORS if any
                                      sensor */
    int band_class,             /* I: Match any of the band classification 
                                      OR'ed together in this parameter */
    int band_class_exclusion,   /* I: Band classification exclusion(s) an OR
                                      of classifications to exclude or the
                                      value 0 for no exclusions */
    int *band_number,           /* O: Pointer to an integer band number array 
                                      to populate */
    int size,                   /* I: Size of the band number array */
    int *number_of_bands        /* O: Pointer to the number of bands returned 
                                      in the band number array */
)
{
    int band_index;
    int band_count = 0;        /* number of bands matching critera */

    /* Check whether the Satellite Attributes haven't been initialized */
    if (!curr_satellite)
    {
        IAS_LOG_ERROR("Satellite Attributes haven't been initialized");
        return ERROR;
    }         

    /* Loop to get all the bands */
    for (band_index = 0; band_index < curr_satellite->total_bands; band_index++)
    {
        const IAS_BAND_ATTRIBUTES *band_attr 
                = &curr_satellite->band_attributes[band_index];

        if (((sensor_id == IAS_MAX_SENSORS) 
             || (band_attr->sensor_id == sensor_id))
            && (band_attr->band_classification & band_class)
            && ((band_attr->band_classification & band_class_exclusion) == 0 ))
        {
            if (band_count < size)
            {
                band_number[band_count] = band_attr->band_number;
                band_count++;
            }
            else
                return ERROR;
        }
    }

    *number_of_bands = band_count;

    return SUCCESS;
}

/*************************************************************************

NAME: ias_sat_attr_band_classification_matches

PURPOSE: Checks whether the band number passed in has a classification that
    exactly matches the provided classification.

RETURNS: TRUE if it matches or FALSE if it doesn't match (or an error happens)

**************************************************************************/
int ias_sat_attr_band_classification_matches
(
    int band_number,        /* I: band number to check */
    IAS_BAND_CLASSIFICATION band_classification /* I: exact classification 
                                                      to match */
)
{
    int band_index;

    /* Check whether the Satellite Attributes haven't been initialized */
    if (!curr_satellite)
    {
        IAS_LOG_ERROR("Satellite Attributes haven't been initialized");
        return FALSE;
    }         

    /* get the band index */
    band_index = ias_sat_attr_convert_band_number_to_index(band_number);
    if (band_index == ERROR)
    {
        IAS_LOG_ERROR("Unable to convert band number %d to an index",
                      band_number);
        return FALSE;
    }

    if (curr_satellite->band_attributes[band_index].band_classification
        == band_classification)
    {
        return TRUE;
    }
    else
        return FALSE;
}

/*************************************************************************

NAME: ias_sat_attr_convert_band_number_to_name

PURPOSE: Convert band number to band name

RETURNS: Constant pointer to the band name, or NULL if the conversion fails.

**************************************************************************/
const char *ias_sat_attr_convert_band_number_to_name
(
    int band_number     /* I: band number to convert to a name */
)
{
    int band_index = ias_sat_attr_convert_band_number_to_index(band_number);

    if (band_index == ERROR)
    {
        IAS_LOG_ERROR("Failed to convert the band number to an index");
        return NULL;
    }

    return curr_satellite->band_attributes[band_index].band_name;
}

/*************************************************************************

NAME: ias_sat_attr_convert_band_name_to_number

PURPOSE: Convert band name to band number

RETURNS: Integer band number if successful, or ERROR if the conversion
         fails

**************************************************************************/
int ias_sat_attr_convert_band_name_to_number
(
    const char *band_name              /* I: Band name to search on */
)
{
    int band_index;                    /* Band loop counter */
    int band_number = ERROR;           /* Band number to return,
                                          initialized to the ERROR
                                          code to start with*/


    /* Check whether the Satellite Attributes haven't been initialized */
    if (!curr_satellite)
    {
        IAS_LOG_ERROR("Satellite Attributes haven't been initialized");
        return ERROR;
    }

    /* Search the band_attributes list for the specified band name */
    for (band_index = 0; band_index < curr_satellite->total_bands;
         band_index++)
    {
        if (strcmp(band_name,
                   curr_satellite->band_attributes[band_index].band_name)
                       == 0)
        {
            band_number
                = curr_satellite->band_attributes[band_index].band_number;
            break;
        }
    }

    /* Failed to find a band number for the specified name */
    if (band_number == ERROR)
    {
        IAS_LOG_ERROR("Band name %s does not have a known band number",
            band_name);
    }

    return band_number;
}

/*************************************************************************

NAME: ias_sat_attr_convert_band_index_to_number

PURPOSE: Convert band index to band number

RETURNS: Integer of the band number, or ERROR if the conversion fails.

**************************************************************************/
int ias_sat_attr_convert_band_index_to_number
(
    int band_index      /* I: band index to convert to a number */
)
{
    /* Check whether the Satellite Attributes haven't been initialized */
    if (!curr_satellite)
    {
        IAS_LOG_ERROR("Satellite Attributes haven't been initialized");
        return ERROR;
    }         

    if ((band_index < 0) || (band_index > curr_satellite->total_bands - 1))
    {
        IAS_LOG_ERROR("Invalid band index");        
        return ERROR;
    }

    return curr_satellite->band_attributes[band_index].band_number;
}

/*************************************************************************

NAME: ias_sat_attr_convert_band_number_to_index

PURPOSE: Convert band number to band index

RETURNS: Integer of the band index, or ERROR if the conversion fails.

**************************************************************************/
int ias_sat_attr_convert_band_number_to_index
(
    int band_number     /* I: band number to convert to an index */
)
{
    /* Check whether the Satellite Attributes haven't been initialized */
    if (!curr_satellite)
    {
        IAS_LOG_ERROR("Satellite Attributes haven't been initialized");
        return ERROR;
    }         

    if ((band_number < 1) || (band_number > curr_satellite->total_bands))
    {
        IAS_LOG_ERROR("Invalid band number: %d", band_number);        
        return ERROR;
    }

    /* Assume that the band number is band_index + 1 as an optimization.  It is
       safe to do that inside this library, but nothing outside the library
       should. */
    return curr_satellite->band_attributes[band_number-1].band_index;
}

/*************************************************************************

NAME: ias_sat_attr_get_band_number_from_type

PURPOSE: Given a band type, find the first band number that matches it.

RETURNS: band number if success, ERROR if it fails to find a matching band.

**************************************************************************/
int ias_sat_attr_get_band_number_from_type
(
    IAS_BAND_TYPE band_type     /* I: band type to convert to a band number */
)
{
    int band_index;
    IAS_BAND_ATTRIBUTES *band;

    /* Check whether the Satellite Attributes haven't been initialized */
    if (!curr_satellite)
    {
        IAS_LOG_ERROR("Satellite Attributes haven't been initialized");
        return ERROR;
    }         

    band = curr_satellite->band_attributes;

    /* search the bands to find a matching one */
    for (band_index = 0; band_index < curr_satellite->total_bands; band_index++)
    {
        if (band[band_index].band_type == band_type)
        {
            return band[band_index].band_number;
        }
    }

    IAS_LOG_ERROR("Unrecognized band type: %d", band_type);
    return ERROR;
}

/*************************************************************************

NAME: ias_sat_attr_get_scas_per_band

PURPOSE: Get the number of SCAs for a band

RETURNS: number of SCAs if success, or ERROR if it fails to find a matching
         band 

**************************************************************************/
int ias_sat_attr_get_scas_per_band
(
    int band_number             /* I: Band number */
)
{
    int band_index = ias_sat_attr_convert_band_number_to_index(band_number);

    if (band_index == ERROR)
    {
        IAS_LOG_ERROR("Failed to convert the band number to an index");
        return ERROR;
    }

    return curr_satellite->band_attributes[band_index].scas;
}

/*************************************************************************

NAME: ias_sat_attr_get_detectors_per_sca

PURPOSE: Get the number of detectors in an SCA for a band

RETURNS: number of detectors if success, or ERROR if it fails to find 
         a matching band 

**************************************************************************/
int ias_sat_attr_get_detectors_per_sca
(
    int band_number             /* I: Band number */
)
{
    int band_index = ias_sat_attr_convert_band_number_to_index(band_number);

    if (band_index == ERROR)
    {
        IAS_LOG_ERROR("Failed to convert the band number to an index");
        return ERROR;
    }

    return curr_satellite->band_attributes[band_index].detectors_per_sca;
}

/*************************************************************************

NAME: ias_sat_attr_get_lines_per_frame

PURPOSE: Get the number of imagery lines per frame collected for the requested
         band.

RETURNS: number of lines per frame if success, or ERROR if it fails to find 
         a matching band 

**************************************************************************/
int ias_sat_attr_get_lines_per_frame
(
    int band_number             /* I: Band number */
)
{
    int band_index = ias_sat_attr_convert_band_number_to_index(band_number);

    if (band_index == ERROR)
    {
        IAS_LOG_ERROR("Failed to convert the band number %d to an index",
            band_number);
        return ERROR;
    }

    return curr_satellite->band_attributes[band_index].lines_per_frame;
}

/*************************************************************************

NAME: ias_sat_attr_get_sensor_name

PURPOSE: Get the name for a specific sensor ID

RETURNS: A pointer to a string for the sensor name.  If the sensor ID is
    invalid, "Unknown" is returned.

**************************************************************************/
const char *ias_sat_attr_get_sensor_name
(
    int sensor_id       /* I: sensor ID to get the name of */
)
{
    static const char *oli_name = "OLI";
    static const char *tirs_name = "TIRS";
    static const char *unknown = IAS_SENSOR_NAME_UNKNOWN;

    switch (sensor_id)
    {
        case IAS_OLI:
            return oli_name;
        case IAS_TIRS:
            return tirs_name;
        default:
            return unknown;
    }
}

/*************************************************************************

NAME: ias_sat_attr_convert_band_number_to_sensor_band_index

PURPOSE: Given a sensor ID and a band number, find the index of that band
    number within just the bands for that sensor ID.

RETURNS: calculated band index for sensor

NOTES:
    - This routine isn't used by the IAS code.  It has been included here for
      ingest deal with the TIRS primary and secondary bands a bit more
      easily.
************************************************************************/
int ias_sat_attr_convert_band_number_to_sensor_band_index
(
    int sensor_id,
    int band_number
)
{
    int band_index;
    int sensor_band_index = 0;

    IAS_BAND_ATTRIBUTES *bands;

    /* Check whether the Satellite Attributes haven't been initialized */
    if (!curr_satellite)
    {
        IAS_LOG_ERROR("Satellite Attributes haven't been initialized");
        return ERROR;
    }

    bands = curr_satellite->band_attributes;

    for (band_index = 0; band_index < curr_satellite->total_bands; band_index++)
    {
        if (bands[band_index].sensor_id == sensor_id)
        {

            if (bands[band_index].band_number == band_number)
                return sensor_band_index;
            sensor_band_index++;
        }
    }

    IAS_LOG_ERROR("Invalid sensor ID (%d) and band number combination (%d)",
                  sensor_id, band_number);
    return ERROR;
}

/*************************************************************************

NAME: ias_sat_attr_get_satellite_name

PURPOSE: Get the name for the configured satellite

RETURNS: A pointer to a string for the satellite name.  If the satellite ID
         is invalid, "Unknown" is returned.

**************************************************************************/
const char *ias_sat_attr_get_satellite_name()
{
    static const char *unknown = IAS_SATELLITE_NAME_UNKNOWN;

    /* Return "unknown" if satellite attributes are not initialized. */
    if (!curr_satellite)
        return unknown;

    return curr_satellite->satellite_name;
}

/**************************************************************************

 NAME:      ias_sat_attr_get_spectral_type_from_band_number

 PURPOSE:  Gets the spectral type for the specified band number

 RETURNS:  Spectral type value if successful, otherwise returns an ERROR

***************************************************************************/
IAS_SPECTRAL_TYPE ias_sat_attr_get_spectral_type_from_band_number
(
    int band_number                    /* I: Current 1-based band number */
)
{
    IAS_BAND_ATTRIBUTES *band_attributes = NULL;
    int band_index;

    /* Has the library been initialized? */
    if (!curr_satellite)
    {
        IAS_LOG_ERROR("Satellite Attributes have not been initialized");
        return ERROR;
    }
    band_attributes = curr_satellite->band_attributes;

    /* Get the band attributes array index for the current band number */
    band_index = ias_sat_attr_convert_band_number_to_index(band_number);
    if (band_index == ERROR)
    {
        IAS_LOG_ERROR("Obtaining correct band index from specified band "
                      "number");
        return ERROR;
    }

    /* We've got a proper index, now return the spectral type */
    return band_attributes[band_index].spectral_type;
}

/****************************************************************************
NAME: ias_sat_attr_get_sensor_max_normal_detectors

PURPOSE: Gets the maximum number of "normal" band imaging detectors for the
         specified sensor ID

RETURNS: Maximum number of detectors if successful, ERROR on failure
*****************************************************************************/
int ias_sat_attr_get_sensor_max_normal_detectors
(
    int sensor_id                     /* I: Sensor ID number for OLI OR
                                         TIRS */
)
{
    int normal_band_list[IAS_MAX_NBANDS];
    int band_count;
    int max_detectors = -1;
    int index;

    if (ias_sat_attr_get_sensor_band_numbers(sensor_id, IAS_NORMAL_BAND, 0,
            normal_band_list, IAS_MAX_NBANDS, &band_count) == ERROR)
    {
        IAS_LOG_ERROR("Failed to get the list of normal bands");
        return ERROR;
    }

    /* find the maximum of all the bands */
    for (index = 0; index < band_count; index++)
    {
        IAS_BAND_ATTRIBUTES *band;
        int band_number = normal_band_list[index];
        int band_index;

        band_index = ias_sat_attr_convert_band_number_to_index(band_number);
        if (band_index == ERROR)
        {
            IAS_LOG_ERROR("Unabled to convert band number %d to an index",
                          band_number);
            return ERROR;
        }

        band = &curr_satellite->band_attributes[band_index];
        if (band->detectors_per_sca > max_detectors)
            max_detectors = band->detectors_per_sca;
    }

    if (max_detectors == -1)
    {
        IAS_LOG_ERROR("No normal bands found for sensor ID %d", sensor_id);
        return ERROR;
    }

    return max_detectors;
}


/****************************************************************************
 NAME:            ias_sat_attr_get_band_name_from_type_and_class

 PURPOSE:  Given a band's spectral type and classification, return its
           name if it exists

 RETURNS:  Integer status code of SUCCESS or ERROR

****************************************************************************/
int ias_sat_attr_get_band_name_from_type_and_class
(
    IAS_BAND_TYPE band_type,                       /* I: Band type */
    IAS_BAND_CLASSIFICATION band_classification,   /* I: Classification */
    char *band_name                                /* O: Band name */
)
{
    int band_index = 0;

    /* Make sure the Satellite Attributes Library has been initialized */
    if (!curr_satellite)
    {
        IAS_LOG_ERROR("Satellite Attributes library not initialized");
        return ERROR;
    }

    /* Search all of the currently defined bands */
    for (band_index = 0;  band_index < IAS_MAX_TOTAL_BANDS; band_index++)
    {

        /* Found a match -- get the name and quit */
        if ((curr_satellite->band_attributes[band_index].band_type
                == band_type)
            && (curr_satellite->band_attributes[band_index]
                    .band_classification == band_classification))
        {
            strcpy(band_name,
                curr_satellite->band_attributes[band_index].band_name);
            return SUCCESS;
        }
    }

    /* Didn't find it, so it's an error */
    IAS_LOG_ERROR("No band name corresponding to specified band type and "
        "classification found");
    return ERROR;
}


/*********************************************************************
 NAME:     ias_sat_attr_get_band_type_from_band_number

 PURPOSE:  Given the current (1-based) band number, return the
           associated band type

 RETURNS:  Enumerated value representing the band type for the
           specified band
**********************************************************************/
IAS_BAND_TYPE ias_sat_attr_get_band_type_from_band_number
(
    int band_number                   /* I: Current band number */
)
{
    int band_index;

    /* Get the corresponding band index */
    band_index = ias_sat_attr_convert_band_number_to_index(band_number);
    if (band_index == ERROR)
    {
        IAS_LOG_ERROR("Unable to convert band number %d to an index",
            band_number);
        return IAS_UNKNOWN_BAND_TYPE;
    }

    /* Get the band type for this band */
    return curr_satellite->band_attributes[band_index].band_type;

}   /* END -- ias_sat_attr_get_band_type_from_band_number */


/*********************************************************************
 NAME:     ias_sat_attr_get_sensor_id_from_band_number

 PURPOSE:  Given the current (1-based) band number, return the
           associated sensor ID

 RETURNS:  Enumerated value representing the sensor ID for the
           specified band
**********************************************************************/
IAS_SENSOR_ID ias_sat_attr_get_sensor_id_from_band_number
(
    int band_number                  /* I: Current band number */
)
{
    int band_index;

    /* Get the corresponding band index */
    band_index = ias_sat_attr_convert_band_number_to_index(band_number);
    if (band_index == ERROR)
    {
        IAS_LOG_ERROR("Unable to convert band number %d to an index",
            band_number);
        return IAS_INVALID_SENSOR_ID;
    }

    /* Get the attributes for this band */
    return curr_satellite->band_attributes[band_index].sensor_id;
}   /* END -- ias_sat_attr_get_sensor_id_from_band_number */

/*************************************************************************

NAME: ias_sat_attr_get_quantization_cal_min

PURPOSE: Get the QUANTIZE_CAL_MIN value for a specific band number

RETURNS: SUCCESS or ERROR

**************************************************************************/
int ias_sat_attr_get_quantization_cal_min
(
    int band_number,            /* I: Current band number */
    int *qcal_min               /* O: qcal min value */
)
{
    const IAS_BAND_ATTRIBUTES *band_attr
        = ias_sat_attr_get_band_attributes(band_number);

    if (!band_attr)
    {
        IAS_LOG_ERROR("Failed to get the band attributes for band number %d",
            band_number);
        return ERROR;
    }

    *qcal_min = band_attr->qcal_min;

    return SUCCESS;
}

/*************************************************************************

NAME: ias_sat_attr_get_quantization_cal_max

PURPOSE: Get the QUANTIZE_CAL_MAX value for a specific band number

RETURNS: SUCCESS or ERROR

**************************************************************************/
int ias_sat_attr_get_quantization_cal_max
(
    int band_number,            /* I: Current band number */
    int *qcal_max               /* O: qcal max value */
)
{
    const IAS_BAND_ATTRIBUTES *band_attr
        = ias_sat_attr_get_band_attributes(band_number);

    if (!band_attr)
    {
        IAS_LOG_ERROR("Failed to get the band attributes for band number %d",
            band_number);
        return ERROR;
    }

    *qcal_max = band_attr->qcal_max;

    return SUCCESS;
}

/**************************************************************************
 NAME:      ias_sat_attr_get_satellite_sensor_name

 PURPOSE:   Gets the sensor name associated with the current satellite

 RETURNS:   Pointer to a string containing the sensor name if successful,
            NULL pointer if there's an error
***************************************************************************/
const char *ias_sat_attr_get_satellite_sensor_name()
{
   static const char *unknown = IAS_SENSOR_NAME_UNKNOWN;

   /* Return "unknown" if satellite attributes are not initialized. */
   if (!curr_satellite)
   {
       return unknown;
   }

   /* Return the associated sensor name. */
   return curr_satellite->sensor_name;
}

/*************************************************************************

NAME: ias_sat_attr_report_saturation_for_band

PURPOSE: Determines if saturation should be reported for a band

RETURNS: TRUE / FALSE
*************************************************************************/
int ias_sat_attr_report_saturation_for_band
(
    int band_number,                   /* I: Band number */
    int *should_report_saturation_flag /* O: Flag to say if saturation should
                                             be reported */
)
{
    int band_index;
    *should_report_saturation_flag = FALSE;

    /* Has the library been initialized? */
    if (!curr_satellite)
    {
        IAS_LOG_ERROR("Satellite Attributes have not been initialized");
        return ERROR;
    }

    band_index = ias_sat_attr_convert_band_number_to_index(band_number);
    if (band_index == ERROR)
    {
        IAS_LOG_ERROR("Unable to convert band number %d to band index",
            band_number);
        return ERROR;
    }

    /* Check if the band itself can saturate */
    if (curr_satellite->band_attributes[band_index].can_saturate)
    {
        *should_report_saturation_flag = TRUE;
    }
    else
    {
        /* Report saturation for this band if any band of this sensor type can
           can have saturation */
        IAS_SENSOR_ID sensor_id = curr_satellite->band_attributes[
            band_index].sensor_id;
        IAS_BAND_CLASSIFICATION band_class = curr_satellite->band_attributes[
            band_index].band_classification;

        /* Loop through the bands to see if any bands of this sensor type can
           saturate */
        for (band_index = 0; band_index < curr_satellite->total_bands;
             band_index++)
        {
            const IAS_BAND_ATTRIBUTES *band_attr
                    = &curr_satellite->band_attributes[band_index];

            if (band_attr->sensor_id == sensor_id
                && band_attr->band_classification == band_class
                && band_attr->can_saturate)
            {
                *should_report_saturation_flag = TRUE;
                break;
            }
        }
    }

    return SUCCESS;
}
