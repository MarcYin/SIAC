/* Standard Library Includes */
#include <stdlib.h>

/* IAS Library Includes */
#include "ias_logging.h"
#include "ias_angle_gen_private.h"

/*******************************************************************************
Name: ias_angle_gen_initialize

Purpose: Initialize the angle generation metadata structure. 

Note: ephem_count must be initialized to the size of the ephemeris and 
      solar vector structs.

Return:
    Type = Integer
    SUCCESS / ERROR
 ******************************************************************************/
int ias_angle_gen_initialize
(
    IAS_ANGLE_GEN_METADATA *metadata /* I/O: Angle metadata struct */
)           
{
    /* Check the number of points requested */
    if (metadata->ephem_count < 1)
    {
        IAS_LOG_ERROR("Invalid ephemeris count");
        return ERROR;
    }

    /* Allocate the ephemeris */
    metadata->ephemeris = (IAS_ANGLE_GEN_EPHEMERIS *)malloc(
        metadata->ephem_count * sizeof(IAS_ANGLE_GEN_EPHEMERIS));
    if (!metadata->ephemeris)
    {
        IAS_LOG_ERROR("Allocating ephemeris");
        return ERROR;
    }

    /* Allocate the solar vector */
    metadata->solar_vector = (IAS_ANGLE_GEN_EPHEMERIS *)malloc(
        metadata->ephem_count * sizeof(IAS_ANGLE_GEN_EPHEMERIS));
    if (!metadata->solar_vector)
    {
        IAS_LOG_ERROR("Allocating solar vector");
        free(metadata->ephemeris);
        metadata->ephemeris = NULL;
        return ERROR;
    }

    return SUCCESS;
}

/*******************************************************************************
Name: ias_angle_gen_free
  
Purpose: Free up the angle generation metadata structure.

Return:
    Type = void
 ******************************************************************************/
void ias_angle_gen_free
(
    IAS_ANGLE_GEN_METADATA *metadata /* I: Metadata structure */
)           
{
    free(metadata->ephemeris);
    metadata->ephemeris = NULL;

    free(metadata->solar_vector);
    metadata->solar_vector = NULL;
}
