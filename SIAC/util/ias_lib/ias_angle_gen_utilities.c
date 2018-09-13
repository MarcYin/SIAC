/* IAS Library Includes */
#include "ias_logging.h"
#include "ias_const.h"
#include "ias_angle_gen_distro.h"
#include "ias_angle_gen_private.h"

/*******************************************************************************
Name: ias_angle_gen_create_ang_filename

Purpose: Creates the ANG filename from the provided scene id.

Returns: 
    Type = integer
    SUCCESS / ERROR
*******************************************************************************/
int ias_angle_gen_create_ang_filename
(
    const char *base_filename, /* I: Base filename */
    int filename_size,         /* I: Filename size */
    char *filename             /* O: Output file name to use */
)       
{
    int status;  /* Routine return status placeholder */

    status = snprintf(filename, filename_size, "%s_ANG.txt", base_filename);
    if (status < 0 || status >= filename_size)
    {
        IAS_LOG_ERROR("Unable to construct the filename using %s", 
            base_filename);
        return ERROR;
    }

    return SUCCESS;
}

/*******************************************************************************
Name: ias_angle_gen_valid_band_index

Purpose: Checks that the band index is valid
 
Returns: 
    Type = integer
    On valid band index: TRUE
    On invalid band index: FALSE
*******************************************************************************/
double ias_angle_gen_valid_band_index
(
    const IAS_ANGLE_GEN_METADATA *metadata, /* I: Angle metadata */
    int band_index                          /* I: Band index to check */
)     
{
    /* Check that the band index is in the valid band range */
    if (band_index < 0 || band_index >= IAS_MAX_NBANDS)
    {
        IAS_LOG_ERROR("Band index %d is out of bounds", band_index);
        return FALSE;
    }
    
    /* Check that the band is present in the metadata */
    if (!metadata->band_present[band_index])
    {
        IAS_LOG_ERROR("Band index %d not present in metadata", band_index);
        return FALSE;
    }

    return TRUE;
}
