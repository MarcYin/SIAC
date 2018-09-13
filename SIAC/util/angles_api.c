/* IAS Library Includes */
#include "ias_logging.h"
#include "ias_angle_gen_distro.h"  

/* Local Includes */
#include "l8_angles.h"

/*******************************************************************************
Name: calculate_angles

Purpose: Calculate the satellite and solar zenith and azimuth angles.

Note: ias_angle_gen_calculate_angles_rpc will return the angles in radians

Return: SUCCESS / ERROR
 ******************************************************************************/
int calculate_angles
(
    const IAS_ANGLE_GEN_METADATA *metadata, /* I: Angle metadata structure */ 
    int line,                               /* I: L1T line coordinate */
    int samp,                               /* I: L1T sample coordinate */
    int band_index,                         /* I: Spectral band number */
    ANGLE_TYPE angle_type,                  /* I: Type of angles to generate */
    double *sat_angles,                     /* O: Satellite angles (radians) */
    double *sun_angles                      /* O: Solar angles (radians) */
)       
{
    double elev;            /* Elevation always set at 0 */
    int outside_image_flag; /* Outside of image flag */ 

    /* Default the elevation to 0 to ensure the full scene coverage */
    elev = 0;

    /* If angle type is not of solar type then it is either calculating
       both angles or just the satellite angles */
    if (angle_type != AT_SOLAR)
    {
        /* Calculate the satellite viewing angles */
        if (ias_angle_gen_calculate_angles_rpc(metadata, line, samp, &elev,
            band_index, IAS_ANGLE_GEN_SATELLITE, &outside_image_flag, 
            sat_angles) != SUCCESS)
        {
            IAS_LOG_ERROR("Evaluating angles for band index %d", band_index);
            return ERROR;
        }
    }

    /* If angle type is not of satellite type then it is either calculating
       both angles or just the solar angles */
    if (angle_type != AT_SATELLITE)
    {
        /* Calculate the solar angles */
        if (ias_angle_gen_calculate_angles_rpc(metadata, line, samp, &elev, 
            band_index, IAS_ANGLE_GEN_SOLAR, &outside_image_flag, sun_angles)
            != SUCCESS) 
        {
            IAS_LOG_ERROR("Evaluating solar angles for band index %d",
                band_index);
            return ERROR;
        }
    }

    return SUCCESS;
}

/*******************************************************************************
Name: get_active_lines

Purpose: Returns the number of active lines for each band

Return: Active lines array
 ******************************************************************************/
const double *get_active_lines
(
    const IAS_ANGLE_GEN_METADATA *metadata, /* I: Angle metadata structure */ 
    int band_index                          /* I: Band index */
)
{
    return metadata->band_metadata[band_index].active_l1t_corner_lines;
}

/*******************************************************************************
Name: get_active_samples

Purpose: Returns the number of active samples for each band

Return: Active samples array
 ******************************************************************************/
const double *get_active_samples
(
    const IAS_ANGLE_GEN_METADATA *metadata, /* I: Angle metadata structure */ 
    int band_index                          /* I: Band index */
)
{
    return metadata->band_metadata[band_index].active_l1t_corner_samps;
}

/*******************************************************************************
Name: get_frame

Purpose: Retrieve the band frame information.

Return: SUCCESS / ERROR
 ******************************************************************************/
int get_frame
(
    const IAS_ANGLE_GEN_METADATA *metadata, /* I: Angle metadata structure */ 
    int band_index,                         /* I: Band index */
    ANGLES_FRAME *frame                     /* O: Image frame info */
)        
{
    /* Return error if band isn't present so caller can issue a warning */
    if (!metadata->band_present[band_index])
        return ERROR;
    

    /* Establish the output image file frame */
    frame->band_number = metadata->band_metadata[band_index].band_number;
    frame->num_lines = metadata->band_metadata[band_index].l1t_lines;
    frame->num_samps = metadata->band_metadata[band_index].l1t_samps;
    frame->projection = metadata->projection;
    frame->pixel_size = metadata->band_metadata[band_index].pixel_size;
    frame->ul_corner.x = metadata->corners.upleft.x;
    frame->ul_corner.y = metadata->corners.upleft.y;

    return SUCCESS;
}
