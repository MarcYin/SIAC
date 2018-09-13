/* Standard Library Includes */
#include <limits.h>

/* IAS Library Includes */
#include "ias_logging.h"       
#include "ias_angle_gen_distro.h"
#include "ias_miscellaneous.h"

/* Local Defines */ 
#define NUM_OUTPUT_SCAS 2 /* Number of output SCAs */

/*******************************************************************************
Name: ias_angle_gen_write_image

Purpose: Write out the solar illumination or the satellite viewing 
         angle images as a binary file.

Returns: 
    Type = integer
    SUCCESS / ERROR
 ******************************************************************************/
int ias_angle_gen_write_image
(
    const char *image_filename, /* I: Image file name */
    const short *azimuth,       /* I: Array of azimuth angles */
    const short *zenith,        /* I: Array of zenith angles */
    IAS_ANGLE_GEN_TYPE sat_or_sun_type, /* I: Image type to write */
    int band_index,             /* I: Outpue band index */
    int num_lines,              /* I: Number of image lines */
    int num_samps,              /* I: Number of image samples */
    IAS_DBL_XY ul_corner,       /* I: Image upper left corner */
    double pixel_size,          /* I: Image pixel size */
    const IAS_PROJECTION *projection /* I: Image framing information */
)        
{
    FILE *output_file;           /* Output file pointer */
    char ang_filename[PATH_MAX]; /* Output angle file name */
    int count;                   /* Total number of samples */
    int status;                  /* Status placeholder */
    const char *description;     /* Envi header description */
    int band_number;             /* Ouput band number */

    /* Check the input band index */
    if (band_index < 0 || band_index > IAS_MAX_NBANDS)
    {
        IAS_LOG_ERROR("Invalid band index %d", band_index);
        return ERROR;
    }

    band_number = ias_sat_attr_convert_band_index_to_number(band_index);
    if (band_number == ERROR)
    {
        IAS_LOG_ERROR("Converting band index %d to band number", band_index);
        return ERROR;
    }

    /* Calculate the total number of samples */
    count = num_lines * num_samps;
    
    if (sat_or_sun_type == IAS_ANGLE_GEN_SATELLITE)
    {
        description = "View Angle Band File";

        /* Construct satellite view angle output file name */
        status = snprintf(ang_filename, sizeof(ang_filename), 
            "%s_sensor_B%02d.img", image_filename, band_number);
    }
    else
    {
        description = "Sun Angle Band File";

        /* Construct solar angle output file name */
        status = snprintf(ang_filename, sizeof(ang_filename), 
            "%s_solar_B%02d.img", image_filename, band_number);
    }

    if (status < 0 || status >= sizeof(ang_filename))
    {
        IAS_LOG_ERROR("Creating the image filename from filename %s", 
            image_filename);
        return ERROR;
    }

    /* Open output file */
    output_file = fopen(ang_filename, "wb");
    if (!output_file)
    {
        IAS_LOG_ERROR("Opening output view angle band file %s", 
            ang_filename);
        return ERROR;
    }

    /* Write the Azimuth Angle layer */
    if (fwrite(azimuth, sizeof(short), count, output_file) != count)
    {
        IAS_LOG_ERROR("Writing azimuth angle layer for band %d to file %s", 
            band_number, ang_filename);
        fclose(output_file);
        return ERROR;
    }

    /* Write the Zenith Angle layer */
    if (fwrite(zenith, sizeof(short), count, output_file) != count)
    {
        IAS_LOG_ERROR("Writing zenith angle layer for band %d to file %s", 
            band_number, ang_filename);
        fclose(output_file);
        return ERROR;
    }

    /* Close the output image file */
    if (fclose(output_file) != 0)
    {
        IAS_LOG_ERROR("Closing the image file");
        return ERROR;
    }
   
    /* Write the envi header */
    if (ias_misc_write_envi_header(ang_filename, projection, description,
        num_lines, num_samps, NUM_OUTPUT_SCAS, ul_corner.x,ul_corner.y,
        pixel_size, pixel_size, "Azimuth, Zenith", IAS_I2) != SUCCESS)
    {
        IAS_LOG_ERROR("Creating the image header");
        return ERROR;
    }

    return SUCCESS;
}
