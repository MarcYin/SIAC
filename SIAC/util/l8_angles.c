/* Standard Library Includes */
#include <libgen.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* IAS Library Includes */
#include "ias_logging.h"
#include "ias_const.h"
#include "ias_angle_gen_distro.h"
#include "ias_miscellaneous.h"

/* Local Includes */
#include "l8_angles.h"

/* Local defines */                                               
#define USE_MESSAGE "Usage: l8_angles \n" \
"     <MetadataFilename>: (Required) Angle coefficient filename \n\n" \
"     <AngleType>: (Required) The type of angles to generate \n" \
"                             VALID: (BOTH, SATELLITE, SOLAR) \n\n" \
"     <SubsampleFactor>: (Required) Sub-sample factor used when calculating\n" \
"                                   the angles (integer)  \n\n" \
"     -f <FillPixelValue>: (Optional) Fill pixel value to use (short int) \n" \
"                                     units used is degrees scaled by 100\n" \
"                                     Default: 0\n" \
"                                     Range: (-32768:32767)\n\n" \
"     -b <BandList>: (Optional) Band list used to calculate angles for, this\n"\
"                            defaults to all bands 1 - 11. Must be comma \n" \
"                            separated with no spaces in between. \n" \
"                            Example: 1,2,3,4,5,6,7,8,9 \n"

/* Prototypes */
int process_parameters(const int argument_count, char *arguments[],
    L8_ANGLES_PARAMETERS *parameters);

/**************************************************************************
NAME:      l8_angles (main)

PURPOSE:   Initializes the IAS_ANGLE_GEN metdata structure and uses the 
            coefficients to generate the satellite viewing angle and solar
            angle image files.    

RETURNS:   EXIT_SUCCESS or EXIT_FAILURE
***************************************************************************/
int main
(
    int argc, 
    char *argv[]
)
{
    int band_index;                 /* Metadata band index */
    int sub_sample;                 /* Subsampling factor */
    int calc_sat_angles;            /* Local satellite angle flag */
    int calc_solar_angles;          /* Local solar angle flag */
    L8_ANGLES_PARAMETERS parameters;/* Parameters read in from file */
    IAS_ANGLE_GEN_METADATA metadata;/* Angle metadata structure */ 
    char root_filename[PATH_MAX];   /* Root filename */
    char *base_ptr;                 /* Basename pointer */
    double r2d = 4500.0 / atan(1.0);/* Conversion to hundredths of degrees */

    /* If there are not enough arguments, generate a usage message */
    if (argc < 4)
    {
        printf(USE_MESSAGE);
        return EXIT_FAILURE;
    }

    /* Initialize the logging library */
    if (ias_log_initialize("L8 Angles") != SUCCESS)
    {
       IAS_LOG_ERROR("Error initializing logging library");
       return EXIT_FAILURE;
    }

    /* Initialize the satellite attributes */
    if (ias_sat_attr_initialize(IAS_L8) != SUCCESS)
    {
        IAS_LOG_ERROR("Initializing satellite attributes library");
        return EXIT_FAILURE;
    }

    /* Process the arguments */
    if (process_parameters(argc, argv, &parameters) != SUCCESS)
    {
        printf(USE_MESSAGE);
        return EXIT_FAILURE;
    }

    /* Setup local sub sampling factor variable */
    sub_sample = parameters.sub_sample_factor;

    /* Read the metadata file */
    if (ias_angle_gen_read_ang(parameters.metadata_filename, &metadata) 
        != SUCCESS)
    {
        IAS_LOG_ERROR("Reading the metadata file %s", 
            parameters.metadata_filename);
        return EXIT_FAILURE;
    }

    /* Check what angles should be calculated */
    calc_sat_angles = TRUE;
    calc_solar_angles = TRUE;
    if (parameters.angle_type == AT_SOLAR)
    {
        calc_sat_angles = FALSE;
    }
    else if (parameters.angle_type == AT_SATELLITE)
    {
        calc_solar_angles = FALSE;
    }
    
    /* Extract the basename from the file path */
    base_ptr = strrchr(parameters.metadata_filename, '/');
    if (base_ptr)
    {
        base_ptr++; /* Move past the the last forward slash */
    }
    else
    {
        base_ptr = parameters.metadata_filename;
    }

    /* Check that the base buffer won't overflow the target buffer */
    if ((strlen(base_ptr) * sizeof(char)) > sizeof(root_filename))
    {
        IAS_LOG_ERROR("Angle coefficient filename too long");
        return ERROR;
    }
    strcpy(root_filename, base_ptr);

    /* Strip off the extension */
    base_ptr = strrchr(root_filename, '.');
    if (base_ptr)
    {
        *base_ptr = '\0'; /* Strip off the '.' by putting ending char */
    }

    /* Extract the root file name */
    base_ptr = strrchr(root_filename, '_');
    if (base_ptr)
    {
        *base_ptr = '\0'; /* Strip off the '_' by putting ending char */
    }

    /* Process the angles for each band */
    for (band_index = 0; band_index < IAS_MAX_NBANDS; band_index++)
    {
        int line;                       /* Line index */
        int samp;                       /* Sample index */
        int num_lines;                  /* Lines in output angle band */
        int num_samps;                  /* Samps in output angle band */
        int index;                      /* Current sample index*/
        size_t angle_size;              /* Malloc angle size */
        short *sat_zn = NULL;           /* Satellite zenith angle */
        short *sat_az = NULL;           /* Satellite azimuth angle */
        short *sun_zn = NULL;           /* Solar zenith angle */
        short *sun_az = NULL;           /* Solar azimuth angle */
        ANGLES_FRAME frame;             /* Output image frame info */
        IAS_MISC_LINE_EXTENT *trim_lut; /* Image trim lookup table */
        int band_number;                /* Band number */ 
                
        /* Retrieve the band number for current index */
        band_number = ias_sat_attr_convert_band_index_to_number(band_index);
        if (band_number == ERROR)
        {
            IAS_LOG_ERROR("Getting band number for band index %d", band_index);
            ias_angle_gen_free(&metadata);
            return EXIT_FAILURE;
        }

        /* Check if this band should be processed */
        if (!parameters.process_band[band_index])
            continue;
        
        /* Get framing information for this band if return is not successful
           then band is not present in metadata so continue */
        if (get_frame(&metadata, band_index, &frame) != SUCCESS)
        {
            IAS_LOG_WARNING("Band not present in metadata for band number %d",
                band_number);
            continue;
        }

        /* Calculate size of subsampled output image */
        num_lines = (frame.num_lines - 1) / sub_sample + 1;
        num_samps = (frame.num_samps - 1) / sub_sample + 1;
        IAS_LOG_INFO("Processing band number %d with %d lines and %d samples "
            "using %d as sub-sampling factor", band_number, num_lines,
            num_samps, sub_sample);

        /* Calculate the angle sizes */
        angle_size = num_lines * num_samps * sizeof(short);

        /* Allocate the satellite buffers if needed */
        if (calc_sat_angles)
        {
            sat_zn = (short *)malloc(angle_size);
            if (!sat_zn)
            {
                IAS_LOG_ERROR("Allocating satellite zenith angle array for "
                    "band number %d", band_number);
                ias_angle_gen_free(&metadata);
                return EXIT_FAILURE;
            }

            sat_az = (short *)malloc(angle_size);
            if (!sat_az)
            {
                IAS_LOG_ERROR("Allocating satellite azimuth angle array for "
                    "band number %d", band_number);
                free(sat_zn);
                ias_angle_gen_free(&metadata);
                return EXIT_FAILURE;
            }
        }

        /* Allocate the solar buffers if needed */
        if (calc_solar_angles)
        {
            sun_zn = (short *)malloc(angle_size);
            if (!sun_zn)
            {
                IAS_LOG_ERROR("Allocating solar angle array for band "
                    "number %d", band_number);
                free(sat_zn);
                free(sat_az);
                ias_angle_gen_free(&metadata);
                return EXIT_FAILURE;
            }

            sun_az = (short *)malloc(angle_size);
            if (!sun_az)
            {
                IAS_LOG_ERROR("Allocating solar azimuth angle array for band "
                    "number %d", band_number);
                free(sat_zn);
                free(sat_az);
                free(sun_zn);
                ias_angle_gen_free(&metadata);
                return EXIT_FAILURE;
            }
        }

        /* Retrieve the trim look up table to remove the scene crenulation */
        trim_lut = ias_misc_create_output_image_trim_lut(
            get_active_lines(&metadata, band_index), 
            get_active_samples(&metadata, band_index),
            frame.num_lines, frame.num_samps);
        if (!trim_lut)
        {
            IAS_LOG_ERROR("Creating the scene trim lookup table for band "
                "number %d", band_number);
            free(sat_zn);
            free(sat_az);
            free(sun_zn);
            free(sun_az);
            ias_angle_gen_free(&metadata);
            return EXIT_FAILURE;
        }

        /* Loop through the L1T lines and samples */
        for (line = 0, index = 0; line < frame.num_lines; line += sub_sample)
        {
            double sun_angles[2];   /* Solar angles */
            double sat_angles[2];   /* Viewing angles */

            for (samp = 0; samp < frame.num_samps; samp += sub_sample, index++)
            {
                if (samp <= trim_lut[line].start_sample || 
                    samp >= trim_lut[line].end_sample)
                {
                    if (calc_sat_angles)
                    {
                        sat_zn[index] = parameters.background;
                        sat_az[index] = parameters.background;
                    }
                    if (calc_solar_angles)
                    {
                        sun_zn[index] = parameters.background;
                        sun_az[index] = parameters.background;
                    }
                    continue;
                }

                /* Calculate the satellite and solar azimuth and zenith */
                if(calculate_angles(&metadata, line, samp, band_index, 
                    parameters.angle_type, sat_angles, sun_angles) != SUCCESS)
                {
                    IAS_LOG_ERROR("Evaluating angles in band %d", band_number);
                    free(sat_zn);
                    free(sat_az);
                    free(sun_zn);
                    free(sun_az);
                    free(trim_lut);
                    ias_angle_gen_free(&metadata);
                    return EXIT_FAILURE;
                }

                /* Quantize the angles by converting from radians to degrees 
                   and scaling by a factor of 100 so it can be stored in the
                   short integer image */
                if (calc_sat_angles)
                {
                    sat_az[index] = (short)round(r2d 
                        * sat_angles[IAS_ANGLE_GEN_AZIMUTH_INDEX]);
                    sat_zn[index] = (short)round(r2d 
                        * sat_angles[IAS_ANGLE_GEN_ZENITH_INDEX]);
                }
                if (calc_solar_angles)  
                {
                    sun_az[index] = (short)round(r2d 
                        * sun_angles[IAS_ANGLE_GEN_AZIMUTH_INDEX]);
                    sun_zn[index] = (short)round(r2d 
                        * sun_angles[IAS_ANGLE_GEN_ZENITH_INDEX]);
                }
            }
        }

        /* Free the lookup table */
        free(trim_lut);
        trim_lut = NULL;

        frame.num_lines = num_lines;
        frame.num_samps = num_samps;
        frame.pixel_size = frame.pixel_size * sub_sample;

        /* Write out the satellite image file */
        if (calc_sat_angles)
        {
            if (ias_angle_gen_write_image(root_filename, sat_az, sat_zn,
                IAS_ANGLE_GEN_SATELLITE, band_index, frame.num_lines, 
                frame.num_samps, frame.ul_corner, frame.pixel_size, 
                &frame.projection) != SUCCESS)
            {
                IAS_LOG_ERROR("Writing satellite angles for band number %d",
                    band_number);
                free(sat_zn);
                free(sat_az);
                free(sun_zn);
                free(sun_az);
                ias_angle_gen_free(&metadata);
                return EXIT_FAILURE;
            }

            free(sat_zn);
            free(sat_az);
            sat_zn = NULL;
            sat_az = NULL;
        }

        /* Write out the solar image file */
        if (calc_solar_angles)
        {
            if (ias_angle_gen_write_image(root_filename, sun_az, sun_zn,  
                IAS_ANGLE_GEN_SOLAR, band_index, frame.num_lines, 
                frame.num_samps, frame.ul_corner, frame.pixel_size, 
                &frame.projection) != SUCCESS)
            {
                IAS_LOG_ERROR("Writing solar angles for band number %d",
                    band_number);
                free(sun_zn);
                free(sun_az);
                ias_angle_gen_free(&metadata);
                return EXIT_FAILURE;
            }

            free(sun_zn);
            free(sun_az);
            sun_zn = NULL;  
            sun_az = NULL;
        }
    }

    /* Release the metadata */
    ias_angle_gen_free(&metadata);

    return EXIT_SUCCESS;
}

/******************************************************************************
NAME: process_parameters

PURPOSE: This fuction processes the command line parameters.

RETURN VALUE: Type = int
    Value     Description
    -----     -----------
    SUCCESS   The parameters were successfully processed
    ERROR     An error occurred processing the parameters
******************************************************************************/
int process_parameters
(
    const int argument_count,        /* I: Arguments count */
    char *arguments[],         /* I: Arguments pointer */
    L8_ANGLES_PARAMETERS *parameters /* O: Generation parameters */
)
{
    char local_angle_type[10];/* Angle type string read in */
    int local_band_count;     /* Number of bands read */
    int local_band_numbers[IAS_MAX_NBANDS] = {ERROR}; /* Read in band numbers */
    int index;                /* Generic loop variable */
    int band_list_init_value; /* Band list initialize value */
    int band_already_in_list[IAS_MAX_NBANDS];
    long local_background; /* Local Background variable */
    int status;

    /* Copy the angle coefficient filename */
    status = snprintf(parameters->metadata_filename, 
        sizeof(parameters->metadata_filename), "%s", arguments[1]); 
    if (status < 0 || status >= sizeof(parameters->metadata_filename))
    {
        IAS_LOG_ERROR("Copying the metadata filename");
        return ERROR;
    }   

    /* Copy the angle type */
    status = snprintf(local_angle_type, sizeof(local_angle_type), "%s", 
        arguments[2]); 
    if (status < 0 || status >= sizeof(local_angle_type))
    {
        IAS_LOG_ERROR("Copying the angle type");
        return ERROR;
    }  

    /* Convert the angle type from string to enumeration */
    ias_misc_convert_to_uppercase(local_angle_type);
    if (strcmp(local_angle_type, "BOTH") == 0)
    {
        parameters->angle_type = AT_BOTH;
    }
    else if (strcmp(local_angle_type, "SATELLITE") == 0)
    {
        parameters->angle_type = AT_SATELLITE;
    }
    else if (strcmp(local_angle_type, "SOLAR") == 0)
    {
        parameters->angle_type = AT_SOLAR;
    }
    else
    {
        IAS_LOG_ERROR("Invalid angle type specified in the input parameters "
            "file.  Valid values are BOTH, SATELLITE, and SOLAR");
        return ERROR;
    }

    /* Parse the sub sample factor */
    if (sscanf(arguments[3], "%d", &parameters->sub_sample_factor) != 1)
    {
        IAS_LOG_ERROR("Parsing the sub_sample_factor from argument %s",
            arguments[3]);
        return ERROR;
    }

    band_list_init_value = 1;
    local_band_count = 0;
    parameters->background = 0;

    
    /* Parse the fill pixel value and band list */
    for (index = 4; index < argument_count; index++)
    {
        if (strcmp("-f", arguments[index]) == 0)
        {
            index++; /* Shift to the fill argument */

            /* Check that a parameter was provided after option handler */
            if (index >= argument_count)
            {
                IAS_LOG_ERROR("No argument specified for fill pixel");
                return ERROR;
            }

            if (sscanf(arguments[index], "%ld", &local_background) != 1)
            {
                IAS_LOG_ERROR("Parsing the background fill pixel from "
                    "argument %s", arguments[index]);
                return ERROR;
            }
            if (local_background < SHRT_MIN || local_background > SHRT_MAX)
            {
                IAS_LOG_ERROR("Background parameter not in valid range of "
                    "(%hi, %hi)", SHRT_MIN, SHRT_MAX);
                return ERROR;
            }
            parameters->background = (short)local_background;
        }
        /* Check if the band list is passed as an argument */
        else if (strcmp("-b", arguments[index]) == 0)
        {
            int offset = 0;
            int character_count = 0;
            const char *bands; 
            int band_length;
            
            /* Shift to the bands argument */
            index++;            

            /* Check that a parameter was provided after option handler */
            if (index >= argument_count)
            {
                IAS_LOG_ERROR("No argument specified for band numbers");
                return ERROR;
            }
            bands = arguments[index];
            band_length = strlen(arguments[index]);

            /* Read in the band list */
            while (character_count <= band_length)
            {
                if (sscanf(bands, "%d%n", &local_band_numbers[local_band_count],
                    &offset) < 1)
                {
                    IAS_LOG_ERROR("Parsing the band list");
                    return ERROR;
                }

                character_count += offset + 1;
                bands += offset + 1; /* Move pointer ahead the number length 
                                        + 1 (for the comma) */
                local_band_count++;
                if (local_band_count > IAS_MAX_NBANDS)
                {
                    IAS_LOG_ERROR("Too many bands specified");
                    return ERROR;   
                }
            }

            band_list_init_value = 0;
        }
        else
        {
            IAS_LOG_ERROR("Unrecognized option %s", arguments[index]);
            return ERROR;
        }
    }

    /* Initialize band list */
    for (index = 0; index < IAS_MAX_NBANDS; index++)
    {
        parameters->process_band[index] = band_list_init_value;
        band_already_in_list[index] = band_list_init_value;
    }

    for (index = 0; index < local_band_count; index++)
    {
        /* Convert the band number entered into a band index */
        int temp_band_index = ias_sat_attr_convert_band_number_to_index(    
                                local_band_numbers[index]);
        if (temp_band_index == ERROR)
        {
            /* No bands numbers entered, so it is an error */
            IAS_LOG_ERROR("Illegal band number %d", local_band_numbers[index]);
            return ERROR;
        }

        /* Make sure the band isn't already in the list */
        if (temp_band_index != ERROR)
        {
            if (!band_already_in_list[temp_band_index])
            {
                /* Put band in list and flag that it has been added 
                   to the list */
                parameters->process_band[temp_band_index] = 1;
                band_already_in_list[temp_band_index] += 1;
            }
            else if (band_already_in_list[temp_band_index] == 1)
            {
                 /* Only print the warning message once for each band */
                 band_already_in_list[temp_band_index] += 1;

                 IAS_LOG_WARNING("Warning: Band number %d appears multiple "
                    "times in the BAND_LIST, only the first one was used",
                    local_band_numbers[index]);
            }
        }
    }

    return SUCCESS;
}
