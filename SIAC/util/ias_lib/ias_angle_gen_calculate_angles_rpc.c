/* Standard Library Includes */
#include <math.h>

/* IAS Library Includes */
#include "ias_logging.h"
#include "ias_math.h"
#include "ias_angle_gen_private.h"

/*******************************************************************************
Name: calculate_rpc_vector_value

Purpose: Calculates the individual vector value used to evaluate the
         rational polynomial coefficient.
 
Return:
    Type = void
 ******************************************************************************/
static void calculate_rpc_vector_value
(
    double l1t_line,           /* I: L1T line coordinate */
    double l1t_samp,           /* I: L1T sample coordinate */
    double l1r_line,           /* I: L1R line coordinate */
    double l1r_samp,           /* I: L1R sample coordinate */
    double mean_offset,        /* I: Vector mean offset */
    double height,             /* I: Input height */
    const double *numerator,   /* I: Numerator values pointer */
    const double *denominator, /* I: Denominator values pointer */
    double *output_value       /* O: Output vector value */ 
)
{
    double equation_num;    /* Equation numerator */
    double equation_den;    /* Equation denominator */

    /* Calculate the numerator */
    equation_num = numerator[0] + numerator[1] * l1t_line
        + numerator[2] * l1t_samp
        + numerator[3] * height
        + numerator[4] * l1r_line
        + numerator[5] * l1t_line * l1t_line
        + numerator[6] * l1t_samp * l1t_line 
        + numerator[7] * l1t_samp * l1t_samp
        + numerator[8] * l1r_samp * l1r_line * l1r_line
        + numerator[9] * l1r_line * l1r_line * l1r_line;
       
    /* Calculate the denominator */
    equation_den = 1.0 + denominator[0] * l1t_line
        + denominator[1] * l1t_samp
        + denominator[2] * height
        + denominator[3] * l1r_line 
        + denominator[4] * l1t_line * l1t_line 
        + denominator[5] * l1t_line * l1t_samp
        + denominator[6] * l1t_samp * l1t_samp
        + denominator[7] * l1r_samp * l1r_line * l1r_line
        + denominator[8] * l1r_line * l1r_line * l1r_line;
   
    *output_value = mean_offset + (equation_num / equation_den);
}

/*******************************************************************************
Name: ias_angle_gen_calculate_angle_rpc

Purpose: Calculates the satellite viewing and solar illumination zenith and 
         azimuth angles for a specified L1T line/sample and height (from DEM)
         using the rational polynomial coefficients for the current band. If
         the SCA location falls in the SCA overlap area it returns the
         average value of the combined SCAs.

Return: 
    Type = integer
    SUCCESS / ERROR
 ******************************************************************************/
int ias_angle_gen_calculate_angles_rpc
(
    const IAS_ANGLE_GEN_METADATA *metadata, /* I: Metadata structure */
    double l1t_line,        /* I: Output space line coordinate */
    double l1t_samp,        /* I: Output space sample coordinate */
    const double *elev,     /* I: Pointer to input elevation or NULL if mean
                              scene height should be used*/
    int band_index,         /* I: Current band index */
    IAS_ANGLE_GEN_TYPE sat_or_sun_type,     /* I: Angle calculation type */
    int *outside_image_flag,/* O: Flag indicating return was outside image */
    double *angle           /* O: Array containing zenith and azimuth angles */
)       
{
    int nsca_found;     /* Number of SCAS containing the point */
    int sca_index;      /* SCA index */
    double height;      /* Model height */
    double l1r_line[2]; /* Input space (L1R) line coordinate declared size 2
                           so it can support 2 SCAs */
    double l1r_samp[2]; /* Input space (L1R) sample coordinate declared size
                           2 so it can support 2 SCAs */
    const IAS_ANGLE_GEN_BAND *band_ptr;    /* Pointer to current band */
    const IAS_ANGLE_GEN_ANG_RPC *data_ptr; /* Solar or satellite data pointer */

    /* Check that the band index is valid */
    if (!ias_angle_gen_valid_band_index(metadata, band_index))
    {
        IAS_LOG_ERROR("Band index %d is invalid", band_index);
        return ERROR;
    }

    /* Initialize the output zenith and azimuth*/
    angle[IAS_ANGLE_GEN_ZENITH_INDEX] = 0.0;
    angle[IAS_ANGLE_GEN_AZIMUTH_INDEX] = 0.0;
    *outside_image_flag = 0;

    /* Setup the band pointer */
    band_ptr = &metadata->band_metadata[band_index];

    /* Set the height to use */
    height = band_ptr->satellite.mean_height;
    if (elev) 
    {
        height = *elev;
    }

    /* See which SCA(s) the point falls inside */
    nsca_found = ias_angle_gen_find_scas(band_ptr, l1t_line, l1t_samp, elev,
        l1r_line, l1r_samp);

    if (nsca_found < 1)
    {
        *outside_image_flag = 1; /* Set the flag to indicate outside of image */
        return SUCCESS;         /* Outside active image, return zero vectors */
    }

    if (nsca_found > 2)
    {
        IAS_LOG_ERROR("Too many SCAs found locating point in active image");
        return ERROR;
    }

    /* Offset the output space coordinates */
    l1t_line -= band_ptr->satellite.line_terms.l1t_mean_offset;
    l1t_samp -= band_ptr->satellite.samp_terms.l1t_mean_offset;
    height -= band_ptr->satellite.mean_height;

    data_ptr = &band_ptr->solar;
    if (sat_or_sun_type == IAS_ANGLE_GEN_SATELLITE)
    {
        data_ptr = &band_ptr->satellite;
    }

    /* Check the located SCAs */
    for (sca_index = 0; sca_index < nsca_found; sca_index++)
    {
        IAS_VECTOR unit_vector;     /* Normalized vector */
        IAS_VECTOR vector;          /* Viewing vector */
        double l1r_line_from_offset;/* L1R line location using offset */        
        double l1r_samp_from_offset;/* L1R sample location using offset */

        /* Determine the line and sample location using the L1R offset */   
        l1r_line_from_offset = l1r_line[sca_index] 
            - band_ptr->satellite.line_terms.l1r_mean_offset;
        l1r_samp_from_offset = l1r_samp[sca_index] 
            - band_ptr->satellite.samp_terms.l1r_mean_offset;

        /* Calculate the rpc vector */ 
        calculate_rpc_vector_value(l1t_line, l1t_samp, l1r_line_from_offset, 
            l1r_samp_from_offset, data_ptr->mean_offset.x, height, 
            data_ptr->x_terms.numerator, data_ptr->x_terms.denominator, 
            &vector.x);

        calculate_rpc_vector_value(l1t_line, l1t_samp, l1r_line_from_offset, 
            l1r_samp_from_offset, data_ptr->mean_offset.y, height, 
            data_ptr->y_terms.numerator, data_ptr->y_terms.denominator, 
            &vector.y);

        calculate_rpc_vector_value(l1t_line, l1t_samp, l1r_line_from_offset, 
            l1r_samp_from_offset, data_ptr->mean_offset.z, height, 
            data_ptr->z_terms.numerator, data_ptr->z_terms.denominator, 
            &vector.z);

        /* Normalize the vector in case the polynomial fit results in non-unit
           vectors that will fail the following trig functions */
        if (ias_math_compute_unit_vector(&vector, &unit_vector) != SUCCESS)
        {
            IAS_LOG_ERROR("Unable to normalize the rpc vector");
            return ERROR;
        }

        /* Calculate zenith and azimuth angles */
        angle[IAS_ANGLE_GEN_ZENITH_INDEX] += acos(unit_vector.z);
        angle[IAS_ANGLE_GEN_AZIMUTH_INDEX] += atan2(unit_vector.x, 
            unit_vector.y);
    }

    /* Average the angles */
    angle[IAS_ANGLE_GEN_ZENITH_INDEX] /= nsca_found;
    angle[IAS_ANGLE_GEN_AZIMUTH_INDEX] /= nsca_found;

    return SUCCESS;
}
