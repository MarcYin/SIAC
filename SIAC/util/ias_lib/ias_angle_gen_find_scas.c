/* IAS Library Includes */
#include "ias_angle_gen_private.h"

/* Local Defines */
#define SCA_OVERLAP 50 /* Max number of overlapping pixels */

/*******************************************************************************
Name: ias_angle_gen_find_scas

Purpose: Uses the L1T to L1R rational polynomials to determine which SCA, or
         SCAs, the input L1T line/sample/height location falls in, and 
         returns the L1R line/sample coordinates associated with each valid 
         SCA

Note: If height not needed in find pass NULL pointer as height. Also the
      l1r_line and l1r_samp pointers need to have space to for 2 SCA 
      line/sample combinations.

Return: 
    Type = integer
    On success: Number of SCAs found
    On error: N/A
 ******************************************************************************/
int ias_angle_gen_find_scas
(
    const IAS_ANGLE_GEN_BAND *metadata,/* I: Metadata for current band */
    double l1t_line,      /* I: Input L1T line */
    double l1t_samp,      /* I: Input L1T sample */
    const double *height, /* I: Input height, NULL for zero height */
    double *l1r_line,     /* O: Array of output L1R line numbers */
    double *l1r_samp      /* O: Array of output L1R sample numbers */
) 
{
    int sca_index;   /* SCA index */
    int nsca_found;  /* Number of SCAs containing input point */
    int scas_tested; /* Number of SCAs tested */

    /* Initialize the number found */
    nsca_found = 0;
    scas_tested = 0;
    sca_index = metadata->num_scas / 2;

    /* Compute the location for this SCA */
    while (sca_index >= 0 && sca_index < metadata->num_scas)
    {
        double line_offset;          /* Offset value of L1T line */
        double samp_offset;          /* Offset value of L1T sample */
        double height_offset;        /* Offset value of height */
        double line;                 /* Local L1R line */
        double sample;               /* Local L1R sample */
        const IAS_ANGLE_GEN_IMAGE_RPC_TERMS *line_terms;/* Line terms pointer */
        const IAS_ANGLE_GEN_IMAGE_RPC_TERMS *samp_terms;/* Samp terms pointer */

        line_terms = &metadata->sca_metadata[sca_index].line_terms;
        samp_terms = &metadata->sca_metadata[sca_index].samp_terms;
        line_offset = l1t_line - line_terms->l1t_mean_offset;
        samp_offset = l1t_samp - samp_terms->l1t_mean_offset;

        /* Check if height should be used */
        if (height)
        {    
            height_offset = *height 
                - metadata->sca_metadata[sca_index].mean_height;
        }
        else
        {
            height_offset = 0;
        }

        /* Calculate the l1r line and sample location */
        line = (line_terms->numerator[0] + line_terms->numerator[1] 
            * line_offset + line_terms->numerator[2] * samp_offset 
            + line_terms->numerator[3] * height_offset 
            + line_terms->numerator[4] * line_offset * samp_offset) 
            / (1.0 + line_terms->denominator[0] * line_offset 
            + line_terms->denominator[1] * samp_offset 
            + line_terms->denominator[2] * height_offset 
            + line_terms->denominator[3] * line_offset * samp_offset) 
            + line_terms->l1r_mean_offset;

        sample = (samp_terms->numerator[0] + samp_terms->numerator[1] 
            * line_offset + samp_terms->numerator[2] * samp_offset
            + samp_terms->numerator[3] * height_offset 
            + samp_terms->numerator[4] * line_offset * samp_offset) 
            / (1.0 + samp_terms->denominator[0] * line_offset 
            + samp_terms->denominator[1] * samp_offset 
            + samp_terms->denominator[2] * height_offset 
            + samp_terms->denominator[3] * line_offset * samp_offset) 
            + samp_terms->l1r_mean_offset;

        /* See if we are in the right place */
        if (sample >= 0 && sample <= (metadata->l1r_samps - 1))
        {
            scas_tested++;
            if (line >= 0 && line < metadata->l1r_lines)
            {
                l1r_line[nsca_found] = line;
                l1r_samp[nsca_found] = sample + sca_index * metadata->l1r_samps;
                nsca_found++;
            }

            /* Return if more than 1 SCA tested */
            if (scas_tested > 1) 
            {
                return nsca_found;
            }

            if (sample < SCA_OVERLAP)
            { 
                sca_index--;
            }
            else if (sample > (metadata->l1r_samps - SCA_OVERLAP))
            { 
                sca_index++;
            }
            else 
            {
                return nsca_found;
            }
        }
        else
        {
            /* Return if a SCA has been tested */
            if ((scas_tested > 0) || (sca_index == 0 && sample < 0) || 
                (sca_index == (metadata->num_scas - 1) && sample 
                > (metadata->l1r_samps - 1)))
            {
                return nsca_found;
            }
            
            if (sample < 0)
            {
                sca_index = (sample + sca_index * metadata->l1r_samps)
                    / metadata->l1r_samps;
            }
            else
            {
                sca_index = ((sample + sca_index * metadata->l1r_samps) + 1) 
                    / metadata->l1r_samps;
            }

            if (sca_index < 0)
            { 
                sca_index = 0;
            }

            if (sca_index >= metadata->num_scas)
            { 
                sca_index = metadata->num_scas - 1;
            }
        }
    }

    return nsca_found;
}
