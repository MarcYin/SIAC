#include <stdlib.h>
#include "ias_logging.h"
#include "ias_miscellaneous.h"

/******************************************************************************
NAME: ias_misc_create_output_image_trim_lut

PURPOSE: Builds a lookup table that contains the start and end sample of the
    area to trim the output image to eliminate the extra SCA data at the top
    and bottom of the imagery.

RETURNS:
    Pointer to the lookup table created (same number of entries as output lines
    in the grid) or NULL if an error occurred.

ALGORITHM REFERENCES:
    Resample ADD

******************************************************************************/
IAS_MISC_LINE_EXTENT *ias_misc_create_output_image_trim_lut
(
    const double *line,            /* I: Line trim box(ul, ur, lr, ll) */
    const double *samp,            /* I: Sample trim box(ul, ur, lr, ll) */
    int output_lines,              /* I: Number of lines */
    int output_samples             /* I: Number of samples */
)
{
    int line_index;                /* line loop index */
    int index;                     /* Counter */
    double start_sample;           /* Start sample for the output line */
    double end_sample;             /* End sample for the output line */
    IAS_MISC_LINE_EXTENT *trim_lut;/* Trimming lookup table */
    IAS_LINE_SEGMENT output_line;  /* Output line segment */
    IAS_LINE_SEGMENT trim_lines[4];/* 4 line segments making up the trimming
                                      bounding box */

    /* Construct the 4 line segments for the trimming box */
    for (index = 0; index < 4; index++)
    {
        trim_lines[index].x1 = samp[index];
        trim_lines[index].x2 = samp[(index + 1) % 4];
        trim_lines[index].y1 = line[index];
        trim_lines[index].y2 = line[(index + 1) % 4];
    }

    /* Make sure "search" line start and end samples extend beyond the
       bounding area limits */
    start_sample = -1.0;
    end_sample = output_samples;
    for (index = 0; index < 4; index ++)
    {
        if (samp[index] < start_sample)
            start_sample = samp[index] - 1.0;
        if (samp[index] > end_sample)
            end_sample = samp[index] + 1.0;
    }

    /* Set the output line X values */
    output_line.x1 = start_sample;
    output_line.x2 = end_sample;

    /* Allocate the look-up table that will define the starting and ending
       sample locations for imagery in output file.  Sample locations are
       defined for each line in output. */
    trim_lut = (IAS_MISC_LINE_EXTENT *)malloc(sizeof(*trim_lut) * output_lines);
    if (!trim_lut)
    { 
        IAS_LOG_ERROR("Error allocating trimming buffer");
        return NULL;
    }

    /* Initialize the lookup table */
    for (line_index = 0; line_index < output_lines; line_index++)
    {
        trim_lut[line_index].start_sample = 0;
        trim_lut[line_index].end_sample = 0;
    }

    /* Check each line of the output image to see if and where it intersects
       the trimming area */
    for (line_index = 0; line_index < output_lines; line_index++)
    {
        double intersect_line[4];   /* line intersection points */
        double intersect_sample[4]; /* sample intersection points */
        int intersect_count = 0;

        output_line.y1 = line_index;
        output_line.y2 = line_index;

        /* Check the output line against each of the lines in the trimming
           area */
        for (index = 0; index < 4; index++)
        {
            int status;                 /* Intersection status */

            /* Check for intersection of output line and the current line
               segment of the trimming box */ 
            status = ias_math_find_line_segment_intersection(&output_line,
                &trim_lines[index], &intersect_sample[intersect_count],
                &intersect_line[intersect_count]);
            if (status == IAS_LINES_INTERSECT)
            {
                /* Intersection found */
                intersect_count++;
            }
        }

        /* If there are 2 or more intersections, find the minimum and maximum
           values since they are the ones that should be used for trimming the
           image.  It should be pretty rare to get more than 2 intersections,
           but it can happen if one of the corners falls exactly on an output
           image line.  If there are zero or one intersections, the lookup 
           table entry can remain at the original zero values it is initialized
           to.  Note that one intersection could happen if the corner of the
           bounding box falls right on the output line. */
        if (intersect_count >= 2)
        {
            double min_samp = intersect_sample[0];
            double max_samp = intersect_sample[0];
            int index;

            for (index = 1; index < intersect_count; index++)
            {
                if (intersect_sample[index] < min_samp)
                    min_samp = intersect_sample[index];
                if (intersect_sample[index] > max_samp)
                    max_samp = intersect_sample[index];
            }

            trim_lut[line_index].start_sample = min_samp;
            trim_lut[line_index].end_sample   = max_samp;
        }
    }

    return trim_lut;
}
