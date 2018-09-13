/******************************************************************************
NAME: ias_math_find_line_segment_intersection

PURPOSE: Find the intersection of two lines segments.

RETURNS:
    IAS_LINES_INTERSECT, IAS_LINES_PARALLEL, or IAS_LINES_DO_NOT_INTERSECT

******************************************************************************/
#include <math.h>
#include "ias_math.h"

/* Tolerance on whether two lines are parallel */
#define TOL 1.0e-8

int ias_math_find_line_segment_intersection
(
    const IAS_LINE_SEGMENT *line1,  /* I: First line segment */
    const IAS_LINE_SEGMENT *line2,  /* I: Second line segment */
    double *intersect_x,            /* O: X coordinate of intersection */
    double *intersect_y             /* O: Y coordinate of intersection */
)
{
    double line1_x;
    double line1_y;
    double line2_x;
    double line2_y;
    double determinant;

    /* Calculate the differences in X and Y directions for each line segment */
    line1_x = line1->x1 - line1->x2;
    line1_y = line1->y1 - line1->y2;
    line2_x = line2->x1 - line2->x2;
    line2_y = line2->y1 - line2->y2;

    determinant = line2_x * line1_y - line2_y * line1_x;

    if (fabs(determinant) > TOL)
    {
        double x2_diff;
        double y2_diff;
        double s;
        double t;

        x2_diff = line2->x2 - line1->x2;
        y2_diff = line2->y2 - line1->y2;

        /* If lines are not parallel then see if they intersect */
        s = (line2_x * y2_diff - line2_y * x2_diff) / determinant;
        t = (line1_x * y2_diff - line1_y * x2_diff) / determinant;

        if (s < 0.0 || s > 1.0 || t < 0.0 || t > 1.0)
        {
            /* Lines do not intersect */
            return IAS_LINES_DO_NOT_INTERSECT;
        }
        else
        {
            /* Lines intersect so calculate the intersection point */
            *intersect_x = line1->x2 + line1_x * s;
            *intersect_y = line1->y2 + line1_y * s;

            return IAS_LINES_INTERSECT;
        }
    }
    else
    {
        /* Lines are parallel */
        return IAS_LINES_PARALLEL;
    }
}
