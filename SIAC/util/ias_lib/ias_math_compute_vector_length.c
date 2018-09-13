/*******************************************************************************
NAME: ias_math_compute_vector_length

PURPOSE: Find the length of a vector

RETURN VALUE
TYPE = double
Calculated length of the vector

*******************************************************************************/
#include <math.h>
#include "ias_math.h"

double ias_math_compute_vector_length
(
    const IAS_VECTOR *vec        /* I Vector to find the length of */
)
{

    /* Calculate the norm */
    return (sqrt(vec->x * vec->x + vec->y * vec->y + vec->z * vec->z));
}
