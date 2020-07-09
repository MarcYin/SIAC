/*******************************************************************************
NAME: ias_math_compute_unit_vector

PURPOSE: Compute a vector's unit vector

RETURNS: SUCCESS or ERROR

*******************************************************************************/
#include <math.h>
#include "ias_math.h"
#include "ias_logging.h"

int ias_math_compute_unit_vector
(
    const IAS_VECTOR *vec,      /* I: Input vector */
    IAS_VECTOR *unit_vector     /* O: Unit vector of the input vector */
)
{
    double magnitude;       /* Magnitude of the input vector */

    /* Calculate the norm */
    magnitude = ias_math_compute_vector_length(vec);

    if (magnitude == 0.0)
    {
        IAS_LOG_ERROR("Input vector magnitude is zero");
        return ERROR;
    }

    /* Calculate the unit vector */
    unit_vector->x = vec->x / magnitude;
    unit_vector->y = vec->y / magnitude;
    unit_vector->z = vec->z / magnitude;

    return SUCCESS;
}
