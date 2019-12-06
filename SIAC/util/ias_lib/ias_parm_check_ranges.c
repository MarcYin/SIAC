/*************************************************************************

NAME:    ias_parm_check_ranges

PURPOSE: This routine validates value(s) read against the range(s) provided.

RETURNS:
Type=int
Value   Description
-----   -----------
SUCCESS All ranges are valid
ERROR   Value does not valid match range 

Note:   If a parameter is not set or there is not range, this routine 
        will return SUCCESS without any processing and no comments.

**************************************************************************/

#include <stdio.h>
#include <string.h>
#include "ias_const.h"
#include "ias_logging.h"
#include "ias_parm.h"
#include "ias_parm_private.h"

int ias_parm_check_ranges
(
    IAS_PARM_TYPE_UNION value,      /* I: pointer to table of values */
    IAS_PARM_TYPE_UNION_CONST valid_values,
                                    /* I: pointer to table of valid values */
    int count_read,                 /* I: count of values read */
    int valid_count,                /* I: count of valid values */
    IAS_PARM_TYPE parm_type,        /* I: parameter type */
    int array_flag                  /* I: IAS_PARM_ARRAY / IAS_PARM_NOT_ARRAY  
                                          flag for array, necessary for 
                                          reading the string */
)
{
    int index;                      /* current list index */
    int valid_index;                /* current list index for valid checking */
    int low_range;                  /* Low range index for valid list */
    int high_range;                 /* High range index for valid list */
    char **str_ptr;                 /* Temporary pointer to the string value */
    int valid_found = FALSE;        /* Flag to determine if valid is found */

    /* Test if any range or valid checking is necessary */
    if ((valid_count == 0) || (count_read == 0))
    {
        /* No valids to check, return success */
        return SUCCESS;
    }

    /* For each type, check the range or valid list */

    /* Integer and Boolean (which is an Integer with a range of 0 - 1 */
    if (parm_type == IAS_PARM_INT || parm_type == IAS_PARM_BOOLEAN)
    {
        for (index = 0; index < count_read; index++)
        {
            /* If 1 range is set, the count value may be set to 1 indicating
               one pair.  In this case, set the valid count to 2 for 2 values */
            if (valid_count == 1)
                valid_count = 2;

            /* The numerical valid list should be a list of range pairs, an
               odd number would not make sense and will return an error */
            if ((valid_count % 2) != 0)
            {
                IAS_LOG_ERROR("Unexpected valid count %d, "
                                    "expect even number for min/max pairs",
                                    valid_count);
                return ERROR;
            }

            for (valid_index = 0; valid_index < valid_count / 2; valid_index++)
            {
                /* Get the index value to go through the array by twos */
                low_range = valid_index * 2;
                high_range = low_range + 1;

                /* The current database/work order system does not support
                   multiple valid ranges that are not a list of single values 
                   (ie 10 to 10, 30 to 30, 50 to 50).  Return an error if
                   this is not the case. */
                if (valid_count > 2)
                {
                    if (valid_values.type_int[low_range] != 
                        valid_values.type_int[high_range])
                    {
                        IAS_LOG_ERROR("Unsupported mutiple range, "
                                    "only single ranges can be entered");
                        return ERROR;
                    }
                }

                if (value.type_int[index] >= valid_values.type_int[low_range] &&
                    value.type_int[index] <= valid_values.type_int[high_range])
                {
                    valid_found = TRUE;
                    break;
                }

            } /* Valid count end */

            if (valid_found == FALSE)
            {
                IAS_LOG_ERROR("Parameter value %d out of range",
                                            value.type_int[index]);
                return ERROR;
            }

            /* Set valid found flag */
            valid_found = FALSE;
        } /* End loop through the parameter array */
    }
    /* Double */
    else if (parm_type == IAS_PARM_DOUBLE)
    {
        for (index = 0; index < count_read; index++)
        {
            /* If one range is set, the count value may be set to one indicating
            one pair.  In this case, set the valid case to 2 for 2 values */
            if (valid_count == 1)
                valid_count = 2;

            /* The numerical valid list should be a list of range pairs, an
               odd number would not make sense and will return an error */
            if ((valid_count % 2) != 0)
            {
                IAS_LOG_ERROR("Unexpected valid count %d, "
                                    "expect even number for min/max pairs",
                                    valid_count);
                return ERROR;
            }

            for (valid_index = 0; valid_index < valid_count/2; valid_index++)
            {
                /* Get the index value to go through the array by twos */
                low_range = valid_index * 2;
                high_range = low_range + 1;

                /* The current database/work order system does not support
                   multiple valid ranges that are not a list of single values 
                   (ie 10.0 to 10.0, 30.0 to 30.0, 50.0 to 50.0).  
                   Return an error if this is not the case. */
                if (valid_count > 2)
                {
                    if (valid_values.type_double[low_range] != 
                        valid_values.type_double[high_range])
                    {
                        IAS_LOG_ERROR("Unsupported mutiple range, "
                                    "only single ranges can be entered");
                        return ERROR;
                    }
                }

                if (value.type_double[index] >= 
                        valid_values.type_double[low_range] &&
                        value.type_double[index] <=
                        valid_values.type_double[high_range])
                {
                    valid_found = TRUE;
                    break;
                }

            } /* Valid count end */

            if (valid_found == FALSE)
            {
                IAS_LOG_ERROR("Parameter value %lf out of range",
                                                value.type_double[index]);
                return ERROR;
            }
            /* Set valid found flag */
            valid_found = FALSE;
        } /* End loop through the parameter array */
                    
    }

    /* String and String Array */

    else if (parm_type == IAS_PARM_STRING)
    {
        /* parm type is not an array */
        if (!array_flag)
        {
            str_ptr = &value.type_string;
        }
        else /* parm type is array of strings */
        {
            str_ptr = value.type_string_array;
        }
        
        /* Initialize string found flag */
        for (index = 0; index < count_read; index++)
        {
            /* initialize the string found flag */
            valid_found = FALSE;

            /* Test if string is in valid list */
            const char **curr = valid_values.type_string_array;
            for (valid_index = 0; valid_index < valid_count; valid_index++)
            {
                if (strcmp(*str_ptr, *curr) == 0)
                {
                    valid_found = TRUE;
                    break;
                }
                curr++;
            } /* End for loop for valid list */
            /* Check if the string was not found */
            if (valid_found == FALSE)
            {
                IAS_LOG_ERROR("Parameter value %s not in valid list", 
                                                    *str_ptr);
                return ERROR;
            }

            /* Increment the pointer if it is an array of strings */
            if (array_flag)
            {
                str_ptr++;
            }

        } /* End of values list */
    } /* End of String condition */

    return SUCCESS;
}


