/*************************************************************************
NAME:    ias_parm_read

PURPOSE: General purpose routine to read a list of items from an ODL file.

RETURNS:
Type=int
Value    Description
-----    -----------
SUCCESS  All required items successfully read from the ODL file
ERROR    error locating or reading ODL file

**************************************************************************/

#include <stdio.h>
#include <string.h>
#include "ias_parm.h"
#include "ias_parm_private.h"
#include "ias_const.h"
#include "ias_logging.h"
#include "ias_odl.h"

int ias_parm_read
(
    const char *odl_file_name,   /* I: ODL file name to read */
    IAS_PARM_PARAMETER_DEFINITION **odl_list_ptr,
                                 /* I/O : pointer to list of items to read
                                    from the ODL file */
    int list_length              /* I: number of items in the list */
)
{
    int status;                  /* function return status */
    IAS_OBJ_DESC *odl_tree = NULL; /* parsed ODL object */
    int list_index;              /* current list index */
    IAS_PARM_PARAMETER_DEFINITION *pd = NULL;
                                 /* pointer to current item in the list */
    IAS_ODL_TYPE type;           /* this is necessary as current *pd item 
                                    could be an array */
    int range_err_flag = FALSE;  /* initiate range flag to false */
    int namelength = 0;          /* Length of parameter name */
    int size_of_array = 0;       /* number of elements defined in the array */


    /* read the ODL file */
    odl_tree = ias_odl_read_tree(odl_file_name);
    if (odl_tree == NULL)
    {
        IAS_LOG_ERROR("Reading parameters from %s", odl_file_name);
        return ERROR;
    }

    /* loop through the list of parameters to read */
    for (list_index = 0; list_index < list_length; list_index++)
    {
        /* get a pointer to the current item */
        pd = odl_list_ptr[list_index];

        /* Check to make sure the parameter name is less than the maximum
           number of allowed characters for the database (40) */
        namelength = strlen(pd->name);
        if (namelength > IAS_PARM_MAX_DB_NAMELENGTH)
        {
            IAS_LOG_ERROR("Parameter name %s has length %d that is "
                          "longer than the maximum allowed database "
                          "length %d", pd->name, namelength,
                          IAS_PARM_MAX_DB_NAMELENGTH);
            ias_odl_free_tree(odl_tree);
            return ERROR;
        }


        /* map the IAS type to ODL type */
        status = ias_parm_map_odl_type(pd, &type);
        if (status != SUCCESS)
        {
            IAS_LOG_ERROR("Mapping type for ODL parameter \"%s\"", pd->name);
            ias_odl_free_tree(odl_tree);
            return ERROR;
        }

        status = ias_odl_get_field(pd->value.type_void, pd->value_bytes,
            type, odl_tree, NULL, pd->name, &pd->count_read);

        if ((status != SUCCESS) || (pd->count_read < pd->min_count))
        {
            /* an error occurred reading the ODL file or not enough items
               were read.  If it was a required parameter, generate
               an error message. 
            */
            if (pd->is_required == IAS_PARM_REQUIRED)
            {
                /* make sure it was really an error */
                if (status != SUCCESS)
                {
                    /* error reading required item, log an error */
                    IAS_LOG_ERROR("Reading ODL tree");
                }
                /* not an error reading the tree, but parameter was required
                   and was not found */
                IAS_LOG_ERROR("Reading required ODL parameter \"%s\"", 
                    pd->name);
                /* free the ODL tree */
                ias_odl_free_tree(odl_tree);
                return ERROR;
            } /* End parameter is required */
            else if (status == IAS_ODL_NOT_FOUND)
            {
                /* an optional parameter wasn't found and there are default
                   values available so use them */

                /* check to make sure defaults are not greater than size of 
                   storage */
                if ((pd->type == IAS_PARM_INT) ||
                    (pd->type == IAS_PARM_BOOLEAN))
                {
                    int bytes_needed = sizeof(int) * pd->number_defaults;
                    if (bytes_needed > pd->value_bytes)
                    {
                        IAS_LOG_ERROR("Parameter values read exceed "
                            "storage available for ODL parameter \"%s\"", 
                            pd->name);
                        ias_odl_free_tree(odl_tree);
                        return ERROR;
                    }
                }
                else if (pd->type == IAS_PARM_DOUBLE)
                {
                    int bytes_needed = sizeof(double) * pd->number_defaults;
                    if (bytes_needed > pd->value_bytes)
                    {
                        IAS_LOG_ERROR("Parameter values read exceed "
                            "storage available for ODL parameter \"%s\"", 
                            pd->name);
                        ias_odl_free_tree(odl_tree);
                        return ERROR;
                    }
                }
                else if (pd->type == IAS_PARM_STRING)
                {
                    int bytes_needed = sizeof(char) 
                        * pd->number_defaults;
                    if (bytes_needed > pd->value_bytes)
                    {
                        IAS_LOG_ERROR("Parameter values read exceed "
                            "storage available for ODL parameter \"%s\"", 
                            pd->name);
                        ias_odl_free_tree(odl_tree);
                        return ERROR;
                    }
                }

                if (pd->number_defaults > 0)
                {
                    const IAS_PARM_TYPE_UNION_CONST *dv = &pd->default_values;
                    int i; /* counter for default value loop */

                    if ((pd->type == IAS_PARM_INT) || 
                        (pd->type == IAS_PARM_BOOLEAN))
                    {
                        for (i = 0; i < pd->number_defaults; i++)
                        {
                            pd->value.type_int[i] = dv->type_int[i];
                        }
                        pd->count_read = pd->number_defaults;
                    }
                    else if (pd->type == IAS_PARM_DOUBLE)
                    {
                        for (i = 0; i < pd->number_defaults; i++)
                        {
                            pd->value.type_double[i] = dv->type_double[i];
                        }
                        pd->count_read = pd->number_defaults;
                    }
                    else if (pd->type == IAS_PARM_STRING)
                    {
                        if (!pd->is_an_array) /* is not an array */
                        {
                            strcpy(pd->value.type_string, 
                                *dv->type_string_array);
                            pd->count_read = 1;
                        }
                        else /* is an array */
                        {
                            /* Get the number of elements in the array */
                            size_of_array = pd->number_defaults;
                            for (i = 0; i < size_of_array; i++)
                            {
                                pd->value.type_string_array[i] = strdup(
                                    dv->type_string_array[i]);
                            }
                            pd->count_read = pd->number_defaults;
                        }
                    }
                }   /* End defaults condition */
            } /* End parameter is not required */
            else
            {
                /* an error happened that needs to be reported */
                IAS_LOG_ERROR("Reading ODL parameter \"%s\"", pd->name);
                ias_odl_free_tree(odl_tree);
                return ERROR;
            }
        } /* End read was not successful and check required flag */

        /* Call range check to ensure the value is valid */
        status = ias_parm_check_ranges(pd->value, pd->valid_values,
                        pd->count_read, pd->number_valids,
                        pd->type, pd->is_an_array);
        if (status != SUCCESS)
        {
            IAS_LOG_ERROR("Parameter %s is invalid",pd->name);
            ias_odl_free_tree(odl_tree);
            return ERROR;
        }
    } /* End for loop through the parameters */

    /* free the ODL tree */
    ias_odl_free_tree(odl_tree);

    /* If range error flag is set to true then return error */
    if (range_err_flag == TRUE)
    {
        IAS_LOG_ERROR("An error occurred testing parameter validity");
        return ERROR;
    }

    return SUCCESS;
}

