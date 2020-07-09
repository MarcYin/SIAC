/*************************************************************************
NAME:    ias_parm_map_odl_type

PURPOSE: Maps from an IAS type to an ODL type

RETURNS:
Type=int
Value    Description
-----    -----------
SUCCESS  All required items successfully read from the ODL file
ERROR    error locating or reading ODL file

**************************************************************************/

#include "ias_parm.h"
#include "ias_parm_private.h"
#include "ias_const.h"
#include "ias_logging.h"
#include "ias_odl.h"

/* Set up the arrays to map IAS type with ODL type */
static const IAS_ODL_TYPE ias_type_to_odl_type[4] =
        {IAS_ODL_Int, IAS_ODL_Int, IAS_ODL_Double, IAS_ODL_String};
static const IAS_ODL_TYPE ias_type_to_odl_type_array[4] =
        {IAS_ODL_Int, IAS_ODL_Int, IAS_ODL_Double, IAS_ODL_ArrayOfString};

int ias_parm_map_odl_type
(
    IAS_PARM_PARAMETER_DEFINITION *pd, /* I: pointer to current parameter
                                             definition */
    IAS_ODL_TYPE *type                 /* O: ODL type to map to */
)
{
    if (pd->is_an_array)
    {
        /* Verify that the type index is within the array limit
           when mapping from ias type to ODL type */
        if ((pd->type + 1) > (sizeof(ias_type_to_odl_type_array)
                                            / sizeof(IAS_ODL_TYPE)))
        {
            IAS_LOG_ERROR("Mapping from IAS type to ODL type");
            return ERROR;
        }
        *type = ias_type_to_odl_type_array[pd->type];
    }
    else
    {
        /* Verify that the type index is within the array limit
               when mapping from ias type to ODL type */
    if ((pd->type + 1) > (sizeof(ias_type_to_odl_type)
                                            / sizeof(IAS_ODL_TYPE)))
        {
            IAS_LOG_ERROR("Mapping from IAS type to ODL type");
            return ERROR;
        }
        *type = ias_type_to_odl_type[pd->type];
    }

    return SUCCESS;
}
