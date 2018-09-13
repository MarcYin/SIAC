/******************************************************************************

UNIT NAME: ias_odl_free_tree.c

PURPOSE: Close the ODL tree.

******************************************************************************/
#include "lablib3.h"                 /* prototypes for the ODL functions */
#include "ias_odl.h"

void ias_odl_free_tree
(
    IAS_OBJ_DESC *p_lp           /* I: Allocated memory to free */
)
{
    OdlFreeTree(p_lp);
}
