/******************************************************************************

UNIT NAME: ias_odl_read_tree.c

PURPOSE: Open and parse the ODL file.

RETURN VALUE:
    Type = OBJDESC *

Value           Description
--------------- ----------------------------------------------------------------
successful        Pointer to ODL tree 
NULL              Failure

******************************************************************************/

#include <sys/param.h>
#include "lablib3.h"                 /* prototypes for the ODL functions */
#include "ias_odl.h"
#include "ias_logging.h"

extern char ODLErrorMessage[];       /* External Variables */

IAS_OBJ_DESC *ias_odl_read_tree 
(
    const char *p_ODLFile   /* I: file to read */
)
{
    OBJDESC *p_lp= NULL;            /* pointer to struct to populate */
    char ODLPathName[MAXPATHLEN];   /* file path buffer */
    ODLErrorMessage[0] = '\0';      /* error message buffer */

    strncpy(ODLPathName, p_ODLFile, MAXPATHLEN);

    /* open and parse ODL file - error messages are suppressed */
    if ((p_lp = OdlParseLabelFile(ODLPathName, NULL, ODL_NOEXPAND, TRUE)) 
        == NULL)
    {
        IAS_LOG_ERROR("ODL Error on File: %s:  %s", ODLPathName, 
            ODLErrorMessage);
        return NULL;
    }

    /* Check of ODL syntax errors */
    if (ODLErrorMessage[0] != '\0')
    {
        IAS_LOG_ERROR("%s", ODLErrorMessage);
        OdlFreeTree(p_lp);
        return NULL;
    }

    return p_lp;
}
