/******************************************************************************
NAME:           ias_misc_convert_to_uppercase

PURPOSE:        
ias_misc_convert_to_uppercase converts all the lower case characters in a 
string to upper case.

RETURN VALUE:
Type = char *
Value    Description
-----    -----------
pointer  Returns a pointer to the string so it can immediately be used as
         a parameter to another function (such as strcmp)

NOTES:
The original string is converted in place.
******************************************************************************/
#include <ctype.h>              /* toupper prototype */
#include "ias_miscellaneous.h"

char *ias_misc_convert_to_uppercase 
(
    char *string_ptr  /* I/O: pointer to string to convert */
)
{
    char *c_ptr = string_ptr;

    while (*c_ptr != '\0')
    {
        *c_ptr = toupper(*c_ptr);
        c_ptr++;
    }

    return string_ptr;
}
