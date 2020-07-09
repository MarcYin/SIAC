/******************************************************************************

UNIT NAME: ias_odl_get_field.c

PURPOSE: Get the requested ODL field and convert it, if necessary.

RETURN VALUE:
    Type = int

Value                                   Description
--------------                          ----------------------------------------
SUCCESS                                 Found and converted field(s) into 
                                            requested type
IAS_ODL_NOT_ENOUGH_MEMORY_SUPPLIED      Not enough memory passed in
IAS_ODL_NOT_FOUND                       Group/label not found
IAS_ODL_INVALID_DATA_TYPE               Data type mismatch
ERROR                                   Fatal error                       

******************************************************************************/

#include "lablib3.h"
#include "ias_odl.h"
#include "ias_logging.h"

extern char ODLErrorMessage[];       /* External Variables */

/* function prototype */
int convert_string
( 
    void *p_destination,            /* I/O: attribute to convert */
    IAS_ODL_TYPE parm_type,         /* I: ODL data type */
    char *kvalue                    /* I: Value to convert */
);

int ias_odl_get_field
( 
    void *p_MemoryAddr,             /* I/O: List of attributes to retrieve */
    int MemorySize,                 /* I: mem size of attributes */ 
    IAS_ODL_TYPE ValueType,         /* I: ODL data type */
    IAS_OBJ_DESC *p_ODLTree,        /* I: ODL tree */
    const char *p_ClassName,        /* I: Group/Object name */
    const char *p_LabelName,        /* I: Field to get */
    int *p_Count                    /* I: number of values in attribute */
)
{
    OBJDESC *p_lp;              /* Object Descriptor */
    KEYWORD *p_kw;              /* Keyword Name */
    char *p_kwv;                /* Keyword Value */
    char *p_keyword;            /* Copy of the Keyword Value */
    int i;                      /* loop counter */
    char *p_word;               /* word to convert */
    int ret_code = 0;           /* function return value */

    *p_Count = 0;

    if ( (p_LabelName == NULL) || (strlen(p_LabelName) == 0) )
    {
        IAS_LOG_ERROR("Attribute name missing");
        return ERROR;
    }

    if ((p_lp = OdlFindObjDesc(p_ODLTree, p_ClassName, p_LabelName, NULL, 1, 
        ODL_RECURSIVE_DOWN)) == NULL)
    {
        /* it's possible that we have NULL for a group(p_ClassName), that's OK;
           no need to log an error, just return */
        return IAS_ODL_NOT_FOUND;
    }

    if ((p_kw = OdlFindKwd(p_lp, p_LabelName, NULL, 1 ,ODL_RECURSIVE_DOWN)) 
        == NULL)
    {
        if ((long)strlen(ODLErrorMessage) <= 1 )
        {
            IAS_LOG_ERROR("%s", ODLErrorMessage);
            IAS_LOG_ERROR("Keyword '%s' not found", p_LabelName);
        }
        return IAS_ODL_NOT_FOUND;
    }

    if ((p_kwv = OdlGetKwdValue(p_kw)) == NULL)
    {
        if ((long)strlen(ODLErrorMessage) <= 1 )
        {
            IAS_LOG_ERROR("%s", ODLErrorMessage);
            IAS_LOG_ERROR("Keyword %s not found", p_LabelName);
        }
        return IAS_ODL_NOT_FOUND;
    }
    if ((p_keyword= malloc(strlen(p_kwv)+1)) == NULL)
    {
        IAS_LOG_ERROR("Malloc error");
        return ERROR;
    }
    (void)strcpy(p_keyword,p_kwv);

    /* When working with a set or a sequence, all of the newline characters
       should be removed from the incoming string.  The ODL library is placing
       newlines into these strings around character 2025 even though they don't
       exist in the ODL file. */
    if ((OdlGetKwdValueType(p_kw) == ODL_SET) ||
        (OdlGetKwdValueType(p_kw) == ODL_SEQUENCE))
    {
        char * p_newline;

        while ((p_newline = strchr(p_keyword, '\n')) != NULL)
        {
            while (*p_newline != '\0')
            {
                *p_newline = *(p_newline + 1);
                p_newline++;
            }
        }
    }

    if (OdlGetKwdValueType(p_kw) == ODL_SET)
    {
        p_word = strtok(p_keyword,"(,) \"\n");
        while(p_word != NULL)
        {
            if (strlen(p_word))
            {
                switch (ValueType)
                {
                case IAS_ODL_Long :
                    MemorySize -= sizeof(long);
                    if (MemorySize < 0 )
                    {
                        IAS_LOG_ERROR("Input ODL value overflows allocated "
                            "space ");
                        free(p_keyword);
                        return IAS_ODL_NOT_ENOUGH_MEMORY_SUPPLIED;
                    }
                    if ( ( ret_code = convert_string(p_MemoryAddr,
                        ValueType,p_word) ) != SUCCESS )
                    {
                        IAS_LOG_ERROR("Converting %s keyword's value "
                            "%s",p_LabelName, p_word);
                        free(p_keyword);
                        return ret_code;
                    }
                    p_MemoryAddr = (long *)p_MemoryAddr + 1;
                    break;

                case IAS_ODL_Int :
                    MemorySize -= sizeof(int);
                    if (MemorySize < 0 )
                    {
                        IAS_LOG_ERROR("Input ODL value overflows allocated "
                            "space ");
                        free(p_keyword);
                        return IAS_ODL_NOT_ENOUGH_MEMORY_SUPPLIED;
                    }
                    if ( ( ret_code = convert_string(p_MemoryAddr,
                        ValueType,p_word) ) != SUCCESS )
                    {
                        IAS_LOG_ERROR("Converting %s keyword's value "
                            "%s",p_LabelName, p_word);
                        free(p_keyword);
                        return ret_code;
                    }
                    p_MemoryAddr = (int *)p_MemoryAddr + 1;
                    break;

                case IAS_ODL_Float :
                    MemorySize -= sizeof(float);
                    if (MemorySize < 0 )
                    {
                        IAS_LOG_ERROR("Input ODL value overflows allocated "
                            "space ");
                        free(p_keyword);
                        return IAS_ODL_NOT_ENOUGH_MEMORY_SUPPLIED;
                    }
                    if ( ( ret_code = convert_string(p_MemoryAddr,
                        ValueType,p_word) ) != SUCCESS )
                    {
                        IAS_LOG_ERROR("Converting %s keyword's value "
                            "%s",p_LabelName, p_word);
                        free(p_keyword);
                        return ret_code;
                    }
                    p_MemoryAddr = (float *)p_MemoryAddr + 1;
                    break;

                case IAS_ODL_Double :
                case IAS_ODL_Sci_Not :
                    MemorySize -= sizeof(double);
                    if (MemorySize < 0 )
                    {
                        IAS_LOG_ERROR("Input ODL value overflows allocated "
                            "space ");
                        free(p_keyword);
                        return IAS_ODL_NOT_ENOUGH_MEMORY_SUPPLIED;
                    }
                    if ( ( ret_code = convert_string(p_MemoryAddr,
                        ValueType,p_word) ) != SUCCESS )
                    {
                        IAS_LOG_ERROR("Converting %s keyword's value "
                            "%s",p_LabelName, p_word);
                        free(p_keyword);
                        return ret_code;
                    }
                    p_MemoryAddr = (double *)p_MemoryAddr + 1;
                    break;

                case IAS_ODL_ArrayOfString :
                    MemorySize -= sizeof(char *);
                    if (MemorySize < 0 )
                    {
                        IAS_LOG_ERROR("Input ODL value overflows allocated "
                            "space ");
                        free(p_keyword);
                        return  IAS_ODL_NOT_ENOUGH_MEMORY_SUPPLIED;
                    }
                    if ( ( ret_code = convert_string(p_MemoryAddr,
                        ValueType,p_word) ) != SUCCESS )
                    {
                        for( i=0;i<*p_Count;i++)
                            free((char *)p_MemoryAddr + i);
                        IAS_LOG_ERROR("Converting %s keyword's value "
                            "%s",p_LabelName, p_word);
                        free(p_keyword);
                        return ret_code;
                    }
                    p_MemoryAddr = (char *)p_MemoryAddr + sizeof(char *);
                    break;

                default:
                    (void)strncpy( p_MemoryAddr, p_word, MemorySize);
                    IAS_LOG_ERROR("Type IAS_ODL_String is not valid for arrays "
                        "of strings");
                    free(p_keyword);
                    return ERROR;
                }
                *p_Count += 1;
                p_word = strtok(NULL,",() \"\n");
            }
        }
    }
    else if (OdlGetKwdValueType(p_kw) == ODL_SEQUENCE)
    {
        p_word = strtok(p_keyword,"{,} \"\n");
        while(p_word != NULL)
        {
            if (strlen(p_word))
            {
                switch (ValueType)
                {
                case IAS_ODL_ArrayOfString :
                    MemorySize -= sizeof(char *);
                    if (MemorySize < 0 )
                    {
                        IAS_LOG_ERROR("Input ODL value overflows allocated "
                            "space ");
                        free(p_keyword);
                        return  IAS_ODL_NOT_ENOUGH_MEMORY_SUPPLIED;
                    }
                    if ( ( ret_code = convert_string(p_MemoryAddr,
                        ValueType,p_word) ) != SUCCESS )
                    {
                        for( i=0;i<*p_Count;i++)
                            free((char *)p_MemoryAddr + i);
                        IAS_LOG_ERROR("Converting %s keyword's value "
                            "%s", p_LabelName, p_word);
                        free(p_keyword);
                        return ret_code;
                    }
                    p_MemoryAddr = (char *)p_MemoryAddr + sizeof(char *);
                    break;

                default:
                    (void)strncpy( p_MemoryAddr, p_word, MemorySize);
                    IAS_LOG_ERROR("An ODL sequence is of type "
                        "IAS_ODL_ArrayOfString only");
                    free(p_keyword);
                    return ERROR;
                }
                *p_Count += 1;
                p_word = strtok(NULL,",{} \"\n");
            }
        }
    }
    else
    {
        if (ValueType == IAS_ODL_String)
        {
            char * p_start = p_kwv;
            int    len     = strlen(p_kwv);

            /* If we are working with a double quoted string, we'll need to
               remove the double quotes in the string that is passed back.
               Adjust the length of the string and the starting point for
               copying, if necessary. */
            if (p_kwv[0] == '\"')
            {
                len -= 2;
                p_start++;
            }

            if (MemorySize < (len + 1))
            {
                IAS_LOG_ERROR("Input ODL value overflows allocated space");
                (void)strncpy(p_MemoryAddr, p_start, MemorySize);
                ((char *)p_MemoryAddr)[MemorySize-1] = '\0';
                free(p_keyword);
                return IAS_ODL_NOT_ENOUGH_MEMORY_SUPPLIED;
            }

            (void)strncpy(p_MemoryAddr, p_start, len);
            ((char *)p_MemoryAddr)[len] = '\0';
        }
        else
        {
            switch(ValueType)
            {
            case IAS_ODL_Long:
                MemorySize -= sizeof(long);
                break;

            case IAS_ODL_Int:
                MemorySize -= sizeof(int);
                break;

            case IAS_ODL_Float:
                MemorySize -= sizeof(float);
                break;

            case IAS_ODL_Double:
            case IAS_ODL_Sci_Not:
                MemorySize -= sizeof(double);
                break;

            default:
                IAS_LOG_ERROR("Invalid Type specified");
                free(p_keyword);
                return ERROR;
            }

            if (MemorySize < 0 )
            {
                IAS_LOG_ERROR("Input ODL value overflows allocated space ");
                free(p_keyword);
                return ERROR;
            }

            if ( ( ret_code = convert_string(p_MemoryAddr,ValueType,
                p_kwv) ) != SUCCESS )
            {
                IAS_LOG_ERROR("Converting %s keyword's value %s",
                    p_LabelName,p_kwv);
                free(p_keyword);
                return ret_code;
            }
        }
        *p_Count += 1;
    }
    free(p_keyword);
    return SUCCESS;
}


/******************************************************************************

UNIT NAME: convert_string.c

PURPOSE: To convert the string to the new type using the memory address given by
         caller

INVOCATION METHOD:
    status = convert_string(p_destination, parm_type, kvalue)

RETURN VALUE:
    Type = int

Value                           Description
---------------                 ------------------------------------------------
SUCCESS                         Converted the String in to requested type
ERROR                           Allocation error or unknown parm_type
IAS_ODL_INVALID_DATA_TYPE       Invalid data type

******************************************************************************/
int convert_string
( 
    void *p_destination,            /* I/O: converted value location */ 
    IAS_ODL_TYPE parm_type,         /* I; ODL data type */
    char *kvalue)                   /* I: word to convert */
{
    char *p_endptr;
    char *p_a;

    switch(parm_type)
    {
        case IAS_ODL_Long:
            (*(long *)p_destination) = strtol(kvalue,&p_endptr,10);
            if (*p_endptr != '\0')
            {
                return IAS_ODL_INVALID_DATA_TYPE;
            }
            break;

        case IAS_ODL_Int:
            (*(int *)p_destination) = (int)strtol(kvalue,&p_endptr,10);
            if (*p_endptr != '\0')
            {
                return IAS_ODL_INVALID_DATA_TYPE;
            }
            break;

        case IAS_ODL_Float:
            (*(float *)p_destination) = (float)strtod(kvalue,&p_endptr);
            if (*p_endptr != '\0')
            {
                return IAS_ODL_INVALID_DATA_TYPE;
            }
            break;

        case IAS_ODL_Sci_Not:
        case IAS_ODL_Double:
            (*(double *)p_destination) = strtod(kvalue,&p_endptr);
            if (*p_endptr != '\0')
            {
                return IAS_ODL_INVALID_DATA_TYPE;
            }
            break;

        case IAS_ODL_ArrayOfString:
            p_a = malloc(strlen(kvalue)+1);
            if (!p_a)
            {
                IAS_LOG_ERROR("Allocating Memory ... ");
                return ERROR;
            }
            (void)strcpy(p_a,kvalue);
            (*(char **)p_destination) = p_a;
            break;

        default:
            IAS_LOG_ERROR("Invalid Conversion type %s ...", kvalue);
            return ERROR;
    }

    return SUCCESS;
}

