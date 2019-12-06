/*************************************************************************
NAME:    ias_parm_provide_help module

PURPOSE: print out the help, template, and "load table" information

**************************************************************************/

#include <unistd.h>
#include <stdio.h>
#include <string.h>
#include "ias_parm.h"
#include "ias_logging.h"
#include "ias_const.h"

/* The maximum length of the parameter description is tied to the size of the
   description column in the database. */
#define IAS_PARM_MAX_DESC_LENGTH 256

/*************************************************************************
Name:  dump_default

Purpose:  print out the default values
*************************************************************************/
static void dump_default
(
    const IAS_PARM_PARAMETER_DEFINITION *pd
                                 /* I: pointer to current item in the list */
)
{
    int i;
    const IAS_PARM_TYPE_UNION_CONST *dv = &pd->default_values;

    if (pd->number_defaults < 1)
    {
        return;
    }

    printf("    Default = ");

    if ((pd->type == IAS_PARM_INT) || (pd->type == IAS_PARM_BOOLEAN))
    {
        for (i = 0; i < pd->number_defaults; i++)
        {
            if (i > 0)
                printf(", ");
            printf("%d", dv->type_int[i]);
        }
        printf("\n");
    }
    else if (pd->type == IAS_PARM_DOUBLE)
    {
        for (i = 0; i < pd->number_defaults; i++)
        {
            if (i > 0)
                printf(", ");
            printf("%f", dv->type_double[i]);
        }
        printf("\n");
    }
    else if (pd->type == IAS_PARM_STRING)
    {
        for (i = 0; i < pd->number_defaults; i++)
        {
            if (i > 0)
                printf(", ");
            printf("%s", dv->type_string_array[i]);
        }
        printf("\n");
    }
}


/*************************************************************************
Name:  dump_range

Purpose:  print out the range of values

Note:  In the case of an Int or Double with multiple ranges, this program
       is expecting an even number of values that will be used in pairs.
       An odd number is not understandable.  This program will ignore the
       last number in the case of an odd valid count.  This case should not
       happen and will cause an error in other routines.
*************************************************************************/
static void dump_range
(
    const IAS_PARM_PARAMETER_DEFINITION *pd,
                                /* I: pointer to current item in the list */
    int output_comments,        /* I: flag to output comment delimiters */
    int split_long_lines        /* I: flag to split lines longer than 80
                                      characters */
)
{
    int i;                      /* loop control variable */
    int line_length;            /* number of characters printed on line */
    char comment_start[10];
    char comment_end[10];
    char comment_end_no_space[10];

    if (pd->number_valids < 1)
        return;

    if (output_comments)
    {
        strcpy(comment_start, "/*");
        strcpy(comment_end, " */");
        strcpy(comment_end_no_space, "*/");
    }
    else
    {
        strcpy(comment_start, "");
        strcpy(comment_end, "");
    }

    if ((pd->type == IAS_PARM_INT) || (pd->type == IAS_PARM_BOOLEAN))
    {
        if (pd->number_valids <= 2) /* one range */
        {
            printf("%s    Range: %d to %d%s\n", comment_start,
                    pd->valid_values.type_int[0], pd->valid_values.type_int[1],
                    comment_end);
        }
        else if (pd->number_valids > 2) /* Multiple ranges listed */
        {
            line_length = printf("%s    Ranges: ", comment_start);

            /* Step through the ranges to print the values */
            for (i = 0; i < pd->number_valids; i += 2)
            {
                if (i > 0)
                    printf(", ");
                if (split_long_lines && (line_length > 80))
                {
                    printf("%s\n    %s    ", comment_end_no_space,
                            comment_start);
                    line_length = 8 + strlen(comment_start);
                }
                line_length += printf("%d to %d",pd->valid_values.type_int[i],
                            pd->valid_values.type_int[i + 1]);
            }
        printf("%s\n", comment_end);
        }
    }
    else if (pd->type == IAS_PARM_DOUBLE)
    {
        if (pd->number_valids <= 2) /* one range */
        {
            printf("%s    Range: %f to %f%s\n", comment_start,
                            pd->valid_values.type_double[0], 
                            pd->valid_values.type_double[1], comment_end);
        }
        else if (pd->number_valids > 2) /* Multiple ranges listed */
        {
            line_length = printf("%s    Ranges: ", comment_start);

            /* Step through the ranges to print the values */
            for (i = 0; i < pd->number_valids; i += 2)
            {
                if (i > 0)
                    line_length += printf(", ");
                if (split_long_lines && (line_length > 80))
                {
                    printf("%s\n    %s    ", comment_end_no_space,
                            comment_start);
                    line_length = 8 + strlen(comment_start);
                }
                line_length += printf("%f to %f",
                            pd->valid_values.type_double[i],
                            pd->valid_values.type_double[i + 1]);
            }
        printf("%s\n", comment_end);
        }
    }
    else if (pd->type == IAS_PARM_STRING)
    {
        line_length = printf("%s    Legal Values: ", comment_start);
        for (i = 0; i < pd->number_valids; i++)
        {
            if (i > 0)
                line_length += printf(", ");
            if (split_long_lines && (line_length > 80))
            {
                printf("%s\n    %s    ", comment_end_no_space, comment_start);
                line_length = 8 + strlen(comment_start);
            }
            line_length += printf("%s", pd->valid_values.type_string_array[i]);
        }
        printf("%s\n", comment_end);
    }
}


/*************************************************************************
Name:  dump_parameter_definition

Purpose:  print out the parameter definitions
*************************************************************************/
static void dump_parameter_definition
(
    const IAS_PARM_PARAMETER_DEFINITION *pd
                                 /* I: pointer to current item in the list */
)
{
    printf("Name: %s\n", pd->name);
    printf("    Description: %s\n", pd->description);
    if (pd->is_required)
    {
        printf("    Required:  Yes\n");
    }
    else
    {
        printf("    Required:  No\n");
    }

    dump_range(pd, 0, 0);
    dump_default(pd);
}


/*************************************************************************
Name:  dump_template

Purposes:  Print out the template
*************************************************************************/
static void dump_template
(
    const IAS_PARM_PARAMETER_DEFINITION *pd
                                 /* I: pointer to current item in the list */
)
{
    int i;
    const IAS_PARM_TYPE_UNION_CONST *dv = &pd->default_values;

    printf("    /* %s */\n", pd->description);

    if (pd->is_required)
    {
        printf("    /*    Required:  Yes */\n");
    }
    else
    {
        printf("    /*    Required:  No */\n");
    }

    if (pd->number_valids)
    {
        printf("    ");
        dump_range(pd, 1, 1);
    }

    printf("    %s = ", pd->name);

    if (pd->number_defaults > 1)
        printf("(");

    if ((pd->type == IAS_PARM_INT) || (pd->type == IAS_PARM_BOOLEAN))
    {
        for (i = 0; i < pd->number_defaults; i++)
        {
            if (i > 0)
                printf(", ");
            printf("%d", dv->type_int[i]);
        }
    }
    else if (pd->type == IAS_PARM_DOUBLE)
    {
        for (i = 0; i < pd->number_defaults; i++)
        {
            if (i > 0)
                printf(", ");
            printf("%f", dv->type_double[i]);
        }
    }
    else if (pd->type == IAS_PARM_STRING)
    {
        for (i = 0; i < pd->number_defaults; i++)
        {
            if (i > 0)
                printf(", ");
            printf("\"%s\"", dv->type_string_array[i]);
        }
    }

    if (pd->number_defaults > 1)
        printf(")");

    printf("\n");
}

/*************************************************************************
Name:  dump_loadtable_parameters

Purpose:  Print out the loadtable PARAMETERS table

Returns:
    SUCCESS or ERROR
*************************************************************************/
static int dump_loadtable_parameters
(
    const IAS_PARM_PARAMETER_DEFINITION *pd,
                                /* I: pointer to current item in the list */
    int param_number            /* I: If this is zero, print headers */
)
{
    int i;
    int size_of_array = 0;      /* number of elements defined in the array */

    /* Load file for the PARAMETERS table */

    /* Print the heading the first time only */
    if (param_number == 0)
    {
        printf("\nload data\ninfile *\nappend\ninto table \"PARAMETERS\"\n");
        printf("fields terminated by '|'\ntrailing nullcols (\n");
        printf("MODULE_ID\n,PARM_NAME\n,PARM_OCCURRENCE\n,PARM_TYPE\n )\n");
        printf("begindata\n");
    }

    if ((pd->type == IAS_PARM_INT) || (pd->type == IAS_PARM_BOOLEAN))
    {
        if (!pd->is_an_array)
        {
            printf("<SCRIPT>|%s|1|long\n", pd->name);
        }
        else  /* Step through array and write each parameter line */
        {
            /* Get the number of elements in the array */
            size_of_array = (pd->value_bytes / sizeof(int));
            for (i = 0; i < size_of_array; i++)
            {
                printf("<SCRIPT>|%s|%d|long\n", pd->name, i + 1);
            }
        }
    }
    else if (pd->type == IAS_PARM_DOUBLE)
    {
        if (!pd->is_an_array)
        {
            printf("<SCRIPT>|%s|1|double\n", pd->name);
        }
        else  /* Step through array and write each parameter line */
        {
            /* Get the number of elements in the array */
            size_of_array = (pd->value_bytes / sizeof(double));
            for (i = 0; i < size_of_array; i++)
            {
                printf("<SCRIPT>|%s|%d|double\n", pd->name, i + 1);
            }
        }
    }
    else if (pd->type == IAS_PARM_STRING)
    {
        if (!pd->is_an_array)
        {
            printf("<SCRIPT>|%s|1|char\n", pd->name);
        }
        else  /* Step through array and write each parameter line */
        {
            /* Get the number of elements in the array */
            size_of_array = (pd->value_bytes / sizeof(&pd->valid_values));
            for (i = 0; i < size_of_array; i++)
            {
                printf("<SCRIPT>|%s|%d|char\n", pd->name, i + 1);
            }
        }
    }
    else    /* Do not recognize the type */
    {
        IAS_LOG_ERROR("Do not recognize type %d for %s",pd->type, pd->name);
        return ERROR;
    }

    return SUCCESS;
}

/*************************************************************************
Name:       dump_loadtable_parm_edits

Purpose:    print out the loadtable PARM_EDITS table

Returns:
    SUCCESS or ERROR

Note:       When the the number of constrained values is more than 2 for 
            integer and double values, that indicates that this is not a 
            true range at all, but a list of valid numbers.  
            Therefore "range" should really be a min and max of the 
            same number.  If not, we really don't support this constraint.  
            Multiple ranges of numerical values are expected to be something 
            like element 1: min 15, 2: max 15, 3: min 30,
            4: max 30, 5: min 60 6: max 60 ... meaning 15, 30, 60... are
            the only valid integers for that parameter
            If not set up this way, an error will be returned.

*************************************************************************/
static int dump_loadtable_parm_edits
(
    const IAS_PARM_PARAMETER_DEFINITION *pd,
                                /* I: pointer to current item in the list */
    int param_number            /* I: If this is zero, print headers */
)
{
    int i;

    /* Load file for the PARM_EDITS table */

    /* Print the heading the first time only */
    if (param_number == 0)
    {
        printf("\nload data\ninfile *\nappend\ninto table \"PARM_EDITS\"\n");
        printf("fields terminated by '|'\ntrailing nullcols (\n");
        printf("PARM_NAME\n,PARM_DESC\n,MIN_VALUE\n,MAX_VALUE )\n");
        printf("begindata\n");
    }

    /* Print the name and description of the parameter */
    printf("%s|%s|", pd->name, pd->description);

    /* Print the min/max if there is one for the numbers only */
    if ((pd->type == IAS_PARM_INT) || (pd->type == IAS_PARM_BOOLEAN))
    {
        if (pd->number_valids)
        {
            if (pd->number_valids > 2)
            {
                for (i = 0; i < pd->number_valids; i += 2)
                {
                    /* Write an error if the ranges are not the same number
                       See the note above */
                    if (pd->valid_values.type_int[i] !=
                                        pd->valid_values.type_int[i + 1])
                    {
                        IAS_LOG_ERROR("Parameter %s has multiple ranges "
                                "which are not supported",pd->name);
                        return ERROR;
                    }
                } 
                printf("|\n");
            }
            else
                printf("%d|%d\n", pd->valid_values.type_int[0],
                                  pd->valid_values.type_int[1]);
        }
        else
            printf("|\n");
    }
    else if (pd->type == IAS_PARM_DOUBLE)
    {
        if (pd->number_valids)
        {
            if (pd->number_valids > 2)
            {
                for (i = 0; i < pd->number_valids; i += 2)
                {
                    /* Write an error if the ranges are not the same number
                       See the note above */
                    if (pd->valid_values.type_double[i] !=
                                        pd->valid_values.type_double[i + 1])
                    {
                        IAS_LOG_ERROR("Parameter %s has multiple ranges "
                                "which are not supported",pd->name);
                        return ERROR;
                    }
                }
                printf("|\n");
            }
            else
                printf("%f|%f\n", pd->valid_values.type_double[0],
                                  pd->valid_values.type_double[1]);
        }
        else
            printf("|\n");
    }
    /* No min/max is possible for a string.  Add closing pipe */
    else if (pd->type == IAS_PARM_STRING)
    {
        printf("|\n");
    }
    else    /* Do not recognize the type */
    {
        IAS_LOG_ERROR("Do not recognize type %d for %s",pd->type, pd->name);
        return ERROR;
    }

return SUCCESS;
}

/*************************************************************************
Name:  dump_loadtable_parm_list

Returns:
    SUCCESS or ERROR

Purpose:   Print out the loadtable PARM_LIST table

*************************************************************************/
static int dump_loadtable_parm_list
(
    const IAS_PARM_PARAMETER_DEFINITION *pd,
                                /* I: pointer to current item in the list */
    int param_number            /* I: If this is zero, print headers */
)
{
    int i;

    /* Load file for the PARM_LIST table */

    /* Print the heading the first time only */
    if (param_number == 0)
    {
        printf("\nload data\ninfile *\nappend\ninto table \"PARM_LIST\"\n");
        printf("fields terminated by '|'\ntrailing nullcols (\n");
        printf("PARM_NAME\n,OUI_VALUE\n,PROGRAM_VALUE )\n");
        printf("begindata\n");
    }

    if ((pd->type == IAS_PARM_INT) || (pd->type == IAS_PARM_BOOLEAN))
    {
        /* There is a list of numerical values, not a range */
        if (pd->number_valids > 2)
        {
            for (i = 0; i < pd->number_valids; i += 2)
            {
                printf("%s|", pd->name);
                printf("%d|\n",pd->valid_values.type_int[i]);
            }
        }
    }
    else if (pd->type == IAS_PARM_DOUBLE)
    {
        /* There is a list of numerical values, not a range */
        if (pd->number_valids > 2)
        {
            for (i = 0; i < pd->number_valids; i += 2)
            {
                printf("%s|", pd->name);
                printf("%f|\n",pd->valid_values.type_double[i]);
            }
        }
    }
    else if (pd->type == IAS_PARM_STRING)
    {
        if (pd->number_valids)
        {
            for (i = 0; i < pd->number_valids; i++)
            {
                printf("%s|", pd->name);
                printf("%s|\n",pd->valid_values.type_string_array[i]);
            }
        }
    }
    else    /* Do not recognize the type */
    {
        IAS_LOG_ERROR("Do not recognize type %d for %s",pd->type, pd->name);
        return ERROR;
    }

return SUCCESS;
}

/*************************************************************************
Name:       dump_loadtable_def_parms

Returns:
    SUCCESS or ERROR

Purpose:    print out the loadtable DEF_PARMS table

*************************************************************************/
static int dump_loadtable_def_parms
(
    const IAS_PARM_PARAMETER_DEFINITION *pd,
                                /* I: pointer to current item in the list */
    int param_number            /* I: If this is zero, print headers */
)
{
    int count;              /* loop control variable */
    char typestring[7];     /* string name for type */
    int param_count;        /* number of rows for this parameter */

    /* Load file for the DEF_PARMS table */

    /* Print the heading the first time only */
    if (param_number == 0)
    {
        printf("\nload data\ninfile *\nappend\ninto table \"DEF_PARMS\"\n");
        printf("fields terminated by '|'\ntrailing nullcols (\n");
        printf("MODULE_ID\n,SUB_MODULE_ID\n,PARM_NAME\n,PARM_OCCURRENCE\n");
        printf(",GLOBAL\n,PARM_TYPE\n,REQUIRED\n,PARM_VALUE\n )");
        printf("begindata\n");
    }

    if ((pd->type == IAS_PARM_INT) || (pd->type == IAS_PARM_BOOLEAN))
    {
        param_count = (pd->value_bytes / sizeof(int));
        strcpy(typestring, "long");
    }
    else if (pd->type == IAS_PARM_DOUBLE)
    {
        param_count = (pd->value_bytes / sizeof(double));
        strcpy(typestring, "double");
    }
    else if (pd->type == IAS_PARM_STRING)
    {
        param_count = (pd->value_bytes / sizeof(char *));
        strcpy(typestring, "char");
    }
    else    /* Do not recognize the type */
    {
        IAS_LOG_ERROR("Do not recognize type %d for %s",pd->type, pd->name);
        return ERROR;
    }

    /* The parm count is calculated above.  If this is not an array, but is
       a string, this method will not be valid so just set the count to one.
       if it's not a string, return error */

    if (!pd->is_an_array)
    {
        /* For a non-array string, just set the count to 1 */
        if (pd->type == IAS_PARM_STRING)
        {
            param_count = 1;
        }
        /* For the numeric values, the count should be set correctly and
           if not, there was an error in the input value and byte sizes */
        else if (param_count != 1)
        {
            IAS_LOG_ERROR("The return value and number of bytes are "
                "invalid for non array %s count: %d", pd->name, param_count);
            return ERROR;
        }
    }

    /* Loop through each parameter in the array */
    for (count = 0; count < param_count; count++)
    {
        /* Begin writing the line */
        printf("<PROCEDURE>|<SCRIPT>|%s|%d|<GLOBAL>|%s|%d|", pd->name, count+1,
                                                typestring,pd->is_required);

        /* Print the default if there is one */
        if ((pd->number_defaults > 0) && (pd->number_defaults - 1 >= count))
        {
            if ((pd->type == IAS_PARM_INT) || (pd->type == IAS_PARM_BOOLEAN))
                printf("%d\n",pd->default_values.type_int[count]);
            else if (pd->type == IAS_PARM_DOUBLE)
                printf("%f\n",pd->default_values.type_double[count]);
            else if (pd->type == IAS_PARM_STRING)
                printf("%s\n",pd->default_values.type_string_array[count]);
        }
        else
            printf("\n");
    }

    return SUCCESS;
}


/*************************************************************************
Name:  ias_parm_provide_help

Purpose:  print out the template or help info

Returns:
    ERROR if an error happens
    0 if no help was provided
    1 if help was provided
*************************************************************************/
int ias_parm_provide_help
(
    const char *option,          /* I: parameter to show help or template */
    IAS_PARM_PARAMETER_DEFINITION **params,
                                 /* I/O : pointer to list of items to read
                                    from the ODL file */
    int count,                   /* I: number of items */
    IAS_PARAMETER_SOURCE file_source_type /* I: Type of file source, OMF and 
                                    Parameter File */
)
{
    int i;

    /* Verify the parameter descriptions aren't larger than 
       IAS_PARM_MAX_DESC_LENGTH characters */
    for(i = 0; i < count; i++)
    {
        if (strlen(params[i]->description) > IAS_PARM_MAX_DESC_LENGTH)
        {
            IAS_LOG_ERROR("Parameter description larger than %d characters "
                "for %s", IAS_PARM_MAX_DESC_LENGTH, params[i]->name);
            return ERROR;
        }
    }

    /* Output the template data.  this data is in the form of an ODL file */
    if(strcmp(option, "--template") == 0)
    {
        if (file_source_type == IAS_INPUT_PARAMETERS)
            printf("OBJECT = PARAMETERS\n");
        else if (file_source_type == IAS_OMF_PARAMETERS)
            printf("OBJECT = OMF\n");
        else
        {
            IAS_LOG_ERROR("Invalid parameter source");
            return ERROR;
        }

        for(i = 0; i < count; i++)
            dump_template(params[i]);

        if (file_source_type == IAS_INPUT_PARAMETERS)
            printf("END_OBJECT = PARAMETERS\n");
        else if (file_source_type == IAS_OMF_PARAMETERS)
            printf("END_OBJECT = OMF\n");
        printf("END\n");
        return 1;
    }
    /* output from loadtable represents the data to be entered in the database
       to support a script */
    else if(strcmp(option, "--loadtable") == 0)
    {
        /* No action if this is OMF */
        if (file_source_type == IAS_OMF_PARAMETERS)
        {
            return 1;
        }
        /* in the loadtable case, there is no need to check the output type */
        for(i = 0; i < count; i++)
        {
            if (dump_loadtable_parameters(params[i],i) != SUCCESS)
            {
                IAS_LOG_ERROR("An error occurred writing load table output "
                                "for the PARAMETERS table");
                return ERROR;
            }
        }

        printf("\n");
        for(i = 0; i < count; i++)
        {
            if (dump_loadtable_parm_edits(params[i],i) != SUCCESS)
            {
                IAS_LOG_ERROR("An error occurred writing load table output "
                                "for the PARM_EDITS table");
                return ERROR;
            }
        }

        printf("\n");
        for(i = 0; i < count; i++)
        {
            if (dump_loadtable_parm_list(params[i],i) != SUCCESS)
            {
                IAS_LOG_ERROR("An error occurred writing load table output "
                                "for the PARM_LIST table");
                return ERROR;
            }
        }

        printf("\n");
        for(i = 0; i < count; i++)
        {
            if (dump_loadtable_def_parms(params[i],i) != SUCCESS)
            {
                IAS_LOG_ERROR("An error occurred writing load table output "
                                "for the DEF_PARMS table");
                return ERROR;
            }
        }

        printf("\n");
        return 1;
    }
    /* The output from help is in a general format desribing the parameters */
    else if(strcmp(option, "--help") == 0)
    {
        if (file_source_type == IAS_INPUT_PARAMETERS)
            printf("\n### Application Input Parameters ###\n\n");
        else if (file_source_type == IAS_OMF_PARAMETERS)
            printf("\n### OMF Parameters ###\n\n");
        else
        {
            IAS_LOG_ERROR("Invalid parameter source");
            return ERROR;
        }

        for(i = 0; i < count; i++)
            dump_parameter_definition(params[i]);

        printf("\n");
        return 1;
    }

    return 0;
}
