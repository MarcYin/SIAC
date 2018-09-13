/*************************************************************************

NAME: ias_logging.c

PURPOSE: Implements the standard message logging interface 

**************************************************************************/

#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <stdio.h>

#define LOGGING_C
#include "ias_logging.h"
#include "ias_const.h"

#define CHANNELS_LENGTH 500 

/*************************************************************************/
enum IAS_LOG_MESSAGE_LEVEL ias_log_message_level; /* Log meeasge level value */
static FILE *file_ptr = NULL;                     /* Output file pointer */
static char program_name[40];                     /* Current running program 
                                                     name */
static pid_t pid;                                 /* Current running processor
                                                     id */
static char ias_log_channels[CHANNELS_LENGTH+3];  /* List of channels to 
                                                     log (3 additional chars:
                                                     start/stop commas and
                                                     null character */
static int channels_on = 0;                       /* Flag set if channels
                                                     are enabled */
static int use_blacklist = 0;                     /* Flag indicates if the
                                                     channel list contains
                                                     disabled channels */
static int dump_registered = 0;                   /* Flag that indicates
                                                     atexit dump handler was
                                                     registered */
/*************************************************************************/

/*************************************************************************

NAME: channel_enabled

PURPOSE: Checks if the specified channel is enabled. 

RETURNS: FALSE -- channel disabled 
         TRUE  -- channel enabled

**************************************************************************/
static int is_channel_enabled
(
    const char *channel
)
{
    /* if no channels are defined all are enabled */
    if ( strlen(ias_log_channels) == 0 )
        return TRUE;

    /* build the search string (prepend/append commas) */
    char search_str[strlen(channel) + 3];
    sprintf(search_str, ",%s,", channel);

    /* return true if the channel is found in the channel list */
    if ( strstr(ias_log_channels, search_str) != NULL )
    {
        if ( use_blacklist )
            return FALSE;
        else
            return TRUE;
    }

    if ( use_blacklist )
        return TRUE;
    else
        return FALSE;
}

/*************************************************************************

NAME: dump_io

PURPOSE: Dumps the process IO file to the log channel 'IO'

RETURNS: Void

**************************************************************************/
static void dump_io
(
)
{
    /* If dump io is registered then we know that the IO channel is enabled
       and the log level <= IAS_LOG_LEVEL_DEBUG */
    char buffer[1000];
    FILE *fp;

    /* Prepare the filename */
    sprintf(buffer, "/proc/%d/io", pid);
    fp = fopen(buffer, "r");
    if (!fp)
    {
        IAS_LOG_WARNING("Unable to open the %s file - not reporting IO stats",
            buffer);
        return;
    }

    /* Reset the output to stdout incase the user is targeting
       a separate log file since the log file has likely been closed
       before this runs */
    file_ptr = stdout;

    /* Loop over each line */
    while (fgets(buffer, sizeof(buffer), fp))
    {
        int last_char = strlen(buffer) - 1;

        /* Remove the new line if present */
        if (last_char >= 0 && buffer[last_char] == '\n')
        {
            buffer[last_char] = '\0';
        }

        /* Log the IO information */
        ias_log_message_with_channel(IAS_LOG_LEVEL_DEBUG, "IO", __FILE__,
            __LINE__, buffer);
    }

    pclose(fp);
}

/*************************************************************************

NAME: ias_log_initialize

PURPOSE: Initilizes the logging library call

Algorithm References: None

RETURNS: SUCCESS -- successfully initialized
         ERROR -- error in initialization 

NOTES:

**************************************************************************/
int ias_log_initialize
(
    const char *log_program_name /* I: name to output with each log message */
)
{
    const char *log_level;
    const char *log_channels;
 
    /* set the program name */
    strncpy(program_name, log_program_name, sizeof(program_name));
    program_name[sizeof(program_name) - 1] = '\0';

    if (file_ptr == NULL)
        file_ptr = stdout;        /* set the output pointer to stdout */
    pid = getpid();               /* get current process id */

    log_level = getenv("IAS_LOG_LEVEL");
    if ( log_level != NULL)
    {
        if (strcmp(log_level,"DEBUG") == 0)
            ias_log_message_level = IAS_LOG_LEVEL_DEBUG;
        else if (strcmp(log_level,"INFO") == 0)
            ias_log_message_level = IAS_LOG_LEVEL_INFO;
        else if (strcmp(log_level,"WARN") == 0)
            ias_log_message_level = IAS_LOG_LEVEL_WARN;
        else if (strcmp(log_level,"ERROR") == 0)
            ias_log_message_level = IAS_LOG_LEVEL_ERROR;
        else 
        {
            IAS_LOG_ERROR("Environment Variable IAS_LOG_LEVEL needs "
                         "to be 'DEBUG', 'INFO', 'WARN', or 'ERROR'");
            return ERROR;
        }
    }
    else
    {
        ias_log_message_level = IAS_LOG_LEVEL_INFO;
    }

    log_channels = getenv("IAS_LOG_CHANNELS");
    if ( log_channels != NULL )
    {
        if ( strlen(log_channels) > CHANNELS_LENGTH )
        {
            IAS_LOG_ERROR("Environment variable IAS_LOG_CHANNELS exceeds max "
                "length of %d", CHANNELS_LENGTH);
            return ERROR;
        }

        /* if channel list starts with '-' treat the list as a black list 
           rather than a white list */
        if (strlen(log_channels) > 0 && log_channels[0] == '-')
            use_blacklist = 1;
    
        /* The value of IAS_LOG_CHANNELS should be a comma delimited list. By
           appending a comma to the start and end of the value, code can check
           if a channel is enabled by searching for ',channel,'. */
        if ( use_blacklist )
        {
            /* don't copy the '-' character */
            sprintf(ias_log_channels, ",%s,", (log_channels+1));
        }
        else
            sprintf(ias_log_channels, ",%s,", log_channels);
        channels_on = 1;
    }
    else
    {
        ias_log_channels[0] = '\0';
    }

    /* line based buffering to the output file to prevent delay */
    if (setvbuf(file_ptr, NULL, _IOLBF, 0) !=0)  
    {
        IAS_LOG_WARNING("Incorrect type or size of buffer for file_ptr");   
    }

    /* Register the at exit handler if not already registered */
    if (is_channel_enabled("IO") && (IAS_LOG_LEVEL_DEBUG
        <= ias_log_message_level) && !dump_registered)
    {
        if (atexit(dump_io) == 0)
            dump_registered = 1;
    }

    return SUCCESS;
} 

/*************************************************************************

NAME: ias_log_set_output_level

PURPOSE: Sets the logging level to the new level

RETURNS: None

**************************************************************************/
IAS_LOG_MESSAGE_LEVEL ias_log_set_output_level
(
    int new_level   /* I: minimum logging level to output */
)
{
    IAS_LOG_MESSAGE_LEVEL old_level;
 
    old_level = ias_log_message_level;
    ias_log_message_level = new_level;

    return old_level;
}

/*************************************************************************

NAME: ias_log_set_output_target

PURPOSE: Sets the output file pointer to the new value

RETURNS: SUCCESS -- successfully setting target
         ERROR -- error in setting target

**************************************************************************/
int ias_log_set_output_target
(
    FILE *new_fp    /* I: File pointer for output message */ 
)
{
    file_ptr = new_fp; 

    /* Line based buffering to the output file to prevent delay */ 
    if (setvbuf(new_fp, NULL, _IOLBF, 0) !=0) 
    {
        IAS_LOG_WARNING("Incorrect type or size of buffer for new_fp");
    }

    return SUCCESS;
}

/*************************************************************************

NAME: format_time

PURPOSE: Gets current timestamp

RETURNS: SUCCESS -- successfully getting time
         ERROR -- error in getting time

**************************************************************************/
static int format_time
(
    char *stamp,           /* O: timestamp for output */
    int stampsize,         /* I: size of timestamp for input */
    const char *format     /* I: format of output timestamp */
)
{
    time_t ptime;                 /* Time in seconds  */
    struct tm *ltime;             /* Time in local time */

    /* Get the current time */
    ptime = time((time_t *) 0);
    if (ptime == ((time_t) - 1))
    {
        stamp[0] = '\0';
        return ERROR;
    }

    /* Convert the current time to local time */
    ltime = localtime(&ptime);
    if (ltime == NULL)
    {
        stamp[0] = '\0';
        return ERROR;
    }

    /* Generate the timestamp */
    if (strftime(stamp, stampsize, format, ltime) == 0) 
    {
        stamp[0] = '\0';
        return ERROR;
    }

    return SUCCESS;
}

/*************************************************************************

NAME: log_message

PURPOSE: Outputs logging message to the output file pointer

RETURNS: None

**************************************************************************/
static void log_message 
(
    int log_level,            /* I: message level for input */
    const char *filename,     /* I: source code file name for input */
    int line_number,          /* I: source code line number for input */
    const char *format,       /* I: format string for message */
    va_list ap                /* I: format string variables */
) 
{
    char time_stamp[20];
    char temp_string[500];
    static const char *log_level_message[] = {"DEBUG", "INFO", "WARN", "ERROR"};

    /* if file_ptr is not set (ias_log_message is not called), stdout is used */
    if (file_ptr == NULL)      
        file_ptr = stdout;

    /* if pid is not set (ias_log_message is not called), getpid is called to
       get the current processor id */
    if (pid == 0)
        pid = getpid();

    /* if log_level is out of range, log_level is set to be 
       IAS_LOG_LEVEL_ERROR */
    if (log_level < IAS_LOG_LEVEL_DEBUG || log_level > IAS_LOG_LEVEL_ERROR)
    {
        log_level = IAS_LOG_LEVEL_ERROR;
        IAS_LOG_WARNING("Log_level is out of range, should be one of " 
                         "IAS_LOG_LEVEL_DEBUG, IAS_LOG_LEVEL_INFO, "
                         "IAS_LOG_LEVEL_WARN, and IAS_LOG_LEVEL_ERROR");
    }

    /* if log_level is high enough, output the message */
    if (log_level >= ias_log_message_level)
    {
        /* Set arg_ptr to beginning of list of optional arguments */
        vsnprintf(temp_string, sizeof(temp_string), format, ap);

        format_time(time_stamp, sizeof(time_stamp),"%F %H:%M:%S");
        fprintf(file_ptr, "%19s  %s  %7d %-20s  %6d  %s %s\n",
                time_stamp, program_name, pid, filename, 
                line_number, log_level_message[log_level], temp_string);     
    }
}

/*************************************************************************

NAME: ias_log_message

PURPOSE: Outputs logging message to the output file pointer

RETURNS: None

**************************************************************************/
void ias_log_message 
(
    int log_level,            /* I: message level for input */
    const char *filename,     /* I: source code file name for input */
    int line_number,          /* I: source code line number for input */
    const char *format, ...   /* I: format string for message */
) 
{
    /* If channels are enabled, ignore non-channel debug messages */
    if (channels_on && log_level == IAS_LOG_LEVEL_DEBUG)
        return;
    va_list arglist;
    va_start(arglist, format);
    log_message(log_level, filename, line_number, format, arglist);
    va_end(arglist);
}

/*************************************************************************

NAME: ias_log_message_with_channel

PURPOSE: Outputs logging message to the output file pointer if logging is
    enabled for the specified channel

RETURNS: None

**************************************************************************/
void ias_log_message_with_channel
(
    int log_level,            /* I: message level for input */
    const char *channel,      /* I: channel name */
    const char *filename,     /* I: source code file name for input */
    int line_number,          /* I: source code line number for input */
    const char *format, ...   /* I: format string for message */
)
{
    if (is_channel_enabled(channel))
    {
        va_list arglist;
        va_start(arglist, format);
        log_message(log_level, filename, line_number, format, arglist);
        va_end(arglist);
    }
}
 
/*************************************************************************/

