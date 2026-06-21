/*
 * ISDA CDS Standard Model
 *
 * Copyright (C) 2009 International Swaps and Derivatives Association, Inc.
 * Developed and supported in collaboration with Markit
 * 
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the ISDA CDS Standard Model Public License.
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <ctype.h>

#include "cgeneral.h"
#include "cerror.h"
#include "cmemory.h"
#include "cfileio.h" 
#include "macros.h"

#define JPMCDS_ERR_MSG_BUFFER 4096

#if defined(UNIX)
#define DEF_LOG_FILE_PATH "error.log"
#elif defined(DOS)
#define DEF_LOG_FILE_PATH "c:\\error.log"
#else                                   /* DOS is default */
#define DEF_LOG_FILE_PATH "c:\\error.log"
#endif

#define ENV_LOG_PATH "JPMCDS_LOG_FILE"     /* Name of environment variable*/


#define TS_SENSITIVITY 2.0      /* min. interval in sec for time stamps */

static TFile *pFp = NULL;

/* pAppendOnOpen indicates whether the error log should be opened in
 * append mode. If a new empty file is to be opened, this should be
 * set to FALSE.
 */
static TBoolean pAppendOnOpen = FALSE;


/* pFileNameSet indicates whether the file name has been set or not.
 */
static TBoolean pFileNameSet = FALSE;   


/* pLogFilePath contains the current error log file name.
 */
#define MAX_LOG_PATH_LEN 256
static char pLogFilePath[MAX_LOG_PATH_LEN];

/* pWriteMessage is toggled by JpmcdsErrMsgOn and JpmcdsErrMsgOff. It simply
 * determines whether a call to JpmcdsErrMsg has any effect or not.
 */
static TBoolean pWriteMessage = FALSE;  /* Default is OFF */

static JpmcdsErrCallBackFunc *errorUserFunc = NULL;
static TBoolean          errorSendTimeStamp = TRUE;
static void             *errorCallBackData;

static int JpmcdsWriteToLog(TBoolean formatted, char *format, va_list parminfo);
static void JpmcdsAddToRecord(char *buffer);


/*
***************************************************************************
** Sets the error log file name.
***************************************************************************
*/
static void SetFileName(char *fileName);


/*
***************************************************************************
** Gets the error log file name.
***************************************************************************
*/
static char *GetFileName(void);


/*
***************************************************************************
** Gets the default error log file name.
***************************************************************************
*/
static char *GetDefaultFileName(void);


/*
***************************************************************************
** Opens file name.
***************************************************************************
*/
static int FileCreate(char *fileName, TBoolean append);

/* TimeStamp routines.
 */
static int TimeStampFill(time_t aStamp_t);
static int TimeStampRequired(void);

typedef struct record_s 
{
   TBoolean   on;
   int        number;
   int        size;
   int        lastLine;
   char       *alloc;
   char       **buf;
} Record;

static Record record = {FALSE, 0, 0, 0, NULL, NULL};


/*
***************************************************************************
** Tells JpmcdsErrMsg to actually write something (either to log file or to
** standard out, depending on how USE_PRINTF is set.)
***************************************************************************
*/
void JpmcdsErrMsgOn(void)
{
    pWriteMessage = TRUE;
}


/*
***************************************************************************
** Tells JpmcdsErrMsg to not write anything.
***************************************************************************
*/
void JpmcdsErrMsgOff(void)
{
    pWriteMessage = FALSE;
}


/*
***************************************************************************
** Returns the current logging state.
***************************************************************************
*/
EXPORT TBoolean JpmcdsErrMsgStatus(void)
{
    return pWriteMessage;
}


/*
***************************************************************************
** Tells JpmcdsErrMsg to write to the provided file name. Routine opens the 
** file, and returns SUCCESS or FAILURE depending on whether the file was 
** successfully opened or not.
***************************************************************************
*/
EXPORT int JpmcdsErrMsgFileName(char *fileName, TBoolean append)
{
    static char routine[] = "JpmcdsErrMsgFileName";

    if (fileName == (char *)NULL)
        fileName = GetDefaultFileName();

#if defined(UNIX)
    else if( strchr(fileName, '~') != NULL )
    {
        /* The '~' character is specific to the C shell.  Most users
           are too familiar with the C shell to know this. */
        JpmcdsErrMsg( "%s:  File name must not contain the '~' character.\n", routine );
        return FAILURE;
    }
#endif

    /* Make sure we can open the file.
     */
    if (FileCreate(fileName, append) == SUCCESS)
    {
        /* Now close it.
         */
        JpmcdsFclose(pFp);
        pFp = NULL;

        /* If we're appending, we don't want to reopen the file a fresh
         * the next time we call JpmcdsErrMsg. Otherwise, we can reopen it
         * when we call JpmcdsErrMsg the first time.
         */
        pAppendOnOpen = append;

        /* Now that we have successfully opened and closed the file,
         * save the file name for next time.
         */
        SetFileName(fileName);
        return SUCCESS;
    }
    else
    {
        /* Comes here on failure.
         */
        JpmcdsErrMsg("%s: Failed to open file \"%s\".\n", routine, fileName);
        return FAILURE;
    }
}


/*
***************************************************************************
** Returns pointer to privately held error log file name which was set by 
** JpmcdsErrMsgFileName.
***************************************************************************
*/
EXPORT char* JpmcdsErrMsgGetFileName(void)
{
    return GetFileName();
}


/*
***************************************************************************
** Writes an error message with a variable number of arguments to a log file.  
***************************************************************************
*/
void JpmcdsErrMsg(char *format,  ...)
{
    va_list parminfo;
    va_start(parminfo, format);
    JpmcdsErrMsgV (format, parminfo);
    va_end(parminfo);
    return;
}


/*
***************************************************************************
** Based on JpmcdsErrMsg except that this version takes va_list as input 
** instead of ...
**
** Writes an error message with a variable number of arguments to a log file
***************************************************************************
*/
void JpmcdsErrMsgV(char *format,  va_list parminfo)
{
    TBoolean localpWriteMessage;     /* Local copy of pWriteMessage */

    /*
     * We MUST turn off error logging while JpmcdsErrMsg is executing.
     * If we don't, we could potentially land in an infinite loop.
     */
    localpWriteMessage = pWriteMessage;

    if (!localpWriteMessage)
    {
        return;                         /* Message writing not turned on */
    }

    pWriteMessage = FALSE; /* Turn off message writing when executing
                              JpmcdsErrMsg because we don't want an
                              infinite loop of calls to JpmcdsErrMsg
                           */

    if (pFp == NULL)
    {
        char *fileName = GetFileName();

        if (FileCreate(fileName, pAppendOnOpen) == FAILURE)
            goto done;

        /* Just opened file from scratch. From now on, open on append.
         * Save file name for next time.
         */
        if (!pAppendOnOpen)
        {
            pAppendOnOpen = TRUE;
            SetFileName(fileName);
        }
    }
    
    if (TimeStampRequired() != SUCCESS)
        goto done;

    if (JpmcdsWriteToLog(TRUE, format, parminfo) != SUCCESS)
    {
        pWriteMessage = localpWriteMessage; /* Restore the variable */
        return; /* failed */
    }

    /* Close file in between calls, unless user supplied file pointer to us.
     * Note that if pFp is NULL, JpmcdsFClose does nothing.
     */
 done:
    JpmcdsFclose(pFp);
    pFp = NULL;
    pWriteMessage = localpWriteMessage; /* Restore the variable */
    return;
}


/*
***************************************************************************
** Writes an error message to a log file..
***************************************************************************
*/
void JpmcdsErrLogWrite(char *message)
{
    if (!pWriteMessage)
        return;  /* Message writing not turned on */

    if (pFp == NULL)
    {
        char *fileName = GetFileName();

        if (FileCreate(fileName, pAppendOnOpen) == FAILURE)
            goto done;

        /* Just opened file from scratch. From now on, open on append.
         * Save file name for next time.
         */
        if (!pAppendOnOpen)
        {
            pAppendOnOpen = TRUE;
            SetFileName(fileName);
        }
    }
    
    if (TimeStampRequired() != SUCCESS)
        goto done;

    if (JpmcdsWriteToLog(FALSE, message, NULL) != SUCCESS)
        return; /* failed */

    /* Close file in between calls, unless user supplied file pointer to us.
     * Note that if pFp is NULL, JpmcdsFClose does nothing.
     */
 done:
    JpmcdsFclose(pFp);
    pFp = NULL;
    return;
}
 
 
/*
***************************************************************************
** Writes the message "%s: Failed.\n" to the error log and returns FAILURE.
** This is to enable a single function call on failing within a function as
** follows: return (JpmcdsErrMsgFailure (routine)); 
***************************************************************************
*/
int JpmcdsErrMsgFailure
(char *routine /* (I) Name of routine published in error log. */
)
{
    JpmcdsErrMsg ("%s: Failed.\n", routine);
    return FAILURE;
}


/*
***************************************************************************
** Flushes the err message stream. This can be called from within a debugger
** to make sure that the latest output written to the error log is flushed. 
***************************************************************************
*/
int JpmcdsErrMsgFlush(void)
{
    if (pFp)
    {
        if (JpmcdsFflush(pFp) != SUCCESS)
            return JpmcdsErrMsgFailure("JpmcdsErrMsgFlush");
    }
    return SUCCESS;
}


/*
***************************************************************************
** Generic routine to open a file.  Note that if pFp is not NULL, it will 
** be overwritten if this routine is successful.
***************************************************************************
*/
static int FileCreate
     (char    *fileName,                /* (I) */
      TBoolean append)                  /* (I) */
{
    TMODE   tMode;
    TFile  *tfp = NULL;              /* Temporary */

    if (append == TRUE)
    {
        tMode = JPMCDS_APPEND;
    }
    else
    {
        tMode = JPMCDS_FWRITE;
    }


    /* Open the file. Only overwrite pFp if the open operation is successful.
     * If the open operation fails when we're trying to append, it may 
     * mean that the file has been deleted.  In this case, we try again,
     * using JPMCDS_FWRITE instead of JPMCDS_APPEND.
     */
    tfp = JpmcdsFopen(fileName, tMode);
    if (tfp == NULL)
    {
        if (append)             /* Now try again, using JPMCDS_FWRITE */
        {
            tfp = JpmcdsFopen(fileName, JPMCDS_FWRITE);
            if (tfp == NULL)
                return FAILURE; /* Now we finally give up */
        }
        else
        {
            return FAILURE;
        }
    }

    /* Close old file before resetting pFp. Otherwise memory is lost.
     */
    if (pFp != NULL)
        JpmcdsFclose(pFp);

    /* Here we have successfully opened the file.
     */
    pFp = tfp;

    return SUCCESS;
}


/*
***************************************************************************
** Determines if a timestamp should be output with the current error message.
**
** If timestamp required, it gets created in the beginning of the  buffer 
** argument.
***************************************************************************
*/
static int TimeStampRequired (void)
{
    static time_t lastStamp_t =0;
    time_t currStamp_t;

    time(&currStamp_t);

    /*
     *  Removed difftime so that it will work with gcc.
     *  Straight subtraction should be ok in this case.  Any
     *  floating point imprecision will simply mean that it
     *  will take a little longer for the timestamp to 
     *  appear.
     */
    if ((currStamp_t - lastStamp_t) > TS_SENSITIVITY)
    {
        lastStamp_t = currStamp_t;
        if (TimeStampFill(currStamp_t) != SUCCESS)
        {
            return FAILURE;
        }
    }

    return SUCCESS;
}



enum offsets 
{ /* offsets for starting and end delimiting characters of ctime() output */
    S_MONTH =  4,
    E_MONTH =  7,
    S_YEAR  = 20, 
    E_YEAR  = 24, 
    S_DAY   =  8,
    E_DAY   = 10,
    S_TIME  = 11,
    E_TIME  = 19
};


/*
***************************************************************************
** Writes the timestamp to the file.
***************************************************************************
*/
static int TimeStampFill(time_t aStamp_t)
{
    char buffer[JPMCDS_ERR_MSG_BUFFER];

   /*  Format of the string returned by ctime() with offsets
    *  0  3   7 9          2   4
    *  Wed Jan 02 02:03:55 1980\n\0
    */
    char *sb;  /* buffer containing ctime() output */

    /* first line of time stamp */
    if (JpmcdsWriteToLog(FALSE, "\n--------------------------\n", NULL) != SUCCESS)
        return FAILURE;

    /* second line of time stamp */
    if ((sb = ctime(&aStamp_t)) != NULL)
    {
        sb[E_MONTH] = sb[E_YEAR] = sb[E_DAY] = sb[E_TIME] = '\0';
        if (JpmcdsSprintf(buffer, JPMCDS_ERR_MSG_BUFFER,
            "-- %s-%s-%s %s --\n", sb + S_DAY, sb + S_MONTH,
            sb + S_YEAR, sb + S_TIME) == FAILURE)
        {
            return FAILURE;
        }
        if (JpmcdsWriteToLog(FALSE, buffer, NULL) != SUCCESS)
            return FAILURE;
    }
    else
    {
        if (JpmcdsWriteToLog(FALSE, "-- **** WRONG TIME **** --\n", NULL) != SUCCESS)
            return FAILURE;
    }                                                

    if (JpmcdsWriteToLog(FALSE, "--------------------------\n", NULL) != SUCCESS)
        return FAILURE;

    return SUCCESS;
}


/*
***************************************************************************
** Writes an error message to the log or to the user function.
***************************************************************************
*/
static int JpmcdsWriteToLog(TBoolean formatted, char *format, va_list parminfo)
{
    char buffer[JPMCDS_ERR_MSG_BUFFER];
    char *bufp = buffer; 

    if (formatted)
    {
        vsprintf(buffer, format, parminfo);
    }
    else
    {
        bufp = format;
    }

    if (record.on == TRUE)
    {
        JpmcdsAddToRecord(bufp);
    }

    if (errorUserFunc == NULL)
    { 
        if (JpmcdsFputs(bufp, pFp) == FAILURE)
        {
            goto done;
        }
    }
    else
    {
        if (errorUserFunc(bufp, errorCallBackData) == TRUE)
        {
            if (JpmcdsFputs(bufp, pFp) == FAILURE)
            {
                goto done;
            }
        }
    }

    return SUCCESS;

 done:
    JpmcdsFclose(pFp);
    pFp = NULL;
    return FAILURE;
}


/*
***************************************************************************
** Turns on the error message record facility.
***************************************************************************
*/
EXPORT int JpmcdsErrMsgEnableRecord(
    int numberOfMessages,   /* (I) Number of messages to save. */
    int messageSize)         /* (I) Maximum size of each message. */
{
    int   i;
    char *all;

    if (record.on)         /* Already on! */
       return SUCCESS;

    record.on = TRUE;
    record.number = numberOfMessages;
    record.size = messageSize;
    record.lastLine = -1;
    /*
    ** The extra pointer is allocated and always set to NULL to provide a "NULL"
    ** terminator of the array of pointers in JpmcdsErrGetMsgRecord
    */
    record.buf = NEW_ARRAY (char*, record.number+1);
    if (record.buf == NULL)
        return FAILURE;

    all = NEW_ARRAY (char, record.number * messageSize);
    if (all == NULL)
        return FAILURE;
    record.alloc = all;
    for (i = 0; i < record.number; i++)
    {
        record.buf[i] = all;
        all += messageSize;
    }
    record.buf[record.number] = NULL; /* not strictly necessary! */

    return SUCCESS;
}


/*
***************************************************************************
** Turns off the error message record facility.
***************************************************************************
*/
EXPORT int JpmcdsErrMsgDisableRecord(void)
{
    if (record.on)                    /* only bother if recording is on */
    {    
        record.on = FALSE;

        FREE (record.alloc);
        record.alloc = NULL; 

        FREE(record.buf);
        record.buf = NULL;
    } 

    return SUCCESS;
}


/*
***************************************************************************
** Returns whether the error message record has been enabled.
**
** Effectively error message record is enabled via the add-in function
** ERR_LOG(1) and disabled by ERR_LOG(0)
**
** Thus we can use this function to determine in a pure add-in based
** interface whether the user has requested error logging.
**
** This enables a special interface error handler to determine whether
** the user has requested a full error log.
***************************************************************************
*/
TBoolean JpmcdsErrMsgRecordEnabled(void)
{
    return (record.on);
}


/*
***************************************************************************
** Turns on the error message record facility.
***************************************************************************
*/
EXPORT char** JpmcdsErrGetMsgRecord(void)
{
    static char routine[] = "JpmcdsErrGetMsgRecord";
    char   **temp;
    int    i;
    int    j;

    if (record.on == FALSE)
        return NULL;

    /*
     *  Determine if there are any error messages yet.
     */
    if (record.lastLine == -1)
        return record.buf;

    /*
     *  First rearrange output so that it is ordered right.
     */
    temp = (char **) JpmcdsMallocSafe((record.number + 1) * sizeof(char*));
    if (temp == NULL)
    {
        JpmcdsErrMsg("%s: JpmcdsMallocSafe failed.", routine);
        return NULL;
    }

    j = 0;
    for (i = record.lastLine + 1; i < record.number; i++)
    {
        temp[j++] = record.buf[i];
    }

    for (i = 0; i <= record.lastLine; i++)
    {
        temp[j++] = record.buf[i];
    }

    record.lastLine = j;
    if (record.lastLine >= record.number)
    {
        record.lastLine = record.number - 1; 
    }

    FREE(record.buf);
    record.buf = temp;
    
    return record.buf;
}


/*
***************************************************************************
** Adds a string to the record list.
***************************************************************************
*/
static void JpmcdsAddToRecord(char *buffer)
{
    int bufLen;

    bufLen = strlen(buffer);
    record.lastLine++;
    if (record.lastLine == record.number)
    {
        record.lastLine = 0; 
    }

    if (bufLen < record.size)
    {
        strcpy (record.buf[record.lastLine], buffer);
    }
    else
    {
        strncpy (record.buf[record.lastLine], buffer, record.size-4);
        strcpy (record.buf[record.lastLine]+record.size-4, "...");
    }
}


/*
***************************************************************************
** Sets the error log file name.
***************************************************************************
*/
static void SetFileName(char *fileName)
{
    pFileNameSet = TRUE;
    /* SetFileName may be passed the result from GetFileName.
     * Here we ensure that this can be done safely */
    if (pLogFilePath != fileName)
    {
        strcpy(pLogFilePath, /* <- */ fileName);
    }
}


/*
***************************************************************************
** Gets the error log file name.
***************************************************************************
*/
static char *GetFileName(void)
{
    if (pFileNameSet)
        return pLogFilePath;
    else
    {
        SetFileName(GetDefaultFileName());
        return pLogFilePath;
    }
}


/*
***************************************************************************
** Gets the default file name, and writes into the memory pointed to by 
** filename.
***************************************************************************
*/
static char* GetDefaultFileName(void)
{
    static char logFilePath[MAX_LOG_PATH_LEN];

    /* Try reading environment variable to get log file path name. */
    char *envp = getenv(ENV_LOG_PATH);
#ifdef UNIX
    char *home = getenv("HOME");
#endif

    if (envp == NULL)       /* If environ. var not set, use def. */
    {
#ifdef UNIX
        if (home != NULL)
        {
            strcpy(logFilePath, home);
            strcpy(&logFilePath[strlen(logFilePath)], "/");
            strcpy(&logFilePath[strlen(logFilePath)], DEF_LOG_FILE_PATH);
        }
        else
        {
#endif
            strcpy(logFilePath, DEF_LOG_FILE_PATH);
#ifdef UNIX
        }
#endif
    }
    else
        strcpy(logFilePath, envp);

    return logFilePath;
}

/*
***************************************************************************
** Allows the user to supply a callback.
***************************************************************************
*/
EXPORT void  JpmcdsErrMsgAddCallback
   (JpmcdsErrCallBackFunc *userFunc, 
    TBoolean               sendTimeStamp,
    void                  *callBackData)
{
    errorUserFunc = userFunc;
    errorSendTimeStamp = sendTimeStamp;
    errorCallBackData = callBackData;
}

/*
***************************************************************************
** Returns callback information currently in use.
***************************************************************************
*/
EXPORT void JpmcdsErrMsgGetCallback
   (JpmcdsErrCallBackFunc **userFunc,        /* (O) */
    TBoolean               *sendTimeStamp,   /* (O) */
    void                  **callBackData)    /* (O) */
{
    *userFunc = errorUserFunc;
    *sendTimeStamp = errorSendTimeStamp;
    *callBackData = errorCallBackData;
}
