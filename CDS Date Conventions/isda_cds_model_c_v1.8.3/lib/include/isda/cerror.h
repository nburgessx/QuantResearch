/*
 * ISDA CDS Standard Model
 *
 * Copyright (C) 2009 International Swaps and Derivatives Association, Inc.
 * Developed and supported in collaboration with Markit
 * 
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the ISDA CDS Standard Model Public License.
 */

#ifndef CERROR_H
#define CERROR_H

#include "cgeneral.h"
#include <stdarg.h>
#include <stdio.h>

#ifdef __cplusplus
extern "C"
{
#endif


/*t
 */
typedef TBoolean (JpmcdsErrCallBackFunc)(char *string, void *callBackData);


/*f
***************************************************************************
** Writes an error message with a variable number of arguments to a log 
** file.
***************************************************************************
*/
void JpmcdsErrMsg(
    char *format,    /* (I) printf style format string. */
    ...)             /* (I) Variable arguments. */
#if __GNUC__ > 1
__attribute__ ((format(printf, 1, 2)))
#endif
;


/*f
***************************************************************************
** Same as JpmcdsErrMsg but using va_list for its arguments instead of ...
***************************************************************************
*/
void JpmcdsErrMsgV(
    char    *format,  /* (I) printf style format string */
    va_list args      /* (I) Variable arguments */
);


/*f
***************************************************************************
** Writes an error message to a log file.
***************************************************************************
*/
void JpmcdsErrLogWrite(char *message);


/*f
***************************************************************************
** Tells JpmcdsErrMsg to actually write something (either to
** log file, or to standard out, depending on how USE_PRINTF is set.)
***************************************************************************
*/
void JpmcdsErrMsgOn(void);


/*f
***************************************************************************
** Tells JpmcdsErrMsg NOT to actually write something (either to
** log file, or to standard out, depending on how USE_PRINTF is set.)
***************************************************************************
*/
void JpmcdsErrMsgOff(void);


/*f
***************************************************************************
** Writes the message "%s: Failed.\n" to the error log and returns FAILURE.
** This is to enable a single function call on failing within a function as
** follows: return (JpmcdsErrMsgFailure (routine)); 
***************************************************************************
*/
int JpmcdsErrMsgFailure
(char *routine /* (I) Name of routine published in error log. */
);


/*f
***************************************************************************
** Flushes the err message stream. This can be called from within a debugger
** to make sure that the latest output written to the error log is flushed.
***************************************************************************
*/
int JpmcdsErrMsgFlush(void);


/*f
***************************************************************************
** Returns the current logging state.
***************************************************************************
*/
EXPORT TBoolean JpmcdsErrMsgStatus(void);


/*f
***************************************************************************
** Tells JpmcdsErrMsg to write to the provided file name. Routine
** opens the file, and returns SUCCESS or FAILURE depending on whether the
* file was successfully opened or not.
***************************************************************************
*/
EXPORT int JpmcdsErrMsgFileName(
    char *fileName,  /* (I)  File name to use. */
    TBoolean append  /* (I)  Append flag. */
);


/*f
***************************************************************************
** Returns pointer to privately held error log file name which was set by 
** JpmcdsErrMsgFileName.
***************************************************************************
*/
EXPORT char* JpmcdsErrMsgGetFileName(void);


/*f
***************************************************************************
** Turns on the error message record facility.
***************************************************************************
*/
EXPORT char** JpmcdsErrGetMsgRecord(void);


/*f
***************************************************************************
** Turns on the error message record facility.
***************************************************************************
*/
EXPORT int JpmcdsErrMsgEnableRecord(
    int numberOfMessages,     /* (I) Number of messages to save. */
    int messageSize);         /* (I) Maximum size of each message. */


/*f
***************************************************************************
** Turns off the error message record facility.
***************************************************************************
*/
EXPORT int JpmcdsErrMsgDisableRecord(void);


/*f
***************************************************************************
** Returns whether the error message record has been enabled.
**
** Effectively error message record is enabled via the add-in function
** ErrorLog(1) and disabled by ErrorLog(0).
**
** Thus we can use this function to determine in a pure add-in based
** interface whether the user has requested error logging.
**
** This enables a special interface error handler to determine whether
** the user has requested a full error log.
***************************************************************************
*/
TBoolean JpmcdsErrMsgRecordEnabled(void);

/*f
 * Allows the user to supply a callback for error handling.
 *
 * After you add a callback, then any call to JpmcdsErrMsg which would otherwise
 * result in a message being written to file (not that if error messaging is
 * off then messages are not being written to file) will first invoke a call
 * to the callback routine.
 *
 * The callback routine returns TRUE to indicate that the message should still
 * be logged to file, and FALSE to indicate that the message should not be
 * logged to file.
 *
 * The callbackData provided in this function call will be made available
 * to the callback routine.
 */
EXPORT void  JpmcdsErrMsgAddCallback(
    JpmcdsErrCallBackFunc *userFunc,       /* (I) User callback function. */
    TBoolean               sendTimeStamp,  /* (I) Send timestamp flag. */
    void                  *callbackData    /* (I) User supplied data passed to cb.*/
);

/*f
 * Returns callback information currently in use.
 */
EXPORT void  JpmcdsErrMsgGetCallback
   (JpmcdsErrCallBackFunc **userFunc,        /* (O) */
    TBoolean               *sendTimeStamp,   /* (O) */
    void                  **callBackData);   /* (O) */


#ifdef __cplusplus
}
#endif

#endif    /* CERROR_H */

