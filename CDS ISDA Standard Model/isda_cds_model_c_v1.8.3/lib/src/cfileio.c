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
#include <stdarg.h>
#include <string.h>

#include "cgeneral.h"
#include "cmemory.h"
#include "cfileio.h"
#include "lprintf.h"
#include "cerror.h"
#include "macros.h"

#include <errno.h>

#include <sys/types.h>
#include <sys/stat.h>

#if defined(WIN32) /* Windows */
#include <io.h>
#include <wtypes.h>
#include <winbase.h>
#include <direct.h>
#elif defined(UNIX)
#include <dirent.h>
#include <unistd.h>
#endif


int JpmcdsLvfscanf(TFile *, char *, va_list);


/*
***************************************************************************
** Opens file.
***************************************************************************
*/
TFile* JpmcdsFopen(
    char *filename,   /* (I) Filename. */
    TMODE mode)    /* (I) File access type. */
{
     static char routine[] = "JpmcdsFopen";
     TFile   *tFile = NULL;

     if (filename == NULL || *filename == '\0')
     {
         JpmcdsErrMsg ("%s: Unspecified file name\n", routine);
         return NULL;
     }

     tFile = (TFile *) JpmcdsMallocSafe(sizeof(TFile));
     if (tFile == NULL)
     {
         return NULL;
     }

     tFile->type = TFILE_FILE;
     tFile->lastChar = 0;
     tFile->hasLastChar = FALSE;

     tFile->filename = NEW_ARRAY(char, strlen(filename)+2);
     if (tFile->filename == NULL)
     {
         FREE(tFile);
         return NULL;
     }
     strcpy(tFile->filename, filename);

     /*
      * Actually open the file as the last step. Makes sure to
      * clean up memory in case the file cannot be opened.
      */

     switch(mode)
     {
     case JPMCDS_FWRITE:
         tFile->fileP = (void *)fopen(filename, "w");
         break;
     case JPMCDS_APPEND:
     case JPMCDS_FAPPEND:
         tFile->fileP = (void *)fopen(filename, "a");
         break;
     case JPMCDS_FREAD:
         tFile->fileP = (void *)fopen(filename, "r");
         break;
     default:
         JpmcdsErrMsg("%s:  Unknown mode %d\n", routine, mode);
         FREE(tFile->filename);
         FREE(tFile);
         return NULL;
     }
     if (tFile->fileP == NULL)
     {
         JpmcdsErrMsg("%s:  Error using file %s, errno %d\n", routine, filename,
                   errno);
         FREE(tFile->filename);
         FREE(tFile);
         return NULL;
     }

     return tFile;
}


/*
***************************************************************************
** Closes file.
***************************************************************************
*/
int JpmcdsFclose(
    TFile *tFile)     /* (I) File handle. */
{
     static char routine[] = "JpmcdsFclose";
     int         status = FAILURE; /* Until successful */

     if (tFile == NULL)
         return SUCCESS;

     if (fclose((FILE *)tFile->fileP) != 0)
         goto done;

     status = SUCCESS;

 done:
     if (status == FAILURE)
     {
         JpmcdsErrMsg("%s:  Error closing file %s.\n", routine, tFile->filename);
     }
     FREE(tFile->filename);
     FREE(tFile);
     return status;
}


/*
***************************************************************************
** Flushes file buffer.
***************************************************************************
*/
int JpmcdsFflush(
    TFile *tFile)     /* (I) File handle. */
{
     static char   routine[] = "JpmcdsFflush";

     if (fflush((FILE *)tFile->fileP) != 0)
     {
         JpmcdsErrMsg("%s:  Error flushing file %s.\n", routine, tFile->filename);
         return FAILURE;
     }
     return SUCCESS;
}


/*
***************************************************************************
** Gets a string from a file.
***************************************************************************
*/
int JpmcdsFgets(
    char *buffer,     /* (I) Buffer to store the string. */
    int numBytes,     /* (I) The number of bytes to get. */
    TFile *tFile)     /* (I) File handle. */
{
     static char routine[] = "JpmcdsFgets";
     char *tmpBuffer;

     tmpBuffer = buffer;
     if (tFile->hasLastChar == TRUE) {
          buffer[0] = tFile->lastChar;
          tmpBuffer = &buffer[1];
          tFile->lastChar = 0;
          tFile->hasLastChar = FALSE;
          numBytes--;
     }

    if (fgets(tmpBuffer, numBytes, (FILE *)tFile->fileP) == NULL)
         return FAILURE;

    return SUCCESS;
}


/*
***************************************************************************
** Puts a string to a file.
***************************************************************************
*/
int JpmcdsFputs(
    char *output,     /* (I) Strings which contain data. */
    TFile *tFile)     /* (I) File handle. */
{
    static char   routine[] = "JpmcdsFputs";

    if (fputs(output, (FILE *)tFile->fileP) == EOF)
    {
        JpmcdsErrMsg("%s:  Error writing file %s.\n", routine, tFile->filename);
        return FAILURE;
    }
    return SUCCESS;
}


/*
***************************************************************************
** Reads characters from a file.
***************************************************************************
*/
int JpmcdsFread(
    char *ptr,          /* (I) Pointer to memory to store data. */
    int numItems,     /* (I) Number of bytes to read. */
    TFile *tFile)     /* (I) File handle. */
{
     static char   routine[] = "JpmcdsFread";
     int     numRead;

    if (tFile->hasLastChar == TRUE)
    {
        JpmcdsErrMsg("%s:  Error mixing read and scanf.\n", routine);
        return FAILURE;
    }

    numRead = fread(ptr, 1, numItems, (FILE *)tFile->fileP);
    if (numRead < 1)
         return FAILURE;

    return numRead;
}


/*
***************************************************************************
** Puts a string to a file.
***************************************************************************
*/
int JpmcdsFwrite(
    char *ptr,        /* (I) Pointer to data to write. */
    int numItems,     /* (I) Number of bytes to write. */
    TFile *tFile)     /* (I) File handle. */
{
     static char   routine[] = "JpmcdsFwrite";
     int     numWrite;

     numWrite = fwrite(ptr, 1, numItems, (FILE *)tFile->fileP);
     if (numWrite < 1)
     {
         JpmcdsErrMsg("%s:  Error writing file %s.\n", routine, tFile->filename);
         return FAILURE;
     }
     return numWrite;
}


/*
***************************************************************************
** Gets a character from a file.
***************************************************************************
*/
int JpmcdsFgetc(
    TFile *tFile)     /* (I) File handle. */
{
     int     c;
     char    tempChar;
     if (tFile->hasLastChar == TRUE)
     {
          tempChar = tFile->lastChar;
          tFile->lastChar = 0;
          tFile->hasLastChar = FALSE;
          return tempChar;
     }

    c = getc((FILE *)tFile->fileP);
    if (c == -1)
    {
        return FAILURE;
    }

    return c;
}


/*
***************************************************************************
** Puts a character to a file.
***************************************************************************
*/
int JpmcdsFputc(
    char c,     /* (I) Character to write. */
    TFile *tFile)     /* (I) File handle. */
{
     static char   routine[] = "JpmcdsFputc";
     int     rc;

     rc = putc(c, (FILE *)tFile->fileP);
     if (rc == EOF)
     {
         JpmcdsErrMsg("%s:  Error writing file %s.\n", routine, tFile->filename);
         return FAILURE;
     }

     return rc;
}


/*
***************************************************************************
** Prints variable and format information to a string.
***************************************************************************
*/
int JpmcdsSprintf(
    char *buffer,     /* (I) Buffer to write data to. */
    int bufLen,       /* (I) Maximum length of the buffer. */
    char *fmt,     /* (I) Format string. */
    ...)        /* (I) List of arguments for format string. */
{
     va_list args;
     TFile   tFile;
     int     bytesSent;

     tFile.type = TFILE_STRING;
     tFile.charPtr = buffer;
     tFile.size = bufLen;
     tFile.curSize = 0;
     tFile.lastChar = 0;
     tFile.hasLastChar = FALSE;

     va_start(args, fmt);
     bytesSent = JpmcdsVfprintf(&tFile, fmt, args);
     va_end(args);

     if (bytesSent < 1)
          return FAILURE;

     return bytesSent;
}


/*
***************************************************************************
** Get information from a string.
***************************************************************************
*/
int JpmcdsSscanf(
    char *buffer,     /* (I) File handle. */
    char *fmt,     /* (I) Format string. */
    ...)        /* (I) List of arguments for format string. */
{
     va_list args;
     TFile   tFile;
     int     numParse;

     tFile.type = TFILE_STRING;
     tFile.charPtr = buffer;
     tFile.curSize = 0;
     tFile.lastChar = 0;
     tFile.hasLastChar = FALSE;

     va_start(args, fmt);
     numParse = JpmcdsLvfscanf(&tFile, fmt, args);
     va_end(args);

     if (numParse == EOF)
          return FAILURE;  /* Note that EOF=FAILURE */

     return numParse;
}
