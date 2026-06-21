/*
 * ISDA CDS Standard Model
 *
 * Copyright (C) 2009 International Swaps and Derivatives Association, Inc.
 * Developed and supported in collaboration with Markit
 * 
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the ISDA CDS Standard Model Public License.
 */

#ifndef CFILEIO_H
#define CFILEIO_H

#include <stdarg.h>
#include "cgeneral.h"

#ifdef __cplusplus
extern "C"
{
#endif
typedef struct
{
    void    *fileP;
    char    *filename;  
    int      type;
    char    *charPtr;
    int      size;
    int      curSize;
    int      lastChar;                  /* Make this int instead of char */
    TBoolean hasLastChar;               /* to avoid problems with char/int */
    int      stringBufSize;             /* conversion. */
} TFile;  




#ifdef __cplusplus
}
#endif

/* lprintf.h needs to use a TFile and so also includes this file. Therefore
 * define TFile before including lprintf.h
 */
#include "lprintf.h"

#ifdef __cplusplus
extern "C"
{
#endif


#define TFILE_STRING 1
#define TFILE_FILE   2

/*
 *  Added JPMCDS_FAPPEND.  Trying to phase out JPMCDS_APPEND.
 */
typedef enum {
    JPMCDS_FWRITE = 1,
    JPMCDS_APPEND = 2,
    JPMCDS_FREAD = 3,
    JPMCDS_FAPPEND =4
} TMODE;

#ifndef JPMCDS_PATH_SEPARATOR

#if defined(UNIX)
#define JPMCDS_PATH_SEPARATOR    "/"
#else 
/* Assumed DOS or NT */
#define JPMCDS_PATH_SEPARATOR     "\\"
#endif
 
#endif /* JPMCDS_PATH_SEPARATOR */

#define JPMCDS_MAX_PATH  256


/*f
***************************************************************************
** Opens file.
***************************************************************************
*/
TFile* JpmcdsFopen(
    char *filename,   /* (I) Filename. */
    TMODE mode);      /* (I) File access type. */


/*f
***************************************************************************
** Closes file.
***************************************************************************
*/
int JpmcdsFclose(
    TFile *tFile);    /* (I) File handle. */


/*f
***************************************************************************
** Flushes file buffer.
***************************************************************************
*/
int JpmcdsFflush(
    TFile *tFile);    /* (I) File handle. */


/*f
***************************************************************************
** Gets a string from a file.
***************************************************************************
*/
int JpmcdsFgets(
    char  *buffer,       /* (I) Buffer to store the string. */
    int    numBytes,     /* (I) The number of bytes to get. */
    TFile *tFile);       /* (I) File handle. */


/*f
***************************************************************************
** Puts a string to a file.
***************************************************************************
*/
int JpmcdsFputs(
    char  *output,       /* (I) Strings which contain data. */
    TFile *tFile);       /* (I) File handle. */


/*f
***************************************************************************
** Reads characters from a file.
***************************************************************************
*/
int JpmcdsFread(
    char  *ptr,          /* (I) Pointer to memory to store data. */
    int    numItems,     /* (I) Number of bytes to read. */
    TFile *tFile);       /* (I) File handle. */


/*f
***************************************************************************
** Puts a string to a file.
***************************************************************************
*/
int JpmcdsFwrite(
    char  *ptr,          /* (I) Pointer to data to write. */
    int    numItems,     /* (I) Number of bytes to write. */
    TFile *tFile);       /* (I) File handle. */


/*f
***************************************************************************
** Gets a character from a file.
***************************************************************************
*/
int JpmcdsFgetc(
    TFile *tFile);       /* (I) File handle. */


/*f
***************************************************************************
** Puts a character to a file.
***************************************************************************
*/
int JpmcdsFputc(
    char   c,            /* (I) Character to write. */
    TFile *tFile);       /* (I) File handle. */


/*f
***************************************************************************
** Prints variable and format information to a string.
***************************************************************************
*/
int JpmcdsSprintf(
    char *buffer,        /* (I) Buffer to write data to. */
    int   bufLen,        /* (I) Maximum length of the buffer. */
    char *fmt,           /* (I) Format string. */
    ...)                 /* (I) List of arguments for format string. */
#if __GNUC__ > 1
__attribute__ ((format(printf, 3, 4)))
#endif
;


/*f
***************************************************************************
** Get information from a string.
***************************************************************************
*/
int JpmcdsSscanf(
    char *buffer,        /* (I) File handle. */
    char *fmt,           /* (I) Format string. */
    ...)                 /* (I) List of arguments for format string. */
#if __GNUC__ > 1
__attribute__ ((format(scanf, 2, 3)))
#endif
;


#ifdef __cplusplus
}
#endif

#endif    /* CFILEIO_H */

