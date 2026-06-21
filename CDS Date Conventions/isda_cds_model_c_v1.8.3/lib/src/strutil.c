/*
 * ISDA CDS Standard Model
 *
 * Copyright (C) 2009 International Swaps and Derivatives Association, Inc.
 * Developed and supported in collaboration with Markit
 * 
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the ISDA CDS Standard Model Public License.
 */

#include "strutil.h"
#include <ctype.h>
#include <string.h>
#include <stdarg.h>
#include "cgeneral.h"
#include "cerror.h"
#include "cmemory.h"
#include "macros.h"


#define TOUPPER(_c) ( (islower(_c)) ? COND_TOUPPER(_c) : (_c) )
#define COND_TOUPPER(_c)    ( (_c)-'a'+'A' )


/*
***************************************************************************
** Duplicates a string using our memory allocation.
***************************************************************************
*/
char* JpmcdsStringDuplicate(          /* duplicate a string */
   char    *in)                            /* (I) input string to duplicate */
{
   static char routine[] = "JpmcdsStringDuplicate";

   char   *out;                        /* duplicated string */
   size_t  len;                        /* length of string we want duped */

   if (in == NULL)
       return NULL;

   len = strlen(in) + 1;
   out = NEW_ARRAY(char,len);

   if (out == NULL)
   {
       JpmcdsErrMsg("%s: out of memory\n",routine);
       return NULL;
   }

   (void) strcpy(out,in);
   return out;
}


/*
***************************************************************************
** Removes any initial whitespace from an input string and converts all 
** characters to uppercase.  Does not allocate any memory.
**
** Needed for efficiency reasons when parsing user input.
***************************************************************************
*/
int JpmcdsStringPreprocess
(
 char   *inputString,    /* (I) String to preprocess. */
 int     n,              /* (I) Max. # of non-whitespace chars to examine. */
 char   *outputString    /* (O) Uppercase w/ no initial whitespace. */
)
{
    int i;

    if( inputString == NULL ||
        outputString == NULL )
    {
        JpmcdsErrMsg( "JpmcdsStringPreprocess: Input or output string is NULL.\n" );
        return FAILURE;
    }

    /* Initialize the output buffer. */
    (void)memset( outputString, 0, n+1 );

    /* Remove initial whitespace. */
    while( *inputString && isspace((int)*inputString) )
    {
        inputString++;
    }
    if( *inputString == 0 )
    {
        /* The output string will be empty, which is OK. */
        return SUCCESS;
    }

    /* Construct outputString from inputString. */
    (void)strncpy( outputString, inputString, n );
    for( i=0; i<n && outputString[i]; i++ )
    {
        outputString[i] = (char)TOUPPER( (int)outputString[i] );
    }

    return SUCCESS;
}


/*
***************************************************************************
** Simple string parser.
**
** The idea is that the input string is parsed into its constituent
** components using a number of delimiters.
**
** If a delimiter is missing in the input string, then this is not an
** error. Instead all remaining constituent strings are returned as NULL.
**
** You need to provide the same number of pointers to char* as the number
** delimiters provided.
**
** The input string is modified - all the delimiters found are replaced
** with the '\0' to terminate the C-string, and the pointers returned are
** within the input string. No new memory is assigned.
**
** Therefore if you wish to preserve the input string, then you should
** take a copy of it before calling this function.
**
** Example:
**   Suppose you have a comma delimited string with three constituents.
**   Your code would be as follows (ignoring error checking):
**
**   char *str1;
**   char *str2;
**   JpmcdsStringParser (str, ",,", &str1, &str2);
**
** Suppose str = "abc,def,ghi"
**
** Then after this call we would have:
**   str = "abc"
**   str1 = "def"
**   str2 = "ghi"
**
***************************************************************************
*/
int JpmcdsStringParser
(char *str,        /* (I/O) Input string - modified by this function */
 char *delimiters, /* (I) Delimiters */
 ...               /* (O) Need to provide strlen(delimiters) char** */
)
{
    static char routine[] = "JpmcdsStringParser";
    int         status    = FAILURE;

    va_list     ap;
    int         len;
    int         idx;
    char      **subString;
    char       *buf;        /* Current pointer */
    char       *delimiter;

    va_start (ap, delimiters);

    if (str == NULL || delimiters == NULL)
    {
        JpmcdsErrMsg ("%s: NULL inputs.\n", routine);
        goto done;
    }

    len = strlen(delimiters);

    buf = str;
    for (idx = 0; idx < len; ++idx)
    {
        subString = va_arg (ap, char**);
        if (buf != NULL)
        {
            delimiter = strchr (buf, delimiters[idx]);
            if (delimiter == NULL)
            {
                buf = NULL;
            }
            else
            {
                *delimiter = '\0'; /* Terminate previous component */
                buf = delimiter + 1;
            }
        }
        *subString = buf;
    }

    status = SUCCESS;

done:

    va_end (ap);

    if (status != SUCCESS)
        JpmcdsErrMsgFailure (routine);

    return status;
}


/*
***************************************************************************
** Splits a string into components given a separator.
**
** For example, bond.type with a separator of '.' would be split into
** "bond" and "type".
**
** Returns the number of components, and an array of strings. You must free
** the array of strings using FREE() or JpmcdsFreeSafe().
**
** Note that the split array returned is an array of (numItems+1) with the
** final item being NULL. This means that you can iterate this array either
** using numItems or by testing the string pointer.
***************************************************************************
*/
int JpmcdsStringSplit
(char   *str,       /* (I) Input string */
 char    separator, /* (I) Separator - note this is a character not a string */
 size_t *numItems,  /* (O) Number of items found */
 char ***split      /* (O) Array of strings found within input. User must
                      FREE() this once they have finished with it. */
)
{
    static char routine[] = "JpmcdsStringSplit";
    int         status    = FAILURE;

    size_t      myNumItems = 0;
    size_t      len;
    char        c;
    char      **mySplit = NULL;
    char       *myString;
    char       *ptr;
    size_t      pos;

    if (split != NULL) *split = NULL;
    if (numItems != NULL) *numItems = 0;
    if (str == NULL || numItems == NULL || split == NULL)
    {
        JpmcdsErrMsg ("%s: NULL inputs\n", routine);
        goto done;
    }

    len = 0;
    myNumItems = 1;
    ptr = str;
    while ((c = *ptr) != '\0')
    {
        ++len;
        if (c == separator)
            ++myNumItems;
        ++ptr;
    }

    /*
    ** mySplit allocates all the memory that we need.
    ** This consists of the array of char*, plus a copy of the input string.
    ** Within the input string, we will replace instances of separator with
    ** the end of the string.
    */
    mySplit  = (char**)(JpmcdsMallocSafe(sizeof(char*)*(myNumItems+1)+len+1));
    if (mySplit == NULL)
        goto done;
    myString = ((char*)mySplit) + sizeof(char*)*(myNumItems+1);
    strcpy (myString, str);

    ptr = myString;
    pos = 0;

    mySplit[pos] = ptr;
    while (*ptr != '\0')
    {
        if (*ptr == separator)
        {
            *ptr = '\0';
            ++pos;
            mySplit[pos] = ptr+1;
        }
        ++ptr;
    }

    ASSERT(pos+1 == myNumItems);

    *split    = mySplit;
    *numItems = myNumItems;
    mySplit   = NULL;
    status    = SUCCESS;

 done:

    FREE(mySplit);
    if (status != SUCCESS)
        JpmcdsErrMsgFailure (routine);
    return status;
}
