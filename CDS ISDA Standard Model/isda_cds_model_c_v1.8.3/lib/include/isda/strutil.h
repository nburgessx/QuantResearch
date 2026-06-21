/*
 * ISDA CDS Standard Model
 *
 * Copyright (C) 2009 International Swaps and Derivatives Association, Inc.
 * Developed and supported in collaboration with Markit
 * 
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the ISDA CDS Standard Model Public License.
 */

#ifndef STRUTIL_H
#define STRUTIL_H

#include "cgeneral.h"

#ifdef __cplusplus
extern "C"
{
#endif


/*f
***************************************************************************
** Duplicates a string using our memory allocation.
***************************************************************************
*/
char* JpmcdsStringDuplicate(
   char    *in);          /* (I) input string to duplicate */


/*f
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
);


/*f
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
** 
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
);


/*f
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
);


#ifdef __cplusplus
}
#endif

#endif    /* STRUTIL_H */




