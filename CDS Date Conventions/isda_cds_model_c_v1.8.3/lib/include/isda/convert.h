/*
 * ISDA CDS Standard Model
 *
 * Copyright (C) 2009 International Swaps and Derivatives Association, Inc.
 * Developed and supported in collaboration with Markit
 * 
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the ISDA CDS Standard Model Public License.
 */

#ifndef CONVERT_H
#define CONVERT_H

#include "cgeneral.h"
#include "cdate.h"

#ifdef __cplusplus
extern "C"
{
#endif




/*f
***************************************************************************
** Formats a TDate for printing. The print format is YYYYMMDD. 
** Can be called eight times from the same print statement, but not more.
***************************************************************************
*/
char* JpmcdsFormatDate(TDate date);


/*f
***************************************************************************
** Formats a TDateInterval.
** Can be called twice from the same print statement, but not more.
***************************************************************************
*/
char* JpmcdsFormatDateInterval
    (TDateInterval *interval);  /* (I) */


/*f
***************************************************************************
** Converts a date in string format to a TDate.
** The following formats are accepted:
**       "YYYYMMDD", "MM/DD/YYYY", "Mth-DD-YYYY", and "DD-Mth-YYYY".
***************************************************************************
*/
int JpmcdsStringToDate
    (char  *mmddyyyy,           /* (I) Month/Day/Year format */
     TDate *tdate);             /* (O) TDate */


/*f
***************************************************************************
** Converts a string representation ("2Y", "1S", etc) to a TDateInterval.
** Number optional, default is 1, hence, e.g., "M" (= "1M") is valid.
***************************************************************************
*/
int JpmcdsStringToDateInterval
    (char          *input,      /* (I) String w/ 1A, 3M, 4D, etc */
     char          *label,      /* (I) Label-for JpmcdsErr Msg only */
     TDateInterval *interval);  /* (O) Value read from file */


/*f
***************************************************************************
** Makes a TDateInterval
***************************************************************************
*/
int JpmcdsMakeDateInterval
    (int            numPeriods, /* (I) Number of periods (Can be negative) */
     char           periodType, /* (I) Period type (D,W,M,S,Q,A,I,J,K,E) */
     TDateInterval *interval);  /* (O) */


#ifdef __cplusplus
}
#endif

#endif    /* CONVERT_H */





