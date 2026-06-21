/*
 * ISDA CDS Standard Model
 *
 * Copyright (C) 2009 International Swaps and Derivatives Association, Inc.
 * Developed and supported in collaboration with Markit
 * 
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the ISDA CDS Standard Model Public License.
 */

#ifndef DATE_SUP_H
#define DATE_SUP_H

#include "cgeneral.h"
#include "cdate.h"

#ifdef __cplusplus
extern "C"
{
#endif


/*f
***************************************************************************
** Returns an array of TDates.
** Returns SUCCESS/FAILURE.
***************************************************************************
*/
int JpmcdsMakeTDateArray(
    TDate baseDate,                     /* (I) Base date */
    TDateInterval *interval,            /* (I) Date increment */
    int startIdx,                       /* (I) # intervals to start at */
                                        /* 0=start @ basedate */
                                        /* 1=start @ baseDate + interval, etc*/
    int arrayIncrement,                 /* (I) For array, not dates */
                                        /* Usually +1 or -1 */
    int ndates,                         /* (I) number of dates to make */
    TDate *dateArray);                  /* (O) array of dates */


/*f
***************************************************************************
** Counts # TDateIntervals in a range of dates. Note that fromDate
** can be either later than or earlier than toDate.
**
** In order to make this routine speedier we start from a low estimate
** of result value that is calculated using division (speedy approach).
** Based on the low estimate we calculate the results.
**
** Note that we allways calculate currDate based on fromDate rather than
** lastDate.  We MUST do it because 
** date + (n*Interval) = ((date + Interval)+Interval) ... )+Interval
** is not always true
***************************************************************************
*/
int JpmcdsCountDates(
    TDate fromDate,                     /* (I) Date to count from */
    TDate toDate,                       /* (I) Date to count to */
    TDateInterval *interval,            /* (I) Interval to count */
    int *numIntervals,                  /* (O) Answer (Quotient) */
    int *extraDays);                    /* (O) Days left over(remainder) */


/*f
***************************************************************************
** Converts a frequency (# times/year) to a TDateInterval
** Returns SUCCESS/FAILURE.
***************************************************************************
*/
int JpmcdsFreq2TDateInterval
    (
     long freq,                      /* (I) # times per year */
     TDateInterval *interval);       /* (O) */


/*f
***************************************************************************
** Converts a TDateInterval to a frequency
***************************************************************************
*/
int JpmcdsDateIntervalToFreq
    (
     TDateInterval *interval,        /* (I) */
     double *freq);                  /* (O) # times per year */


/*f
***************************************************************************
** Converts a TDateInterval to # years. Note that if the TDateInterval is a
** month type (A,S,Q), the routine uses 30/360 to compute the year fraction. 
** If it is a day type (D,W), it uses Act/365F.
***************************************************************************
*/
int JpmcdsDateIntervalToYears
    (
     TDateInterval *interval,        /* (I) */
     double *years);                 /* (O) # Years */


/*f
***************************************************************************
** Calculates a date from another date and offset represented by interval 
** and index: newDate = oldDate + index*interval
***************************************************************************
*/
int JpmcdsDateFromDateAndOffset
    (
     TDate oldDate,                     /* (I) */
     TDateInterval *interval,           /* (I) */
     int index,                         /* (I) */
     TDate *newDate);                   /* (O) */


#ifdef __cplusplus
}
#endif

#endif    /* DATE_SUP_H */
