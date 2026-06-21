/*
 * ISDA CDS Standard Model
 *
 * Copyright (C) 2009 International Swaps and Derivatives Association, Inc.
 * Developed and supported in collaboration with Markit
 * 
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the ISDA CDS Standard Model Public License.
 */

#include "date_sup.h"
#include <ctype.h>
#include <string.h>
#include <math.h>
#include "cerror.h"
#include "cgeneral.h"
#include "ldate.h"
#include "macros.h"
#include "mdydate.h"
#include "convert.h"
#include "dateconv.h"
#include "macros.h"
#include "cdate.h"


/*
***************************************************************************
** Checks validity of TDateInterval.
**
** The following TDateIntervals are considered to be invalid:
** - Zero intervals ( = 0 days)
** - intervals, whose signum is opposite to that of toDate - fromDate
***************************************************************************
*/
static TBoolean checkDateInterval(
                        TDateInterval *ptrDateInterval,
                        TDate fromDate,
                        TDate toDate,
                        char *caller);


/*m
** Macro: determines if a year has a leap day in it.
** Pass in the year number. Such as 1998. used by JpmcdsIntervalBetweenDates
*/
#define EXACTLY_DIVIDES(a,b)    (((a) % (b)) == 0)
#define YEAR_HAS_A_LEAP_DAY(y) (EXACTLY_DIVIDES(y,400) || \
 (EXACTLY_DIVIDES(y,4) && !EXACTLY_DIVIDES(y,100)))


/*
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
    TDate *dateArray)                   /* (O) array of dates */
{
    int idx;
    TDateInterval offsetInterval;
    offsetInterval.prd_typ = interval->prd_typ;

    for (idx = 0; idx < ndates; idx++)
    {
        offsetInterval.prd = interval->prd * (startIdx + idx);

        if (JpmcdsDtFwdAny(baseDate, &offsetInterval,
            &dateArray[idx*arrayIncrement]) == FAILURE)
            return FAILURE;
    }

    return SUCCESS;
}


/*
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
    int *extraDays)                     /* (O) Days left over(remainder) */
{
    static char routine[]="JpmcdsCountDates";
    int status = FAILURE;     /* Until proven successful */
    double intervalYears;     /* TDateInterval expressed in years */
    double fromToYears;       /* interval between to & from dates in years */
    int lowNumIntervals;      /* low estimation of *numIntervals */
    int index;                /* running approximation of *numIntervals */
    TDate currDate;           /* date index times interval away from fromDate*/
    TDate lastDate;           /* last date fitting guaranteed to be between */

   /* Make sure interval has the right sign.
    */
    if (checkDateInterval(interval,fromDate,toDate,routine) == FAILURE)
        goto done;

    /* Convert TDateInterval to years
     */
    if (JpmcdsDateIntervalToYears(interval, &intervalYears) == FAILURE)
        goto done;

    /* Convert toDate-fromDate to years
     */
    fromToYears = (double)(toDate-fromDate)/JPMCDS_DAYS_PER_YEAR;

   /* Compute & round the quotient to the lower integer and then subtract TWO
    * to make sure that corresponding lowToDate is INSIDE
    */
    lowNumIntervals = MAX(0,(int)floor(fabs(fromToYears/intervalYears))-2);

   /* Initialize the LOOP by index of intervals
    */
   index = lowNumIntervals;
   if (JpmcdsDateFromDateAndOffset(fromDate,interval,index,&currDate) == FAILURE)
       goto done;
   lastDate = currDate;

   /* Keep advancing currDate while currDate is between fromDate and toDate.
    */
   while ( IS_BETWEEN (currDate,fromDate,toDate) )
   {
       ++index;
       lastDate = currDate;
       if (JpmcdsDateFromDateAndOffset(fromDate,interval,index,&currDate)
           == FAILURE)
           goto done;
   }

   *numIntervals=index-1; /* step back inside (fromDate,toDate) interval */
   if (*numIntervals < lowNumIntervals)
   {
      JpmcdsErrMsg("%s: Failed with parameters:\n"
                "\tfromDate: %s  toDate: %s interval: %s\n",
                routine,
                JpmcdsFormatDate(fromDate), JpmcdsFormatDate(toDate),
                JpmcdsFormatDateInterval(interval));
      goto done;
   }
    *extraDays = (int) ABS(toDate - lastDate);

    status = SUCCESS;

done:
    if (status == FAILURE)
        JpmcdsErrMsg("%s: Failed.\n", routine);

    return status;
}



#define DAYS_PER_MONTH 30


/*
***************************************************************************
** Converts a frequency (# times/year) to a TDateInterval
** Returns SUCCESS/FAILURE.
***************************************************************************
*/
int JpmcdsFreq2TDateInterval(
        long freq,                      /* (I) */
        TDateInterval *interval)        /* (O) */
{
    static char routine[]="JpmcdsFreq2TDateInterval";

    if (freq > 0 && freq <= 12)
    {
        SET_TDATE_INTERVAL((*interval), JPMCDS_MONTHS_PER_YEAR / (int)freq, 'M');
        return SUCCESS ;
    }
    else
    {
        JpmcdsErrMsg("%s: bogus frequency %ld detected.\n", routine, freq);
        return FAILURE ;
    }
}


/*
***************************************************************************
** Converts a TDateInterval to a frequency
***************************************************************************
*/
int JpmcdsDateIntervalToFreq(TDateInterval *interval, /* (I) */
                          double *freq) /* (O) # times per year */
{
    double years;
    static char routine[]="JpmcdsDateIntervalToFreq";

    if (JpmcdsDateIntervalToYears(interval, &years) == FAILURE)
    {
        JpmcdsErrMsg("%s: Failed.\n", routine);
        return FAILURE ;
    }
    if (years > 0.)
    {
        *freq = 1./years;
        return SUCCESS ;
    }
    else
    {
        JpmcdsErrMsg("%s: interval is zero.\n", routine);
        return FAILURE ;
    }
}


/*
***************************************************************************
** Converts a TDateInterval to # years. Note that if the TDateInterval is a
** month type (A,S,Q), the routine uses 30/360 to compute the year fraction. 
** If it is a day type (D,W), it uses Act/365F.
***************************************************************************
*/
int JpmcdsDateIntervalToYears
    (TDateInterval *interval,           /* (I) */
     double *years)                     /* (O) # times per year */
{
    static char routine[]="JpmcdsDateIntervalToYears";

    switch(toupper(interval->prd_typ))
    {
        case 'A':
        case 'Y':
            *years = (double)interval->prd; /* 30/360 */
            break;
        case 'S':
            *years = (double)(interval->prd)/2.; /* 30/360 */
            break;
        case 'Q':
        case 'I':
        case 'K':
        case 'L':
            *years = (double)(interval->prd)/4.; /* 30/360 */
            break;
        case 'M':
        case 'E':
        case 'F':
        case 'G':
        case 'H':
        case 'J':
        case 'T':
            *years = (double)(interval->prd)/JPMCDS_MONTHS_PER_YEAR; /* 30/360 */
            break;
        case 'W':
            *years = (double)interval->prd *
                     (double)JPMCDS_DAYS_PER_WEEK/JPMCDS_DAYS_PER_YEAR; /* Act/365F*/
            break;
        case 'D':
            *years = (double)interval->prd/JPMCDS_DAYS_PER_YEAR; /* Act/365F */
            break;
        case 'U':
            *years = (double)interval->prd *
                     (double)JPMCDS_DAYS_PER_LUNAR_MONTH/JPMCDS_DAYS_PER_YEAR; /* Act/365F*/
        default:
            JpmcdsErrMsg("%s: unknown interval type %c.\n", routine, interval->prd_typ);
            return FAILURE;
    }

    return SUCCESS;
}


/*
***************************************************************************
** Calculates a date from another date and offset represented by interval 
** and index: newDate = oldDate + index*interval
***************************************************************************
*/
int JpmcdsDateFromDateAndOffset (TDate oldDate,
                               TDateInterval *pInterval,
                               int index,
                               TDate *pNewDate)
{
    TDateInterval compoundInterval = *pInterval;

    compoundInterval.prd *= index;

    return JpmcdsDtFwdAny (oldDate,&compoundInterval,pNewDate);
}


/*
***************************************************************************
** Checks validity of TDateInterval.
**
** The following TDateIntervals are considered to be invalid:
** - Zero intervals ( = 0 days)
** - intervals, whose signum is opposite to that of toDate - fromDate
***************************************************************************
*/
static TBoolean checkDateInterval (
          TDateInterval *interval,
          TDate fromDate,
          TDate toDate,
          char *caller)
{
   if (interval->prd == 0)
   {
      JpmcdsErrMsg("%s: Zero TDateInterval is not accepted.\n",caller);
      return FAILURE;
   }

   if ((toDate - fromDate) * interval->prd  < 0)
   {
      JpmcdsErrMsg("%s: cannot count from %s to %s with interval %s.\n",
         caller, JpmcdsFormatDate(fromDate), JpmcdsFormatDate(toDate),
         JpmcdsFormatDateInterval(interval));
      return FAILURE;
   }
   return SUCCESS;
}
