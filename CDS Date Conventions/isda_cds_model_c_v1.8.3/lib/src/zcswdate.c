/*
 * ISDA CDS Standard Model
 *
 * Copyright (C) 2009 International Swaps and Derivatives Association, Inc.
 * Developed and supported in collaboration with Markit
 * 
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the ISDA CDS Standard Model Public License.
 */

#include <math.h>
#include "cgeneral.h"
#include "macros.h"
#include "badday.h"
#include "cerror.h"
#include "convert.h"
#include "date_sup.h"
#include "datelist.h"
#include "ldate.h"
#include "busday.h"
#include "zcswdate.h"


/*
***************************************************************************
** First finds out if a date is on cycle, counting from the value date. 
** Then computes the previous coupon date.
***************************************************************************
*/
static int SetPrevDateAndOnCycle
     (TDate          valueDate,         /* (I) Value date */
      TDate          origDate,          /* (I) Unadjusted date  */
      int            freq,              /* (I) Time between canon dates */
      TBoolean      *onCycle,           /* (O) On cycle from valuedate */
      TDate         *prevDate);         /* (O) Previous coupon date */


/*
***************************************************************************
** Uses swap maturity dates adjusted for bad days and a TBadDayList to set 
** up a TSwapDates.
***************************************************************************
*/
TSwapDates *JpmcdsSwapDatesNewFromAdjusted
     (TDate          valueDate,           /* (I) Must be a good day. */
      int            freq,                /* (I) Used to set previous */
      TDate         *adjustedDates,       /* (I) Adjusted swap mat dates */
      int            numDates,            /* (I) Length of adjustedDates */
      TBadDayList   *badDayList)          /* (I) Maps good to bad & vice-vsa*/
{
    static char routine[]="JpmcdsSwapDatesNewFromAdjusted";
    int    status = FAILURE;
    int    idx;                 
 
    TSwapDates *sd = JpmcdsSwapDatesNewEmpty(numDates);
    if (sd == NULL)
        goto done;
    
    /* For each swap maturity date... */
    for (idx=0; idx < numDates; idx++)
    {
        sd->adjusted[idx] = adjustedDates[idx];
        sd->original[idx] = JpmcdsGood2BadBadDayList(badDayList, adjustedDates[idx]);

        /* Now set sd->onCycle[idx] and sd->previous[idx]. */
        if (SetPrevDateAndOnCycle
            (valueDate, sd->original[idx], freq,
             &sd->onCycle[idx], &sd->previous[idx]) == FAILURE)
            goto done;
    }

    status = SUCCESS;

  done:
    if (status == FAILURE)
    {
        JpmcdsErrMsg("%s: Failed.\n", routine);
        JpmcdsSwapDatesFree(sd);
        sd = NULL;
    }

    return sd;
}


/*
***************************************************************************
** Uses un-adjusted swap maturity dates, a holiday file, and a bad day 
** convention to set up a TSwapDates.
***************************************************************************
*/
TSwapDates *JpmcdsSwapDatesNewFromOriginal
     (TDate          valueDate,         /* (I) Must be a good day. */
      int            freq,              /* (I) Used to set previous */
      TDate         *originalDates,     /* (I) UN-Adjusted swap mat dates */
      int            numDates,          /* (I) Length of originalDates */
      TBadDayList   *badDayList,        /* (I) Maps good to bad & vice-vsa*/
      long           badDayConv,        /* (I) See JpmcdsIsBusinessDay */
      char          *holidayFile)       /* (I) List of holidays */
{
    static char routine[]="JpmcdsSwapDatesNewFromOriginal";
    int    status = FAILURE;
    int    idx;                 
 
    TSwapDates *sd = JpmcdsSwapDatesNewEmpty(numDates);
    if (sd == NULL)
        goto done;
    
    /* For each swap maturity date... */
    for (idx=0; idx < numDates; idx++)
    {
        /* Set original (UNadjusted) and adjusted dates.
         */
        sd->original[idx] = originalDates[idx];
        if (JpmcdsZCAdjustDate(originalDates[idx], 
                            badDayList, badDayConv, holidayFile, 
                            &sd->adjusted[idx]) == FAILURE)
            goto done;

        /* Now set sd->onCycle[idx] and sd->previous[idx]. */
        if (SetPrevDateAndOnCycle
            (valueDate, originalDates[idx], freq,
             &sd->onCycle[idx], &sd->previous[idx]) == FAILURE)
            goto done;
    }

    status = SUCCESS;

  done:
    if (status == FAILURE)
    {
        JpmcdsErrMsg("%s: Failed.\n", routine);
        JpmcdsSwapDatesFree(sd);
        sd = NULL;
    }

    return sd;
}


/*
***************************************************************************
** Constructor for TSwapDates.
***************************************************************************
*/
TSwapDates *JpmcdsSwapDatesNewEmpty(int numDates)
{
    static char routine[]="JpmcdsSwapDatesNewEmpty";
    int    status = FAILURE;
    TSwapDates *sd = NEW(TSwapDates);
    if (sd == NULL)
        goto done;

    sd->numDates = numDates;

    if (numDates > 0)
    {
        sd->adjusted = NEW_ARRAY(TDate, numDates);
        sd->original = NEW_ARRAY(TDate, numDates);
        sd->previous = NEW_ARRAY(TDate, numDates);
        sd->onCycle = NEW_ARRAY(TBoolean, numDates);

        if (sd->adjusted == NULL ||
            sd->original == NULL ||
            sd->previous == NULL ||
            sd->onCycle == NULL)
            goto done;
    }
    else
    {
        sd->adjusted = NULL;
        sd->original = NULL;
        sd->previous = NULL;
        sd->onCycle = NULL;
    }
    status = SUCCESS;

  done:
    if (status == FAILURE)
    {
        JpmcdsErrMsg("%s: Failed.\n", routine);
        JpmcdsSwapDatesFree(sd);
        sd = NULL;
    }
    return sd;
}


/*
***************************************************************************
** Destructor for TSwapDates.
***************************************************************************
*/
void JpmcdsSwapDatesFree(TSwapDates *sd)
{
    if (sd != NULL)
    {
        FREE(sd->adjusted);
        FREE(sd->original);
        FREE(sd->previous);
        FREE(sd->onCycle);
        FREE(sd);
    }
}


/*
***************************************************************************
** First finds out if a date is on cycle, counting from the value date. 
** Then computes the previous coupon date.
***************************************************************************
*/
static int SetPrevDateAndOnCycle
     (TDate          valueDate,         /* (I) Value date */
      TDate          origDate,          /* (I) UNadjusted date  */
      int            freq,              /* (I) Used to set previous */
      TBoolean      *onCycle,           /* (O) On cycle from valuedate */
      TDate         *prevDate)          /* (O) Previous coupon date */
{
    static char routine[]="SetPrevDateAndOnCycle";
    int    status = FAILURE;

    TDateInterval interval;             /* Translated from freq */
    int           numIntervals;         /* # coupon intervals to swap mat */
    int           extraDays;            /* number of extra days (stub) */
    TMonthDayYear valueMDY;
    TMonthDayYear origMDY;

    if (JpmcdsFreq2TDateInterval(freq, &interval) == FAILURE)
        goto done;

    if (JpmcdsDateToMDY(valueDate, &valueMDY) == FAILURE)
        goto done;

    if (JpmcdsDateToMDY(origDate, &origMDY) == FAILURE)
        goto done;

    if ((valueMDY.day <= 28) &&
        (origMDY.day <= 28))
    {
        /* 
         * We assume we can only be on cycle if date is not on or after 
         * the 29th of the month.
         */

        /* Find out if unadjusted date is on cycle.
         */
        if (JpmcdsCountDates(valueDate, origDate, &interval, &numIntervals, &extraDays)==FAILURE)
            goto done;

        *onCycle = (extraDays == 0);
    }
    else
    {
        *onCycle = FALSE;
    }

    /* Now compute prevDate. If on cycle, count forwards
     * from value date. If off cycle, count backwards from matDate.
     */
    if (*onCycle)                       /* On cycle */
    {
        if (JpmcdsDateFromDateAndOffset
            (valueDate, &interval, numIntervals-1,
             prevDate) == FAILURE)
            goto done;
    }
    else                                /* Off cycle */
    {
        if (JpmcdsDateFromDateAndOffset
            (origDate, &interval, -1,
             prevDate) == FAILURE)
            goto done;
    }
    status = SUCCESS;

  done:
    if (status == FAILURE)
        JpmcdsErrMsg("%s: Failed.\n", routine);

    return status;
}


/*
***************************************************************************
** Adjusts date to fall on a business day, if necesssary.
***************************************************************************
*/
int JpmcdsZCAdjustDate
   (TDate       date,           /* (I) Date to be adjusted */
    TBadDayList *badDayList,    /* (I) */
    long         badDayConv,    /* (I) */
    char        *holidayFile,   /* (I) */
    TDate       *adjDate)       /* (O) */
{
    static char *routine="JpmcdsZCAdjustDate";
    int          status = FAILURE;

    if (badDayList != NULL)
    {
        *adjDate = JpmcdsBad2GoodBadDayList(badDayList, date);
        return SUCCESS;
    }
    else
    {
        if (JpmcdsBusinessDay(date, badDayConv, holidayFile, adjDate) == FAILURE)
            goto done;
    }

    status = SUCCESS;

  done:
    if (status == FAILURE)
        JpmcdsErrMsg("%s: Failed.\n", routine);

    return status;
}
