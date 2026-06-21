/*
 * ISDA CDS Standard Model
 *
 * Copyright (C) 2009 International Swaps and Derivatives Association, Inc.
 * Developed and supported in collaboration with Markit
 * 
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the ISDA CDS Standard Model Public License.
 */

#include "datelist.h"
#include <string.h>
#include <memory.h>
#include "cheaders.h"
#include "bastypes.h"
#include "yearfrac.h"
#include "date_sup.h"
#include "convert.h"
#include "ldate.h"
#include "macros.h"


/*
***************************************************************************
** Makes an array of dates from startDate, MaturityDate, & interval.
** If (maturityDate-startDate) divided by interval is not an integer,
** there is a stub. If stubAtEnd is set, the stub is placed at the end;
** otherwise, it is placed at the beginning.
** Unlike JpmcdsNewDateListExtended, the startDate and maturityDate are always
** included, and are the first and last dates respectively.
** Assuming there is no stub,  dates created are of the form:
** baseDate + idx*Interval, where startIdx <= idx <= Time2Maturity/Interval
** Returns a new DateList on success, and NULL on failure.
***************************************************************************
*/
TDateList* JpmcdsNewDateList(
    TDate startDate,                    /* (I) Start Date */
    TDate maturityDate,                 /* (I) Maturity Date */
    TDateInterval *interval,            /* (I) Increment */
    TBoolean stubAtEnd)                 /* (I) T=Stub at end; F=Stub at beg. */
{
    TDateList *dateList = NULL;
    TDateInterval intval = *interval;   /* Make local copy */
    int numIntervals;
    int extraDays;
    int numDates;
    static char routine[]="JpmcdsNewDateList";

    if (stubAtEnd)                      /* Count forwards from startDate */
    {
        if (JpmcdsCountDates(startDate, maturityDate, &intval,
            &numIntervals, &extraDays) == FAILURE)
            goto done;
    }
    else                                /* Count backwards from matDate */
    {
        intval.prd = -intval.prd;
        if (JpmcdsCountDates(maturityDate, startDate, &intval,
            &numIntervals, &extraDays) == FAILURE)
            goto done;
    }

    if (extraDays > 0)
        numDates = numIntervals + 2;
    else
        numDates = numIntervals + 1;

    dateList = JpmcdsNewEmptyDateList(numDates);
    if (dateList == NULL)
        goto done;
    
    /* Fill in the dates
     */
    if (stubAtEnd)                      /* Stub at end */
    {
        if (JpmcdsMakeTDateArray(startDate, &intval, 0, 1,
                              numDates-1, dateList->fArray) == FAILURE)
            goto done;
        dateList->fArray[numDates-1] = maturityDate;
    }
    else                                /* Stub at beginning */
    {
        if(JpmcdsMakeTDateArray(maturityDate, &intval, 0, -1, numDates-1,
                             &dateList->fArray[numDates-1]) == FAILURE)
            goto done;
        dateList->fArray[0] = startDate;        
    }

    /* Succeeded. */
    return dateList;  

    
done:
    JpmcdsErrMsg("%s: Failed.\n", routine);
    JpmcdsFreeDateList(dateList);
    return NULL;
}


/*
***************************************************************************
** Makes an array of dates using startDate, MaturityDate, & interval.
** If time2Maturity/interval is not an integer, there is a stub. If stubAtEnd
** is set, the stub is placed at the end; otherwise, it is placed at the
** beginning.
** If there is a stub at the END, the MaturityDate is NOT included
** in the list; but the NEXT date is.
** If there is a stub at the BEGINNING, the startDate is NOT included
** in the list; but the PREVIOUS date is.
** This is necessary for Actual/Actual stub calculation.
** Assuming there is no stub,  dates created are of the form:
** baseDate + idx*Interval, where startIdx <= idx <= Time2Maturity/Interval
** Returns a new DateList on success, and NULL on failure.
***************************************************************************
*/
TDateList* JpmcdsNewDateListExtended(
    TDate startDate,                    /* (I) Start Date */
    TDate maturityDate,                 /* (I) Maturity Date */
    TDateInterval *interval,            /* (I) Increment */
    TBoolean stubAtEnd)                 /* (I) T=Stub at end; F=Stub at beg. */
{
    static char routine[]="JpmcdsNewDateListExtended";
    TDateList *dateList = NULL;
    TDateInterval intval = *interval;   /* Make local copy */
    int numIntervals;
    int extraDays;
    int numDates;


    if (stubAtEnd)                      /* Count forwards from startDate */
    {
        if (JpmcdsCountDates(startDate, maturityDate, &intval,
            &numIntervals, &extraDays) == FAILURE)
            goto done;
    }
    else                                /* Count backwards from matDate */
    {
        intval.prd = -intval.prd;
        if (JpmcdsCountDates(maturityDate, startDate, &intval,
            &numIntervals, &extraDays) == FAILURE)
            goto done;
    }

    if (extraDays > 0)
        numDates = numIntervals + 2;
    else
        numDates = numIntervals + 1;

    dateList = JpmcdsNewEmptyDateList(numDates);
    if (dateList == NULL)
        goto done;
    
    /* Fill in the dates
     */
    if (stubAtEnd)                      /* Stub at end */
    {
        if (JpmcdsMakeTDateArray(startDate, &intval, 0, 1,
                              numDates, dateList->fArray) == FAILURE)
            goto done;
    }
    else                                /* Stub at beginning */
    {
        if(JpmcdsMakeTDateArray(maturityDate, &intval, 0, -1, numDates,
                             &dateList->fArray[numDates-1]) == FAILURE)
            goto done;
    }

    /* Success */
    return dateList;

done:
    JpmcdsErrMsg("%s: Failed.\n", routine);
    JpmcdsFreeDateList(dateList);
    return NULL;
}


/*
***************************************************************************
** Makes an array of dates using startDate, MaturityDate, rollDate & 
** interval. If rollDate is non-zero and stubAtEnd is False, then the rollDate
** is used to generate cash flow dates. First starting from roll-date we see if
** we can count forward an integral number of intervals to the maturity date.
** If this is so then we can generate our datelist by counting forward from
** the rolldate.
**
** Returns a new DateList on success, and NULL on failure.
***************************************************************************
*/
TDateList* JpmcdsNewDateListExtendedRoll(
    TDate startDate,                    /* (I) Start Date */
    TDate maturityDate,                 /* (I) Maturity Date */
    TDate rollDate,                     /* (I) Date for roll information */
    TDateInterval *interval,            /* (I) Increment */
    TBoolean stubAtEnd)                 /* (I) T=Stub at end; F=Stub at beg. */
{
    static char routine[]="JpmcdsNewDateListExtendedRoll";
    int   numIntervals;
    int   extraDays;
    TDate firstDate;


    if (rollDate == 0 || stubAtEnd == TRUE)
    {
        return JpmcdsNewDateListExtended(startDate, maturityDate, interval, stubAtEnd);
    }

    if (JpmcdsCountDates(startDate, maturityDate, interval, &numIntervals, &extraDays) != SUCCESS)
        goto done;
    
    if (extraDays == 0)
    {
        /* We can use the roll date - so we count forward from this day */
        if (startDate < rollDate)
        {
            JpmcdsErrMsg("%s: Start date cannot be before roll date.\n", routine);
            goto done;
        }

        /* We roll forward until we are on the startDate or if there are
         * a non-integral number of periods so that we are on the flow
         * date preceeding the start date
         */
        if (JpmcdsCountDates(rollDate, startDate, interval, &numIntervals, &extraDays) != SUCCESS)
            goto done;

        if (extraDays == 0)
        {
            firstDate = startDate;
        }
        else
        {
            if (JpmcdsDateFromDateAndOffset(rollDate, interval, numIntervals, &firstDate) != SUCCESS)
                goto done;
        }

        return JpmcdsNewDateListExtended(firstDate, maturityDate, interval, TRUE);
    }
    else
    {
        return JpmcdsNewDateListExtended(startDate, maturityDate, interval, stubAtEnd);
    }

done:
    JpmcdsErrMsg("%s: Failed.\n", routine);
    return NULL;
}


/*
***************************************************************************
** Allocates a new TDateList and copies supplied dates into it from 
** supplied TCurve.
**
** This provides a convenient way to isolate the dates of a TCurve.
***************************************************************************
*/
TDateList* JpmcdsNewDateListFromTCurve
    (TCurve     *curve)         /* (I) Given curve */
{
    static char routine[]="JpmcdsNewDateListFromTCurve";
    int idx;
    TDateList *newDateList = JpmcdsNewEmptyDateList(curve->fNumItems);
    if (newDateList == (TDateList *)NULL)
        goto done;
    
    for (idx=0; idx < curve->fNumItems; idx++)
    {
        newDateList->fArray[idx] = curve->fArray[idx].fDate;
    }

done:
    if (newDateList == (TDateList *)NULL)
        JpmcdsErrMsg("%s: Failed.\n", routine);
    
    return newDateList;
}


/*
***************************************************************************
** Allocates a new TDateList by calling JpmcdsNewDateList, and then removing 
** the startDate.
***************************************************************************
*/
TDateList* JpmcdsNewPayDates
    (TDate          startDate,          /* (I) This date is not included */
     TDate          matDate,            /* (I) */
     TDateInterval *payInterval,        /* (I) */
     TBoolean       stubAtEnd)          /* (I) */
{
    static char routine[]="JpmcdsNewPayDates";
    int idx;
    TDateList *payDates;
    
    payDates = JpmcdsNewDateList(startDate, matDate, payInterval,
                              stubAtEnd);
    
    if (payDates == (TDateList *)NULL)
        goto done;
    
    /* Now remove startDate, and move all dates back by one.
     */
    for (idx=0; idx < payDates->fNumItems-1; idx++)
    {
        payDates->fArray[idx] = payDates->fArray[idx+1];
    }
    payDates->fNumItems--;

done:
    if (payDates == (TDateList *)NULL)
        JpmcdsErrMsg("%s: Failed.\n", routine);
    
    return payDates;
}


/*
***************************************************************************
** Creates a new business-day adjusted dateList from an UN-adjusted
** datelist, using the supplied busisness day adjustment convention.
***************************************************************************
*/
TDateList* JpmcdsNewDateListBusDayAdj
(
    TDateList *dateList,        /* (I) Unadjusted date list  */ 
    long       badDayConv,      /* (I) Business day adjustment convention */
    char      *holidayFile      /* (I) Holiday date file */
)
{
    static char routine[]="JpmcdsNewDateListBusDayAdj";
    TDateList *adjDateList;
    int status;

    adjDateList = JpmcdsNewEmptyDateList( dateList->fNumItems);
    if(adjDateList == (TDateList *)NULL) 
        goto error;
    

    (void)memcpy((char *)adjDateList->fArray, (char *)dateList->fArray,
                 (size_t)(dateList->fNumItems * sizeof(TDate)));

    status = JpmcdsDateListBusDayAdj(adjDateList, badDayConv, holidayFile);

    if (status == SUCCESS)
        return adjDateList;

error:
   JpmcdsErrMsg(" %s: Failed.\n", routine);  
   JpmcdsFreeDateList(adjDateList);
   return (TDateList *) NULL;
}


/*
***************************************************************************
** Adjusts a datelist (in-place) according to the supplied business day 
** adjustment convention.
***************************************************************************
*/
int JpmcdsDateListBusDayAdj
(
    TDateList *dateList,        /* (I/O) Unadjusted date list  */ 
    long       badDayConv,      /* (I)   Business day adjustment convention */
    char      *holidayFile      /* (I)   Holiday date file */
)
{
    static char routine[]="JpmcdsDateListBusDayAdj";
    long idx;

    for( idx=0; idx < dateList->fNumItems; idx++)
    {
        if( JpmcdsBusinessDay(dateList->fArray[idx], badDayConv, 
                           holidayFile, &dateList->fArray[idx]) == FAILURE)
        {
            JpmcdsErrMsg(" %s: Failed.\n", routine);  
            return FAILURE;
        }
    }

    return SUCCESS;
}


