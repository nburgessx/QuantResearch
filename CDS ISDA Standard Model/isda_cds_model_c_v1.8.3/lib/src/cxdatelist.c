/*
 * ISDA CDS Standard Model
 *
 * Copyright (C) 2009 International Swaps and Derivatives Association, Inc.
 * Developed and supported in collaboration with Markit
 * 
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the ISDA CDS Standard Model Public License.
 */

#include "cxdatelist.h"
#include <assert.h>
#include "cxbsearch.h"
#include "macros.h"
#include "cerror.h"
#include "ldate.h"


static int qsortDateCompare(const void *date1, const void *date2)
{
    TDate d1 = *(TDate*)date1;
    TDate d2 = *(TDate*)date2;

    if (d1 == d2) return 0;
    if (d1 < d2)  return -1;
    return 1;
}


/*
***************************************************************************
** Adds dates to a TDateList. 
**
** If the original date list and date list to be added are sorted, then
** the resulting date list will be sorted and duplicate dates will be
** removed. Sorting assumes ascending order ([0] < [1] etc).
**
** If either of the inputs are not sorted, then the resulting date list
** will not be sorted, and some duplicates may remain.
**
** For efficiency, we do not automatically try to sort the resulting
** date list for unsorted inputs. Sorting the date list each time appears
** to be a huge performance issue in some algorithms (where the input
** dates would all be sorted anyway).
**
** Note that if dl is NULL, then this will create a new date list from
** the given dates.
**
** Note that if numItems=0, this will copy the given date list.
***************************************************************************
*/
TDateList* JpmcdsDateListAddDates
(TDateList *dl,          /* (I) Initial date list            */
 int         numItems,    /* (I) Number of dates to add       */
 TDate      *array)      /* (I) [numItems] Dates to be added */
{
    static char routine[] = "JpmcdsDateListAddDates";
    int         status    = FAILURE;

    TDateList tmp = {0, NULL};
    TDateList *result = NULL;

    REQUIRE (numItems >= 0);
    REQUIRE (dl == NULL || dl->fNumItems >= 0);

    if (dl == NULL)
    {
        result = JpmcdsNewDateListFromDates (array, numItems);
    }
    else if (numItems <= 0)
    {
        result = JpmcdsCopyDateList (dl);
    }
    else if (dl->fNumItems == 0 && numItems == 0)
    {
        result = JpmcdsNewDateListFromDates (NULL, 0);
    }
    else
    {
        int totalItems = dl->fNumItems + numItems;
        int i = 0;
        int j = 0;
        int k = 0;

        result = JpmcdsNewEmptyDateList (totalItems);
        if (result == NULL)
            goto done;

        while (i < dl->fNumItems && j < numItems)
        {
            if (dl->fArray[i] == array[j])
            {
                /* exclude duplicates */
                ++j;
                --totalItems;
            }
            else if (dl->fArray[i] < array[j])
            {
                result->fArray[k] = dl->fArray[i];
                ++i;
                ++k;
            }
            else
            {
                assert (dl->fArray[i] > array[j]);
                result->fArray[k] = array[j];
                ++j;
                ++k;
            }
        }

        if (i < dl->fNumItems)
        {
            int n = dl->fNumItems - i;
            COPY_ARRAY (result->fArray+k, dl->fArray+i, TDate, n);
            k += n;
        }

        if (j < numItems)
        {
            int n = numItems - j;
            COPY_ARRAY (result->fArray+k, array+j, TDate, n);
            k += n;
        }

        assert (k == totalItems);
        result->fNumItems = totalItems;
    }
    if (result == NULL)
        goto done;

    status = SUCCESS;

 done:

    if (status != SUCCESS)
    {
        JpmcdsErrMsgFailure (routine);
        JpmcdsFreeDateList (result);
        result = NULL;
    }

    FREE(tmp.fArray);
    
    return result;
}


/*
***************************************************************************
** Adds dates to a TDateList and frees the input date list.
**
** If the original date list and date list to be added are sorted, then
** the resulting date list will be sorted and duplicate dates will be
** removed. Sorting assumes ascending order ([0] < [1] etc).
**
** If either of the inputs are not sorted, then the resulting date list
** will not be sorted, and some duplicates may remain.
**
** For efficiency, we do not automatically try to sort the resulting
** date list for unsorted inputs. Sorting the date list each time appears
** to be a huge performance issue in some algorithms (where the input
** dates would all be sorted anyway).
**
** Note that if dl is NULL, then this will create a new date list from
** the given dates.
**
** Note that if numItems=0, this will copy the given date list.
**
** The input date list is FREE'd by this routine. Thus if you have an
** algorithm which involves building up a datelist gradually, you can
** do something like this:
**
**  TDateList* dl = NULL;
**  ...
**  dl = JpmcdsDateListAddDatesFreeOld (dl, numItems, array);
**  if (dl == NULL) goto done;
**  ..
**  dl = JpmcdsDateListAddDatesFreeOld (dl, numItems, array);
**  if (dl == NULL) goto done;
**  ..
**  etc.
**
** with the point being that you don't have to worry about the original
** date list at each step since this routine frees it for you.
***************************************************************************
*/
TDateList* JpmcdsDateListAddDatesFreeOld
(TDateList *dl,         /* (I/O) Initial date list - gets freed */
 int         numItems,   /* (I) Number of dates to add           */
 TDate     *array)      /* (I) [numItems] Dates to be added     */
{
    static char routine[] = "JpmcdsDateListAddDatesFreeOld";
    
    TDateList *output;

    output = JpmcdsDateListAddDates (dl, numItems, array);
    JpmcdsFreeDateList (dl);
    
    if (output == NULL)
        JpmcdsErrMsgFailure(routine);

    return output;
}


/*
***************************************************************************
** Truncates a date list at the specified date. The resulting date list
** will contain all dates previous to (or following) the specified date.
** Dates in the datelist which match the specified date may be optionally
** included.
**
** The datelist may optionally be modified in place or a new copy is
** returned.
**
** The input date list must be sorted.
***************************************************************************
*/
TDateList* JpmcdsDateListTruncate
(TDateList *dateList,       /* (I/O) Date list to be modified in place */
 TDate      truncationDate, /* (I) Date on which to perform trunctation */
 TBoolean      inclusive,      /* (I) TRUE=include truncation date if in list*/
 TBoolean      excludeBefore,  /* (I) TRUE=exclude dates before truncation date*/
 TBoolean      inPlace         /* (I) TRUE=modify date list in place */
)
{
    static char  routine[] = "JpmcdsDateListTruncate";
    TDateList  *truncated = NULL;
    int          numItems;
    int          size;
    int          offset;
    int          truncatePt;
    int          i;

    REQUIRE (dateList != NULL);

    /*  First we find the truncation point in the date list and the size
     *  of the new date list
     */

    /* perhaps we could do this more efficiently with a binary search -
       but this is the code which was imported and therefore should work */
    numItems = dateList->fNumItems;
    if (excludeBefore)
    {
        truncatePt = 0;
        for (i = 0; i < numItems; i++)
        {
            if (dateList->fArray[i] > truncationDate)
            {
                truncatePt = i;
                break;
            }
            if (inclusive && dateList->fArray[i] == truncationDate)
            {
                truncatePt = i;
                break;
            }
        }
        size = numItems - truncatePt;
        offset = truncatePt;
    }
    else
    {
        truncatePt = numItems - 1;
        for (i = numItems - 1; i > 0; i--)
        {
            if (dateList->fArray[i] < truncationDate)
            {
                truncatePt = i;
                break;
            }
            if (inclusive && dateList->fArray[i] == truncationDate)
            {
                truncatePt = i;
                break;
            }
        }
        size = truncatePt + 1;
        offset = 0;
    }

    /* Next we get a pointer to the datelist where we store the result */
    if (inPlace)
    {
        truncated = dateList;
    }
    else
    {
        truncated = JpmcdsNewEmptyDateList (size);
        if (truncated == NULL)
            goto done;
    }

    /* finally we populate the result */
    if (!inPlace || offset != 0)
    {
        for (i = 0; i < size; i++)
        {
            truncated->fArray[i] = dateList->fArray[i + offset];
        }
    }
    truncated->fNumItems = size;

 done:

    if (truncated == NULL)
        JpmcdsErrMsgFailure(routine);

    return truncated;
}


/*
***************************************************************************
** Makes a date list from a given start date to a given end date with dates
** seperated by a given interval.
**
** Use the stub parameter to determine whether the stub appears at the
** start or the end of the date list, and whether the stub is long or
** short.
**
** The start date and end date are always both in the date list.
** The end date must be strictly after the start date.
** The date interval must be positive.
***************************************************************************
*/
TDateList* JpmcdsDateListMakeRegular
(TDate          startDate,  /* (I) Start date */
 TDate          endDate,    /* (I) End date */
 TDateInterval *interval,   /* (I) Date interval */
 TStubMethod   *stubType)   /* (I) Stub type */
{
    static char routine[] = "JpmcdsDateListMakeRegular";
    int         status = FAILURE;

    TDateList *dl = NULL;
    TDate      tmpDates[100];
    int        i;
    int        numIntervals;
    int        numTmpDates = sizeof(tmpDates) / sizeof(TDate);
    int        totalDates = 0;
    TDate      date;
    TDateInterval multiInterval;

    REQUIRE (interval != NULL);
    REQUIRE (interval->prd > 0);
    REQUIRE (endDate > startDate);

    /* we calculate tmpDates in blocks of 100 and add to the datelist */
    if (!stubType->stubAtEnd)
    {
        /* front stub - so we start at endDate and work backwards */
        numIntervals = 0;
        i            = numTmpDates;
        date         = endDate;
        while (date > startDate)
        {
            if (i == 0)
            {
                dl = JpmcdsDateListAddDatesFreeOld (dl, numTmpDates, tmpDates);
                if (dl == NULL)
                    goto done;
                i = numTmpDates;
            }
            --i;
            --numIntervals;
            ++totalDates;
            assert (i >= 0);
            tmpDates[i] = date;

            SET_TDATE_INTERVAL(multiInterval, interval->prd * numIntervals, interval->prd_typ);
            if (JpmcdsDtFwdAny(endDate, &multiInterval, &date) != SUCCESS)
                goto done;
        }
        assert (totalDates > 0);
        assert (date <= startDate);
        if (date == startDate || totalDates == 1 || !stubType->longStub)
        {
            /* don't change existing tmpDates[] but need to add startDate */
            if (i == 0)
            {
                dl = JpmcdsDateListAddDatesFreeOld (dl, numTmpDates, tmpDates);
                if (dl == NULL)
                    goto done;
                i = numTmpDates;
            }
            --i;
            ++totalDates;
            tmpDates[i] = startDate;
        }
        else
        {
            assert (!stubType->stubAtEnd && stubType->longStub);
            assert (date < startDate);
            /* the existing date in tmpDates[] should be changed to be
               the start date */
            tmpDates[i] = startDate;
        }
        /* now add from tmpDates[i] to tmpDates[numTmpDates-1] to date list */
        dl = JpmcdsDateListAddDatesFreeOld (dl, numTmpDates-i, tmpDates+i);
        if (dl == NULL)
            goto done;
    }
    else
    {
        /* back stub - so we start at startDate and work forwards */
        numIntervals = 0;
        i            = -1;
        date         = startDate;
        while (date < endDate)
        {
            ++i;
            ++totalDates;
            if (i == numTmpDates)
            {
                dl = JpmcdsDateListAddDatesFreeOld (dl, numTmpDates, tmpDates);
                if (dl == NULL)
                    goto done;
                i = 0;
            }
            ++numIntervals;
            assert (i < numTmpDates);
            tmpDates[i] = date;

            SET_TDATE_INTERVAL(multiInterval, interval->prd * numIntervals, interval->prd_typ);
            if (JpmcdsDtFwdAny(startDate, &multiInterval, &date) != SUCCESS)
                goto done;
        }
        assert (totalDates > 0);
        assert (date >= endDate);
        if (date == endDate || totalDates == 1 || stubType->stubAtEnd && !stubType->longStub)
        {
            /* don't change existing tmpDates[] but need to add endDate */
            ++i;
            ++totalDates;
            if (i == numTmpDates)
            {
                dl = JpmcdsDateListAddDatesFreeOld (dl, numTmpDates, tmpDates);
                if (dl == NULL)
                    goto done;
                i = 0;
            }
            tmpDates[i] = endDate;
        }
        else
        {
            assert (stubType->stubAtEnd && stubType->longStub);
            assert (date > endDate);
            /* the existing date in tmpDates[] should be changed to be
               the end date */
            tmpDates[i] = endDate;
        }
        /* now add from tmpDates[0] to tmpDates[i] to the date list */
        dl = JpmcdsDateListAddDatesFreeOld (dl, i+1, tmpDates);
        if (dl == NULL)
            goto done;
    }
    ASSERT (totalDates >= 2);
    ASSERT (dl->fNumItems == totalDates);
    status = SUCCESS;

 done:

    if (status != SUCCESS)
    {
        JpmcdsFreeDateList (dl);
        dl = NULL;
        JpmcdsErrMsgFailure (routine);
    }

    return dl;
}

