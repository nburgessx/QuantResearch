/*
 * ISDA CDS Standard Model
 *
 * Copyright (C) 2009 International Swaps and Derivatives Association, Inc.
 * Developed and supported in collaboration with Markit
 * 
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the ISDA CDS Standard Model Public License.
 */

#include "dtlist.h"
#include <string.h>
#include "bsearch.h"
#include "cerror.h"
#include "convert.h"
#include "macros.h"


/*
***************************************************************************
** Merge the two input date lists into a single list.
**
** If allowDuplicates is TRUE all dates in both lists will be included, 
** if FALSE only unique dates will be in the results list
***************************************************************************
*/
static TDateList *mergeDateLists
(TDateList  *dateList1,        /* (I) */
 TDateList  *dateList2,        /* (I) */
 TBoolean    allowDuplicates); /* (I) TRUE if duplicates in merged list */


/*
***************************************************************************
** Creates a new empty TDateList.
**
** Note that it is permissible to create a TDateList with 0 dates.
***************************************************************************
*/
TDateList* JpmcdsNewEmptyDateList(int numDates)
{
    static char routine[]="JpmcdsNewEmptyDateList";
    TDateList *dateList = NULL;

    if (numDates < 0)
    {
        JpmcdsErrMsg("%s: # dates (%d) must be >= 0.\n", routine, numDates);
        goto done;
    }

    dateList = NEW(TDateList);
    if (dateList == NULL)
        goto done;

    dateList->fNumItems = numDates;

    if (numDates > 0)
    {
        dateList->fArray = NEW_ARRAY(TDate, numDates);
        if (dateList->fArray == NULL)
            goto done;
    }
    else                                    /* numItems=0 */
    {
        dateList->fArray = NULL;
    }

    /* Success */
    return dateList;

done:
    JpmcdsErrMsg("%s: Failed.\n",routine);
    JpmcdsFreeDateList(dateList);
    return NULL;
}


/*
***************************************************************************
** Extends a date list to the given number of dates. This involves
** re-allocating the array within the TDateList structure.
**
** Requesting a smaller number of dates returns FAILURE and leaves the
** dateList unchanged.
***************************************************************************
*/
int JpmcdsExtendDateList
(TDateList *dateList,    /* (I/O) Date list to extend. */
 int        numDates     /* (I) New size of the date list. */
)
{
    static char routine[] = "JpmcdsExtendDateList";
    int         status    = FAILURE;

    if (numDates < dateList->fNumItems)
    {
        JpmcdsErrMsg ("%s: Number of dates requested (%d) smaller than existing "
                   "size of date list (%d).\n",
                   routine, numDates, dateList->fNumItems);
        goto done;
    }

    if (numDates > dateList->fNumItems)
    {
        TDate *newArray = NEW_ARRAY(TDate, numDates);
        if (newArray == NULL)
            goto done;

        COPY_ARRAY (newArray, dateList->fArray, TDate, dateList->fNumItems);
        FREE_ARRAY (dateList->fArray);

        dateList->fNumItems = numDates;
        dateList->fArray    = newArray;
    }

    status = SUCCESS;

done:

    if (status != SUCCESS)
        JpmcdsErrMsgFailure (routine);

    return status;
}


/*f
***************************************************************************
** Frees a DateList
***************************************************************************
*/
void JpmcdsFreeDateList(TDateList *theDL)
{
    if (theDL != NULL)
    {
        FREE(theDL->fArray);
        FREE(theDL);
    }
}


/*f
***************************************************************************
** Makes a copy of a date list.
**
** Returns NULL on error, or if the original datelist was NULL.
***************************************************************************
*/
TDateList* JpmcdsCopyDateList
    (TDateList   *a)                    /* (I) Date list to copy */
{
    static char routine[]="JpmcdsCopyDateList";
    TDateList *b = NULL;

    if (a != NULL)
    {
        b = JpmcdsNewDateListFromDates(a->fArray,a->fNumItems);
        if (b==NULL)
            goto done;
    }
    else
    {
        JpmcdsErrMsg("%s: Cannot copy NULL TDateList.\n", routine);
    }

  done:
    if (b == NULL)
        JpmcdsErrMsg("%s: Failed.\n", routine);

    return b;
}


/*
***************************************************************************
** Allocates a new TDateList and copies supplied dates into it.
***************************************************************************
*/
TDateList* JpmcdsNewDateListFromDates
    (TDate      *dates,         /* (I) */
     int         numDates)      /* (I) */
{
    static char routine[]="JpmcdsNewDateListFromDates";
    int idx;
    TDateList *newDateList = JpmcdsNewEmptyDateList(numDates);
    if (newDateList == (TDateList *)NULL)
        goto done;

    for (idx=0; idx < numDates; idx++)
    {
        newDateList->fArray[idx] = dates[idx];
    }

done:
    if (newDateList == (TDateList *)NULL)
        JpmcdsErrMsg("%s: Failed.\n", routine);

    return newDateList;
}


/*f
***************************************************************************
** Merges two date lists such that the resultant date list represents
** the union of the original lists (no duplicates).  Resultant list
** will be in sorted order.
**
** Note:
** 1. Both input lists must be in sorted order
** 2. Each original date list must not contain duplicates
** 3. If either of the datelists is NULL, it just makes a copy of
**    the other and returns it. If both are NULL, it returns NULL.
***************************************************************************
*/
TDateList* JpmcdsMergeDateLists
    (TDateList *dateList1,      /* (I) */
     TDateList *dateList2)      /* (I) */
{

    TDateList *outputDates = NULL;
    static char routine[] = "JpmcdsMergeDateLists";

    outputDates = mergeDateLists (dateList1, dateList2, FALSE);
    if (outputDates == (TDateList *)NULL)
        JpmcdsErrMsg ("%s: Failed to Merge Date Lists.\n", routine);

    return outputDates;
}


/*
***************************************************************************
** Sorts a datelist.
** 
** Returns NULL on FAILURE, and the sorted DateList on success.
***************************************************************************
*/
static int compareTDatesFcn(const void *date1, const void *date2);

TDateList* JpmcdsSortDateList
     (TDateList *dl)           /* (I) */
{
    static char routine[]="JpmcdsSortDateList";
    TDateList *distinctDL = NULL;
    TDateList *sortedDL = NULL;

    sortedDL = JpmcdsNewDateListFromDates(dl->fArray,
                                       dl->fNumItems);
    if (sortedDL == NULL)
        goto done;

    qsort (sortedDL->fArray,
           (size_t)sortedDL->fNumItems,
           (size_t)sizeof(TDate),
           compareTDatesFcn);

    distinctDL = JpmcdsNewDateListDistinct(sortedDL);

  done:
    JpmcdsFreeDateList(sortedDL);
    if (distinctDL == NULL)
        JpmcdsErrMsg("%s: Failed.\n", routine);

    return distinctDL;
}


static int compareTDatesFcn(const void *date1, const void *date2)
{

    TDate  result;
    TDate *date1p = (TDate *)date1;
    TDate *date2p = (TDate *)date2;

    result = *date1p - *date2p;
    if ((long)result < 0)
    {
        return (-1);
    }
    else if ((long)result > 0)
    {
        return (1);
    }
    return (0);
}


/*
***************************************************************************
** Merge the two input date lists into a single list.
**
** If allowDuplicates is TRUE all dates in both lists will be included, 
** if FALSE only unique dates will be in the results list
***************************************************************************
*/
static TDateList *mergeDateLists
(TDateList  *dateList1,       /* (I) */
 TDateList  *dateList2,       /* (I) */
 TBoolean    allowDuplicates) /* (I) TRUE if duplicates in merged list */
{

    TDateList  *outputDates = NULL;

    if (dateList1 == NULL)
    {
        outputDates = JpmcdsCopyDateList(dateList2);
    }
    else if (dateList2 == NULL)
    {
        outputDates = JpmcdsCopyDateList(dateList1);
    }
    else
    {
        outputDates = JpmcdsJpmcdsDateListAddDates (dateList1,
                                           dateList2->fNumItems,
                                           dateList2->fArray,
                                           allowDuplicates);
        if (outputDates == (TDateList *)NULL)
        {
            goto done;
        }
    }

 done:

    /* don't print error message, will be done by calling routine */
    return (outputDates);
}                                     /* End mergeDateLists */


/*
***************************************************************************
** Adds dates to a date list, and returns a new date list. Both the list of
** dates and the date list must be sorted. Keeps duplicates dates depending
** on the allowDuplicate flag (FALSE = remove duplicates, TRUE = keep them).
***************************************************************************
*/
TDateList* JpmcdsJpmcdsDateListAddDates
(TDateList *dateList,        /* (I) Original date list. */
 int        numDates,        /* (I) Number of dates to add. */
 TDate     *addDates,        /* (I) Dates to be added. */
 TBoolean   allowDuplicates  /* (I) Keep duplicate dates (or not) */
)
{
    static char routine[] = "JpmcdsJpmcdsDateListAddDates";
    int         status  = FALSE;

    int         indx1   = 0;
    int         indx2   = 0;
    int         indxOut;
    TDateList  *outputDates = NULL;

    if (dateList == NULL ||
        (numDates > 0 && addDates == NULL))
    {
        JpmcdsErrMsg ("%s: NULL inputs.\n", routine);
        goto done;
    }

    if (dateList->fNumItems < 0 || numDates < 0)
    {
        JpmcdsErrMsg ("%s: Negative number of dates specified.\n", routine);
        goto done;
    }

    /* This may be bigger than it has to be - don't worry about waste. */
    outputDates = JpmcdsNewEmptyDateList(dateList->fNumItems + numDates);
    if (outputDates == NULL)
    {
        goto done;
    }

    /* Now merge the two lists - repeat until one list is exhausted */
    for (indxOut = 0;
         indx1 < dateList->fNumItems && indx2 < numDates;
         indxOut++)
    {

        /* Three cases to consider: List1 <, >, = List 2 */
        if (dateList->fArray[indx1] < addDates[indx2])
        {

            /* Take from List 1 */
            outputDates->fArray[indxOut] = dateList->fArray[indx1++];
        }
        else if (dateList->fArray[indx1] > addDates[indx2])
        {

            /* Take from List 2 */
            outputDates->fArray[indxOut] = addDates[indx2++];
        }
        else
        {

            /* The two dates are equal, grab 1st & bump both indices */
            outputDates->fArray[indxOut] = dateList->fArray[indx1++];
            if (allowDuplicates)
            {
                indxOut++;
                outputDates->fArray[indxOut] = addDates[indx2++];
            }
            else
            {
                indx2++;
            }
        }
    } /* endfor - exhaust one of the lists */

    /*
    ** Now tack remaining dates onto list.  No need to check for which
    ** list was exhausted - either or both 'for' loops will be null
    */
    for ( ; indx1 < dateList->fNumItems; indxOut++)
    {
        outputDates->fArray[indxOut] = dateList->fArray[indx1++];
    }
    for ( ; indx2 < numDates; indxOut++)
    {
        outputDates->fArray[indxOut] = addDates[indx2++];
    }

    /* Update the item count */
    outputDates->fNumItems = indxOut;

    status = SUCCESS;

done:

    if (status != SUCCESS)
    {
        JpmcdsFreeDateList (outputDates);
        outputDates = NULL;
        JpmcdsErrMsgFailure (routine);
    }
    return outputDates;
}


/*
***************************************************************************
** Removes one date list from another date list, and returns a new date list.
** The input date lists may be unsorted. Output date list is sorted.
** Note the algorithm used may be inefficient if subDates is large.
***************************************************************************
*/
TDateList* JpmcdsSubtractDateLists
(TDateList *dateList,        /* (I) Original date list                     */
 TDateList *subDates         /* (I) date list to remove from original      */
)
{
    static char routine[] = "JpmcdsSubtractDateLists";
    int         status  = FALSE;

    int         indx1   = 0;
    int         indx2   = 0;
    int         indxOut = 0;
    long        maxDate = 0;
    TDateList  *unsortedDates = NULL;    /* before sorting                      */
    TDateList  *sortedDates = NULL;      /* after sorting and what we send back */

    if (dateList == NULL || subDates == NULL)
    {
        JpmcdsErrMsg ("%s: NULL inputs.\n", routine);
        goto done;
    }

    if (dateList->fNumItems < 0 || subDates->fNumItems < 0)
    {
        JpmcdsErrMsg ("%s: Negative number of dates specified.\n", routine);
        goto done;
    }

    /* This may be bigger than it has to be */
    unsortedDates = JpmcdsNewEmptyDateList(dateList->fNumItems);
    if (unsortedDates == NULL)
    {
        goto done;
    }

    /* store a date if in original dateList and not in subDates */
    for (indx1 =  indxOut = 0; indx1 < dateList->fNumItems; indx1++ )
    {
        for (indx2 = 0; indx2 < subDates->fNumItems; indx2++)
        {
            /* we have a match therefore do not store */
            if ( dateList->fArray[indx1] == subDates->fArray[indx2] )
                break;

            /* no match therefore store */
            if(indx2 == (subDates->fNumItems -1))
            {
                unsortedDates->fArray[indxOut] = dateList->fArray[indx1];
                indxOut++;
            }
        }

    }

    /* save the count */
    unsortedDates->fNumItems = indxOut;

    /* sort the dates */
    sortedDates = JpmcdsSortDateList(unsortedDates);
    if(sortedDates == NULL)
    {
         goto done;
    }

    /* now that dates are sorted truncate any surplus */
    for (indxOut = 0; indxOut < sortedDates->fNumItems; indxOut++)
    {
        /* find higest date */
        if(sortedDates->fArray[indxOut] > maxDate)
            maxDate = sortedDates->fArray[indxOut];
    }

    /* do the truncation in place */
    JpmcdsJpmcdsDateListTruncate(sortedDates, maxDate, FALSE, FALSE, TRUE);
    if(sortedDates == NULL)
    {
         goto done;
    }

    /* free up memory */
    JpmcdsFreeDateList (unsortedDates); unsortedDates = NULL;

    /* pat yourself on the back for success */
    status = SUCCESS;

done:

    if (status != SUCCESS)
    {
        JpmcdsFreeDateList (unsortedDates); unsortedDates = NULL;
        JpmcdsFreeDateList (sortedDates);     sortedDates = NULL;
        JpmcdsErrMsgFailure (routine);
    }
    return sortedDates;
}


/*
***************************************************************************
** Takes away the excess same date elements in a DateList
**
** Especially, useful to clean up a DateList after going through bad-day adj.
***************************************************************************
*/
TDateList* JpmcdsNewDateListDistinct
(
    TDateList *theDateList      /* (I) the input DateList to be cleaned */
)
{
    static char routine[]="JpmcdsNewDateListDistinct";
    TDateList *distinctDL = NULL;
    int i;                              /* Old index */
    int distIdx;                        /* Distinct index */

    distinctDL = JpmcdsNewEmptyDateList(theDateList->fNumItems);
    if (distinctDL == (TDateList *)NULL)
        goto done;

    /* Make sure it handles the caase where theDateList->fNumItems = 0.
     */
    for(i=0, distIdx = 0; i < distinctDL->fNumItems; i++)
    {
        if(i == 0 || theDateList->fArray[i] != theDateList->fArray[i-1] )
        {
           distinctDL->fArray[distIdx++] = theDateList->fArray[i];
        }
    }
    distinctDL->fNumItems = distIdx;

  done:
    if (distinctDL == NULL)
        JpmcdsErrMsg("%s: Failed.\n", routine);

    return distinctDL;
}


/*f
***************************************************************************
** Truncates a date list at the specified date. The resulting date list will
** contain all dates previous to (or following) the specified date. Dates in
** the datelist which match the specified date may be optionally included.
** The datelist may optionally be modified in place or a new copy is returned.
**
** The input date list must be sorted.
***************************************************************************
*/
TDateList* JpmcdsJpmcdsDateListTruncate
(TDateList *dateList,       /* (I/O) Date list to be modified in place */
 TDate      truncationDate, /* (I) Date on which to perform trunctation */
 TBoolean   inclusive,      /* (I) TRUE=include truncation date if in list*/
 TBoolean   excludeBefore,  /* (I) TRUE=exclude dates before truncation date */
 TBoolean   inPlace)        /* (I) TRUE=modify date list in place */
{
    static char  routine[] = "JpmcdsJpmcdsDateListTruncate";
    TDateList   *truncated = NULL;
    int          numItems;
    int          size;
    int          offset;
    int          truncatePt;
    int          i;

    if (dateList == NULL)
    {
        JpmcdsErrMsg ("%s: Null inputs.\n", routine);
        goto done;
    }

    /*  First we find the truncation point in the date list and the size
     *  of the new date list
     */
    numItems = dateList->fNumItems;
    if (excludeBefore)
    {
        truncatePt = numItems;
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
        truncatePt = -1;
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
        truncated = JpmcdsNewEmptyDateList(size);
        if (truncated == NULL)
        {
            goto done;
        }
    }

    /* finally we populate the result */
    if (inPlace == FALSE || offset != 0)
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
