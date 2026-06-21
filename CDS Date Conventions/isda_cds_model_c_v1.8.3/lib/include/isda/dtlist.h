/*
 * ISDA CDS Standard Model
 *
 * Copyright (C) 2009 International Swaps and Derivatives Association, Inc.
 * Developed and supported in collaboration with Markit
 * 
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the ISDA CDS Standard Model Public License.
 */

#ifndef DTLIST_H
#define DTLIST_H

#include "cgeneral.h"
#include "cdate.h"

#ifdef __cplusplus
extern "C"
{
#endif


/*f
***************************************************************************
** Creates a new empty TDateList.
**
** Note that it is permissible to create a TDateList with 0 dates.
***************************************************************************
*/
TDateList* JpmcdsNewEmptyDateList(int numDates);


/*f
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
);


/*f
***************************************************************************
** Allocates a new TDateList and copies supplied dates into it.
***************************************************************************
*/
TDateList* JpmcdsNewDateListFromDates
    (TDate      *dates,                 /* (I) */
     int         numDates);             /* (I) */


/*f
***************************************************************************
** Frees a DateList
***************************************************************************
*/
void JpmcdsFreeDateList
     (TDateList *);                      /* (I)  */


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
TDateList* JpmcdsMergeDateLists(
               TDateList *dateList1,     /* (I) */
               TDateList *dateList2);    /* (I) */


/*f
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
);


/*f
***************************************************************************
** Removes one date list from another date list, and returns a new date list.
** The input date lists may be unsorted. Output date list is sorted.
** Note the algorithm used may be inefficient if subDates is large.
***************************************************************************
*/
TDateList* JpmcdsSubtractDateLists
(TDateList *dateList,        /* (I) Original date list                     */
 TDateList *subDates         /* (I) date list to remove from original      */
);


/*f
***************************************************************************
** Makes a copy of a date list.
**
** Returns NULL on error, or if the original datelist was NULL.
***************************************************************************
*/
TDateList* JpmcdsCopyDateList(   /* duplicate a datelist */
      TDateList   *dl);                        /* (I) date list to copy */

/*f
***************************************************************************
** Sorts a datelist.
** 
** Returns NULL on FAILURE, and the sorted DateList on success.
***************************************************************************
*/
TDateList* JpmcdsSortDateList
     (TDateList *dl);           /* (I) */


/*f
***************************************************************************
** Takes away the excess same date elements in a DateList
**
** Especially, useful to clean up a DateList after going through bad-day adj.
***************************************************************************
*/
TDateList* JpmcdsNewDateListDistinct
(
    TDateList *theDateList      /* (I) Input DateList to be cleaned */
);


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
 TBoolean   inPlace         /* (I) TRUE=modify date list in place */
);


#ifdef __cplusplus
}
#endif

#endif    /* DTLIST_H */
