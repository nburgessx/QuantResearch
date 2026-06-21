/*
 * ISDA CDS Standard Model
 *
 * Copyright (C) 2009 International Swaps and Derivatives Association, Inc.
 * Developed and supported in collaboration with Markit
 * 
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the ISDA CDS Standard Model Public License.
 */

#ifndef CX_DATELIST_H
#define CX_DATELIST_H

#include "datelist.h"


#ifdef __cplusplus
extern "C"
{
#endif


/*f
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
 TDate      *array);     /* (I) [numItems] Dates to be added */


/*f
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
 TDate     *array);     /* (I) [numItems] Dates to be added     */


/*f
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
 TBoolean   inclusive,      /* (I) TRUE=include truncation date if in list*/
 TBoolean   excludeBefore,  /* (I) TRUE=exclude dates before truncation date*/
 TBoolean   inPlace         /* (I) TRUE=modify date list in place */
);


/*f
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
 TStubMethod   *stubType);  /* (I) Stub type */


#ifdef __cplusplus
}
#endif

#endif

