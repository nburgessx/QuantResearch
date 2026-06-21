/*
 * ISDA CDS Standard Model
 *
 * Copyright (C) 2009 International Swaps and Derivatives Association, Inc.
 * Developed and supported in collaboration with Markit
 * 
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the ISDA CDS Standard Model Public License.
 */

#ifndef DATELIST_H
#define DATELIST_H

#include "cgeneral.h"
#include "bastypes.h"
#include "dtlist.h"

#ifdef __cplusplus
extern "C"
{
#endif


/*f
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
    TBoolean stubAtEnd);                /* (I) T=Stub at end; F=Stub at beg. */


/*f
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
    TBoolean stubAtEnd);                /* (I) T=Stub at end; F=Stub at beg. */


/*f
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
TDateList *   JpmcdsNewDateListExtendedRoll(
    TDate startDate,                    /* (I) Start Date */
    TDate maturityDate,                 /* (I) Maturity Date */
    TDate rollDate,                     /* (I) Date for roll information */
    TDateInterval *interval,            /* (I) Increment */
    TBoolean stubAtEnd);                 /* (I) T=Stub at end; F=Stub at beg. */


/*f
***************************************************************************
** Allocates a new TDateList and copies supplied dates into it from 
** supplied TCurve.
**
** This provides a convenient way to isolate the dates of a TCurve.
***************************************************************************
*/
TDateList* JpmcdsNewDateListFromTCurve
    (TCurve     *curve);        /* (I) Given curve */


/*f
***************************************************************************
** Allocates a new TDateList by calling JpmcdsNewDateList, and then removing 
** the startDate.
***************************************************************************
*/
TDateList* JpmcdsNewPayDates
    (TDate          startDate,          /* (I) This date is not included */
     TDate          matDate,            /* (I) */
     TDateInterval *payInterval,        /* (I) */
     TBoolean       stubAtEnd);         /* (I) */


/*f
***************************************************************************
** Creates a new business-day adjusted dateList from an UN-adjusted
** datelist, using the supplied busisness day adjustment convention.
***************************************************************************
*/
TDateList* JpmcdsNewDateListBusDayAdj
(
    TDateList *dateList,        /* (I) Unadjusted date list  */ 
    long       badDayConv,      /* (I) See JpmcdsBusinessDay */
    char      *holidayFile      /* (I) Holiday date file */
);


/*f
***************************************************************************
** Adjusts a datelist (in-place) according to the supplied business day 
** adjustment convention.
***************************************************************************
*/
int JpmcdsDateListBusDayAdj
(
    TDateList *dateList,        /* (I/O) Date list to adjust (in-place) */ 
    long       badDayConv,      /* (I) See JpmcdsBusinessDay */
    char      *holidayFile      /* (I) Holiday date file */
);


#ifdef __cplusplus
}
#endif

#endif    /* DATELIST_H */

