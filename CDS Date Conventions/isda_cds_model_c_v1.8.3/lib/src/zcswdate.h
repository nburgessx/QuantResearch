/*
 * ISDA CDS Standard Model
 *
 * Copyright (C) 2009 International Swaps and Derivatives Association, Inc.
 * Developed and supported in collaboration with Markit
 * 
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the ISDA CDS Standard Model Public License.
 */

#ifndef ZCSWDATE_H
#define ZCSWDATE_H

#include "cgeneral.h"
#include "bastypes.h"
#include "badday.h"

#ifdef __cplusplus
extern "C"
{
#endif

typedef struct _swapDates
{
    TDate    *adjusted;                 /* Adjusted dates */
    TDate    *original;                 /* UNadjusted dates */
    TDate    *previous;                 /* Previous coupon dates */
    TBoolean *onCycle;                  /* T=Count fwd from Val date */
                                        /* F=Count bckwd from Mat date */
    int       numDates;                 /* Length of all 4 arrays. */
} TSwapDates;


/*f
***************************************************************************
** Constructor for TSwapDates.
***************************************************************************
*/
TSwapDates* JpmcdsSwapDatesNewEmpty(int numDates);


/*f
***************************************************************************
** Uses swap maturity dates adjusted for bad days and a TBadDayList to set 
** up a TSwapDates.
***************************************************************************
*/
TSwapDates* JpmcdsSwapDatesNewFromAdjusted
     (TDate          valueDate,         /* (I) Must be a good day. */
      int            freq,              /* (I) Used to set prevDates */
      TDate         *adjustedDates,     /* (I) Adjusted swap mat dates */
      int            numDates,          /* (I) Length of adjustedDates */
      TBadDayList   *badDayList);       /* (I) Maps good to bad & vice-vsa*/


/*f
***************************************************************************
** Uses un-adjusted swap maturity dates, a holiday file, and a bad day 
** convention to set up a TSwapDates.
***************************************************************************
*/
TSwapDates* JpmcdsSwapDatesNewFromOriginal
     (TDate          valueDate,         /* (I) Must be a good day. */
      int            freq,              /* (I) Used to set prevDates */
      TDate         *originalDates,     /* (I) UN-Adjusted swap mat dates */
      int            numDates,          /* (I) Length of originalDates */
      TBadDayList   *badDayList,        /* (I) Maps good to bad & vice-vsa*/
      long           badDayConv,        /* (I) See JpmcdsIsBusinessDay */
      char          *holidayFile);      /* (I) List of holidays */


/*f
***************************************************************************
** Destructor for TSwapDates.
***************************************************************************
*/
void JpmcdsSwapDatesFree(TSwapDates *sd);


/*f
***************************************************************************
** Adjusts date to fall on a business day, if necesssary.
***************************************************************************
*/
int JpmcdsZCAdjustDate
   (TDate       date,           /* (I) Date to be adjusted */
    TBadDayList *badDayList,    /* (I) */
    long         badDayConv,    /* (I) */
    char        *holidayFile,   /* (I) */
    TDate       *adjDate);      /* (O) */



#ifdef __cplusplus
}
#endif

#endif    /* ZCSWDATE_H */

