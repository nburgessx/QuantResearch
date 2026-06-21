/*
 * ISDA CDS Standard Model
 *
 * Copyright (C) 2009 International Swaps and Derivatives Association, Inc.
 * Developed and supported in collaboration with Markit
 * 
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the ISDA CDS Standard Model Public License.
 */

#ifndef SCHEDULE_H
#define SCHEDULE_H

#include "bastypes.h"
#include "convert.h"
#include "ldate.h"

#ifdef __cplusplus
extern "C"
{
#endif


typedef enum
{
    JPMCDS_SINGLE_REFIX,
} TSwapType;


typedef enum
{
    JPMCDS_NO_STUB,
    JPMCDS_BACK_STUB,
    JPMCDS_FRONT_STUB,
    JPMCDS_TWO_STUB
} TStubType;

 

/*f
***************************************************************************
** Creates a date schedule for one stream of a swap.
***************************************************************************
*/
TCouponDateList* JpmcdsNewCouponDatesSwap(
    TDate           startDate,          /* (I) Effective date of stream     */
    TDate           matDate,            /* (I) Maturity date of stream      */
    TDateInterval   *couponInterval,    /* (I) Time between payments        */
    TBoolean        adjustLastAccDate,  /* (I) Maturity date adjusted       */
    TBoolean        inArrears,          /* (I) Arrears setting flag         */
    int             payOffsetDays,      /* (I) Days from acc end to payment */
    int             resetOffsetDays,    /* (I) Days from accrual to reset   */
    TBoolean        stubAtEnd,          /* (I) TRUE=back stub,FALSE=front stub*/
    TBoolean        longStub,           /* (I) TRUE=long stub,FALSE=short stub*/
    TDate           firstRollDate,      /* (I) Date of first cashflow       */
    TDate           lastRollDate,       /* (I) Acc start of final period    */
    TDate           fullFirstCoupDate,  /* (I) Acc start of first period    */
    long            payBadDayConv,      /* (I) Bad day convention for payment*/
    long            accBadDayConv,      /* (I) Bad day convention for accrual*/
    long            resetBadDayConv,    /* (I) Bad day convention for resets*/
    char            *holidayFile);      /* (I) See JpmcdsBusinessDay           */


/*f
***************************************************************************
** Creates a new empty TCouponDateList.
***************************************************************************
*/
TCouponDateList* JpmcdsNewEmptyCouponDateList(int numPeriods);


/*f
***************************************************************************
** Frees a TCouponDateList.
***************************************************************************
*/
void JpmcdsFreeCouponDateList
       (TCouponDateList *theSchedule); /* (I)  */


/*f
***************************************************************************
** Checks the details for a stream schedule.
***************************************************************************
*/
int JpmcdsCouponDateListCheck(TCouponDateList *schedule);


/*f
***************************************************************************
** Generate a list of unadjusted coupon dates for the stream.
***************************************************************************
*/
int JpmcdsGenerateUnadjDates(
    TDate           startDate,
    TDate           matDate,
    TDateInterval   *couponInterval,
    TBoolean        stubAtEnd,
    TBoolean        longStub,
    TDate           firstRollDate,
    TDate           lastRollDate,
    TDateList       **dateList,
    long            *stubInfo);


#ifdef __cplusplus
}
#endif

#endif    /* SCHEDULE_H */

