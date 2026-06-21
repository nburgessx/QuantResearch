/*
 * ISDA CDS Standard Model
 *
 * Copyright (C) 2009 International Swaps and Derivatives Association, Inc.
 * Developed and supported in collaboration with Markit
 * 
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the ISDA CDS Standard Model Public License.
 */

#ifndef BASTYPES_H
#define BASTYPES_H

#include "cgeneral.h"
#include "cdate.h"
#include "mdydate.h"
#include "cfinanci.h"
#include "stub.h"

#ifdef __cplusplus
extern "C"
{
#endif


/*t
 */
typedef struct
{
    TDate  fDate;
    double fRate;
} TRatePt;


/*t
 * Holds a zero-coupon rate curve or clean spread curve.
 */
typedef struct _TCurve
{
    int       fNumItems;     /* Number of TRatePts in fArray */
    TRatePt  *fArray;        /* Dates & rates */
    TDate     fBaseDate;     /* Discount date */ 
    double    fBasis;        /* Number compounding periods / year */
    long      fDayCountConv; /* How the year fraction is computed */
} TCurve;


/*t
 * Defines a single fixed cashflow.
 */
typedef struct
{
    TDate  fDate;
    double fAmount;
} TCashFlow;


/*t
 * CashFlowList is a list of CashFlows (date and amount)
 */
typedef struct
{
    int        fNumItems;
    TCashFlow *fArray;
} TCashFlowList;


/*t
 * Defines dates required for one floating payment.
 */
typedef struct _TCouponDates
{
    TDate    resetDate;
    TDate    payDate;
    TDate    accrueStartDate;
    TDate    accrueEndDate;
    TDate    exDividendDate; 
} TCouponDates;


/*t
 * Defines dates needed for a list of floating payments.
 */
typedef struct _TCouponDateList
{
    int             fNumItems;
    TCouponDates   *fArray;
    long            stubLocation;
} TCouponDateList;


/*t 
 * TFloatRate defines a floating rate. For example, a 5 year semi-annual
 * swap rate, LIBOR or a compounding rate.
 * Note that payInterval is ONLY used if rateType = JPMCDS_SIMPLE_BASIS. 
 * This setting is used to compute LIBOR and par swap rates.
 * In order to compute compounding rates, set rateType = JPMCDS_ANNUAL_BASIS (1)
 * or 2 for semi-annual compounding, etc.
 * See the documentation for JpmcdsDayCountFraction for a list of possible 
 * day count convention constants. Note that there is no difference
 * between a TFloatRate which has a payInterval *equal* to the 
 * matInterval (for the same instrument) one which has a payInterval
 * which is *greater* than the matInterval, and one which has payInterval = 0.
 * This is because the instruments *always* have a payment at maturity, 
 * and also payments at intervals of payInterval, for dates on or 
 * before the maturity date. 
 * Simple LIBOR, etc (a zero coupon rate) is expressed by making
 * the payInterval = matInterval. 
 * Note that the spread is added to the rate *AFTER* multiplying by
 * the weight. 
 * Note that the holiday file and bad day convention in the
 * TDateAdjIntvl are used to adjust dates used to compute the rate. 
 */
typedef struct _TFloatRate
{
    TDateInterval   matInterval;        /* Time to maturity of rate */
    TDateInterval   payInterval;        /* Time between payments for rate */
    long            dayCountConv;       /* Day count convention of rate */
    TDateAdjIntvl   spotOffset;         /* From reset to rate effective date */
    double          spread;             /* Added to the rate  */
    long            rateType;           /* JPMCDS_SIMPLE_BASIS, JPMCDS_ANNUAL_BASIS*/
    double          weight;             /* Multiplied by rate */
} TFloatRate;                           



#ifdef __cplusplus
}
#endif

#endif    /* BASTYPES_H */


