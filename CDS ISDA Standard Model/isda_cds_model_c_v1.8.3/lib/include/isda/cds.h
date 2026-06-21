/*
 * ISDA CDS Standard Model
 *
 * Copyright (C) 2009 International Swaps and Derivatives Association, Inc.
 * Developed and supported in collaboration with Markit
 * 
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the ISDA CDS Standard Model Public License.
 */

#ifndef CDS_H
#define CDS_H

#include "cx.h"


#ifdef __cplusplus
extern "C"
{
#endif

/*f
***************************************************************************
** Makes a contingent leg for a vanilla CDS.
**
** The CDS starts at startDate and ends at endDate. The last date is always
** protected - the start date is only protected if protectStart=True.
***************************************************************************
*/
TContingentLeg* JpmcdsCdsContingentLegMake(
    /** Date when protection begins. Either at start or end of day (depends
        on protectStart) */
    TDate     startDate,
    /** Date when protection ends (end of day) */
    TDate     endDate,
    /** Notional value protected */
    double    notional,
    /** Should protection include the start date */
    TBoolean  protectStart);


/*f
***************************************************************************
** Computes the PV for a contingent leg for a vanilla CDS.
**
** Risk starts at the end of today. The PV is computed for a given value date.
** The CDS starts at startDate and ends at endDate. The last date is always
** protected - the start date is only protected if protectStart=True.
***************************************************************************
*/
int JpmcdsCdsContingentLegPV(
    /** Risk starts at the end of today */
    TDate             today,
    /** Date for which the PV is calculated and cash settled */
    TDate             valueDate,
    /** Date when protection begins. Either at start or end of day (depends
        on protectStart) */
    TDate             startDate,
    /** Date when protection ends (end of day) */
    TDate             endDate,
    /** Notional value protected */
    double            notional, 
    /** Interest rate discount curve - assumes flat forward interpolation */
    TCurve           *discCurve,
    /** Credit clean spread curve */
    TCurve           *spreadCurve,
    /** Assumed recovery rate in case of default */
    double            recoveryRate,
    /** True => protection includes start date */
    TBoolean          protectStart,
    /** Output - the present value is returned */
    double           *pv);


/*f
***************************************************************************
** Makes a fixed fee leg for a vanilla CDS.
**
** The CDS starts at startDate and ends at endDate. The last date is always
** protected - the start date is only protected if protectStart=True.
**
** Interest accrues for the same number of days as there is protection.
** Thus if protectStart=True you get one extra day of accrued interest in
** comparison with an interest rate swap. This extra day is assumed to be
** the last day of the CDS and means that the last period is one day longer
** than for an interest rate swap.
***************************************************************************
*/
TFeeLeg* JpmcdsCdsFeeLegMake(
    /** Date when protection begins. Either at start or end of day (depends
        on protectStart) */
    TDate           startDate,
    /** Date when protection ends (end of day) */
    TDate           endDate,
    /** Should accrued interest be paid on default. Usually set to TRUE */
    TBoolean        payAccOnDefault,
    /** Interval between coupon payments. Can be NULL when 3M is assumed */
    TDateInterval  *couponInterval,
    /** If the startDate and endDate are not on cycle, then this parameter
        determines location of coupon dates. */
    TStubMethod    *stubType,
    /** Notional value protected */
    double          notional, 
    /** Fixed coupon rate (a.k.a. spread) for the fee leg */
    double          couponRate,
    /** Day count convention for coupon payment. Normal is ACT_360 */
    long            paymentDcc,
    /** Bad day convention for adjusting coupon payment dates. */
    long            badDayConv,
    /** Calendar used when adjusting coupon dates. Can be NULL which equals
        a calendar with no holidays and including weekends. */
    char           *calendar,
    /** Should protection include the start date */
    TBoolean        protectStart);


/*f
***************************************************************************
** Computes the PV for a fixed fee leg for a vanilla CDS.
**
** Risk starts at the end of today. The PV is computed for a given value date.
** The CDS starts at startDate and ends at endDate. The last date is always
** protected - the start date is only protected if protectStart=True.
**
** Interest accrues for the same number of days as there is protection.
** Thus if protectStart=True you get one extra day of accrued interest in
** comparison with an interest rate swap. This extra day is assumed to be
** the last day of the CDS and means that the last period is one day longer
** than for an interest rate swap.
***************************************************************************
*/
int JpmcdsCdsFeeLegPV(
    /** Risk starts at the end of today */
    TDate           today,
    /** Date for which the PV is calculated and cash settled */
    TDate           valueDate,
    /** Date when step-in becomes effective */
    TDate           stepinDate,
    /** Date when protection begins. Either at start or end of day (depends
        on protectStart) */
    TDate           startDate,
    /** Date when protection ends (end of day) */
    TDate           endDate,
    /** Should accrued interest be paid on default. Usually set to TRUE */
    TBoolean        payAccOnDefault,
    /** Interval between coupon payments. Can be NULL when 3M is assumed */
    TDateInterval  *couponInterval,
    /** If the startDate and endDate are not on cycle, then this parameter
        determines location of coupon dates. */
    TStubMethod    *stubType,
    /** Notional value protected */
    double          notional, 
    /** Fixed coupon rate (a.k.a. spread) for the fee leg */
    double          couponRate,
    /** Day count convention for coupon payment. Normal is ACT_360 */
    long            paymentDcc,
    /** Bad day convention for adjusting coupon payment dates. */
    long            badDayConv,
    /** Calendar used when adjusting coupon dates. Can be NULL which equals
        a calendar with no holidays and including weekends. */
    char           *calendar,
    /** Interest rate discount curve - assumes flat forward interpolation */
    TCurve         *discCurve,
    /** Credit clean spread curve */
    TCurve         *spreadCurve,
    /** Should protection include the start date */
    TBoolean        protectStart,
    /** Should the present value be computed as a clean price (removing
        accrued interest) */
    TBoolean        isPriceClean,
    /** Output - the present value is returned */
    double         *pv);


/*f
***************************************************************************
** Computes the price (a.k.a. upfront charge) for a vanilla CDS
**
** Risk starts at the end of today. The PV is computed for a given value date.
** The CDS starts at startDate and ends at endDate. The last date is always
** protected - the start date is also protected because internally protectStart
** is set to True.
**
** Interest accrues for the same number of days as there is protection.
** Since protectStart is set to True you get one extra day of accrued interest in
** comparison with an interest rate swap. This extra day is assumed to be
** the last day of the CDS and means that the last period is one day longer
** than for an interest rate swap.
***************************************************************************
*/
EXPORT int JpmcdsCdsPrice(
    /** Risk starts at the end of today */
    TDate           today,
    /** Date for which the PV is calculated and cash settled */
    TDate           valueDate,
    /** Date when step-in becomes effective */
    TDate           stepinDate,
    /** Date when protection begins. Either at start or end of day (depends
        on protectStart) */
    TDate           startDate,
    /** Date when protection ends (end of day) */
    TDate           endDate,
    /** Fixed coupon rate (a.k.a. spread) for the fee leg */
    double          couponRate,
    /** Should accrued interest be paid on default. Usually set to TRUE */
    TBoolean        payAccOnDefault,
    /** Interval between coupon payments. Can be NULL when 3M is assumed */
    TDateInterval  *couponInterval,
    /** If the startDate and endDate are not on cycle, then this parameter
        determines location of coupon dates. */
    TStubMethod    *stubType,
    /** Day count convention for coupon payment. Normal is ACT_360 */
    long            paymentDcc,
    /** Bad day convention for adjusting coupon payment dates. */
    long            badDayConv,
    /** Calendar used when adjusting coupon dates. Can be NULL which equals
        a calendar with no holidays and including weekends. */
    char           *calendar,
    /** Interest rate discount curve - assumes flat forward interpolation */
    TCurve         *discCurve,
    /** Credit clean spread curve */
    TCurve         *spreadCurve,
    /** Assumed recovery rate in case of default */
    double          recoveryRate,
    /** Is the price expressed as a clean price (removing accrued interest) */
    TBoolean        isPriceClean,
    /** Output - price (a.k.a. upfront charge) for the CDS is returned 
        (see also isPriceClean) */
    double         *price);

  
/*f
***************************************************************************
** Computes the benchmark par spreads 
**
** Risk starts at the end of today. The PV is computed for a given value date.
** The CDS starts at startDate and ends at endDate. The last date is always
** protected - the start date is also protected because internally protectStart
** is set to True.
**
** Interest accrues for the same number of days as there is protection.
** Since protectStart is set to True you get one extra day of accrued interest in
** comparison with an interest rate swap. This extra day is assumed to be
** the last day of the CDS and means that the last period is one day longer
** than for an interest rate swap.
***************************************************************************
*/
EXPORT int JpmcdsCdsParSpreads(
    /** Risk starts at the end of today */
    TDate           today,
    /** Date when step-in becomes effective  */
    TDate           stepinDate,
    /** Date when protection begins. Either at start or end of day (depends
        on protectStart) */
    TDate           startDate,
    /** Number of benchmark dates */
    long            nbEndDates,
    /** Date when protection ends (end of day), no bad day adjustment */
    TDate          *endDates,
    /** Should accrued interest be paid on default. Usually set to TRUE */
    TBoolean        payAccOnDefault,
    /** Interval between coupon payments. Can be NULL when 3M is assumed */
    TDateInterval  *couponInterval,
    /** If the startDate and endDate are not on cycle, then this parameter
        determines location of coupon dates. */
    TStubMethod    *stubType,
    /** Day count convention for coupon payment. Normal is ACT_360 */
    long            paymentDcc,
    /** Bad day convention for adjusting coupon payment dates. */
    long            badDayConv,
    /** Calendar used when adjusting coupon dates. Can be NULL which equals
        a calendar with no holidays and including weekends. */
    char           *calendar,
    /** Interest rate discount curve - assumes flat forward interpolation */
    TCurve         *discCurve,
    /** Credit clean spread curve */
    TCurve         *spreadCurve,
    /** Assumed recovery rate in case of default */
    double          recoveryRate,
    /** Output - par spreads for the CDS are returned (see also isPriceClean) */
    double         *parSpread);


/*f
***************************************************************************
** Computes the non-contingent cash flows for a fee leg. These are the
** cash flows you will receive if there is no default.

** Interest accrues for the same number of days as there is protection.
** Since protectStart is set to True you get one extra day of accrued interest in
** comparison with an interest rate swap. This extra day is assumed to be
** the last day of the CDS and means that the last period is one day longer
** than for an interest rate swap.
***************************************************************************
*/
EXPORT TCashFlowList* JpmcdsCdsFeeLegFlows(
    /** Date when protection begins. Either at start or end of day (depends
        on protectStart) */
    TDate           startDate,
    /** Date when protection ends for each benchmark (end of day).*/
    TDate           endDate,
    /** Interval between coupon payments. Can be NULL when 3M is assumed */
    TDateInterval  *dateInterval,
    /** If the startDate and endDate are not on cycle, then this parameter
        determines location of coupon dates. */
    TStubMethod    *stubType,
    /** Notional of the fee leg */
    double          notional,
    /** Fixed coupon rate (a.k.a. spread) for the fee leg */
    double          couponRate,
    /** Day count convention for coupon payment. Normal is ACT_360 */
    long            paymentDcc,
    /** Bad day convention for adjusting coupon payment dates. */
    long            badDayConv,
    /** Calendar used when adjusting coupon dates. Can be NULL which equals
        a calendar with no holidays and including weekends. */
    char           *calendar);


/*f
***************************************************************************
** Bootstraps a clean spread curve from par spread inputs.
**
** Risk starts at the end of today. The PV is computed for a given value date.
** The CDS starts at startDate and ends at endDate. The last date is always
** protected - the start date is also protected because internally protectStart
** is set to True.
**
** Interest accrues for the same number of days as there is protection.
** Since protectStart is set to True you get one extra day of accrued interest in
** comparison with an interest rate swap. This extra day is assumed to be
** the last day of the CDS and means that the last period is one day longer
** than for an interest rate swap.
***************************************************************************
*/
EXPORT TCurve* JpmcdsCleanSpreadCurve(
    /** Risk starts at the end of today */
    TDate           today,
    /** Interest rate discount curve - assumes flat forward interpolation */
    TCurve         *discCurve,
    /** Effective date of the benchmark CDS */
    TDate           startDate,
    /** Step in date of the benchmark CDS */
    TDate           stepinDate,
    /** Date when payment should be make */
    TDate           cashSettleDate,
    /** Number of benchmark dates */
    long            nbDate,
    /** Dates when protection ends for each benchmark (end of day).
        Array of size nbDate */
    TDate          *endDates,
    /** Coupon rates for each benchmark instrument. Array of size nbDate */
    double         *couponRates,
    /** Flags to denote that we include particular benchmarks. This makes it
        easy for the user to include or exclude benchmarks on a one-by-one
        basis. Can be NULL if all are included. Otherwise an array of size
        nbDate. */
    TBoolean       *includes,
    /** Recovery rate in case of default */
    double          recoveryRate,
    /** Should accrued interest be paid on default. Usually set to TRUE */
    TBoolean        payAccOnDefault,
    /** Interval between coupon payments. Can be NULL when 3M is assumed */
    TDateInterval  *couponInterval,
    /** Day count convention for coupon payment. Normal is ACT_360 */
    long            paymentDcc,
    /** If the startDate and endDate are not on cycle, then this parameter
        determines location of coupon dates. */
    TStubMethod    *stubType,
    /** Bad day convention for adjusting coupon payment dates. */
    long            badDayConv,
    /** Calendar used when adjusting coupon dates. Can be NULL which equals
        a calendar with no holidays and including weekends. */
    char           *calendar
);

#ifdef __cplusplus
}
#endif

#endif

