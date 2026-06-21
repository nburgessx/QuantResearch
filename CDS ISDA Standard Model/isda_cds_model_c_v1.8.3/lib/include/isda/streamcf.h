/*
 * ISDA CDS Standard Model
 *
 * Copyright (C) 2009 International Swaps and Derivatives Association, Inc.
 * Developed and supported in collaboration with Markit
 * 
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the ISDA CDS Standard Model Public License.
 */

#ifndef STREAMCF_H
#define STREAMCF_H

#include "bastypes.h"
#include "ldate.h"
#include "schedule.h"
#include "cashflow.h"


#ifdef __cplusplus
extern "C"
{
#endif

typedef struct _TStubRates
{
    double          stubRate;
    TDateInterval   stubInterval;
} TStubRates;


typedef struct _TStreamFloat
{
    TCouponDateList  *dl;
    TSwapType        swapType;
    double           principal;
    long             payDayCountConv;
    long             flags;
    TStubRates       firstCoupon;
    TStubRates       finalCoupon;
    double           averSoFar;
    double           compSoFar;
    TFloatRate       *floatingRate;
} TStreamFloat;

typedef struct _TStreamFixed
{
    TCouponDateList  *dl;
    TSwapType        swapType;
    double           principal;
    long             payDayCountConv;
    long             flags;
    double           firstCouponRate;
    double           finalCouponRate;
    double           fixedRate;
} TStreamFixed;


/*f
***************************************************************************
** Calculates the cashflows for a floating stream of a swap.
**
** Principals are included as per specification in the stream->flags field.
***************************************************************************
*/
TCashFlowList* JpmcdsNewStreamFloatCFLGen(
     TCurve             *zeroCurve,          /* (I) Zero Curve              */
     long                zcInterpType,       /* (I) Zero curve interp type  */
     TStreamFloat       *stream,             /* (I) See JpmcdsNewStreamFloat   */
     long               rateBadDayConv,      /* (I) Bad Day Conv for index  */
     char               *holidayFile);       /* (I) Holiday file            */


/*f
***************************************************************************
** Calculates the cashflows for a fixed stream of a swap.
**
** Principals are included as per specification in the stream->flags field.
***************************************************************************
*/
TCashFlowList* JpmcdsNewStreamFixedCFL(
     TDate              valueDate,           /* (I) Value Date              */
     TStreamFixed       *stream);            /* (I) See JpmcdsNewStreamFixed   */


/*f
***************************************************************************
** Given one stream of a swap, this function estimates the rates for each 
** coupon.
***************************************************************************
*/
int JpmcdsGetFixedCouponRates(
     TDate              valueDate,           /* (I) Value Date              */
     TStreamFixed       *stream,             /* (I) See JpmcdsMakeTStram       */
     double             *couponRates);       /* (O) Coupon Rates            */


/*f
***************************************************************************
** Makes a new TStreamFloat.
***************************************************************************
*/
TStreamFloat* JpmcdsNewStreamFloat(
     TCouponDateList   *dl,                  /* (I) See JpmcdsSwapMakeSchedule */
     TSwapType          swapType,            /* (I) The type of swap        */
     double             principal,           /* (I) Principal amount        */
     long               payDayCountConv,     /* (I) See JpmcdsDayCountConvention*/
     long               flags,               /* (I) Principal exchanges     */
     TStubRates        *firstCoupon,         /* (I) Data for first coupon   */
     TStubRates        *finalCoupon,         /* (I) Data for final coupon   */
     double             averSoFar,           /* (I) Known avg rate portion  */
     double             compSoFar,           /* (I) Known compounded rate   */
     TFloatRate        *floatingRate);       /* (I) Defn of floating index  */


/*f
***************************************************************************
** Makes a new TStreamFixed.
***************************************************************************
*/
TStreamFixed* JpmcdsNewStreamFixed(
     TCouponDateList   *dl,                  /* (I) See JpmcdsSwapMakeSchedule */
     TSwapType          swapType,            /* (I) The type of swap        */
     double             principal,           /* (I) Principal amount        */
     long               payDayCountConv,     /* (I) See JpmcdsDayCountConvention*/
     long               flags,               /* (I) Principal exchanges     */
     double             firstCoupon,         /* (I) Data for first coupon   */
     double             finalCoupon,         /* (I) Data for final coupon   */
     double             fixedRate);          /* (I) Fixed Coupon Rate       */


/*f
***************************************************************************
** Creates a new empty TStreamFloat.
***************************************************************************
*/
TStreamFloat* JpmcdsNewEmptyTStreamFloat(int  numPeriods);


/*f
***************************************************************************
** Creates a new empty TStreamFixed.
***************************************************************************
*/
TStreamFixed* JpmcdsNewEmptyTStreamFixed(int  numPeriods);


/*f
***************************************************************************
** Frees a floating stream. 
***************************************************************************
*/
void JpmcdsFreeTStreamFloat(TStreamFloat *theStream);


/*f
***************************************************************************
** Frees a fixed stream. 
***************************************************************************
*/
void JpmcdsFreeTStreamFixed(TStreamFixed *theStream);


/*f
***************************************************************************
** Same as JpmcdsFloatingCashFlow, except user can choose the zero curve 
** interpolation type. 
***************************************************************************
*/
int JpmcdsFloatingCashFlowGen
   (TCurve         *zeroCurve,           /* (I) Zero curve                  */
    long            zcInterpType,        /* (I) Zero curve interp type      */
    TDate          indexStartDate,       /* (I) Start date of rate          */
    TDate          indexMatDate,         /* (I) Maturity date of rate       */
    TDateInterval  payInterval,          /* (I) Payment interval of rate    */
    TCouponDates   *couponDates,         /* (I) AccStart,AccEnd,Pay,Reset   */
    long           rateDayCountConv,     /* (I) Market convention of rate   */
    long           payDayCountConv,      /* (I) Convention for payments     */
    double         notional,             /* (I) Notional amount             */
    double         spread,               /* (I) Added to the coupon rate    */
    double         *cashflow);           /* (O) Coupon payment              */

     
#ifdef __cplusplus
}
#endif

#endif    /* STREAMCF_H */
    
