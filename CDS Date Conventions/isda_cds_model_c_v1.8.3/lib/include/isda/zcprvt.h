/*
 * ISDA CDS Standard Model
 *
 * Copyright (C) 2009 International Swaps and Derivatives Association, Inc.
 * Developed and supported in collaboration with Markit
 * 
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the ISDA CDS Standard Model Public License.
 */

#ifndef ZCPRVT_H
#define ZCPRVT_H

#include "cgeneral.h"
#include "bastypes.h"
#include "badday.h"
#include "tcurve.h"
#include "interp.h"

#ifdef __cplusplus
extern "C"
{
#endif

typedef struct {                       /* stores curve of zero rates */
   TDate           valueDate;          /* base date all rates start at */
   double         *rate;               /* array of rates */
   TDate          *date;               /* array of dates */
   double         *discount;           /* array of discounts */
   int             numItems;           /* number of valid rates in curve */
   int             numAlloc;           /* number of items allowed in arrays */
   long            basis;              /* compounding basis: usually annual */
   long            dayCountConv;       /* day counting convention */
} ZCurve;


typedef struct {                       /* data for forward stub methods */
   double   userRate;                  /* user rate: annualized rate */
   long     userBasis;                 /* user rate: compounding basis */
   long     userDayCountConv;          /* user rate: day counting conv */
   TDate   *interpDates;               /* array of TDates for interping */
   double  *interpRates;               /* array of rates for interping */
   int      interpN;                   /* number of elements for interping*/
   double   rate;                      /* euro-money rate (simple annual) */
} TStubData;


typedef struct TInterpDataTag
{                                       /* interpolation data type */
   TDateInterval   interval;            /* linear FWDs: length of forward */
   long            dayCountConv;        /*    ""   ""   day counting conv */
   long            basis;               /*    ""   ""   compounding basis */
   TBadDayList    *badDayList;          /*    ""   ""   bad day list */
   TDateInterval   addPointsIntvls[5];  /* intervals to add points between
                                         * start,1y,2y,5y,10y,end */
   TBoolean        enableGeneration;    /* enable generation of zc points 
                                         * at coupon dates */
} TInterpData;

/*f
***************************************************************************
** Constructs a TCurve from a ZCurve.
***************************************************************************
*/
TCurve* JpmcdsZCToTCurve(
   ZCurve  *zc);                           /* (I) ZCurve to copy */


/*f
***************************************************************************
** Constructs a ZCurve from a TCurve.
***************************************************************************
*/
ZCurve* JpmcdsZCFromTCurve(
   TCurve  *tc);                           /* (I) TCurve to copy */


/*f
***************************************************************************
** Constructs a ZCurve, which is basically an array of <rate,date> pairs, 
** where the rates are all zero-coupons with the same compounding and day 
** count convention.  This allows a bunch of different types of interest rate
** instruments to be combined into one yield curve in the same terms.
**
***************************************************************************
*/
ZCurve* JpmcdsZCMake(
   TDate           valueDate,              /* (I) base date for rates */
   int             n,                      /* (I) expected number of rates */
   long            basis,                  /* (I) compounding period */
   long            dayCountConv);          /* (I) day count convention */


/*f
***************************************************************************
** Frees memory associated with a zero curve.
***************************************************************************
*/
void JpmcdsZCFree(
   ZCurve                 *zc);            /* zero curve to free */


/*f
***************************************************************************
** Adds points to a zero curve from another zero curve, but only those 
** dates before first date of the other curve.
**
** ZCurve is updated to include the new points.
**
** Note:curve value dates MUST be the same, basis and dayCounts may differ.
***************************************************************************
*/
int JpmcdsZCAddPrefixCurve(
   ZCurve                     *zc1,        /* (M) zCurve to added to */
   ZCurve                     *zc2);       /* (I) zCurve to add from */


/*f
***************************************************************************
** Adds points to a zero curve from another zero curve, but only those 
** dates after last date of the other curve.
**
** ZCurve is updated to include the new points.
**
** Note: curve value dates MUST be the same, basis and dayCounts may differ.
***************************************************************************
*/
int JpmcdsZCAddSuffixCurve(
   ZCurve                     *zc1,        /* (M) zCurve to added to */
   ZCurve                     *zc2);       /* (I) zCurve to add from */


/*f
***************************************************************************
** Adds a single swap instrument to a ZCurve.
**
** ZCurve is updated to contain swap information (with points at each of 
** the coupon dates, as well as the maturity date).
***************************************************************************
*/
int JpmcdsZCAddSwap(
   ZCurve              *zc,                /* (M) ZCurve to add to          */
   TCurve              *discZC,            /* (I) ZC used for discounting   */
   double              price,              /* (I) Par price: usually 1.0    */
   TDate               matDate,            /* (I) Unadjusted Maturity date  */
   TBoolean            onCycle,            /* (I) If matDate on cycle       */
   double              rate,               /* (I) Coupon rate (0.06 for 6%) */
   int                 fixedSwapFreq,      /* (I) Fixed coupon frequency    */
   int                 floatSwapFreq,      /* (I) Floating coupon frequency */
   long                fixDayCountConv,    /* (I) Fixed leg convention      */
   long                floatDayCountConv,  /* (I) Floating leg convention   */
   long                interpType,         /* (I) Interpolation method type */
   TInterpData         *interpData,        /* (I) Interpolation data        */
   TBadDayList         *badDayList,        /* (I) List of bad-business days */
   TBadDayAndStubPos   badDayAndStubPos,   /* (I) See JpmcdsBusinessDay        */
   char                *holidayFile);      /* (I) See JpmcdsBusinessDay        */


/*f
***************************************************************************
** Adds a strip of swap instruments to a ZCurve.
**
** The zero interpolation kicks in whenever a point is calculated to
** determine the discount factor for a cash flow, deep inside the objective
** function called by the root finder, which calls JpmcdsZCInterpolate() to 
** determine the appropriate discount factor.  
**
** In this routine, there are two mutually exclusive ways to specify 
** bad days: badDayList, and badDayConv/holidayFile. If badDayList is used, 
** bad business day lists can be used to make sure swap dates near to
** canonical dates (i.e. integral number of years from start date) are
** properly used for coupon dates of later instruments.  Say the 2-year
** point has been adjusted because of holidays to be 3 days later.  The
** coupon at 2-years for the 3 (and subsequent) swap instruments must be
** adjusted to that date.  To do this, the badDayList contains a list of
** <good,bad> day pairs, which are used to adjust dates to/from bad
** business days.
**
** INPUTS: ZCurve to add swaps to, maturity dates and rates for them, how many
** there are, coupon frequency: annual or semi-annual, dayCountConvention,
** interpolation method and data, bad business day list. 
** WARNING: In this routine, there are two mutually exclusive ways to specify
** bad days: badDayList, and badDayConv/holidayFile. 
** If badDayList is non-NULL, the maturity dates passed in must be
** adjusted. If the badDayConv is not JPMCDS_BAD_DAY_NONE, the maturity dates
** passed in must be UN-adjusted. 
**
** OUTPUTS: Extended zero curve. 
***************************************************************************
*/
int JpmcdsZCAddSwaps(
 ZCurve          *zc,               /* (M) ZCurve to add to        */
 TCurve          *discZC,           /* (I) Discount zero curve     */
 TDate           *inDates,          /* (I) MatDates; adj. if badDayList!=0*/
 double          *inRates,          /* (I) Fixed rates (0.06 for 6%) */
 int              numSwaps,         /* (I) # instruments to add    */
 int              fixedSwapFreq,    /* (I) Fixed leg freq          */
 int              floatSwapFreq,    /* (I) Floating leg freq       */
 long             fixDayCountConv,  /* (I) Convention for fixed leg*/
 long             floatDayCountConv,/* (I) Convention for float leg*/
 long             interpType,       /* (I) Zero interpolation method */
 TInterpData     *interpData,       /* (I) Zero interpolation data */
 TBadDayList     *badDayList,       /* (I) Bad-good day pairs; if non-NULL */
                                    /*     maturity dates must be adjusted */
 TBadDayAndStubPos badDayAndStubPos,/* (I) See JpmcdsBusinessDay      */
 char            *holidayFile);     /* (I) See JpmcdsBusinessDay      */


/*f
***************************************************************************
** Makes a date list for all coupons associated w/ a swap instrument.
**
** Only glitch is possible inclusion of a stub date, which is necessary if the
** maturity date isn't an integral number of frequency intervals away, e.g.
** a swap date 5 years and 1 month from the value date, which would have a
** stub date 1 month from now, followed by coupons every year from then.
**
** Returns date list associated with the swap's coupons, NULL if error found.
***************************************************************************
*/
TDateList* JpmcdsZCGetSwapCouponDL
  (TDate           valueDate,           /* (I) Value date */
   TDate           matDate,             /* (I) Unadjusted maturity date */
   TBoolean        stubAtEnd,           /* (I) If matDate onCycle from vdate */
   TDateInterval  *interval,            /* (I) Payment interval */
   TBadDayList    *badDayList,          /* (I) Bad day adjustment list */
   long            badDayConv,          /* (I) */
   char           *holidayFile);        /* (I) */


/*f
***************************************************************************
** Makes a cash flow list for a swap instrument.
**
** return s cash flow list associated with the swap, NULL if an error found.
***************************************************************************
*/
TCashFlowList* JpmcdsZCGetSwapCFL
  (TDate           valueDate,           /* (I) Value date                  */
   TDate           matDate,             /* (I) Unadjusted maturity date    */
   TBoolean        stubAtEnd,           /* (I) If matDate on-cycle         */
   double          rate,                /* (I) Coupon rate                 */
   TDateInterval  *interval,            /* (I) Coupon payment interval     */
   long            dayCountConv,        /* (I) Day count convention        */
   TBadDayList    *badDayList,          /* (I) Bad day adjustment list     */
   long            badDayConv,          /* (I) Bad day convention for rate */
   char           *holidayFile);        /* (I) Name of holiday file        */


/*f
***************************************************************************
** Adds simple-interest money market bond to ZCurve
**
** ZCurve is updated to include information from money market instruments
***************************************************************************
*/
int JpmcdsZCAddMoneyMarket(
   ZCurve                 *zc,             /* (M) ZCurve to add to */
   TDate                   dates[],        /* (I) maturity dates of instrs */
   double                  rates[],        /* (I) rates (e.g. 0.06 for 6%) */
   int                     n,              /* (I) number of instruments */
   long                    dayCountConv);  /* (I) day-count-convention */


/*f
***************************************************************************
** Adds points into a zero curve to match all points in a cashflow list 
** after a certain date. Note that only the *dates* in the cashflow list
** are used.
**
** Zero curve is modified.
***************************************************************************
*/
int JpmcdsZCAddCFLPoints(
   ZCurve                 *zc,          /* (M) ZCurve to add info to */
   TCashFlowList          *cfl,         /* (I) Use dates only for pts to add */
   TDate                   lastDate,    /* (I) Only add pts after this date*/
   long                    interpType,  /* (I) Interpolation method's type */
   TInterpData            *interpData); /* (I) Interpolation method's data */


/*f
***************************************************************************
** Adds information represented by a list-of-cash-flows to a zero curve.  
** Any cash flows which are already covered by the zero curve are 
** discounted at rates derived from the zero curve.  Cash flows beyond the
** zero curve imply discount factors, which are added to the zero curve.  If
** there is more than one such cash flow, several points are added to the curve,
** which are calculated by using an interative root-finding secant method,
** where the discount factor for the last cash flow is guessed (and the other
** discount factors are implied by interpolation) where the current price = 
** net present value of all the cash flows.
**
** The zero curve is updated to reflect the cash flows.  A point is added for 
** every cash flow, if not already in ZCurve list.  For linear forwards, all 
** interpolated points are returned, e.g. for 1 month forwards in an annual 
** market, monthly points will be returned, not just yearly. 
**
** Notes: date may be set for non-linear-forward interpolation methods to a 
** date to be added to the zero curve.  This allows production of a curve 
** with "nice" dates.
***************************************************************************
*/
int JpmcdsZCAddCashFlowList(
   ZCurve                 *zc,          /* (M) ZCurve to add info to */
   TCashFlowList          *cfl,         /* (I) Cash flows to add to ZCurve */
   double                  price,       /* (I) Current price of cash-flows */
   TDate                   date,        /* (I) Date to add last point at */
   long                    interpType,  /* (I) Interpolation method's type */
   TInterpData            *interpData); /* (I) Interpolation method's data */


/*f
***************************************************************************
** Adds a general zero rate to a ZCurve.
**
** ZCurve is updated to include the rate.
***************************************************************************
*/
int JpmcdsZCAddGenRate(
   ZCurve         *zc,                     /* (M) ZCurve to add rate to */
   TDate           date,                   /* (I) date associated with rate */
   double          rate,                   /* (I) rate, in ZCurve terms */
   long            basis,                  /* (I) compounding basis of rate*/
   long            dayCountConv);          /* (I) day counting convention */


/*f
***************************************************************************
** Adds a discount factor (at a specified date) to a ZCurve.
**
** Modified ZCurve, and error indicator.
***************************************************************************
*/
int JpmcdsZCAddDiscountFactor(
   ZCurve                  *zc,            /* (M) ZCurve to add discount to */
   TDate                    date,          /* (I) date assoc'd w/ disc factor */
   double                   disc);         /* (I) discount factor to add */


/*f
***************************************************************************
** Adds a zero rate to a ZCurve.
** The input rate must be in the proper ZCurve terms (zc->basis and zc->dayCountConv).
**
** ZCurve is updated to include the rate.
***************************************************************************
*/
int JpmcdsZCAddRate(
   ZCurve         *zc,                     /* (M) ZCurve to add rate to */
   TDate           date,                   /* (I) date associated with rate */
   double          rate);                  /* (I) rate, in ZCurve terms */


/*f
***************************************************************************
** Adds a zero rate and discount to a ZCurve.  The input rate must 
** The input rate must be in the proper ZCurve terms (zc->basis and zc->dayCountConv).
**
** ZCurve is updated to include the rate.
***************************************************************************
*/
int JpmcdsZCAddRateAndDiscount(
   ZCurve         *zc,                     /* (M) ZCurve to add rate to */
   TDate           date,                   /* (I) date associated with rate */
   double          rate,                   /* (I) rate, in ZCurve terms */
   double          disc);                  /* (I) discount factor to add */


/*f
***************************************************************************
** Calculates a simple interest money market rate from a zero curve.
***************************************************************************
*/
int JpmcdsZCMoneyMarketRate(
   ZCurve                 *zc,             /* (I) ZCurve */
   TDate                   date,           /* (I) date to start forward */
   long                    dayCountConv,   /* (I) day counting convention */
   long                    interpType,     /* (I) interpolation type to use */
   TInterpData            *interpData,     /* (I) interpolation method to use */
   double                 *rate);          /* (O) money market rate */


/*f
***************************************************************************
** Calculates a forward rate from a zero curve.
***************************************************************************
*/
int JpmcdsZCForwardRate(
   ZCurve                 *zc,             /* (I) ZCurve */
   TDate                   date1,          /* (I) date to start forward */
   TDate                   date2,          /* (I) date to end forward */
   long                    basis,          /* (I) compounding basis */
   long                    dayCountConv,   /* (I) day counting convention */
   long                    interpType,     /* (I) interpolation type to use */
   TInterpData            *interpData,     /* (I) interpolation method to use */
   double                 *rate);          /* (O) forward rate */


/*
***************************************************************************
** Calculates a discount factor for a forward between two dates.
***************************************************************************
*/
int JpmcdsZCFwdDiscountFactor(
   ZCurve                 *zc,             /* (I) ZCurve */
   TDate                   date1,          /* (I) date to start forward */
   TDate                   date2,          /* (I) date to end forward */
   long                    interpType,     /* (I) interpolation type to use */
   TInterpData            *interpData,     /* (I) interpolation method to use */
   double                 *discount);      /* (O) forward discount factor */


/*f
***************************************************************************
** Calculates net-present-value of a cash-flow list.
***************************************************************************
*/
int JpmcdsJpmcdsZCPresentValueCFL(
   ZCurve                 *zc,             /* (I) ZCurve for discounting */
   TCashFlowList          *cfl,            /* (I) cash flow list to value */
   int                     iLo,            /* (I) index of 1st c.f. to value */
   int                     iHi,            /* (I) index of last c.f. to value */
   long                    interpType,     /* (I) type of interpolation */
   TInterpData            *interpData,     /* (I) data for interpolation */
   double                 *pv);            /* (O) present value of cash flows */


/*f
***************************************************************************
** Calculates net-present-value of a cash-flow (a payment at a given date 
** in the future).
***************************************************************************
*/
int JpmcdsZCPresentValue(
   ZCurve                     *zc,         /* (I) ZCurve used for discounting */
   double                      price,      /* (I) amount of payment */
   TDate                       date,       /* (I) date of payment */
   long                        interpType, /* (I) type of interpolation */
   TInterpData                *interpData, /* (I) data for interpolation */
   double                     *pv);        /* (I) present value of price */


/*f
***************************************************************************
** Calculates discount factor for a date.
***************************************************************************
*/
int JpmcdsZCDiscountFactor(
   ZCurve                 *zc,             /* (I) ZCurve */
   TDate                   date,           /* (I) date to calculate d.f. for */
   long                    interpType,     /* (I) interpolation type to use */
   TInterpData            *interpData,     /* (I) interpolation method to use */
   double                 *discount);      /* (O) discount factor */


/*f
***************************************************************************
** Converts zc-style rate into a discount factor.
**
** Note: As a special optimization, the normal ZCurve case is hard coded for
** speed (since this routine is critical to performance of building zero-curves).
** That is, annually-compounded rates w/ 360 or 365-fixed year are calculated
** directly.
***************************************************************************
*/
int JpmcdsZCComputeDiscount(
   ZCurve         *zc,                     /* (I) zero curve */
   TDate           date,                   /* (I) date assoc'd w/ rate */
   double          rate,                   /* (I) rate in ZC basis & dayConv */
   double         *discount);              /* (O) discount factor */


/*f
***************************************************************************
** Calculates an interpolated rate for a date.
**
** Note: Piece wise interpolation allows different areas of the zero curve
** to be interpolated using different methods.  Basically an array of 
** <date,interpolationStuff> is given, where the interpolationStuff is used
** for any date before the given date.  The code allows the interpolationStuff
** to be another pieceWise interpolation type, although the utility of this
** is unknown.
***************************************************************************
*/
int JpmcdsZCInterpolate(
   ZCurve               *zc,               /* (I) ZCurve to interpolate on */
   TDate                 date,             /* (I) date to interpolate to */
   long                  interpType,       /* (I) type of interpolation's use */
   TInterpData          *interpData,       /* (I) data for interp method */
   double               *rate);            /* (O) rate to return */


/*f
***************************************************************************
** Finds an exact date in a ZC.
** Returns index of matching date if found, -1 if not.
***************************************************************************
*/
int JpmcdsZCFindDateExact(
   ZCurve         *zc,                     /* (I) ZCurve to search */
   TDate           date);                  /* (I) date to search for */


/*
***************************************************************************
** Finds index of entry closest to an input date.
** Returns index closest to date.
***************************************************************************
*/
int JpmcdsZCFindDateClosest(
   ZCurve         *zc,                     /* (I) ZCurve to search */
   TDate           date);                  /* (I) date to search for */


/*f
***************************************************************************
** Models the swap points which are to be added to the zero curve using the 
** more detailed swap structures.
***************************************************************************
*/
int JpmcdsZCValueFixFltSwap(
 ZCurve           *zc,                 /* (M) ZCurve to add to           */
 TCurve           *discZC,             /* (I) Discount zero curve        */
 double            price,              /* (I) Par price: usually 1.0     */
 TDate             matDate,            /* (I) Maturity date              */
 double            rate,               /* (I) Coupon rate (0.06 for 6%)  */
 int               fixedSwapFreq,      /* (I) Fixed coupon frequency     */
 int               floatSwapFreq,      /* (I) Floating coupon frequency  */
 long              fixDayCountConv,    /* (I) Fixed leg convention       */
 long              floatDayCountConv,  /* (I) Floating leg convention    */
 long              interpType,         /* (I) Interpolation method       */
 TInterpData      *interpData,         /* (I) Interpolation data         */
 TBadDayAndStubPos badDayAndStubPos,   /* (I) See JpmcdsBusinessDay      */
 char             *holidayFile);       /* (I) See JpmcdsBusinessDay      */


#ifdef __cplusplus
}
#endif

#endif    /* ZCPRVT_H */
