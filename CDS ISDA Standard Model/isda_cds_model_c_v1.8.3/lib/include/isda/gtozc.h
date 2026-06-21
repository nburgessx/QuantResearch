/*
 * ISDA CDS Standard Model
 *
 * Copyright (C) 2009 International Swaps and Derivatives Association, Inc.
 * Developed and supported in collaboration with Markit
 * 
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the ISDA CDS Standard Model Public License.
 */

#ifndef GTOZC_H
#define GTOZC_H

#include "bastypes.h"
#include "cheaders.h"
#include "ldate.h"
#include "tcurve.h"

#ifdef __cplusplus
extern "C"
{
#endif


/*f
***************************************************************************
** Adds cash points to a given zero curve.
***************************************************************************
*/
TCurve* JpmcdsZCCash(
   TCurve        *zeroCurve,         /* (I) Zero curve to add to            */
   TDate         *dates,             /* (I) Array of cash dates             */
   double        *rates,             /* (I) Array of cash rates (0.06 for 6%)*/
   int           numRates,           /* (I) Number of cash rates to add     */
   long          dayCountConv);      /* (I) See JpmcdsDayCountFraction         */


/*f
***************************************************************************
** Adds swap points to a given zero curve.
***************************************************************************
*/
TCurve* JpmcdsZCSwaps(
 TCurve          *zeroCurve,         /* (I) Zero curve to add swap points to*/
 TCurve          *discZC,            /* (I) Zero curve used for discounting */
 TDate           *dates,             /* (I) Unadjusted swap maturity dates  */
 double          *rates,             /* (I) Swap par fixed rates (0.06=6%)  */
 int             numSwaps,           /* (I) Len of dates,rates,prices       */
 int             fixedSwapFreq,      /* (I) Fixed leg coupon frequency      */
 int             floatSwapFreq,      /* (I) Floating leg coupon frequency   */
 long            fixDayCountConv,    /* (I) See JpmcdsDayCountConvention       */
 long            floatDayCountConv,  /* (I) See JpmcdsDayCountConvention       */
 char            fwdLength,          /* (I) For fwd smoothing, length of fwds*/
 TBadDayAndStubPos badDayAndStubPos, /* (I) Bad day conv + stub pos         */
 char            *holidayFile);      /* (I) See JpmcdsBusinessDay              */


#ifdef __cplusplus
}
#endif

#endif    /* GTOZC_H */
