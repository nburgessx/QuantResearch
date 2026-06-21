/*
 * ISDA CDS Standard Model
 *
 * Copyright (C) 2009 International Swaps and Derivatives Association, Inc.
 * Developed and supported in collaboration with Markit
 * 
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the ISDA CDS Standard Model Public License.
 */

#include <assert.h>
#include <ctype.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "cxbsearch.h"
#include "cxdatelist.h"
#include "macros.h"
#include "cxzerocurve.h"

#include "datelist.h"
#include "ldate.h"
#include "tcurve.h"
#include "zr2fwd.h"
#include "convert.h"
#include "cerror.h"


#define NaN sqrt(-1.0)

static int zcInterpRate (TCurve*, TDate, long, long, double*);
static int zcRateCC (TCurve*, int, double*);


/*
***************************************************************************
** Calculates the zero price for a given start date and maturity date.
** Returns NaN for errors.
***************************************************************************
*/
double JpmcdsForwardZeroPrice
(TCurve* zeroCurve,
 TDate   startDate,
 TDate   maturityDate)
{
    double startPrice    = JpmcdsZeroPrice(zeroCurve, startDate);
    double maturityPrice = JpmcdsZeroPrice(zeroCurve, maturityDate);
    return maturityPrice / startPrice;
}


/*
***************************************************************************
** Calculates the zero price for a given date. Returns NaN for errors.
***************************************************************************
*/
EXPORT double JpmcdsZeroPrice
(TCurve* zeroCurve,
 TDate   date)
{
    double      zeroPrice = 0.0;
    double      rate;
    double      time;

    rate = JpmcdsZeroRate (zeroCurve, date);

    /* 
    ** rate is continuously compounded calculated between valueDate of
    ** the zeroCurve and the required date
    */
    time = (date - zeroCurve->fBaseDate) / 365.0;
    zeroPrice = exp(-rate * time);
    return zeroPrice;
}


/*
***************************************************************************
** Calculates the zero rate for a given date using ACT/365F and continously
** compounded rates.
**
** We do not trust the ALIB to do this because it does not extrapolate
** beyond the last date of the curve.
***************************************************************************
*/
double JpmcdsZeroRate
(TCurve* zeroCurve,
 TDate   date)
{
    static char routine[] = "JpmcdsZeroRate";
    int         status    = FAILURE;

    long        exact;
    long        lo;
    long        hi;
    double      rate = 0.0;

    REQUIRE (zeroCurve != NULL);
    REQUIRE (zeroCurve->fNumItems > 0);
    REQUIRE (zeroCurve->fArray != NULL);

    if (JpmcdsBinarySearchLong (date,
                            &zeroCurve->fArray[0].fDate,
                            sizeof(TRatePt),
                            zeroCurve->fNumItems,
                            &exact,
                            &lo,
                            &hi) != SUCCESS) 
        goto done;

    if (exact >= 0)
    {
        /* date found in zeroDates */
        if (zcRateCC (zeroCurve, exact, &rate) != SUCCESS)
            goto done;
    }
    else if (lo < 0)
    {
        /* date before start of zeroDates */
        if (zcRateCC (zeroCurve, 0, &rate) != SUCCESS)
            goto done;
    }
    else if (hi >= zeroCurve->fNumItems)
    {
        /* date after end of zeroDates */
        if (zeroCurve->fNumItems == 1)
        {
            if (zcRateCC (zeroCurve, 0, &rate) != SUCCESS)
                goto done;
        }
        else
        {
            /* extrapolate using last flat segment of the curve */
            lo = zeroCurve->fNumItems-2;
            hi = zeroCurve->fNumItems-1;
            if (zcInterpRate (zeroCurve, date, lo, hi, &rate) != SUCCESS)
                goto done;
        }
    }
    else
    {
        /* date between start and end of zeroDates */
        if (zcInterpRate (zeroCurve, date, lo, hi, &rate) != SUCCESS)
            goto done;
    }

    status = SUCCESS;

 done:

    if (status != SUCCESS)
    {
        JpmcdsErrMsgFailure (routine);
        return NaN;
    }

    return rate;
}


/*
***************************************************************************
** Interpolates a rate segment of a zero curve expressed with continuously
** compounded rates using flat forwards.
**
** Always returns a continously compounded ACT/365F rate.
***************************************************************************
*/
static int zcInterpRate 
(TCurve* zc, TDate date, long lo, long hi, double *rate)
{
    static char routine[] = "zcInterpRate";
    int         status    = FAILURE;

    long   t1;
    long   t2;
    long   t;
    double zt;
    double z1t1;
    double z2t2;
    double z1;
    double z2;

    t1   = zc->fArray[lo].fDate - zc->fBaseDate;
    t2   = zc->fArray[hi].fDate - zc->fBaseDate;
    t    = date - zc->fBaseDate;

    assert (t > t1);
    assert (t2 > t1);

    if (zcRateCC (zc, lo, &z1) != SUCCESS)
        goto done;

    if (zcRateCC (zc, hi, &z2) != SUCCESS)
        goto done;

    /* rates are continuously compounded, i.e. exp(-rt) */
    /* flat forwards => (zt) should be linear in t */
    z1t1 = z1 * t1;
    z2t2 = z2 * t2;
    if (t == 0)
    {
        /* If the date equals the base date, then the zero rate
           is essentially undefined and irrelevant - so let us
           get the rate for the following day which is in thr
           right ballpark at any rate. */
        /* An exception to this approach is when t2 = 0 as well.
           In this case, rate = z2 */
        if (t2 == 0)
        {
            *rate = z2;
            goto success;
        }
        t = 1;
    }
    zt   = z1t1 + (z2t2 - z1t1) * (double)(t - t1) / (double)(t2 - t1);
    *rate = zt / t;

success:
    status = SUCCESS;

 done:

    if (status != SUCCESS)
        JpmcdsErrMsgFailure (routine);

    return status;
}


/*
***************************************************************************
** Gets a continuously compounded rate from a TCurve for element idx.
***************************************************************************
*/
static int zcRateCC 
(TCurve *tc,
 int     idx,
 double *ccRate)
{
    int status;

    status =  JpmcdsConvertCompoundRate(tc->fArray[idx].fRate,
                                     tc->fBasis,
                                     tc->fDayCountConv,
                                     JPMCDS_CONTINUOUS_BASIS,
                                     JPMCDS_ACT_365F,
                                     ccRate);
    return status;
}


/*
***************************************************************************
** Converts a compound rate from one frequency to another.
** Can also convert between ACT-style day count conventions.
***************************************************************************
*/
int JpmcdsConvertCompoundRate 
(double  inRate,
 double  inBasis,
 long    inDayCountConv,
 double  outBasis,
 long    outDayCountConv,
 double *outRate)
{
    static char routine[] = "JpmcdsConvertCompoundRate";
    int         status    = FAILURE;

    double      ccRate;

    /* this routine is a hotspot and was taking too long for the case where we
       do nothing */

    if (IS_EQUAL(inBasis,outBasis))
    {
        if (inDayCountConv == outDayCountConv)
        {
            *outRate = inRate;
        }
        else if (inDayCountConv == JPMCDS_ACT_365F && outDayCountConv == JPMCDS_ACT_360)
        {
            *outRate = inRate * 360.0/365.0;
        }
        else if (inDayCountConv == JPMCDS_ACT_360 && outDayCountConv == JPMCDS_ACT_365F)
        {
            *outRate = inRate * 365.0/360.0;
        }
        else
        {
            JpmcdsErrMsg ("%s: Can only convert between ACT/360 and ACT/365F day count "
                       "conventions\n", routine);
            goto done;
        }
    }
    else
    {
        double dayFactor;

        if (inDayCountConv == outDayCountConv)
        {
            dayFactor = 1.0;
        }
        else if (inDayCountConv == JPMCDS_ACT_365F && outDayCountConv == JPMCDS_ACT_360)
        {
            dayFactor = 360.0/365.0;
        }
        else if (inDayCountConv == JPMCDS_ACT_360 && outDayCountConv == JPMCDS_ACT_365F)
        {
            dayFactor = 365.0/360.0;
        }
        else
        {
            JpmcdsErrMsg ("%s: Can only convert between ACT/360 and ACT/365F day count "
                       "conventions\n", routine);
            goto done;
        }

        /* convert inRate to ccRate, then convert to outRate */
        if (IS_EQUAL(inBasis, JPMCDS_CONTINUOUS_BASIS))
        {
            ccRate = inRate * dayFactor;
        }
        else if (inBasis >= 1.0 && inBasis <= 365.0)
        {
            ccRate = dayFactor * inBasis * log (1.0 + inRate / inBasis);
        }
        else
        {
            JpmcdsErrMsg ("%s: Input basis %f is not a compounding frequency\n",
                     routine, inBasis);
            goto done;
        }

        if (IS_EQUAL(outBasis, JPMCDS_CONTINUOUS_BASIS))
        {
            *outRate = ccRate;
        }
        else if (outBasis >= 1.0 && outBasis <= 365.0)
        {
            *outRate = outBasis * (exp (ccRate/outBasis) - 1.0);
        }
        else
        {
            JpmcdsErrMsg ("%s: Output basis %f is not a compounding frequency\n",
                     routine, outBasis);
            goto done;
        }
    }
    status = SUCCESS;

 done:

    if (status != SUCCESS)
        JpmcdsErrMsgFailure (routine);

    return status;
}
