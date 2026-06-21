/*
 * ISDA CDS Standard Model
 *
 * Copyright (C) 2009 International Swaps and Derivatives Association, Inc.
 * Developed and supported in collaboration with Markit
 * 
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the ISDA CDS Standard Model Public License.
 */

#include <ctype.h>
#include <math.h>
#include "bastypes.h"
#include "cerror.h"
#include "cgeneral.h"
#include "interp.h"
#include "ldate.h"
#include "lintrp.h" 
#include "macros.h" 
#include "convert.h"
#include "cfinanci.h"
#include "bsearch.h"
#include "strutil.h"
#include "tcurve.h"    
#include "zr2fwd.h"                                              
#include "yearfrac.h"


/*
***************************************************************************
** Interpolates rates using Flat forwards.
**
** Note: if extrapolation is necessary before the curve, then no calculation 
** is done, simply the beginning point is returned.  Extrapolation does occur
** after the curve though...
***************************************************************************
*/
static int InterpRateFlatFwds(
   double  desiredDate,                /* (I) Date to interpolate at */
   TCurve *curve,                      /* (I) Curve with rates */
   double *rate);                      /* (O) Interpolated rate */


/*
***************************************************************************
** Linearly Interpolates for a TCurve which actually stores discount factors.
**
** Note: if extrapolation is necessary before the curve, then no calculation 
** is done, simply the beginning point is returned.  extrapolation does occur
** after the curve though...
***************************************************************************
*/
static int InterpRateDiscountFactors( /* Used for basis=JPMCDS_DISCOUNT_FACTOR */
   double  desiredDate,                /* (I) Date to interpolate at */
   TCurve *curve,                      /* (I) Curve with rates */
   double *rate);                      /* (O) Interpolated rate */

/*f
***************************************************************************
** Finds a date in a curve for interpolation.
***************************************************************************
*/
static int InterpRateFind(     
   double  date,                       /* (I) Date to interpolate at */
   TCurve *tc,                         /* (I) Curve with rates */
   int    *lo,                         /* (O) lower bracketing point */
   int    *hi);                        /* (O) upper bracketing point */


/*
***************************************************************************
** Interpolates a TCurve at a TDate using the specified interp method. 
***************************************************************************
*/
int JpmcdsInterpRate(
   TDate   desiredDate,                /* (I) Date to interpolate at */
   TCurve *curve,                      /* (I) Curve with rates */
   long    interpType,                 /* (I) Type of interpolation to use */
   double *rate)                       /* (O) Interpolated rate */
{
    static char routine[]="JpmcdsInterpRate";
    int         status = FAILURE;

    if (JpmcdsInterpRateExact((double)desiredDate, curve, interpType, rate) == FAILURE)
        goto done;

    status = SUCCESS;

done:
    if (status != SUCCESS)
        JpmcdsErrMsgFailure(routine);

    return status;
}


/*
***************************************************************************
** Using curve dates & rates, routine interpolates to find a new rate at 
** the desired date. Only difference between this and JpmcdsInterpRate is that 
** this routine's version of desired date is a double, enabling it to handle 
** continuous time, which is sometimes used in trees.
***************************************************************************
*/
int   JpmcdsInterpRateExact( 
   double  desiredDate,                /* (I) Date to interpolate at */
   TCurve *curve,                      /* (I) Curve with rates */
   long    interpType,                 /* (I) Type of interpolation to use */
   double *rate)                       /* (O) Interpolated rate */
{
    static char routine[]="JpmcdsInterpRateExact";
    int status = FAILURE;              /* Until proved a success */
    int order = 0;
    int curveLen = curve->fNumItems;

    /*
    ** Since we are using the general interpolation functions, but we
    ** know that the TCurve should have been constructed properly,
    ** we will turn off order checking in the binary search routines.
    **
    ** This might well produce a dramatic improvement in all the ALIB
    ** zero curve analytics!
    */
    TBoolean checkOrder = JpmcdsBSearchCheckOrder (FALSE);

    /* Check for at least one point in TCurve.
     */
    if (curveLen < 1)
    {
        JpmcdsErrMsg("%s: Number of points (%d) in TCurve < 1.\n", routine, curveLen);
        goto done;
    }


    /* Do not extrapolate, regardless of what the interpolation type is.
     * Just use the first or last point.  Naturally, discount factors
     * can't be flat. 
     */
    if (curve->fBasis != JPMCDS_DISCOUNT_FACTOR)
    {
       if (desiredDate > curve->fArray[curveLen-1].fDate)
       {
           *rate = curve->fArray[curveLen-1].fRate;
           return SUCCESS;
       }
       else if (desiredDate < curve->fArray[0].fDate)
       {
           *rate = curve->fArray[0].fRate;
           return SUCCESS;
       }
    }

    switch (interpType)
    {
    case JPMCDS_FLAT_FORWARDS:
        if (InterpRateFlatFwds(desiredDate,curve,rate) == FAILURE)
            goto done;
        break;
        
    case JPMCDS_LINEAR_INTERP:
        if (curve->fBasis==JPMCDS_DISCOUNT_FACTOR)    /* special handling */
        {
            if (InterpRateDiscountFactors(desiredDate,curve,rate) == FAILURE)
                goto done;
            break;
        }
        if (JpmcdsLinInterpLongPoint1(&curve->fArray[0].fDate,
                                   sizeof(TRatePt), curve->fNumItems,
                                   &curve->fArray[0].fRate,
                                   sizeof(TRatePt),
                                   desiredDate,
                                   (TMetricDoubleFunc)NULL,
                                   rate) == FAILURE)
            goto done;
        break;
        
    default:
        JpmcdsErrMsg("%s: Unknown interpolation type %ld.\n", routine, interpType);
        goto done;
    }
    
    status = SUCCESS;


done:
    /* Re-set the global check order flag. */
    JpmcdsBSearchCheckOrder (checkOrder);

    if (status == FAILURE)
        JpmcdsErrMsg("%s: Failed.\n", routine);

    return status;
}


/*
***************************************************************************
** Interpolates rates using Flat forwards.
**
** Note: if extrapolation is necessary before the curve, then no calculation 
** is done, simply the beginning point is returned.  Extrapolation does occur
** after the curve though...
***************************************************************************
*/
static int InterpRateFlatFwds(
   double  desiredDate,                /* (I) Date to interpolate at */
   TCurve *curve,                      /* (I) Curve with rates */
   double *rate)                       /* (O) Interpolated rate */
{
    int         status = FAILURE;
    static char routine[] = "InterpRateFlatFwds";
    double      discLo;                 /* discount factor for lo */
    double      discHi;                 /* discount factor for hi */
    double      discDate;               /* discount factor for date */
    int         hi;                     /* index of known point above desired*/
    int         lo;                     /* index of known point below desired*/
    TDate       dateHi;                 /* date at hi */
    TDate       dateLo;                 /* date at lo */
    double      rateHi;                 /* rate at hi */
    double      rateLo;                 /* rate at lo */
    long        diffh;
    long        diffl;
    double      doubleDiffh;
    double      doubleDiffl;
    TBoolean    isFractionalDate;       /* is desiredDate not an integer */
    TDate       desiredTDate;           /* cast of desiredDate to TDate */

    desiredTDate = (TDate) desiredDate;
    isFractionalDate = ! ARE_ALMOST_EQUAL(desiredDate, (double)desiredTDate);

    if (InterpRateFind(desiredDate,curve,&lo,&hi)==FAILURE)
        goto done;

    dateHi = curve->fArray[hi].fDate;
    dateLo = curve->fArray[lo].fDate;
    rateHi = curve->fArray[hi].fRate;
    rateLo = curve->fArray[lo].fRate;

    if (dateLo == desiredDate ||                 /* no need to interpolate */
        (dateLo > desiredDate && curve->fBasis!=JPMCDS_DISCOUNT_FACTOR))
    {
        *rate = rateLo;
        return SUCCESS;
    }

    if (dateHi == desiredDate)
    {
        *rate = rateHi;
        return SUCCESS;
    }
    /* special case: basedate between lo&hi*/
    if (dateLo<=curve->fBaseDate && curve->fBaseDate<=dateHi &&
        curve->fBasis!=JPMCDS_DISCOUNT_FACTOR)
    {
        if (desiredDate < curve->fBaseDate)
        {
            *rate = rateLo;             /* rates are flat between lo and hi */
        }                               /* iff base date is between them. */
        else                            /* This ensures that basedate will */
        {                               /* be associated with a "kink" in the*/
            *rate = rateHi;             /* forwards, and discount(base)== 1.0*/
        }
        return SUCCESS;
    }

    if (JpmcdsRateToDiscount(rateLo, 
                          curve->fBaseDate, 
                          dateLo,
                          curve->fDayCountConv, 
                          (long) curve->fBasis,
                          &discLo) == FAILURE)
        goto done;

    if (JpmcdsRateToDiscount(rateHi, 
                          curve->fBaseDate, 
                          dateHi,
                          curve->fDayCountConv, 
                          (long) curve->fBasis,
                          &discHi) == FAILURE)
        goto done;


    if (discLo == 0.0)
    {
        JpmcdsErrMsg("%s: Zero discount factor\n",routine);
        goto done;
    }


    switch (curve->fDayCountConv) 
    {
      case JPMCDS_B30_360:         /* Fractional dates not supported */
      case JPMCDS_B30E_360:
      case JPMCDS_ACT_365:
        
        if (isFractionalDate)
        {
            JpmcdsErrMsg ("%s: Unsupported curve daycount convention %s for "
                       "fractional date %f\n",
                       routine,
                       JpmcdsFormatDayCountConv (curve->fDayCountConv),
                       desiredDate);
            goto done;
        }
        if (JpmcdsDaysDiff(dateLo, desiredTDate, curve->fDayCountConv, &diffh) != SUCCESS)
            goto done;
        
        if (JpmcdsDaysDiff(dateLo, dateHi, curve->fDayCountConv, &diffl) != SUCCESS)
            goto done;

        doubleDiffh = (double)diffh;
        doubleDiffl = (double)diffl;
        break;
        
        
      case JPMCDS_ACT_360:         /* Fractional dates supported */
      case JPMCDS_ACT_365F:
        doubleDiffl = dateHi - dateLo;
        doubleDiffh = desiredDate - dateLo;
        break;
        
    default:
        JpmcdsErrMsg("%s: Unsupported curve daycount convention %s\n",routine,
                  JpmcdsFormatDayCountConv(curve->fDayCountConv));
        goto done;
    }
    
    if (IS_ALMOST_ZERO(doubleDiffl) || doubleDiffl < 0.0)
    {
        JpmcdsErrMsg("%s: Unexpected 0 difference between curve dates.\n", routine);
        goto done;
    }

/*
** This formula computes the discount to desiredDate using the
** assumption of a flat forward rate between dateLo and dateHi
*/
    discDate = discLo * pow(discHi/discLo, doubleDiffh/doubleDiffl);

/*
** Now convert the discount to desiredDate into a rate. This is complicated
** by the possibility that the desiredDate is not a date, but is in fact
** a double. A double only makes sense for ACT/365F and ACT/360.
*/
    if (isFractionalDate)
    {
        /* Calculate year fraction by hand and use DiscountToRateYearFrac */

        double desiredDateYearFrac;

        switch (curve->fDayCountConv)
        {
        case JPMCDS_ACT_365F:
            desiredDateYearFrac = (desiredDate - curve->fBaseDate) / 365.0;
            break;
        case JPMCDS_ACT_360:
            desiredDateYearFrac = (desiredDate - curve->fBaseDate) / 360.0;
            break;
        default:
            PROGRAM_BUG();
            JpmcdsErrMsg ("%s: Unexpected day count convention %s\n",
                       routine,
                       JpmcdsFormatDayCountConv(curve->fDayCountConv));
            goto done;
        }
        if (JpmcdsDiscountToRateYearFrac (discDate,
                                       desiredDateYearFrac,
                                       curve->fBasis,
                                       rate) != SUCCESS)
        {
            goto done;
        }
    }
    else
    {
        /* Not a fractional date, so use DiscountToRate */
        if (JpmcdsDiscountToRate(discDate,
                              curve->fBaseDate,
                              desiredTDate,
                              curve->fDayCountConv,
                              (long) curve->fBasis,
                              rate) == FAILURE)
            goto done;
    }

    status = SUCCESS;

done:
    if (status == FAILURE)
        JpmcdsErrMsg("%s: Failed.\n",routine);

    return status;
}


/*
***************************************************************************
** Linearly Interpolates for a TCurve which actually stores discount factors.
**
** Method:
**          Find points in curve which surround desired date.
**          If desired date is in curve, return it's associated value.
**          Otherwise, get discount factors from curve.
**          Convert them to annually-compounded, actual/365Fixed rates.
**          Linearly interpolate (or use constant extrapolation on RATES).
**          Convert resulting rate back into a discount factor.
**          Return that as the result.
**
** Note: if extrapolation is necessary before the curve, then no calculation 
** is done, simply the beginning point is returned.  extrapolation does occur
** after the curve though...
**
** Confusing bit: The TCurve actually stores discount factors. This routine
** has local variables discXXX which refer to the discount factors,
** The rateXXX variables actually store Annually Compounded, Act/365Fixed
** rates.  That is how the interpolation is actually done.  The result
** of this routine is called rate, but is actually a discount factor.
***************************************************************************
*/
static int InterpRateDiscountFactors(
   double  desiredDate,                /* (I) Date to interpolate at */
   TCurve *curve,                      /* (I) Curve with rates */
   double *rate)                       /* (O) Interpolated rate */
{
    int         status = FAILURE;
    static char routine[] = "InterpRateDiscountFactors";
    double      discLo;                 /* discount factor for lo */
    double      discHi;                 /* discount factor for hi */
    int         hi;                     /* index of known point above desired*/
    int         lo;                     /* index of known point below desired*/
    TDate       dateHi;                 /* date at hi */
    TDate       dateLo;                 /* date at lo */
    double      rateHi;                 /* rate at hi */
    double      rateLo;                 /* rate at lo */
    double      rateDate;               /* rate at desiredDate */
    TBoolean    isFractionalDate;       /* is desiredDate not an integer */
    TDate       desiredTDate;           /* cast of desiredDate to TDate */

    desiredTDate = (TDate) desiredDate;
    isFractionalDate = ! ARE_ALMOST_EQUAL(desiredDate, (double)desiredTDate);

    if (InterpRateFind(desiredDate,curve,&lo,&hi)==FAILURE)
        goto done;

    dateHi = curve->fArray[hi].fDate;
    dateLo = curve->fArray[lo].fDate;
    discHi = curve->fArray[hi].fRate;   /* NOTE: Rates in TCurve are stored as*/
    discLo = curve->fArray[lo].fRate;   /* discount FACTORS. */

    if (dateLo == desiredDate)          /* no need to interpolate */
    {
        *rate = discLo;
        return SUCCESS;
    }

    if (dateHi == desiredDate)
    {
        *rate = discHi;
        return SUCCESS;
    }

    if (JpmcdsDiscountToRate(discLo, 
                          curve->fBaseDate, 
                          dateLo,
                          JPMCDS_ACT_365F,
                          JPMCDS_ANNUAL_BASIS,
                          &rateLo) == FAILURE)  goto done;
    
    if (JpmcdsDiscountToRate(discHi, 
                          curve->fBaseDate, 
                          dateHi,
                          JPMCDS_ACT_365F,
                          JPMCDS_ANNUAL_BASIS,
                          &rateHi) == FAILURE)  goto done;

    if (dateLo > desiredDate)      /* extrapolate at a constant RATE */
    {
       rateDate = rateLo;
    } 
    else if (desiredDate > dateHi) /* ditto */
    {
       rateDate = rateHi;
    }
    else                           /* use linear interpolation: */
    {
       rateDate = ((rateHi-rateLo)/(dateHi-dateLo))*(desiredDate-dateLo)+rateLo;
    }
/*
** Now convert the rate to desiredDate into a discount. This is complicated
** by the possibility that the desiredDate is not a date, but is in fact
** a double. A double only makes sense for ACT/365F and ACT/360.
*/
    if (isFractionalDate)
    {
        /* Calculate year fraction by hand and use DiscountToRateYearFrac */

        double desiredDateYearFrac;

        switch (curve->fDayCountConv)
        {
        case JPMCDS_ACT_365F:
            desiredDateYearFrac = (desiredDate - curve->fBaseDate) / 365.0;
            break;
        case JPMCDS_ACT_360:
            desiredDateYearFrac = (desiredDate - curve->fBaseDate) / 360.0;
            break;
        default:
            PROGRAM_BUG();
            JpmcdsErrMsg ("%s: Unexpected day count convention %s\n",
                       routine,
                       JpmcdsFormatDayCountConv(curve->fDayCountConv));
            goto done;
        }
        if (JpmcdsRateToDiscountYearFrac (rateDate,
                                       desiredDateYearFrac,
                                       JPMCDS_ANNUAL_BASIS,
                                       rate) != SUCCESS)
        {
            goto done;
        }
    }
    else
    {
        /* Not a fractional date, so use RateToDiscount */
        if (JpmcdsRateToDiscount(rateDate,
                              curve->fBaseDate,
                              desiredTDate,
                              JPMCDS_ACT_365F,
                              JPMCDS_ANNUAL_BASIS,
                              rate) == FAILURE)
            goto done;
    }

    status = SUCCESS;

done:
    if (status == FAILURE)
        JpmcdsErrMsg("%s: Failed.\n",routine);

    return status;
}


/*
***************************************************************************
** Finds a date in a curve for interpolation.
***************************************************************************
*/
static int InterpRateFind(     
   double  date,                       /* (I) Date to interpolate at */
   TCurve *tc,                         /* (I) Curve with rates */
   int    *lo,                         /* (O) lower bracketing point */
   int    *hi)                         /* (O) upper bracketing point */
{
    int status = FAILURE;
    static char routine[] = "InterpRateFind";

   if (JpmcdsBSearchLong(date, 
                      &tc->fArray[0].fDate, 
                      sizeof(TRatePt), 
                      tc->fNumItems,
                      lo, 
                      hi) == FAILURE)
       goto done;
    status = SUCCESS;

  done:
    if (status == FAILURE)
        JpmcdsErrMsg("%s: Failed.\n", routine);

    return status;
}
