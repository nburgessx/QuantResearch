/*
 * ISDA CDS Standard Model
 *
 * Copyright (C) 2009 International Swaps and Derivatives Association, Inc.
 * Developed and supported in collaboration with Markit
 * 
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the ISDA CDS Standard Model Public License.
 */

#include "tcurve.h" 
#include "cheaders.h"
#include "bastypes.h"
#include "datelist.h"
#include "dtlist.h"
#include "convert.h"
#include "interp.h"
#include "lintrp.h"
#include "yearfrac.h"
#include "ldate.h"
#include "macros.h"
#include <math.h>


/*
***************************************************************************
** Make sure the TCurve is valid.
**
** This routine guarantees we can generate positive discount factors.
***************************************************************************
*/
static int checkTCurve(
    TCurve *tc,                /* (I) TCurve to check */
    char   *name,              /* (I) Name of TCurve for error msgs */
    char   *routine);          /* (I) For error messages. */


/*
***************************************************************************
** Allocates and sets a TCurve.
**
** Checks parameters for validity. Returns pointer to allocated TCurve. 
** Only allocates memory for the array when it is returning a valid TCurve.
***************************************************************************
*/
TCurve* JpmcdsNewTCurve
    (TDate  baseDate,           /* (I) Base date */
     int    numPts,             /* (I) Length of dates & rates */
     double basis,              /* (I) Compounding periods/year */
     long   dayCountConv)       /* (I) Day count convention */
{
    static char routine[] = "JpmcdsNewTCurve";
    int         status    = FAILURE;
    TCurve     *tc        = NULL;

    /* Check numPts
    */
    if (numPts < 0)
    {
        JpmcdsErrMsg ("%s: Attempt to create curve with %d points.\n",
            routine, (int)numPts);
        goto done;
    }


    /* Check basis
    */
    if (basis < -2.0)
    {
        JpmcdsErrMsg ("%s: Attempt to create curve with basis %f.\n",
                routine, basis);
        goto done;
    }

    /* Allocate the TCurve and set parameters.
    */
    tc = NEW(TCurve);
    if (tc == NULL)
        goto done;

    /* Set params. */
    tc->fDayCountConv = dayCountConv;
    tc->fBasis = basis;
    tc->fBaseDate = baseDate;
    tc->fNumItems = numPts;
    if (numPts > 0)
    {
        tc->fArray = NEW_ARRAY(TRatePt, numPts);
        if( tc->fArray == NULL )
        {
            goto done;
        }
    }
    else
    {
        tc->fArray = NULL;
    }

    status = SUCCESS; 
 
done: 

    if ( status == FAILURE )
    {
        JpmcdsFreeTCurve(tc);
        JpmcdsErrMsg("%s: Failed.\n", routine );
        return NULL;
    }

    return tc;
}


/*
***************************************************************************
** Creates a new TCurve, and fills it in with rates, dates, yearfraction 
** function, and basis.
** 
** If it encounters a TDate = 0, it assumes that this is a terminator. 
** This allows spreadsheet users to select a big range of dates, but use 
** less than the whole range by terminating with a 1 or 0.
***************************************************************************
*/
EXPORT TCurve* JpmcdsMakeTCurve(
    TDate   baseDate,                    /* (I) Discount date  */
    TDate  *dates,                       /* (I) Dates in curve */
    double *rates,                       /* (I) Associated rates */
    int     numPts,                      /* (I) Len of dates and rates arrays */
    double  basis,                       /* (I) Compounding periods/year */
    long    dayCountConv)                /* (I) */
{
    static char routine[]= "JpmcdsMakeTCurve";
    int         status   = FAILURE;      /* Until successful */
    TCurve     *curve    = NULL;

    /* Set up curve struct from Dates and Rates. */
    curve = JpmcdsMakeTCurveNoRateCheck(baseDate,
                                     dates,
                                     rates,
                                     numPts,
                                     basis,
                                     dayCountConv);
    if (curve == NULL)
        goto done;

    /* Make sure the TCurve is valid. This routine will guarantee we can
     * generate positive discount factors.
     */
    if (checkTCurve(curve, "TCurve", routine) != SUCCESS)
        goto done;

    status = SUCCESS;

 done:
    if ( status == FAILURE )
    {
        JpmcdsFreeTCurve(curve);
        curve = NULL;
        JpmcdsErrMsg("%s: Failed.\n", routine );
    }

    return curve;
}


/*
***************************************************************************
** Creates a new TCurve, and fills it in with rates, dates, yearfraction 
** function, and basis.
** 
** If it encounters a TDate = 0, it assumes that this is a terminator. 
** This allows spreadsheet users to select a big range of dates, but use 
** less than the whole range by terminating with a 1 or 0.
**
** It is different from JpmcdsMakeTCurve in that allows any values for the 
** rates. JpmcdsMakeTCurve checks that the rates are all greater than -1. 
***************************************************************************
*/
TCurve* JpmcdsMakeTCurveNoRateCheck(
    TDate   baseDate,                   /* (I) Discount date  */
    TDate  *dates,                      /* (I) Dates in curve */
    double *rates,                      /* (I) Associated rates */
    int     numPts,                     /* (I) Len of dates and rates arrays */
    double  basis,                      /* (I) Compounding periods/year */
    long    dayCountConv)               /* (I) */
{
    static char routine[]= "JpmcdsMakeTCurveNoRateCheck";
    int         status   = FAILURE;
    int         idx;
    TCurve     *curve    = NULL;

    /* Set up curve struct from Dates and Rates. */
    curve = JpmcdsNewTCurve(baseDate, numPts, basis,
                         dayCountConv);

    if (curve == NULL)
        goto done;

    for (idx = 0; idx < numPts; idx++)
    {
        if (dates[idx] == 0)            /* Check for NULL termination */
        {
            curve->fNumItems = idx;     /* Extra locations simply not used */
            break;
        }
        
        /* Make sure dates are in increasing order. Otherwise
         * interpolation doesn't work.
         */
        if (idx > 0 && dates[idx-1] >= dates[idx])
        {
            JpmcdsErrMsg("%s: Date %d (%s) >= date %d (%s) in TCurve.\n",
                      routine,
                      idx, JpmcdsFormatDate(dates[idx-1]),
                      idx+1, JpmcdsFormatDate(dates[idx]));
            goto done;
        }
        else
        {
            curve->fArray[idx].fDate = dates[idx];
            curve->fArray[idx].fRate = rates[idx];
        }
    }

    status = SUCCESS;

done:
    if ( status == FAILURE )
    {
        JpmcdsFreeTCurve(curve);
        curve = NULL;
        JpmcdsErrMsg("%s: Failed.\n", routine );
    }

    return curve;
}


/*
***************************************************************************
** Creates an array of dates from the input curve and returns a pointer to 
** the array on success and NULL on failure.
***************************************************************************
*/
TDate* JpmcdsDatesFromCurve
    (TCurve   *curve)    /* (I) TCurve */
{
    static char routine[] = "JpmcdsDatesFromCurve";

    int       i;
    TDate   *dates = NULL;

    if (curve->fNumItems == 0)
        return NULL; /* no dates in curve is not an error */

    dates = NEW_ARRAY(TDate, curve->fNumItems);

    if (dates==NULL)
       goto done;

    for (i=0; i<curve->fNumItems; i++)
        dates[i] = curve->fArray[i].fDate;

 done:
    if (dates==NULL)
    {
       JpmcdsErrMsg("%s Failed.\n", routine);
    }

    return dates;
}


/*
***************************************************************************
** Calculate present value of a cash flow list.
***************************************************************************
*/
int JpmcdsCashFlowPV
    (TCashFlowList *cfl,                /* (I) Cash Flow List */
     TCurve        *zeroCurve,          /* (I) Zero Curve */
     long           interpType,         /* (I) Interp method */
     double        *presValue)          /* (O) Present Value */
{
    int status = FAILURE;
    static char routine[]="JpmcdsCashFlowPV";
    int idx;
    double discountFactor;
    double pv = 0;

    for (idx=0;idx < cfl->fNumItems; idx++)
    {
        if (JpmcdsDiscountDate(cfl->fArray[idx].fDate,
                            zeroCurve,
                            interpType,
                            &discountFactor) == FAILURE)
            goto done;

        pv += (discountFactor * cfl->fArray[idx].fAmount);
    }
    
    status = SUCCESS;

done:
    if (status == FAILURE)
        JpmcdsErrMsg("%s: Failed.\n", routine);
    else
        *presValue = pv;
        
    return status;
}


/*
***************************************************************************
** Calculate discount factor for specified date and zero curve.
***************************************************************************
*/
int JpmcdsDiscountDate
     (TDate   date,             /* (I) Date */
      TCurve *zc,               /* (I) Zero Curve */
      long    interpType,       /* (I) See interp.h */
      double *discount)         /* (O) Present Value */
{
    static char routine[]="JpmcdsDiscountDate";
    double interpRate;
    double yearsToDate;
    int status = FAILURE;

    if (zc == NULL || discount == NULL)
    {
        JpmcdsErrMsg ("%s: NULL inputs\n", routine);
        goto done;
    }

    /* Interpolate  to get the zero-coupon rate. */
    if (JpmcdsInterpRate (date, zc, interpType, &interpRate) == FAILURE)
        goto done;


    if (JpmcdsDayCountFraction(zc->fBaseDate, date, zc->fDayCountConv, &yearsToDate) == FAILURE)
        goto done;

    status  = JpmcdsRateToDiscountYearFrac(interpRate,yearsToDate,zc->fBasis, discount);
    
done:
    if (status == FAILURE)
        JpmcdsErrMsg("%s: Failed.\n", routine);

    return status;
}


/*
***************************************************************************
** Frees a TCurve.
***************************************************************************
*/
void JpmcdsFreeTCurve(TCurve *theCurve)
{
    if (theCurve != (TCurve *)NULL)
    {
        FREE(theCurve->fArray);
        FREE(theCurve);
    }
}


/*
***************************************************************************
** Creates a copy of the curve.
***************************************************************************
*/
TCurve* JpmcdsCopyCurve(TCurve *aCurve) /* (I) */
{
    static char routine[]="JpmcdsCopyCurve";
    TCurve * newCurve= NULL;
    int index;

    if (aCurve == NULL)
    {
        JpmcdsErrMsg("%s: Cannot copy NULL TCurve.\n", routine);
        goto done;
    }

    newCurve = JpmcdsNewTCurve(aCurve->fBaseDate,
                            aCurve->fNumItems,
                            aCurve->fBasis,
                            aCurve->fDayCountConv);
    if (newCurve == NULL)
        goto done;
    
    for (index = 0; index < aCurve->fNumItems; index++)
    {
        newCurve->fArray[index] = aCurve->fArray[index];
    }

done:
    if (newCurve == NULL)
        JpmcdsErrMsg("%s: Failed.\n", routine);

    return newCurve;
}


/*
***************************************************************************
** Computes the first valid date for a zero curve.
**
** For an encapsulated curve, this is a method of the curve. Otherwise
** this is the minimum of the base date of the curve and the first date
** in the date array.
**
** The date returned is the first date for which we have genuine
** information regarding rates and discount factors. It is possible that
** using the curve for dates before this date will give answers, but they
** might not be based on real information.
***************************************************************************
*/
int JpmcdsCurveFirstDate
(TCurve *curve,    /* (I) */
 TDate  *firstDate /* (O) */
)
{
    static char routine[] = "JpmcdsCurveFirstDate";

    if (curve == NULL || firstDate == NULL)
    {
        JpmcdsErrMsg ("%s: NULL inputs.\n", routine);
        return FAILURE;
    }

    if (curve->fNumItems == 0)
    {
        *firstDate = curve->fBaseDate;
    }
    else
    {
        *firstDate = MIN (curve->fBaseDate, curve->fArray[0].fDate);
    }
    return SUCCESS;
}


/*
***************************************************************************
** Make sure the TCurve is valid.
**
** This routine guarantees we can generate positive discount factors.
***************************************************************************
*/
static int checkTCurve(
    TCurve *tc,                /* (I) TCurve to check */
    char   *name,              /* (I) Name of TCurve for error msgs */
    char   *routine)           /* (I) For error messages. */
{
    int idx;
    int status = SUCCESS;       /* Until proven a failure */

    if (tc == NULL)
    {
       JpmcdsErrMsg("%s: Null %s received.\n", routine, name);
       return FAILURE;
    }

    /* Check number of points.
     */
    if (tc->fNumItems < 0)
    {
        JpmcdsErrMsg("%s: %s with %d points received.\n", routine, name, tc->fNumItems);
        status = FAILURE;
    }

    /* Check basis
    */
    if (tc->fBasis < -2.0)
    {
        JpmcdsErrMsg ("%s: Attempt to create %s with basis %f.\n", routine, name, tc->fBasis);
        status = FAILURE;
    }

    /* Check day count convention.
     */
    if (JpmcdsDayCountValid(routine, tc->fDayCountConv) == FAILURE)
    {
        JpmcdsErrMsg("%s: Bad day count convention in %s.\n", routine, name);
        status = FAILURE;
    }

    /* Check dates and rates for each point.
     */
    for (idx = 0; idx < tc->fNumItems; idx++)
    {
        /* Make sure dates are in increasing order. 
         * Otherwise interpolation doesn't work.
         */
        if (idx > 0 && 
            tc->fArray[idx-1].fDate >= tc->fArray[idx].fDate)
        {
            JpmcdsErrMsg("%s: Date %d (%s) >= date %d (%s) in %s.\n",
                      routine,
                      idx, JpmcdsFormatDate(tc->fArray[idx-1].fDate),
                      idx+1, JpmcdsFormatDate(tc->fArray[idx].fDate), name);
            status = FAILURE;
        }

        /* Make sure rates in TCurve are valid (i.e. produce positive
         * discount factors.
         */       
        if (JpmcdsRateValid(routine,
                         tc->fArray[idx].fRate,
                         tc->fBaseDate,
                         tc->fArray[idx].fDate,
                         tc->fDayCountConv,
                         tc->fBasis) != SUCCESS)
        {
            JpmcdsErrMsg("%s: Rate %d (%f) at %s in %s implies "
                      "non-positive discount factor.\n",
                      routine, idx+1, tc->fArray[idx].fRate,
                      JpmcdsFormatDate(tc->fArray[idx].fDate), 
                      name);
            status = FAILURE;
        }
    }

    return status;
}

