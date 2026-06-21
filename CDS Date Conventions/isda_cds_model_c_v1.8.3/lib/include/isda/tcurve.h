/*
 * ISDA CDS Standard Model
 *
 * Copyright (C) 2009 International Swaps and Derivatives Association, Inc.
 * Developed and supported in collaboration with Markit
 * 
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the ISDA CDS Standard Model Public License.
 */

#ifndef TCURVE_H
#define TCURVE_H

#include "bastypes.h"
#include "interp.h"

#ifdef __cplusplus
extern "C"
{
#endif


/*f
***************************************************************************
** Allocates and sets a TCurve.
**
** Checks parameters for validity. Returns pointer to allocated TCurve. 
** Only allocates memory for the array when it is returning a valid TCurve.
***************************************************************************
*/
TCurve* JpmcdsNewTCurve(
    TDate  baseDate,            /* (I) Value date */
    int    numPts,              /* (I) Length of dates & rates */
    double basis,               /* (I) Basis. See JpmcdsRateToDiscount */
    long   dayCountConv);       /* (I) See JpmcdsDayCountFraction */


/*f
***************************************************************************
** Creates a new TCurve, and fills it in with rates, dates, yearfraction 
** function, and basis.
** 
** If it encounters a TDate = 0, it assumes that this is a terminator. 
** This allows spreadsheet users to select a big range of dates, but use 
** less than the whole range by terminating with a 1 or 0.
***************************************************************************
*/
EXPORT TCurve* JpmcdsMakeTCurve
    (TDate   baseDate,          /* (I) Value date */
     TDate  *dates,             /* (I) Array of dates for curve */
     double *rates,             /* (I) Array of rates for dates of curve */
     int     numPts,            /* (I) # elements in dates & rates arrays */
     double  basis,             /* (I) Basis. See JpmcdsRateToDiscount */
     long    dayCountConv);     /* (I) See JpmcdsDayCountFraction */


/*f
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
    TDate  baseDate,                    /* (I) Discount date  */
    TDate *dates,                       /* (I) Dates in curve */
    double*rates,                       /* (I) Associated rates */
    int    numPts,                      /* (I) Len of dates and rates arrays */
    double basis,                       /* (I) Compounding periods/year */
    long   dayCountConv);               /* (I) */


/*f
***************************************************************************
** Creates an array of dates from the input curve and returns a pointer to 
** the array on success and NULL on failure.
***************************************************************************
*/
TDate* JpmcdsDatesFromCurve
    (TCurve   *curve);   /* (I) Curve to get dates from */


/*f
***************************************************************************
** Calculate present value of a cash flow list.
***************************************************************************
*/
int JpmcdsCashFlowPV
    (TCashFlowList *cashFlows,          /* (I) Arbitrary cashflows */
     TCurve        *zeroCurve,          /* (I) Zero-coupon curve */
     long           interpType,         /* (I) For zcurve */
     double         *thePV);            /* (O) Present value of cashflows */


/*f
***************************************************************************
** Calculate discount factor for specified date and zero curve.
***************************************************************************
*/
int JpmcdsDiscountDate
    (TDate    date,                /* (I) Date */
     TCurve   *zeroCurve,          /* (I) Zero Curve */
     long     interpType,          /* (I) For zcurve. */
     double   *discount);          /* (O) Present Value */


/*f
***************************************************************************
** Frees a TCurve.
***************************************************************************
*/
void JpmcdsFreeTCurve
    (TCurve *);                    /* (I) TCurve to be deallocated */


/*f
***************************************************************************
** Creates a copy of the curve.
***************************************************************************
*/
TCurve* JpmcdsCopyCurve
    (TCurve *aCurve);               /* (I) */ 
 

/*f
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
);


#ifdef __cplusplus
}
#endif

#endif    /* TCURVE_H */

