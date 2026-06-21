/*
 * ISDA CDS Standard Model
 *
 * Copyright (C) 2009 International Swaps and Derivatives Association, Inc.
 * Developed and supported in collaboration with Markit
 * 
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the ISDA CDS Standard Model Public License.
 */

#ifndef ZR2FWD_H
#define ZR2FWD_H

#include "bastypes.h"
#include "ldate.h"
#include "cfinanci.h"
#include "interp.h"


#ifdef __cplusplus
extern "C"
{
#endif


/*f
***************************************************************************
** Computes the zero-coupon forward rate from startDate to the maturity 
** date using a zero-coupon curve.
**
** Note the routine does not check to make sure that the start date comes 
** before the maturity date, since zero curve creation requires the 
** calculation of forwards which go backwards (JpmcdsDiscountToRate should 
** handle that check anyway).
***************************************************************************
*/
int JpmcdsForwardFromZCurve(
    TCurve *zeroCurve,                  /* (I) Zero Curve */
    long    interpType,                 /* (I) */
    TDate   startDate,                  /* (I) Start Date */
    TDate   endDate,                    /* (I) End Date */
    long    dayCountConv,               /* (I) See JpmcdsDayCountFraction */
    long    basis,                      /* (I) See JpmcdsDiscountToRate */
    double *rate);                      /* (O) Rate from start to end */


#ifdef __cplusplus
}
#endif

#endif    /* ZR2FWD_H */
