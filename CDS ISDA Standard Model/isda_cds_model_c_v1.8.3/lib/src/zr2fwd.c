/*
 * ISDA CDS Standard Model
 *
 * Copyright (C) 2009 International Swaps and Derivatives Association, Inc.
 * Developed and supported in collaboration with Markit
 * 
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the ISDA CDS Standard Model Public License.
 */

#include "zr2fwd.h"
#include "cgeneral.h"
#include "bastypes.h"
#include "cfinanci.h"
#include "tcurve.h"
#include "cerror.h"


/*
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
    TCurve *zc,                         /* (I) Zero Curve */
    long interpType,                    /* (I) Type of interpolation */
    TDate startDate,                    /* (I) Start Date */
    TDate endDate,                      /* (I) End Date */
    long dayCountConv,                  /* (I) See ldate.h */
    long basis,                         /* (I) See cfinanci.h */
    double *rate)                       /* (O) rate from start to end */
{
    static char routine[]="JpmcdsForwardFromZCurve";
    int         status = FAILURE;               /* Until proven successful */
    double      startDiscount;
    double      endDiscount;
    double      forwardDiscount;

    if (JpmcdsDiscountDate(startDate, zc, interpType, &startDiscount) == FAILURE)
        goto done;

    if (JpmcdsDiscountDate(endDate, zc, interpType, &endDiscount) == FAILURE)
        goto done;
    
    forwardDiscount = endDiscount/startDiscount;

    if (JpmcdsDiscountToRate(forwardDiscount, startDate, endDate,
                          dayCountConv, basis, rate) == FAILURE)
        goto done;

    status = SUCCESS;

done:
    if (status == FAILURE)
        JpmcdsErrMsg("%s: failed.\n", routine);

    return status;
}
