/*
 * ISDA CDS Standard Model
 *
 * Copyright (C) 2009 International Swaps and Derivatives Association, Inc.
 * Developed and supported in collaboration with Markit
 * 
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the ISDA CDS Standard Model Public License.
 */

#include <string.h>
#include <math.h>
#include <float.h>
#include "zr2coup.h"
#include "cgeneral.h"
#include "bastypes.h"
#include "cerror.h"
#include "cashflow.h"
#include "datelist.h"
#include "ldate.h"
#include "date_sup.h"
#include "convert.h"
#include "tcurve.h"
#include "macros.h"


/*
***************************************************************************
** Calculates the coupon rate for an instrument, such that cashflows 
** discounted by the zero-rates equal the present value.
**
** Uses:
** 1) zero-coupon rates in the zc structure,
** 2) cashFlowList, which contains coupons payments for an instrument based
**    on an annual rate of 1 (100%), (in other words,  fraction of annual
**    coupon which should be paid at the given date)
** 3) the present Value of the instrument
** 
** It is assumed that the last payment represents a coupon plus principal.
***************************************************************************
*/
static int JpmcdsCalcCoupon(
  TCurve        *zc, 
  long           interpType, 
  TCashFlowList *cfl,
  double         presentValue, 
  double        *couponRate);


/*
***************************************************************************
** A convenience function for JpmcdsZerosToCouponsPointAdj with no bad-day 
** adjustment.
**
** Calls JpmcdsZerosToCouponsPoint with badDayConv = JPMCDS_BAD_DAY_NONE.
***************************************************************************
*/
int JpmcdsZerosToCouponsPoint
   (TCurve        *zc,                  /* (I) Zero-coupon rate curve */
    long           interpType,          /* (I) */
    TDate          startDate,           /* (I) Date instrument begins at */
    TDateInterval *interval,            /* (I) Time between payments */
    TDate          maturityDate,        /* (I) Date instrument matures at*/
    long           fixedDCC,            /* (I) See JpmcdsDayCountFraction */
    TBoolean       stubAtEnd,           /* (I) TRUE=Back, FALSE=Front*/
    double        *couponRate)          /* (O) Rate for this instrument */
{
    static char routine[]="JpmcdsZerosToCouponsPoint";
    int         status = FAILURE; 

    /* Turn off bad day adjustment. */
    if (JpmcdsZerosToCouponsPointAdj
        (zc, interpType, 
         startDate, interval, maturityDate, fixedDCC,
         stubAtEnd, 
         JPMCDS_BAD_DAY_NONE, JPMCDS_BAD_DAY_NONE, "NO_WEEKENDS",
         couponRate) == FAILURE)
        goto done;

    status = SUCCESS;

done:
    if (status == FAILURE)
        JpmcdsErrMsg("%s: Failed.\n", routine);

    return status;
}


/*
***************************************************************************
** Calculates a par swap rate for the swap starting at startDate,
** maturing at maturityDate, with fixed day count convention 
** fixedDayCountConv and fixed payments occuring at time 
** intervals defined by interval. In other words, the routine calculates 
** the fixed rate such that the present value of the fixed and 
** floating sides are equal. The floating side is assumed to be at par.
**
** Bad day adjustment is performed on the startDate, maturityDate,
** and all fixed coupon payment dates. 
**
** The stubAtEnd argument only has an effect if startDate and maturityDate 
** are not on cycle. 
***************************************************************************
*/
int JpmcdsZerosToCouponsPointAdj
   (TCurve        *zc,                /* (I) Zero-coupon rate curve */
    long           interpType,        /* (I) */
    TDate          startDate,         /* (I) Date instrument begins at */
    TDateInterval *interval,          /* (I) Time between payments */
    TDate          maturityDate,      /* (I) Date instrument matures at*/
    long           fixedDCC,          /* (I) See JpmcdsDayCountFraction */
    TBoolean       stubAtEnd,         /* (I) TRUE=Back, FALSE=Front*/
    long           accBadDayConv,     /* (I) See JpmcdsBusinessDay for methods*/
    long           payBadDayConv,     /* (I) See JpmcdsBusinessDay for methods*/
    char          *holidayFile,       /* (I) Used to adjust swap dates */
    double        *couponRate)        /* (O) Rate for this instrument */
{
    static char routine[]="JpmcdsZerosToCouponsPointAdj";
    int            status = FAILURE; 
    TCashFlowList *cfl = (TCashFlowList *)NULL;
    TDate          startDateAdj;      /* Adjusted for bad days */
    double         startDatePV;       /* Present value of $1 at startDate */
    TDate          zcFirstDate;       /* First valid date in curve */


    /* Move startDate to a good business day. */
    if (JpmcdsBusinessDay(startDate, accBadDayConv, holidayFile, &startDateAdj) == FAILURE)
        goto done;

    if (JpmcdsCurveFirstDate (zc, &zcFirstDate) != SUCCESS)
        goto done;

    if (startDateAdj < zcFirstDate)
    {
        JpmcdsErrMsg("%s: Start date (%s) <  first valid zero date (%s).\n",
                  routine, JpmcdsFormatDate(startDateAdj),
                  JpmcdsFormatDate(zcFirstDate));

        if (accBadDayConv != JPMCDS_BAD_DAY_NONE && startDate != startDateAdj)
        {
            JpmcdsErrMsg
                ("%s: Bad day adjustment moved start date from %s to %s.\n",
                 routine, JpmcdsFormatDate(startDate), 
                 JpmcdsFormatDate(startDateAdj));
        }
        goto done;
    }

    /* Make cash flow list. Don't pass startDateAdj, because the 
     * routine will adjust that date itself.
     */
    cfl = JpmcdsMakeCFL(1. /* coupon*/, startDate, interval, maturityDate,
                     fixedDCC, stubAtEnd,
                     0 /* No initial payment, no final repayment */,
                     accBadDayConv, payBadDayConv, holidayFile);
     
    if (cfl == (TCashFlowList *)NULL)
        goto done;

    /* Get present value of $1 at startDateAdj. */
    if (JpmcdsDiscountDate(startDateAdj, zc, interpType, &startDatePV) == FAILURE)
        goto done;

    /* Compute coupon for the given zero-coupon rates.
    * Note that because the amount field for the 1st cashFlow is 0,
    * (which corresponds to the value date,) it is not taken into
    * account.
    */
    status = JpmcdsCalcCoupon(zc, interpType, cfl, startDatePV, couponRate);
        
done:
    if (status == FAILURE)
        JpmcdsErrMsg("%s: Failed.\n", routine);
    JpmcdsFreeCFL(cfl);    

    return status;
}


/*
***************************************************************************
** Calculates the coupon rate for an instrument, such that cashflows 
** discounted by the zero-rates equal the present value.
**
** Uses:
** 1) zero-coupon rates in the zc structure,
** 2) cashFlowList, which contains coupons payments for an instrument based
**    on an annual rate of 1 (100%), (in other words,  fraction of annual
**    coupon which should be paid at the given date)
** 3) the present Value of the instrument
** 
** It is assumed that the last payment represents a coupon plus principal.
***************************************************************************
*/
static int JpmcdsCalcCoupon(
    TCurve *zc,                         /* (I) Zero-coupon rate curve */
    long interpType,                    /* (I) Interpolation */
    TCashFlowList *cfl,                 /* (I) Cash Flow list */
    double presentValue,                /* (I) Present value */
    double *couponRate)                 /* (O) Annual coupon rate */
{
    static char routine[]="JpmcdsCalcCoupon";
    int         status = FAILURE;

    int num_pts = cfl->fNumItems;
    double coupons_PV;                  /* Present value of coupons */
    double last_PV;                     /* Present value of final payment */
    
    if (JpmcdsCashFlowPV(cfl, zc, interpType, &coupons_PV) == FAILURE)
        goto done;

    if (coupons_PV == 0.)
    {
        JpmcdsErrMsg("JpmcdsCalcCoupon: coupons have zero present value.\n");
        goto done;
    }

    if (JpmcdsDiscountDate(cfl->fArray[num_pts-1].fDate, zc, interpType, &last_PV) == FAILURE)
        goto done;

    *couponRate = (presentValue - last_PV) / coupons_PV;

    status = SUCCESS;

  done:
    if (status == FAILURE)
        JpmcdsErrMsg("%s: Failed.\n", routine);

    return status;
}
