/*
 * ISDA CDS Standard Model
 *
 * Copyright (C) 2009 International Swaps and Derivatives Association, Inc.
 * Developed and supported in collaboration with Markit
 * 
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the ISDA CDS Standard Model Public License.
 */

#include "cds.h"
#include "contingentleg.h"
#include "feeleg.h"
#include "cxzerocurve.h"
#include "cxdatelist.h"
#include "rtbrent.h"
#include "convert.h"
#include "tcurve.h"
#include "ldate.h"
#include "macros.h"
#include "cerror.h"


typedef struct
{
    int             i;
    TDate           stepinDate;
    TDate           cashSettleDate;
    TCurve         *discountCurve;
    TCurve         *cdsCurve;
    double          recoveryRate;
    TContingentLeg *cl;
    TFeeLeg        *fl;
} CDS_BOOTSTRAP_CONTEXT;


/* static function declarations */
static int cdsBootstrapPointFunction
(double   cleanSpread,
 void    *data,
 double  *pv);


static TCurve* CdsBootstrap
(TDate           today,           /* (I) Used as credit curve base date     */
 TCurve         *discountCurve,   /* (I) Risk-free discount curve           */
 TDate           startDate,       /* (I) Start of CDS for accrual and risk  */
 TDate           stepinDate,      /* (I) Stepin date for accrual etc        */
 TDate           cashSettleDate,  /* (I) Pay date                           */
 long            nbDate,          /* (I) Number of benchmark dates          */
 TDate          *endDates,        /* (I) Maturity dates of CDS to bootstrap */
 double         *couponRates,     /* (I) Coupons (e.g. 0.05 = 5% = 500bp)   */ 
 double          recoveryRate,    /* (I) Recovery rate                      */
 TBoolean        payAccOnDefault, /* (I) Pay accrued on default             */
 TDateInterval  *couponInterval,  /* (I) Interval between fee payments      */
 long            paymentDCC,      /* (I) DCC for fee payments and accrual   */
 TStubMethod    *stubType,        /* (I) Stub type for fee leg              */
 long            badDayConv,
 char           *calendar);

/**
***************************************************************************
** Routine for converting the compounding basis of a credit curve.
**
** The curve is amended in place.
***************************************************************************
*/
int CreditCurveConvertRateType(TCurve *cc, long newRateType);


/*
***************************************************************************
** The main bootstrap routine.
***************************************************************************
*/
EXPORT TCurve* JpmcdsCleanSpreadCurve
(TDate              today,           /* (I) Used as credit curve base date       */
 TCurve            *discountCurve,   /* (I) Risk-free discount curve             */
 TDate              startDate,       /* (I) Start of CDS for accrual and risk    */
 TDate              stepinDate,      /* (I) Stepin date                          */
 TDate              cashSettleDate,  /* (I) Pay date                             */
 long               nbDate,          /* (I) Number of benchmark dates            */
 TDate             *endDates,        /* (I) Maturity dates of CDS to bootstrap   */
 double            *couponRates,     /* (I) CouponRates (e.g. 0.05 = 5% = 500bp) */ 
 TBoolean          *includes,        /* (I) Include this date. Can be NULL if    
                                        all are included.                        */
 double             recoveryRate,    /* (I) Recovery rate                        */
 TBoolean           payAccOnDefault, /* (I) Pay accrued on default               */
 TDateInterval     *couponInterval,  /* (I) Interval between fee payments        */
 long               paymentDCC,      /* (I) DCC for fee payments and accrual     */
 TStubMethod       *stubType,        /* (I) Stub type for fee leg                */
 long               badDayConv,
 char              *calendar
)
{
    static char routine[] = "JpmcdsCleanSpreadCurve";
    TCurve *out = NULL;

    TDate           *includeEndDates = NULL;
    double          *includeCouponRates = NULL;

    TDateInterval ivl3M;

    SET_TDATE_INTERVAL(ivl3M,3,'M');
    if (couponInterval == NULL)
        couponInterval = &ivl3M;

    /* put common requirements between two routines here */
    REQUIRE (discountCurve != NULL);
    REQUIRE (nbDate > 0);
    REQUIRE (endDates != NULL);
    REQUIRE (couponRates != NULL);

    if (includes != NULL)
    {
        /* need to pick and choose which names appear */
        int  i;
        long nbInclude = 0;
        long j;
        for (i = 0; i < nbDate; ++i) 
        {
            if (includes[i]) ++nbInclude;
        }
        REQUIRE (nbInclude > 0);

        includeEndDates    = NEW_ARRAY(TDate, nbInclude);
        includeCouponRates = NEW_ARRAY(double, nbInclude);
        
        j = 0;
        for (i = 0; i < nbDate; ++i)
        {
            if (includes[i])
            {
                includeEndDates[j]    = endDates[i];
                includeCouponRates[j] = couponRates[i];
                ++j;
            }
        }
        ASSERT (j == nbInclude);
        
        nbDate      = nbInclude;
        endDates    = includeEndDates;
        couponRates = includeCouponRates;
    }

    out = CdsBootstrap (today,
                        discountCurve,
                        startDate,
                        stepinDate,
                        cashSettleDate,
                        nbDate,
                        endDates,
                        couponRates,
                        recoveryRate,
                        payAccOnDefault,
                        couponInterval,
                        paymentDCC,
                        stubType,
                        badDayConv,
                        calendar);

 done:
    FREE(includeEndDates);
    FREE(includeCouponRates);
    if (out == NULL)
        JpmcdsErrMsgFailure (routine);

    return out;
}


/*
***************************************************************************
** This is the CDS bootstrap routine.
**
** Very little attempt has been made at extreme optimisation - this is quite
** a basic bootstrap routine which simply calls the underlying CDS pricer
** for each benchmark instrument while it changes the CDS zero rate at the
** maturity date of the benchmark instrument.
**
***************************************************************************
*/
static TCurve* CdsBootstrap
(TDate             today,           /* (I) Used as credit curve base date     */
 TCurve           *discountCurve,   /* (I) Risk-free discount curve           */
 TDate             startDate,       /* (I) Start of CDS for accrual and risk  */
 TDate             stepinDate,      /* (I) Stepin date                        */
 TDate             cashSettleDate,  /* (I) Pay date                           */
 long              nbDate,          /* (I) Number of benchmark dates          */
 TDate            *endDates,        /* (I) Maturity dates of CDS to bootstrap */
 double           *couponRates,     /* (I) CouponRates (e.g. 0.05 = 5% = 500bp)   */ 
 double            recoveryRate,    /* (I) Recovery rate                      */
 TBoolean          payAccOnDefault, /* (I) Pay accrued on default         */
 TDateInterval    *couponInterval,  /* (I) Interval between fee payments      */
 long              paymentDCC,      /* (I) DCC for fee payments and accrual   */
 TStubMethod      *stubType,        /* (I) Stub type for fee leg              */
 long              badDayConv,
 char             *calendar)
{
    static char routine[] = "CdsBootstrap";
    int         status    = FAILURE;

    TCurve     *cdsCurve = NULL;
    int         i;
    
    CDS_BOOTSTRAP_CONTEXT context;
    TContingentLeg *cl = NULL;
    TFeeLeg        *fl = NULL;
    double          settleDiscount = 0.0;
    TBoolean        protectStart = TRUE;

    /* we work with a continuously compounded curve since that is faster -
       but we will convert to annual compounded since that is traditional */
    cdsCurve = JpmcdsMakeTCurve (today,
                                 endDates,
                                 couponRates,
                                 nbDate,
                                 JPMCDS_CONTINUOUS_BASIS,
                                 JPMCDS_ACT_365F);
    if (cdsCurve == NULL)
        goto done;

    context.discountCurve = discountCurve;
    context.cdsCurve      = cdsCurve;
    context.recoveryRate  = recoveryRate;
    context.stepinDate    = stepinDate;
    context.cashSettleDate = cashSettleDate;
    
    for (i = 0; i < nbDate; ++i)
    {
        double guess;
        double spread;

        guess = couponRates[i] / (1.0 - recoveryRate);

        cl = JpmcdsCdsContingentLegMake (MAX(today, startDate),
                                         endDates[i],
                                         1.0,
                                         protectStart);
        if (cl == NULL)
            goto done;

        fl = JpmcdsCdsFeeLegMake(startDate,
                                 endDates[i],
                                 payAccOnDefault,
                                 couponInterval,
                                 stubType,
                                 1.0,
                                 couponRates[i],
                                 paymentDCC,
                                 badDayConv,
                                 calendar,
                                 protectStart);
        if (fl == NULL)
            goto done;

        context.i  = i;
        context.cl = cl;
        context.fl = fl;

        if (JpmcdsRootFindBrent ((TObjectFunc)cdsBootstrapPointFunction,
                                 (void*) &context,
                                 0.0,    /* boundLo */
                                 1e10,   /* boundHi */
                                 100,    /* numIterations */
                                 guess,
                                 0.0005, /* initialXstep */
                                 0,      /* initialFDeriv */
                                 1e-10,  /* xacc */
                                 1e-10,  /* facc */
                                 &spread) != SUCCESS)
        {
            JpmcdsErrMsg ("%s: Could not add CDS maturity %s spread %.2fbp\n",
                          routine,
                          JpmcdsFormatDate(endDates[i]),
                          1e4 * couponRates[i]);
            goto done;
        }
        cdsCurve->fArray[i].fRate = spread;

        FREE(cl);
        JpmcdsFeeLegFree (fl);
        cl = NULL;
        fl = NULL;

        /** check if forward hazard rate is negative */
        if(i > 0)
        {
            double fwdPrice = JpmcdsForwardZeroPrice(cdsCurve, 
                                                     endDates[i-1],
                                                     endDates[i]);

            if(fwdPrice > 1)
            {
                JpmcdsErrMsg ("%s: Negative forward hazard rate at maturity %s with "
                              "spread %.2fbp\n",
                              routine,
                              JpmcdsFormatDate(endDates[i]),
                              1e4 * couponRates[i]);
                goto done;
            }
        }
    }

    if (CreditCurveConvertRateType (cdsCurve, JPMCDS_ANNUAL_BASIS) != SUCCESS)
        goto done;

    status = SUCCESS;

 done:

    if (status != SUCCESS)
    {
        JpmcdsFreeTCurve(cdsCurve);
        cdsCurve = NULL;
        JpmcdsErrMsgFailure (routine);
    }

    FREE(cl);
    JpmcdsFeeLegFree (fl);
        
    return cdsCurve;
}


/*
***************************************************************************
** Objective function for root-solver.
** Returns the PV as at the startDate of the CDS curve.
***************************************************************************
*/
static int cdsBootstrapPointFunction
(double   cleanSpread,
 void    *data,
 double  *pv)
{
    static char routine[] = "cdsBootstrapPointFunction";
    int         status    = FAILURE;

    CDS_BOOTSTRAP_CONTEXT *context = (CDS_BOOTSTRAP_CONTEXT*)data;

    int             i             = context->i;
    TCurve         *discountCurve = context->discountCurve;
    TCurve         *cdsCurve      = context->cdsCurve;
    double          recoveryRate  = context->recoveryRate;
    TContingentLeg *cl            = context->cl;
    TFeeLeg        *fl            = context->fl;
    TDate           cdsBaseDate   = cdsCurve->fBaseDate;
    TDate           stepinDate    = context->stepinDate;
    TDate           cashSettleDate = context->cashSettleDate;
    TBoolean        isPriceClean  = 1;

    double          pvC; /* PV of contingent leg */
    double          pvF; /* PV of fee leg */

    cdsCurve->fArray[i].fRate = cleanSpread;

    if (JpmcdsContingentLegPV (cl,
                               cdsBaseDate,
                               cashSettleDate,
                               stepinDate,
                               discountCurve,
                               cdsCurve,
                               recoveryRate,
                               &pvC) != SUCCESS)
        goto done;
                              
    if (JpmcdsFeeLegPV(fl,
                       cdsBaseDate,
                       stepinDate,
                       cashSettleDate,
                       discountCurve,
                       cdsCurve,
                       isPriceClean,
                       &pvF) != SUCCESS)
        goto done;

    /* Note: price is discounted to cdsBaseDate */
    *pv = pvC - pvF;
    status = SUCCESS;

 done:

    if (status != SUCCESS)
        JpmcdsErrMsgFailure (routine);

    return status;
}


/*
***************************************************************************
** Routine for converting the compounding basis of a credit curve.
** The curve is amended in place.
***************************************************************************
*/
int CreditCurveConvertRateType(TCurve *cc, long newRateType)
{
    static char routine[] = "CreditCurveConvertRateType";
    int         status    = FAILURE;

    int i;

    REQUIRE(cc != NULL);

    if (IS_EQUAL(newRateType, cc->fBasis))
    {
        /* do nothing */
    }
    else
    {
        for (i = 0; i < cc->fNumItems; ++i)
        {
            if (JpmcdsConvertCompoundRate(cc->fArray[i].fRate,
                                          cc->fBasis,
                                          cc->fDayCountConv,
                                          newRateType,
                                          cc->fDayCountConv,
                                          &cc->fArray[i].fRate) != SUCCESS)
                goto done;
        }
        cc->fBasis = newRateType;
    }
    
    status = SUCCESS;

 done:

    if (status != SUCCESS)
        JpmcdsErrMsgFailure(routine);

    return status;
}
