/*
 * ISDA CDS Standard Model
 *
 * Copyright (C) 2009 International Swaps and Derivatives Association, Inc.
 * Developed and supported in collaboration with Markit
 * 
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the ISDA CDS Standard Model Public License.
 */


#include "feeleg.h"
#include <math.h>
#include "timeline.h"
#include "macros.h"
#include "cxdatelist.h"
#include "cxzerocurve.h"
#include "cxbsearch.h"
#include "ldate.h"
#include "cerror.h"
#include "cashflow.h"
#include "convert.h"



/*
***************************************************************************
** Calculates the PV of a single fee payment.
** Uses a pre-calculated timeline for efficiency.
***************************************************************************
*/
static int FeePaymentPVWithTimeLine
(TAccrualPayConv accrualPayConv,
 TDate           today,
 TDate           stepinDate,
 TDate           accStartDate,
 TDate           accEndDate,
 TDate           payDate,
 long            accrueDCC,
 double          notional,
 double          couponRate,
 TCurve         *discCurve,
 TCurve         *spreadCurve,
 TDateList      *tl,
 TBoolean        obsStartOfDay,
 double         *pv);


/*
***************************************************************************
** Calculates the accrued interest as of today
***************************************************************************
*/
static int FeeLegAI
(TFeeLeg      *fl,
 TDate         today,
 double       *ai);


/*
***************************************************************************
** Calculates the PV of a fee leg with fixed fee payments.
***************************************************************************
*/
int JpmcdsFeeLegPV
(TFeeLeg      *fl,
 TDate         today,
 TDate         stepinDate,
 TDate         valueDate,
 TCurve       *discCurve,
 TCurve       *spreadCurve,
 TBoolean      payAccruedAtStart,
 double       *pv)
{
    static char routine[] = "JpmcdsFeeLegPV";
    int         status    = FAILURE;
    int         i;
    double      myPv;
    double      valueDatePv;
    TDateList  *tl = NULL;
    TDate       matDate;

    REQUIRE (spreadCurve != NULL);

    REQUIRE (fl != NULL);
    REQUIRE (discCurve != NULL);
    REQUIRE (spreadCurve != NULL);
    REQUIRE (pv != NULL);
    REQUIRE (valueDate >= today);
    REQUIRE (stepinDate >= today);

    myPv = 0.0;

    if (fl->nbDates > 1)
    {
        /* it is more efficient to compute the timeLine just the once
           and truncate it for each payment */
        TDate startDate = fl->accStartDates[0];
        TDate endDate   = fl->accEndDates[fl->nbDates-1];
        
        if (!JpmcdsZeroPrice(spreadCurve, endDate) > JPMCDS_LOG0_THRESHOLD) 
        {
            JpmcdsErrMsg("%s: Default hazard rate is not defined at maturity %s!", 
                    routine, JpmcdsFormatDate(endDate));
            goto done;
        }

        if (!JpmcdsZeroPrice(discCurve, endDate) > JPMCDS_LOG0_THRESHOLD)
        {
            JpmcdsErrMsg("%s: Zero discount rate is not defined at maturity %s!",
                    routine, JpmcdsFormatDate(endDate));
            goto done;
        }

        tl = JpmcdsRiskyTimeLine(startDate,
                                 endDate,
                                 discCurve,
                                 spreadCurve);

        if (tl == NULL)
            goto done;
    }

    matDate = (fl->obsStartOfDay == TRUE ? 
               fl->accEndDates[fl->nbDates - 1] - 1 :
               fl->accEndDates[fl->nbDates - 1]);

    if(today > matDate || stepinDate > matDate)
    {
        status = SUCCESS;
        *pv = 0;
        goto done;
    }

    for (i = 0; i < fl->nbDates; ++i)
    {
        double thisPv = 0;
        
        if (FeePaymentPVWithTimeLine (fl->accrualPayConv,
                                      today,
                                      stepinDate,
                                      fl->accStartDates[i],
                                      fl->accEndDates[i],
                                      fl->payDates[i],
                                      fl->dcc,
                                      fl->notional,
                                      fl->couponRate,
                                      discCurve,
                                      spreadCurve,
                                      tl,
                                      fl->obsStartOfDay,
                                      &thisPv) != SUCCESS)
            goto done;

        myPv += thisPv;
    }

    valueDatePv = JpmcdsForwardZeroPrice (discCurve, today, valueDate);

    *pv = myPv / valueDatePv;
    
    if(payAccruedAtStart) /* clean price */
    {
        double ai;
        if(FeeLegAI(fl, stepinDate, &ai) == FAILURE)
        {
           JpmcdsErrMsg ("%s: accrued interest calculation failed.\n", routine);
           goto done;
        }
        *pv -= ai;
    }

    status = SUCCESS;

 done:

    if (status != SUCCESS)
        JpmcdsErrMsgFailure (routine);

    JpmcdsFreeDateList (tl);

    return status;
}


/*
***************************************************************************
** Calculates the PV of a single fee payment.
** Uses a pre-calculated timeline for efficiency.
**
** Calling this function repeatedly with sensible inputs and summing the
** result would give the same answer as using JpmcdsFeeLegPV directly
** with one difference - the pv returned by this function is for today,
** whereas the pv returned by JpmcdsFeeLegPV is for valueDate.
**
** The conversion is a matter of dividing by the discount factor between
** today and valueDate using the risk-free curve.
***************************************************************************
*/
static int FeePaymentPVWithTimeLine
(TAccrualPayConv accrualPayConv,
 TDate           today,
 TDate           stepinDate,
 TDate           accStartDate,
 TDate           accEndDate,
 TDate           payDate,
 long            accrueDCC,
 double          notional,
 double          couponRate,
 TCurve         *discCurve,
 TCurve         *spreadCurve,
 TDateList      *tl,
 TBoolean        obsStartOfDay,
 double         *pv)
{
    static char routine[] = "FeePaymentPVWithTimeLine";
    int         status    = FAILURE;

    double myPv = 0.0;

    /*
     * Because survival is calculated at the end of the day, then if
     * we observe survival at the start of the day, we need to subtract
     * one from the date.
     */
    int    obsOffset = obsStartOfDay ? -1 : 0;

    REQUIRE (discCurve != NULL);
    REQUIRE (spreadCurve != NULL);
    REQUIRE (pv != NULL);

    if(accEndDate <= stepinDate)
    {
        *pv = 0;
        return SUCCESS;
    }

    switch (accrualPayConv)
    {
    case ACCRUAL_PAY_NONE:
    {
        /* fee leg pays at pay date if it has survived to accrual end date */
        double accTime;
        double amount;
        double survival;
        double discount;
        
        if (JpmcdsDayCountFraction(accStartDate, accEndDate, accrueDCC, &accTime) != SUCCESS)
            goto done;
        
        amount   = notional * couponRate * accTime;
        survival = JpmcdsForwardZeroPrice(spreadCurve, today, accEndDate + obsOffset);
        discount = JpmcdsForwardZeroPrice(discCurve, today, payDate);
        myPv = amount * survival * discount;
        break;
    }
    case ACCRUAL_PAY_ALL:
    {
        /* fee leg pays accrual on default - otherwise it pays at pay date
           if it has survived to accrual end date */
        double accTime;
        double amount;
        double survival;
        double discount;
        double accrual;
        
        if (JpmcdsDayCountFraction(accStartDate, accEndDate, accrueDCC, &accTime) != SUCCESS)
            goto done;

        amount   = notional * couponRate * accTime;
        survival = JpmcdsForwardZeroPrice(spreadCurve, today, accEndDate + obsOffset);
        discount = JpmcdsForwardZeroPrice(discCurve, today, payDate);
        myPv = amount * survival * discount;
        
        /* also need to calculate accrual PV */
        
        if (JpmcdsAccrualOnDefaultPVWithTimeLine(today,
                                                 stepinDate + obsOffset,
                                                 accStartDate + obsOffset,
                                                 accEndDate + obsOffset,
                                                 amount,
                                                 discCurve,
                                                 spreadCurve,
                                                 tl,
                                                 &accrual) != SUCCESS)
            goto done;
        
        myPv += accrual;
        break;
    }
    default:
        JpmcdsErrMsg("%s: Invalid accrual payment type %d\n", routine, (int)accrualPayConv);
        goto done;
    }

    status = SUCCESS;
    *pv = myPv;

 done:

    if (status != SUCCESS)
        JpmcdsErrMsgFailure (routine);

    return status;
}


/*
***************************************************************************
** Calculates the PV of the accruals which occur on default.
** Uses a pre-calculated timeline for efficiency.
***************************************************************************
*/
int JpmcdsAccrualOnDefaultPVWithTimeLine
(TDate           today,
 TDate           stepinDate,
 TDate           startDate,
 TDate           endDate,
 double          amount,
 TCurve         *discCurve,
 TCurve         *spreadCurve,
 TDateList      *criticalDates,
 double         *pv)
{
    static char routine[] = "JpmcdsAccrualOnDefaultPVWithTimeLine";
    int         status    = FAILURE;

    double  myPv = 0.0;
    int     i;

    double t;
    double s0;
    double s1;
    double df0;
    double df1;
    double accRate;
    TDate  subStartDate;

    TDateList  *tl = NULL;

    REQUIRE (endDate > startDate);
    REQUIRE (discCurve != NULL);
    REQUIRE (spreadCurve != NULL);
    REQUIRE (pv != NULL);
    
    /*
    ** Timeline is points on the spreadCurve between startDate and endDate,
    ** combined with points from the discCurve, plus
    ** the startDate and endDate.
    */
    if (criticalDates != NULL)
    {
        tl = JpmcdsTruncateTimeLine (criticalDates, startDate, endDate);
    }
    else
    {
        tl = JpmcdsRiskyTimeLine (startDate, endDate, discCurve, spreadCurve);
    }
    if (tl == NULL)
        goto done;

    /* the integration - we can assume flat forwards between points on
       the timeline - this is true for both curves 

       we are integrating -Zt dS/dt where Z is the discount factor and
       S is the survival probability and t is the accrual time

       assuming flat forwards on each part of the integration, this is an
       exact integral
    */
    
    subStartDate = MAX(stepinDate, startDate);
    t       = (double)(endDate-startDate)/365.0;
    accRate = amount/t;
    s0      = JpmcdsForwardZeroPrice(spreadCurve, today, subStartDate);
    df0     = JpmcdsForwardZeroPrice(discCurve, today, MAX(today, subStartDate));

    for (i = 1; i < tl->fNumItems; ++i)
    {
        double lambda;
        double fwdRate;
        double thisPv = 0.;
        double t0;
        double t1;
        double lambdafwdRate;
        if(tl->fArray[i] <= stepinDate)
            continue;

        s1  = JpmcdsForwardZeroPrice(spreadCurve, today, tl->fArray[i]);
        df1 = JpmcdsForwardZeroPrice(discCurve, today, tl->fArray[i]);

        t0  = (double)(subStartDate + 0.5 - startDate)/365.0;
        t1  = (double)(tl->fArray[i] + 0.5- startDate)/365.0;
        t   = t1-t0;

		/*************************Markit Proposed Fix***************************************************
		 *
		 * Some of the division of original ISDA model can be removed
		 * lambda  = log(s0 / s1) / t = (log(s0) - log(s1)) / t
		 * fwdRate = log(df0 / df) / t = (log(df0) - log(df1)) / t
		 * Divisions by t can be absorbed by later formulas as well.
		 */
        lambda  = log(s0) - log(s1);
        fwdRate = log(df0) - log(df1);
        lambdafwdRate = lambda + fwdRate + 1.0e-50;


		/**
		 * If lambdafwdRate is extremely small, original calculation generates big noise on computer 
		 * due to the small denominators.
		 * In this case, Talyor expansion is employed to remove lambdafwdRate from denomintors 
		 * so that numerical noise is signicantly reduced. 
		 */
		if (fabs(lambdafwdRate) > 1e-4)
		{
			/*This is the original formula which contains an integral*/
			thisPv  = lambda * accRate * s0 * df0 * ( \
				(t0 + t/(lambdafwdRate))/(lambdafwdRate) - \
				(t1 + t/(lambdafwdRate))/(lambdafwdRate) * \
				s1/s0 * df1/df0);
			
			/* This is the accrual on default formula fix
			thisPv  = lambda * accRate * s0 * df0 * t * ( \
				1.0 / lambdafwdRate / lambdafwdRate - \
				(1.0 + 1.0 / lambdafwdRate) / lambdafwdRate * \
				s1 / s0 * df1 / df0);
			*/
		}
		else
		{
			/* 
			This is the numerical fix corresponding to the original formula
			*/
			const double lambdaAccRate = lambda * s0 * df0 * accRate * 0.5;
			const double thisPv1 = lambdaAccRate * (t0 + t1);

			const double lambdaAccRateLamdaFwdRate = lambdaAccRate * lambdafwdRate / 3.;
			const double thisPv2 = -lambdaAccRateLamdaFwdRate * (t0 + 2. * t1);

			const double lambdaAccRateLamdaFwdRate2 = lambdaAccRateLamdaFwdRate * lambdafwdRate * .25;
			const double thisPv3 = lambdaAccRateLamdaFwdRate2 * (t0 + 3. * t1);

			const double lambdaAccRateLamdaFwdRate3 = lambdaAccRateLamdaFwdRate2 * lambdafwdRate * .2;
			const double thisPv4 = -lambdaAccRateLamdaFwdRate3 * (t0 + 4. * t1);

			const double lambdaAccRateLamdaFwdRate4 = lambdaAccRateLamdaFwdRate3 * lambdafwdRate / 6.;
			const double thisPv5 = lambdaAccRateLamdaFwdRate4 * (t0 + 5. * t1);
			
			/* This is the numerical fix along with accrual on default model fix
			const double lambdaAccRate = lambda * s0 * df0 * accRate * t;
			const double thisPv1 = lambdaAccRate * 0.5;

			const double lambdaAccRateLamdaFwdRate = lambdaAccRate * lambdafwdRate;
			const double thisPv2 = -lambdaAccRateLamdaFwdRate / 3.;

			const double lambdaAccRateLamdaFwdRate2 = lambdaAccRateLamdaFwdRate * lambdafwdRate;
			const double thisPv3 = lambdaAccRateLamdaFwdRate2 * .125;

			const double lambdaAccRateLamdaFwdRate3 = lambdaAccRateLamdaFwdRate2 * lambdafwdRate;
			const double thisPv4 = -lambdaAccRateLamdaFwdRate3 / 30.;

			const double lambdaAccRateLamdaFwdRate4 = lambdaAccRateLamdaFwdRate3 * lambdafwdRate;
			const double thisPv5 = lambdaAccRateLamdaFwdRate4 / 144.;
			*/
			
			thisPv += thisPv1;
			thisPv += thisPv2;
			thisPv += thisPv3;
			thisPv += thisPv4;
			thisPv += thisPv5;

		}

        myPv += thisPv;
        s0  = s1;
        df0 = df1;
        subStartDate = tl->fArray[i];
    }

    status = SUCCESS;
    *pv = myPv;
        
 done:

    if (status != SUCCESS)
        JpmcdsErrMsgFailure (routine);

    JpmcdsFreeDateList (tl);
    return status;
}


/*
***************************************************************************
** Computes the non-contingent cash flows for a fee leg.
***************************************************************************
*/
TCashFlowList* JpmcdsFeeLegFlows
(TFeeLeg      *fl)
{
    static char routine[] = "JpmcdsFeeLegFlows";
    int         status    = FAILURE;

    TCashFlowList *cfl = NULL;
    int i;

    cfl = JpmcdsNewEmptyCFL(fl->nbDates);
    if (cfl == NULL)
        goto done;

    for (i = 0; i < fl->nbDates; ++i)
    {
        double amount;
        double time;

        if(JpmcdsDayCountFraction(fl->accStartDates[i],
                               fl->accEndDates[i],
                               fl->dcc,
                               &time) != SUCCESS)
        {
            goto done;
        }

        amount = time * fl->couponRate * fl->notional;

        cfl->fArray[i].fDate   = fl->payDates[i];
        cfl->fArray[i].fAmount = amount;
    }

    status = SUCCESS;

 done:

    if (status != SUCCESS)
    {
        JpmcdsFreeCFL(cfl);
        cfl = NULL;
        JpmcdsErrMsgFailure(routine);
    }

    return cfl;
}



/*
***************************************************************************
** Calculates the accrued interest as of today
***************************************************************************
*/
static int FeeLegAI
(TFeeLeg      *fl,
 TDate         today,
 double       *ai)
{
    static char routine[] = "FeeLegAI";
    int         status    = FAILURE;

    long        exact, lo, hi;

    if(today <= fl->accStartDates[0] ||
       today >= fl->accEndDates[fl->nbDates -1])
    {
        *ai = 0;
        status = SUCCESS;
        goto done;
    }

    if(JpmcdsBinarySearchLong(today,
                          fl->accStartDates,
                          sizeof(TDate),
                          fl->nbDates,
                          &exact,
                          &lo,
                          &hi) != SUCCESS)
        goto done;

    if(exact >= 0)
    {
        /* today is on accrual start date */
        *ai = 0;
    }
    else
    {
        /* calculate ai, today is bracketed in some accrual period */
        double accrual;
        if(JpmcdsDayCountFraction(fl->accStartDates[lo], today, fl->dcc, &accrual) != SUCCESS)
            goto done;

        accrual *= fl->couponRate * fl->notional;
        *ai = accrual;
    }
    status = SUCCESS;

 done:
    if (status != SUCCESS)
        JpmcdsErrMsgFailure(routine);

    return status;
}
