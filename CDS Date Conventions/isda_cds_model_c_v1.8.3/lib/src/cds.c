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
#include "cxdatelist.h"
#include "macros.h"
#include "cxzerocurve.h"
#include "busday.h"
#include "dtlist.h"
#include "cerror.h"


/*
***************************************************************************
** Makes a contingent leg for a vanilla CDS
**
** Protection starts either at the beginning of startDate (protectStart=True)
** or at the end of startDate.
**
** Protection ends at the end of endDate.
**
** Notional is the amount of notional protected.
***************************************************************************
*/
TContingentLeg* JpmcdsCdsContingentLegMake
(TDate     startDate,
 TDate     endDate,  
 double    notional, 
 TBoolean  protectStart)
{
    static char routine[] = "JpmcdsCdsContingentLegMake";
    int         status    = FAILURE;
    
    TContingentLeg *cl = NULL;
    
    cl = JpmcdsContingentLegMake(startDate,
                                 endDate,
                                 notional,
                                 PROT_PAY_DEF,
                                 protectStart);
    if (cl == NULL)
        goto done;

    status = SUCCESS;

 done:

    if (status != SUCCESS)
    {
        FREE(cl);
        cl = NULL;
        JpmcdsErrMsgFailure(routine);
    }

    return cl;
}


/*
***************************************************************************
** Computes the PV for a contingent leg for a vanilla CDS.
***************************************************************************
*/
int JpmcdsCdsContingentLegPV
(TDate             today,
 TDate             valueDate,
 TDate             startDate,
 TDate             endDate,
 double            notional,
 TCurve           *discCurve,
 TCurve           *spreadCurve,
 double            recoveryRate,
 TBoolean          protectStart,
 double           *pv)
{
    static char routine[] = "JpmcdsCdsContingentLegPV";
    int         status    = FAILURE;

    TContingentLeg *cl = NULL;

    cl = JpmcdsCdsContingentLegMake (startDate, endDate, notional, protectStart);
    if (cl == NULL)
        goto done;

    if (JpmcdsContingentLegPV (cl, today, valueDate, startDate, discCurve, spreadCurve, 
                               recoveryRate, pv) != SUCCESS)
        goto done;

    status = SUCCESS;

 done:

    if (status != SUCCESS)
        JpmcdsErrMsgFailure (routine);

    FREE(cl);

    return status;
}


/*
***************************************************************************
** Makes a fixed fee leg for a vanilla CDS.
** Note that you are protected on both startDate and endDate.
***************************************************************************
*/
TFeeLeg* JpmcdsCdsFeeLegMake
(TDate           startDate,
 TDate           endDate,
 TBoolean        payAccOnDefault,
 TDateInterval  *dateInterval,
 TStubMethod    *stubType,
 double          notional,
 double          couponRate,
 long            paymentDcc,
 long            badDayConv,
 char           *calendar,
 TBoolean        protectStart)
{
    static char routine[] = "JpmcdsCdsFeeLegMake";
    int         status    = FAILURE;

    TDateList    *dl = NULL;
    TFeeLeg      *fl = NULL;
    TDateInterval ivl3M;
    TDate         prevDate;
    TDate         prevDateAdj;

    int i;

    SET_TDATE_INTERVAL(ivl3M,3,'M');
    if (dateInterval == NULL)
        dateInterval = &ivl3M;

    if(protectStart)
        REQUIRE (endDate >= startDate);
    else
        REQUIRE (endDate > startDate);

    if(protectStart && endDate == startDate)
    {
        TDate dates[2];
        dates[0] = startDate;
        dates[1] = endDate;
        dl = JpmcdsNewDateListFromDates(dates, 2);
    } else
    {
        dl = JpmcdsDateListMakeRegular (startDate, endDate, dateInterval, stubType);
    }
    if (dl == NULL)
        goto done;

    /* the datelist includes both start date and end date */
    /* therefore it has one more element than the fee leg requires */

    fl = JpmcdsFeeLegMakeEmpty (dl->fNumItems-1);
    if (fl == NULL)
        goto done;

    if (payAccOnDefault)
    {
        fl->accrualPayConv = ACCRUAL_PAY_ALL;
    }
    else
    {
        fl->accrualPayConv = ACCRUAL_PAY_NONE;
        /* and we will assume that it observes at end of the period */
    }
    fl->dcc = paymentDcc;

    prevDate = dl->fArray[0];
    prevDateAdj = prevDate; /* first date is not bad day adjusted */

    for (i = 0; i < fl->nbDates; ++i)
    {
        TDate nextDate = dl->fArray[i+1];
        TDate nextDateAdj;

        if (JpmcdsBusinessDay (nextDate, badDayConv, calendar, &nextDateAdj) != SUCCESS)
            goto done;

        fl->accStartDates[i] = prevDateAdj;
        fl->accEndDates[i]   = nextDateAdj;
        fl->payDates[i]      = nextDateAdj;

        prevDate    = nextDate;
        prevDateAdj = nextDateAdj;
    }

    fl->notional    = notional;
    fl->couponRate  = couponRate;

    /* the last accrual date is not adjusted */
    /* also we may have one extra day of accrued interest */
    if (protectStart)
    {
        fl->accEndDates[fl->nbDates-1] = prevDate+1;
        fl->obsStartOfDay              = TRUE;
    }
    else
    {
        fl->accEndDates[fl->nbDates-1] = prevDate;
        fl->obsStartOfDay              = FALSE;
    }
    status = SUCCESS;

 done:

    if (status != SUCCESS)
    {
        JpmcdsErrMsgFailure(routine);
        JpmcdsFeeLegFree (fl);
        fl = NULL;
    }

    JpmcdsFreeDateList (dl);

    return fl;
}


/*
***************************************************************************
** Computes the PV for a fee leg for a vanilla CDS.
** Note that you are protected on both start date and end date.
***************************************************************************
*/
int JpmcdsCdsFeeLegPV
(TDate           today,
 TDate           valueDate,
 TDate           stepinDate,
 TDate           startDate,
 TDate           endDate,
 TBoolean        payAccOnDefault,
 TDateInterval  *dateInterval,
 TStubMethod    *stubType,
 double          notional,
 double          couponRate,
 long            paymentDcc,
 long            badDayConv,
 char           *calendar,
 TCurve         *discCurve,
 TCurve         *spreadCurve,
 TBoolean        protectStart,
 TBoolean        cleanPrice,
 double         *pv)
{
    static char routine[] = "JpmcdsCdsFeeLegPV";
    int         status    = FAILURE;

    TFeeLeg *fl = NULL;

    fl = JpmcdsCdsFeeLegMake (startDate, endDate, payAccOnDefault,
                              dateInterval, stubType, notional,
                              couponRate, paymentDcc, badDayConv,
                              calendar, protectStart);
    if (fl == NULL)
        goto done;

    if (JpmcdsFeeLegPV (fl, today, stepinDate, valueDate, discCurve, spreadCurve, 
                        cleanPrice, pv) != SUCCESS)
        goto done;

    status = SUCCESS;

 done:

    if (status != SUCCESS)
        JpmcdsErrMsgFailure (routine);

    JpmcdsFeeLegFree (fl);
    return status;
}


/*
***************************************************************************
** Computes the price for a vanilla CDS.
** Note that you are protected on both start date and end date.
***************************************************************************
*/
EXPORT int JpmcdsCdsPrice
(TDate             today,
 TDate             settleDate,
 TDate             stepinDate,
 TDate             startDate,
 TDate             endDate,
 double            couponRate,
 TBoolean          payAccOnDefault,
 TDateInterval    *dateInterval,
 TStubMethod      *stubType,
 long              paymentDcc,
 long              badDayConv,
 char             *calendar,
 TCurve           *discCurve,
 TCurve           *spreadCurve,
 double            recoveryRate,
 TBoolean          isPriceClean,
 double           *price)
{
    static char routine[] = "JpmcdsCdsPrice";
    int         status    = FAILURE;

    double      feeLegPV = 0;
    double      contingentLegPV = 0;
    TDate       valueDate;
    TBoolean    protectStart = TRUE;

    REQUIRE(price != NULL);
    REQUIRE(stepinDate >= today);

    /* all other requirements can be handled by the routines we call */

    valueDate = settleDate;

    if (JpmcdsCdsFeeLegPV(today,
                          valueDate,
                          stepinDate,
                          startDate,
                          endDate,
                          payAccOnDefault,
                          dateInterval,
                          stubType,
                          1.0, /* notional */
                          couponRate,
                          paymentDcc,
                          badDayConv,
                          calendar,
                          discCurve,
                          spreadCurve,
                          protectStart,
                          isPriceClean,
                          &feeLegPV) != SUCCESS)
        goto done;

    if (MAX(stepinDate, startDate) <= endDate) 
    {
        if (JpmcdsCdsContingentLegPV(today,
                                     valueDate,
                                     MAX(stepinDate, startDate),
                                     endDate,
                                     1.0, /* notional */
                                     discCurve,
                                     spreadCurve,
                                     recoveryRate,
                                     protectStart,
                                     &contingentLegPV) != SUCCESS)
            goto done;
    }

    *price = contingentLegPV - feeLegPV;
    status = SUCCESS;

 done:

    if (status != SUCCESS)
        JpmcdsErrMsgFailure (routine);

    return status;
}


/*
***************************************************************************
** Computes the par spread for a vanilla CDS which produces a zero price.
***************************************************************************
*/
EXPORT int JpmcdsCdsParSpreads(
    TDate           today,
    TDate           stepinDate,
    TDate           startDate,
    long            nbEndDates,
    TDate          *endDates,
    TBoolean        payAccOnDefault,
    TDateInterval  *couponInterval,
    TStubMethod    *stubType,
    long            paymentDcc,
    long            badDayConv,
    char           *calendar,
    TCurve         *discCurve,
    TCurve         *spreadCurve,
    double          recoveryRate,
    double         *parSpread)
{
    static char routine[] = "JpmcdsCdsParSpread";
    int         status    = FAILURE;

    double      feeLegPV;
    double      contingentLegPV;
    long        i;
    TBoolean    isPriceClean = 1;
    TBoolean    protectStart = TRUE;

    REQUIRE(parSpread != NULL);
    REQUIRE(nbEndDates >= 1);
    REQUIRE(stepinDate >= today);
    /* all other requirements can be handled by the routines we call */

    for(i = 0; i < nbEndDates; ++i)
    {
        if (JpmcdsCdsFeeLegPV(today,
                              stepinDate, /* valueDate */
                              stepinDate,
                              startDate,
                              endDates[i],
                              payAccOnDefault,
                              couponInterval,
                              stubType,
                              1.0, /* notional */
                              1.0, /* couponRate */
                              paymentDcc,
                              badDayConv,
                              calendar,
                              discCurve,
                              spreadCurve,
                              protectStart,
                              isPriceClean,
                              &feeLegPV) != SUCCESS)
            goto done;

        if (JpmcdsCdsContingentLegPV(today,
                                     stepinDate,
                                     MAX(stepinDate, startDate),
                                     endDates[i],
                                     1.0, /* notional */
                                     discCurve,
                                     spreadCurve,
                                     recoveryRate,
                                     protectStart,
                                     &contingentLegPV) != SUCCESS)
            goto done;
        
        parSpread[i] = contingentLegPV / feeLegPV;
    }

    status = SUCCESS;

 done:

    if (status != SUCCESS)
        JpmcdsErrMsgFailure (routine);

    return status;
}


/*
***************************************************************************
** Computes the non-contingent cash flows for a fee leg. These are the 
** cash flows you will receive if there is no default.
***************************************************************************
*/
EXPORT TCashFlowList *JpmcdsCdsFeeLegFlows
(TDate           startDate,
 TDate           endDate,
 TDateInterval  *dateInterval,
 TStubMethod    *stubType,
 double          notional,
 double          couponRate,
 long            paymentDcc,
 long            badDayConv,
 char           *calendar)
{
    static char routine[] = "JpmcdsCdsFeeLegFlows";

    TCashFlowList *cfl = NULL;
    TFeeLeg       *fl  = NULL;
    TBoolean       protectStart = TRUE;
    TBoolean       payAccOnDefault = 1;
     
    fl = JpmcdsCdsFeeLegMake(startDate,
                             endDate,
                             payAccOnDefault,
                             dateInterval,
                             stubType,
                             notional,
                             couponRate,
                             paymentDcc,
                             badDayConv,
                             calendar,
                             protectStart);

    if (fl == NULL)
        goto done;

    cfl = JpmcdsFeeLegFlows (fl);

 done:

    JpmcdsFeeLegFree (fl);
    if (cfl == NULL)
        JpmcdsErrMsgFailure (routine);

    return cfl;
}
