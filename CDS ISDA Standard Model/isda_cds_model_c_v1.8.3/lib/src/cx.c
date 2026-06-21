/*
 * ISDA CDS Standard Model
 *
 * Copyright (C) 2009 International Swaps and Derivatives Association, Inc.
 * Developed and supported in collaboration with Markit
 * 
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the ISDA CDS Standard Model Public License.
 */

#include "cx.h"
#include <ctype.h>
#include "convert.h"
#include "tcurve.h"
#include "macros.h"
#include "cerror.h"


/*
***************************************************************************
** Constructor for TContingentLeg
***************************************************************************
*/
TContingentLeg* JpmcdsContingentLegMake(
TDate           startDate,           /* (I) */
TDate           endDate,             /* (I) */
double          notional,            /* (I) */
TProtPayConv    payType,             /* (I) */
TBoolean        protectStart         /* (I) */
)
{
    static char routine[] = "JpmcdsContingentLegMake";
    int status = FAILURE;

    TContingentLeg* p = NULL;

    p = NEW(TContingentLeg);
    if (p==NULL) goto done;

    /* p->startDate is defined as giving protection from end of startDate.
       So if we want to protect on the start date, we need to move this
       date forward by one. */
    if (protectStart)
    {
        p->startDate = startDate-1;
    }
    else
    {
        p->startDate = startDate;
    }

    p->protectStart = protectStart;
    p->endDate    = endDate;
    p->notional   = notional;
    p->payType    = payType;

    REQUIRE (p->endDate > p->startDate);

    status = SUCCESS;

done:
    if (status != SUCCESS)
    {
        JpmcdsErrMsg("%s: Failed!\n", routine);
        FREE(p);
        p=NULL;
    }

    return p;
}

/*
***************************************************************************
** Constructor for TFeeLeg
***************************************************************************
*/
TFeeLeg* JpmcdsFeeLegMake(
int             nbDates,             /* (I) */
TDate*          accStartDates,       /* (I) [nbDates] */
TDate*          accEndDates,         /* (I) [nbDates] */
TDate*          payDates,            /* (I) [nbDates] */
double          notional,           /* (I) */
double          couponRate,         /* (I) */
long            dcc,                 /* (I) */
TAccrualPayConv accrualPayConv,      /* (I) */
TBoolean        obsStartOfDay        /* (I) */
)
{
    static char routine[] = "JpmcdsFeeLegMake";
    int status = FAILURE;

    TFeeLeg* p = NULL;


    p = NEW(TFeeLeg);
    if (p==NULL) goto done;

    p->nbDates         = nbDates;
    if (p->nbDates > 0 && accStartDates != NULL)
    {
        p->accStartDates = NEW_ARRAY(TDate, p->nbDates);
        if (p->accStartDates == NULL) goto done;
        COPY_ARRAY (p->accStartDates, accStartDates, TDate, p->nbDates);
    }
    else
    {
        p->accStartDates = NULL;
    }

    if (p->nbDates > 0 && accEndDates != NULL)
    {
        p->accEndDates = NEW_ARRAY(TDate, p->nbDates);
        if (p->accEndDates == NULL) goto done;
        COPY_ARRAY (p->accEndDates, accEndDates, TDate, p->nbDates);
    }
    else
    {
        p->accEndDates = NULL;
    }

    if (p->nbDates > 0 && payDates != NULL)
    {
        p->payDates = NEW_ARRAY(TDate, p->nbDates);
        if (p->payDates == NULL) goto done;
        COPY_ARRAY (p->payDates, payDates, TDate, p->nbDates);
    }
    else
    {
        p->payDates = NULL;
    }

    p->notional        = notional;
    p->couponRate      = couponRate;
    p->dcc             = dcc;
    p->accrualPayConv  = accrualPayConv;
    p->obsStartOfDay   = obsStartOfDay;

    status = SUCCESS;

done:
    if (status != SUCCESS)
    {
        JpmcdsErrMsg("%s: Failed!\n", routine);
        JpmcdsFeeLegFree(p);
        p=NULL;
    }

    return p;
}


/*
***************************************************************************
** Memory allocator for TFeeLeg
***************************************************************************
*/
TFeeLeg* JpmcdsFeeLegMakeEmpty(
int             nbDates              /* (I) */
)
{
    static char routine[] = "JpmcdsFeeLegMakeEmpty";
    int status = FAILURE;

    TFeeLeg* p = NULL;


    p = NEW(TFeeLeg);
    if (p==NULL) goto done;

    p->nbDates         = nbDates;

    if (p->nbDates > 0)
    {
        p->accStartDates = NEW_ARRAY(TDate, p->nbDates);
        if (p->accStartDates == NULL) goto done;
    }
    else
    {
        p->accStartDates = NULL;
    }

    if (p->nbDates > 0)
    {
        p->accEndDates = NEW_ARRAY(TDate, p->nbDates);
        if (p->accEndDates == NULL) goto done;
    }
    else
    {
        p->accEndDates = NULL;
    }

    if (p->nbDates > 0)
    {
        p->payDates = NEW_ARRAY(TDate, p->nbDates);
        if (p->payDates == NULL) goto done;
    }
    else
    {
        p->payDates = NULL;
    }

    status = SUCCESS;

done:
    if (status != SUCCESS)
    {
        JpmcdsErrMsg("%s: Failed!\n", routine);
        JpmcdsFeeLegFree(p);
        p=NULL;
    }

    return p;
}


/*
***************************************************************************
** Destructor for TFeeLeg
***************************************************************************
*/
void JpmcdsFeeLegFree(TFeeLeg *p)
{
    if (p != NULL)
    {
        FREE(p->accStartDates);
        FREE(p->accEndDates);
        FREE(p->payDates);
        FREE(p);
    }
}
