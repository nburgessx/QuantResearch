/*
 * ISDA CDS Standard Model
 *
 * Copyright (C) 2009 International Swaps and Derivatives Association, Inc.
 * Developed and supported in collaboration with Markit
 * 
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the ISDA CDS Standard Model Public License.
 */

#include "defaulted.h"
#include "macros.h"
#include "cds.h"
#include "cerror.h"

EXPORT int JpmcdsDefaultAccrual(
        TDate           tradeDate, 
        TDate           edd, 
        TDate           startDate, 
        TDate           endDate, 
        TDateInterval  *couponInterval,
        TStubMethod    *stubType,
        double          notional,
        double          couponRate,
        long            paymentDcc,
        long            badDayConv,
        char           *calendar,
        double         *accrualDays,
        double         *defaultAccrual)
{
    static char         routine[]       = "JpmcdsDefaultedCDS";
    int                 status          = FAILURE;
    TFeeLeg            *fl              = NULL;

    int i;
    TDate accrualStartDateAdj;
    TDate accrualEndDateAdj;
    
    *accrualDays = 0.;
    *defaultAccrual = 0.;

    if (tradeDate < edd)
    {
        status = SUCCESS;
        goto done;
    }
    
    if (edd < startDate)
        goto done;
    
    fl = JpmcdsCdsFeeLegMake(
            startDate,
            endDate,
            TRUE,
            couponInterval,
            stubType,
            notional,
            couponRate,
            paymentDcc,
            badDayConv,
            calendar,
            TRUE);
    
    if (NULL == fl) 
        goto done;

    i = 0;
    while (i < fl->nbDates)
    {
        if (fl->accStartDates[i] <= tradeDate && tradeDate < fl->accEndDates[i])
        {
            *accrualDays = edd + 1.0 - fl->accStartDates[i];
            if(JpmcdsDayCountFraction(
                    fl->accStartDates[i], 
                    edd + 1, 
                    paymentDcc, 
                    defaultAccrual) != SUCCESS)
                goto done;
            *defaultAccrual *= couponRate * notional;
            status = SUCCESS;
            goto done;
        }
        ++i;
    }
    
    done:
        if (status != SUCCESS)
        {
            JpmcdsErrMsgFailure(routine);
        }

        JpmcdsFeeLegFree(fl);

        return status;
}
