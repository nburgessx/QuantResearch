/*
 * ISDA CDS Standard Model
 *
 * Copyright (C) 2009 International Swaps and Derivatives Association, Inc.
 * Developed and supported in collaboration with Markit
 * 
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the ISDA CDS Standard Model Public License.
 */

#include "cashflow.h"
#include <limits.h>
#include <string.h>
#include <math.h>
#include "cheaders.h"
#include "bastypes.h"
#include "bsearch.h"
#include "convert.h"
#include "datelist.h"
#include "stub.h"
#include "ldate.h"
#include "busday.h"
#include "macros.h"
#include "tcurve.h"
#include "date_sup.h"


/*
***************************************************************************
** Merges two cash flow lists.
***************************************************************************
*/
TCashFlowList* JpmcdsMergeCFL(
   TCashFlowList *a,                       /* (I) first cash flow list */
   TCashFlowList *b)                       /* (I) second cash flow list */
{
   static char routine[] = "JpmcdsMergeCFL";
   TCashFlowList *cfl = NULL;              /* cash flow list to return */
   int            i;                       /* where to add next cf in cfl */
   int            ia;                      /* iterates over cfs in a */
   int            ib;                      /* iterates over cfs in b */

   if (a == NULL || b == NULL)
   {
       JpmcdsErrMsg ("%s: NULL inputs\n", routine);
       goto done;
   }

   cfl = JpmcdsNewEmptyCFL(a->fNumItems + b->fNumItems);
   if (cfl == NULL)
       goto done;

   i = 0;
   ia = 0;
   ib = 0;
   while (ia<a->fNumItems && ib<b->fNumItems)
   {
       TDate aDate = a->fArray[ia].fDate;
       TDate bDate = b->fArray[ib].fDate;

       if (aDate < bDate)
       {
           cfl->fArray[i++] = a->fArray[ia++];
       }
       else if (aDate > bDate)
       {
           cfl->fArray[i++] = b->fArray[ib++];
       }
       else
       {
           cfl->fArray[i].fDate = aDate;
           cfl->fArray[i++].fAmount = a->fArray[ia++].fAmount + 
                                      b->fArray[ib++].fAmount;
       }
   }

   while (ia<a->fNumItems)             /* add any from A left over */
   {
       cfl->fArray[i++] = a->fArray[ia++];
   }

   while (ib<b->fNumItems)             /* add any from B left over */
   {
       cfl->fArray[i++] = b->fArray[ib++];
   }

   cfl->fNumItems = i;

done:

   if (cfl == NULL)
       JpmcdsErrMsg("%s: FAILED\n",routine);

   return cfl;
}


/*
***************************************************************************
** Makes a new empty TCashFlowList.
***************************************************************************
*/
TCashFlowList* JpmcdsNewEmptyCFL(int numItems)
{
    static char    routine[]="JpmcdsNewEmptyCFL";
    TCashFlowList *cfl = NULL;
    
    if (numItems < 0)
    {
        JpmcdsErrMsg("%s: # cashflows (%d) must be >= 0.\n", routine, numItems);
        goto done;
    }
    
    cfl = NEW(TCashFlowList);
    if (cfl == NULL)
        goto done;
    
    cfl->fNumItems = numItems;
    if (numItems > 0)
    {
        cfl->fArray = NEW_ARRAY (TCashFlow, numItems);
        if (cfl->fArray == NULL)
            goto done;
    }
    else
    {
        cfl->fArray = NULL;
    }
    /* Success
     */
    return cfl;


  done:
    JpmcdsErrMsg("%s: Failed.\n", routine);
    JpmcdsFreeCFL(cfl);
    return NULL;
}

/*
***************************************************************************
** Makes a new TCashFlowList.
***************************************************************************
*/
TCashFlowList* JpmcdsNewCFL(TDate  *dates,    /* (I) Dates */
                         double *amounts,  /* (I) Amounts */
                         int     numItems) /* (I) Length of arrays */
{
    int            idx;
    TCashFlowList *cfl;

    cfl = JpmcdsNewEmptyCFL(numItems);
    if( cfl == NULL )
        return NULL;

    for (idx = 0; idx < numItems; idx++)
    {
        cfl->fArray[idx].fDate = dates[idx];
        cfl->fArray[idx].fAmount = amounts[idx];
    }
    return cfl;
}


/*
***************************************************************************
** Makes a new  cash flow list.
**
** Depending on the value of flags, an initial paymemt of value 1 is 
** subtracted and/or  a final payment of 1 is added. If 
** JPMCDS_PRESTART_ZERO_PAYMENT is set, and there is a front stub, the first 
** date will be the payment date that would have occured before the start 
** date. In this case, JPMCDS_SUBTRACT_INITIAL cannot be set.
***************************************************************************
*/
TCashFlowList* JpmcdsMakeCFL
    (double couponRate,                 /* (I) Annual coupon Rate  */
     TDate startDate,                   /* (I) Date instrument begins at */
     TDateInterval *interval,           /* (I) Time between payments */
     TDate maturityDate,                /* (I) Date instrument matures at*/
     long couponDayCountConv,           /* (I) For coupon calculation */
     long stubPlacement,                /* (I) [Front / Back] + [Short (default) / Long] */
     long flags,                        /* (I) What type (See cashflow.h) */
     long accrualBadDayConv,            /* (I) Accrual bad day convention */
     long payBadDayConv,                /* (I) Pay bad day convention */
     char *holidayFile)                 /* (I) Holiday file */
{
    static char routine[]="JpmcdsMakeCFLRoll";

    TDateList *dl = NULL;
    TCashFlowList *cfl = NULL;
    
    if (maturityDate <= startDate)
    {
        JpmcdsErrMsg("%s: maturityDate(%s) <= startDate(%s)\n",
            routine, JpmcdsFormatDate(maturityDate), JpmcdsFormatDate(startDate));
        return NULL;
    }

    if ((flags & JPMCDS_SUBTRACT_INITIAL || flags & JPMCDS_KEEP_START_DATE) && flags & JPMCDS_PRESTART_ZERO_PAYMENT)
    {
        JpmcdsErrMsg("%s: cannot subtract initial payment at prestart date.\n", routine);
        JpmcdsErrMsg("\tAlso cannot keep both the start and prestart dates.\n");
        return NULL;        
    }

    if (flags & JPMCDS_POST_MATURITY_ZERO_PAYMENT)
    {
        if (stubPlacement & JPMCDS_STUB_POSN_MASK)
        {
            JpmcdsErrMsg("%s: no post maturity zero payment w/ stub at end.\n", routine);
            return NULL;        
        }
        else if (flags & JPMCDS_ADD_FINAL)
        {
            JpmcdsErrMsg("%s: no post maturity zero payment w/ final payment.\n", routine);
            return NULL;        
        }
        else                            /* Advance maturityDate by interval */
        {
            if (JpmcdsDtFwdAny(maturityDate, interval, &maturityDate) == FAILURE)
                goto done;
        }
    }
    
   /* Set up the TDateList of coupon dates; if stub at end, last date
    * could be later than maturityDate. If there is a stub at the start, 1st date
    * could be earlier than startDate.
    */
    dl = JpmcdsNewDateListExtendedRoll(startDate,
                                maturityDate,
                                startDate /* roll date */,
                                interval,
                                stubPlacement & JPMCDS_STUB_POSN_MASK);
    if (dl == (TDateList *)NULL)
        goto done;

    /* Make cash flow list based on an annual rate of couponRate.
    * The startDate is included, with a cashflow of 0. 
    */
    cfl = JpmcdsNewCFLFromDL3(dl, startDate, maturityDate, couponRate,
                           couponDayCountConv, stubPlacement & JPMCDS_STUB_LENGTH_LONG,
                           accrualBadDayConv, payBadDayConv, holidayFile, holidayFile,
                           (TBoolean)(flags & JPMCDS_PRESTART_ZERO_PAYMENT));

    JpmcdsFreeDateList(dl);                /* Free memory */
                          
    if (cfl == NULL)
        goto done;
    if( cfl->fNumItems == 0 )
        return cfl;

    /* Initial payment of principal */
    if (flags & JPMCDS_SUBTRACT_INITIAL)
    {
        cfl->fArray[0].fAmount-= 1.;
    }
    
    /* Add on the final repayment of principal */
    if (flags & JPMCDS_ADD_FINAL)
    {
        cfl->fArray[cfl->fNumItems-1].fAmount += 1.;
    }

    if (flags & JPMCDS_POST_MATURITY_ZERO_PAYMENT)
    {
        cfl->fArray[cfl->fNumItems-1].fAmount = 0.;
    }
    
    /* Dont want zero-cashflows around, so move
     * all cashflows back by one, unless JPMCDS_PRESTART_ZERO_PAYMENT
     * or JPMCDS_KEEP_START_DATE is set. Testing equality to zero
     * is no problem since no calculations were done.
     */
    if (!(flags & JPMCDS_PRESTART_ZERO_PAYMENT || flags & JPMCDS_KEEP_START_DATE) && cfl->fArray[0].fAmount == 0.)
    {
        int idx;
        for (idx=0; idx < cfl->fNumItems-1; idx++)
        {
            cfl->fArray[idx] = cfl->fArray[idx+1];
        }
        cfl->fNumItems--;
    }

done:
    if (cfl == NULL)
        JpmcdsErrMsg("%s: Failed.\n", routine);
    
    return cfl;
}


/*
***************************************************************************
** Makes a cashflow list from a DateList. Cashflows (for non stub payments) 
** are equal to the rate times the fraction of the year between  each date 
** and the previous date. The year fraction is calculated by the year 
** fraction function passed in. For stub payments, Actual/Actual is used.
**
** If there is a stub at the BEGINNING, then the first element of the 
** dateList is assumed to PRECEED startDate.
**
** If there is a stub at the END, then the last element of the dateList is 
** assumed to COME AFTER maturityDate.
***************************************************************************
*/
TCashFlowList* JpmcdsNewCFLFromDL3
    (TDateList *dates,                  /* (I) Dates */
     TDate startDate,                   /* (I) startDate */
     TDate maturityDate,                /* (I) MaturityDate */
     double rate,                       /* (I) Multiplied by yearFrac */
     long couponDayCountConv,           /* (I) Coupon calc  */
     TBoolean longStub,                 /* (I) Long stub */
     long accrualBadDayConv,            /* (I) Accrual bad day convention */
     long payBadDayConv,                /* (I) Pay bad day convention */
     char *accrualHolidayFile,          /* (I) Accrual Holiday file. */
     char *payHolidayFile,              /* (I) Payment Holiday file. */ 
     TBoolean usePreStartDate)          /* (I) Use prestart date if stub */
{
    static char routine[]="JpmcdsNewCFLFromDL3";
    int idx;
    int startIdx;                       /* First index for regular proc */
    int endIdx;                         /* Last index for regular proc */   
    int numDates = dates->fNumItems;
    TDate accrualStart, payStart;       /* Accrue and pay start dates */
    TDate accrualMaturity, payMaturity; /* Accrue and pay maturity dates */    
    TDateList *accrualDates = NULL;
    TDateList *payDates = NULL;
    TBoolean frontStub, backStub;
    TDateList *datesCopy = NULL;
    TCashFlowList *cfl = NULL;
    int status = FAILURE;

    if ( numDates == 0 )
    {
        cfl = JpmcdsNewEmptyCFL(0);
        if ( cfl == NULL )
            goto done;

        status = SUCCESS;
        goto done;
    }

    frontStub = ( startDate > dates->fArray[0] );
    backStub = ( maturityDate < dates->fArray[numDates-1] );

    /* For long stubs, simply coalesce first/last two coupon periods
     * (i.e. short stub + adjacent regular period) into one...
     */
    if ( longStub && (frontStub || backStub) && numDates > 2 )
    {
        datesCopy = JpmcdsNewEmptyDateList(numDates - 1);
        if ( datesCopy == NULL )
            goto done;

        if ( frontStub )
        {
            datesCopy->fArray[0] = dates->fArray[0];
            memcpy(&datesCopy->fArray[1], &dates->fArray[2], (numDates - 2) * sizeof(TDate));
            dates = datesCopy; /* simply overwrite pointer, hiding original */
            --numDates;
        }

        if (backStub)
        {
            memcpy(&datesCopy->fArray[0], &dates->fArray[0], (numDates - 2) * sizeof(TDate));
            datesCopy->fArray[numDates - 2] = dates->fArray[numDates - 1];
            dates = datesCopy; /* simply overwrite pointer, hiding original */
            --numDates;
        }
    }

    cfl = JpmcdsNewEmptyCFL(numDates);
    if( cfl == NULL )
        goto done;

    /* Get accrualDates
     */
    if (accrualBadDayConv == JPMCDS_BAD_DAY_NONE)
    {
        accrualDates = dates;            /* efficiency */
        accrualStart = startDate;
        accrualMaturity = maturityDate;
    }
    else
    {
        /* Adjust for bad days */
        accrualDates = JpmcdsNewDateListBusDayAdj(dates, accrualBadDayConv,
                                               accrualHolidayFile);
        if (accrualDates == NULL)
            goto done;
        
        if (JpmcdsBusinessDay(startDate, accrualBadDayConv, accrualHolidayFile,
                           &accrualStart) == FAILURE)
            goto done;

        if (JpmcdsBusinessDay(maturityDate, accrualBadDayConv, accrualHolidayFile,
                           &accrualMaturity) == FAILURE)
            goto done;
    }

    /* Get payDates
     */
    if (payBadDayConv == JPMCDS_BAD_DAY_NONE)
    {
        payDates = dates;            /* efficiency */
        payStart = startDate;
        payMaturity = maturityDate;
    }
    else
    {
        /* Adjust for bad days */
        payDates = JpmcdsNewDateListBusDayAdj(dates, payBadDayConv, payHolidayFile);
        if (payDates == NULL)
            goto done;

        if (JpmcdsBusinessDay(startDate, payBadDayConv, payHolidayFile,
                           &payStart) == FAILURE)
            goto done;

        if (JpmcdsBusinessDay(maturityDate, payBadDayConv, payHolidayFile,
                           &payMaturity) == FAILURE)
            goto done;
    }

    if (usePreStartDate)                /* Useful to JpmcdsTrinomialModel */
        cfl->fArray[0].fDate = payDates->fArray[0];
    else      
        cfl->fArray[0].fDate = payStart;
    
    cfl->fArray[0].fAmount = 0.;
    cfl->fArray[numDates-1].fDate = payMaturity;

    for (idx=1; idx < numDates-1; idx++)
        cfl->fArray[idx].fDate = payDates->fArray[idx];
    
    startIdx = 1;
    endIdx = numDates-1;
    if ( frontStub )
    {
        startIdx = 2;                   /* For non-stub processing */
        if (JpmcdsStubPayment(accrualDates->fArray[0], /*(I)Coupon start*/
                           accrualDates->fArray[1], /*(I)Coupon end */
                           accrualStart,            /*(I)Accrue start*/
                           accrualDates->fArray[1], /*(I)Accrue End */
                           rate,                            /*(I) */
                           couponDayCountConv,              /*(I)  */
                           &cfl->fArray[1].fAmount)         /*(O) */
            == FAILURE)
            goto done;
    }
    
    if ( backStub )
    {
        endIdx = numDates-2;            /* For non-stub processing */
        if (JpmcdsStubPayment(accrualDates->fArray[numDates-2], /*(I)Coup start*/
                           accrualDates->fArray[numDates-1], /*(I)Coupon end*/
                           accrualDates->fArray[numDates-2], /*(I)Accru start*/
                           accrualMaturity,                  /*(I)Accrue end*/
                           rate,                             /*(I)*/
                           couponDayCountConv,               /*(I)*/
                           &cfl->fArray[numDates-1].fAmount) /*(O)*/
                           == FAILURE)
            goto done;
    }

    
    /* Compute non-stub cashflows */
    for (idx=startIdx; idx <= endIdx; idx++)
    {
        double yearFrac;

        if (JpmcdsDayCountFraction(accrualDates->fArray[idx-1],
                                accrualDates->fArray[idx],
                                couponDayCountConv, &yearFrac) == FAILURE)
            goto done;

        cfl->fArray[idx].fAmount = rate * yearFrac;
    }

    status = SUCCESS;

done:

    JpmcdsFreeDateList(datesCopy);

    if (payDates != dates)
        JpmcdsFreeDateList(payDates);
    
    if (accrualDates != dates)
        JpmcdsFreeDateList(accrualDates);

    if (status == FAILURE)
    {
        JpmcdsFreeCFL(cfl);
        return NULL;
    }
    return cfl;
}


/*
***************************************************************************
** Free cash flow list.
***************************************************************************
*/
void JpmcdsFreeCFL(TCashFlowList *theCashFlowList)
{
    if (theCashFlowList != (TCashFlowList *)NULL)
    {
        FREE(theCashFlowList->fArray);
        FREE(theCashFlowList);
    }
}
