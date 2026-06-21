/*
 * ISDA CDS Standard Model
 *
 * Copyright (C) 2009 International Swaps and Derivatives Association, Inc.
 * Developed and supported in collaboration with Markit
 * 
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the ISDA CDS Standard Model Public License.
 */

#include <math.h>
#include "cgeneral.h"
#include "badday.h"
#include "cashflow.h"
#include "cerror.h"
#include "cmemory.h"
#include "convert.h"
#include "datelist.h"
#include "date_sup.h"
#include "ldate.h" 
#include "macros.h"
#include "tcurve.h"
#include "zcprvt.h"
#include "zcswdate.h"



/*
***************************************************************************
** Makes a cash flow list for a swap instrument.
**
** return s cash flow list associated with the swap, NULL if an error found.
***************************************************************************
*/
TCashFlowList* JpmcdsZCGetSwapCFL
  (TDate           valueDate,           /* (I) Value date                  */
   TDate           matDate,             /* (I) Unadjusted maturity date    */
   TBoolean        stubAtEnd,           /* (I) If matDate on-cycle         */
   double          rate,                /* (I) Coupon rate                 */
   TDateInterval  *interval,            /* (I) Coupon payment interval     */
   long            dayCountConv,        /* (I) Day count convention        */
   TBadDayList    *badDayList,          /* (I) Bad day adjustment list     */
   long            badDayConv,          /* (I) Bad day convention for rate */
   char           *holidayFile)         /* (I) Name of holiday file        */
{
    static char routine[]="JpmcdsZCGetSwapCFL";
    int status = FAILURE;               /* Until proven successful */

    TCashFlowList  *cfl = NULL;         /* cash-flow-list to return */
    TDateList      *dl = NULL;          /* list of coupon dates */
    TDate           prevDate;           /* prev date added to cash-flow list */
    int             i;                  /* loops over coupon dates */

    /* Bizarre rate==0.0 case 
     */
    if (rate==0.0)                      
    {
        cfl = JpmcdsNewEmptyCFL(1);        /* simply cash flow of 1.0 at date */
        if (cfl == NULL)
            goto done;
        cfl->fArray[0].fAmount = 1;
        if (JpmcdsZCAdjustDate(matDate, badDayList, badDayConv, holidayFile,
                            &cfl->fArray[0].fDate) == FAILURE)
            goto done;
        
        return cfl;
    }
    
     /* Get dateList with adjusted coupon dates & mat date 
      */
    dl = JpmcdsZCGetSwapCouponDL(valueDate,   
                              matDate, stubAtEnd, interval,
                              badDayList, badDayConv, holidayFile);
   if (dl == NULL)
       goto done;

   cfl = JpmcdsNewEmptyCFL(dl->fNumItems);    /* make cash flow list w/ coupons */
   if (cfl == NULL)
       goto done;


   prevDate = valueDate;
   for (i=0; i<dl->fNumItems; i++)
   {
       TDate  cDate = dl->fArray[i];       /* coupon date */
       double yearFraction;                /* fraction of year for coupon #i */
       if (JpmcdsDayCountFraction
           (prevDate, cDate, dayCountConv, &yearFraction) == FAILURE)
           goto done;

       cfl->fArray[i].fAmount = rate * yearFraction;
       cfl->fArray[i].fDate = cDate;      /* store date */
       prevDate = cDate;
   }

   cfl->fArray[cfl->fNumItems-1].fAmount += 1.0;   /* add principal */

   status = SUCCESS;


 done:
   JpmcdsFreeDateList(dl);                /* Free date list memory */
   if (status == FAILURE)
   {
       JpmcdsErrMsg("%s: Failed.\n", routine);
       JpmcdsFreeCFL(cfl);
       cfl = NULL;
   }
   return cfl;
}


/*
***************************************************************************
** Makes a date list for all coupons associated w/ a swap instrument.
**
** Only glitch is possible inclusion of a stub date, which is necessary if the
** maturity date isn't an integral number of frequency intervals away, e.g.
** a swap date 5 years and 1 month from the value date, which would have a
** stub date 1 month from now, followed by coupons every year from then.
**
** Returns date list associated with the swap's coupons, NULL if error found.
***************************************************************************
*/
TDateList *JpmcdsZCGetSwapCouponDL
  (TDate           valueDate,           /* (I) Value date */
   TDate           matDate,             /* (I) Unadjusted maturity date */
   TBoolean        stubAtEnd,           /* (I) If matDate onCycle from vdate */
   TDateInterval  *interval,            /* (I) Coupon payment interval */
   TBadDayList    *badDayList,          /* (I) Bad day adjustment list */
   long            badDayConv,          /* (I) */
   char           *holidayFile)         /* (I) */
{
    static char routine[] = "JpmcdsZCGetSwapCouponDL";
    int             status = FAILURE;   /* Until proven successful */
    TDateList      *dl = NULL;          /* List of coupon dates */

    /* If the matDate is on-cycle, then stub is at end, because we
     * are counting forward from the value date.
     */
    dl = JpmcdsNewPayDates(valueDate, matDate, interval, stubAtEnd);
    if (dl == NULL)
        goto done;


    /* Now adjust for bad days.
     */
    if (badDayConv != JPMCDS_BAD_DAY_NONE)
    {
        if (JpmcdsDateListBusDayAdj(dl,badDayConv,holidayFile)==FAILURE)
            goto done;
    }

    JpmcdsFixBadDLBadDayList(badDayList,dl);
    status = SUCCESS;

done:
    if (status == FAILURE)
    {
        JpmcdsErrMsg("%s: Failed.\n", routine);
        JpmcdsFreeDateList(dl);
        dl = NULL;
    }

    return dl;
}



