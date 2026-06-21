/*
 * ISDA CDS Standard Model
 *
 * Copyright (C) 2009 International Swaps and Derivatives Association, Inc.
 * Developed and supported in collaboration with Markit
 * 
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the ISDA CDS Standard Model Public License.
 */

#ifndef CASHFLOW_H
#define CASHFLOW_H

#include "bastypes.h"   
#include "stub.h"
#include "busday.h"
#include "ldate.h"

#ifdef __cplusplus
extern "C"
{
#endif


/*f
***************************************************************************
** Makes a new empty TCashFlowList.
***************************************************************************
*/
TCashFlowList* JpmcdsNewEmptyCFL(int numItems);


/*f
***************************************************************************
** Makes a cashflow list from a DateList. Cashflows (for non stub payments) 
** are equal to the rate times the fraction of the year between each date 
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
     TBoolean usePreStartDate);         /* (I) Use prestart date if stub */

 
/*f
***************************************************************************
** Makes a new TCashFlowList.
***************************************************************************
*/
TCashFlowList* JpmcdsNewCFL
    (TDate  *dates,      /* (I) Dates  */
     double *amounts,    /* (I) Amounts */
     int     numItems);  /* (I) Length of both arrays */


#define JPMCDS_ADD_FINAL 1                 /* Add 1. on maturityDate */
#define JPMCDS_SUBTRACT_INITIAL 2          /* Subtract 1. at startDate */
#define JPMCDS_PRESTART_ZERO_PAYMENT 4     /* If stub, leave prestart  */
                                        /* zero payment. */
#define JPMCDS_POST_MATURITY_ZERO_PAYMENT 8 /* Add final zero payment  */
                                        /* Used for libor-in arrears */
#define JPMCDS_KEEP_START_DATE 16          /* Keep start date, even if  */
                                        /* JPMCDS_SUBTRACT_INITIAL not set */


/*f
***************************************************************************
** Makes a new  cash flow list.
**
** Depending on the value of flags, an initial paymemt of value 1 is 
** subtracted and/or a final payment of 1 is added. If 
** JPMCDS_PRESTART_ZERO_PAYMENT is set, and there is a front stub, the first 
** date will be the payment date that would have occured before the start 
**  date. In this case, JPMCDS_SUBTRACT_INITIAL cannot be set.
***************************************************************************
*/
TCashFlowList *  JpmcdsMakeCFL
    (double         couponRate,     /* (I) Annual coupon Rate  */
     TDate          startDate,      /* (I) Date instrument begins at */
     TDateInterval *interval,       /* (I) Time between payments */
     TDate          maturityDate,   /* (I) Date instrument matures at*/
     long           dayCountConv,   /* (I) See JpmcdsDayCountFraction */
     long           stubPlacement,  /* (I) [Front / Back] + [Short (default) / Long] */
     long           flags,          /* (I) JPMCDS_ADD_FINAL/JPMCDS_SUBTRACT_INITIAL*/
     long           accBadDayConv,  /* (I) Accrual conv (See JpmcdsBusinessDay) */
     long           payBadDayConv,  /* (I) Payment conv (See JpmcdsBusinessDay) */
     char          *holidayFile);   /* (I) See JpmcdsBusinessDay */


/*f
***************************************************************************
** Merges two cash flow lists. The input cash flow lists cannot be NULL.
***************************************************************************
*/
TCashFlowList *  JpmcdsMergeCFL(
   TCashFlowList *a,                       /* (I) first cash flow list */
   TCashFlowList *b);                      /* (I) second cash flow list */


/*f
***************************************************************************
** Free cash flow list.
***************************************************************************
*/
void  JpmcdsFreeCFL
      (TCashFlowList *);                /* (I) */

#ifdef __cplusplus
}
#endif

#endif    /* CASHFLOW_H */
