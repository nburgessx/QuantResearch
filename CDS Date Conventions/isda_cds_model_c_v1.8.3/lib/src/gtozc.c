/*
 * ISDA CDS Standard Model
 *
 * Copyright (C) 2009 International Swaps and Derivatives Association, Inc.
 * Developed and supported in collaboration with Markit
 * 
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the ISDA CDS Standard Model Public License.
 */

#include <stdlib.h>
#include "bastypes.h"
#include "cgeneral.h"
#include "cerror.h"
#include "macros.h"
#include "convert.h"
#include "ldate.h"
#include "interp.h"
#include "gtozc.h"
#include "tcurve.h"
#include "zcprvt.h"
#include "yearfrac.h"

#define ZC_DEFAULT_BASIS 1
#define ZC_DEFAULT_DAYCNT JPMCDS_ACT_365F

/*
***************************************************************************
** Check input arguments for adding cash to a given zero curve.
***************************************************************************
*/
static int CheckZCCashInputs(
        TCurve    *zeroCurve,
        TDate     *dateArr,
        double    *rateArr,
        int       numRates,
        long      dayCountConv);


/*
***************************************************************************
** Check input arguments for adding swaps to a given zero curve.
***************************************************************************
*/
static int CheckZCSwapsInputs(
        TCurve     *zeroCurve,
        TCurve     *discZC,
        TDate      *dates,
        double     *rates,
        int        numSwaps,
        int        fixedSwapFreq,
        int        floatSwapFreq,
        long       fixDayCountConv,
        long       floatDayCountConv,
        char       fwdLength,
        long       badDayConv,
        char       *holidayFile);


/*
***************************************************************************
** Adds cash points to a given zero curve.
***************************************************************************
*/
TCurve* JpmcdsZCCash(
   TCurve        *zeroCurve,         /* (I) Zero curve to add to             */
   TDate         *dates,             /* (I) Array of cash dates              */
   double        *rates,             /* (I) Array of cash rates              */
   int           numRates,           /* (I) Number of cash rates to add      */
   long          dayCountConv)       /* (I) See JpmcdsDayCountConvention        */
{
   static char routine[] = "JpmcdsZCCash";

   int       status = FAILURE;
   TCurve    *tcCash = NULL;
   ZCurve    *zcCash = NULL;
   ZCurve    *zcStub = NULL;
  
   if (numRates == 0)    /* nothing to add */
   {
      tcCash = JpmcdsCopyCurve (zeroCurve);
      return tcCash;
   }

   if (CheckZCCashInputs(zeroCurve, dates, rates, numRates, dayCountConv) == FAILURE)
      goto done;

   /* convert the input TCurve into a ZCurve so we can call all the usual
      zero curve routines */
   if (zeroCurve->fNumItems > 0)      /* only want it if we need it */
   {
       zcStub = JpmcdsZCFromTCurve(zeroCurve);

       if (zcStub == NULL)
       {
           JpmcdsErrMsg("%s: failed to make ZCurve structure from input TCurve.\n", routine);
           goto done;
       }
   }

   zcCash = JpmcdsZCMake(zeroCurve->fBaseDate, numRates + 2, ZC_DEFAULT_BASIS, ZC_DEFAULT_DAYCNT);
   if (zcCash == NULL)
   {
       JpmcdsErrMsg("%s: couldn't make Zero Curve for cash.\n", routine);
       goto done;
   }

   if (JpmcdsZCAddMoneyMarket(zcCash, dates, rates, numRates, dayCountConv) == FAILURE)
   {
       JpmcdsErrMsg("%s: Adding cash instruments failed.\n", routine);
       goto done;
   }

   /* combine the stub curve with the cash part */
   if (zcStub != NULL)
   {
       if (JpmcdsZCAddPrefixCurve(zcCash, zcStub) == FAILURE)
       {
           JpmcdsErrMsg("%s: Error adding Cash and stub zero curve together.\n", routine);
           goto done;
       }

       if (JpmcdsZCAddSuffixCurve(zcCash, zcStub) == FAILURE)
       {
           JpmcdsErrMsg("%s: Error adding cash and stub zero curve together(2).\n", routine);
           goto done;
       }
   }

   tcCash = JpmcdsZCToTCurve(zcCash);
   if (tcCash == NULL)
   {
       JpmcdsErrMsg("%s: couldn't make a TCurve copy of Cash zero curve.\n", routine);
       goto done;
   }

   status = SUCCESS;
   goto done;

   done: 
       if (status == FAILURE)
       {
           JpmcdsFreeTCurve(tcCash);
           tcCash = NULL;
           JpmcdsErrMsg("%s: Failed.\n", routine);
       }
   
   JpmcdsZCFree(zcCash);
   JpmcdsZCFree(zcStub);
   return tcCash;
}


/*
***************************************************************************
** Adds swap points to a given zero curve.
***************************************************************************
*/
TCurve* JpmcdsZCSwaps(
 TCurve          *zeroCurve,         /* (I) Zero curve to add swap points to*/
 TCurve          *discZC,            /* (I) Zero curve used for discounting */
 TDate           *dates,             /* (I) Unadjusted swap maturity dates  */
 double          *rates,             /* (I) Swap par fixed rates (0.06=6%)  */
 int             numSwaps,           /* (I) Len of dates,rates,prices       */
 int             fixedSwapFreq,      /* (I) Fixed leg coupon frequency      */
 int             floatSwapFreq,      /* (I) Floating leg coupon frequency   */
 long            fixDayCountConv,    /* (I) See JpmcdsDayCountConvention       */
 long            floatDayCountConv,  /* (I) See JpmcdsDayCountConvention       */
 char            fwdLength,          /* (I) For fwd smoothing, length of fwds*/
 TBadDayAndStubPos badDayAndStubPos, /* (I) Bad day conv + stub pos         */
 char            *holidayFile)       /* (I) See JpmcdsBusinessDay              */
{
   static char routine[] = "JpmcdsZCSwaps";

   int          status=FAILURE;      /* until proven successful              */
   int          offset;              /* offset in array of swap rates        */
   TInterpData  *interpData=NULL;    /* interpolation data for zero rates    */
   TDate        lastStubDate;
   TCurve       *tcSwaps = NULL;
   ZCurve       *zcSwaps = NULL;
   TBadDayList  *badDayList = NULL;
   long         badDayConv;
   TStubPos     stubPos;

   
   if (numSwaps == 0)       /* nothing to add */
   {
       tcSwaps = JpmcdsCopyCurve (zeroCurve);
       return tcSwaps;
   }

   if (JpmcdsBadDayAndStubPosSplit(badDayAndStubPos, &badDayConv, &stubPos) != SUCCESS)
        goto done;

   /* Check inputs */
   if (CheckZCSwapsInputs(
             zeroCurve,
             discZC,
             dates,
             rates,
             numSwaps,
             fixedSwapFreq,
             floatSwapFreq,
             fixDayCountConv,
             floatDayCountConv,
             fwdLength,
             badDayConv,
             holidayFile)==FAILURE)
   {
      goto done;
   }


   /* convert the input TCurve into a ZCurve so we can call all the usual zero curve routines */
   zcSwaps = JpmcdsZCFromTCurve(zeroCurve);
   if (zcSwaps == NULL)
   {
       JpmcdsErrMsg("%s: failed to make ZCurve structure from input TCurve.\n", routine);
       goto done;
   }

   /* only want to add swap points not already covered by stub curve,
      don't want to replace anything in the existing curve */

   if (zeroCurve->fNumItems < 1)
   {
       lastStubDate = zeroCurve->fBaseDate;
   }
   else
   {
       lastStubDate = zeroCurve->fArray[zeroCurve->fNumItems-1].fDate;
   }

   offset = 0;

   while (numSwaps > 0 && dates[offset] < lastStubDate)
   {
       offset++;
       numSwaps--;    
   }

   if (numSwaps > 0)
   {
       if (JpmcdsZCAddSwaps(zcSwaps,
                      discZC,
                      &dates[offset],
                      &rates[offset],
                      numSwaps,
                      fixedSwapFreq,
                      floatSwapFreq,
                      fixDayCountConv,
                      floatDayCountConv,
                      JPMCDS_FLAT_FORWARDS,
                      interpData,
                      badDayList,
                      badDayAndStubPos,
                      holidayFile)==FAILURE)
       {
           goto done;
       }
   }
 
   tcSwaps = JpmcdsZCToTCurve(zcSwaps);
   if (tcSwaps == NULL)
   {
       JpmcdsErrMsg("%s: couldn't make a TCurve copy of swaps zero curve.\n", routine);
       goto done;
   }

   status = SUCCESS;
   goto done;
    
 done:
   FREE(interpData);
   if (status == FAILURE)
   {
       JpmcdsFreeTCurve(tcSwaps);
       tcSwaps = NULL;
       JpmcdsErrMsg("%s: Failed.\n", routine);
   }

   JpmcdsZCFree(zcSwaps);
   return tcSwaps;
}


/*
***************************************************************************
** Check input arguments for adding cash to a given zero curve.
***************************************************************************
*/
static int CheckZCCashInputs(
        TCurve    *zeroCurve,
        TDate     *dateArr,
        double    *rateArr,
        int       numRates,
        long      dayCountConv)
{
    int   i;
    int   status = SUCCESS;     /* Until proven a failure */
    static char routine[] = "CheckZCCashInputs";

    /* Check whether zeroCurve exists 
     */
    if(zeroCurve == NULL)
    {
        JpmcdsErrMsg("%s: ZeroCurve is NULL\n", routine);
        status = FAILURE;
    }

    if (numRates < 0)
    {
        JpmcdsErrMsg("%s: Number of cash points (%d) < 0.\n", 
                  routine, numRates);
        status = FAILURE;
    }

    if (JpmcdsDayCountValid(routine,dayCountConv)==FAILURE)
    {
        JpmcdsErrMsg("%s: Bad day count convention (%ld).\n", 
                  routine, dayCountConv);
        status = FAILURE;
    }

    for (i=0 ;i<numRates ;i++ )
    {
       if (i>0 && (dateArr[i-1] > dateArr[i]))
       {
           JpmcdsErrMsg("%s: date #%d (%s) is after date #%d (%s).\n",
                      routine,
                      i,
                      JpmcdsFormatDate(dateArr[i-1]),
                      i+1,
                      JpmcdsFormatDate(dateArr[i]));
           status = FAILURE;
       }
        
       if (rateArr[i] <= -1.0 )
       {
           JpmcdsErrMsg("%s: rate #%d (%f) <= -1.\n", 
                     routine, i+1, rateArr[i]);
           status = FAILURE;
       }

       if (dateArr[i] <= zeroCurve->fBaseDate)
       {
           JpmcdsErrMsg("%s: Input cash date (%s) on or before value date (%s).\n",
                     routine, 
                     JpmcdsFormatDate(dateArr[i]),
                     JpmcdsFormatDate(zeroCurve->fBaseDate));
           status = FAILURE;
       }
    }
   

    return status;
}


/*
***************************************************************************
** Check input arguments for adding swaps to a given zero curve.
***************************************************************************
*/
static int CheckZCSwapsInputs(
        TCurve     *zeroCurve,
        TCurve     *discZC,
        TDate      *dates,
        double     *rates,
        int        numSwaps,
        int        fixedSwapFreq,
        int        floatSwapFreq,
        long       fixDayCountConv,
        long       floatDayCountConv,
        char       fwdLength,
        long       badDayConv,
        char       *holidayFile)
{
   int         i;
   TBoolean    isBusDay;
   int         status = SUCCESS;   /* Until proven a failure */
   static char routine[] = "CheckZCSwapsInputs";


   /* Check whether zeroCurve exists 
    */
    if(zeroCurve == NULL)
    {
        JpmcdsErrMsg("%s: zeroCurve is NULL.\n",routine);
        status = FAILURE;
    }

   if (numSwaps < 0)
   {
       JpmcdsErrMsg("%s: Number of swaps (%d) < 0.\n", routine, numSwaps);
       status = FAILURE;
   }

   if (JpmcdsDayCountValid(routine,fixDayCountConv)==FAILURE)
   {
       JpmcdsErrMsg("%s: Bad day count convention (%ld) on fixed leg.\n", 
                 routine, fixDayCountConv);
       status = FAILURE;
   }

   if (JpmcdsDayCountValid(routine, floatDayCountConv)==FAILURE)
   {
       JpmcdsErrMsg("%s: Bad day count convention (%ld) on floating leg.\n", 
                 routine, floatDayCountConv);
       status = FAILURE;
   }

   if (fixedSwapFreq != 0) 
   {
       if (12%fixedSwapFreq != 0)
       {
           JpmcdsErrMsg("%s: Invalid swap frequency (%d) on fixed leg.\n", 
                     routine, fixedSwapFreq);
           status = FAILURE;
       }
   }

   if (floatSwapFreq != 0)
   {
       if (12%floatSwapFreq != 0)
       {
           JpmcdsErrMsg("%s: Invalid swap frequency (%d) on floating leg.\n", 
                     routine, floatSwapFreq);
           status = FAILURE;
       }
   }

   if (JpmcdsBadDayConvValid(routine, badDayConv)==FAILURE)
   {
       JpmcdsErrMsg("%s: Invalid bad day convention (%c).\n", 
                 routine, (char)badDayConv);
       status = FAILURE;
   }

   for (i=0;i<numSwaps ;i++ )
   {
       if (rates[i] <= -1 )
       {
           JpmcdsErrMsg("%s: rate #%d (%f) <= -1.\n",
                     routine, i+1, rates[i]);
           status = FAILURE;
       }
   }

   for (i=0;i<numSwaps-1 ;i++ )
   {
       if (dates[i+1] <= dates[i])
       {
           JpmcdsErrMsg("%s: date #%d (%s) is after date #%d (%s).\n",
                      routine,
                      i+1, JpmcdsFormatDate(dates[i]),
                      i+2, JpmcdsFormatDate(dates[i+1]));
           status = FAILURE;
       }
   }

   for (i=0;i<numSwaps ;i++ )
   {
       if (dates[i] < zeroCurve->fBaseDate)
       {
           JpmcdsErrMsg("%s: swap #%d (%s) before value date (%s).\n", 
                     routine, i+1, 
                     JpmcdsFormatDate(dates[i]),
                     JpmcdsFormatDate(zeroCurve->fBaseDate));
           status = FAILURE;
       }
   }


   /* Check whether valueDate is a business date 
    */
   if(JpmcdsIsBusinessDay(zeroCurve->fBaseDate, holidayFile, &isBusDay)==FAILURE)
   {
       JpmcdsErrMsg("%s: JpmcdsIsBusinessDay failed for valueDate.\n",routine);
       status = FAILURE;
       goto done;               /* Cant check isBusDay now. */
   }

   if(isBusDay == FALSE)
   {
       JpmcdsErrMsg("%s: value date %s is not a business day\n",
                 routine, JpmcdsFormatDate(zeroCurve->fBaseDate));
       status = FAILURE;
   }

                 
 done:
   return status;
}
