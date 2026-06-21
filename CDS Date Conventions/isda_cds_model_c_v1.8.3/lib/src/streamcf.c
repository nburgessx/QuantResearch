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
#include <string.h>

#include "bastypes.h"
#include "cashflow.h"
#include "cerror.h"
#include "cmemory.h"
#include "dateadj.h"
#include "fltrate.h"
#include "interp.h"
#include "macros.h"
#include "schedule.h"
#include "streamcf.h"
#include "stub.h"
#include "zr2fwd.h"
#include "zr2coup.h"


/*
***************************************************************************
** Gets the payment for a back stub on a floating leg.
***************************************************************************
*/
static int GetBackStub
     (TCurve        *zeroCurve,      /* (I) Zero Curve          */
      long           zcInterpType,   /* (I) Zero Curve interp   */
      TStreamFloat  *stream,         /* (I) One stream of swap  */
      long          rateBadDayConv,  /* (I) Bad day conv for rate*/
      char          *holidayFile,    /* (I) Holiday file        */
      double        *stubAmount);    /* (O) Back stub rate      */


/*
***************************************************************************
** Gets the payment for a front stub on a floating leg.
***************************************************************************
*/
static int GetFrontStub
    (TCurve        *zeroCurve,     /* (I) Zero Curve        */
     long           zcInterpType,  /* (I) Zero curve interp type  */
     TStreamFloat  *stream,        /* (I) one stream of swap*/
     long          rateBadDayConv, /* (I) bad day conv      */
     char          *holidayFile,   /* (I) holiday file      */
     double        *stubAmount);   /* (O) front stub rate   */


/*
***************************************************************************
** Gets stub rates for a back stub on a floating leg.
***************************************************************************
*/
static int GetFloatStubRate
    (TCurve            *zeroCurve,      /* (I) Zero Curve           */
     long               zcInterpType,   /* (I) Zero curve interp    */
     TDateAdjIntvl     *spotOffset,     /* (I) Period to spot       */
     double             averSoFar,      /* (I) Known avg rate portion*/
     double             compSoFar,      /* (I) Known comp rate portn*/
     long               dayCountConv,   /* (I) See JpmcdsDayCountConv  */
     TCouponDates       coupon,         /* (I) Coupon forming stub  */
     TSwapType          swapType,       /* (I) The type of swap     */
     TDateInterval      stubInterval,   /* (I) Maturity of stub rate*/
     long               rateBadDayConv, /* (I) Bad day conv for rate*/
     char              *holidayFile,    /* (I) See JpmcdsBusinessDay   */
     double            *stubRate);      /* (O) Rate for stub period */


/*
***************************************************************************
** Checks the details for a fixed stream of a swap.
***************************************************************************
*/
static int CheckFixedStreamDetails(
     TStreamFixed      *stream,     /* (I) Stream of swap to be checked      */
     TDate             valueDate,   /* (I) Value Date                        */
     TBoolean          pricing);     /* (I) If pricing, no fixed rate expected*/


/*
***************************************************************************
** Checks the input arguments to JpmcdsFloatingCashFlow.
***************************************************************************
*/
static int CheckFloatingCFArgs(
            TDate         valueDate,
            TDate         indexStartDate,
            TDate         indexMatDate,
            long          rateDayCountConv,
            long          payDayCountConv,
            double        spread);


/*
***************************************************************************
** Calculates the cashflows for a floating stream of a swap.
**
** Principals are included as per specification in the stream->flags field.
***************************************************************************
*/
TCashFlowList* JpmcdsNewStreamFloatCFLGen
    (TCurve             *zeroCurve,          /* (I) Zero Curve              */
     long                zcInterpType,       /* (I) Zero curve interp type  */
     TStreamFloat       *stream,             /* (I) See JpmcdsNewStreamFloat   */
     long               rateBadDayConv,      /* (I) Bad Day Conv for index  */
     char               *holidayFile)        /* (I) Holiday file            */
{
   int           i;                         /* loops over coupon period     */
   int           firstCoupon;               /* where we start looping       */
   int           lastCoupon;                /* where we finish looping      */
   int           status = FAILURE;          /* until proven successful      */
   TBoolean      stubCompound = FALSE;      /* TRUE for tam, FALSE for t4m  */
   TDate         indexStartDate;            /* start date of fwd rate       */
   TDate         indexMatDate;              /* mat date of fwd rate         */
   TDateInterval payInterval;               /* Payment interval             */
   double        *amounts = NULL;           /* cashflows for coupons        */
   TCashFlowList *cfl = NULL;               /* cfl to return                */

   static char routine[] = "JpmcdsNewStreamFloatCFLGen";
   
   amounts = NEW_ARRAY(double, stream->dl->fNumItems);
   if (amounts == NULL)
       goto done;

   if (stream->flags & JPMCDS_SUBTRACT_INITIAL)
   {
        cfl = JpmcdsNewEmptyCFL(stream->dl->fNumItems+1);
   }
   else
   {
       cfl = JpmcdsNewEmptyCFL(stream->dl->fNumItems);
   }

   if (cfl==NULL)
       goto done;

   /* set up how we have to loop over the regular coupons */

    switch (stream->dl->stubLocation)
    {
      case JPMCDS_NO_STUB:
        firstCoupon = 0;
        lastCoupon = stream->dl->fNumItems;
        break;
        
      case JPMCDS_BACK_STUB:
        if(GetBackStub
               (zeroCurve,
                zcInterpType,
                stream,
                rateBadDayConv,
                holidayFile,
                &amounts[stream->dl->fNumItems-1]) == FAILURE)
            {
                JpmcdsErrMsg("%s: Failed to calculate back stub rate.\n", routine);
                goto done;
            }
            
        firstCoupon = 0;
        lastCoupon = stream->dl->fNumItems-1;
        break;

      case JPMCDS_FRONT_STUB:
        if(GetFrontStub
               (zeroCurve,
                zcInterpType,
                stream,
                rateBadDayConv,
                holidayFile,
                &amounts[0]) == FAILURE)
            {
                JpmcdsErrMsg("%s: Failed to calculate front stub rate.\n", routine);
                goto done;
            }
            
            firstCoupon = 1;
            lastCoupon = stream->dl->fNumItems;
        break;
        
        case JPMCDS_TWO_STUB:
        if(GetFrontStub
               (zeroCurve,
                zcInterpType,
                stream,
                rateBadDayConv,
                holidayFile,
                &amounts[0]) == FAILURE)
            {
                JpmcdsErrMsg("%s: Failed to calculate front stub rate.\n", routine);
                goto done;
            }
            
            if(GetBackStub
               (zeroCurve,
                zcInterpType,
                stream,
                rateBadDayConv,
                holidayFile,
                &amounts[stream->dl->fNumItems-1]) == FAILURE)
            {
                JpmcdsErrMsg("%s: Failed to calculate back stub rate.\n", routine);
                goto done;
            }
            
        firstCoupon = 1;
        lastCoupon = stream->dl->fNumItems-1;
        break;
        
      default:
            JpmcdsErrMsg("%s: unknown stub location (%ld).\n", routine, stream->dl->stubLocation);
            goto done;
    } /* switch */

   for (i=firstCoupon; i<lastCoupon ;i++)
   {
       switch (stream->swapType)
       {
         case JPMCDS_SINGLE_REFIX:
           indexStartDate = stream->dl->fArray[i].resetDate;
           
           if (JpmcdsDtFwdAdj
               (indexStartDate,
                &stream->floatingRate->spotOffset,
                &indexStartDate) == FAILURE)
           {
               JpmcdsErrMsg("%s: Error when finding index Date\n", routine);
               goto done;
           }
           
           if (JpmcdsDtFwdAny(indexStartDate,
                           &stream->floatingRate->matInterval,
                           &indexMatDate)==FAILURE)
           {
               JpmcdsErrMsg("%s: error finding index maturity date(%d)\n", routine, i);
               goto done;
           }
           
           if (JpmcdsBusinessDay
               (indexMatDate,
                rateBadDayConv,
                holidayFile,
                &indexMatDate) == FAILURE)
           {
               JpmcdsErrMsg("%s: error adjusting index mat date (%d)\n", routine, i);
               goto done;
           }

           /* If the maturity and payInterval are the same, then
            * make sure that the holiday adjustments which affected
            * the indexStartDate and the index maturity date are 
            * reflected in the payment interval. Otherwise we may
            * get a par fixed rate with a very small stub, when in
            * fact we just want a simple LIBOR (zero-coupon) type 
            * rate. 
            */
           if (TDATE_INTERVALS_EQUAL(stream->floatingRate->payInterval,
                                     stream->floatingRate->matInterval))
           {
               if (JpmcdsMakeDateInterval
                   (indexMatDate-indexStartDate,'D', &payInterval) 
                   == FAILURE)
                   goto done;
           }
           else
           {
               payInterval = stream->floatingRate->payInterval;
           }
          if (JpmcdsFloatingCashFlowGen
               (zeroCurve,
                zcInterpType,
                indexStartDate,
                indexMatDate,
                payInterval,
                &stream->dl->fArray[i],
                stream->floatingRate->dayCountConv,
                stream->payDayCountConv,
                stream->principal,
                stream->floatingRate->spread,
                &amounts[i]) == FAILURE)
           {
               JpmcdsErrMsg("%s: Failed to get floating cash flow.\n", routine);
               goto done;
           }
           break;
            
         default:
           JpmcdsErrMsg("%s: Unknown swap type (%d).\n", routine, i);
           goto done;
       } /* end of switch */
   }     /* end of i loop over c'flows */


   for (i=0;i<stream->dl->fNumItems ;i++ )
   {
       if (stream->flags & JPMCDS_SUBTRACT_INITIAL)
       {
           if (i==0)
           {
               cfl->fArray[i].fAmount = -stream->principal;
               cfl->fArray[i].fDate = stream->dl->fArray[i].accrueStartDate;
           }
           cfl->fArray[i+1].fAmount = amounts[i];
           cfl->fArray[i+1].fDate = stream->dl->fArray[i].payDate;
       }
       else
       {
           cfl->fArray[i].fAmount = amounts[i];
           cfl->fArray[i].fDate = stream->dl->fArray[i].payDate;
       }
   }

   if (stream->flags & JPMCDS_ADD_FINAL)
   {
       cfl->fArray[cfl->fNumItems-1].fAmount += stream->principal;
   }
   
   status = SUCCESS;

 done:
   if (status == FAILURE)
   {
       JpmcdsFreeCFL(cfl);
       cfl = NULL;
       JpmcdsErrMsg("%s: Failed.\n", routine);
   }
   
   FREE(amounts);
   return cfl;
}


/*
***************************************************************************
** Calculates the cashflows for a fixed stream of a swap.
**
** Principals are included as per specification in the stream->flags field.
***************************************************************************
*/
TCashFlowList* JpmcdsNewStreamFixedCFL
    (TDate              valueDate,           /* (I) Value Date              */
     TStreamFixed       *stream)             /* (I) See JpmcdsMakeTStream      */
{
    int              i;
    int              j;
    int              status = FAILURE;
    double           rate;
    double           yearFrac;
    double           *couponRates = NULL;
    TDate            accStartDate;
    TDate            accEndDate;
    TCashFlowList    *cfl = NULL;
    
    static char routine[] = "JpmcdsNewStreamFixedCFL";
    
    couponRates = NEW_ARRAY(double, stream->dl->fNumItems);
    if (couponRates == NULL)
        goto done;

    if(JpmcdsGetFixedCouponRates(valueDate, stream, couponRates) == FAILURE)
    {
        JpmcdsErrMsg("%s: Failed to calculate coupon rates.\n", routine);
        goto done;
    }

    if (stream->flags & JPMCDS_SUBTRACT_INITIAL)
    {
        cfl = JpmcdsNewEmptyCFL(stream->dl->fNumItems+1);
        if (cfl== NULL)
            goto done;

        
        cfl->fNumItems = 0;
        cfl->fArray[0].fDate = stream->dl->fArray[0].accrueStartDate;
        cfl->fArray[0].fAmount = -stream->principal;
        ++cfl->fNumItems;
    }
    else
    {
        cfl = JpmcdsNewEmptyCFL(stream->dl->fNumItems);
        if(cfl==NULL)
            goto done;

        cfl->fNumItems = 0;
    }
 
    for (i=0, j=cfl->fNumItems; i<stream->dl->fNumItems; i++, j++)
    {
         accStartDate = stream->dl->fArray[i].accrueStartDate;
         accEndDate = stream->dl->fArray[i].accrueEndDate;
         
         if(JpmcdsDayCountFraction(
                 accStartDate,
                 accEndDate,
                 stream->payDayCountConv,
                 &yearFrac) == FAILURE)
         {
             JpmcdsErrMsg("%s: error getting day cnt fract for coupon %d\n", routine, i);
             goto done;
         }

         rate = couponRates[i];
         
         cfl->fArray[j].fDate = stream->dl->fArray[i].payDate;
          cfl->fArray[j].fAmount =
           stream->principal * rate * yearFrac;
         ++cfl->fNumItems;
           
         if (i==stream->dl->fNumItems-1 &&
             stream->flags & JPMCDS_ADD_FINAL)
         {
             cfl->fArray[j].fAmount += stream->principal;
         }
    }
    
    status = SUCCESS;

done:
    if (status == FAILURE)
    {
        JpmcdsFreeCFL(cfl);
        cfl = NULL;
        JpmcdsErrMsg("%s: Failed.\n",routine);
    }

    FREE(couponRates);
    return cfl;
}


/*
***************************************************************************
** Same as JpmcdsFloatingCashFlow, except user can choose the zero curve 
** interpolation type. 
***************************************************************************
*/
int JpmcdsFloatingCashFlowGen
   (TCurve         *zeroCurve,           /* (I) Zero curve                  */
    long            zcInterpType,        /* (I) Zero curve interp type      */
    TDate          indexStartDate,       /* (I) Start date of rate          */
    TDate          indexMatDate,         /* (I) Maturity date of rate       */
    TDateInterval  payInterval,          /* (I) Payment interval of rate    */
    TCouponDates   *couponDates,         /* (I) AccStart,AccEnd,Pay,Reset   */
    long           rateDayCountConv,     /* (I) Market convention of rate   */
    long           payDayCountConv,      /* (I) Convention for payments     */
    double         notional,             /* (I) Notional amount             */
    double         spread,               /* (I) Added to the coupon rate    */
    double         *cashflow)            /* (O) Coupon payment              */
{
    int           status = FAILURE;       /* Until successful */
    static char   routine[] = "JpmcdsFloatingCashFlowGen";
    double        couponRate;
    double        yearFrac;
    TFloatRate    floatRateDef; /* Used for convexity adjustment */
    TDateInterval matInterval;  /* Goes into TFloatRate  */

    if (CheckFloatingCFArgs(
                            zeroCurve->fBaseDate,
                            indexStartDate,
                            indexMatDate,
                            rateDayCountConv,
                            payDayCountConv,
                            spread) == FAILURE)
    {
        JpmcdsErrMsg("%s: Checking inputs uncovered error.\n", routine);
        goto done;
    }


    /* Set up TFloatRate, floatRateDef. */
    if (JpmcdsMakeDateInterval(indexMatDate-indexStartDate, 'D', &matInterval) == FAILURE)
        goto done;

    if (JpmcdsFloatRateSet(&matInterval, &payInterval, rateDayCountConv,
                        0       /*SpotOffsetDays*/, 0./*Spread*/, 1./*Weight*/,
                        &floatRateDef) == FAILURE)
        goto done;
   

    if (JpmcdsZerosToCouponsPoint
        (zeroCurve, zcInterpType,
         indexStartDate, &payInterval, indexMatDate,
         rateDayCountConv,
         FALSE, /* Stub at front */
         &couponRate) == FAILURE)
    {
        goto done;
    }

    if (JpmcdsDayCountFraction
        (couponDates->accrueStartDate,
         couponDates->accrueEndDate,
         payDayCountConv,
         &yearFrac) == FAILURE)
    {
        goto done;
    }

    *cashflow = notional * yearFrac * (couponRate + spread);
    status = SUCCESS;

  done:
    if (status == FAILURE)
        JpmcdsErrMsg("%s: Failed.\n", routine);

    return status;
}


/*
***************************************************************************
** Checks the input arguments to JpmcdsFloatingCashFlow.
***************************************************************************
*/
static int CheckFloatingCFArgs
    (TDate         valueDate,
     TDate         indexStartDate,
     TDate         indexMatDate,
     long          rateDayCountConv,
     long          payDayCountConv,
     double        spread)
{
   static char routine[] = "CheckFloatingCFArgs";

   if (indexStartDate < valueDate)
   {
       JpmcdsErrMsg("%s: rate start (%s) is before value date (%s).\n",
                  routine,
                  JpmcdsFormatDate(indexStartDate),
                  JpmcdsFormatDate(valueDate));
       return FAILURE;
   }

   if (indexStartDate > indexMatDate)
   {
       JpmcdsErrMsg("%s: rate maturity (%s) is before rate start (%s).\n",
                  routine,
                  JpmcdsFormatDate(indexMatDate),
                  JpmcdsFormatDate(indexStartDate));
       return FAILURE;
   }

   if (JpmcdsDayCountValid(routine, rateDayCountConv) == FAILURE)
   {
       JpmcdsErrMsg("%s: error when checking rate day count convention.\n", routine);
       return FAILURE;
   }

   if (JpmcdsDayCountValid(routine, payDayCountConv) == FAILURE)
   {
       JpmcdsErrMsg("%s: error when checking pay day count convention.\n", routine);
       return FAILURE;
   }

   return SUCCESS;
}


/*
***************************************************************************
** Given one stream of a swap, this function estimates the rates for each 
** coupon.
***************************************************************************
*/
int JpmcdsGetFixedCouponRates(
     TDate              valueDate,           /* (I) Value Date              */
     TStreamFixed       *stream,             /* (I) See JpmcdsMakeTStram       */
     double             *couponRates)        /* (O) Coupon Rates            */

{
    int        i;                         /* loops over coupon period     */
    int        firstCoupon;               /* where we start looping       */
    int        lastCoupon;                /* where we finish looping      */
    

    static char routine[] = "JpmcdsGetFixedCouponRates";
    
    if (CheckFixedStreamDetails(stream, valueDate, FALSE)==FAILURE)
    {
        JpmcdsErrMsg("%s: checking stream details uncovered error.\n", routine);
        return FAILURE;
    }
    
    /* set up how we loop over the regular coupons */
    switch (stream->dl->stubLocation)
    {
        case JPMCDS_NO_STUB:
        firstCoupon = 0;
        lastCoupon = stream->dl->fNumItems;
        break;
        
    case JPMCDS_BACK_STUB:
        couponRates[stream->dl->fNumItems-1] = stream->finalCouponRate;
        firstCoupon = 0;
        lastCoupon = stream->dl->fNumItems-1;
        break;

    case JPMCDS_FRONT_STUB:
        couponRates[0] = stream->firstCouponRate;
        firstCoupon = 1;
        lastCoupon = stream->dl->fNumItems;
        break;
        
    case JPMCDS_TWO_STUB:
        couponRates[stream->dl->fNumItems-1] = stream->finalCouponRate;
        couponRates[0] = stream->firstCouponRate;
        firstCoupon = 1;
        lastCoupon = stream->dl->fNumItems-1;
        break;
        
    default:
        JpmcdsErrMsg("%s: unknown stub location (%ld).\n", routine, stream->dl->stubLocation);
        return FAILURE;
    }

   for (i=firstCoupon; i<lastCoupon; i++)
   {
        couponRates[i] = stream->fixedRate;
   }
    
   return SUCCESS;
}


/*
***************************************************************************
** Makes a new TStreamFloat.
***************************************************************************
*/
TStreamFloat* JpmcdsNewStreamFloat(
     TCouponDateList    *dl,                 /* (I) See JpmcdsSwapMakeSchedule */
     TSwapType          swapType,            /* (I) The type of swap        */
     double             principal,           /* (I) Principal amount        */
     long               payDayCountConv,     /* (I) See JpmcdsDayCountConvention*/
     long               flags,               /* (I) Principal exchanges     */
     TStubRates         *firstCoupon,        /* (I) Data for first coupon   */
     TStubRates         *finalCoupon,        /* (I) Data for final coupon   */
     double             averSoFar,           /* (I) Known avg rate portion*/
     double             compSoFar,           /* (I) Known compounded rate   */
     TFloatRate         *floatingRate)       /* (I) Defn of floating index  */
{
   int                  i;
   TStreamFloat *stream = NULL;
   int           status = FAILURE;
   static char   routine[] = "JpmcdsNewStreamFloat";
    
   if (swapType > JPMCDS_SINGLE_REFIX)
   {
       JpmcdsErrMsg("%s: unknown swap Type (%d).\n", routine, swapType);
       return NULL;
   }
   if (JpmcdsDayCountValid(routine,payDayCountConv) == FAILURE)
   {
       JpmcdsErrMsg("%s: Failed checking day count convention.\n", routine);
       return NULL;
   }

   stream = JpmcdsNewEmptyTStreamFloat(dl->fNumItems);
    
   if (stream == NULL)
   {
       JpmcdsErrMsg("%s: Failed to generate new TStreamFloat.\n", routine);
       return NULL;
    }
    
   stream->dl->stubLocation = dl->stubLocation;
   for (i=0; i<dl->fNumItems; i++)
   {
       stream->dl->fArray[i] = dl->fArray[i];
   } 
   
   stream->swapType = swapType;
   stream->principal = principal;
   stream->payDayCountConv = payDayCountConv;
   stream->flags = flags;
   stream->firstCoupon = *firstCoupon;
   stream->finalCoupon = *finalCoupon;
   stream->averSoFar = averSoFar;
   stream->compSoFar = compSoFar;

   stream->floatingRate = NEW(TFloatRate);
   if (stream->floatingRate == NULL) 
       goto done;
   *(stream->floatingRate) = *floatingRate;

   status = SUCCESS;

done:
   if (status == FAILURE)
   {
       JpmcdsErrMsg("%s: Failed.\n", routine);
       JpmcdsFreeTStreamFloat(stream);
       stream = NULL;
   }

   return stream;  
}


/*
***************************************************************************
** Makes a new TStreamFixed.
***************************************************************************
*/
TStreamFixed* JpmcdsNewStreamFixed(
     TCouponDateList    *dl,                 /* (I) See JpmcdsSwapMakeSchedule */
     TSwapType          swapType,            /* (I) The type of swap        */
     double             principal,           /* (I) Principal amount        */
     long               payDayCountConv,     /* (I) See JpmcdsDayCountConvention*/
     long               flags,               /* (I) Principal exchanges     */
     double             firstCouponRate,     /* (I) Data for first coupon   */
     double             finalCouponRate,     /* (I) Data for final coupon   */
     double             fixedRate)           /* (I) Fixed Coupon Rate       */
{
   int                  i;
   TStreamFixed     *stream = NULL; 
    
   static char routine[] = "JpmcdsNewStreamFixed";
    
   if (swapType > JPMCDS_SINGLE_REFIX)
   {
       JpmcdsErrMsg("%s: unknown swap type.\n", routine);
       return NULL;
   }

   if (JpmcdsDayCountValid(routine,payDayCountConv)==FAILURE)
   {
      JpmcdsErrMsg("%s: failed on checking pay day count convention.\n", routine);
      return NULL;
   }
   
   stream = JpmcdsNewEmptyTStreamFixed(dl->fNumItems);
    
   if (stream == NULL)
   {
       JpmcdsErrMsg("%s: failed to generate new TStreamFixed.\n", routine);
       return NULL;
   }

   stream->dl->stubLocation = dl->stubLocation;
   for (i=0; i<dl->fNumItems; i++)
   {
       stream->dl->fArray[i] = dl->fArray[i];
   } 

   stream->swapType = swapType;
   stream->principal = principal;
   stream->payDayCountConv = payDayCountConv;
   stream->flags = flags;
   stream->firstCouponRate = firstCouponRate;
   stream->finalCouponRate = finalCouponRate;
   stream->fixedRate = fixedRate;
    
   return stream;
}


/*
***************************************************************************
** Creates a new empty TStreamFloat.
***************************************************************************
*/
TStreamFloat* JpmcdsNewEmptyTStreamFloat(int numPeriods)
{
    static char routine[] = "JpmcdsNewEmptyTStreamFloat";
    
    TStreamFloat    *stream = NULL;
    
    if (numPeriods < 0)
    {
         JpmcdsErrMsg("%s: number of coupon periods (%d) must be >= 0.\n", routine, numPeriods);
         return NULL;
    }
    
    stream = NEW(TStreamFloat);
    
    if (stream == NULL)
    {
        JpmcdsErrMsg("%s: out of memory when allocating TStreamFloat.\n", routine);
        return NULL;
    }
    
    stream->dl = JpmcdsNewEmptyCouponDateList(numPeriods);
    
    if (stream->dl == NULL)
    {
        FREE(stream);
        JpmcdsErrMsg("%s: out of memory when allocating TStreamFloat.\n", routine);
        return NULL;
    }
    
    return stream;
}


/*
***************************************************************************
** Creates a new empty TStreamFixed.
***************************************************************************
*/
TStreamFixed* JpmcdsNewEmptyTStreamFixed(int  numPeriods)
{
    static char routine[] = "JpmcdsNewEmptyTStreamFixed";
    
    TStreamFixed    *stream = NULL;
    
    if (numPeriods < 0)
    {
         JpmcdsErrMsg("%s: number of coupon periods (%d) must be >= 0.\n", routine, numPeriods);
         return NULL;
    }
    
    stream = NEW(TStreamFixed);
    
    if (stream == NULL)
    {
        JpmcdsErrMsg("%s: out of memory when allocating TStreamFixed.\n",routine);
        return NULL;
    }
    
    stream->dl = JpmcdsNewEmptyCouponDateList(numPeriods);
    
    if (stream->dl == NULL)
    {
        FREE(stream);
        JpmcdsErrMsg("%s: out of memory when allocating TStreamFixed.\n",routine);
        return NULL;
    }
    
    return stream;
}


/*
***************************************************************************
** Frees a floating stream. 
***************************************************************************
*/
void JpmcdsFreeTStreamFloat(TStreamFloat *theStream)
{
    if (theStream != NULL)
    {
        JpmcdsFreeCouponDateList(theStream->dl);
        FREE(theStream->floatingRate);
        FREE(theStream);
    }
}


/*
***************************************************************************
** Frees a fixed stream. 
***************************************************************************
*/
void JpmcdsFreeTStreamFixed(TStreamFixed *theStream)
{
    if (theStream != NULL)
    {
        JpmcdsFreeCouponDateList(theStream->dl);
        FREE(theStream);
    }
}


/*
***************************************************************************
** Gets the payment for a front stub on a floating leg.
***************************************************************************
*/
static int GetFrontStub
    (TCurve        *zeroCurve,     /* (I) Zero Curve        */
     long           zcInterpType,  /* (I) Zero curve interp type  */
     TStreamFloat  *stream,        /* (I) one stream of swap*/
     long          rateBadDayConv, /* (I) bad day conv      */
     char          *holidayFile,   /* (I) holiday file      */
     double        *stubAmount)    /* (O) front stub rate   */
{
    double    yearFrac;
    double    stubRate;

    static char routine[] = "GetFrontStub";

    if (GetFloatStubRate
        (zeroCurve,
         zcInterpType,
         &stream->floatingRate->spotOffset,
         stream->averSoFar,
         stream->compSoFar,
         stream->floatingRate->dayCountConv,
         stream->dl->fArray[0],
         stream->swapType,
         stream->firstCoupon.stubInterval,
         rateBadDayConv,
         holidayFile,
         &stubRate) == FAILURE)
    {
        JpmcdsErrMsg("%s: failed to calculate floating front stub.\n", routine);
        return FAILURE;
    }

    /* Now add floating spread as done for all the other floating c'flows */
    /* $$$  NB: This can be improved when the floating rate for the stubs */
    /* is changed to be described by a TFloatRate (LB).                   */
    stubRate += stream->floatingRate->spread;

    if (JpmcdsDayCountFraction(
           stream->dl->fArray[0].accrueStartDate,
           stream->dl->fArray[0].accrueEndDate,
           stream->payDayCountConv,
           &yearFrac) == FAILURE)
    {
        JpmcdsErrMsg("%s: error getting day cnt fraction.\n", routine);
        return FAILURE;
    }
    
    *stubAmount = 
        stubRate * stream->principal * yearFrac;
   
   return SUCCESS;
}


/*
***************************************************************************
** Gets the payment for a back stub on a floating leg.
***************************************************************************
*/
static int GetBackStub
     (TCurve        *zeroCurve,      /* (I) Zero Curve          */
      long           zcInterpType,   /* (I) Zero Curve interp   */
      TStreamFloat  *stream,         /* (I) One stream of swap  */
      long           rateBadDayConv, /* (I) Bad day conv for rate*/
      char          *holidayFile,    /* (I) Holiday file        */
      double        *stubAmount)     /* (O) Back stub rate      */
{
    double    yearFrac;
    double    stubRate;

    static char routine[] = "GetBackStub";
    
    if(GetFloatStubRate
       (zeroCurve,
        zcInterpType,
        &stream->floatingRate->spotOffset,
        stream->averSoFar,
        stream->compSoFar,
        stream->floatingRate->dayCountConv,
        stream->dl->fArray[stream->dl->fNumItems-1],
        stream->swapType,
        stream->finalCoupon.stubInterval,
        rateBadDayConv,
        holidayFile,
        &stubRate) == FAILURE)
    {
        JpmcdsErrMsg("%s: Failed to calculate floating back stub.\n", routine);
        return FAILURE;
    }
    
    /* Now add floating spread as done for all other floating c'flows */
    /* $$$  NB: This can be improved when the floating rate for stubs */
    /* is changed to be described by a TFloatRate (LB).               */
    stubRate += stream->floatingRate->spread;

    if (JpmcdsDayCountFraction
        (stream->dl->fArray[stream->dl->fNumItems-1].accrueStartDate,
         stream->dl->fArray[stream->dl->fNumItems-1].accrueEndDate,
         stream->payDayCountConv,
         &yearFrac) == FAILURE)
    {
        JpmcdsErrMsg("%s: error getting day cnt fraction.\n", routine);
        return FAILURE;
    }

    *stubAmount = 
        stubRate * stream->principal * yearFrac;
    
    return SUCCESS;
}


/*
***************************************************************************
** Gets stub rates for a back stub on a floating leg.
***************************************************************************
*/
static int GetFloatStubRate
    (TCurve            *zeroCurve,      /* (I) Zero Curve           */
     long               zcInterpType,   /* (I) Zero curve interp    */
     TDateAdjIntvl     *spotOffset,     /* (I) Period to spot       */
     double             averSoFar,      /* (I) Known avg rate portion*/
     double             compSoFar,      /* (I) Known comp rate portn*/
     long               dayCountConv,   /* (I) See JpmcdsDayCountConv  */
     TCouponDates       coupon,         /* (I) Coupon forming stub  */
     TSwapType          swapType,       /* (I) The type of swap     */
     TDateInterval      stubInterval,   /* (I) Maturity of stub rate*/
     long               rateBadDayConv, /* (I) Bad day conv for rate*/
     char              *holidayFile,    /* (I) See JpmcdsBusinessDay   */
     double            *stubRate)       /* (O) Rate for stub period */
{   
    static char   routine[] = "GetFloatStubRate";
    int           status = FAILURE;     /* Until proven successful */
    double        rate;
    
    TDate         startDate;
    TDate         matDate;    

    /* use interpolated stub */
    startDate = coupon.accrueStartDate;
    matDate = coupon.accrueEndDate;

    switch (swapType)
    {
      case JPMCDS_SINGLE_REFIX:
        if (JpmcdsForwardFromZCurve
            (zeroCurve,
             zcInterpType,
             startDate,
             matDate,
             dayCountConv,
             JPMCDS_SIMPLE_BASIS,
             &rate) == FAILURE)
        {
            JpmcdsErrMsg("%s: Failed on libor floating back stub rate.\n", routine);
            goto done;
        }
        break;
    
      default:
        JpmcdsErrMsg("%s: unknown swap type. Can't calculate stub rate.\n", routine);
        goto done;
    }

    status = SUCCESS;


 done:
   if (status == FAILURE)
       JpmcdsErrMsg("%s: Failed to calculate back stub rate.\n", routine);

    *stubRate = rate;

    return status;
}


/*
***************************************************************************
** Checks the details for a fixed stream of a swap.
***************************************************************************
*/
static int CheckFixedStreamDetails(
     TStreamFixed      *stream,     /* (I) Stream of swap to be checked      */
     TDate             valueDate,   /* (I) Value Date                        */
     TBoolean          pricing)     /* (I) If pricing, no fixed rate expected*/
{
    int status = SUCCESS;
    static char routine[] = "CheckFixedStreamDetails";
    
    if(JpmcdsCouponDateListCheck(stream->dl)==FAILURE)
    {
       JpmcdsErrMsg("%s: Stream has incorrect schedule\n", routine);
       status = FAILURE;
    }
     
    if (stream->swapType > JPMCDS_SINGLE_REFIX)
    {
        JpmcdsErrMsg("%s: Unknown swap type.\n", routine);
        status = FAILURE;
    }

    if (JpmcdsDayCountValid(routine,stream->payDayCountConv)==FAILURE)
    {
       JpmcdsErrMsg("%s: Invalid pay day count convention.\n", routine);
       status = FAILURE;
    }
    
    if (stream->principal == 0.0)
    {
        JpmcdsErrMsg("%s: Stream has principal = 0.\n", routine);
        status = FAILURE;
    }
    
    if (pricing == FALSE &&
        (stream->fixedRate < -1.0 ||
         stream->fixedRate > JPMCDS_MAX_RATE))
    {
        JpmcdsErrMsg("%s: stream has bad fixed rate (%f).\n", routine,
                   stream->fixedRate);
        status = FAILURE;
    }    
    
    return status;    
}
