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
#include <float.h>
#include <memory.h>
#include "cgeneral.h"
#include "bastypes.h"
#include "badday.h"
#include "cashflow.h"
#include "cerror.h"
#include "cmemory.h"
#include "rtbrent.h"
#include "convert.h"
#include "datelist.h"
#include "date_sup.h"
#include "ldate.h" 
#include "lintrp.h"
#include "lintrp.h"
#include "macros.h"
#include "schedule.h"
#include "streamcf.h"
#include "stub.h"
#include "tcurve.h"
#include "zcprvt.h"
#include "zcswdate.h"



/* Defaults for JpmcdsNewCouponDatesSwap. */
#define PAY_OFFSET             0
#define RESET_OFFSET           0
#define IN_ARREARS             FALSE

#define LONG_STUB              FALSE
#define ADJUST_LAST_ACC_DATE   TRUE
#define FIRST_ROLL_DATE        0L
#define LAST_ROLL_DATE         0L
#define FULL_FIRST_COUPON_DATE 0L

/* Define constants for JpmcdsRootFindBrent */
#define   INITIAL_X_STEP      (1E-8)
#define   INITIAL_F_DERIV     (0)
#define   X_TOLERANCE         (DBL_MAX)
#define   F_TOLERANCE         (1E-11)
#define   MAX_ITERATIONS      (50)
#define   LOWER_BOUND         (-.99)
#define   UPPER_BOUND         (100.)



/*
***************************************************************************
** Optimization, adds a swap to a zero curve assuming the one added just 
** before has the same coupon dates (w/ one additional payment).  This is
** much quicker than generating a cash-flow-list for the swap and adding those.
** And it often happens that this routine can be used.
**
** ZCurve is updated to contain swap information.
**
** Assumptions:
**   Naturally it assumed that the new swap is equal to the old one with 
**   an additional coupon.  Also assume that the previous couponRate!=0.0,
**   and that the last discount factor in the zero curve corresponds to the date
**   of the previous swap and was placed there based on the previous swap.
***************************************************************************
*/
static int AddSwapFromPrevious(           /* adds swap based on one just added */
   ZCurve                      *zc,        /* (M) ZCurve to add to */
   TDate                        dateNew,   /* (I) date of new swap */
   double                       rateNew,   /* (I) rate of new swap */
   double                       priceNew,  /* (I) price of new swap */
   TDate                        dateOld,   /* (I) date of previous swap */
   double                       rateOld,   /* (I) rate of previous swap */
   double                       priceOld,  /* (I) price of previous swap */
   long                         dayCount); /* (I) day counting convention */


/*
***************************************************************************
** Objective function for root-finder to use.
***************************************************************************
*/
static int JpmcdsObjFuncPVtheSwap(
        double    rateGuess,          /* (I) variable x - guess of rate     */
        void      *p,                 /* (I) data ptr (TPfunctionParams)  */
        double    *f);                /* (O) f(x) == value to make zero   */


/*
***************************************************************************
** Calculates the PV of a swap.
***************************************************************************
*/
static int CalcPV
    (TCurve           *zeroCurve,       /* (I) Zero curve for estim cfls  */ 
     long              interpType,      /* (I) Type of zero interpolation */
     TCashFlowList    *fixedFlows,      /* (I) CFLs of fixed leg          */
     TStreamFloat     *floatStream,     /* (I) Floating leg               */
     double            price,           /* (I) Price - 1.0==par           */
     long              rateBadDayConv,  /* (I) Bad day convention         */
     char             *holidayFile,     /* (I) Name of holiday file       */
     TCurve           *discZC,          /* (I) Discount zero curve        */
     double           *dealPV);         /* (O) PV of the swap             */


typedef struct    /* parameters for passing into the solver */
{
   TCashFlowList    *fixedCfls;
   TCurve           *zeroCurve;
   long             interpType;
   TStreamFloat     *floatStream;
   double           price;
   long             rateBadDayConv;
   char             *holidayFile;
   TCurve           *discZC;
}TPfunctionParams;


/*
***************************************************************************
** Adds a strip of swap instruments to a ZCurve.
**
** The zero interpolation kicks in whenever a point is calculated to
** determine the discount factor for a cash flow, deep inside the objective
** function called by the root finder, which calls JpmcdsZCInterpolate() to 
** determine the appropriate discount factor.  
**
** In this routine, there are two mutually exclusive ways to specify 
** bad days: badDayList, and badDayConv/holidayFile. If badDayList is used, 
** bad business day lists can be used to make sure swap dates near to
** canonical dates (i.e. integral number of years from start date) are
** properly used for coupon dates of later instruments.  Say the 2-year
** point has been adjusted because of holidays to be 3 days later.  The
** coupon at 2-years for the 3 (and subsequent) swap instruments must be
** adjusted to that date.  To do this, the badDayList contains a list of
** <good,bad> day pairs, which are used to adjust dates to/from bad
** business days.
**
** INPUTS: ZCurve to add swaps to, maturity dates and rates for them, how many
** there are, coupon frequency: annual or semi-annual, dayCountConvention,
** interpolation method and data, bad business day list. 
** WARNING: In this routine, there are two mutually exclusive ways to specify
** bad days: badDayList, and badDayConv/holidayFile. 
** If badDayList is non-NULL, the maturity dates passed in must be
** adjusted. If the badDayConv is not JPMCDS_BAD_DAY_NONE, the maturity dates
** passed in must be UN-adjusted. 
**
** OUTPUTS: Extended zero curve. 
***************************************************************************
*/
int  JpmcdsZCAddSwaps(            /* adds a strip of swaps to ZCurve */
	ZCurve          *zc,               /* (M) ZCurve to add to        */
	TCurve          *discZC,           /* (I) Discount zero curve     */
	TDate           *inDates,          /* (I) MatDates; adj. if badDayList!=0*/
	double          *inRates,          /* (I) Fixed rates (0.06 for 6%) */
	int              numSwaps,         /* (I) # instruments to add    */
	int              fixedSwapFreq,    /* (I) Fixed leg freq          */
	int              floatSwapFreq,    /* (I) Floating leg freq       */
	long             fixDayCountConv,  /* (I) Convention for fixed leg*/
	long             floatDayCountConv,/* (I) Convention for float leg*/
	long             interpType,       /* (I) Zero interpolation method */
	TInterpData     *interpData,       /* (I) Zero interpolation data */
	TBadDayList     *badDayList,       /* (I) Bad-good day pairs; if non-NULL */
	/*     maturity dates must be adjusted */
	TBadDayAndStubPos badDayAndStubPos,/* (I) See JpmcdsBusinessDay      */
	char            *holidayFile)      /* (I) See JpmcdsBusinessDay      */
{
	static char    routine[] = "JpmcdsZCAddSwaps";
	int            status = FAILURE;       /* Until proven successful */

	TBoolean       oneAlreadyAdded;        /* If a swap was already added */
	int            i;                      /* loops over all instruments */

	TSwapDates    *swapDates = NULL;         /* Dates from caller */
	double        *swapRates = NULL;

	TSwapDates    *newDates = NULL;        /* Swap dates w/ synthetics */
	double        *newRates = NULL;        /* Swap rates w/ synthetics */

	TBoolean       useFastZC = FALSE;      /* whether to use swap zc */
	long           badDayConv;
	TStubPos       stubPos;


	if (zc == NULL)     /* need a ZCurve to start with */
	{
		JpmcdsErrMsg("%s: input zero curve must contain data.\n", routine);
		goto done;
	}

	if (JpmcdsBadDayAndStubPosSplit(badDayAndStubPos, &badDayConv, &stubPos) != SUCCESS)
		goto done;

	if (badDayList != NULL && badDayConv != JPMCDS_BAD_DAY_NONE)
	{
		JpmcdsErrMsg("%s: Bad days can be defined either by badDayList or\n"
			"\tbadDayConv, but not both.\n", routine);
		goto done;
	}

	/* Set up TSwapDates with input swap maturity dates. */
	if (badDayList != NULL)    /* This means dates are adjusted already */
	{
		swapDates = JpmcdsSwapDatesNewFromAdjusted
			(zc->valueDate, fixedSwapFreq, inDates, numSwaps, badDayList);
		if (swapDates == NULL)
			goto done;
	}
	else                         /* Dates NOT adjusted */
	{
		swapDates = JpmcdsSwapDatesNewFromOriginal
			(zc->valueDate, fixedSwapFreq, inDates, numSwaps,
			badDayList, badDayConv, holidayFile);
		if (swapDates == NULL)
			goto done;
	}

	/* Just to keep date/rate names the same. */
	swapRates = inRates;

	oneAlreadyAdded = FALSE;

	/* Add individual swap instruments  */
	for (i = 0; i < swapDates->numDates; i++)
	{
		/* Add those beyond stub zero curve
		 */
		if (zc->numItems > 0 && swapDates->adjusted[i] > zc->date[zc->numItems - 1])
		{
			/* Check if optimization okay. Note linear forwards
			 * not OK because they must include intermed. fwds.
			 */
			if (oneAlreadyAdded                                      &&
				discZC == NULL                                      &&
				swapRates[i - 1] != 0.0                              &&
				swapDates->adjusted[i - 1] == zc->date[zc->numItems - 1] &&
				swapDates->previous[i] == swapDates->original[i - 1] &&
				swapDates->onCycle[i] &&
				interpType != JPMCDS_LINEAR_FORWARDS)
			{
				/* Optimization: compute from last
				 */
				if (AddSwapFromPrevious
					(zc,
					swapDates->adjusted[i], swapRates[i],
					1.0,
					swapDates->adjusted[i - 1], swapRates[i - 1],
					1.0,
					fixDayCountConv) == FAILURE)
				{
					goto done;
				}
			}
			else             /* No efficiency */
			{
				if (JpmcdsZCAddSwap
					(zc, discZC,
					1.0,
					swapDates->original[i], swapDates->onCycle[i],
					swapRates[i],
					fixedSwapFreq, floatSwapFreq,
					fixDayCountConv, floatDayCountConv,
					interpType, interpData,
					badDayList,
					badDayAndStubPos, holidayFile) == FAILURE)
				{
					goto done;
				}
				oneAlreadyAdded = TRUE;
			}
		} /* if (zc->numItems > 0 && swapDates->adjusted[i] > zc->date[zc->numItems - 1]) */
		else {
			if (JpmcdsZCAddSwap
				(zc, discZC,
				1.0,
				swapDates->original[i], swapDates->onCycle[i],
				swapRates[i],
				fixedSwapFreq, floatSwapFreq,
				fixDayCountConv, floatDayCountConv,
				interpType, interpData,
				badDayList,
				badDayAndStubPos, holidayFile) == FAILURE)
			{
				goto done;
			}
			oneAlreadyAdded = TRUE;
		} /* else */
	}/* for (i=0;  i < swapDates->numDates; i++)  */

   status = SUCCESS;
   
 done:
   JpmcdsSwapDatesFree(swapDates);
   FREE_ARRAY(newRates);

   if (status == FAILURE)
       JpmcdsErrMsg("%s: Failed.\n", routine);

   return status;
}


/*
***************************************************************************
** Optimization, adds a swap to a zero curve assuming the one added just 
** before has the same coupon dates (w/ one additional payment).  This is
** much quicker than generating a cash-flow-list for the swap and adding those.
** And it often happens that this routine can be used.
**
** ZCurve is updated to contain swap information.
**
** How it works:
**    Assuming present value of principal + coupon payments == current price,
**        price = discount[n] + sum of coupon[i]*discount[i]
**    replacing for coupon:
**        coupon[i] = yearFraction(date[i] - date[i-1]) * couponRate
**        price = discount[n] + sum of discount[i]*yearFraction(i)*couponRate
**        price = discount[n] + couponRate * (sum of discount[i]*yearFraction[i])
**
**  We assume that the new swap shares the same yearFraction[], with previous.
**  Also, all discounts except the last (which we're trying to find) are known.
**
**  So, from the previous swap we back out sum of discount[i]*yearFraction[i],
**  which is done via:  sumDY = (priceOld - discount[n]) / couponRateOld
**
**  So, 
**        priceNew = discount[n] + couponRate(discount[n]*yearFraction[n]+sumDY)
**        priceNew = discount[n] * (1 + couponRate*yearFraction[n]) + cRate*sumDY
**  Solving for discount[n]:
**        discount[n] = (priceNew - couponRate*sumDY) / (1+couponRate*yearF[n])
**
** Assumptions:
**   Naturally it assumed that the new swap is equal to the old one with 
**   an additional coupon.  Also assume that the previous couponRate!=0.0,
**   and that the last discount factor in the zero curve corresponds to the date
**   of the previous swap and was placed there based on the previous swap.
***************************************************************************
*/
static int AddSwapFromPrevious(           /* adds swap based on one just added */
   ZCurve                      *zc,        /* (M) ZCurve to add to */
   TDate                        dateNew,   /* (I) date of new swap */
   double                       rateNew,   /* (I) rate of new swap */
   double                       priceNew,  /* (I) price of new swap */
   TDate                        dateOld,   /* (I) date of previous swap */
   double                       rateOld,   /* (I) rate of previous swap */
   double                       priceOld,  /* (I) price of previous swap */
   long                         dayCount)  /* (I) day counting convention */
{
    static char routine[]="AddSwapFromPrevious";
    int         status = FAILURE;       /* Until successfull  */
    double      sumDY;                  /* sum of discount[i]*yearFrac[i] */
    double      yf;                     /* year fraction dateNew-dateOld */
    double      discount;               /* discount calculated from swapNew */
    double      divisor;                /* a factor used to divide */

    /* Compute year fraction  */
    if (JpmcdsDayCountFraction(dateOld, dateNew, dayCount, &yf) == FAILURE)
        goto done;

    /* Compute divisor, and make sure it isn't = 0. */
    divisor = 1.0 + rateNew*yf;
    if (IS_ALMOST_ZERO(divisor))
    {
        JpmcdsErrMsg("%s: Rate (%f) implies zero discount factor.\n", routine, rateNew);
        goto done;
    }

    /* Here we assume that calling routine has already checked that
     * rateOld != 0.
     */
    sumDY = (priceOld - zc->discount[zc->numItems-1]) / rateOld;

    discount = (priceNew - rateNew*sumDY) / divisor;
   
    if (discount <= 0.0)
    {
        /*
         * The discount may be <= 0. This is the classical case where a 
         * consistent zero curve cannot be built.
         */
        JpmcdsErrMsg("%s: Implied discount factor (%f) on %s is <= 0.0.\n", 
                  routine, discount, JpmcdsFormatDate(dateNew));
        JpmcdsErrMsg("%s: The swap rates may be inconsistent with "
                  "one another.\n", routine);
        goto done;
    }

    /* Add the discount factor to the zero curve. */
    if (JpmcdsZCAddDiscountFactor (zc, dateNew, discount) == FAILURE)
        goto done;

    status = SUCCESS;

  done:
    if (status == FAILURE)
    {
       JpmcdsErrMsg("%s: Failed for swap at %s(adj), rate=%f price=%f.\n",
                 routine, JpmcdsFormatDate(dateNew), rateNew, priceNew);
    }
    return status;
}


/*
***************************************************************************
** Adds a single swap instrument to a ZCurve.
**
** ZCurve is updated to contain swap information (with points at each of 
** the coupon dates, as well as the maturity date).
***************************************************************************
*/
int JpmcdsZCAddSwap(                             /* adds a swap to ZCurve */
   ZCurve              *zc,                /* (M) ZCurve to add to          */
   TCurve              *discZC,            /* (I) ZC used for discounting */
   double              price,              /* (I) Par price: usually 1.0    */
   TDate               matDate,            /* (I) Unadjusted Maturity date  */
   TBoolean            onCycle,            /* (I) If matDate on cycle       */
   double              rate,               /* (I) Coupon rate (0.06 for 6%) */
   int                 fixedSwapFreq,      /* (I) Fixed coupon frequency    */
   int                 floatSwapFreq,      /* (I) Floating coupon frequency */
   long                fixDayCountConv,    /* (I) Fixed leg convention      */
   long                floatDayCountConv,  /* (I) Floating leg convention   */
   long                interpType,         /* (I) Interpolation method type */
   TInterpData         *interpData,        /* (I) Interpolation data        */
   TBadDayList         *badDayList,        /* (I) List of bad-business days */
   TBadDayAndStubPos   badDayAndStubPos,   /* (I) Bad day conv + stub pos   */
   char                *holidayFile)       /* (I) See JpmcdsBusinessDay        */
{
   static char    routine[] = "JpmcdsZCAddSwap";
   int            status = FAILURE;  
   TDate          adjMatDate;          /* Mat date adjusted to be good date */
   TCashFlowList *cfl = NULL;          /* Cash flow list representing swap*/
   TDateInterval  ivl;
   long           badDayConv;
   TStubPos       stubPos;

   if (JpmcdsBadDayAndStubPosSplit(badDayAndStubPos, &badDayConv, &stubPos) != SUCCESS)
   {
       goto done;
   }

   /* If floating side at par. */
   if (discZC == NULL)
   {
       TBoolean isEndStub;
       /* Here we must pass in the unadjusted matDate, so that the 
        * cashflow dates are correctly generated. 
        */
       if (JpmcdsFreq2TDateInterval (fixedSwapFreq, &ivl) != SUCCESS)
           goto done;

       if (onCycle)
       {
           isEndStub = TRUE;
       }
       else
       {
           if (JpmcdsIsEndStub(zc->valueDate,
                            matDate,
                            &ivl,
                            stubPos,
                            &isEndStub) != SUCCESS)
           {
               goto done;
           }
       }
       cfl = JpmcdsZCGetSwapCFL(zc->valueDate,
                                matDate, isEndStub,
                                rate,
                               &ivl,
                                fixDayCountConv,
                                badDayList, badDayConv, holidayFile);
       if (cfl==NULL)
           goto done;

       /* Adjust matDate to be a good business day, store it in adjMatDate.
        */
       if (JpmcdsZCAdjustDate(matDate, badDayList, badDayConv, holidayFile,
                           &adjMatDate) == FAILURE)
           goto done;


       /* Add rate implied by cashflow list to the zeroCurve.
        */
       if (JpmcdsZCAddCashFlowList(zc, cfl, price, adjMatDate,
                             interpType,interpData) == FAILURE)
           goto done;


       /* Add zero rates at all dates in cfl.
        * This guarantees that regardless of the interpolation used when
        * pricing one of the original benchmark swaps with the zero curve
        * created by this routine, the swap will be priced to par.
        * Since JpmcdsZCInterpolate does not support linear forwards,
        * we switch to linear interp here if linear fwds was specified.
        */
       if (JpmcdsZCAddCFLPoints
           (zc, cfl, 0, 
            interpType == JPMCDS_LINEAR_FORWARDS ? JPMCDS_LINEAR_INTERP:interpType, 
            interpData) == FAILURE)
           goto done;
   }
   else         /* Need to value floating side */
   {
       if (JpmcdsZCValueFixFltSwap
           (zc, discZC,
            price, matDate, rate,
            fixedSwapFreq, floatSwapFreq,
            fixDayCountConv, floatDayCountConv,
            interpType, interpData,
            badDayAndStubPos, holidayFile) == FAILURE)
       {
           goto done;
       }
   } /* ValueFloating */

   status = SUCCESS;

 done:
   JpmcdsFreeCFL(cfl);
   if (status == FAILURE)
   {
       JpmcdsErrMsg("%s: Failed for swap at %s(unadj), rate=%f price=%f.\n",
                 routine, JpmcdsFormatDate(matDate), rate, price);
   }

   return status;
}


/*
***************************************************************************
** Models the swap points which are to be added to the zero curve using the 
** more detailed swap structures.
***************************************************************************
*/
int JpmcdsZCValueFixFltSwap
(ZCurve           *zc,                 /* (M) ZCurve to add to           */
 TCurve           *discZC,             /* (I) Discount zero curve        */
 double            price,              /* (I) Par price: usually 1.0     */
 TDate             matDate,            /* (I) Maturity date              */
 double            rate,               /* (I) Coupon rate (0.06 for 6%)  */
 int               fixedSwapFreq,      /* (I) Fixed coupon frequency     */
 int               floatSwapFreq,      /* (I) Floating coupon frequency  */
 long              fixDayCountConv,    /* (I) Fixed leg convention       */
 long              floatDayCountConv,  /* (I) Floating leg convention    */
 long              interpType,         /* (I) Interpolation method       */
 TInterpData      *interpData,         /* (I) Interpolation data         */
 TBadDayAndStubPos badDayAndStubPos,   /* (I) See JpmcdsBusinessDay         */
 char             *holidayFile)       /* (I) See JpmcdsBusinessDay         */

{
   static char routine[] = "JpmcdsZCValueFixFltSwap";
   int                status = FAILURE;

   TCouponDateList    *floatDL = NULL;
   TCouponDateList    *fixedDL = NULL;
   TStreamFloat       *floatStream = NULL;
   TStreamFixed       *fixedStream = NULL;
   TDateInterval      floatInterval;
   TDateInterval      fixedInterval;
   TDateInterval      payInterval;
   TStubRates         firstCoupon;
   TStubRates         finalCoupon;
   TFloatRate         floatingRate;
   TCashFlowList      *fixedCFL = NULL;
   TCashFlowList      *floatCFL = NULL;
   TCashFlowList      *swapCFL = NULL;
   TPfunctionParams   funcParams;
   TCurve             *tc = NULL;
   TDate              indexMatDate;
   TDate              minDate;          /* Min of pay and index mat date  */
   TDate              maxDate;          /* Max of pay and index mat date */
   long               badDayConv;
   TStubPos           stubPos;
   TBoolean           isEndStub;

   if (JpmcdsBadDayAndStubPosSplit(badDayAndStubPos,
                                &badDayConv,
                                &stubPos) != SUCCESS)
   {
      goto done;
   }

   /* Convert the frequency to date intervals */
   if (JpmcdsFreq2TDateInterval((long) fixedSwapFreq, &fixedInterval)==FAILURE)
   {
       PROGRAM_BUG();
       goto done;
   }
   if (JpmcdsFreq2TDateInterval((long) floatSwapFreq, &floatInterval)==FAILURE)
   {
       PROGRAM_BUG();
       goto done;
   }

   /* make sure that the pay interval that is contained in the floating rate
      is bigger than the maturity interval or otherwise JpmcdsZerosToCouponsPoint
      will not calculate a basic *IBOR rate */
   if (JpmcdsMakeDateInterval(2,'A',&payInterval)==FAILURE)
   {
       PROGRAM_BUG();
       goto done;
   }

   /* set up the stub and floating rate structures */
   firstCoupon.stubRate = 0.0;
   firstCoupon.stubInterval = floatInterval;

   finalCoupon.stubRate = 0.0;
   finalCoupon.stubInterval = floatInterval;

   floatingRate.matInterval = floatInterval;
   floatingRate.payInterval = payInterval;
   floatingRate.dayCountConv = floatDayCountConv;
   JPMCDS_SET_ADJ_INTERVAL_DAYS(floatingRate.spotOffset, RESET_OFFSET)
   floatingRate.weight = 1.0;
   floatingRate.spread = 0.0;


   if (JpmcdsIsEndStub(zc->valueDate,
                    matDate,
                    &fixedInterval,
                    stubPos,
                    &isEndStub) != SUCCESS)
   {
       goto done;
   } 
   fixedDL = JpmcdsNewCouponDatesSwap(
                zc->valueDate,
                matDate,
                &fixedInterval,
                ADJUST_LAST_ACC_DATE,
                IN_ARREARS,
                PAY_OFFSET,
                RESET_OFFSET,
                isEndStub,
                LONG_STUB,
                FIRST_ROLL_DATE,
                LAST_ROLL_DATE,
                FULL_FIRST_COUPON_DATE,
                badDayConv,
                badDayConv,
                badDayConv,
                holidayFile);
   if (fixedDL == NULL)
       goto done;


   /* Make notional of fixed side = +1.
    */
   fixedStream = JpmcdsNewStreamFixed(
                   fixedDL,
                   JPMCDS_SINGLE_REFIX,
                   1.0,
                   fixDayCountConv, 
                   0,
                   rate,
                   rate,
                   rate);
   if (fixedStream == NULL)
       goto done;

   fixedCFL = JpmcdsNewStreamFixedCFL(zc->valueDate, fixedStream);
   if (fixedCFL == NULL)
       goto done;

   if (JpmcdsIsEndStub(zc->valueDate,
                    matDate,
                    &floatInterval,
                    stubPos,
                    &isEndStub) != SUCCESS)
   {
       goto done;
   } 
   floatDL = JpmcdsNewCouponDatesSwap(
                zc->valueDate,
                matDate,
                &floatInterval,
                ADJUST_LAST_ACC_DATE,
                IN_ARREARS,
                PAY_OFFSET,
                RESET_OFFSET,
                isEndStub,
                LONG_STUB,
                FIRST_ROLL_DATE,
                LAST_ROLL_DATE,
                FULL_FIRST_COUPON_DATE,
                badDayConv,
                badDayConv,
                badDayConv,
                holidayFile);
   if (floatDL == NULL)
       goto done;


   /* Make notional of floating side -1.
    */
   floatStream = JpmcdsNewStreamFloat(
                   floatDL,
                   JPMCDS_SINGLE_REFIX,
                   -1.0,
                   floatDayCountConv,
                   0,
                   &firstCoupon,
                   &finalCoupon,
                   0.0,         /* Average so far */
                   1.0,         /* Compounded so far */
                   &floatingRate);

   if (floatStream == NULL)
       goto done;

   /* find the furthest reaching date */

   if (JpmcdsDtFwdAny
       (floatDL->fArray[floatDL->fNumItems-1].resetDate,
        &floatInterval,
        &indexMatDate) == FAILURE)
       goto done;

   if (JpmcdsBusinessDay(indexMatDate,
                      badDayConv,
                      holidayFile,
                      &indexMatDate) == FAILURE)
       goto done;


   /* Pricing the floating side depends on both the payDate and the
    * the index maturity date. We want the max of these two, so
    * that when we interpolate the other, the answer is not changed
    * by adding the zero for the *next* swap. We put *rate* in for
    * now, and then replace this with the correct value below. 
    */
   maxDate = MAX(indexMatDate, floatDL->fArray[floatDL->fNumItems-1].payDate);
   minDate = MIN(indexMatDate, floatDL->fArray[floatDL->fNumItems-1].payDate);
   if (JpmcdsZCAddRate(zc,maxDate,rate) == FAILURE)
       goto done;


   /* Need to make a temporary TCurve copy of the ZCurve, since all the 
    * swap routines only take TCurve structures 
    */

   tc = JpmcdsZCToTCurve(zc);     
   if (tc == NULL)
       goto done;
   
   funcParams.fixedCfls = fixedCFL;
   funcParams.zeroCurve = tc;
   funcParams.floatStream = floatStream;
   funcParams.interpType = interpType;
   funcParams.price = price;
   funcParams.rateBadDayConv = badDayConv;
   funcParams.holidayFile = holidayFile;
   funcParams.discZC = discZC;

   /* Use the par swap rate as the first zero rate guess. Use that
    * to set the bounds. Use the fact that JpmcdsRootFind starts iterating
    * at average between rlow and rhigh. Set up rlow and rhigh so that
    * that average is "rate".
    */
   if (JpmcdsRootFindBrent(JpmcdsObjFuncPVtheSwap, 
                        &funcParams,
                        LOWER_BOUND, UPPER_BOUND,
                        MAX_ITERATIONS, rate,
                        INITIAL_X_STEP, INITIAL_F_DERIV,
                        X_TOLERANCE, F_TOLERANCE,
                        &zc->rate[zc->numItems-1]) == FAILURE)
   {
       JpmcdsErrMsg("%s: Root finder failed.\n", routine);
       goto done;
   }
   
   if (JpmcdsZCComputeDiscount(zc,
                         zc->date[zc->numItems-1],
                         zc->rate[zc->numItems-1],
                         &zc->discount[zc->numItems-1])==FAILURE)
   {
       goto done;
   }

   /* We are only generating the cashflows for the dates-we don't
    * care about the values of the cashflows. We use the dates only
    * to insert points in the zero curve so that the benchmarks will
    * price out to par. Note that there could be small errors if
    * we are adjusting for bad days, and if the interpolation method
    * used by the pricing software does not match the zero curve 
    * interpolation method, because the index maturity dates do not
    * necessarily match the pay dates.
    */
   floatCFL = JpmcdsNewStreamFloatCFLGen(tc, interpType, floatStream, badDayConv, holidayFile);
   if (floatCFL == NULL)
       goto done;

   swapCFL = JpmcdsMergeCFL(fixedCFL, floatCFL);
   if (swapCFL == NULL)
       goto done;


   /* Add zero rates at all dates in swapCFL.
    * This guarantees that regardless of the interpolation used when
    * pricing one of the original benchmark swaps with the zero curve
    * created by this routine, the swap will be priced to par.
    */
   if (JpmcdsZCAddCFLPoints(zc,
                      swapCFL,
                      0,     /* Add point after this one*/
                      interpType,
                      interpData) == FAILURE)
   {
       goto done;
   }

   /* If we are adjusting for bad days, the index maturity date 
    * is not necessarily the same as the pay date. Add another
    * zero point to help guarantee that the benchmarks price out to par.
    */
   if (minDate != maxDate)
   {
       double zeroRate;

       if (JpmcdsZCInterpolate(zc, minDate, interpType, interpData, &zeroRate) == FAILURE)
           goto  done;

       if (JpmcdsZCAddRate(zc, minDate, zeroRate) == FAILURE)
           goto done;
   }

   status = SUCCESS;

 done:
   if (status == FAILURE)
   {
       JpmcdsErrMsg("%s: Failed.\n", routine);
   }
   
   JpmcdsFreeCouponDateList(floatDL);
   JpmcdsFreeCouponDateList(fixedDL);
   JpmcdsFreeTStreamFloat(floatStream);
   JpmcdsFreeTStreamFixed(fixedStream);
   JpmcdsFreeCFL(fixedCFL);
   JpmcdsFreeCFL(floatCFL);
   JpmcdsFreeCFL(swapCFL);
   JpmcdsFreeTCurve(tc);
   return status;
}


/*
***************************************************************************
** Objective function for root-finder to use.
***************************************************************************
*/
static int JpmcdsObjFuncPVtheSwap(
        double    rateGuess,          /* (I) variable x - guess of rate     */
        void      *p,                 /* (I) data ptr (TPfunctionParams)  */
        double    *f)                 /* (O) f(x) == value to make zero   */
{
   TPfunctionParams *params = ( TPfunctionParams * ) p;

   params->zeroCurve->fArray[params->zeroCurve->fNumItems-1].fRate = rateGuess;

   if (CalcPV(
              params->zeroCurve,
              params->interpType,
              params->fixedCfls,
              params->floatStream,
              params->price,
              params->rateBadDayConv,
              params->holidayFile,
              params->discZC,
              f) == FAILURE)
   {
       *f = -1.0;
       return FAILURE;
   }
   else
       return SUCCESS;
}


/*
***************************************************************************
** Calculates the PV of a swap.
***************************************************************************
*/
static int CalcPV
    (TCurve           *zeroCurve,       /* (I) Zero curve for estim cfls  */ 
     long              interpType,      /* (I) Type of zero interpolation */
     TCashFlowList    *fixedFlows,      /* (I) CFLs of fixed leg          */
     TStreamFloat     *floatStream,     /* (I) Floating leg               */
     double            price,           /* (I) Price - 1.0==par           */
     long              rateBadDayConv,  /* (I) Bad day convention         */
     char             *holidayFile,     /* (I) Name of holiday file       */
     TCurve           *discZC,          /* (I) Discount zero curve        */
     double           *dealPV)          /* (O) PV of the swap             */
{
    static char routine[]="CalcPV";
    int status = FAILURE;

    TCashFlowList *floatFlows = NULL; 
    double floatPV = 0.0;
    double fixedPV = 0.0;

    floatFlows = JpmcdsNewStreamFloatCFLGen (zeroCurve, interpType, floatStream, rateBadDayConv, holidayFile);
    if (floatFlows == NULL)
        goto done;

    if(discZC == NULL)
    {
        discZC = zeroCurve;
    }

    if (JpmcdsCashFlowPV(floatFlows, discZC, interpType, &floatPV) == FAILURE)
        goto done;

    if (JpmcdsCashFlowPV(fixedFlows, discZC, interpType, &fixedPV) == FAILURE)
        goto done;

    *dealPV =  fixedPV + floatPV + price - 1.;


    status = SUCCESS;

  done:
    if (status == FAILURE)
        JpmcdsErrMsg("%s: Failed.\n", routine);

    JpmcdsFreeCFL(floatFlows);
    return status;
}
