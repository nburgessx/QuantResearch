/*
 * ISDA CDS Standard Model
 *
 * Copyright (C) 2009 International Swaps and Derivatives Association, Inc.
 * Developed and supported in collaboration with Markit
 * 
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the ISDA CDS Standard Model Public License.
 */

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h> 
#include "badday.h"
#include "cerror.h"
#include "cmemory.h"
#include "rtbrent.h"
#include "convert.h"
#include "datelist.h"
#include "date_sup.h"
#include "ldate.h"
#include "lintrp.h"
#include "interp.h"
#include "macros.h"
#include "bsearch.h"
#include "tcurve.h"
#include "yearfrac.h"
#include "zcprvt.h"


#define JpmcdsBasisValid(x)  ((x)>=1 && (x)<=5001 ? SUCCESS : FAILURE)

#define MIN_ZERO_RATE  -10000.00       /* minimum zero rate to place in curve*/
#define MAX_ZERO_RATE   10000.00       /* maximum zero rate to place in curve*/

/* Define constants for JpmcdsRootFindBrent */
#define INITIAL_X_STEP          (1E-8)
#define INITIAL_F_DERIV         (0)
#define X_TOLERANCE             (DBL_MAX)
#define F_TOLERANCE             (1E-10)
#define MAX_ITERATIONS          (50)
#define LOWER_BOUND             (-DBL_MAX/1000)
#define UPPER_BOUND             (DBL_MAX/1000)

typedef struct {                           /* data of objective-functions */
   ZCurve         *zc;                     /* ZCurve being fit */
   TCashFlowList  *cfl;                    /* cash flows specifying fit */
   int             firstUncovered;         /* index of first unconvered c.f. */
   double          pvUnCovered;            /* PV of uncovered c.f.s*/
   long            interpType;             /* interpolation type */
   TInterpData    *interpData;             /* data for interpolation method */
   int             zcIndex;                /* index in ZC of rate being fit */
   double         *discount;               /* lin fwds only: array of discnt*/
   double          offset;                 /*  ""   "" : const for linear eqn*/
   double          startingDiscount;       /*  ""   "" : discount to 1st fwd */
   TDateList      *fwdDL;                  /*  ""   "" : list of forward dts*/
} TObjectiveData;


/*
***************************************************************************
** Objective function for use with JpmcdsRootFindBrent.  
**
** Basically, JpmcdsZCAddCashFlowList() uses this function to iteratively solve
** for a rate to add to the zero curve, so that when the uncovered cash 
** flows are discounted (note: they are now covered by the curve), their 
** present value == the instrument price - npv of covered cash flows. 
***************************************************************************
*/
static int JpmcdsObjFunctionRate(            /* objective func for root finder */
   double              rate,               /* (I) newest rate guess */
   TObjectiveData     *data,               /* (I) representing problem solved*/
   double             *result);            /* (O) objective function value */


/*
***************************************************************************
** Constructs a TCurve from a ZCurve.
***************************************************************************
*/
TCurve* JpmcdsZCToTCurve(        /* construct TCurve from a ZCurve */
   ZCurve  *zc)                            /* (I) ZCurve to copy */
{
   static char routine[] = "JpmcdsZCToTCurve";

   TCurve *tc = JpmcdsMakeTCurve(zc->valueDate, /* make TCurve copy of ZCurve */
                          zc->date,
                          zc->rate,
                          zc->numItems,
                          (double)zc->basis,
                          zc->dayCountConv);
   if (tc == NULL)
       JpmcdsErrMsg("%s: couldn't make TCurve copy of zero curve\n",routine);

   return tc;
}


/*
***************************************************************************
** Constructs a ZCurve from a TCurve.
***************************************************************************
*/
ZCurve* JpmcdsZCFromTCurve(        /* construct ZCurve from a TCurve */
   TCurve  *tc)                            /* (I) TCurve to copy */
{
   static char routine[] = "JpmcdsZCFromTCurve";

   ZCurve *zc = NULL;

   zc = JpmcdsZCMake( tc->fBaseDate,
                tc->fNumItems + 32, /* Allow room for expansion */
                (long) tc->fBasis,
                tc->fDayCountConv);

   if (zc!=NULL)
   {
       int     i;                      /* loops over elements in TCurve */
       for (i=0; i<tc->fNumItems; i++)
       {
           if (JpmcdsZCAddRate(zc,
                            tc->fArray[i].fDate,
                            tc->fArray[i].fRate)==FAILURE)
           {
               JpmcdsZCFree(zc);
               JpmcdsErrMsg("%s: couldn't add rate\n",routine);
               return NULL;
           }
       }
   }

   return zc;
}
 

/*
***************************************************************************
** Constructs a ZCurve, which is basically an array of <rate,date> pairs, 
** where the rates are all zero-coupons with the same compounding and day 
** count convention.  This allows a bunch of different types of interest rate
** instruments to be combined into one yield curve in the same terms.
**
***************************************************************************
*/
ZCurve* JpmcdsZCMake(              /* constructor for ZCurve */
   TDate           valueDate,              /* (I) base date for rates */
   int             n,                      /* (I) expected number of rates */
   long            basis,                  /* (I) compounding period */
   long            dayCountConv)           /* (I) day count convention */
{
   ZCurve *zc;                         /* new ZCurve */
   double *rates;                      /* allocated rates      array for ZC */
   double *disc;                       /* allocated discounts  array for ZC */
   TDate  *dates;                      /* allocated dates      array for ZC */

   if (n<=0 || 
       JpmcdsBasisValid(basis)==FAILURE || 
       JpmcdsDayCountValid("JpmcdsZCMake",dayCountConv)==FAILURE)
   {
       JpmcdsErrMsg("JpmcdsZCMake: bad input: n:%d  basis:%ld  dayCountConv:%s\n",
                   n,
                   basis,
                   JpmcdsFormatDayCountConv(dayCountConv));
       return NULL;
   }

   zc =    NEW(ZCurve);
   rates = NEW_ARRAY(double,n);
   disc  = NEW_ARRAY(double,n);
   dates = NEW_ARRAY(TDate,n);

   if (zc==NULL || rates==NULL || dates==NULL || disc==NULL)
   {
       FREE(zc);
       FREE_ARRAY(rates);
       FREE_ARRAY(disc);
       FREE_ARRAY(dates);
       JpmcdsErrMsg("JpmcdsZCMake: out of memory\n");
       return NULL;
   }

   zc->valueDate    = valueDate;       /* copy inputs to new structure */
   zc->rate         = rates;
   zc->date         = dates;
   zc->discount     = disc;
   zc->numItems     = 0;
   zc->numAlloc     = n;
   zc->basis        = basis;
   zc->dayCountConv = dayCountConv;

   return zc;
}


/*
***************************************************************************
** Frees memory associated with a zero curve.
***************************************************************************
*/
void JpmcdsZCFree(              /* frees a zero curve */
   ZCurve                 *zc)             /* zero curve to free */
{
   if (zc != NULL)
   {
       FREE_ARRAY(zc->rate);
       FREE_ARRAY(zc->discount);
       FREE_ARRAY(zc->date);
       FREE(zc);
   }
}


/*
***************************************************************************
** Adds points to a zero curve from another zero curve, but only those 
** dates before first date of the other curve.
**
** ZCurve is updated to include the new points.
**
** Note: curve value dates MUST be the same, basis and dayCounts may differ.
***************************************************************************
*/
int JpmcdsZCAddPrefixCurve(        /* adds prefix from other ZCurve */
   ZCurve                     *zc1,        /* (M) zCurve to added to */
   ZCurve                     *zc2)        /* (I) zCurve to add from */
{
   TDate   firstDate;                  /* first date of zc */
   int     i;                          /* loops over zc2 elements */

   if (zc1->numItems<=0)               /* if nothing in zCurve, add everything*/
   {
       if (zc2->numItems<=0)           /* easy to add nothing to nothing */
       {
           return SUCCESS;
       }
       firstDate = zc2->date[zc2->numItems-1];
   }
   else
   {
       firstDate = zc1->date[0];
   }

   for (i=0; i<zc2->numItems && zc2->date[i]<firstDate; i++)
   {
       if (JpmcdsZCAddGenRate(zc1,
                           zc2->date[i],
                           zc2->rate[i],
                           zc2->basis,
                           zc2->dayCountConv)==FAILURE)
       {
           JpmcdsErrMsg("ZCAddPrefix Error adding rate to curve\n");
           return FAILURE;
       }
   }

   return SUCCESS;
}


/*
***************************************************************************
** Adds points to a zero curve from another zero curve, but only those 
** dates after last date of the other curve.
**
** ZCurve is updated to include the new points.
**
** Note: curve value dates MUST be the same, basis and dayCounts may differ.
***************************************************************************
*/
int JpmcdsZCAddSuffixCurve(        /* adds suffix from other ZCurve */
   ZCurve                     *zc1,        /* (M) zCurve to added to */
   ZCurve                     *zc2)        /* (I) zCurve to add from */
{
   TDate   lastDate;                   /* last date of zc */
   int     i;                          /* loops over zc2 elements */

   if (zc1->numItems<=0)               /* if nothing in zCurve, add everything*/
   {
       if (zc2->numItems<=0)           /* easy to add nothing to nothing */
       {
           return SUCCESS;
       }
       lastDate = zc2->date[0] - 1;
   }
   else
   {
       lastDate = zc1->date[zc1->numItems-1];
   }

   for (i=zc2->numItems-1; i>=0 && zc2->date[i]>lastDate; i--)
   {
       if (JpmcdsZCAddGenRate(zc1,
                           zc2->date[i],
                           zc2->rate[i],
                           zc2->basis,
                           zc2->dayCountConv)==FAILURE)
       {
           JpmcdsErrMsg("ZCAddSuffix Error adding rate to curve\n");
           return FAILURE;
       }
   }

   return SUCCESS;
}


/*
***************************************************************************
** Adds simple-interest money market bond to ZCurve
**
** ZCurve is updated to include information from money market instruments
***************************************************************************
*/
int JpmcdsZCAddMoneyMarket(        /* adds simple-interest MMKT to ZC */
   ZCurve                 *zc,             /* (M) ZCurve to add to */
   TDate                   dates[],        /* (I) maturity dates of instrs */
   double                  rates[],        /* (I) rates (e.g. 0.06 for 6%) */
   int                     n,              /* (I) number of instruments */
   long                    dayCountConv)   /* (I) day-count-convention */
{
   int     i;                          /* loops over instruments */

   if (JpmcdsDayCountValid("JpmcdsZCAddMoneyMarket",dayCountConv)==FAILURE)
       return FAILURE;

   for (i=0; i<n; i++)
   {
       if (JpmcdsZCAddGenRate(zc,
                           dates[i],
                           rates[i],
                           JPMCDS_SIMPLE_BASIS,
                           dayCountConv)==FAILURE)
       {
           JpmcdsErrMsg("JpmcdsZCAddMoneyMarket: error adding rate[%d]: %6.4f%% %s\n",
                       i,rates[i],JpmcdsFormatDate(dates[i]));
           return FAILURE;
       }
   }

   return SUCCESS;
}


/*
***************************************************************************
** Adds information represented by a list-of-cash-flows to a zero curve.  
** Any cash flows which are already covered by the zero curve are 
** discounted at rates derived from the zero curve.  Cash flows beyond the
** zero curve imply discount factors, which are added to the zero curve.  If
** there is more than one such cash flow, several points are added to the curve,
** which are calculated by using an interative root-finding secant method,
** where the discount factor for the last cash flow is guessed (and the other
** discount factors are implied by interpolation) where the current price = 
** net present value of all the cash flows.
**
** The zero curve is updated to reflect the cash flows.  A point is added for 
** every cash flow, if not already in ZCurve list.  For linear forwards, all 
** interpolated points are returned, e.g. for 1 month forwards in an annual 
** market, monthly points will be returned, not just yearly. 
**
** Notes: date may be set for non-linear-forward interpolation methods to a 
** date to be added to the zero curve.  This allows production of a curve 
** with "nice" dates.
***************************************************************************
*/
int JpmcdsZCAddCashFlowList(       /* adds cashFlowList info to a ZCurve */
   ZCurve                 *zc,          /* (M) ZCurve to add info to */
   TCashFlowList          *cfl,         /* (I) Cash flows to add to ZCurve */
   double                  price,       /* (I) Current price of cash-flows */
   TDate                   date,        /* (I) Date to add last point at */
   long                    interpType,  /* (I) Interpolation method's type */
   TInterpData            *interpData)  /* (I) Interpolation method's data */
{
    static char routine[]="JpmcdsZCAddCashFlowList";
    int         status = FAILURE;       /* Until proven successful */
    double      sumNPV = 0.0;           /* sum of n.p.vs of uncovered c.f.s */
    int         firstUncovered;         /* index in cfl of 1st uncovered c.f.*/

    if (date==0L)                       /* add at last c.f. if not set up */
    {
        date = cfl->fArray[ cfl->fNumItems-1].fDate;
    }

    if (zc->numItems<=0)                /* if no points in ZCurve ... */
    {
        firstUncovered = 0;
        /* if no data, default to lin interp*/
        interpType = JPMCDS_LINEAR_INTERP;
    }
    else
    {
        TDate lastZCDate = zc->date[ zc->numItems-1 ]; /* last date in ZCurve*/

        if (date<=lastZCDate)
        {
            JpmcdsErrMsg("%s: date to add already covered.\n", routine);
            goto done;
        }

        firstUncovered = cfl->fNumItems - 1;
        if (firstUncovered<0 || 
            cfl->fArray[firstUncovered].fDate < lastZCDate)
        {
            JpmcdsErrMsg("%s: no cash flows in list beyond"
                      " ZeroCurve.  Therefore, nothing to add.\n", routine);
            goto done;
        }
        while (firstUncovered>=0 && 
               cfl->fArray[firstUncovered].fDate>lastZCDate)
        {
            firstUncovered--;           /* decrement until covered */
        }
        firstUncovered++;               /* move up to first not covered */
    }


    if (firstUncovered > 0)             /* calc n.p.v. of covered cash flows */
    {
        if (JpmcdsJpmcdsZCPresentValueCFL( zc, cfl, 0, firstUncovered-1, interpType, interpData, &sumNPV) == FAILURE)
            goto done;
    }

    if (firstUncovered == cfl->fNumItems-1 /* if only one c.f. not covered */
        && cfl->fArray[firstUncovered].fDate==date) /* and c.f. is date to add */
    {
        double discFactor;              /* discount implied by last c.f. */
        double pvOfLast = price - sumNPV; /* present value of last c.f. */
        double futOfLast;               /* future value of last c.f. */
        futOfLast = cfl->fArray[firstUncovered].fAmount;

        if (futOfLast==0.0 || pvOfLast == 0.) 
        {
            JpmcdsErrMsg("%s: Last cash flow has zero value.\n", routine);
            goto done;
        }
        discFactor = pvOfLast / futOfLast;
        if (JpmcdsZCAddDiscountFactor(zc, date, discFactor) == FAILURE) 
            goto done;
    }
    else if(firstUncovered < cfl->fNumItems)
    {                                  
        double          rate;       /* rate to use for guess (and result)*/
        TObjectiveData  objData;    /* data for objective function */
       
        if (zc->numItems<=0)        /* make a guess based on zc */
        {
            rate = 0.06;
        } 
        else 
        {                           /* linearly extrapolate guess */
            if (JpmcdsZCInterpolate(zc, date, JPMCDS_LINEAR_INTERP, NULL, &rate) == FAILURE)
                rate = 0.06;
            rate = MAX(0.01, MIN(1.00, rate));
        }
        if (JpmcdsZCAddRate(zc, date, rate) == FAILURE)
            goto done;

        objData.zc               = zc; /* zero curve */
        objData.cfl              = cfl; /* cash flow list */
        objData.firstUncovered   = firstUncovered; /* 1st uncovered c.f.*/
        objData.pvUnCovered      = price - sumNPV; /* p.v.of uncoverd cfs*/
        objData.interpType       = interpType; /* interpolation type */
        objData.interpData       = interpData; /* interpolation data */
        objData.zcIndex          = zc->numItems-1; /*index into ZC ofrate*/
        objData.discount         = NULL; /* not used for interp*/
        objData.offset           = 0.0; /* ditto */
        objData.startingDiscount = 0.0; /* ditto */
        objData.fwdDL            = NULL; /* ditto */
       
        if (JpmcdsRootFindBrent((TObjectFunc) JpmcdsObjFunctionRate,
                             &objData,
                             MIN_ZERO_RATE, MAX_ZERO_RATE,
                             MAX_ITERATIONS, rate,
                             INITIAL_X_STEP, INITIAL_F_DERIV,
                             X_TOLERANCE, F_TOLERANCE,
                             &rate) == FAILURE)
        {
            goto done;
        }
        zc->rate[ objData.zcIndex ] = rate;
    }

    /* Ensure all dates used to discount cash flows are placed in zero curve.
     * This will make sure they can be priced appropriately by the curve.  Note
     * that non-linear interpolation might actually change when points are 
     * added, so they are calculated before any are added, so they
     * won't change. 
     */
    if (firstUncovered < cfl->fNumItems)
    {
        if(JpmcdsZCAddCFLPoints(zc,
                             cfl,
                             cfl->fArray[firstUncovered].fDate-1,
                             interpType,
                             interpData) == FAILURE)
        {
            goto done;
        }
    }

    status = SUCCESS;

  done:
    if (status == FAILURE)
        JpmcdsErrMsg("%s: Failed for CFL maturing %s.\n", routine, JpmcdsFormatDate(cfl->fArray[cfl->fNumItems-1].fDate));

    return status;
}


/*
***************************************************************************
** Adds points into a zero curve to match all points in a cashflow list 
** after a certain date. Note that only the *dates* in the cashflow list
** are used.
**
** Zero curve is modified.
***************************************************************************
*/
int  JpmcdsZCAddCFLPoints(          /* adds points in a CFL to zcurve */
   ZCurve                 *zc,          /* (M) ZCurve to add info to */
   TCashFlowList          *cfl,         /* (I) Use dates only for pts to add */
   TDate                   lastDate,    /* (I) Only add pts after this date*/
   long                    interpType,  /* (I) Interpolation method's type */
   TInterpData            *interpData)  /* (I) Interpolation method's data */
{
   static char routine[] = "JpmcdsZCAddCFLPoints";
   int     status = FAILURE;           /* status to return */
   int     n = 0;                      /* counts added points */
   int     i;                          /* iterates over cashflows */
   TDate  *dates = NEW_ARRAY(TDate,  cfl->fNumItems);
   double *rates = NEW_ARRAY(double, cfl->fNumItems);
   if (dates == NULL || rates == NULL)
       goto done;

   /* Ensure all uncovered cash flow dates in curve. Note that we perform
    * all interpolations *first*, before inserting any points back in
    * curve, because some interp methods are such that new points will
    * affect the values of remaining points to be interpolated.
    */
   for (i=0; i<cfl->fNumItems; i++)
   {
       TDate d = cfl->fArray[i].fDate;
       if (lastDate < d && 
           JpmcdsZCFindDateExact(zc,d)<0) /* only add if not already in */
       {
           double rate;
           if (JpmcdsZCInterpolate(zc,
                             d,
                             interpType,
                             interpData,
                            &rate) == FAILURE)
           {
               goto done;
           }
           if (rate >= 0.0) 
           {
               dates[n] = d;
               rates[n++] = rate;
           }
       }
   }

   while (--n>=0)
   {
       if (JpmcdsZCAddRate(zc,dates[n],rates[n]) == FAILURE)
           goto done;
   }
   
   status = SUCCESS;

 done:
   FREE_ARRAY(dates);
   FREE_ARRAY(rates);
   if (status==FAILURE)
       JpmcdsErrMsg("%s: Failed.\n",routine);

   return status;
}


/*
***************************************************************************
** Objective function for use with JpmcdsRootFindBrent.  
**
** Basically, JpmcdsZCAddCashFlowList() uses this function to iteratively solve
** for a rate to add to the zero curve, so that when the uncovered cash 
** flows are discounted (note: they are now covered by the curve), their 
** present value == the instrument price - npv of covered cash flows. 
***************************************************************************
*/
static int JpmcdsObjFunctionRate( /* objective function for root finder */
   double              rate,    /* (I) newest rate guess */
   TObjectiveData     *data,    /* (I) representing problem solved */
   double             *result)  /* (O) objective function value */
{
   double      sumPV = 0.0;         /* net present value of uncovered c.fs */
   static char routine[]="JpmcdsObjFunctionRate";
   int         status = FAILURE; /* Until proven successful */


   if (rate<=MIN_ZERO_RATE || MAX_ZERO_RATE<=rate)
   {
       JpmcdsErrMsg("%s: Rate guess out of limits.\n", routine);
       goto done;
   }

   data->zc->rate[ data->zcIndex ] = rate; /* put rate in zero curve */
   if (JpmcdsZCComputeDiscount(data->zc,
                            data->zc->date[ data->zcIndex],
                            rate,
                           &data->zc->discount[ data->zcIndex ]) == FAILURE)
   {
       goto done;
   }

   if (JpmcdsJpmcdsZCPresentValueCFL(data->zc,
                            data->cfl,
                            data->firstUncovered,
                            data->cfl->fNumItems-1,
                            data->interpType,
                            data->interpData,
                           &sumPV) == FAILURE)
   {
       goto done;
   }

   status = SUCCESS;
   *result = data->pvUnCovered - sumPV;   

 done:
   if (status == FAILURE)
       JpmcdsErrMsg("%s: Failed.\n", routine);

   return status;
}


/*
***************************************************************************
** Adds a general zero rate to a ZCurve.
**
** ZCurve is updated to include the rate.
***************************************************************************
*/
int JpmcdsZCAddGenRate(            /* adds rate to ZCurve */
   ZCurve         *zc,                  /* (M) ZCurve to add rate to */
   TDate           date,                /* (I) date associated with rate */
   double          rate,                /* (I) rate, in ZCurve terms */
   long            basis,               /* (I) compounding basis of rate*/
   long            dayCountConv)        /* (I) day counting convention */
{
    static char routine[]="JpmcdsZCAddGenRate";
    int         status = FAILURE;       /* Until proved succesfull */
    double  df;                         /* discount factor assoc w/ rate */

   if (basis==zc->basis && dayCountConv==zc->dayCountConv)
   {
        if (JpmcdsZCAddRate(zc,date,rate) == FAILURE)
            goto done;
   }

   if (JpmcdsRateToDiscount(rate,         /* convert to discount factor */
                         zc->valueDate,
                         date,
                         dayCountConv,
                         basis,
                        &df) == FAILURE)
       goto done;


    if (JpmcdsZCAddDiscountFactor(zc,date,df) == FAILURE)
        goto done;

    status = SUCCESS;

  done:
    if (status == FAILURE)
        JpmcdsErrMsg("%s: Failed.\n", routine);

   return status;
}


/*
***************************************************************************
** Adds a discount factor (at a specified date) to a ZCurve.
**
** Modified ZCurve, and error indicator.
***************************************************************************
*/
int JpmcdsZCAddDiscountFactor(                /* add discount factor to ZCurve */
   ZCurve                  *zc,         /* (M) ZCurve to add discount to */
   TDate                    date,       /* (I) date assoc'd w/ disc factor */
   double                   disc)       /* (I) discount factor to add */
{
    static char routine[]="JpmcdsZCAddDiscountFactor";
    int         status = FAILURE;       /* Until proven successful */
    double      rate;


    /* Convert discount factor to rate 
     */
    if (JpmcdsDiscountToRate(disc, zc->valueDate,
                          date, zc->dayCountConv, zc->basis,
                          &rate) == FAILURE)
        goto done;

    status = JpmcdsZCAddRateAndDiscount(zc, date, rate, disc);

  done:
    if (status == FAILURE)
        JpmcdsErrMsg("%s: Failed.\n", routine);

   return status;
}


/*
***************************************************************************
** Adds a zero rate to a ZCurve.
** The input rate must be in the proper ZCurve terms (zc->basis and zc->dayCountConv).
**
** ZCurve is updated to include the rate.
***************************************************************************
*/
int JpmcdsZCAddRate(                          /* adds rate to ZCurve */
   ZCurve         *zc,                  /* (M) ZCurve to add rate to */
   TDate           date,                /* (I) date associated with rate */
   double          rate)                /* (I) rate, in ZCurve terms */
{
    static char routine[]="JpmcdsZCAddRate";
    int         status = FAILURE;
    double      discount;               /* discount associated with rate */

    if (JpmcdsZCComputeDiscount(zc,date,rate,&discount) == FAILURE)
        goto done;

    status = JpmcdsZCAddRateAndDiscount(zc, date, rate, discount);

  done:
    if (status == FAILURE)
        JpmcdsErrMsg("%s: Failed.\n", routine);

    return status;
}


/*
***************************************************************************
** Adds a zero rate and discount to a ZCurve.  The input rate must 
** The input rate must be in the proper ZCurve terms (zc->basis and zc->dayCountConv).
**
** ZCurve is updated to include the rate.
***************************************************************************
*/
int JpmcdsZCAddRateAndDiscount(                  /* adds rate and discount factor to ZC*/
   ZCurve         *zc,                     /* (M) ZCurve to add rate to */
   TDate           date,                   /* (I) date associated with rate */
   double          rate,                   /* (I) rate, in ZCurve terms */
   double          disc)                   /* (I) discount factor to add */
{
    static char routine[] = "JpmcdsZCAddRateAndDiscount";

    if (zc->numItems > 0 &&
        date <= zc->date[zc->numItems-1])
    {                             /* make sure date not already in list */
        int i = JpmcdsZCFindDateExact(zc,date);
        if (i>=0)
        {
            if (fabs(rate - zc->rate[i]) < 0.0000001) /* if same, no error */
            {
                return SUCCESS;
            } 
            else
            {
                JpmcdsErrMsg("%s: Date %s (rate:%6.4f%%) already in curve"
                          " at index %d\n (rate:%6.4f%%)\n",
                          routine,
                          JpmcdsFormatDate(date),
                          rate*100.0,
                          i,
                          zc->rate[i]);
                return FAILURE;
            }
        }
    }

    if (zc->numAlloc <= zc->numItems)
    { /* re-allocate arrays */
        int n = zc->numItems + 32;      /* number of items to allocate */
        TDate  *dateNew = NEW_ARRAY(TDate,n);
        double *rateNew = NEW_ARRAY(double,n);
        double *discNew = NEW_ARRAY(double,n);
        if (dateNew==NULL || rateNew==NULL || discNew==NULL)
        {
            FREE_ARRAY(dateNew);
            FREE_ARRAY(rateNew);
            FREE_ARRAY(discNew);
            JpmcdsErrMsg("%s: out of memory\n",routine);
            return FAILURE;
        }
        memcpy ((char*)dateNew,
                (char*)(zc->date),
                sizeof(TDate)*(zc->numItems));
        memcpy ((char*)rateNew,
                (char*)(zc->rate), 
                sizeof(double)*zc->numItems);
        memcpy ((char*)discNew, 
                (char*)(zc->discount),
                sizeof(double)*(zc->numItems));
        FREE_ARRAY(zc->date);           /* free old arrays */
        FREE_ARRAY(zc->rate);
        FREE_ARRAY(zc->discount);
        zc->date =     dateNew;         /* use new arrays */
        zc->rate =     rateNew;
        zc->discount = discNew;
    } /* insert new data point...*/
    if (zc->numItems==0 || zc->date[zc->numItems-1]<date)
    {
        zc->date[     zc->numItems] = date; /* insert at end */
        zc->rate[     zc->numItems] = rate;
        zc->discount[ zc->numItems] = disc;
    }
    else
    {
        int i = zc->numItems - 1;       /* loops backwards, until before date */

        while (i>=0)
        {
            if (zc->date[i]<date)       /* stop when before date */
            {
                break;
            } /* shift down */
            zc->date[     i+1] = zc->date[i];
            zc->rate[     i+1] = zc->rate[i];
            zc->discount[ i+1] = zc->discount[i];
            i--;
        }
        zc->date[     i+1] = date;
        zc->rate[     i+1] = rate;
        zc->discount[ i+1] = disc;
    }
    zc->numItems++;
    return SUCCESS;
}


/*
***************************************************************************
** Calculates a simple interest money market rate from a zero curve.
***************************************************************************
*/
int JpmcdsZCMoneyMarketRate(                     /* gets money market rate from ZCurve */
   ZCurve                 *zc,             /* (I) ZCurve */
   TDate                   date,           /* (I) date to start forward */
   long                    dayCountConv,   /* (I) day counting convention */
   long                    interpType,     /* (I) interpolation type to use */
   TInterpData            *interpData,     /* (I) interpolation method to use */
   double                 *rate)           /* (O) money market rate */
{
   return JpmcdsZCForwardRate(zc,
                           zc->valueDate,
                           date,
                           JPMCDS_SIMPLE_BASIS,
                           dayCountConv,
                           interpType,
                           interpData,
                           rate);
}


/*
***************************************************************************
** Calculates a forward rate from a zero curve.
***************************************************************************
*/
int JpmcdsZCForwardRate(                         /* gets forward rate from ZCurve */
   ZCurve                 *zc,             /* (I) ZCurve */
   TDate                   date1,          /* (I) date to start forward */
   TDate                   date2,          /* (I) date to end forward */
   long                    basis,          /* (I) compounding basis */
   long                    dayCountConv,   /* (I) day counting convention */
   long                    interpType,     /* (I) interpolation type to use */
   TInterpData            *interpData,     /* (I) interpolation method to use*/
   double                 *rate)           /* (O) forward rate */
{
   double discount;                    /* discount factor over forward */

   if (JpmcdsZCFwdDiscountFactor(zc,
                              date1,
                              date2,
                              interpType,
                              interpData,
                             &discount)==FAILURE)
   {
       JpmcdsErrMsg("JpmcdsZCForwardRate: couldn't calculate discount factor\n");
       return FAILURE;
   }

   if (discount<=0.0)                  /* bad discount factor */
   {
       JpmcdsErrMsg("JpmcdsZCForwardRate: bad discount for forward\n");
       return FAILURE;
   }
   else
   {
       if (JpmcdsDiscountToRate( discount,
                              date1,
                              date2,
                              dayCountConv,
                              basis,
                              rate)==FAILURE)
       {
           JpmcdsErrMsg("JpmcdsZCForwardRate: couldn't convert discount to rate\n");
           return FAILURE;
       }
   }

   return SUCCESS;
}


/*
***************************************************************************
** Calculates a discount factor for a forward between two dates.
***************************************************************************
*/
int JpmcdsZCFwdDiscountFactor(                   /* discount factor for given forward */
   ZCurve                 *zc,             /* (I) ZCurve */
   TDate                   date1,          /* (I) date to start forward */
   TDate                   date2,          /* (I) date to end forward */
   long                    interpType,     /* (I) interpolation type to use */
   TInterpData            *interpData,     /* (I) interpolation method to use */
   double                 *discount)       /* (O) forward discount factor */
{
   double discount1;                   /* discount factor valueDate to date1 */
   double discount2;                   /* discount factor valueDate to date2 */

   if (JpmcdsZCDiscountFactor(zc,date1,interpType,interpData,&discount1)==FAILURE)
   {
       JpmcdsErrMsg("JpmcdsZCFwdDiscountFactor: couldn't calc discount1\n");
       return FAILURE;
   }

   if (JpmcdsZCDiscountFactor(zc,date2,interpType,interpData,&discount2)==FAILURE)
   {
       JpmcdsErrMsg("JpmcdsZCFwdDiscountFactor: couldn't calc discount2\n");
       return FAILURE;
   }

   *discount = 1.0;
   if (discount1<=0.0)                 /* bad discount factor */
   {
       JpmcdsErrMsg("JpmcdsZCFwdDiscountFactor: bad discount at start of fwd\n");
   }
   else
   {
       *discount = discount2 / discount1;
   }

   return SUCCESS;
}


/*
***************************************************************************
** Calculates net-present-value of a cash-flow list.
**
** Note: This embeds some logic in JpmcdsZCDiscountFactor() to optimize speed.
** Basically it avoids searching for matching dates in the zero curve.  This can
** usually works well, since a swap's coupon dates have been added by prior
** instruments.
***************************************************************************
*/
int JpmcdsJpmcdsZCPresentValueCFL(                     /* present value of a Cash-Flow-List */
   ZCurve                 *zc,             /* (I) ZCurve for discounting */
   TCashFlowList          *cfl,            /* (I) cash flow list to value */
   int                     iLo,            /* (I) index of 1st c.f. to value*/
   int                     iHi,            /* (I) index of last c.f. to value*/
   long                    interpType,     /* (I) type of interpolation */
   TInterpData            *interpData,     /* (I) data for interpolation */
   double                 *outPV)          /* (O) present value of cash flows*/
{
    static char routine[] = "JpmcdsJpmcdsZCPresentValueCFL";
    int         status = FAILURE;       /* Until proven successsful */

    double      sumPV = 0.0;            /* sum of p.v. of cash flows */
    int         i;                      /* loops over cash flows[iLo..iHi] */
    int         j = 0;                  /* loops over ZC entries */

    if(iLo < 0)
    {
        JpmcdsErrMsg("%s: iLo = %d is negative\n",routine,iLo);
        return FAILURE;
    }
    if(iLo > iHi)
    {
        JpmcdsErrMsg("%s: iLo(%d) > iHi(%d)\n",routine,iLo,iHi);
        return FAILURE;
    }
    if(cfl->fNumItems <= iHi)
    {
        JpmcdsErrMsg("%s: cfl->fNumItems(%d) <= iHi(%d)\n", routine,cfl->fNumItems,iHi);
        return FAILURE;
    }

    for (i=iLo; i<=iHi; i++)
    {
        double amt  = cfl->fArray[i].fAmount; /* short hand for c.f. amount */
        TDate  date = cfl->fArray[i].fDate; /* short hand for c.f. date */
        double pv;                      /* present value of c.f. */

        while (j<zc->numItems && zc->date[j]<date) /* push j up to c.f. date */
        {
            j++;
        }

        if (j<zc->numItems && zc->date[j]==date) /* found exact date in zc */
        {
            pv = zc->discount[j]*amt;
        }
        else
        {
            if (JpmcdsZCPresentValue(zc,amt,date,interpType,interpData, &pv) == FAILURE)
                goto done;
        }
        sumPV += pv;
    }
    *outPV = sumPV;
    status = SUCCESS;

  done:
    if (status == FAILURE)
        JpmcdsErrMsg("%s: Failed.\n", routine);

    return status;
}


/*
***************************************************************************
** Calculates net-present-value of a cash-flow (a payment at a given date 
** in the future).
***************************************************************************
*/
int JpmcdsZCPresentValue(                        /* present-value of future payment */
   ZCurve                     *zc,         /* (I) ZCurve used for discounting*/
   double                      price,      /* (I) amount of payment */
   TDate                       date,       /* (I) date of payment */
   long                        interpType, /* (I) type of interpolation */
   TInterpData                *interpData, /* (I) data for interpolation */
   double                     *pv)         /* (I) present value of price */
{
    static char routine[]="JpmcdsZCPresentValue";
    int         status = FAILURE;

    double  discFactor;

    if (JpmcdsZCDiscountFactor(zc,date,interpType,interpData,&discFactor) == FAILURE)
        goto done;

    status = SUCCESS;

    *pv = price * discFactor;

  done:
    if (status == FAILURE)
        JpmcdsErrMsg("%s: Failed.\n", routine);

    return status;
}


/*
***************************************************************************
** Calculates discount factor for a date.
***************************************************************************
*/
int JpmcdsZCDiscountFactor              /* discount factor for given date */
   (ZCurve       *zc,             /* (I) ZCurve */
    TDate         date,           /* (I) date to calculate d.f. for */
    long          interpType,     /* (I) interpolation type to use */
    TInterpData  *interpData,     /* (I) interpolation method to use */
    double       *discount)       /* (O) discount factor */
{
    static char routine[]="JpmcdsZCDiscountFactor";
    int         status = FAILURE;
    double      rate = 0;         /* Zero rate to date */

    if (date == zc->valueDate)    /* special case for valueDate */
    {
        *discount = 1.0;
        return SUCCESS;
    }
   
    if (JpmcdsZCInterpolate(zc,date,interpType,interpData,&rate) == FAILURE)
        goto done;

    /* Calculate discount factor from rate
     */
    if (JpmcdsZCComputeDiscount(zc, date, rate, discount) == FAILURE)
        goto done;

    status = SUCCESS;

  done:
    if (status == FAILURE)
        JpmcdsErrMsg("%s: Failed at %s, rate=%f.\n", routine, JpmcdsFormatDate(date), rate);

    return status;
}


/*
***************************************************************************
** Calculates an interpolated rate for a date.
**
** Note: Piece wise interpolation allows different areas of the zero curve
** to be interpolated using different methods.  Basically an array of 
** <date,interpolationStuff> is given, where the interpolationStuff is used
** for any date before the given date.  The code allows the interpolationStuff
** to be another pieceWise interpolation type, although the utility of this
** is unknown.
***************************************************************************
*/
int JpmcdsZCInterpolate(                         /* interpolation on a ZCurve */
   ZCurve               *zc,               /* (I) ZCurve to interpolate on */
   TDate                 date,             /* (I) date to interpolate to */
   long                  interpTypeIn,     /* (I) type of interpolation's use*/
   TInterpData          *interpDataIn,     /* (I) data for interp method */
   double               *rate)             /* (O) rate to return */
{
   int             status = FAILURE;       /* result of interpolation */
   TMetricDoubleFunc metricFP;             /* day counting metric function */
   int             lo;                     /* low index in ZCurve for linear */
   int             hi;                     /* high index ""     "" */
   long            interpType = interpTypeIn;  /* interpolation type */
   TInterpData    *interpData = interpDataIn;  /* interpolation data */
   char           routine[] = "JpmcdsZCInterpolate"; /* KM */
   double          yearsFromBaseDate = (date - zc->valueDate) / 365.25;

   if (zc->numItems < 1)
   {
       JpmcdsErrMsg("%s: No points in zero curve.\n", routine);
       goto done;
   }

   /* Do *flat* extrapolation only when going backwards. This
    * is done so that swaps which have payments before the beginning
    * of the stub zero curve will still value to par. This can happen
    * very easily if there are swaps with front stubs.
    * We still permit forward non-flat extrapolation.
    */
   if (date <= zc->date[0])
   {
       *rate = zc->rate[0];
       return SUCCESS;
   }

   if (date == zc->valueDate)
   {
       ++date;
       /* cannot calculate rate for value date, so get the value at the
       ** next date which is the next best thing
       */
   }

   /* Find indices in zero curve which bracket date.
    */
   if (JpmcdsBSearchLongFast((double)date,
                          zc->date, 
                          sizeof(TDate), 
                          zc->numItems,
                          &lo,
                          &hi) == FAILURE)
   {
       goto done;
   }

   if (zc->date[lo] == date)             /* if exact match, don't interpolate */
   {
       *rate = zc->rate[lo];
       return SUCCESS;
   }

   if (zc->date[hi] == date)             /* if exact match, don't interpolate */
   {
       *rate = zc->rate[hi];
       return SUCCESS;
   }

   metricFP = JpmcdsDayCountToMetricFunc(zc->dayCountConv);

   switch (interpType)
   {
   case JPMCDS_LINEAR_INTERP:
       {
           long hi_lo = zc->date[hi] - zc->date[lo];
           long dt_lo = date         - zc->date[lo];
           *rate = zc->rate[lo];
           if (hi_lo!=0)
           {
               *rate += ((zc->rate[hi] - zc->rate[lo]) / hi_lo) * dt_lo;
           }
       }
       break;

   case JPMCDS_FLAT_FORWARDS:
       {
           double  discLo;                  /* discount factor for lo */
           double  discHi;                  /* discount factor for hi */
           double  discDate;                /* discount factor for date */
           long    hi_lo = zc->date[hi] - zc->date[lo];
           long    dt_lo = date         - zc->date[lo];
           
           
           if (hi_lo==0)
           {
               *rate = zc->rate[lo];
           }
           else
           {
               if (JpmcdsZCComputeDiscount
                   (zc, zc->date[lo], zc->rate[lo], &discLo) == FAILURE)
                   goto done;
               
               
               if (JpmcdsZCComputeDiscount
                   (zc, zc->date[hi], zc->rate[hi], &discHi) == FAILURE)
                   goto done;
               
               if (discLo==0.0)
               {
                   JpmcdsErrMsg("%s: zero discount factor.\n", routine);
                   goto done;
               }
               
               
               discDate = discLo * pow(discHi/discLo, dt_lo / (double) hi_lo);
               
               if (JpmcdsDiscountToRate
                   (discDate,
                    zc->valueDate, date,
                    zc->dayCountConv, zc->basis,
                    rate) == FAILURE)
                   goto done;
           } /* else */
       }
       break;

   default:
       JpmcdsErrMsg("%s: Bad interpolation type (%ld).\n", routine, interpType);
       status = FAILURE;
       goto done;
       break;
   }

   /* Made it through. */
   status = SUCCESS;

 done:
   if (status == FAILURE)
       JpmcdsErrMsg("%s: Failed.\n", routine);

   return status;
}


/*
***************************************************************************
** Converts zc-style rate into a discount factor.
**
** Note: As a special optimization, the normal ZCurve case is hard coded for
** speed (since this routine is critical to performance of building zero-curves).
** That is, annually-compounded rates w/ 360 or 365-fixed year are calculated
** directly.
***************************************************************************
*/
int JpmcdsZCComputeDiscount(                     /* calculate discount from ZC rate */
   ZCurve               *zc,               /* (I) ZCurve assoc'd with rate */
   TDate                 date,             /* (I) date assoc'd with rate */
   double                rate,             /* (I) rate in ZC basis & dayCount*/
   double               *discount)         /* (O) discount factor */
{
   if (zc->basis == 1              && 
       rate      >= -1.0           &&
       date      >= zc->valueDate  &&
       (zc->dayCountConv == JPMCDS_ACT_365F || zc->dayCountConv == JPMCDS_ACT_360)) 
   {
       *discount = pow(1+rate, 
                       (zc->valueDate - date) / 
                       (zc->dayCountConv == JPMCDS_ACT_360 ? 360.0 : 365.0));
       return SUCCESS;
   }

   if (JpmcdsRateToDiscount(rate,
                         zc->valueDate,
                         date,
                         zc->dayCountConv,
                         zc->basis,
                         discount) == FAILURE)
   {
      JpmcdsErrMsg("JpmcdsZCComputeDiscount: couldn't calculate discount\n");
      return FAILURE;
   }

   return SUCCESS;
}


/*
***************************************************************************
** Finds an exact date in a ZC.
** Returns index of matching date if found, -1 if not.
***************************************************************************
*/
int JpmcdsZCFindDateExact(                       /* finds an exact date in a ZCurve */
   ZCurve         *zc,                     /* (I) ZCurve to search */
   TDate           date)                   /* (I) date to search for */
{
   int i = JpmcdsZCFindDateClosest(zc,date);  /* search for closest date */

   if (i>=0 && zc->date[i]!=date)
   {
       i = -1;
   }

   return i;
}


/*
***************************************************************************
** Finds index of entry closest to an input date.
** Returns index closest to date.
***************************************************************************
*/
int JpmcdsZCFindDateClosest(                     /* finds a closest date in a ZCurve */
   ZCurve         *zc,                     /* (I) ZCurve to search */
   TDate           date)                   /* (I) date to search for */
{
   int   lo;                            /* index near date, but below */
   int   hi;                            /* index near date, but above */
   int   diffLo;                        /* # days away from lo date */
   int   diffHi;                        /* # days away from hi date */
   
   if (zc == NULL || zc->numItems <= 0)  /* if no entries, return -1 */
       return -1;

   if (JpmcdsBSearchLongFast((double)date,
                          zc->date,
                          sizeof(TDate),
                          zc->numItems,
                          &lo,
                          &hi) == FAILURE)
   {
       JpmcdsErrMsg("JpmcdsZCFindDateClosest: bSearch failed\n");
       return -1;
   }

   diffLo = (int)ABS(date - zc->date[lo]);
   diffHi = (int)ABS(date - zc->date[hi]);
    
   return (diffLo <= diffHi) ? lo : hi;
}

