/*
 * ISDA CDS Standard Model
 *
 * Copyright (C) 2009 International Swaps and Derivatives Association, Inc.
 * Developed and supported in collaboration with Markit
 * 
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the ISDA CDS Standard Model Public License.
 */

#include "cfinanci.h"
#include <math.h>

#ifdef __cplusplus
extern "C"
{
#endif

#ifdef __cplusplus
}
#endif

#include "cgeneral.h"
#include "bastypes.h"
#include "cerror.h"
#include "convert.h"
#include "ldate.h"
#include "macros.h"


/*
***************************************************************************
** Convert from a given discount factor to a rate.
***************************************************************************
*/
int JpmcdsDiscountToRate
    (double         discount,         /* (I) Discount factor */
     TDate          startDate,        /* (I) Start date */
     TDate          endDate,          /* (I) End date */
     long           rateDayCountConv, /* (I) Day count convention for rate */
     long           rateBasis,        /* (I) Basis for the rate */
     double          *rate)           /* (O) output rate */
{
    static char routine[] = "JpmcdsDiscountToRate" ;
    int status = FAILURE;               /* Until proven successful */

    double      rateYF ;

    if (discount <= 0.0)
    {
        JpmcdsErrMsg ("%s: Bad discount factor %f.\n",routine,discount) ;
        goto done;
    }

    if (rateBasis == JPMCDS_DISCOUNT_FACTOR)
    {
        if (startDate == endDate)
        {
            /*
             * In this case, the discount factor should be 1.0.
             */
            if (! ARE_ALMOST_EQUAL(discount, 1.0))
            {
                JpmcdsErrMsg("%s: Start date (%s) equals end date (%s), "
                          "but input disount factor (%f) <> 1.0\n",
                          routine,
                          JpmcdsFormatDate(startDate),
                          JpmcdsFormatDate(endDate),
                          discount);
                goto done;
            }
            *rate = 1.0;
            return SUCCESS;
        }
        else
        {
            *rate = discount;
            return SUCCESS;
        }
    }

    if (startDate == endDate)
    {
        JpmcdsErrMsg("%s: Start date %s cannot equal end date.\n",
                  routine, 
                  JpmcdsFormatDate(startDate));
        goto done;
    }

    if (rateBasis < JPMCDS_SIMPLE_BASIS)
    {
        JpmcdsErrMsg ("%s: Basis (%ld) < %d.\n",
                   routine, rateBasis, JPMCDS_SIMPLE_BASIS) ;
        goto done;
    }
    
    /* get year fractions for rate */
    if (JpmcdsDayCountFraction(startDate,
                            endDate,
                            rateDayCountConv,
                            &rateYF) == FAILURE)
    {
        goto done;
    }

    if (JpmcdsDiscountToRateYearFrac(discount,
                                  rateYF,
                                  (double) rateBasis,
                                  rate) != SUCCESS)
    {
        goto done;
    }

    status = SUCCESS;

done:
    if (status == FAILURE)
    {
        JpmcdsErrMsg("%s: Failed.\n", routine);
    }

    return status;
}


/*
***************************************************************************
** Convert from a given rate to a discount factor.
***************************************************************************
*/
int JpmcdsRateToDiscount
    (double         rate,             /* (I) Rate */
     TDate          startDate,        /* (I) Start date */
     TDate          endDate,          /* (I) End date */
     long           rateDayCountConv, /* (I) Day count convention for rate */
     long           rateBasis,        /* (I) Basis for the rate */
     double        *discount)         /* (O) Discount factor */
{
    int status = FAILURE;               /* Until proven successful */
    static char routine[] = "JpmcdsRateToDiscount" ;
    double      rateYF ;

    if (rateBasis == JPMCDS_DISCOUNT_FACTOR)
    {
       if (rate <= 0.0)
       {
           JpmcdsErrMsg("%s: Bad rate (discount factor) %f\n",routine,rate);
           goto done;
       }
       *discount = rate;
       return SUCCESS;
    }
       
    if (rateBasis < JPMCDS_SIMPLE_BASIS)
    {
        JpmcdsErrMsg ("%s: Basis (%ld) < %d.\n",routine,rateBasis,
                   JPMCDS_SIMPLE_BASIS);
        goto done;                      /* Failed */
    }

    /* get year fractions for rate */
    if (JpmcdsDayCountFraction(startDate, 
                            endDate, 
                            rateDayCountConv, 
                            &rateYF) == FAILURE)
    {
        goto done;                      /* Failed */
    }

    if (JpmcdsRateToDiscountYearFrac(rate, 
                                  rateYF, 
                                  (double)rateBasis,
                                  discount) != SUCCESS)
    {
        goto done;
    }

    status = SUCCESS;
    
done:
    if (status == FAILURE)
    {
        JpmcdsErrMsg("%s: Failed.\n", routine);
    }    
    
    return status;
}


/*
***************************************************************************
** Convert from a discount factor to a rate using a year fraction.
**
** Calculate an interest rate from a discount factor, given the daycount 
** and basis.
***************************************************************************
*/
int JpmcdsDiscountToRateYearFrac(
    double  discount,                  /* (I) Discount factor */
    double  yearFraction,              /* (I) Year fraction */ 
    double  basis,                     /* (I) Basis for the rate */
    double *rate )                     /* (O) Output rate */
{
    static char routine[]="JpmcdsDiscountToRateYearFrac";
    
    if (discount <= 0.0)
    {
        JpmcdsErrMsg("%s: Discount %f <= 0.\n", routine, discount);
        *rate = 0.0;
        return FAILURE;
    }

    switch (long (basis))
    {

    case JPMCDS_SIMPLE_BASIS:
        if (IS_ALMOST_ZERO(yearFraction))
        {
            JpmcdsErrMsg("%s: Bad discount:%f or year fraction:%f\n",
                      routine, discount, yearFraction) ;
            *rate = 0;
            return FAILURE;
        }
        *rate = (1.0/discount - 1.0)/yearFraction ;
        break;

    case JPMCDS_DISCOUNT_RATE:
        if (IS_ALMOST_ZERO(yearFraction))
        {
            JpmcdsErrMsg("%s: Undefined discount rate\n",routine);
            *rate = 0.0;
            return FAILURE;
        }
        *rate = (1.0 - discount) / yearFraction;
        break;

    case JPMCDS_CONTINUOUS_BASIS:
        if (IS_ALMOST_ZERO(yearFraction))
        {
            JpmcdsErrMsg("%s: Div by zero - yr frac %f.\n",
                      routine, yearFraction) ;
            *rate = 0.0;
            return FAILURE;
        }
        *rate = (-log(discount)/yearFraction) ;
        break;

    case JPMCDS_DISCOUNT_FACTOR:
        *rate = discount;
        break;

    default:
        if (IS_ALMOST_ZERO(yearFraction))
        {
            JpmcdsErrMsg ("%s: Div by zero - yr frac %f.\n",
                       routine, yearFraction) ;
            *rate = 0.0;
            return FAILURE;
        }
        /* We prefer to be stodgy and not do the log/exp thing
         * here for performance.
         */
        *rate = basis * (pow(discount, -1.0/(basis*yearFraction)) - 1.0);
        break;

    }

    return SUCCESS;
}


/*
***************************************************************************
** Convert from a rate to a discount factor using a year fraction.
**
** Calculate discount factor from an interest rate, given the daycount 
** and basis.
***************************************************************************
*/
int JpmcdsRateToDiscountYearFrac(
    double  rate,           /* (I) The rate */
    double  yearFraction,   /* (I) Year fraction */ 
    double  basis,          /* (I) Basis for the rate */
    double *discount )      /* (O) Output discount rate */
{
    static char routine[]="JpmcdsRateToDiscountYearFrac";

    switch ((long)basis)
    {
    case JPMCDS_SIMPLE_BASIS:
        {
            /* Simple basis is 0 so don't worry bout dividing by it
             */
            double denom = 1.0 + rate * yearFraction;
            if (denom <= 0.0 ||
                IS_ALMOST_ZERO(denom))
            {
                JpmcdsErrMsg("%s: Invalid simple interest rate:%f\n",
                          routine, rate);
                *discount = 0.0;
                goto done;
            }
            *discount = 1.0 / denom;
        }
        break;

    case JPMCDS_DISCOUNT_RATE:
        if (IS_ALMOST_ZERO(yearFraction))
        {
            *discount = 1.0;
        }
        else
        {
            *discount = 1.0 - rate * yearFraction;
            if (*discount <= 0.0)
            {
                JpmcdsErrMsg("%s: Invalid discount rate:%f\n",
                          routine, rate);
                *discount = 0.0;
                goto done;
            }
        }
        break;

    case JPMCDS_CONTINUOUS_BASIS:
        *discount = exp(-rate*yearFraction);
        break;

    case JPMCDS_DISCOUNT_FACTOR:
        *discount = rate;
        break;
       
    default:
        {
            double tmp = 1.0 + rate / basis;
            /* Since pow(x,y) is not defined when x < 0 and y is not an integer,
             * check before calling it.
             */
            if (tmp <= 0.0 ||
                IS_ALMOST_ZERO(tmp))
            {
                JpmcdsErrMsg("%s: Bad rate: %f.\n", routine, rate);
                *discount = 0.0;
                goto done;
            }
            else
            {
                /* We prefer to be stodgy and not do the log/exp thing
                 * here for performance.
                 */
                *discount = pow( tmp, -basis*yearFraction);
            }
        }
        break;
    }
    
    return SUCCESS;

 done:
    JpmcdsErrMsg("%s: Failed.\n", routine);
    return FAILURE;
}


/*
***************************************************************************
** Checks whether a rate is valid between two dates.
***************************************************************************
*/
int JpmcdsRateValid(
    char  *routine,           /* (I) Routine name to print */
    double rate,              /* (I) Rate to validate */
    TDate  startDate,         /* (I) Starting date */
    TDate  endDate,           /* (I) Ending date */
    long   rateDayCountConv,  /* (I) Day count convention */
    double rateBasis)         /* (I) Compounding basis */
{
    int status = FAILURE;

    double yearFraction;

    switch ((long) rateBasis) 
    {
    case JPMCDS_SIMPLE_BASIS:
    case JPMCDS_DISCOUNT_RATE:
        /* Only need year fraction in these cases.
         */
        if (JpmcdsDayCountFraction(startDate, 
                                endDate, 
                                rateDayCountConv,
                                &yearFraction) != SUCCESS)
        {
            goto done;
        }
        break;

    default:
        yearFraction = 1.0;
    }

    if (JpmcdsRateValidYearFrac(routine,
                             rate,
                             yearFraction,
                             rateBasis) != SUCCESS)
    {
        goto done;
    }


    status = SUCCESS;

 done:
    return status;
}


/*
***************************************************************************
** Checks whether a rate is valid between two dates.
***************************************************************************
*/
int JpmcdsRateValidYearFrac(
    char  *routine,           /* (I) Routine name to print */
    double rate,              /* (I) Rate to validate */
    double yearFraction,      /* (I) Fraction of year */
    double basis)             /* (I) Compounding basis */
{
    int status = FAILURE;

    switch ((long) basis)
    {
    case JPMCDS_SIMPLE_BASIS:
        if (rate * yearFraction <= -1.0)
        {
            JpmcdsErrMsg("%s: Simple Rate (%f) * Year Fraction (%f) must "
                      "be > -1.0.\n",
                      routine, rate, yearFraction);
            goto done;
        }
        break;
        
    case JPMCDS_DISCOUNT_RATE:
        if (rate * yearFraction >= 1.0)
        {
            JpmcdsErrMsg("%s: Discount Rate (%f) * Year Fraction (%f) must "
                      "be < 1.0.\n",
                      routine, rate, yearFraction);
            goto done;
        }
        break;

    case JPMCDS_CONTINUOUS_BASIS:
        /* Any rate is valid */
        status = SUCCESS;
        goto done;

    case JPMCDS_DISCOUNT_FACTOR:
        if (rate <= 0.0)
        {
            JpmcdsErrMsg("%s: Discount factor (%f) must be > 0.0.\n",
                      routine, rate);
            goto done;
        }
        break;

    default:
        if (rate <= -basis)
        {
            JpmcdsErrMsg("%s: Rate (%f) must be greater than -basis (%f).\n",
                      routine, rate, -basis);
            goto done;
        }
        break;
    }

    status = SUCCESS;

 done:
    return status;
}
