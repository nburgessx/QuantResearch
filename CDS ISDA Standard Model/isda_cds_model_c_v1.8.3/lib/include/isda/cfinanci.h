/*
 * ISDA CDS Standard Model
 *
 * Copyright (C) 2009 International Swaps and Derivatives Association, Inc.
 * Developed and supported in collaboration with Markit
 * 
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the ISDA CDS Standard Model Public License.
 */

#ifndef CFINANCI_H
#define CFINANCI_H

#include "cgeneral.h"
#include "cdate.h"

#ifdef __cplusplus
extern "C"
{
#endif


/*------- Calculating discount factors----------*/
#define JPMCDS_CONTINUOUS_BASIS 5000     /* Means continuous compounding */
#define JPMCDS_DISCOUNT_RATE    512      /* Discount RATE (not factor) K's b-Day*/
#define JPMCDS_SIMPLE_BASIS 0            /* Means simple interest */
#define JPMCDS_ANNUAL_BASIS 1            /* Annual compounding */
#define JPMCDS_DISCOUNT_FACTOR  -2       /* Rate is stored as a discount factor */


/*f
***************************************************************************
** Convert from a given rate to a discount factor.
***************************************************************************
*/
int JpmcdsRateToDiscount
    (double         rate,               /* (I) Rate */
     TDate          startDate,          /* (I) Start date */
     TDate          endDate,            /* (I) End date */
     long           rateDayCountConv,   /* (I) See JpmcdsDayCountFraction */
     long           rateType,           /* (I) Type of rate */
     double        *discount);          /* (O) Output discount factor */


/*f
***************************************************************************
** Convert from a given discount factor to a rate.
***************************************************************************
*/
int JpmcdsDiscountToRate
    (double         discount,           /* (I) Discount factor */
     TDate          startDate,          /* (I) Start date */
     TDate          endDate,            /* (I) End date */
     long           rateDayCountConv,   /* (I) See JpmcdsDayCountFraction */
     long           rateType,           /* (I) Type of rate */
     double        *rate);              /* (O) Output rate */


/*f
***************************************************************************
** Convert from a discount factor to a rate using a year fraction.
**
** Calculate an interest rate from a discount factor, given the daycount 
** and basis.
***************************************************************************
*/
int JpmcdsDiscountToRateYearFrac
    (double discount,                   /* (I) Discount factor */
     double yearFraction,               /* (I) See JpmcdsDayCountFraction */
     double basis,                      /* (I) Basis for the rate */
     double *rate);                     /* (O) Output rate */


/*f
***************************************************************************
** Convert from a rate to a discount factor using a year fraction.
**
** Calculate discount factor from an interest rate, given the daycount 
** and basis.
***************************************************************************
*/
int JpmcdsRateToDiscountYearFrac
    (double rate,                       /* (I) Rate */
     double yearFraction,               /* (I) See JpmcdsDayCountFraction */
     double basis,                      /* (I) Basis for the rate */
     double *discount );                /* (O) Output discount rate */


/*f
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
    double rateBasis);        /* (I) Compounding basis */


/*f
***************************************************************************
** Checks whether a rate is valid between two dates.
***************************************************************************
*/
int JpmcdsRateValidYearFrac(
    char  *routine,           /* (I) Routine name to print */
    double rate,              /* (I) Rate to validate */
    double yearFraction,      /* (I) Fraction of year */
    double basis);            /* (I) Compounding basis */


#ifdef __cplusplus
}
#endif

#endif    /* CFINANCI_H */

