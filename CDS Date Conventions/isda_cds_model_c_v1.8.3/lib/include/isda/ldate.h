/*
 * ISDA CDS Standard Model
 *
 * Copyright (C) 2009 International Swaps and Derivatives Association, Inc.
 * Developed and supported in collaboration with Markit
 * 
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the ISDA CDS Standard Model Public License.
 */

#ifndef LDATE_H
#define LDATE_H

#include "cdate.h"
#include "dateconv.h"
#include "busday.h"
#include "dateadj.h"

#ifdef __cplusplus
extern "C"
{
#endif


/* DAY COUNT FRACTION METHODS. Note that JPMCDS_ACT_ACT is used in
 * two ways in this system. It is used as a day count convention
 * in which the denominator is a weighted average of 365 and 366,
 * depending on the number of days in the years covered in the period.
 * It is also used by the bond library to calculate the accrued-
 * which depends on the fraction of the coupon period elapsed.
 * Note that JPMCDS_B30_360_FIXED is not a day count convention like
 * the others. It can be used to generate bond coupon payments.
 * These bond coupon payments have the property that they are 
 * always the same. JPMCDS_B30_360 does not have this property in some
 * special cases. 
 */
#define JPMCDS_BAD_DCC       -2L                /* Dummmy Day Count Conv for 
                                             ** Brazilian bods*/
#define JPMCDS_ACT_365        1L               /* Actual/365 */
#define JPMCDS_ACT_365F       2L               /* Actual/365 Fixed */
#define JPMCDS_ACT_360        3L               /* Actual/360 */
#define JPMCDS_B30_360        4L               /* 30/360 */
#define JPMCDS_B30E_360       5L               /* 30E/360 */
#define JPMCDS_ACT_ACT        JPMCDS_ACT_365
#define JPMCDS_EFFECTIVE_RATE 8L             /* Effective rate - YF always 1 */

#define JPMCDS_DCC_LAST       JPMCDS_EFFECTIVE_RATE


/*f
***************************************************************************
** Moves forward by the specified interval.
***************************************************************************
*/
int JpmcdsDtFwdAny(
    TDate          startDate,                  /* (I) date            */
    TDateInterval *interval,                   /* (I) dateInterval    */
    TDate         *sumDate);                   /* (O) date + dateInterval */


/*f
***************************************************************************
** Calculates Day Count Fractions for following methods:
** 30/360, 30E/360, ACT/365, ACT/365F, ACT/360, etc.
***************************************************************************
*/
int JpmcdsDayCountFraction (
    TDate    beginDate,  /* (I)                                              */
    TDate    endDate,    /* (I)                                              */
    long     method,     /* (I) Specifies day count method:                  */
                         /*    JPMCDS_B30_360   - 30/360 USA NASD               */
                         /*    JPMCDS_B30E_360  - 30E/360 AIBD RULE 220         */
                         /*    JPMCDS_ACT_365   - ACT/365 or ACT/ACT            */
                         /*    JPMCDS_ACT_365F  - ACT/365 Fixed                 */
                         /*    JPMCDS_ACT_360   - ACT/360                       */
    double  *yearFrac);  /* (O) Fraction of year between two dates           */


/*f
***************************************************************************
** Calculates difference between two dates.
**
** Methods currently covered:
** - Bond basis: 30, 30E
** - Actual
**
** This function calculates the numerator for day count fraction
***************************************************************************
*/
int JpmcdsDaysDiff (
    TDate    beginDate,  /* (I)                                              */
    TDate    endDate,    /* (I)                                              */
    long     method,     /* (I) Specifies day count method:                  */
                         /*    JPMCDS_B30_360   - 30/360 USA NASD               */
                         /*    JPMCDS_B30E_360  - 30E/360 AIBD RULE 220         */
                         /*    JPMCDS_ACT_365   - ACT/365 or ACT/ACT            */
                         /*    JPMCDS_ACT_365F  - ACT/365 Fixed                 */
                         /*    JPMCDS_ACT_360   - ACT/360                       */
    long *numDays);      /* (O)  Number days between two dates               */


/*f
***************************************************************************
** Test if date falls in a leap year.
***************************************************************************
*/
int JpmcdsIsLeap(
    TDate  date,                  /* (I)                             */
    long  *isLeap);               /* (O) 0 - if not a leap year      */
                                  /*     1 - if is a leap year       */

/*f
***************************************************************************
** Given a date in TDate, retrieve year as a four digit year.
***************************************************************************
*/
int JpmcdsYear(TDate , int *);


/*f
***************************************************************************
** Calculates the number of days left in a year.
***************************************************************************
*/
int JpmcdsDaysLeftThisYear(TDate , int, long *);


/*f
***************************************************************************
** Given a date in TDate, retrieve beginning of the year in TDate format
***************************************************************************
*/
int JpmcdsYearStart(TDate , TDate *);


/*f
***************************************************************************
** Moves forward by the specified TDateInterval, and then adjust to find a
** good day. Basically a call to JpmcdsDtFwdAny followed by a call to
** JpmcdsBusinessDay.
***************************************************************************
*/
int  JpmcdsDateFwdThenAdjust
(TDate           date,               /* (I) Start date */
 TDateInterval  *interval,           /* (I) Interval to advance by */
 long            badDayMethod,       /* (I) JPMCDS_BAD_DAY_XYZ */
 char           *holidayFile,        /* (I) Holiday file to use */
 TDate          *advAdjustedDate);   /* (O) Advanced adjusted date */


#ifdef __cplusplus
}
#endif

#endif    /* LDATE_H */

