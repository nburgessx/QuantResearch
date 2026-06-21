/*
 * ISDA CDS Standard Model
 *
 * Copyright (C) 2009 International Swaps and Derivatives Association, Inc.
 * Developed and supported in collaboration with Markit
 * 
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the ISDA CDS Standard Model Public License.
 */

#include "ldate.h"
#include <math.h>

#ifdef __cplusplus
extern "C"
{
#endif

#include <string.h>
#include <ctype.h>

#ifdef __cplusplus
}
#endif

#include "cgeneral.h"
#include "bastypes.h"
#include "cdate.h"
#include "cmemory.h"
#include "cerror.h"
#include "strutil.h"
#include "macros.h"
#include "date_sup.h"
#include "convert.h"


static int dateToMDYFast
    (TDate         startDate,           /* (I) */
     TMonthDayYear *mdy);               /* (O) */


/*
***************************************************************************
** Test if date falls in a leap year.
***************************************************************************
*/
int JpmcdsIsLeap
   (TDate DATE,                 /* (I) DATE */
    long *IS_LEAP               /* (O) 0 - not a leap year
                                      1 - is a leap year */
    )
{
    TMonthDayYear temp;
    int year;


    /*  Returns zero if argument date is not a leap year     */
    /*  otherwise returns 1.                                 */


    *IS_LEAP=0;     /* assume not a leap year */


    if (JpmcdsDateToMDY(DATE,&temp) == FAILURE)
        return FAILURE;

    year=(int)temp.year;
    /*LINTED*/
    if (year%4 != 0)
        ; /* not divisible by 4 */
    else
        if (year%100 != 0)
            *IS_LEAP=1; /* divisible by 4 but not 100 */
        else
            /*LINTED*/
            if (year%400 != 0)
                ; /* divisible by 100 but not 400 */
            else
                *IS_LEAP=1;       /* divisible by 400, so is a leap year */

    return SUCCESS;
}


/*
***************************************************************************
** Moves forward by the specified interval.
***************************************************************************
*/
int JpmcdsDtFwdAny
    (TDate           startDate,      /* (I) date */
     TDateInterval   *interval,      /* (I) dateInterval */
     TDate           *sumDate        /* (O) date + dateInterval */
     )
{
    static char          routine[]="JpmcdsDtFwdAny";
    int                  status = FAILURE; /* Until proven successful */
    TMonthDayYear        mdy;
    TDateInterval        intval;
    char                 upperPrdTyp;


    if (startDate <= 0)
    {
        JpmcdsErrMsg("%s: startDate must be > 0.\n",  routine);
        goto done;
    }


    upperPrdTyp = (char) toupper(interval->prd_typ);
    switch (upperPrdTyp)
    {
      case 'M':                     /* MONTHly increments */
      case 'A':                     /* ANNUAL increments */
      case 'Y':                     /* YEARly increments */
      case 'S':                     /* SEMIANNUALL increments */
      case 'Q':                     /* QUARTERly increments */

        /* First reduce all of these types to monthly. Note that we shouldn't
         * really need code to handle A,S & Q, since JpmcdsMakeDateInterval
         * converts all of them to 'M' anyway. This code is left here for
         * people who have set up their own TDateIntervals without calling
         * JpmcdsMakeDateInterval or one of the supplied macros. The month
         * type is checked first for efficiency reasons.
         */
        intval.flag = 0;
        intval.prd_typ = 'M';           /* months */

        if (upperPrdTyp == 'M')
            intval.prd = interval->prd;
        else if (upperPrdTyp == 'A' || interval->prd_typ == 'Y')
            intval.prd = interval->prd * JPMCDS_MONTHS_PER_YEAR;
        else if (upperPrdTyp == 'S')
            intval.prd = interval->prd * JPMCDS_MONTHS_PER_SEMI;
        else
            intval.prd = interval->prd * JPMCDS_MONTHS_PER_QUARTER;


        if (dateToMDYFast(startDate, &mdy) == FAILURE)
            goto done;


        mdy.month += intval.prd;
        if (JpmcdsNormalizeMDY(&mdy) == FAILURE)
            goto done;

        if (JpmcdsMDYToDate(&mdy,sumDate) == FAILURE)
            goto done;
        break;


      case 'D':                         /* DAYly increments */
        *sumDate=startDate+interval->prd;
        break;


      case 'W':
        *sumDate=startDate+interval->prd*JPMCDS_DAYS_PER_WEEK;
        break;                          /* WEEKly increments */

      default:
        if (isalpha((int)interval->prd_typ))
        {
             JpmcdsErrMsg("%s:  Period type %c unknown.\n",
                       routine, interval->prd_typ);
        }
        else
        {
             JpmcdsErrMsg("%s:  Period type %d unknown.\n",
                       routine, (int)interval->prd_typ);
        }

        goto done;
    }

    status = SUCCESS;


  done:
    if (status == FAILURE)
        JpmcdsErrMsg("%s: Failed.\n", routine);

    return status;
}


static int dateToMDYFast
    (TDate         startDate,           /* (I) */
     TMonthDayYear *mdy)                /* (O) */
{
    static TDate         lastStartDate; /* For efficiency */
    static TMonthDayYear lastMDY;       /* For efficiency */

    /* Check if we've already done this before to avoid
     * calling JpmcdsDateToMDY if not necessary.
     */
    if (startDate == lastStartDate)
    {
        *mdy = lastMDY;                 /* Same as last time */
    }
    else                                /* Need to calculate */
    {
        if (JpmcdsDateToMDY(startDate,mdy) == FAILURE)
            return FAILURE;
        else
        {
            lastStartDate = startDate;  /* Save for next time */
            lastMDY = *mdy;
        }
    }
    return SUCCESS;
}


/*
***************************************************************************
** Calculates Day Count Fractions for following methods:
** 30/360, 30E/360, ACT/365, ACT/365F, ACT/360, etc.
***************************************************************************
*/
int JpmcdsDayCountFraction
     (TDate   date1,            /* (I) Start date */
      TDate   date2,            /* (I) End date   */
      long    method,           /* (I) method: B30/360
                                                   B30E/360
                                                   ACT/365 (ACT/ACT)
                                                   ACT/365 FIXED
                                                   ACT/360
                                                   NL/365 FIXED (no leap)
                                                   B30E/360 (Italy)
                                                   30/360SIA
                                        */
       double *result           /* (O) Day count fraction*/
                        )
{
    static char routine[] = "JpmcdsDayCountFraction";

    int status = FAILURE;       /* Until successful */
    TDate cur_date;
    TDate temp;
    int start_year,end_year;
    double sign = 1.0;
    int i;
    long act_days,leap_days,non_leap_days,days_left;
    long isleap;

    /*
     * This duplicates some code below and JpmcdsDaysDiff() and is
     * done for performance reasons.
     */
    if (method == JPMCDS_ACT_365F)
    {
        *result = (date2-date1)/365.;
        return SUCCESS;
    }
    else if (method == JPMCDS_ACT_360)
    {
        *result = (date2-date1)/360.;
        return SUCCESS;
    }

    leap_days = 0;
    non_leap_days = 0;

    /* Check if same date
     */
    if (date1 == date2)
    {
        *result = 0;
        return SUCCESS;
    }

    /* Check the date order
     */
    if (date1 > date2)
    {
        sign = -1.0;

        /* reverse order */
        temp = date1;
        date1 = date2;
        date2 = temp;
    }

    if (method == JPMCDS_EFFECTIVE_RATE)
    {
        /*
         * Effective rates have a year fraction of 1.0 or -1.0, depending
         * on the order of the dates (note: if the dates are the same,
         * the year fraction is 0 (handled above).
         */
        *result = sign;
        return SUCCESS;
    }

    /* Calculate the number of days between two dates
     */
    if (JpmcdsDaysDiff(date1,date2,method,&act_days) == FAILURE)
        goto done;

    switch (method)
    {
    case JPMCDS_B30_360:
    case JPMCDS_B30E_360:
    case JPMCDS_ACT_360:
        *result = act_days/360.0;
        break;

    case JPMCDS_ACT_365F:
        *result = act_days/365.0;
        break;

    case JPMCDS_ACT_365:

        /* weighted average of leap days divided by 366 +
           weighted average of non-leap days divided by 365 */

        /* handle first year */
        if (JpmcdsDaysLeftThisYear(date1,(int)method,&days_left) == FAILURE)
            goto done;

        if (JpmcdsIsLeap(date1,&isleap) == FAILURE)
            goto done;


        if (isleap == 1)
            leap_days += MIN(act_days,days_left);
        else
            non_leap_days += MIN(act_days,days_left);
        /**********************/


        /* loop through the years */
        if(JpmcdsYear(date1, &start_year) == FAILURE)
            goto done;

        if(JpmcdsYear(date2, &end_year) == FAILURE)
            goto done;

        cur_date = date1;

        /* loop through full years */
        for (i=start_year+1; i<end_year; i++)
        {
            /* check if previous year is a leap */
            if (isleap == 1)
                cur_date += 366; /* next year in days*/
            else
                cur_date += 365; /* next year in days */

            /* check if new year is a leap */
            if (JpmcdsIsLeap(cur_date, &isleap) == FAILURE)
                goto done;

            if (isleap == 1)
                leap_days += 366;
            else
                non_leap_days += 365;
        }

        /* handle last year */
        if (start_year != end_year)
        {
            if (JpmcdsIsLeap(date2, &isleap) == FAILURE)
                goto done;

            if (JpmcdsYearStart(date2, &cur_date) == FAILURE)
                goto done;

            /* number of days elapsed in the last year */
            if (JpmcdsDaysDiff(cur_date, date2, method, &days_left) == FAILURE)
                goto done;

            if (isleap == 1)
                leap_days += days_left;
            else
                non_leap_days += days_left;
        }

        /* calculate final day count fraction */

        /* ISDA interpretation */
        *result = leap_days/366.0 + non_leap_days/365.0;
        break;

    default:
        JpmcdsErrMsg("Invalid method (%ld).\n", method);
        goto done;

    }
    status = SUCCESS;

    *result *= sign;


  done:
    if (status == FAILURE)
        JpmcdsErrMsg("%s: Failed.\n", routine);

    return status;
}


/*
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
int JpmcdsDaysDiff
    (TDate date1,          /* (I) start date */
     TDate date2,          /* (I) end date   */
     long method,          /* (I) method: B30/360
                              B30E/360
                              ACT/365
                              ACT/365 FIXED
                              ACT/360
                           */
     long *result          /* (O) Numerator for day count frac */
     )
{
    static char routine[]="JpmcdsDaysDiff";
    int           status = FAILURE;
    TMonthDayYear date1MDY, date2MDY;
    TDate         temp;
    long          Y1,Y2,M1,M2,D1,D2;
    int           negative;

    negative = 0;

    if (date1 > date2)
    {
        negative = 1;

    /* reverse order */
        temp = date1;
        date1 = date2;
        date2 = temp;
    }

    switch (method)
    {
      case JPMCDS_B30_360:

        if (JpmcdsDateToMDY(date1, &date1MDY) == FAILURE)
            goto done;

        if (JpmcdsDateToMDY(date2, &date2MDY) == FAILURE)
            goto done;

        Y1 = date1MDY.year;
        M1 = date1MDY.month;
        D1 = date1MDY.day;

        Y2 = date2MDY.year;
        M2 = date2MDY.month;
        D2 = date2MDY.day;

        /* D1=31 => change D1 to 30 */
        if (D1 == 31)
            D1 = 30;

        /* D2=31 and D1 is 30 or 31 => change D2 to 30 */
        if (D2 == 31 && D1 == 30)
            D2 = 30;

        *result = (Y2-Y1)*360 + (M2-M1)*30 + (D2-D1);
        break;

      case JPMCDS_B30E_360:

        if (JpmcdsDateToMDY(date1, &date1MDY) == FAILURE)
            goto done;

        if (JpmcdsDateToMDY(date2, &date2MDY) == FAILURE)
            goto done;

        Y1 = date1MDY.year;
        M1 = date1MDY.month;
        D1 = date1MDY.day;

        Y2 = date2MDY.year;
        M2 = date2MDY.month;
        D2 = date2MDY.day;

        /* D1=31 => change D1 to 30 */
        if (D1 == 31)
            D1 = 30;

        /* D2=31 => change D2 to 30 */
        if (D2 == 31)
            D2 = 30;

        *result = (Y2-Y1)*360 + (M2-M1)*30 + (D2-D1);
        break;

      case JPMCDS_ACT_365:
      case JPMCDS_ACT_365F:
      case JPMCDS_ACT_360:
        *result = date2 - date1;
        break;

      default:
        JpmcdsErrMsg("%s: Invalid method (%ld).\n", routine, method);
        goto done;

    }

    /* negative difference ??? */
    if (negative == 1)
        *result *= -1;
    status = SUCCESS;

 done:
   if (status == FAILURE)
        JpmcdsErrMsg("%s: Failed.\n", routine);

    return status;
}


/*
***************************************************************************
** Given a date in TDate, retrieve beginning of the year in TDate format
***************************************************************************
*/
int JpmcdsYearStart(
                TDate  DATE,            /* input date */
                TDate  *result          /* output           */
                )
{
    TMonthDayYear  next_year;

    if (JpmcdsDateToMDY(DATE,&next_year) == FAILURE)
        return FAILURE;

    next_year.month=1;
    next_year.day=1;

    if (JpmcdsMDYToDate(&next_year,result) == FAILURE)
        return FAILURE;

    return SUCCESS;
}


/*
***************************************************************************
** Given a date in TDate, retrieve year as a four digit year.
***************************************************************************
*/
int JpmcdsYear(
                TDate DATE,         /* input date */
                int *result         /* output           */
                )
{

    TMonthDayYear  temp;

    if (JpmcdsDateToMDY(DATE,&temp) == FAILURE)
        return FAILURE;

    *result=(int)temp.year;
    return 0;
}


/*
***************************************************************************
** Calculates the number of days left in a year.
***************************************************************************
*/
int JpmcdsDaysLeftThisYear(
                        TDate DATE,         /* input date */
                        int method,         /* day count method */
                        long *result        /* output           */
                        )
{
long year;
TDate next_year_days;
TMonthDayYear temp, next_year;

    if (JpmcdsDateToMDY(DATE,&temp) == FAILURE)
        return FAILURE;

    year=temp.year;
    next_year.year=year+1;
    next_year.month=1;
    next_year.day=1;

    if (JpmcdsMDYToDate(&next_year,&next_year_days) == FAILURE)
        return FAILURE;

    return JpmcdsDaysDiff(DATE,next_year_days,method,result);
}


/*
***************************************************************************
** Find minimum in an array of integers.
***************************************************************************
*/
static int mymin (int *values)
{
  int num, minimum_val,i;
  num=values[0];

  minimum_val=values[1];

  for (i=2;i<=num;i++)
    if (values[i]<minimum_val)
      minimum_val=values[i];

  return minimum_val;
}


/*
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
 TDate          *advAdjustedDate)    /* (O) Advanced adjusted date */
{
    static char  routine[] = "JpmcdsDateFwdThenAdjust";
    int   status = FAILURE;      /* Until successful */
    TDate advancedDate;

    /* Advance the date */
    if (JpmcdsDtFwdAny (date, interval, &advancedDate) == FAILURE)
        goto done;

    /* Adjust to find a good business day */
    if (JpmcdsBusinessDay (advancedDate, badDayMethod,
                        holidayFile, advAdjustedDate) == FAILURE)
        goto done;

    status = SUCCESS;

done:
    if (status == FAILURE)
        JpmcdsErrMsg ("%s: Failed.\n", routine);

    return status;
}

