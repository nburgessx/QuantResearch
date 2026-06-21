/*
 * ISDA CDS Standard Model
 *
 * Copyright (C) 2009 International Swaps and Derivatives Association, Inc.
 * Developed and supported in collaboration with Markit
 * 
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the ISDA CDS Standard Model Public License.
 */

#ifndef CDATE_H
#define CDATE_H

#include "cgeneral.h"

#ifdef __cplusplus
extern "C"
{
#endif

typedef long int TDate;

/*
 * (double) value of TDate with any fractional part
 * being a fraction of elapsed day.
 */
typedef double TDateTime;

/* Date constants */
#define MAX_DAYS_PER_MONTH         31
#define JPMCDS_MONTHS_PER_YEAR        12
#define JPMCDS_MONTHS_PER_SEMI         6
#define JPMCDS_MONTHS_PER_QUARTER      3
#define JPMCDS_DAYS_PER_WEEK           7
#define JPMCDS_DAYS_PER_YEAR         365
#define JPMCDS_YEARS_PER_DAY          (1./365.)
#define JPMCDS_DAYS_PER_LUNAR_MONTH   28
#define JPMCDS_WEEKS_PER_LUNAR_MONTH   4

/* --------- Definition of a date interval -------------------- */
typedef struct
{
    int prd;        /* number of periods from offset date                  */
    char prd_typ;   /* type of periods                                     */
                    /* D - day; M - month; W - week                        */
                    /* Q - 3 months; S - 6 months                          */
                    /* A - 12 months; Y - 12 months                        */
                    /* I - quarterly IMM period                            */
                    /* F - flexible end of month                           */
                    /* E - end of month                                    */
                    /* J - monthly IMM period                              */
                    /* K - quarterly Australian futures period             */
                    /* L - quarterly New Zealand (kiwi) futures period     */
                    /* T - equity derivatives expiry - 3rd Friday monthly  */
                    /* U - Lunar (i.e. 28 Day) period                      */

    int flag;       /* 0 - offset is value date
                       -1 - offset is the previous date in the date array
                       x - any other number is index into array of intervals.
                           the date at that location is an offset */
} TDateInterval;


/*
 * SET_TDATE_INTERVAL macro.
 */
#define SET_TDATE_INTERVAL(di,periods,period_type)   \
do {                                                 \
(di).prd = (periods);                                \
(di).prd_typ = (period_type);                        \
(di).flag = 0;                                       \
} while(0)

#define TDATE_INTERVALS_EQUAL(d1,d2)            \
   ((d1).prd == (d2).prd && (d1).prd_typ == (d2).prd_typ)

#ifdef __cplusplus
}
#endif

/*t
 *  TDateList is a list of dates.
 */
typedef struct
{
    int    fNumItems;
    TDate *fArray;
} TDateList;

/*t
 * TDateAdjIntvl is a time interval expressed in business
 * or calendar days or week days.
 */
#define JPMCDS_DATE_ADJ_TYPE_CALENDAR 0
#define JPMCDS_DATE_ADJ_TYPE_BUSINESS 1
#define JPMCDS_DATE_ADJ_TYPE_WEEKDAY  2

/*t
** TDateAdjIntvl is a holiday-adjusted date interval.
**
** Note that isBusDays used to be TBoolean. Now it can have three
** possible values - the third value (JPMCDS_DATE_ADJ_TYPE_WEEKDAY)
** indicates that we move forward a number of weekdays and then
** adjust for bad days.
**
** The interval can be daily based or monthly based. If it is
** monthly based and the isBusDays flag is not calendar, then
** this indicates that business end of month adjustment applies,
** i.e. if the start date is the last business day of the month,
** then the end date is the chosen to be the last day of the month
** (which can subsequently be adjusted via the badDayConv to be
** the last business day of the month, but this is not done
** automatically).
*/
typedef struct
{
    TDateInterval interval;    /* Must be in days if isBusDays=T */
    int           isBusDays;   /* see JPMCDS_DATE_ADJ_TYPE_... constants */
    char         *holidayFile; /* Holiday file specification */
    long          badDayConv;  /* Only applies if isBusDays=F */
} TDateAdjIntvl;


#define JPMCDS_SET_ADJ_INTERVAL_DAYS(adjIntval, days) \
adjIntval.interval.prd = days;             \
adjIntval.interval.prd_typ = 'D';          \
adjIntval.isBusDays = JPMCDS_DATE_ADJ_TYPE_CALENDAR;  \
adjIntval.holidayFile = NULL;              \
adjIntval.badDayConv = JPMCDS_BAD_DAY_NONE;

#endif    /* CDATE_H */
