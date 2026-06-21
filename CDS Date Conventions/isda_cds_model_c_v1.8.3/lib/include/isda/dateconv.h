/*
 * ISDA CDS Standard Model
 *
 * Copyright (C) 2009 International Swaps and Derivatives Association, Inc.
 * Developed and supported in collaboration with Markit
 * 
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the ISDA CDS Standard Model Public License.
 */

#ifndef DATECONV_H
#define DATECONV_H

#include "cgeneral.h"
#include "cdate.h"
#include "mdydate.h"

#ifdef __cplusplus
extern "C"
{
#endif

#define JPMCDS_SUNDAY     0                /* constants for day-of-week */
#define JPMCDS_MONDAY     1
#define JPMCDS_TUESDAY    2
#define JPMCDS_WEDNESDAY  3
#define JPMCDS_THURSDAY   4
#define JPMCDS_FRIDAY     5
#define JPMCDS_SATURDAY   6

#define JPMCDS_TDATE_BASE_YEAR 1601

/* A generally useful macro. */
#define JPMCDS_IS_LEAP(year) (                                              \
                            (((year)%4 == 0) && ((year)%100 != 0)) || \
                            ((year)%400 == 0)                            \
                          )


/*f
***************************************************************************
** Converts year, month, day to a date and returns the date. This is useful
** for initializing dates within test routines. Returns FAILURE for invalid
** dates.
***************************************************************************
*/
TDate JpmcdsDate
(long year,  /* (I) Year */
 long month, /* (I) Month */
 long day    /* (I) Day */
);


/*f
***************************************************************************
** Converts TDate to Month, Day, Year.
**
** NOTE: routine was written to avoid divides especially as well as 
** multiplies. Whenever possible, the type int is used-which is presumably 
** the fastest type on a given architecture.
***************************************************************************
*/
int JpmcdsDateToMDY
    (TDate date,                        /* (I) TDate format */
     TMonthDayYear *mdyDate);           /* (O) Date in mm/dd/yyyy format */


/*f
***************************************************************************
** Converts Month, Day, Year to TDate.
**
** NOTE: routine was written to avoid divides (especially) and multiplies.
***************************************************************************
*/
int JpmcdsMDYToDate
    (TMonthDayYear *mdyDate,            /* (I) Date in mm/dd/yyyy format */
     TDate *date);                      /* (O) TDate format */


/*f
***************************************************************************
** Converts TDate to day of week (0-6)
** 0 is Sunday, 1 is Monday.... 6 is Saturday
***************************************************************************
*/
int JpmcdsDayOfWeek
    (TDate date,                        /* (I) Date */
     long *dayOfWeek);                  /* (O) 0=Sunday, 6=Saturday */


/*f
***************************************************************************
** Normalizes a month/day/year. If month is out of range, it is brought 
** into range, and the years are incremented or decremented as appropriate.
** If day belongs to a month/year combination which does not exist, (such 
** as April 31) the day is reduced so that it becomes valid.(to April 30).
***************************************************************************
*/
int JpmcdsNormalizeMDY
    (TMonthDayYear *mdy);               /* (I/O) */


/*f
***************************************************************************
** Returns # days in month.
***************************************************************************
*/
int JpmcdsDaysInMonth
    (long year,                         /* (I) Year */
     long month);                       /* (I) Month (1-12) */


/*f
***************************************************************************
** Moves a date to the end of the current month.
***************************************************************************
*/
int JpmcdsDateToEOM
(TDate  inDate,  /* (I) Input date */
 TDate *outDate  /* (O) Output date - at end of month, in the same month */
);


#ifdef __cplusplus
}
#endif

#endif    /* DATECONV_H */

