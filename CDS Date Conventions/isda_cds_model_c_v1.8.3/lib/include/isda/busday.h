/*
 * ISDA CDS Standard Model
 *
 * Copyright (C) 2009 International Swaps and Derivatives Association, Inc.
 * Developed and supported in collaboration with Markit
 * 
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the ISDA CDS Standard Model Public License.
 */

#ifndef BUSDAY_H
#define BUSDAY_H

#include "cgeneral.h"
#include "cdate.h"
#include "buscache.h"

#ifdef __cplusplus
extern "C"
{
#endif

/* Bad Day conventions
 */
#define JPMCDS_BAD_DAY_FOLLOW   ((long)'F')   /* Use following good day */
#define JPMCDS_BAD_DAY_PREVIOUS ((long)'P')   /* Use previous good day */
#define JPMCDS_BAD_DAY_NONE     ((long)'N')   /* Dont check for bad days */
#define JPMCDS_BAD_DAY_MODIFIED ((long)'M')   /* Try following, then previous*/
                                           /* if a month boundary is crossed.*/
#define JPMCDS_NUM_BAD_DAY_CONVS 4         /* # bad day conventions */

/*m
***************************************************************************
** Determines whether a day is at a weekend as given by the weekends flag.
***************************************************************************
*/
#define JPMCDS_IS_WEEKEND(date, weekends) ((1 << ((date) % 7)) & (weekends))
#define JPMCDS_IS_WEEKDAY(date, weekends) (!(JPMCDS_IS_WEEKEND((date),(weekends))))


/*
***************************************************************************
** This file is broken into three types of functions.
**
** Note that the holiday cache manipulation functions are now in
** buscache.h.
**
** 1. The original interface. These are functions that use a string for
**    the holiday file name, and the software gets the holiday list from
**    the holiday cache.
**
** 2. Holiday list manipulation functions. These are functions that
**    manipulate the holiday list structure but perform no significant
**    analytics.
**
** 3. The new interface. These are functions that use a holiday list
**    structure and perform the actual holiday adjustments.
***************************************************************************
*/


/*
***************************************************************************
** 1. The original interface. These are functions that use a string for
**    the holiday file name, and the software gets the holiday list from
**    the holiday cache.
***************************************************************************
*/

/*f
***************************************************************************
** Calculates a business date being at offset business days
** from the original date
***************************************************************************
*/
int  JpmcdsDateFromBusDaysOffset
(TDate     fromDate,       /* (I) input date                 */
 long      offset,         /* (I) number of business days    */
 char     *holidayFile,    /* (I) holiday file specification */
 TDate    *result);        /* (O) resulting business date    */


/*f
***************************************************************************
** Using business day conventions (Following, Preceding, and Modified
** Following) calculates the next bisiness day if the input date is not a
** business day. Input and Output dates are represented as TDate types.
***************************************************************************
*/
int  JpmcdsBusinessDay
(TDate   date,                       /* (I) Input Date */
 long    method,                     /* (I) See ldate.h */
 char   *holidayFile,                /* (I) Filename w/ Holidays */
 TDate  *outDate);                   /* (O) output date */

/*f
***************************************************************************
** Validates a bad day convention.
***************************************************************************
*/
int JpmcdsBadDayConvValid(char *routine, long badDayConv);


/*f
***************************************************************************
** Indicates whether a date is a business day.
***************************************************************************
*/
int JpmcdsIsBusinessDay
(TDate      date,                   /* (I) Input Date */
 char      *name,                   /* (I) Filename w/ Holidays */
 TBoolean  *isBusinessDay);         /* (O) */


/*f
***************************************************************************
** Computes the last business day of the month.
***************************************************************************
*/
int JpmcdsDateToBusinessEOM
(TDate  inDate,  /* (I) Date to be tested */
 char  *name,    /* (I) Name of holiday file */
 TDate *outDate  /* (O) TRUE or FALSE returned */
);


/*f
***************************************************************************
** Loads a holiday cache entry from a disk file, removing old one if it
** exists.
***************************************************************************
*/
EXPORT int JpmcdsHolidayLoadFromDisk
(char  *name,                        /* (I) name associated with holidays */
 char  *filename);                   /* (I) filename to load */


/*
***************************************************************************
** 2. The new interface. These are functions that use a holiday list
**    structure and perform the actual holiday adjustments.
***************************************************************************
*/


/*f
***************************************************************************
** Equivalent to JpmcdsIsBusinessDay but uses THolidayList structure in place
** of holiday file name.
**
** Indicates whether a date is a business day.
***************************************************************************
*/
int JpmcdsHolidayListIsBusinessDay
(TDate         date,              /* (I) Input Date */
 THolidayList *hl,                /* (I) Holiday list structure */
 TBoolean     *isBusinessDay);    /* (O) TRUE if a business day */


/*f
***************************************************************************
** Equivalent to JpmcdsIsHoliday, but uses a THolidayList structure instead of
** a holiday file name (and output type is TBoolean rather than long).
**
** Does not take week-ends into account.
***************************************************************************
*/
int JpmcdsHolidayListIsHoliday
(TDate         date,       /* (I) Arbitrary date          */
 THolidayList *hl,         /* (I) Holiday list structure  */
 TBoolean     *isHoliday); /* (O) 0 = not a holiday, 1 = is a holiday */

/*f
***************************************************************************
** Converts date to a valid business date using the holiday list and bad
** day convention passed in.
**
** Returns SUCCESS if convention OK, FAILURE otherwise.
**
** Equivalent to JpmcdsBusinessDay except using a THolidayList structure
** instead of a holiday file. Should be used if efficiency is an issue.
***************************************************************************
*/
int JpmcdsHolidayListBusinessDay
(TDate         date,       /* (I) Arbitrary date. */
 long          badDayConv, /* (I) Bad day convention for adjusting
                              non-business days. Use one of the following:
                              JPMCDS_BAD_DAY_FOLLOW
                              JPMCDS_BAD_DAY_PREVIOUS
                              JPMCDS_BAD_DAY_MODIFIED
                              JPMCDS_BAD_DAY_NONE */
 THolidayList *hl,         /* (I) Holiday file list. */
 TDate        *outDate);   /* (O) Valid business day. */

/*
***************************************************************************
** FUNCTION: JpmcdsMultiHolidayListBusinessDay
**
** (See JpmcdsHolidayListBusinessDay.)
***************************************************************************
*/
int JpmcdsMultiHolidayListBusinessDay
(TDate         date,             /* (I) Arbitrary date. */
 long          badDayConv,       /* (I) Bad day convention for adjusting
                                    non-business days. Use one of the following:
                                    JPMCDS_BAD_DAY_FOLLOW
                                    JPMCDS_BAD_DAY_PREVIOUS
                                    JPMCDS_BAD_DAY_MODIFIED
                                    JPMCDS_BAD_DAY_NONE */
 int           numHolidayLists,  /* (I) Size of holidayLists[] array. */
 THolidayList *holidayLists[],   /* (I) Holiday lists. */
 TDate        *outDate);         /* (O) Valid business day. */

/*f
***************************************************************************
** Calculates the number of business days between two dates (FROM & TO).
**
** Algorithm:
** 1. if FROM = TO, the result is 0
** 2. if FROM < TO, the result is the number of business days in the CLOSED
**    interval of [FROM+1,TO]
** 3. if FROM > TO, the result is negated number of business days in the
**    CLOSED interval of [TO,FROM-1]
**
** Equivalent to JpmcdsBusinessDaysDiff but uses a THolidayList structure
** instead of a holiday file name. Should be used if efficiency is an
** issue.
***************************************************************************
*/
int JpmcdsHolidayListBusinessDaysDiff
    (TDate         fromDate,     /* (I) Earlier date       */
     TDate         toDate,       /* (I) Later date         */
     THolidayList *hl,           /* (I) Holiday list structure  */
     long         *result);      /* (O) Result value       */


/*f
***************************************************************************
** Adds a number of business days to the original date.
**
** Equivalent to JpmcdsDateFromBusDaysOffset, except using a THolidayList
** structure instead of a holiday file name. Should be used if efficiency
** is an issue.
***************************************************************************
*/
int JpmcdsHolidayListAddBusinessDays
    (TDate         fromDate,     /* (I) Input date                 */
     long          numBusDays,   /* (I) Number of business days    */
     THolidayList *hl,           /* (I) Holiday list structure     */
     TDate        *resultDate);  /* (O) Resulting business date    */


/*f
***************************************************************************
** Computes the last business day of the month.
***************************************************************************
*/
int JpmcdsHolidayListDateToBusinessEOM
(TDate         inDate,  /* (I) Date to be tested */
 THolidayList *hl,    /* (I) Name of holiday file */
 TDate        *outDate  /* (O) TRUE or FALSE returned */
);



#ifdef __cplusplus
}
#endif

#endif    /* BUSDAY_H */

