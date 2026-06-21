/*
 * ISDA CDS Standard Model
 *
 * Copyright (C) 2009 International Swaps and Derivatives Association, Inc.
 * Developed and supported in collaboration with Markit
 * 
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the ISDA CDS Standard Model Public License.
 */

#include "busday.h"
#include <ctype.h>
#include "bsearch.h"
#include "cerror.h"
#include "dateconv.h"
#include "dtlist.h"
#include "macros.h"
#include "cfileio.h"


/*---------------------------------------------------------------------------
 *                  MACROS used in this file.
 *---------------------------------------------------------------------------
 */

#if defined(SIGN)
#   undef   SIGN
#endif
#define SIGN(x)                      ((x) < 0 ? (-1) : (1))


#define IS_EMPTY_DATELIST(pdatelist) ( (pdatelist)->fNumItems > 0  ? \
                                       FALSE : \
                                       TRUE)

#if defined(DATE_BUFFER_LEN)
#   undef   DATE_BUFFER_LEN
#endif
#define DATE_BUFFER_LEN         20     /* maximum length of string date to read in */

#if defined(HOLC_BUFFER_LEN)
#   undef   HOLC_BUFFER_LEN
#endif
#define HOLC_BUFFER_LEN         20     /* maximum length of holiday centre to read in */

#if defined(INITIAL_ARRAY_SIZE)
#   undef   INITIAL_ARRAY_SIZE
#endif
#define INITIAL_ARRAY_SIZE      2000    /* starting size of array for holding dates for a single holiday centre */

#if defined(INCREMENT_ARRAY_SIZE)
#   undef   INCREMENT_ARRAY_SIZE
#endif
#define INCREMENT_ARRAY_SIZE      64     /* increment for increasing array size */


/*
***************************************************************************
** Private function declarations
***************************************************************************
*/
static int   holidayCompare (const void *e1, const void *e2);

static int getNextBusDateMulti
(TDate          startDate,        /* (I) starting date. */
 long           direction,        /* (I) +1=forwards, -1=backwards */
 int            numHolidayLists,  /* (I) Size of holidayLists[] array. */
 THolidayList  *holidayLists[],   /* (I) Holiday lists. */
 TDate         *nextDate);        /* (O) next business day */

static int   getNextBusDate (
    TDate          startDate,    /* (I) starting date. */
    long           direction,    /* (I) +1=forwards, -1=backwards */
    THolidayList * hl,           /* (I) holiday list */
    TDate        * nextDate      /* (O) next business day */
    );

static int   forwardNonStandardWeekends (
    TDate          fromDate,       /* (I) start date */
    long           numBusDaysLeft, /* (I) abs. num. bus. days */
    long           intervalSign,   /* (I) +1=forwards, -1=backwards */
    long           weekends,       /* (I) Weekends flags */
    long           busDaysPerWeek, /* (I) num. bus. days per week */
    TDateList    * hlDatelist,     /* (I) holiday dateList */
    TDate        * resultDate      /* (O) forwarded date. */
);


static int   holAdjustStandardWeekends (
    TDate           fromDate,       /* (I) start date */
    TDateList     * hlDatelist,     /* (I) holiday dateList */
    TDate         * resultDate      /* (I/O) forwarded date. */
);

static long  calcNumWeekdayHolidays (
    TDate              toDate,     /* (I) End date.   */
    long               startHolIdx,/* (I) Starting holiday index. */
    long               direction,  /* (I) 1=forwards, 2=backwards*/
    long               weekends,   /* (I) Weekends flags */
    TDateList        * hlDatelist, /* (I) Holiday datelist. */
    long             * endHolIdx   /* (O) Idx where hol=toDate, +1 */
);

static int   findFirstHolidayIdx (
    TDate             date,         /* (I) Date to search for. */
    TDateList       * hlDatelist,   /* (I) Holiday date list. */
    long              direction,    /* (I) +1=forwards, -1=backwards. */
    long            * holIdx,       /* (O) Index of first holiday. */
    TBoolean        * doneSearching /* (O) TRUE=date is beyond list bounds. */
);

static int readHolidaysAll
(char                 *fileName,              /* (I) file to read                 */
 const TBoolean       supplementaryDates,     /* (I) are these extra dates or not */
 const TBoolean       addSupplementaryDates   /* (I) supplemantary dates if any are add (TRUE) or subtract (FALSE) */
);

static int stringToDate(
    const char* str,            /* (I) string to convert */
    long  *date                 /* (O) date              */
);

/*---------------------------------------------------------------------------
 *                              Static Data.
 *---------------------------------------------------------------------------
 */

/*
 * Table containing the offset to add to a date to skip a given number of
 * business days.  Only valid for less than 5 business days, and when saturday
 * and sunday are both holidays.
 *
 * So:
 * offsetTable[dayOfWeek(date)][nDays] is the number of days to add to date
 * to get a new business day.
 *
 * Note: The indexes are set up as TDate mod 7, which is slightly different
 * than JpmcdsDayOfWeek(), which adds one to TDate.
 */
static long  offsetTable[7][5] =
                                    {/*n:  0  1  2  3  4 */

                                        {  0, 1, 2, 3, 4 },  /* Monday */
                                        {  0, 1, 2, 3, 6 },  /* Tuesday */
                                        {  0, 1, 2, 5, 6 },  /* Wednesday */
                                        {  0, 1, 4, 5, 6 },  /* Thursday */
                                        {  0, 3, 4, 5, 6 },  /* Friday */
                                        { -1, 2, 3, 4, 5 },  /* Saturday */
                                        { -2, 1, 2, 3, 4 }   /* Sunday */
                                    };
/*---------------------------------------------------------------------------
 *                              End Static Data.
 *---------------------------------------------------------------------------
 */


/*
***************************************************************************
** This file is broken into two types of functions.
**
** Note that holiday cache and non-analytics holiday list functions have
** been moved to buscache.c.
**
** 1. The original interface. These are functions that use a string for
**    the holiday file name, and the software gets the holiday list from
**    the holiday cache.
**
** 2. The new interface. These are functions that use a holiday list
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

/*
***************************************************************************
** Calculates a business date being at offset business days
** from the original date
***************************************************************************
*/
int  JpmcdsDateFromBusDaysOffset
(TDate     fromDate,       /* (I) input date                 */
 long      offset,         /* (I) number of business days    */
 char     *holidayFile,    /* (I) holiday file specification */
 TDate    *result)         /* (O) resulting business date    */
{
    static char          routine[]="JpmcdsDateFromBusDaysOffset";
    int                  status = FAILURE;

    THolidayList       * hl = NULL;

    /*
    ** Search for the holiday file in the cache.
    */
    hl = JpmcdsHolidayListFromCache (holidayFile);
    if (hl == NULL)
        goto done;

    /* call the real routine, and return status */
    if (JpmcdsHolidayListAddBusinessDays (fromDate, offset, hl, result) != SUCCESS)
        goto done;

    status = SUCCESS;

done:
    if (status != SUCCESS)
        JpmcdsErrMsg ("%s: Failed.\n", routine);

    return status;
}


/*
***************************************************************************
** Using business day conventions (Following, Preceding, and Modified
** Following) calculates the next bisiness day if the input date is not a
** business day. Input and Output dates are represented as TDate types.
***************************************************************************
*/
int  JpmcdsBusinessDay
(TDate   date,                       /* Input Date */
 long    method,                     /* See ldate.h */
 char   *holidayFile,                /* Filename w/ Holidays */
 TDate  *outDate)                    /* output date */
{
    static char   routine[] = "JpmcdsBusinessDay";
    THolidayList *hl = NULL;
    int           status = FAILURE;

    /* determine whether we should do anything */
    if (method == JPMCDS_BAD_DAY_NONE)
    {
        *outDate = date;
        status = SUCCESS;
        goto done; /* success */
    }

    /* if method is other than NONE, pass to hl routine */
    hl = JpmcdsHolidayListFromCache (holidayFile);
    if (hl == NULL)
        goto done;

    /* call the real routine, and return status */
    if (JpmcdsHolidayListBusinessDay (date, method, hl, outDate) != SUCCESS)
        goto done;

    status = SUCCESS;

done:

    if (status != SUCCESS)
        JpmcdsErrMsg ("%s: Failed.\n", routine);

    return status;
}


/*
***************************************************************************
** Validates a bad day convention.
***************************************************************************
*/
int  JpmcdsBadDayConvValid
(char *routine,
 long  badDayConv)
{
    switch (badDayConv)
    {
    case JPMCDS_BAD_DAY_FOLLOW:
    case JPMCDS_BAD_DAY_PREVIOUS:
    case JPMCDS_BAD_DAY_NONE:
    case JPMCDS_BAD_DAY_MODIFIED:
        return SUCCESS;

    default:
        /*
        ** There's a good chance that the user entered a lower-case
        ** version of a valid bad day convention.
        */
        switch (toupper((int)badDayConv))
        {
        case JPMCDS_BAD_DAY_FOLLOW:
        case JPMCDS_BAD_DAY_PREVIOUS:
        case JPMCDS_BAD_DAY_NONE:
        case JPMCDS_BAD_DAY_MODIFIED:
            JpmcdsErrMsg ("%s: Bad day convention %ld must be uppercase\n",
                              routine, badDayConv);
            break;
        default:
            JpmcdsErrMsg ("%s: unknown bad day convention %ld\n",
                       routine, badDayConv);
            break;
        }
        return FAILURE;
    }
}




/*
***************************************************************************
** Indicates whether a date is a business day.
***************************************************************************
*/
int  JpmcdsIsBusinessDay
(TDate      date,                   /* (I) Input Date */
 char      *name,                   /* (I) Filename w/ Holidays */
 TBoolean  *isBusinessDay)          /* (O) */
{
    static char   routine[] = "JpmcdsIsBusinessDay";
    THolidayList *hl = NULL;
    int           status = FAILURE;

    hl = JpmcdsHolidayListFromCache (name);
    if (hl == NULL)
        goto done;

    /* call the real routine, and return status */
    if (JpmcdsHolidayListIsBusinessDay (date, hl, isBusinessDay) != SUCCESS)
        goto done;

    status = SUCCESS;

done:

    if (status != SUCCESS)
        JpmcdsErrMsg ("%s: Failed.\n", routine);

    return status;
}


/*
***************************************************************************
** Computes the last business day of the month.
***************************************************************************
*/
int JpmcdsDateToBusinessEOM
(TDate  inDate,  /* (I) Date to be tested */
 char  *name,    /* (I) Name of holiday file */
 TDate *outDate  /* (O) TRUE or FALSE returned */
)
{
    static char   routine[] = "JpmcdsDateToBusinessEOM";
    int           status = FAILURE;

    THolidayList *hl;

    hl = JpmcdsHolidayListFromCache (name);
    if (hl == NULL)
        goto done;

    if (JpmcdsHolidayListDateToBusinessEOM (inDate, hl, outDate) != SUCCESS)
        goto done;

    status = SUCCESS;

done:

    if (status != SUCCESS)
        JpmcdsErrMsgFailure(routine);

    return status;
}


/*
***************************************************************************
** Loads a holiday cache entry from a disk file, removing old one if it
** exists.
***************************************************************************
*/
EXPORT int JpmcdsHolidayLoadFromDisk
(char  *name,                        /* (I) name associated with holidays */
 char  *filename)                    /* (I) filename to load */
{
    static char routine[] = "JpmcdsHolidayLoadFromDisk";
    int         status    = FAILURE;

    THolidayList *hl;

    hl = JpmcdsHolidayListRead (filename);
    if (hl == NULL)
        goto done;

    if (JpmcdsHolidayListAddToCache (name, hl) != SUCCESS)
        goto done;

    status = SUCCESS;

done:

    if (status != SUCCESS)
        JpmcdsErrMsg ("%s: Failed.\n", routine);

    return status;
}


/*
***************************************************************************
** 2. The new interface. These are functions that use a holiday list
**    structure and perform the actual holiday adjustments.
***************************************************************************
*/


/*
***************************************************************************
** Equivalent to JpmcdsIsBusinessDay but uses THolidayList structure in place
** of holiday file name.
**
** Indicates whether a date is a business day.
***************************************************************************
*/
int  JpmcdsHolidayListIsBusinessDay
(TDate         date,              /* (I) Input Date */
 THolidayList *hl,                /* (I) Holiday list structure */
 TBoolean     *isBusinessDay)      /* (O) TRUE if a business day */
{
    static char     routine[] = "JpmcdsHolidayListIsBusinessDay";
    TBoolean        aHoliday = FALSE;

    /* Must have a holiday list */
    if (hl == NULL)
    {
        JpmcdsErrMsg ("%s: hl is NULL.\n", routine);
        return FAILURE;
    }

    /* First check for week-ends */
    if (JPMCDS_IS_WEEKEND (date, hl->weekends))
    {
        *isBusinessDay = FALSE;
        return SUCCESS;
    }

    /* Check if it is a holiday. */
    if (JpmcdsHolidayListIsHoliday (date, hl, &aHoliday) == FAILURE)
    {
        JpmcdsErrMsg("%s: Failed.\n", routine);
        return FAILURE;
    }

    *isBusinessDay = !aHoliday;
    return SUCCESS;
}



/*
***************************************************************************
** Equivalent to JpmcdsIsHoliday, but uses a THolidayList structure instead of
** a holiday file name (and output type is TBoolean rather than long).
**
** Does not take week-ends into account.
***************************************************************************
*/
int  JpmcdsHolidayListIsHoliday
(TDate         date,      /* (I) Arbitrary date          */
 THolidayList *hl,        /* (I) Holiday list structure  */
 TBoolean     *isHoliday) /* (O) 0 = not a holiday, 1 = is a holiday */
{
    static char   routine[] = "JpmcdsHolidayListIsHoliday";

    TDate       * dates;
    TDate       * bRet;

    *isHoliday = FALSE;

    if (hl == NULL || hl->dateList == NULL)
    {
        JpmcdsErrMsg("%s: hl is NULL.\n", routine);
        return FAILURE;
    }

    if (IS_EMPTY_DATELIST( hl->dateList ))
    {
        return SUCCESS;
    }

    dates = hl->dateList->fArray ;
    if (dates == NULL)
    {
        /* Should only be possible if numItems == zero. */
        PROGRAM_BUG();
        return FAILURE;
    }

    /*
    ** Do a binary search of the date list to see if date
    ** is in the list.
    */
    bRet = (TDate *) bsearch ((void *) &date,
                              (void *) &dates[0],
                              (size_t) hl->dateList->fNumItems,
                              (size_t) sizeof(TDate),
                              holidayCompare);

    if (bRet != NULL)                  /* found == holiday */
    {
        *isHoliday = TRUE;
    }


    return SUCCESS;
}


/****************************************************************************
 * Comparison function for the binary search done in JpmcdsIsHoliday and
 * findFirstHolidayIdx. This one is for forwards comparisons.
 ****************************************************************************/
static int holidayCompare (const void *e1, const void *e2)
{
    TDate   dateDiff;

    dateDiff = *((TDate *)e1) - *((TDate *)e2);

    /*
    ** change date difference to 1, -1, or 0.
    ** Makes 16 bit windows happy, but it slows
    ** things down noticibly when this function
    ** is called more than 1 million times over time.
    */
    if (dateDiff < 0)
        return (-1);
    else if (dateDiff > 0)
        return (1);
    else
        return (0);
}

/*
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
int  JpmcdsHolidayListBusinessDay
(TDate         date,       /* (I) Arbitrary date. */
 long          badDayConv, /* (I) Bad day convention for adjusting
                              non-business days. Use one of the following:
                              JPMCDS_BAD_DAY_FOLLOW
                              JPMCDS_BAD_DAY_PREVIOUS
                              JPMCDS_BAD_DAY_MODIFIED
                              JPMCDS_BAD_DAY_NONE */
 THolidayList *hl,         /* (I) Holiday file list. */
 TDate        *outDate)    /* (O) Valid business day. */

{
    static char          routine[] = "JpmcdsHolidayListBusinessDay";
    int                  status;

    /* just pass call down to JpmcdsMultiHolidayListBusinessDay */

    status = JpmcdsMultiHolidayListBusinessDay(date, badDayConv, 1, &hl, outDate);

    if (status != SUCCESS)
        JpmcdsErrMsg("%s: Failed.\n", routine);

    return status;
}

/*
***************************************************************************
** (See JpmcdsHolidayListBusinessDay.)
***************************************************************************
*/
int  JpmcdsMultiHolidayListBusinessDay
(TDate         date,             /* (I) Arbitrary date. */
 long          badDayConv,       /* (I) Bad day convention for adjusting
                                    non-business days. Use one of the following:
                                    JPMCDS_BAD_DAY_FOLLOW
                                    JPMCDS_BAD_DAY_PREVIOUS
                                    JPMCDS_BAD_DAY_MODIFIED
                                    JPMCDS_BAD_DAY_NONE */
 int           numHolidayLists,  /* (I) Size of holidayLists[] array. */
 THolidayList *holidayLists[],   /* (I) Holiday lists. */
 TDate        *outDate)          /* (O) Valid business day. */

{
    static char          routine[] = "JpmcdsMultiHolidayListBusinessDay";
    int                  status = FAILURE;

    TMonthDayYear        mdy1, mdy2;
    TDate                nextDate = 0L;
    long                 intervalSign = 1;

    if (holidayLists == NULL)
    {
        JpmcdsErrMsg ("%s: hl is NULL.\n", routine);
        return FAILURE;
    }

    *outDate = date;

    switch (badDayConv)
    {
    case JPMCDS_BAD_DAY_NONE:
        /* Doesn't matter if business day or not. */
        return SUCCESS;

    case JPMCDS_BAD_DAY_FOLLOW:
        intervalSign = +1;
        if (getNextBusDateMulti (date, intervalSign, numHolidayLists, holidayLists, &nextDate)
            != SUCCESS)
            goto done;
        break;

    case JPMCDS_BAD_DAY_PREVIOUS:
        intervalSign = -1;
        if (getNextBusDateMulti (date, intervalSign, numHolidayLists, holidayLists, &nextDate)
            != SUCCESS)
            goto done;
        break;

    case JPMCDS_BAD_DAY_MODIFIED:
        /*
        ** Go forwards first. If you wind up in a different
        ** month, then go backwards.
        */
        intervalSign = +1;
        if (getNextBusDateMulti (date, intervalSign, numHolidayLists, holidayLists, &nextDate)
            != SUCCESS)
            goto done;

        if (JpmcdsDateToMDY (nextDate, &mdy1) != SUCCESS)
            goto done;

        if (JpmcdsDateToMDY (date, &mdy2) != SUCCESS)
            goto done;

        if (mdy1.month != mdy2.month)
        {
            /* Go backwards. */
            if (getNextBusDateMulti (date, -intervalSign, numHolidayLists, holidayLists, &nextDate)
                != SUCCESS)
                goto done;
        }
        break;

    default:
        JpmcdsErrMsg ("%s: Unrecognized bad day convention = %ld.\n", routine, badDayConv);
        goto done;
    }

    *outDate = nextDate;
    status = SUCCESS;

done:
    if (status != SUCCESS)
        JpmcdsErrMsg("%s: Failed.\n", routine);

    return status;
}

/*
***************************************************************************
** (See getNextBusDate.)
***************************************************************************
*/
static int getNextBusDateMulti
(TDate          startDate,        /* (I) starting date. */
 long           direction,        /* (I) +1=forwards, -1=backwards */
 int            numHolidayLists,  /* (I) Size of holidayLists[] array. */
 THolidayList  *holidayLists[],   /* (I) Holiday lists. */
 TDate         *nextDate          /* (O) next business day */
)
{
    static char          routine[] = "getNextBusDateMulti";
    int                  status = FAILURE;
    int                  n;
    TDate                adjDate;

    *nextDate = startDate;

    /* this loop executes getNextBusDate for each market in sequence; whenever weekends/
       holidays force an adjustment, the process recommences at the first market; eventually,
       a date will be reached which is a business day in all markets, whereupon the loop
       will execute to completion

       we can be certain that all of the skipped dates were weekends/holidays in at least
       one of the markets, i.e. output will be the _first_ suitable day */

    n = 0;
    while (n < numHolidayLists)
    {
        if ( getNextBusDate(*nextDate, direction, holidayLists[n], &adjDate) != SUCCESS )
        {
            goto done;
        }

        if ( adjDate != *nextDate )
        {
            *nextDate = adjDate;

            /* date has changed to acommodate this market; now
               have to recommence testing for all other markets */

            /* in fact, this date does not need retesting for the current market,
               and we exploit this in the special case below, which avoids a redundant
               duplicate test in "single holiday list" usage */

            n = ( n == 0 ) ? 1 : 0;

            /* ^-- if adjusted for first market, now move on to second; in all
               other cases, recommence testing from scratch */
        }
        else
        {
            n++; /* next market */
        }
    }

    status = SUCCESS;

done:
    if (status != SUCCESS)
        JpmcdsErrMsg("%s: Failed.\n", routine);

    return status;
}

/*
***************************************************************************
** Compute the date of next business day, given a start date and
** a holiday list. The function currently always returns SUCCESS.
***************************************************************************
*/
static int getNextBusDate
(TDate          startDate,    /* (I) starting date. */
 long           direction,    /* (I) +1=forwards, -1=backwards */
 THolidayList * hl,           /* (I) holiday list */
 TDate        * nextDate      /* (O) next business day */
)
{
    TBoolean               doneSearchingList;
    TBoolean               aHoliday;
    long                   holIdx;

    long                   numHols = hl->dateList->fNumItems;
    TDate                * holArray = hl->dateList->fArray;
    TDate                  curDate = startDate;

    if (IS_EMPTY_DATELIST(hl->dateList))
    {
        doneSearchingList = TRUE;
    }
    else
    {
        if (findFirstHolidayIdx (startDate, hl->dateList,
                                 direction, &holIdx, &doneSearchingList)
            != SUCCESS)
        {
            JpmcdsErrMsg ("getNextBusDate: Failed.\n");
            return FAILURE;
        }
    }

    /*
    ** holIdx should now point to the first holiday on or before
    ** (after) the start date. Now loop though each day until
    ** you reach the first business day.
    */

    aHoliday = TRUE;

    do
    {   /* Check if it's a holiday first. */
        if (!doneSearchingList && curDate == holArray[ holIdx ])
        {
            holIdx  += direction;
            curDate += direction;
            doneSearchingList =
              (TBoolean)(holIdx < 0 || holIdx >= numHols);
        }
        else
        {
          /* Day is not a holiday. Check if it's a week-end day. */
            if (JPMCDS_IS_WEEKEND (curDate, hl->weekends))
            {
                 curDate += direction;
            }
            else
            {
                 aHoliday = FALSE;
            }
        }
    }
    while ( aHoliday == TRUE );

    *nextDate = curDate;

    return SUCCESS;
}


/*
*****************************************************************************
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
****************************************************************************
*/
int  JpmcdsHolidayListBusinessDaysDiff
    (TDate         fromDate,     /* (I) Earlier date       */
     TDate         toDate,       /* (I) Later date         */
     THolidayList *hl,           /* (I) Holiday list structure  */
     long         *result)       /* (O) Result value       */
{

    static char           routine[] = "JpmcdsHolidayListBusinessDaysDiff";
    int                   busDaysPerWeek = -1;

    long                  nrExtraDays;
    long                  numWeeks = 0L;

    int                   signum = 1;  /* assume going forward by default */
    long                  numHolidays = 0L;
    TDate                 curDate = fromDate;

/*
 * The following tables are used to calculate the number of business days
 * between weekdays.
 * If (toDate > fromDate) use fwdDiffTable.
 * eg NumDays = fwdDiffTable[ fromDate % 7 ][ toDate % 7 ]
 *
 * If (fromDate > toDate) use -bwdDiffTable.
 * eg NumDays = bwdDiffTable[ toDate % 7 ][ fromDate % 7 ]
 *
 * The values of the table were chosen so that the current ALIB behaviour
 * was left unchanged. Note that because fwdDiffTable is not identical
 * to (-1 x bwdDiffTable), swapping fromDate and toDate in the call to this
 * function does not simply reverse the sign of the result.
 */
    static long           fwdDiffTable[7][7] =

                              {/*   Mo Tu We Th Fr Sa Su*/

                                  {  0, 1, 2, 3, 4, 4, 4 },  /* Monday */
                                  {  4, 0, 1, 2, 3, 3, 3 },  /* Tuesday */
                                  {  3, 4, 0, 1, 2, 2, 2 },  /* Wednesday */
                                  {  2, 3, 4, 0, 1, 1, 1 },  /* Thursday */
                                  {  1, 2, 3, 4, 0, 0, 0 },  /* Friday */
                                  {  1, 2, 3, 4, 5, 0, 0 },  /* Saturday */
                                  {  1, 2, 3, 4, 5, 5, 0 }   /* Sunday */
                              };

    static long           bwdDiffTable[7][7] =

                              {/*   Mo Tu We Th Fr Sa Su*/

                                  {  0,-1,-2,-3,-4,-5,-5 },  /* Monday */
                                  { -4, 0,-1,-2,-3,-4,-4 },  /* Tuesday */
                                  { -3,-4, 0,-1,-2,-3,-3 },  /* Wednesday */
                                  { -2,-3,-4, 0,-1,-2,-2 },  /* Thursday */
                                  { -1,-2,-3,-4, 0,-1,-1 },  /* Friday */
                                  { -0,-1,-2,-3,-4, 0, 0 },  /* Saturday */
                                  { -0,-1,-2,-3,-4,-5, 0 }   /* Sunday */
                              };

    if (hl == NULL)
    {
        JpmcdsErrMsg ("%s: hl is NULL.\n", routine);
        return FAILURE;
    }

    *result = 0L; /* TDate EQ long */

    if (fromDate == toDate)
        return SUCCESS;

    /* Set the direction. */
    if (toDate < fromDate)
    {
        signum = -1;
    }

    /*
    ** First, compute the date difference adjusting for weekends only.
    */


    /*
    ** Get the number of business days per week. In-line for speed.
    */
    switch (hl->weekends)
    {
    case JPMCDS_WEEKEND_STANDARD:
        busDaysPerWeek = 5;
        numWeeks = (toDate - fromDate) / JPMCDS_DAYS_PER_WEEK;
        curDate += JPMCDS_DAYS_PER_WEEK * numWeeks;
        if (curDate > toDate)  /* backwards */
        {
            nrExtraDays = bwdDiffTable [ toDate %7 ][ fromDate % 7 ];
        }
        else
        {
            nrExtraDays = fwdDiffTable [ fromDate %7 ][ toDate % 7 ];
        }
        *result = numWeeks * busDaysPerWeek + nrExtraDays;
        break;
    case JPMCDS_WEEKEND_NO_WEEKENDS:
        *result = toDate - fromDate;
        break;
    default:
        busDaysPerWeek = 7;
        if (hl->weekends & JPMCDS_WEEKEND_SUNDAY) busDaysPerWeek--;
        if (hl->weekends & JPMCDS_WEEKEND_MONDAY) busDaysPerWeek--;
        if (hl->weekends & JPMCDS_WEEKEND_TUESDAY) busDaysPerWeek--;
        if (hl->weekends & JPMCDS_WEEKEND_WEDNESDAY) busDaysPerWeek--;
        if (hl->weekends & JPMCDS_WEEKEND_THURSDAY) busDaysPerWeek--;
        if (hl->weekends & JPMCDS_WEEKEND_FRIDAY) busDaysPerWeek--;
        if (hl->weekends & JPMCDS_WEEKEND_SATURDAY) busDaysPerWeek--;
        numWeeks    = ABS((toDate - fromDate) / JPMCDS_DAYS_PER_WEEK);
        curDate    += (JPMCDS_DAYS_PER_WEEK * numWeeks * signum);
        nrExtraDays = 0L;

        while ( ABS(curDate - toDate) > 0 )
        {
            curDate += signum;
            if (JPMCDS_IS_WEEKDAY( curDate, hl->weekends ))
            {
                nrExtraDays++;
            }
        }

        (*result) = ((numWeeks * busDaysPerWeek) + nrExtraDays) * signum;
        break;
    }


    /*
    ** Now count the number of weekday holidays
    ** and adjust the date difference by the result. The
    ** call to calcNumWeekdayHolidays searches the
    ** holiday list at most once.
    */

    if (! IS_EMPTY_DATELIST( hl->dateList ))
    {
        TBoolean      doneSearching = FALSE;
        long          holIdx = -1;

        if (findFirstHolidayIdx (fromDate + signum, hl->dateList,
                                 signum, &holIdx, &doneSearching)
            != SUCCESS)
            return FAILURE;

        if (!doneSearching)
        {
            numHolidays     = calcNumWeekdayHolidays (toDate, holIdx,
                                          signum, hl->weekends,
                                          hl->dateList, &holIdx);
        }
    }

    (*result) -= numHolidays * signum;

    return SUCCESS;
}


/*
*****************************************************************************
** Adds a number of business days to the original date.
**
** Equivalent to JpmcdsDateFromBusDaysOffset, except using a THolidayList
** structure instead of a holiday file name. Should be used if efficiency
** is an issue.
****************************************************************************
*/
int  JpmcdsHolidayListAddBusinessDays
    (TDate         fromDate,     /* (I) Input date                 */
     long          numBusDays,   /* (I) Number of business days    */
     THolidayList *hl,           /* (I) Holiday list structure     */
     TDate        *resultDate)   /* (O) Resulting business date    */
{

    static char        routine[] = "JpmcdsHolidayListAddBusinessDays";
    int                status = FAILURE;

    long               intervalSign   = SIGN( numBusDays );
    long               numBusDaysLeft = ABS( numBusDays );
    long               busDaysPerWeek = -1;

    if (hl == NULL)
    {
        JpmcdsErrMsg ("%s: hl is NULL.\n", routine);
        goto done;
    }

    /*
    ** Get the number of business days per week. In-line for speed.
    */
    switch (hl->weekends)
    {
    case JPMCDS_WEEKEND_STANDARD:
        busDaysPerWeek = 5;
        break;
    case JPMCDS_WEEKEND_NO_WEEKENDS:
        busDaysPerWeek = 7;
        break;
    default:
        busDaysPerWeek = 7;
        if (hl->weekends & JPMCDS_WEEKEND_SUNDAY) busDaysPerWeek--;
        if (hl->weekends & JPMCDS_WEEKEND_MONDAY) busDaysPerWeek--;
        if (hl->weekends & JPMCDS_WEEKEND_TUESDAY) busDaysPerWeek--;
        if (hl->weekends & JPMCDS_WEEKEND_WEDNESDAY) busDaysPerWeek--;
        if (hl->weekends & JPMCDS_WEEKEND_THURSDAY) busDaysPerWeek--;
        if (hl->weekends & JPMCDS_WEEKEND_FRIDAY) busDaysPerWeek--;
        if (hl->weekends & JPMCDS_WEEKEND_SATURDAY) busDaysPerWeek--;
        break;
    }

    if (IS_EMPTY_DATELIST (hl->dateList) &&
        hl->weekends == JPMCDS_WEEKEND_NO_WEEKENDS)
    {
        /*
        ** No adjustments at all.
        */
        *resultDate = fromDate + (intervalSign * numBusDaysLeft);
    }
    else if (intervalSign == 1 && hl->weekends == JPMCDS_WEEKEND_STANDARD)
    {
        /* First, check for special case. The code for this case
         * was provided by Tom Robbins-Milne. It is much faster
         * than the general algorithm.
         */

        if (numBusDaysLeft == 0)
        {
            *resultDate = fromDate;
        }
        else
        {
           /*
            * Calculate result if no holidays are involved.
            *
            * Move forward over a week for every 5 business days.
            * Use a table for moving 0..4 bus days from each day of the week.
            */
            *resultDate = fromDate +  7 * (numBusDaysLeft / 5) +
                             offsetTable[ fromDate % 7 ][ numBusDaysLeft % 5 ];

           /*
            * Now adjust for any holidays
            */
            if (hl->dateList->fNumItems > 0)
            {
                if (holAdjustStandardWeekends (fromDate,
                                               hl->dateList,
                                               resultDate) != SUCCESS)
                    goto done;
            }
        }
    }
    else
    {
        if (forwardNonStandardWeekends (fromDate,
                                        numBusDaysLeft,
                                        intervalSign,
                                        hl->weekends,
                                        busDaysPerWeek,
                                        hl->dateList,
                                        resultDate) != SUCCESS)
            goto done;
    }

    status = SUCCESS;

done:
    if (status == FAILURE)
        JpmcdsErrMsg("%s: Failed.\n", routine);

    return status;
}


/*---------------------------------------------------------------------------
 * Start of private functions used by function JpmcdsHolidayListAddBusinessDays.
 *---------------------------------------------------------------------------
 */

/*
***************************************************************************
** This function does the following:
**
** If going forward in time, it finds the index of the 1st holiday date
** which is >= the input date.
**
** If going backwards in time, it finds the index of the first of holiday
** which is <= the input date.
**
** The function calls JpmcdsBSearchLong to do the search. It returns SUCCESS
** or FAILURE.
**
** Note that the function DOES NOT check if the input holiday list is
** empty or NULL.
***************************************************************************
*/
static int  findFirstHolidayIdx (
    TDate             date,         /* (I) Date to search for. */
    TDateList       * hlDatelist,   /* (I) Holiday date list. */
    long              direction,    /* (I) +1=forwards, -1=backwards. */
    long            * holIdx,       /* (O) Index of first holiday. */
    TBoolean        * doneSearching /* (O) TRUE=date is beyond list bounds. */
)
{
    int       idxLo;
    int       idxHi;

    long      numHols = hlDatelist->fNumItems;
    TDate   * holArray = hlDatelist->fArray;

    /* Check if the date is outside array bounds. */
    (*doneSearching) =
        (TBoolean) ( (direction ==  1 && holArray[numHols-1] < date) ||
                     (direction == -1 && holArray[0]         > date) );

    if (!(*doneSearching))
    {
        int      maxIdx;
        int      minIdx;

        /*
        ** Do a binary search of the date list to find the
        ** bracketing low and high index values.
        */
        JpmcdsBSearchCheckOrder(0);     /* Turn order checking OFF. */
        if (JpmcdsBSearchLong( (double) date,
                                     holArray,
                                     sizeof(TDate),
                                     numHols,
                                    &idxHi,
                                    &idxLo) != SUCCESS)
        {
            JpmcdsErrMsg("findFirstHolidayIdx: Failed.\n");
            return FAILURE;
        }
        JpmcdsBSearchCheckOrder(1);     /* Turn order checking back ON. */

        /* Clean up after JpmcdsBSearchLong because it does not
           work as advertised. */
        minIdx  = MIN( idxHi, idxLo );
        maxIdx  = MAX( idxHi, idxLo );

        if (date >=  holArray[maxIdx])
        {
            *holIdx = (long)maxIdx;
            if (direction == 1 && date > holArray[maxIdx])
            {
                (*doneSearching) = TRUE;
            }
        }
        else if (date <= holArray[minIdx])
        {
            *holIdx = (long)minIdx;
            if (direction == -1 && date < holArray[minIdx])
            {
                (*doneSearching) = TRUE;
            }
        }
        else
        {
            if (direction == 1)
                *holIdx = (long)idxLo;
            else
                *holIdx = (long)idxHi;
        }
    }
    return SUCCESS;
}


/*
***************************************************************************
** This function handles the special case where the number of business days
** to move forward by is very large.
***************************************************************************
*/
static int forwardNonStandardWeekends (
    TDate          fromDate,       /* (I) start date */
    long           numBusDaysLeft, /* (I) abs. num. bus. days */
    long           direction,      /* (I) +1=forwards, -1=backwards */
    long           weekends,       /* (I) Weekends flags */
    long           busDaysPerWeek, /* (I) num. bus. days per week */
    TDateList    * hlDatelist,     /* (I) holiday dateList */
    TDate        * resultDate      /* (O) forwarded date. */
)
{
    static char routine[] = "forwardNonStandardWeekends";
    int         status = FAILURE;

    long               numWeeks;

    long               numHolidays = 0L;
    TDate              curDate = fromDate;
    long               holIdx = -1;
    TBoolean           doneSearchingList = FALSE;

    TDate            * holArray = hlDatelist->fArray;
    long               numHols  = hlDatelist->fNumItems;

    /*
    ** First, adjust for weekends only, pretending there are no holidays.
    */
    numWeeks  = MAX(((numBusDaysLeft / busDaysPerWeek) - 1), 0);
    curDate  += (JPMCDS_DAYS_PER_WEEK * numWeeks * direction);

    if (numHols > 0)
    {
        /*
        ** Search the holiday list for the first holiday
        ** strictly after (if going forward in time) or strictly
        ** before (if going backward in time) the start date. Note
        ** that the holiday list is assumed to be sorted in increasing
        ** order.
        */
        if (findFirstHolidayIdx (fromDate + direction,
                                 hlDatelist,
                                 direction,
                                 &holIdx,
                                 &doneSearchingList) != SUCCESS)
        {
            goto done;
        }
    }
    else
    {
        doneSearchingList = TRUE; /* No holidays - no searching */
    }

    /*
    ** Count the number of weekday holidays starting from the
    ** current holiday index and going up to the curDate.
    ** Count weekday holidays only because holidays occurring
    ** on week-end days have been skipped by the previous calculation.
    */
    if (!doneSearchingList)
    {
        numHolidays = calcNumWeekdayHolidays (curDate,
                                              holIdx,
                                              direction,
                                              weekends,
                                              hlDatelist,
                                              &holIdx);
        doneSearchingList =
            (TBoolean) (holIdx  <  0  ||  holIdx  >=  numHols);
    }
    else
    {
        numHolidays = 0;
    }
    numBusDaysLeft -= ((busDaysPerWeek * numWeeks) - numHolidays);

    /*
    ** Now search one day at a time starting one day beyond
    ** the current date.
    */
    while (numBusDaysLeft > 0)
    {
        curDate += direction;

        /*
        ** Check if the day is a holiday first. If it is, don't decrement
        ** numBusDaysLeft.
        */
        if (!doneSearchingList && curDate  ==  holArray[ holIdx ])
        {
            holIdx  += direction;
            doneSearchingList =
                (TBoolean) (holIdx  <  0  ||  holIdx  >=  numHols);
        }
        else  /* Not a holiday. */
        {
            /*
            ** If the day is a weekday, then decrement numBusDaysLeft,
            ** otherwise continue looping.
            */
            if (JPMCDS_IS_WEEKDAY( curDate, weekends ))
            {
                numBusDaysLeft--;
            }
        }
    }

    (*resultDate) = curDate;

    status = SUCCESS;

done:

    if (status != SUCCESS)
        JpmcdsErrMsg ("%s: Failed.\n", routine);

    return status;
}


/*
***************************************************************************
** This function adjusts a forwarded date to take account of any holidays 
** which might have occured during the forward calculation.
**
** The algorithm assumes the following conditions:
**
** (1) Saturdays and Sundays are both holidays.
** (2) The direction is forwards in time.
**
** The function is optimized for speed.
***************************************************************************
*/
static int holAdjustStandardWeekends (
    TDate           fromDate,       /* (I) start date */
    TDateList     * hlDatelist,     /* (I) holiday dateList */
    TDate         * resultDate      /* (I/O) forwarded date. */
)
{
    static char routine[] = "holAdjustStandardWeekends";
    int         status    = FAILURE;

    TDate       dateEnd;      /* Date to return */
    int         i, idxLo;

    int         numHols  = hlDatelist->fNumItems;
    TDate       *holDates = hlDatelist->fArray;

    dateEnd = *resultDate;

    /*
    ** If there are any holidays, then we want to find the first holiday
    ** strictly after fromDate. There are some special cases.
    **
    ** 1. There are no holidays or fromDate >= last holiday.
    **
    **    In this case no further adjustments are necessary.
    **
    ** 2. fromDate < first holiday.
    **
    **    JpmcdsBSearchLong returns i=1 in this case instead of i=0.
    **
    ** 3. fromDate is a holiday.
    **
    **    JpmcdsBSearchLong could return i such that holDates[i] = fromDate.
    */
    if (numHols <= 0 ||
        fromDate >= holDates[numHols-1])
    {
        /* Special case 1 - no further adjustment is necessary. */
    }
    else
    {
        if (fromDate <= holDates[0])
        {
            i = 0;
        }
        else
        {
        /*
         * Handle any holidays there are. Use JpmcdsBSearchLong() to find date in
         * the datelist.
         */
            JpmcdsBSearchCheckOrder(0);
            if (JpmcdsBSearchLong( (double) fromDate,
                                holDates,
                                sizeof(TDate),
                                numHols,
                                &idxLo,
                                &i) == FAILURE)
            {
                goto done;
            }
            JpmcdsBSearchCheckOrder(1);
        }

        /*
         * Don't count holidays that are on the original date.
         */
        if (holDates[i] == fromDate)
        {
            i++;
        }

        while (i < numHols && holDates[i] <= dateEnd)
        {
            dateEnd += offsetTable[ dateEnd % 7 ][ 1 ];
            i++;
        }
    }

    (*resultDate) = dateEnd;

    status = SUCCESS;

done:

    if (status != SUCCESS)
        JpmcdsErrMsg ("%s: Failed.\n", routine);

    return status;
}


/*
***************************************************************************
** Compute the number of week-day holidays starting at the input holiday 
** index and continuing until the holiday list is at or beyond toDate.
**
** Week-day holidays are holidays which occur on a non-weekend day.
***************************************************************************
*/
static long calcNumWeekdayHolidays (
    TDate              toDate,     /* (I) End date.   */
    long               startHolIdx,/* (I) Starting holiday index. */
    long               direction,  /* (I) 1=forwards, 2=backwards*/
    long               weekends,   /* (I) Weekend flags */
    TDateList        * hlDatelist, /* (I) Holiday datelist. */
    long             * endHolIdx   /* (O) Idx where hol=toDate, +1 */
    )
{
    long       k;

    long       numHols = 0L;
    TDate    * dateArray = hlDatelist->fArray;
    long       numItems = hlDatelist->fNumItems;

    if (direction == 1)     /* forward */
    {
        for (k = startHolIdx;
             (k < numItems && dateArray[k] <= toDate);
             k++)
        {
            numHols++;
        }
        *endHolIdx = k;
        if (k < numItems   &&
            dateArray[k] == toDate)
        {
               (*endHolIdx)++;
        }
    }
    else                   /* backward */
    {
        for (k = startHolIdx;
             (k >= 0 && dateArray[k] >= toDate);
             k--)
        {
             numHols++;
        }
        *endHolIdx = k;
        if (k >= 0  &&
            dateArray[k] == toDate)
        {
              (*endHolIdx)--;
        }
    }
    return (numHols); /* SUCCESS */
}



/*---------------------------------------------------------------------------
 * End of private functions used by function JpmcdsHolidayListAddBusinessDays.
 *---------------------------------------------------------------------------
 */


/*
***************************************************************************
** Computes the last business day of the month.
***************************************************************************
*/
int JpmcdsHolidayListDateToBusinessEOM
(TDate         inDate,  /* (I) Date to be tested */
 THolidayList *hl,    /* (I) Name of holiday file */
 TDate        *outDate  /* (O) TRUE or FALSE returned */
)
{
    /*
    ** Calculate the last day of the month.
    ** Adjust backwards for holidays.
    */

    static char routine[] = "JpmcdsHolidayListDateToBusinessEOM";
    int         status    = FAILURE;

    if (JpmcdsDateToEOM (inDate, outDate) != SUCCESS)
        goto done;

    if (JpmcdsHolidayListBusinessDay (*outDate,
                                   JPMCDS_BAD_DAY_PREVIOUS,
                                   hl,
                                   outDate) != SUCCESS)
        goto done;

    status = SUCCESS;

done:

    if (status != SUCCESS)
        JpmcdsErrMsgFailure (routine);

    return status;
}
