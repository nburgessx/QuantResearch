/*
 * ISDA CDS Standard Model
 *
 * Copyright (C) 2009 International Swaps and Derivatives Association, Inc.
 * Developed and supported in collaboration with Markit
 * 
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the ISDA CDS Standard Model Public License.
 */

#include <ctype.h>
#include <string.h>

#include "convert.h"
#include "cgeneral.h"
#include "cerror.h"
#include "macros.h"
#include "cfileio.h"
#include "cdate.h"
#include "dateconv.h"

#define CROSSOVER_YEAR 60     /* If double digit year less than this, assume year belongs to 21 century. */


/*
***************************************************************************
** Helper function to convert 3-character representation of month to a 
** long integer.
***************************************************************************
*/
static int NStringToMonth(char *cp, long *monthN);


/*
***************************************************************************
** Converts a date in string format to a TDate.
** The following formats are accepted:
**       "YYYYMMDD", "MM/DD/YYYY", "Mth-DD-YYYY", and "DD-Mth-YYYY".
***************************************************************************
*/
int JpmcdsStringToDate
    (char *mmddyyyy,            /* (I) */
     TDate *tdate)              /* (O) */
{
    char           dateCopy[32];
    TMonthDayYear  mdy;
    char          *cp;
    static char    routine[]="JpmcdsStringToDate";
    TBoolean       isMonth;

    /*
    ** use strncpy to prevent accidental use of long strings
    ** use sizeof() to allow the size of dateCopy to be increased above
    */
    (void)strncpy(dateCopy, mmddyyyy, sizeof(dateCopy)-1);         /* Don't alter original */
    dateCopy[sizeof(dateCopy)-1] = '\0';

    /* Support the international YYYYMMDD format.  This is needed
       for holiday files, regression test input and output files,
       etc. */
    if( strpbrk(dateCopy, "-/. ") == NULL )   /* Check for delimiter chars */
    {
        char smallBuf[5];

        if( strlen(dateCopy) != 8 )
            goto done;
        cp = dateCopy;

        /* Create the year. */
        smallBuf[0] = *cp++;
        smallBuf[1] = *cp++;
        smallBuf[2] = *cp++;
        smallBuf[3] = *cp++;
        smallBuf[4] = '\0';
        mdy.year = atol(smallBuf);
        if( mdy.year == 0 )
            goto done;
        /* No adjustments are needed for the year. */

        /* Create the month. */
        smallBuf[0] = *cp++;
        smallBuf[1] = *cp++;
        smallBuf[2] = '\0';
        mdy.month = atol(smallBuf);
        if( mdy.month == 0 )
            goto done;

        /* Create the day. */
        smallBuf[0] = *cp++;
        smallBuf[1] = *cp++;
        mdy.day = atol(smallBuf);
        if( mdy.day == 0 )
            goto done;
    }
    else                                   /* Parse delimited format */
    {
        /* This block parses the following formats:
              MM-DD-YYYY, Mth-DD-YYYY, DD-Mth-YYYY   */
        /* Get month (or possibly day) */
        isMonth = FALSE;
        cp = strtok(dateCopy, "-/. ");
        if (cp == NULL)
            goto done;
        if (JpmcdsSscanf(cp,"%ld",&mdy.month) == 0)  /* If not integer try string */
        {
            if (NStringToMonth(cp,&mdy.month) == FAILURE)
                goto done;

            isMonth = TRUE;    /* Make sure input doesn't contain 2 months */
        }

        /* Get day (or possibly month) */
        cp = strtok(NULL, "-/. ");
        if (cp == NULL)
            goto done;
        if (JpmcdsSscanf(cp,"%ld",&mdy.day) == 0)   /* If not integer try string */
        {
            mdy.day = mdy.month;              /* Month was really day */
            if (isMonth == TRUE || NStringToMonth(cp,&mdy.month) == FAILURE)
                goto done;
        }

        /* Get year and check */
        cp = strtok(NULL, "-/. ");
        if (cp == NULL)
            goto done;
        mdy.year = atol(cp);
        if (mdy.year > 99 && mdy.year < JPMCDS_TDATE_BASE_YEAR)
        {
            JpmcdsErrMsg("%s: Year %ld out of range.\n", routine, mdy.year);
            return FAILURE;
        }

        if (mdy.year < CROSSOVER_YEAR)  /* Assume its past year 2000 */
            mdy.year += 2000;           /* To get years since 1900 */
        else if (mdy.year < JPMCDS_TDATE_BASE_YEAR)
            mdy.year += 1900;           /* Else assume we got only 2 digits */
    }

    return JpmcdsMDYToDate(&mdy, tdate);

done:
    return FAILURE;
}


/*
***************************************************************************
** Formats a TDate for printing. The print format is YYYYMMDD. 
** Can be called eight times from the same print statement, but not more.
***************************************************************************
*/
char* JpmcdsFormatDate(TDate date) /* (I) */
{
    static int ibuf;
#define MAX_STR_LEN 16
#define MAX_AT_ONCE 8                  /* Must be a power of 2 */
    static char format[MAX_AT_ONCE][MAX_STR_LEN];
    TMonthDayYear mdy;
    ibuf = (ibuf+1)&(MAX_AT_ONCE-1); /* Toggle buffers */

    if (JpmcdsDateToMDY(date, &mdy) == FAILURE)
        sprintf(&format[ibuf][0], "%s", "bad date");
    else
    {
        if (mdy.month < 10 && mdy.day < 10)
        {
            sprintf(&format[ibuf][0], "%ld0%ld0%ld",
                    mdy.year, mdy.month, mdy.day );
        }
        else if (mdy.month < 10 && mdy.day >= 10)
        {
            sprintf(&format[ibuf][0], "%ld0%ld%ld",
                    mdy.year, mdy.month, mdy.day );
        }
        else if (mdy.month >= 10 && mdy.day < 10)
        {
            sprintf(&format[ibuf][0], "%ld%ld0%ld",
                    mdy.year, mdy.month, mdy.day );
        }
        else   /* month && day >= 10 */
        {
            sprintf(&format[ibuf][0], "%ld%ld%ld",
                    mdy.year, mdy.month, mdy.day );
        }
    }
    return &format[ibuf][0];
}


/*
***************************************************************************
** Formats a TDateInterval.
** Can be called twice from the same print statement, but not more.
***************************************************************************
*/
char* JpmcdsFormatDateInterval(TDateInterval *interval) /* (I) */
{
    static int ibuf;
#define MAX_STR_LEN 16
    static char format[2][MAX_STR_LEN];
    char periodType;
    int numPeriods;

    ibuf = !ibuf;                       /* Toggle buffers */

    if (interval == NULL)
    {
        format[ibuf][0] = 0; /* empty string */
        return &format[ibuf][0];
    }

    switch(interval->prd_typ)
    {
        case 'M':                       /* Months */
            if (interval->prd % JPMCDS_MONTHS_PER_YEAR == 0)
            {
                periodType = 'A';
                numPeriods = interval->prd/JPMCDS_MONTHS_PER_YEAR;
            }
            else if (interval->prd % JPMCDS_MONTHS_PER_SEMI == 0)
            {
                periodType = 'S';
                numPeriods = interval->prd/JPMCDS_MONTHS_PER_SEMI;
            }
            else if (interval->prd % JPMCDS_MONTHS_PER_QUARTER == 0)
            {
                periodType = 'Q';
                numPeriods = interval->prd/JPMCDS_MONTHS_PER_QUARTER;
            }
            else
            {
                periodType = 'M';
                numPeriods = interval->prd;
            }
            break;

        case 'D':                       /* Days */
            if (interval->prd % JPMCDS_DAYS_PER_WEEK == 0)
            {
                periodType = 'W';
                numPeriods = interval->prd/JPMCDS_DAYS_PER_WEEK;
            }
            else
            {
                periodType = 'D';
                numPeriods = interval->prd;
            }

        default:
            periodType = interval->prd_typ;
            numPeriods = interval->prd;
    }

    sprintf(&format[ibuf][0], "%d%c", numPeriods, periodType);
    return &format[ibuf][0];
}


/*
***************************************************************************
** Converts a string representation ("2Y", "1S", etc) to a TDateInterval.
** Number optional, default is 1, hence, e.g., "M" (= "1M") is valid.
***************************************************************************
*/
int JpmcdsStringToDateInterval
    (char *input,          /* (I) String w/ 1Y, 3M, 4D, etc */
     char *label,          /* (I) Label-for JpmcdsErrMsg only */
     TDateInterval *interval) /* (O) Value read from file */
{
    static char routine[]="JpmcdsStringToDateInterval";
    int status;
    char *inp=input;                    /* Pointer to input */
#define MAX_DIGITS 32
    char numberBuff[MAX_DIGITS];
    char *nump = numberBuff;            /* Pointer to number */
    char periodType;
    int numPeriods;

    /* Copy sign,if any, to numberBuff */
    if (*inp == '-' || *inp == '+')
    {
        *nump++ = *inp++;
    }

    /* Copy digits, if any, to numberBuff */
    while (isdigit((int)*inp))
        *nump++ = *inp++;
    *nump = '\0';                       /* Null terminate */

    if (inp != input)                 /* Found some digits */
        numPeriods = atoi(numberBuff);
    else                                /* Found on digits */
    {
        numPeriods = 1;
    }

    periodType = (char)toupper(*inp);   /* To upper case */

    status = JpmcdsMakeDateInterval(numPeriods, periodType, interval);

    if (status == FAILURE)
        JpmcdsErrMsg("%s: Failed interpreting %s.\n", routine, label);

    return status;
}


/*
***************************************************************************
** Makes a TDateInterval
***************************************************************************
*/
int JpmcdsMakeDateInterval
    (
     int numPeriods,                    /* (I) Number of periods */
     char periodType,                   /* (I) Period type */
     TDateInterval *interval)           /* (O) Value read from file */
{
    static char routine[]="JpmcdsMakeDateInterval";

    interval->flag = 0;

    switch(toupper(periodType))
    {
        case 'A':                       /* Year */
        case 'Y':                       /* Year */
            interval->prd_typ = 'M';
            interval->prd = numPeriods * JPMCDS_MONTHS_PER_YEAR;
            break;

        case 'S':                       /* Semi annual */
            interval->prd_typ = 'M';
            interval->prd = numPeriods * JPMCDS_MONTHS_PER_SEMI;
            break;

        case 'Q':                       /* Quarter */
            interval->prd_typ = 'M';
            interval->prd = numPeriods * JPMCDS_MONTHS_PER_QUARTER;
            break;

        case 'W':                       /* Week */
            interval->prd_typ = 'D';
            interval->prd = numPeriods * JPMCDS_DAYS_PER_WEEK;
            break;

        case 'U':                       /* Lunar period - 28 Days */
            interval->prd_typ = 'D';
            interval->prd = numPeriods * JPMCDS_DAYS_PER_LUNAR_MONTH;
            break;

        case 'D':                       /* Day                                 */
        case 'M':                       /* Normal Month                        */
        case 'F':                       /* Flexible End of month               */
        case 'E':                       /* End of month unconditional          */
        case 'B':                       /* Beginning of month unconditional          */
        case 'G':                       /* 29th of month                       */
        case 'H':                       /* 30th of month                       */
        case 'I':                       /* IMM date                            */
        case 'J':                       /* monthly IMM date                    */
        case 'K':                       /* Aussie quarterly IMM date           */
        case 'L':                       /* Kiwi quarterly IMM date             */
        case 'T':                       /* equity derivative expiry 3rd Friday */
        case 'V':                       /* virtual IMM dates (fortnightly)     */
            interval->prd_typ = (char) toupper(periodType);
            interval->prd = numPeriods;
            break;

        default:
            if (isprint((int)periodType))    /* If periodType is printable */
            {
                JpmcdsErrMsg("%s: Interval type %c not valid.\n", routine, periodType);
            }
            else
            {
                JpmcdsErrMsg("%s: Interval type (unprintable) not valid.\n", routine);
            }
            return FAILURE;
    }

    return SUCCESS;
}


/*
***************************************************************************
** Helper function to convert 3-character representation of month to a 
** long integer.
***************************************************************************
*/
static int NStringToMonth(char *cp, long *monthN)
{
        cp[0] = tolower(cp[0]);
        cp[1] = tolower(cp[1]);
        cp[2] = tolower(cp[2]);
        if (strcmp(cp,"jan") == 0)
        {
            *monthN = 1;
        }
        else if (strcmp(cp,"feb") == 0)
        {
            *monthN = 2;
        }
        else if (strcmp(cp,"mar") == 0)
        {
            *monthN = 3;
        }
        else if (strcmp(cp,"apr") == 0)
        {
            *monthN = 4;
        }
        else if (strcmp(cp,"may") == 0)
        {
            *monthN = 5;
        }
        else if (strcmp(cp,"jun") == 0)
        {
            *monthN = 6;
        }
        else if (strcmp(cp,"jul") == 0)
        {
            *monthN = 7;
        }
        else if (strcmp(cp,"aug") == 0)
        {
            *monthN = 8;
        }
        else if (strcmp(cp,"sep") == 0)
        {
            *monthN = 9;
        }
        else if (strcmp(cp,"oct") == 0)
        {
            *monthN = 10;
        }
        else if (strcmp(cp,"nov") == 0)
        {
            *monthN = 11;
        }
        else if (strcmp(cp,"dec") == 0)
        {
            *monthN = 12;
        }
        else
        {
            return FAILURE;
        }
        return SUCCESS;
}

