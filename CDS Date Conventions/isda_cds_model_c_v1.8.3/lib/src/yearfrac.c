/*
 * ISDA CDS Standard Model
 *
 * Copyright (C) 2009 International Swaps and Derivatives Association, Inc.
 * Developed and supported in collaboration with Markit
 * 
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the ISDA CDS Standard Model Public License.
 */

#include "yearfrac.h"  
#include <ctype.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "cgeneral.h"
#include "ldate.h"
#include "cdate.h"
#include "cerror.h"
#include "bastypes.h"
#include "strutil.h"
#include "macros.h"

/* Centralized place for strings for day count conventions. */
#define JPMCDS_ACT_365_STR       "Actual/365"
#define JPMCDS_ACT_365F_STR      "Actual/365F"
#define JPMCDS_ACT_360_STR       "Actual/360"
#define JPMCDS_B30_360_STR       "30/360"
#define JPMCDS_B30E_360_STR      "30E/360"
#define JPMCDS_EFFECTIVE_RATE_STR "Effective"

/* The following daycount convention can be encoded and decoded but is not
 * a valid DCC for any calculations */
#define JPMCDS_DEFAULT_STR       "Default"


/*
***************************************************************************
** Metric function which returns actual # days between 2 dates.
***************************************************************************
*/
static double actualDaysDiff(double startDate, double endDate)
{
    return (endDate - startDate);
}


/*
***************************************************************************
** Metric function which returns actual # bond days between 2 dates.
***************************************************************************
*/
static double bondDaysDiff(double startDate, double endDate)
{
    long days;

    (void)JpmcdsDaysDiff((TDate)startDate, (TDate)endDate, JPMCDS_B30_360, &days);
    return ((double)days);
}


/*
***************************************************************************
** Metric function which returns actual # bond E days between 2 dates.
***************************************************************************
*/
static double bondEDaysDiff(double startDate, double endDate)
{
    long days;

    (void)JpmcdsDaysDiff((TDate)startDate, (TDate)endDate, JPMCDS_B30E_360, &days);
    return ((double)days);
}


/*
***************************************************************************
** Returns metric function to use for interpolation routines for a given 
** day count convention.
***************************************************************************
*/
TMetricDoubleFunc JpmcdsDayCountToMetricFunc(long dayCountConv)
{
    static char routine[]="JpmcdsDayCountToMetricFunc";

    switch(dayCountConv)
    {
        case JPMCDS_ACT_365F:
        case JPMCDS_ACT_365:
        case JPMCDS_ACT_360:
            return (TMetricDoubleFunc) actualDaysDiff;
            /*NOTREACHED*/
            
        case JPMCDS_B30_360:
            return (TMetricDoubleFunc) bondDaysDiff;
            /*NOTREACHED*/
            
        case JPMCDS_B30E_360:
            return (TMetricDoubleFunc) bondEDaysDiff;          
            /*NOTREACHED*/

        default:
            JpmcdsErrMsg("%s: Day count convention %ld unknown.\n", routine, dayCountConv);
            return (TMetricDoubleFunc) NULL;
    }
}


/*
***************************************************************************
** Maps day count conventions to (private) strings.
***************************************************************************
*/
char* JpmcdsFormatDayCountConv(long dayCountConv)
{
    static char routine[]="JpmcdsFormatDayCountConv";

    switch(dayCountConv)
    {
    case JPMCDS_DEFAULT:
        return JPMCDS_DEFAULT_STR;
        /*NOTREACHED*/

    case JPMCDS_ACT_365F:
        return JPMCDS_ACT_365F_STR;
        /*NOTREACHED*/

    case JPMCDS_ACT_365:
        return  JPMCDS_ACT_365_STR;
        /*NOTREACHED*/

    case JPMCDS_ACT_360:
        return JPMCDS_ACT_360_STR;
        /*NOTREACHED*/

    case JPMCDS_B30_360:
        return JPMCDS_B30_360_STR;
        /*NOTREACHED*/

    case JPMCDS_B30E_360:
        return JPMCDS_B30E_360_STR;
        /*NOTREACHED*/
 
    case JPMCDS_EFFECTIVE_RATE:
        return JPMCDS_EFFECTIVE_RATE_STR;
        /*NOTREACHED*/

    default:
        goto done;
    }

 done:
    JpmcdsErrMsg("%s: Day count convention %ld unknown.\n", routine, dayCountConv);
    return "unrecognized";
    /*NOTREACHED*/
}


/*
***************************************************************************
** Converts String to Day Count Convention.
**
** This function is fairly permissive in terms of what it will accept.
** For example, "Actual/360", "A", etc.
** Strings of the form "30/Actual" are disallowed, however.
***************************************************************************
*/
int JpmcdsStringToDayCountConv(char *dayCountString, long *type)
{
    static char  routine[]="JpmcdsStringToDayCountConv";
#define MAX_DCC_CHARS 32
    static char  privDccString[MAX_DCC_CHARS + 1];
    char *str = NULL;

    /* Remove white space, convert to capitals.
    */
    if( JpmcdsStringPreprocess(dayCountString,
                            MAX_DCC_CHARS,
                            privDccString) == FAILURE )
    {
        goto error;
    }

    
    if (strstr(privDccString, "A"))           /* It's an Actual type. */
    {
        if (strstr(privDccString, "/A"))      /* Must come before "365". */
        {
            *type = JPMCDS_ACT_365;              /* Act/Act */
            return SUCCESS;
        }
        if (strstr(privDccString, "365"))
        {  
            if (strstr(privDccString, "F"))
            {  
                *type = JPMCDS_ACT_365F;         /* Act/365F */
                return SUCCESS;
            }
            *type = JPMCDS_ACT_365;              /* Act/365 */
            return SUCCESS;
        }
        if (strstr(privDccString, "360"))
        {  
            *type = JPMCDS_ACT_360;
            return SUCCESS;
        }
    }
    /* The following two main cases are disjoint, e.g., "30/Actual" is illegal. */
    else if ((strstr(privDccString, "30") ||    /* It's a 30/Bond type */
             strstr(privDccString, "B"))  &&
             strstr(privDccString, "360"))
    {
        if (strstr(privDccString, "E"))
        {
            *type = JPMCDS_B30E_360;
            return SUCCESS;
        }
        *type = JPMCDS_B30_360;
        return SUCCESS;
    }
    else if (strstr(privDccString, "EFF"))
    {
        /* Effective rate */
        *type = JPMCDS_EFFECTIVE_RATE;
        return SUCCESS;
    }
    
    if (strcmp(privDccString, "DEFAULT") == 0)
    {
        /* DCC which can be encoded/decoded but not usable for calculations */
        *type = JPMCDS_DEFAULT;
        return SUCCESS;
    }

    if (strcmp(privDccString, "NONE") == 0)
    {
        /* DCC which can be encoded/decoded but not usable for calculations */
        *type = 0;
        return SUCCESS;
    }

    JpmcdsErrMsg("%s: Unrecognized day count convention \"%s\"\n", routine, dayCountString);
    
 error:
    FREE(str);
    JpmcdsErrMsg("%s: Failed.\n", routine);
    return FAILURE;
}


/*
***************************************************************************
** Returns SUCCESS if day count convention is valid, FAILURE otherwise.
**
** If routine != NULL, then populates the error log on FAILURE. Otherwise
** if routine == NULL, then the error log is not populated on FAILURE.
***************************************************************************
*/
int JpmcdsDayCountValid(char *routine, long dayCountConv)
{
    switch(dayCountConv)
    {
    case JPMCDS_ACT_365F:
    case JPMCDS_ACT_365:
    case JPMCDS_ACT_360:
    case JPMCDS_B30_360:
    case JPMCDS_B30E_360:
    case JPMCDS_EFFECTIVE_RATE:
        return SUCCESS ;
        /*NOTREACHED*/

    default:
        return FAILURE;
    }
}

