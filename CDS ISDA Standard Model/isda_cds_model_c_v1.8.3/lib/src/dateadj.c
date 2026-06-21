/*
 * ISDA CDS Standard Model
 *
 * Copyright (C) 2009 International Swaps and Derivatives Association, Inc.
 * Developed and supported in collaboration with Markit
 * 
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the ISDA CDS Standard Model Public License.
 */

#include "dateadj.h"
#include <ctype.h>
#include <string.h>
#include "convert.h"
#include "cerror.h"
#include "ldate.h"
#include "macros.h"
#include "busday.h"
#include "strutil.h"


/*f
***************************************************************************
** Calculates a date being at a business or calendar days offset from the 
** input date.
***************************************************************************
*/
int JpmcdsDtFwdAdj
   (TDate          startDate,    /* (I) input date       */
    TDateAdjIntvl *adjIvl,       /* (I) in business or calendar days  */
    TDate         *result)       /* (O) resulting date   */
{
    int         status = FAILURE;
    static char routine[]="JpmcdsDtFwdAdj";

    switch (adjIvl->isBusDays)
    {
    case JPMCDS_DATE_ADJ_TYPE_BUSINESS:  /* offset by business days */
        if (adjIvl->interval.prd_typ == 'D')
        {
            if (JpmcdsDateFromBusDaysOffset
                (startDate,
                 adjIvl->interval.prd, /* Number days */
                 adjIvl->holidayFile,
                 result) != SUCCESS)
            {
                goto done;
            }
        }
        else
        {
            /*
            ** Non-daily interval and offset by business days.
            ** This means that if the start date is the end of the business
            ** month, then the end date must also be the end of the month.
            ** Actually we do adjust the end date by the bad day convention,
            ** but in the context of generating a set of cash flow dates,
            ** this is often not done when computing the maturity date.
            */
            TDate businessEOM;
            TDate endDate;
            
            if (JpmcdsDateToBusinessEOM (startDate,
                                      adjIvl->holidayFile,
                                      &businessEOM) != SUCCESS)
            {
                goto done;
            }

            if (JpmcdsDtFwdAny (startDate,
                             &adjIvl->interval,
                             &endDate) != SUCCESS)
            {
                goto done;
            }

            if (startDate == businessEOM)
            {
                if (JpmcdsDateToEOM (endDate, &endDate) != SUCCESS)
                    goto done;
            }

            if (JpmcdsBusinessDay (endDate, adjIvl->badDayConv,
                                adjIvl->holidayFile, result) != SUCCESS)
            {
                goto done;
            }
        }
        break;

    case JPMCDS_DATE_ADJ_TYPE_CALENDAR:
        /* offset by calendar period */
        if (JpmcdsDateFwdThenAdjust
            (startDate, &adjIvl->interval,
             adjIvl->badDayConv, adjIvl->holidayFile, result) != SUCCESS)
        {
            goto done;
        }
        break;

    case JPMCDS_DATE_ADJ_TYPE_WEEKDAY:
        /* Adjust interval as weekdays, and then bad day adjust the result. */
        if (adjIvl->interval.prd_typ == 'D')
        {
            if (JpmcdsDateFromBusDaysOffset
                (startDate,
                 adjIvl->interval.prd, /* Number days */
                 "NONE", /* A holiday file with week-ends only */
                 result) != SUCCESS)
            {
                goto done;
            }
            if (JpmcdsBusinessDay (*result, 
                                adjIvl->badDayConv, 
                                adjIvl->holidayFile,
                                result) != SUCCESS)
            {
                goto done;
            }
        }
        else
        {
            JpmcdsErrMsg ("%s: Weekday adjustment requires date interval "
                       "%s to be in days.\n",
                       routine, JpmcdsFormatDateInterval (&adjIvl->interval));
            goto done;
        }
        break;

    default:
        JpmcdsErrMsg ("%s: Invalid date adjustment type %d.\n", routine, adjIvl->isBusDays);
        goto done;
    }
   
    status = SUCCESS;

done:
    
    if (status != SUCCESS)
        JpmcdsErrMsg("%s: Failed.\n", routine);

    return status;
}  
