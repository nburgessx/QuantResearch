/*
 * ISDA CDS Standard Model
 *
 * Copyright (C) 2009 International Swaps and Derivatives Association, Inc.
 * Developed and supported in collaboration with Markit
 * 
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the ISDA CDS Standard Model Public License.
 */

#include "timeline.h"
#include "cxdatelist.h"
#include "macros.h"
#include "cerror.h"
#include "datelist.h"
#include "tcurve.h"


/*
***************************************************************************
** Returns a timeline for use with risky integrations assuming flat
** forward curves.
**
** The timeline will contain the following dates:
** - startDate
** - endDate
** - all the points in the discount curve
** - all the points in the risky curve
** - nothing before startDate
** - nothing after endDate
***************************************************************************
*/
TDateList* JpmcdsRiskyTimeLine
(TDate             startDate,
 TDate             endDate,
 TCurve*           discCurve,
 TCurve           *spreadCurve)
{
    static char routine[] = "JpmcdsRiskyTimeLine";

    TDateList *tl  = NULL;
    TDate     *dates = NULL;

    REQUIRE (discCurve != NULL);
    REQUIRE (spreadCurve != NULL);
    REQUIRE (endDate > startDate);
    REQUIRE (spreadCurve != NULL);

    /*
    ** Timeline is points on the spreadCurve between startDate and endDate,
    ** plus the startDate and endDate, plus the critical dates.
    */
    tl = JpmcdsNewDateListFromTCurve(discCurve);
    if (tl == NULL) goto done;

    dates = JpmcdsDatesFromCurve(spreadCurve);
    tl = JpmcdsDateListAddDatesFreeOld(tl, spreadCurve->fNumItems, dates);
    if (tl == NULL) goto done;
    tl = JpmcdsDateListAddDatesFreeOld(tl, 1, &startDate);
    if (tl == NULL) goto done;
    tl = JpmcdsDateListAddDatesFreeOld(tl, 1, &endDate);
    if (tl == NULL) goto done;

    /* remove dates strictly before startDate and strictly after endDate */
    tl = JpmcdsDateListTruncate(tl, startDate, TRUE, TRUE, TRUE);
    tl = JpmcdsDateListTruncate(tl, endDate, TRUE, FALSE, TRUE);
    
 done:

    if (tl == NULL)
        JpmcdsErrMsgFailure(routine);

    FREE (dates);

    return tl;
}


/*
***************************************************************************
** Truncate timeline.
**
** Truncates a timeline so that it will contain the following:
** - startDate
** - endDate
** - criticalDates (can be NULL)
** - nothing before startDate
** - nothing after endDate
***************************************************************************
*/
TDateList* JpmcdsTruncateTimeLine
(TDateList* criticalDates,
 TDate      startDate,
 TDate      endDate)
{
    static char routine[] = "JpmcdsTruncateTimeLine";
    
    TDateList *tl = NULL;
    TDate      startEndDate[2];

    REQUIRE (endDate > startDate);

    startEndDate[0] = startDate;
    startEndDate[1] = endDate;

    tl = JpmcdsDateListAddDates (criticalDates, 2, startEndDate);
    if (tl == NULL)
        goto done;

    /* remove dates strictly before startDate and strictly after endDate */
    tl = JpmcdsDateListTruncate (tl, startDate, TRUE, TRUE, TRUE);
    tl = JpmcdsDateListTruncate (tl, endDate, TRUE, FALSE, TRUE);
    
 done:

    if (tl == NULL)
        JpmcdsErrMsgFailure (routine);
    
    return tl;
}
