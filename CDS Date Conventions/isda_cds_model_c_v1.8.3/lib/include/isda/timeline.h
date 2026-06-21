/*
 * ISDA CDS Standard Model
 *
 * Copyright (C) 2009 International Swaps and Derivatives Association, Inc.
 * Developed and supported in collaboration with Markit
 * 
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the ISDA CDS Standard Model Public License.
 */

#ifndef TIMELINE_H
#define TIMELINE_H

#include "cx.h"

#ifdef __cplusplus
extern "C"
{
#endif

/*f
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
 TCurve           *discCurve,
 TCurve           *riskyCurve);


/*f
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
 TDate      endDate);

#ifdef __cplusplus
}
#endif

#endif



