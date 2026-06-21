/*
 * ISDA CDS Standard Model
 *
 * Copyright (C) 2009 International Swaps and Derivatives Association, Inc.
 * Developed and supported in collaboration with Markit
 * 
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the ISDA CDS Standard Model Public License.
 */

#include "cgeneral.h"
#include "bastypes.h"
#include "interp.h"
#include "lintrp.h"
#include "tcurve.h"
#include "yearfrac.h"
#include "macros.h"
#include "cerror.h"


/*
***************************************************************************
** Using oldCurve dates & rates, routine interpolates linearly to find new 
** rates at the desiredDates; sets up a new TCurve struct, and uses the 
** interpolated rates and desired Dates to fill in the new TCurve.
**
** Returns pointer to a new TCurve on success, NULL on failure.
***************************************************************************
*/
TCurve* JpmcdsLinInterpCurve
(
    TDateList *desiredDates,            /* (I) Dates to interpolate at */
    TCurve *oldCurve,                   /* (I) Old curve */
    void *interpData                    /* (I) Interpolation data-not used */
)
{
    static char routine[] = "JpmcdsLinInterpCurve";
    int idx;
    TCurve *interpCurve = NULL;
    TDate *modDates = NULL;


    /* Make a new curve to store the new rates in.
     */
    interpCurve = JpmcdsNewTCurve(oldCurve->fBaseDate, desiredDates->fNumItems,
                               oldCurve->fBasis, oldCurve->fDayCountConv);
    if (interpCurve == (TCurve *)NULL)
        goto done;

    modDates = NEW_ARRAY(TDate, desiredDates->fNumItems);
    if (modDates == NULL)
    {
        JpmcdsFreeTCurve(interpCurve);
        interpCurve = NULL;
        goto done;
    }


    /* Make sure that we don't extrapolate. Before the first date
     * in oldCurve, we just use the first rate in oldCurve. After 
     * the last date in oldCurve, we just use the last rate in oldCurve.
     */
    for (idx = 0; idx < desiredDates->fNumItems; idx++)
    {
        modDates[idx] = MIN(desiredDates->fArray[idx],
                            oldCurve->fArray[oldCurve->fNumItems-1].fDate);
        modDates[idx] = MAX(modDates[idx],
                            oldCurve->fArray[0].fDate);
        interpCurve->fArray[idx].fDate = desiredDates->fArray[idx];
    }


    /* Interpolate */
    if (JpmcdsLinInterpLongArray1
        (&oldCurve->fArray[0].fDate, sizeof(TRatePt), oldCurve->fNumItems,
         &oldCurve->fArray[0].fRate, sizeof(TRatePt), 
         modDates, sizeof(TDate), desiredDates->fNumItems,
         JpmcdsDayCountToMetricFunc(oldCurve->fDayCountConv),
         sizeof(TRatePt), &interpCurve->fArray[0].fRate) == FAILURE)
    {
        JpmcdsFreeTCurve(interpCurve);
        interpCurve = NULL;
        goto done;
    }

  done:
    FREE_ARRAY(modDates);
    if (interpCurve == NULL)
        JpmcdsErrMsg("%s: Failed.\n", routine);

    return interpCurve;
}

