/*
 * ISDA CDS Standard Model
 *
 * Copyright (C) 2009 International Swaps and Derivatives Association, Inc.
 * Developed and supported in collaboration with Markit
 * 
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the ISDA CDS Standard Model Public License.
 */

#ifndef FEE_LEG_H
#define FEE_LEG_H

#include "cx.h"

#ifdef __cplusplus
extern "C"
{
#endif


/*f
***************************************************************************
** Calculates the PV of a fee leg with fixed fee payments.
***************************************************************************
*/
int JpmcdsFeeLegPV
(TFeeLeg       *fl,
 TDate          today,
 TDate          stepinDate,
 TDate          valueDate,
 TCurve        *discCurve,
 TCurve        *spreadCurve,
 TBoolean       payAccruedAtStart,
 double        *pv);


/*f
***************************************************************************
** Calculates the PV of the accruals which occur on default with delay.
***************************************************************************
*/
int JpmcdsAccrualOnDefaultPVWithTimeLine
(TDate           today,
 TDate           stepinDate,
 TDate           startDate,
 TDate           endDate,
 double          amount,
 TCurve         *discCurve,
 TCurve         *spreadCurve,
 TDateList      *tl,
 double         *pv);


/*f
***************************************************************************
** Computes the non-contingent cash flows for a fee leg.
***************************************************************************
*/
TCashFlowList* JpmcdsFeeLegFlows
(TFeeLeg      *fl);

#ifdef __cplusplus
}
#endif

#endif

