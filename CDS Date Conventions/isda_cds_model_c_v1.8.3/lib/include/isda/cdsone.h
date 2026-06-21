/*
 * ISDA CDS Standard Model
 *
 * Copyright (C) 2009 International Swaps and Derivatives Association, Inc.
 * Developed and supported in collaboration with Markit
 * 
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the ISDA CDS Standard Model Public License.
 */

#ifndef CDSONE_H
#define CDSONE_H

#include "cx.h"
#include "stub.h"

#ifdef __cplusplus
extern "C"
{
#endif

/*f
***************************************************************************
** Computes the upfront charge for a flat spread par curve.
***************************************************************************
*/
EXPORT int JpmcdsCdsoneUpfrontCharge
(TDate           today,
 TDate           valueDate,
 TDate           benchmarkStartDate,  /* start date of benchmark CDS for
                                      ** internal clean spread bootstrapping */
 TDate           stepinDate,
 TDate           startDate,           /* CDS start date, can be in the past */
 TDate           endDate,
 double          couponRate,
 TBoolean        payAccruedOnDefault,
 TDateInterval  *dateInterval,
 TStubMethod    *stubType,
 long            accrueDCC,
 long            badDayConv,
 char           *calendar,
 TCurve         *discCurve,
 double          oneSpread,
 double          recoveryRate,
 TBoolean        payAccruedAtStart,
 double         *upfrontCharge);


/*f
***************************************************************************
** Computes the flat spread required to match the upfront charge.
***************************************************************************
*/
EXPORT int JpmcdsCdsoneSpread
(TDate           today,
 TDate           valueDate,
 TDate           benchmarkStartDate,  /* start date of benchmark CDS for
                                      ** internal clean spread bootstrapping */
 TDate           stepinDate,
 TDate           startDate,           /* CDS start date, can be in the past */
 TDate           endDate,
 double          couponRate,
 TBoolean        payAccruedOnDefault,
 TDateInterval  *dateInterval,
 TStubMethod    *stubType,
 long            accrueDCC,
 long            badDayConv,
 char           *calendar,
 TCurve         *discCurve,
 double          upfrontCharge,
 double          recoveryRate,
 TBoolean        payAccruedAtStart,
 double         *oneSpread);

#ifdef __cplusplus
}
#endif

#endif

