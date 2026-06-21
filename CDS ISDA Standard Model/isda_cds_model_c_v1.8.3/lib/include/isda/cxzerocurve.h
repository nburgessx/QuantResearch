/*
 * ISDA CDS Standard Model
 *
 * Copyright (C) 2009 International Swaps and Derivatives Association, Inc.
 * Developed and supported in collaboration with Markit
 * 
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the ISDA CDS Standard Model Public License.
 */

#ifndef CX_ZEROCURVE_H
#define CX_ZEROCURVE_H

#include "tcurve.h"

#ifdef __cplusplus
extern "C"
{
#endif


/*f
***************************************************************************
** Calculates the zero price for a given date. Returns NaN for errors.
***************************************************************************
*/
EXPORT double JpmcdsZeroPrice
(TCurve* zeroCurve,
 TDate   date);


/*f
***************************************************************************
** Calculates the zero price for a given start date and maturity date.
** Returns NaN for errors.
***************************************************************************
*/
double JpmcdsForwardZeroPrice
(TCurve* zeroCurve,
 TDate   startDate,
 TDate   maturityDate);


/*f
***************************************************************************
** Calculates the zero rate for a given date using ACT/365F and continously
** compounded rates.
***************************************************************************
*/
double JpmcdsZeroRate
(TCurve* zeroCurve,
 TDate   date);


/*f
***************************************************************************
** Converts a compound rate from one frequency to another.
** Can also convert between ACT-style day count conventions.
***************************************************************************
*/
int JpmcdsConvertCompoundRate 
(double  inRate,
 double  inBasis,
 long    inDayCountConv,
 double  outBasis,
 long    outDayCountConv,
 double *outRate);


#ifdef __cplusplus
}
#endif

#endif /* CX_ZERO_CURVE_D_H */

