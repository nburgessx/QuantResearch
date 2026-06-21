/*
 * ISDA CDS Standard Model
 *
 * Copyright (C) 2009 International Swaps and Derivatives Association, Inc.
 * Developed and supported in collaboration with Markit
 * 
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the ISDA CDS Standard Model Public License.
 */

#ifndef INTERP_H
#define INTERP_H

#include "bastypes.h"

#ifdef __cplusplus
extern "C"
{
#endif

/*t
 */
typedef TCurve * (TCurveInterpFunc)(TDateList *desiredDates,
                                    TCurve    *oldCurve,
                                    void      *interpData);


/*f
***************************************************************************
** Interpolates a TCurve at a TDate using the specified interp method. 
***************************************************************************
*/
int JpmcdsInterpRate(
    TDate   desiredDate,                /* (I) Date to interpolate at */
    TCurve *curve,                      /* (I) Curve with rates */
    long    interpType,                 /* (I) Type of interpolation to use */
    double *rate);                      /* (O) Interpolated rate */


/*f
***************************************************************************
** Performs interpolation at "exact" or potentially fractional dates.
***************************************************************************
*/
int JpmcdsInterpRateExact( 
   double  desiredDate,                /* (I) Date to interpolate at */
   TCurve *curve,                      /* (I) Curve with rates */
   long    interpType,                 /* (I) Type of interpolation to use */
   double *rate);                      /* (O) Interpolated rate */


#define JPMCDS_LINEAR_INTERP          0L

/* Interpolation methods which can be applied to a zero curve only. */
#define JPMCDS_LINEAR_FORWARDS      123L     /* Linear forward interpolation  */
#define JPMCDS_FLAT_FORWARDS        124L     /* Flat forward interpolation  */


/*f 
***************************************************************************
** Interpolates a TCurve at a set of dates using linear interp, and produces
** a new TCurve.
***************************************************************************
*/
TCurveInterpFunc JpmcdsLinInterpCurve;


#ifdef __cplusplus
}
#endif

#endif    /* INTERP_H */




