/*
 * ISDA CDS Standard Model
 *
 * Copyright (C) 2009 International Swaps and Derivatives Association, Inc.
 * Developed and supported in collaboration with Markit
 * 
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the ISDA CDS Standard Model Public License.
 */

#ifndef LINTRP_H
#define LINTRP_H

#include "cgeneral.h"
#include "metric.h"


#ifdef __cplusplus
extern "C"
{
#endif


int JpmcdsLinInterpLongPoint1 (
    long             *x,                /* (I) Ordered Array of X values */
    int               xskip,            /* (I) # bytes between x values */
    int               N,                /* (I) Length of X & F arrays */
    double           *f,                /* (I) Ordered Array of F values */
    int               fskip,            /* (I) # bytes between f values */
    double            xDesired,         /* (I) X for which F is desired */
    TMetricDoubleFunc mfunc,            /* (I) Metric Function */
    double           *fInterp);         /* (O) Interpolated value */

int JpmcdsLinInterpLongArray1 (
    long             *x,                /* (I) Ordered Array of X values */
    int               xskip,            /* (I) # bytes between x values */
    int               N,                /* (I) Length of X & F arrays */
    double           *f,                /* (I) Ordered Array of F values */
    int               fskip,            /* (I) # bytes between f values */
    long             *xDesired,         /* (I) Unordered array desired x */
    int               dskip,            /* (I) # bytes between xDesired vals */
    int               numInterps,       /* (I) Length of xDesired & yInterp */
    TMetricDoubleFunc mfunc,            /* (I) Metric Function */
    int               oskip,            /* (I) # bytes between fInterp vals */
    double           *fInterp);         /* (O) Interpolated F values */


#ifdef __cplusplus
}
#endif

#endif    /* LINTRP_H */
