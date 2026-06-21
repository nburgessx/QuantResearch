/*
 * ISDA CDS Standard Model
 *
 * Copyright (C) 2009 International Swaps and Derivatives Association, Inc.
 * Developed and supported in collaboration with Markit
 * 
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the ISDA CDS Standard Model Public License.
 */

#ifndef FLTRATE_H
#define FLTRATE_H

#include "cgeneral.h"
#include "bastypes.h"
#include "interp.h"

#ifdef __cplusplus
extern "C"
{
#endif


/*f
***************************************************************************
** Sets a TFloatRate. 
***************************************************************************
*/
int JpmcdsFloatRateSet
    (TDateInterval *matInterval,        /* (I) */
     TDateInterval *payInterval,        /* (I) */
     long           dayCountConv,       /* (I) */
     long           spotOffsetDays,     /* (I) */
     double         spread,             /* (I) */
     double         weight,             /* (I) */
     TFloatRate    *floatRate);         /* (O) Floating rate definition */


#ifdef __cplusplus
}
#endif

#endif    /* FLTRATE_H */

