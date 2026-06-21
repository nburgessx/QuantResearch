/*
 * ISDA CDS Standard Model
 *
 * Copyright (C) 2009 International Swaps and Derivatives Association, Inc.
 * Developed and supported in collaboration with Markit
 * 
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the ISDA CDS Standard Model Public License.
 */

#include "bastypes.h"
#include "cgeneral.h"
#include "cmemory.h"
#include "cerror.h"
#include "convert.h"                                    
#include "date_sup.h"
#include "ldate.h"
#include "stub.h"
#include "yearfrac.h"
#include "zr2coup.h"
#include "zr2fwd.h"
#include "fltrate.h"
#include <math.h>


/*
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
     TFloatRate    *floatRate)          /* (O) Floating rate definition */
{
    floatRate->matInterval = *matInterval;
    floatRate->payInterval = *payInterval;
    floatRate->dayCountConv = dayCountConv;
    floatRate->spread = spread;
    floatRate->weight = weight;
    JPMCDS_SET_ADJ_INTERVAL_DAYS(floatRate->spotOffset,spotOffsetDays)
    return SUCCESS;
}
