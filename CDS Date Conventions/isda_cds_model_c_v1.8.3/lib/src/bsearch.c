/*
 * ISDA CDS Standard Model
 *
 * Copyright (C) 2009 International Swaps and Derivatives Association, Inc.
 * Developed and supported in collaboration with Markit
 * 
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the ISDA CDS Standard Model Public License.
 */

#include <math.h>
#include "cgeneral.h"
#include "cerror.h"
#include "bsearch.h"

#define STR(x) # x
#define STRING(x) STR(x)

static TBoolean bSearchCheckOrder_g = TRUE;


/*
***************************************************************************
** Binary search through an array. This can be an array of doubles (set
** skip=sizeof(double)) or an array of some other structure which includes
** a double as one of its fields (skip=sizeof(structure)).
**
** Normally checks that the array is in sequence. This slows the function
** down to a linear search + binary search. You can either call
** JpmcdsBSearchCheckOrder to turn this checking off (and back on again), or
** else call JpmcdsBSearchDoubleFast.
***************************************************************************
*/
int JpmcdsBSearchDouble
    (double xDesired,                   /* (I) X for which Y is desired */
     double *x,                         /* (I) Ordered Array of X values */
     int skip,                          /* (I) # bytes between X values */
                                        /* (=sizeof(X_TYPE) if x is an array)*/
     int N,                             /* (I) # elements of X array */
     int *lo_idx,                       /* (O) Low index for X array */
     int *hi_idx)                       /* (O) High index for X array */
{
    static char routine[] = "JpmcdsBSearchDouble";

    if (bSearchCheckOrder_g &&
        JpmcdsCheckDoubleArrayOrder (routine, x, skip, N) != SUCCESS)
    {
        return FAILURE;
    }

    if (JpmcdsBSearchDoubleFast (xDesired, x, skip, N, lo_idx, hi_idx)
        != SUCCESS)
    {
        JpmcdsErrMsg ("%s: Failed.\n", routine);
        return FAILURE;
    }

    return SUCCESS;
}


#define SEARCH_ROUTINE_NAME JpmcdsBSearchDoubleFast
#define CHECK_ORDER_ROUTINE_NAME JpmcdsCheckDoubleArrayOrder
#define X_TYPE double
#define XD_TYPE double
#include "bsearch.inc"
#undef SEARCH_ROUTINE_NAME
#undef CHECK_ORDER_ROUTINE_NAME
#undef X_TYPE
#undef XD_TYPE


/*
***************************************************************************
** Binary search through an array. This can be an array of longs (set
** skip=sizeof(long)) or an array of some other structure which includes
** a long as one of its fields (skip=sizeof(structure)).
**
** Normally checks that the array is in sequence. This slows the function
** down to a linear search + binary search. You can either call
** JpmcdsBSearchCheckOrder to turn this checking off (and back on again), or
** else call JpmcdsBSearchLongFast.
***************************************************************************
*/
int JpmcdsBSearchLong
    (double xDesired,                   /* (I) X for which Y is desired */
     long *x,                           /* (I) Ordered Array of X values */
     int skip,                          /* (I) # bytes between X values */
                                        /* (=sizeof(X_TYPE) if x is an array)*/
     int N,                             /* (I) # elements of X array */
     int *lo_idx,                       /* (O) Low index for X array */
     int *hi_idx)                       /* (O) High index for X array */
{
    static char routine[] = "JpmcdsBSearchLong";

    if (bSearchCheckOrder_g &&
        JpmcdsCheckLongArrayOrder (routine, x, skip, N) != SUCCESS)
    {
        return FAILURE;
    }

    if (JpmcdsBSearchLongFast (xDesired, x, skip, N, lo_idx, hi_idx)
        != SUCCESS)
    {
        JpmcdsErrMsg ("%s: Failed.\n", routine);
        return FAILURE;
    }

    return SUCCESS;
}


/* Define another binary routine with the following types. */
#define SEARCH_ROUTINE_NAME JpmcdsBSearchLongFast
#define CHECK_ORDER_ROUTINE_NAME JpmcdsCheckLongArrayOrder
#define X_TYPE long
#define XD_TYPE double
#include "bsearch.inc"
#undef SEARCH_ROUTINE_NAME
#undef CHECK_ORDER_ROUTINE_NAME
#undef X_TYPE
#undef XD_TYPE

TBoolean JpmcdsBSearchCheckOrder(
    TBoolean        doCheckOrder)      /* (I) TRUE=> Default Yes check; FALSE=>Assume in order */
{
    TBoolean     checkOrderWas = bSearchCheckOrder_g;

    bSearchCheckOrder_g = doCheckOrder;

    return checkOrderWas;
}

