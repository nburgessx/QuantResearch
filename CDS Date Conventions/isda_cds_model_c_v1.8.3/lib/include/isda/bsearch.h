/*
 * ISDA CDS Standard Model
 *
 * Copyright (C) 2009 International Swaps and Derivatives Association, Inc.
 * Developed and supported in collaboration with Markit
 * 
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the ISDA CDS Standard Model Public License.
 */
 
#ifndef BSEARCH_H
#define BSEARCH_H

#include "cgeneral.h"

#ifdef __cplusplus
extern "C"
{
#endif


    /*f
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
     int *hi_idx);                      /* (O) High index for X array */

/*f
***************************************************************************
** Binary search through an array. This can be an array of doubles (set
** skip=sizeof(double)) or an array of some other structure which includes
** a double as one of its fields (skip=sizeof(structure)).
**
** Assumes that the array is in sequence.
***************************************************************************
*/
int JpmcdsBSearchDoubleFast
    (double xDesired,                   /* (I) X for which Y is desired */
     double *x,                         /* (I) Ordered Array of X values */
     int skip,                          /* (I) # bytes between X values */
                                        /* (=sizeof(X_TYPE) if x is an array)*/
     int N,                             /* (I) # elements of X array */
     int *lo_idx,                       /* (O) Low index for X array */
     int *hi_idx);                      /* (O) High index for X array */


/*f
***************************************************************************
** Validates that an array is in order. This can be an array of doubles
** (set skip=sizeof(double)) or an array of some other structure which
** includes a double as one of its fields (skip=sizeof(structure)).
***************************************************************************
*/
int JpmcdsCheckDoubleArrayOrder
(char *routine,                     /* (I) Routine name for error messages */
 double *x,                         /* (I) Ordered Array of X values */
 int skip,                          /* (I) # bytes between X values */
                                    /* (=sizeof(X_TYPE) if x is an array)*/
 int N);                            /* (I) # elements of X array */


/*f
***************************************************************************
** Binary search through an array. This can be an array of longs (set
** skip=sizeof(long)) or an array of some other structure which includes
** a long as one of its fields (skip=sizeof(structure)).
**
** Normally checks that the array is in sequence. This slows the function
** down to a linear search + binary search. You can either call
** JpmcdsBSearchCheckOrder to turn this checking off (and back on again), or
** else call JpmcdsBSearchDoubleFast.
***************************************************************************
*/
int JpmcdsBSearchLong
    (
     double xDesired,                   /* (I) X for which Y is desired */
     long *x,                           /* (I) Ordered Array of X values */
     int skip,                          /* (I) # bytes between X values */
                                        /* (=sizeof(X_TYPE) if x is an array)*/
     int N,                             /* (I) # elements of X array */
     int *lo_idx,                       /* (O) Low index for X array */
     int *hi_idx);                      /* (O) High index for X array */


/*f
***************************************************************************
** Binary search through an array. This can be an array of longs (set
** skip=sizeof(long)) or an array of some other structure which includes
** a long as one of its fields (skip=sizeof(structure)).
**
** Assumes that the array is in sequence.
***************************************************************************
*/
int JpmcdsBSearchLongFast
    (
     double xDesired,                   /* (I) X for which Y is desired */
     long *x,                           /* (I) Ordered Array of X values */
     int skip,                          /* (I) # bytes between X values */
                                        /* (=sizeof(X_TYPE) if x is an array)*/
     int N,                             /* (I) # elements of X array */
     int *lo_idx,                       /* (O) Low index for X array */
     int *hi_idx);                      /* (O) High index for X array */


/*f
***************************************************************************
** Validates that an array is in order. This can be an array of longs
** (set skip=sizeof(long)) or an array of some other structure which
** includes a long as one of its fields (skip=sizeof(structure)).
***************************************************************************
*/
int JpmcdsCheckLongArrayOrder
(char *routine,                   /* (I) Routine name for error messages */
 long *x,                         /* (I) Ordered Array of X values */
 int skip,                        /* (I) # bytes between X values */
                                  /* (=sizeof(X_TYPE) if x is an array)*/
 int N);                          /* (I) # elements of X array */


/* Returns signed difference beween p1 and p2. (p1-p2). */
typedef double (*TCompareFunc)(const void *p1, const void *p2);


/*f
***************************************************************************
** Sets flag indicating whether to check that input Array of X vals is in order
** The default state is TRUE (check IS performed).
** JpmcdsBSearchCheckOrder() returns state before it is switched to 'doCheckOrder'
** NB This function should be used with caution since it affects behaviour of
** other higher-level functions which may be relying upon the check in 
** JpmcdsBSearchLong, JpmcdsBSearchDouble (according to usage).
***************************************************************************
*/
TBoolean JpmcdsBSearchCheckOrder(
    TBoolean        doCheckOrder);     /* (I) TRUE  => Default Yes check */
                                       /*     FALSE => Assume in order */



#ifdef __cplusplus
}
#endif

#endif    /* BSEARCH_H */
