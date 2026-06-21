/*
 * ISDA CDS Standard Model
 *
 * Copyright (C) 2009 International Swaps and Derivatives Association, Inc.
 * Developed and supported in collaboration with Markit
 * 
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the ISDA CDS Standard Model Public License.
 */

#ifndef CX_BSEARCH_H
#define CX_BSEARCH_H

#include <stdlib.h>

#ifdef __cplusplus
extern "C"
{
#endif


/*f
***************************************************************************
** Searches for a long value within an array.
**
** The array can either be an array of longs, or else an array of something
** bigger than a long which has one of its elements equal to the value that
** we seek.
**
** Returns three values - the exact value (if it exists), the lower bound
** and the upper bound.
**
** The lower bound is defined as the index of the last value which is
** strictly less than the desired value.
**
** The upper bound is defined as the index of the first value which is
** strictly greater than the desired value.
**
** The exact value is the index of a value which is equal to the desired
** value.
**
** In the case that such an index cannot be found, then the value of -1
** is returned for exact and loBound, and arraySize for hiBound.
**
** The idea behind these return values is that hiBound-loBound-1 should
** give the number of exact matches in all cases, and that everything in
** the range [loBound+1, hiBound-1] inclusive should be an exact match.
**
** If the input array is not sorted in ascending order, then the output
** is undefined. No checks are performed on the order of xArray - this
** would in fact defeat the purpose of having a fast binary search if you
** are going to slow it down by doing a linear check of the inputs.
***************************************************************************
*/
int JpmcdsBinarySearchLong
(long    xDesired,  /* (I) Value for which we seek */
 long   *xArray,    /* (I) Array in which we seek - assumes sorted in
                       ascending order.*/
 size_t  skip,      /* (I) Size of elements in the array */
 long    arraySize, /* (I) Size of the array */
 long   *exact,     /* (O) Index of a value which is equal to xDesired,
                       -1 if not found */
 long   *loBound,   /* (O) Index of last value strictly less than xDesired,
                       -1 if no such value exists
                       - can be NULL if you don't care. */
 long   *hiBound);  /* (O) Index of first value strictly greater than xDesired,
                       arraySize if no such value exists. 
                       - can be NULL if you don't care */


#ifdef __cplusplus
}
#endif

#endif
