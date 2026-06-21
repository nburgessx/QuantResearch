/*
 * ISDA CDS Standard Model
 *
 * Copyright (C) 2009 International Swaps and Derivatives Association, Inc.
 * Developed and supported in collaboration with Markit
 * 
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the ISDA CDS Standard Model Public License.
 */

#include "cxbsearch.h"
#include <assert.h>
#include "macros.h"
#include "macros.h"
#include "cerror.h"


/*
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
                       -1 if no such value exists. */
 long   *hiBound)   /* (O) Index of first value strictly greater than xDesired,
                       -1 if no such value exists. */
{
    static char routine[] = "JpmcdsBinarySearchLong";
    int         status    = FAILURE;

    int count;                 /* Used to check # search steps */
    int lo;                    /* Index of low estimate */
    int hi;                    /* Index of high estimate */
    int mid = 0;               /* Index of best estimate */
    char *xp = (char *)xArray; /* Ptr to x array */

    REQUIRE (arraySize > 0);
    REQUIRE (skip >= sizeof(long));
    REQUIRE (exact != NULL);

#undef ELEMENT
#define ELEMENT(idx) *(long*)(xp + skip*(idx))

    /* If we are not in the range then we are done */
    if (xDesired < ELEMENT(0))
    {
        *exact   = -1;
        *loBound = -1;
        *hiBound = 0;
        return SUCCESS;
    }
    else if (xDesired > ELEMENT(arraySize-1))
    {
        *exact   = -1;
        *loBound = arraySize-1;
        *hiBound = arraySize;
        return SUCCESS;
    }

    /* arraySize of 1 we are done */
    if (arraySize == 1)
    {
        assert (xDesired == ELEMENT(0));
        *exact   = 0;
        *loBound = -1;
        *hiBound = arraySize;
        return SUCCESS;
    }
    
    lo = 0;
    hi = arraySize - 2;
    
    /* Do binary search to find pair of x's which surround the desired
     * X value.
     */
    for (count = arraySize+1; count > 0; --count)
    {
        mid = (hi + lo) / 2;
        
        if (xDesired < ELEMENT(mid))
            hi = mid - 1;
        else if (xDesired > ELEMENT(mid + 1))
            lo = mid + 1;
        else
            break;                  /* Done */
    }
    
    if (count == 0)
    {
        JpmcdsErrMsg("%s: x array not in increasing order.n", routine);
        return FAILURE;
    }
    
    /* Protect against a run of x values which are the same.
     * Set two surrounding indices to be lo and hi.
     * Note that there is no danger of running off the end
     * since the only way for x[lo] = x[hi] is for both
     * to be equal to xDesired. But from check at beginning,
     * we know X[N-1] <> xDesired. 
     */

    assert (mid < arraySize);
    assert (xDesired >= ELEMENT(mid));
    assert (xDesired <= ELEMENT(mid+1));

    lo = mid;
    hi = mid+1;
    
    if (ELEMENT(lo) == xDesired)
        *exact = lo;
    else if (ELEMENT(hi) == xDesired)
        *exact = hi;
    else
        *exact = -1;

    if (loBound != NULL)
    {
        while (lo >= 0 && ELEMENT(lo) >= xDesired)
            --lo;
        
        if (lo >= 0)
            *loBound = lo;
        else
            *loBound = -1;
    }

    if (hiBound != NULL)
    {
        while (hi < arraySize && ELEMENT(hi) <= xDesired)
            ++hi;
        
        if (hi < arraySize)
            *hiBound = hi;
        else
            *hiBound = arraySize;
    }
    status = SUCCESS;

 done:
    
    if (status != SUCCESS) JpmcdsErrMsgFailure(routine);
    return status;
}

