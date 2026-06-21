/*
 * ISDA CDS Standard Model
 *
 * Copyright (C) 2009 International Swaps and Derivatives Association, Inc.
 * Developed and supported in collaboration with Markit
 * 
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the ISDA CDS Standard Model Public License.
 */

#ifndef MACROS_H
#define MACROS_H

#include "cmemory.h"
#include <string.h>
#include <float.h>


#ifdef __cplusplus
extern "C"
{
#endif

#ifndef MAX
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#endif

#ifndef MIN
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#endif

#ifndef ABS
#define ABS(x) ((x) > 0 ? (x) : -(x))
#endif

/* This should be correct on every platform that supports the IEEE
 * floating point standard, i.e., everything except DEC VAXes.
 *
 * Note that it is only defined here if it was _not_ defined in float.h
 *
 */
#ifndef DBL_EPSILON
#define DBL_EPSILON     2.2204460492503131E-16
#endif


/* REQUIRE should be used to check input conditions. */
/* REQUIRE assumes existence of 'done' label for error cases. */
/* REQUIRE assumes that the function is named 'routine'. */
#undef REQUIRE
#define REQUIRE(cond) do { if (!(cond))\
{\
    JpmcdsErrMsg("%s: Required condition (%s) fails!\n",routine,#cond);\
    goto done;\
}} while(0)


/* ASSERT should be used to check things should be true from prior code      */
/* ASSERT should not be used to check input conditions - use REQUIRE instead */
/* ASSERT assumes existence of 'done' label for error cases.                 */
/* ASSERT assumes that the function is named 'routine'.                      */
/* ASSERT has same syntax as the system macro assert, but failures are not   */
/*        quite so drastic as assert().                                      */
/* ASSERT is checked in release code, unlike assert() which is only checked  */
/*        in debug code.                                                     */
#undef ASSERT
#define ASSERT(cond) do { if (!(cond))\
{\
    JpmcdsErrMsg("%s: Assertion (%s) fails: %s line %d\n",routine,#cond,__FILE__,__LINE__);\
    goto done;\
}} while(0)

#ifndef PROGRAM_BUG
#define PROGRAM_BUG()   JpmcdsErrMsg("Program bug:%s line %d\n",__FILE__,__LINE__)
#endif

#ifndef JPMCDS_HUGE_RATE
#define JPMCDS_HUGE_RATE  1E10             /* Used for root finding  */
#endif

#ifndef JPMCDS_MAX_RATE
#define JPMCDS_MAX_RATE    100.            /* 10,000%. Used for validation */
#endif
    
#ifndef JPMCDS_LOG0_THRESHOLD
#define JPMCDS_LOG0_THRESHOLD      1E-100  /* Used for avoiding LOG0 error  */
#endif


/* IS_ALMOST_ZERO returns TRUE if and only if 1.0 + x = 1.0, assuming x is a
 * double. In other words, it compares x with 1.0
 */
#define IS_ALMOST_ZERO(x) ((x) < DBL_EPSILON && (x) > -DBL_EPSILON)
#define ARE_ALMOST_EQUAL(x,y) IS_ALMOST_ZERO((x)-(y))
#define IS_BETWEEN(x,a,b) ((a) < (b) ? \
                           ((x) >= (a) && (x) <= (b)) : \
                           ((x) >= (b) && (x) <= (a)))

/*
** With gcc we use -Wfloat-equal.
**
** In some cases though we really want to allow exact float comparisons.
** In these cases use these macros to avoid the warning messages.
*/
#ifdef __GNUC__

#define IS_NOT_EQUAL(x,y) ((x)<(y) || (x)>(y))
#define IS_EQUAL(x,y) !IS_NOT_EQUAL((x),(y))

#else

#define IS_NOT_EQUAL(x,y) ((x)!=(y))
#define IS_EQUAL(x,y) ((x)==(y))

#endif


#ifndef NEW
#define NEW(t)                  (t *) JpmcdsMallocSafe(sizeof(t))
#define NEW_ARRAY(t,n)          (t *) JpmcdsMallocSafe(sizeof(t)*(n))
#define FREE_ARRAY(ptr)         FREE(ptr)
#endif

#ifndef COPY_ARRAY
#define COPY_ARRAY(dst,src,t,n) memcpy((char*)(dst),(char*)(src),(n)*sizeof(t))
#endif

#ifndef FREE
#define FREE(x) JpmcdsFreeSafe ((void *) (x))
#endif


#ifdef __cplusplus
}
#endif

#endif    /* MACROS_H */
