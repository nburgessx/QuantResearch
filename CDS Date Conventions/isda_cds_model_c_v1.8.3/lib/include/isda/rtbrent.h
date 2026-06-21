/*
 * ISDA CDS Standard Model
 *
 * Copyright (C) 2009 International Swaps and Derivatives Association, Inc.
 * Developed and supported in collaboration with Markit
 * 
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the ISDA CDS Standard Model Public License.
 */

#ifndef RTBRENT_H
#define RTBRENT_H

#include "cgeneral.h"

#ifdef __cplusplus
extern "C"
{
#endif


typedef int (*TObjectFunc) (double x, void * para, double *f);


/*f
***************************************************************************
** Finds the root of f(x) = 0 using a combination of secant, bisection and 
** an inverse quadratic interpolation method.
***************************************************************************
*/
int JpmcdsRootFindBrent(
   TObjectFunc funcd,              /* (I) function to be solved */
   void        *data,              /* (I) data to pass into funcd */
   double      boundLo,            /* (I) lower bound on legal X */
   double      boundHi,            /* (I) upper bound on legal X */
   int         numIterations,      /* (I) Maximum number of iterations */
   double      guess,              /* (I) Initial guess */
   double      initalXStep,        /* (I) Size of step in x */
   double      initialFDeriv,      /* (I) Derivative, defaults to 0 */
   double      xacc,               /* (I) X accuracy tolerance */
   double      facc,               /* (I) function accuracy tolerance */
   double      *solution);         /* (O) root found */


#ifdef __cplusplus
}
#endif

#endif    /* RTBRENT_H */
