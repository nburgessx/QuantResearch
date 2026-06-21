/*
 * ISDA CDS Standard Model
 *
 * Copyright (C) 2009 International Swaps and Derivatives Association, Inc.
 * Developed and supported in collaboration with Markit
 * 
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the ISDA CDS Standard Model Public License.
 */

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
int   ARRAY_ROUTINE_NAME 
    (ARG_TYPE         *x,               /* (I) Ordered Array of X values */
     int               xskip,           /* (I) # bytes between x values */
     int               N,               /* (I) Length of X & F arrays */
     double           *f,               /* (I) Ordered Array of F values */
     int               fskip,           /* (I) # bytes between f values */
     ARG_TYPE         *xDesired,        /* (I) Unordered array desired x */
     int               dskip,           /* (I) # bytes between xDesired vals */
     int               numInterps,      /* (I) Length of xDesired & yInterp */
     TMetricDoubleFunc mfunc,           /* (I) Metric func (0=default) */
     int               oskip,           /* (I) # bytes between fInterp vals */
     double           *fInterp)         /* (O) Interpolated F values */
/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
{
    char *xdp = (char *)xDesired;       /* Pointer into xDesired array*/
    char *fp = (char *)fInterp;         /* Pointer into fInterp */
    int status = FAILURE;               /* Until proven successful */
    static char routine[]=STRING(ARRAY_ROUTINE_NAME);

    while (numInterps-- > 0)
    {
        if (POINT_ROUTINE_NAME( x, xskip, N, f, fskip, 
            (double)(*((ARG_TYPE *)xdp)), mfunc, (double *)fp) == FAILURE)
            goto done;
        xdp += dskip;
        fp += oskip;
    }
    status = SUCCESS;

done:
    if (status == FAILURE)
        JpmcdsErrMsg("%s: Failed.\n", routine);
    return status;
}

#define XARRAY(idx)    *((ARG_TYPE *)(xp + (idx)*xskip))
#define FARRAY(idx)    *((double *)(fp + (idx)*fskip))

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
int   POINT_ROUTINE_NAME
    (ARG_TYPE         *x,               /* (I) Ordered Array of X values */
     int               xskip,           /* (I) # bytes between x values */
     int               N,               /* (I) Length of X & F arrays */
     double           *f,               /* (I) Ordered Array of F values */
     int               fskip,           /* (I) # bytes between f values */
     double            xDesired,        /* (I) X for which F is desired */
     TMetricDoubleFunc mfunc,           /* (I) Metric func (0=default) */
     double           *fInterp)         /* (O) Interpolated value */
/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
{
    int lo;                             /* Low index */
    int hi;                             /* High index */
    double step_fraction;
    char *xp = (char *)x;
    char *fp = (char *)f;
    int status = FAILURE;               /* Until proven successful */
    static char routine[]=STRING(POINT_ROUTINE_NAME);

    if (SEARCH_ROUTINE_NAME(xDesired, x, xskip, N, &lo, &hi) == FAILURE)
        goto done;                      /* Failed */


    if (lo == hi)
    {                                   /* Happens if only 1 pt in curve */
        *fInterp = FARRAY(lo);         
    }
    else                                /* More than 1 pt in curve */
    {
        if (mfunc ==  NULL) /* Use default */
            step_fraction = (double)(xDesired-XARRAY(lo)) /
                            (double)(XARRAY(hi) - XARRAY(lo));
        else
            step_fraction = (*mfunc)(xDesired, (double)XARRAY(lo)) /
                            (*mfunc)((double)XARRAY(hi), (double)XARRAY(lo));
        *fInterp = FARRAY(lo) +   step_fraction * (FARRAY(hi) - FARRAY(lo));
    }
    status = SUCCESS;
done:
    if (status == FAILURE)
        JpmcdsErrMsg("%s: Failed.\n", routine);

    return status;
}





/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
int   DISCONTINUOUS_POINT_ROUTINE_NAME
    (ARG_TYPE         *x,               /* (I) Ordered Array of X values */
     int               xskip,           /* (I) # bytes between x values */
     int               N,               /* (I) Length of X & F arrays */
     double           *f,               /* (I) Ordered Array of F values */
     int               fskip,           /* (I) # bytes between f values */
     double            xDesired,        /* (I) X for which F is desired */
     TMetricDoubleFunc mfunc,           /* (I) Metric func (0=default) */
     TBoolean          rightLimit,      /* (I) T=Use right limit;F=Use left.*/
     double           *fInterp)         /* (O) Interpolated value */
/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
{
    int         lo;                             /* Low index */
    int         hi;                             /* High index */
    double      step_fraction;
    char       *xp = (char *)x;
    char       *fp = (char *)f;
    int         status = FAILURE;               /* Until proven successful */
    static char routine[]=STRING(DISCONTINUOUS_POINT_ROUTINE_NAME);
    TBoolean    checkOrder;


    /* Turn off order checking in SEARCH_ROUTINE_NAME, since we may have
     * repeated X values. We do our own checking below, if checkOrder
     * is set to TRUE (the default), which permits repeated X values.
     * (In other words it is OK to have X(i) = X(i+1).)
     * Discontinuities are expressed by having repeated X values,
     * with different associated F values. 
     */
    checkOrder = JpmcdsBSearchCheckOrder(FALSE);

    if (checkOrder)
    {    
        int idx;
        for (idx = 0; idx < N-1; idx++)
        {
            if (XARRAY(idx) > XARRAY(idx+1))
            {
                JpmcdsErrMsg
                ("%s: Domain array element %d (%f) < element %d (%f).\n",
                 STRING(ROUTINE_NAME), 
                 idx+1, (double)XARRAY(idx+1), 
                 idx, (double)XARRAY(idx));           
                 return FAILURE;
            }
        }
    }


    if (SEARCH_ROUTINE_NAME(xDesired, x, xskip, N, &lo, &hi) == FAILURE)
        goto done;                      /* Failed */

    /* Now reset checking to original value.
     */
    JpmcdsBSearchCheckOrder(checkOrder);


    /* If we're not on a discontinuity, perform interpolation as before.
     * Otherwise, we take the lefthand or righthand limit, depending
     * on the value of argument rightLimit.
     */
    if (lo == hi)
    {
        *fInterp = FARRAY(lo);           /* Only one point in curve */
    }
    else if (xDesired != XARRAY(lo))   /* Didn't find exact X point */
    {
        if (mfunc ==  NULL) /* Use default */
            step_fraction = (double)(xDesired-XARRAY(lo)) /
                            (double)(XARRAY(hi) - XARRAY(lo));
        else
            step_fraction = (*mfunc)(xDesired, (double)XARRAY(lo)) /
                            (*mfunc)((double)XARRAY(hi), (double)XARRAY(lo));
        *fInterp = FARRAY(lo) +   step_fraction * (FARRAY(hi) - FARRAY(lo));
    }
    else                                 /* We have found xDesired exactly. */
    {
        /* Make sure we get the min "lo" index which has the same X value.
        */ 
        while (lo > 0 && XARRAY(lo) == XARRAY(lo-1))
        {
            lo--;
        }

       /* Now reset hi to be the biggest index such that XARRAY(hi)
        * equals XARRAY(lo). This means we are ignoring the value
        * hi which comes out of SEARCH_ROUTINE_NAME. 
        */ 
        for (hi=lo; hi < N-1 && XARRAY(hi+1) == XARRAY(lo); )
        {
            hi++;
        }
        if (hi - lo >= 2)                /* Take middle value if 3 x values */
        {                                /* have the same value */
            *fInterp = FARRAY((hi+lo)/2);
        }
        else if (rightLimit)             /* x[lo] = x[hi], Use rightLimit */
        {
            *fInterp = FARRAY(hi);
        }
        else                             /* x[lo] = x[hi], Use leftLimit */
        {
            *fInterp = FARRAY(lo);
        }
    }

    status = SUCCESS;
done:
    if (status == FAILURE)
        JpmcdsErrMsg("%s: Failed.\n", routine);

    return status;
}





#undef XARRAY
#undef FARRAY
