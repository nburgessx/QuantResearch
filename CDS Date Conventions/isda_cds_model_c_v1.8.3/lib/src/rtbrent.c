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
#include <ctype.h> 
#include "macros.h"
#include "cgeneral.h" 
#include "cerror.h" 
#include "rtbrent.h"


/*
***************************************************************************
** Uses the secant method until one of the following things happens. On the 
** right, the setting of foundIt, bracketed, and the return status is given. 
**
**                                      foundIt  bracketed  status
**                                      -------  ---------  ------
**  1) the root is found                TRUE     TRUE       SUCCESS
**  2) root is bracketed                FALSE    TRUE       SUCCESS
**  3) max # steps reached              FALSE    FALSE      SUCCESS
**  4) next x not in [boundLo, boundHi] FALSE    FALSE      SUCCESS
**  5) Objective func returns FAILURE   FALSE    FALSE      FAILURE
**
**  Note that the routine only returns FAILURE if the objective function 
**returns FAILURE. 
***************************************************************************
*/
static int secantMethod 
    (TObjectFunc funcd,        /* (I) Name of function to call */
     void       *data,         /* (I) Data to pass into funcd */
     int        numIterations, /* (I) Maximum number of iterations */
     double     xacc,          /* (I) X accuracy */
     double     facc,          /* (I) F accuracy */
     double     boundLo,       /* (I) Lower bound */
     double     boundHi,       /* (I) Upper bound */
     double     *xPoints,      /* (I/O) Array of x values */
     double     *yPoints,      /* (I/O) Array of y values */
     TBoolean   *foundIt,      /* (O) If solution was found */
     TBoolean   *bracketed,    /* (O) If root was bracketed */
     double     *solution);    /* (O) Root of equation */


/*
***************************************************************************
** Finds the root using a combination of inverse quadratic method and 
** bisection method.            
***************************************************************************
*/
static int brentMethod 
    (TObjectFunc funcd,         /* (I) Name of function to call */
     void       *data,         /* (I) Data to pass into funcd */
     int        numIterations, /* (I) Maximum number of iterations*/
     double     xacc,          /* (I) X accuracy tolerance */
     double     facc,          /* (I) Function accuracy tolerance */
     double     *xPoints,      /* (I) Array of x values */
     double     *yPoints,      /* (I) Array of y values */
     double     *solution);    /* (O) X root found */


#define SWITCH(a,b)  \
{                    \
    double temp = a; \
    a = b;           \
    b = temp;        \
}               

/* 
 * This is the percentage of x=(boundHi - boundLo) used to set x[2] 
 * if x[2] is outside [boundLo, boundHi]. 
 */
#define ONE_PERCENT 0.01


/* file scope as private functions report same function name */
static char routine[] = "JpmcdsRootFindBrent"; 


/*
***************************************************************************
** Finds the root of f(x) = 0 using a combination of secant, bisection and 
** an inverse quadratic interpolation method.
***************************************************************************
*/
int JpmcdsRootFindBrent(
   TObjectFunc funcd,                   /* (I) Function to call */
   void       *data,                    /* (I) Data to pass into funcd */
   double      boundLo,                 /* (I) Lower bound on legal X */
   double      boundHi,                 /* (I) Upper bound on legal X */
   int         numIterations,           /* (I) Maximum number of iterations */
   double      guess,                   /* (I) Initial guess */
   double      initialXStep,            /* (I) Size of step in x */
   double      initialFDeriv,           /* (I) Initial derivative or 0*/
   double      xacc,                    /* (I) X accuracy tolerance */
   double      facc,                    /* (I) Function accuracy tolerance */
   double      *solution)               /* (O) Root found */

{
   double      fLo;                    /* Function evaluated at boundLo */
   double      fHi;                    /* Function evaluated at boundHi */
   TBoolean    bracketed;              /* If root was bracketed by secant */
   TBoolean    foundIt;                /* If root was found by secant */
   double      xPoints[3];             /* Array of x values */
   double      yPoints[3];             /* Array of y values */
   double      boundSpread;            /* Diff between the hi and lo bounds */ 


   xPoints[0] = guess;

   /* Make sure lower bound is below higher bound.
    */
   if (boundLo >= boundHi)
   {
       JpmcdsErrMsg ("%s: Lower bound (%2.6e) >= higher bound (%2.6e).\n", routine, boundLo, boundHi);
       return FAILURE;
   }

   /* Check to make sure the guess is within the bounds 
    */
   if (xPoints[0]<boundLo || xPoints[0]>boundHi)
   {
       JpmcdsErrMsg ("%s: Guess (%2.6e) is out of range [%2.6e,%2.6e].\n", routine, xPoints[0], boundLo, boundHi);
       return FAILURE;
   }

   if ((*funcd)(xPoints[0],data,&yPoints[0])==FAILURE) 
   {
       JpmcdsErrMsg("%s: Supplied function failed at point %2.6e.\n", routine, xPoints[0]);
       return FAILURE;
   }

   /* Check if guess is the root (use bounds for x-accuracy) */
   if (yPoints[0] == 0.0 || 
       (ABS(yPoints[0]) <= facc && (ABS(boundLo-xPoints[0]) <= xacc ||
                                    ABS(boundHi-xPoints[0]) <= xacc)))
   {
       *solution = xPoints[0]; 
       return SUCCESS;
   }

   /* If the initialXStep is 0, set it to ONE_PERCENT of 
    *  of (boundHi - BoundLo). 
    */
   boundSpread  = boundHi - boundLo; 
   if ( initialXStep == 0.0 ) 
   {
       initialXStep = ONE_PERCENT * boundSpread; 
   }
 
   /* Take a step of the size passed in, if the derivative is not
      passed in */
   if (initialFDeriv == 0)
   {
       xPoints[2] = xPoints[0] + initialXStep;
   }
   else
   {
   /* If initial derivative is known, use Newton's Method 
    * it to find the next point.
    */
       xPoints[2] = xPoints[0] - (yPoints[0])/initialFDeriv;
   }

   /* Now check to make sure that xPoints[2] 
    * is within the hi-lo bounds. If it isn't, then adjust it 
    * so that it is. 
    */
   if ( xPoints[2] < boundLo || xPoints[2] > boundHi ) 
   {
       /* Switch the direction of the step */ 
       xPoints[2] = xPoints[0] - initialXStep; 
       if ( xPoints[2] < boundLo )
       {
           /* xPoints[2] is still too small, so we make 
            * it boundLo
            */ 
           xPoints[2] = boundLo; 
       }
       
       if ( xPoints[2] > boundHi )
       {
           /* xPoints[2] is too large, then set it 
            * to boundHi. 
            */ 
           xPoints[2] = boundHi; 
       }

       if ( xPoints[2] == xPoints[0] )
       {
           /* We cannot have xPoints[0] and 
            * xPoints[2] be the same. 
            */ 
           if ( xPoints[2] == boundLo ) 
           {
               xPoints[2] = boundLo + ONE_PERCENT * boundSpread; 
           }
           else 
           {
               xPoints[2] = boundHi - ONE_PERCENT * boundSpread; 
           }
       }
   }
       
   /* Finally, try to call (*funcd) with xPoints[2], to make 
    * that the function can return a value at that point. 
    */
   if ((*funcd)(xPoints[2],data,&yPoints[2])==FAILURE) 
   {
       JpmcdsErrMsg("%s: Supplied function failed at point %2.6e\n.", routine,xPoints[2]);
       return FAILURE;
   }


   /* Check if the new point meets the tolerance requirements */
   if (yPoints[2] == 0.0 || 
       (ABS(yPoints[2]) <= facc  && ABS(xPoints[2]-xPoints[0]) <= xacc))
   {
       *solution = xPoints[2]; 
       return SUCCESS;
   }

   /* Call secant method to find the root, or to get a 
      third point, so that two of the three points bracket the root. */
   if (secantMethod (funcd, data, numIterations,
                     xacc, facc, boundLo, boundHi,
                     xPoints, yPoints, 
                     &foundIt, &bracketed, solution) == FAILURE)
   {
       return FAILURE;
   }
   else if (foundIt)
   {
       return SUCCESS;
   }
   else if (bracketed)
   {
       goto bracketed;
   }


   /* Root was not bracketed, now try at the bounds 
    */
   if ((*funcd)(boundLo,data,&fLo)==FAILURE) 
   {
       JpmcdsErrMsg("%s: Supplied function failed at point %2.6e.\n", routine, boundLo);
       return FAILURE;
   }

   if (fLo == 0.0 || (ABS(fLo) <= facc && ABS(boundLo - xPoints[0]) <= xacc))
   {
       *solution = boundLo; 
       return SUCCESS;
   }
   
   /* If these points bracket the root, assign points so that
      xPoints[0] < xPoints[2] */
   if (yPoints[0]*fLo < 0)
   {
       xPoints[2] = xPoints[0];
       xPoints[0] = boundLo;

       yPoints[2] = yPoints[0];
       yPoints[0] = fLo;
       
   }
   else 
   /* Root is still not bracketed, so try at the upper bound now. */
   {
       if ((*funcd)(boundHi,data,&fHi)==FAILURE)
       {
           JpmcdsErrMsg("%s: Supplied function failed at point %2.6e.\n", routine, boundHi);
           return FAILURE;
       }
       if (fHi == 0.0 || (ABS(fHi) <= facc && ABS(boundHi-xPoints[0]) <= xacc)) 
       {
           *solution = boundHi; 
           return SUCCESS;
       }

       /* If points bracket the root, assign xPoints[2] to boundHi */
       if (yPoints[0]*fHi < 0)
       {
           xPoints[2] = boundHi;
           yPoints[2] = fHi;
       }
       else 
       /* Root could not be bracketed at the bounds. */
       {
           JpmcdsErrMsg("%s: Function values (%2.6e,%2.6e) at bounds\n"
                     "\t(%2.6e, %2.6e) imply no root exists.\n",
                     routine, fLo, fHi, boundLo, boundHi);
           return FAILURE;
       }
   }        

   /* xPoints[0] and xPoints[2] bracket the root, but we need third
      point to do Brent method. Take the midpoint. */
   xPoints[1] = 0.5*(xPoints[0]+xPoints[2]);
   if ((*funcd)(xPoints[1],data,&yPoints[1])==FAILURE) 
   {
       JpmcdsErrMsg("%s: Supplied function failed at point %2.6e.\n", routine, xPoints[1]);
       return FAILURE;
   }
   if (yPoints[1] == 0.0 || 
       (ABS(yPoints[1]) <= facc && ABS(xPoints[1] - xPoints[0]) <= xacc)) 
   {
       *solution = xPoints[1]; 
       return SUCCESS;
   }
   
 bracketed:
   /* xPoints is an array of three points, two of which bracket the root.
      Call brent Method now, to find the root */
   if (brentMethod (funcd, data, numIterations, xacc, facc,
                    xPoints, yPoints, solution)  == FAILURE)
   {
       return FAILURE;
   }
   else
   {
       return SUCCESS;
   }
}


/*
***************************************************************************
** Finds the root using a combination of inverse quadratic method and 
** bisection method.            
***************************************************************************
*/
static int brentMethod 
    (TObjectFunc  funcd,        /* (I) Function to evaluate */
     void        *data,         /* (I) Data to pass into funcd */
     int         numIterations, /* (I) Maximum number of iterations*/
     double      xacc,          /* (I) X accuracy tolerance */
     double      facc,          /* (I) Function accuracy tolerance */
     double      *xPoints,      /* (I) Array of x values */
     double      *yPoints,      /* (I) Array of y values */
     double      *solution)     /* (O) X root found */
{
    int         j;                      /* Index */
    double      ratio;                  /* (x3-x1)/(x2-x1) */
    double      x31;                    /* x3-x1*/
    double      x21;                    /* x2-x1*/
    double      f21;                    /* f2-f1 */
    double      f31;                    /* f3-f1 */
    double      f32;                    /* f3-f2 */
    double      xm;                     /* New point found using Brent method*/
    double      fm;                     /* f(xm) */

    double      x1 = xPoints[0];        /* xN short hand for xPoints[n] */
    double      x2 = xPoints[1];
    double      x3 = xPoints[2];

    double      f1 = yPoints[0];
    double      f2 = yPoints[1];
    double      f3 = yPoints[2];

    for (j=1; j<=numIterations; j++) 
    {
        /* Always want to be sure that f1 and f2 have opposite signs,
         * and f2 and f3 have the same sign.
         */
        if (f2*f1>0.0)
        {   
            SWITCH(x1,x3);
            SWITCH(f1,f3);
        }
        f21 = f2-f1;
        f32 = f3-f2;
        f31 = f3-f1;
        x21 = x2-x1;
        x31 = x3-x1;
        /* Check whether this is suitable for interpolation. When checking
         * for f21,f31,f32 = 0, we don't use IS_ALMOST_ZERO for efficiency
         * reasons. If the objective function has been truncated to a 
         * certain number of digits of accuracy, f21,f31,or f32 could be
         * (were in one case) zero. In this case we need to protect against
         * division by zero. So we use bisection instead of brent.
         */
        ratio = (x3-x1)/(x2-x1);
        if (f3*f31<ratio*f2*f21 || f21 == 0. || f31 == 0. || f32 == 0.) 
        {
            /* This is not suitable, do bisection 
             */
            x3 = x2;
            f3 = f2; 

        }
        else 
        {
            xm = x1 - (f1/f21)*x21 + ((f1*f2)/(f31*f32))*x31 - 
                ((f1*f2)/(f21*f32))*x21;
            if ((*funcd)(xm,data,&fm)==FAILURE)
            {
                JpmcdsErrMsg("%s: Supplied function failed at point %2.6e.\n", routine,xm);
                return FAILURE;
            } 
            if (fm == 0.0 || (ABS(fm) <= facc && ABS(xm-x1) <= xacc))
            {
                *solution=xm;
                return SUCCESS;
            }
            /* If the new point and point1 bracket the root,
               replace point3 with the new point */
            if (fm*f1<0.0)
            {
                x3=xm;
                f3=fm;
            }
            /* If the root is not bracketed, replace point1 with new point,
               and point3 with point2 */
            else
            {
                x1=xm;
                f1=fm;
                x3=x2;
                f3=f2;
            }
        }             
        x2 = 0.5*(x1+x3); 
        if ((*funcd)(x2,data,&f2)==FAILURE)
        {
            JpmcdsErrMsg("%s: Supplied function failed at point %2.6e.\n", routine, x2);
            return FAILURE;
        } 
        if (f2 == 0.0 || (ABS(f2) <= facc && ABS(x2 - x1) <= xacc))
        {
            *solution = x2;
            return SUCCESS;
        }
    }

    JpmcdsErrMsg("%s: Maximum number of iterations exceeded.\n", routine);
    return FAILURE;
}    


/*
***************************************************************************
** Uses the secant method until one of the following things happens. On the 
** right, the setting of foundIt, bracketed, and the return status is given. 
**
**                                      foundIt  bracketed  status
**                                      -------  ---------  ------
**  1) the root is found                TRUE     TRUE       SUCCESS
**  2) root is bracketed                FALSE    TRUE       SUCCESS
**  3) max # steps reached              FALSE    FALSE      SUCCESS
**  4) next x not in [boundLo, boundHi] FALSE    FALSE      SUCCESS
**  5) Objective func returns FAILURE   FALSE    FALSE      FAILURE
**
**  Note that the routine only returns FAILURE if the objective function 
**returns FAILURE. 
***************************************************************************
*/
static int secantMethod
    (TObjectFunc funcd,        /* (I) Function to evaluate */
     void       *data,         /* (I) data to pass into funcd */
     int        numIterations, /* (I) Maximum number of iterations */
     double     xacc,          /* (I) Accuracy of x */
     double     facc,          /* (I) Accuracy of f */
     double     boundLo,       /* (I) Lower bound */
     double     boundHi,       /* (I) Upper bound */
     double     *xPoints,      /* (I/O) Array of x points */
     double     *yPoints,      /* (I/O) Array of y points */
     TBoolean   *foundIt,      /* (O) If solution was found */
     TBoolean   *bracketed,    /* (O) if root was bracketed*/
     double     *solution)     /* (O) Root of function */
{
    int           j=numIterations;      /* Index */
    double        dx;                   /* Delta x used for secant */       

    *foundIt = FALSE;           /* Until solution is found. */
    *bracketed = FALSE;         /* Until bracketed */

    while (j--)
    {
        /* Swap points so that yPoints[0] is smaller than yPoints[2] */
        if (ABS(yPoints[0])>ABS(yPoints[2]))
        {   
            SWITCH(xPoints[0],xPoints[2]);
            SWITCH(yPoints[0],yPoints[2]);
        }

        /* Make sure that you do not divide by a very small value */
        if (ABS(yPoints[0]-yPoints[2]) <= facc)
        {
            if (yPoints[0] - yPoints[2] > 0)
            {
                dx = -yPoints[0]*(xPoints[0] - xPoints[2])/facc;
            }
            else
            {
                dx = yPoints[0]*(xPoints[0] - xPoints[2])/facc;
            }
        }
        else
        {
            dx= (xPoints[2]-xPoints[0])* yPoints[0]/(yPoints[0]-yPoints[2]);
        }
        xPoints[1] = xPoints[0] + dx;

        /* Make sure that the point is within bounds 
         */
        if (xPoints[1] < boundLo || xPoints[1] > boundHi)
        {
            return SUCCESS;     /* Not bracketed, not found */
        }

        if ((*funcd)(xPoints[1],data,&yPoints[1])==FAILURE) 
        {
            JpmcdsErrMsg("%s: Supplied function failed at point %2.6e.\n", routine, xPoints[1]);
            return FAILURE;     /* Not bracketed, not found */
        }
        if (yPoints[1] == 0.0 || 
            (ABS(yPoints[1]) <= facc && ABS(xPoints[1] - xPoints[0]) <= xacc))
        {
            *solution  = xPoints[1];
            *foundIt   = TRUE;
            *bracketed = TRUE;
            return SUCCESS;     /* Found, bracketed */
        }

        if ((yPoints[0] < 0 && yPoints[1] < 0 && yPoints[2] < 0) ||
            (yPoints[0] > 0 && yPoints[1] > 0 && yPoints[2] > 0))
        {
            /* Swap points so that yPoints[0] is always smallest 
             */
            if (ABS(yPoints[0]) > ABS(yPoints[1]))
            {
                xPoints[2] = xPoints[0];
                yPoints[2] = yPoints[0];
                xPoints[0] = xPoints[1];
                yPoints[0] = yPoints[1];
            }
            else 
            {
                xPoints[2] = xPoints[1];
                yPoints[2] = yPoints[1];
            }
            continue;
        }
        else
        { 
            /* Root was bracketed. 
             * Swap points so that yPoints[0]*yPoints[2] < 0 
             */
            if (yPoints[0]*yPoints[2] > 0)
            {
                /* Make sure that you swap so that 
                 * xPoints[0]<xPoints[1]<xPoints[2] 
                 */
                if (xPoints[1] < xPoints[0])
                {
                    SWITCH(xPoints[0], xPoints[1]);
                    SWITCH(yPoints[0], yPoints[1]);
                }
                else
                {
                    SWITCH (xPoints[1], xPoints[2]);
                    SWITCH (yPoints[1], yPoints[2]);
                }
            }
            /* Root was bracketed, but not found.
             */
            *bracketed = TRUE;
            return SUCCESS;
        }
    } /* while */


    /* Root not bracketed or found.
     */
    return SUCCESS;
}
