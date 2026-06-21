/*
 * ISDA CDS Standard Model
 *
 * Copyright (C) 2009 International Swaps and Derivatives Association, Inc.
 * Developed and supported in collaboration with Markit
 * 
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the ISDA CDS Standard Model Public License.
 */

#ifndef CONTINGENT_LEG_H
#define CONTINGENT_LEG_H

#include "cx.h"

#ifdef __cplusplus
extern "C"
{
#endif

/*f
***************************************************************************
** Computes the PV of a contingent leg as a whole.
**
** For each payment period this is the integral of LGD(t) . Z(t) . dS/dt dt
** where S is the survival function and LGD is the loss given default
** function and Z is the discount function. Discounting is calculated at the
** payment date and not at the observation date.
***************************************************************************
*/
int JpmcdsContingentLegPV
(TContingentLeg *cl,               /* (I) Contingent leg                  */
 TDate           today,            /* (I) No observations before today    */
 TDate           valueDate,        /* (I) Value date for discounting      */
 TDate           stepinDate,       /* (I) Step-in date                    */
 TCurve         *discountCurve,    /* (I) Risk-free curve                 */
 TCurve         *spreadCurve,      /* (I) Spread curve                    */
 double          recoveryRate,     /* (I) Recovery rate                   */
 double         *pv);              /* (O) Present value of contingent leg */

#ifdef __cplusplus
}
#endif

#endif

