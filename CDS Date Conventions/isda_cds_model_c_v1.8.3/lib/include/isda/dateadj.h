/*
 * ISDA CDS Standard Model
 *
 * Copyright (C) 2009 International Swaps and Derivatives Association, Inc.
 * Developed and supported in collaboration with Markit
 * 
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the ISDA CDS Standard Model Public License.
 */

#ifndef DATEADJ_H
#define DATEADJ_H

#include "cgeneral.h"
#include "cdate.h"

#ifdef __cplusplus
extern "C"
{
#endif


/*f
***************************************************************************
** Calculates a date being at a business or calendar days offset from the 
** input date.
***************************************************************************
*/
int JpmcdsDtFwdAdj
   (TDate          startDate,    /* (I) Input date       */
    TDateAdjIntvl *interval,     /* (I) Business days or calendar period  */
    TDate         *result);      /* (O) Resulting date   */


#ifdef __cplusplus
}
#endif

#endif    /* DATEADJ_H */
