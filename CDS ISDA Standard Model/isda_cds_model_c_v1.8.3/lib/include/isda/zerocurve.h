/*
 * ISDA CDS Standard Model
 *
 * Copyright (C) 2009 International Swaps and Derivatives Association, Inc.
 * Developed and supported in collaboration with Markit
 * 
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the ISDA CDS Standard Model Public License.
 */

#ifndef ZEROCURVE_H_
#define ZEROCURVE_H_

#include "bastypes.h"


#ifdef __cplusplus
extern "C"
{
#endif


/*f
***************************************************************************
** Build zero curve from money market, and swap instruments.
***************************************************************************
*/
EXPORT TCurve* JpmcdsBuildIRZeroCurve(
    TDate      valueDate,      /* (I) Value date                       */
    char      *instrNames,     /* (I) Array of 'M' or 'S'              */
    TDate     *dates,          /* (I) Array of swaps dates             */
    double    *rates,          /* (I) Array of swap rates              */
    long       nInstr,         /* (I) Number of benchmark instruments  */
    long       mmDCC,          /* (I) DCC of MM instruments            */
    long       fixedSwapFreq,  /* (I) Fixed leg freqency               */
    long       floatSwapFreq,  /* (I) Floating leg freqency            */
    long       fixedSwapDCC,   /* (I) DCC of fixed leg                 */
    long       floatSwapDCC,   /* (I) DCC of floating leg              */
    long       badDayConv,     /* (I) Bad day convention               */
    char      *holidayFile);   /* (I) Holiday file                     */


#ifdef __cplusplus
}
#endif

#endif /*ZEROCURVE_H_*/
