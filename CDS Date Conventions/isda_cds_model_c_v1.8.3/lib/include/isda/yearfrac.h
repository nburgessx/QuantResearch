/*
 * ISDA CDS Standard Model
 *
 * Copyright (C) 2009 International Swaps and Derivatives Association, Inc.
 * Developed and supported in collaboration with Markit
 * 
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the ISDA CDS Standard Model Public License.
 */

#ifndef YEARFRAC_H
#define YEARFRAC_H

#include "cgeneral.h"
#include "bastypes.h"
#include "metric.h"


#ifdef __cplusplus
extern "C"
{
#endif


/*f
***************************************************************************
** Returns metric function to use for interpolation routines for a given 
** day count convention.
***************************************************************************
*/
TMetricDoubleFunc JpmcdsDayCountToMetricFunc(long dayCountConv);


/*f
***************************************************************************
** Maps day count conventions to (private) strings.
***************************************************************************
*/
char* JpmcdsFormatDayCountConv(long dayCountConv);


/*f
***************************************************************************
** Converts String to Day Count Convention.
**
** This function is fairly permissive in terms of what it will accept.
** For example, "Actual/360", "A", etc.
** Strings of the form "30/Actual" are disallowed, however.
***************************************************************************
*/
int JpmcdsStringToDayCountConv(char *dayCountString, long *type);


/*f
***************************************************************************
** Returns SUCCESS if day count convention is valid, FAILURE otherwise.
**
** If routine != NULL, then populates the error log on FAILURE. Otherwise
** if routine IS NULL, then the error log is not populated on FAILURE.
***************************************************************************
*/
int JpmcdsDayCountValid(char *routine, long dayCountConv);


#ifdef __cplusplus
}
#endif

#endif    /* YEARFRAC_H */
