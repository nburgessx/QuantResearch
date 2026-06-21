/*
 * ISDA CDS Standard Model
 *
 * Copyright (C) 2009 International Swaps and Derivatives Association, Inc.
 * Developed and supported in collaboration with Markit
 * 
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the ISDA CDS Standard Model Public License.
 */

#ifndef MDYDATE_H
#define MDYDATE_H

#ifdef __cplusplus
extern "C"
{
#endif

typedef struct _TMDYDate
{
    long month;                         /* In range [1,12] */
    long day;                           /* In range [1-31] */
    long year;                          /* In range [1600-] */
} TMonthDayYear;

#ifdef __cplusplus
}
#endif

#endif    /* MDYDATE_H */
