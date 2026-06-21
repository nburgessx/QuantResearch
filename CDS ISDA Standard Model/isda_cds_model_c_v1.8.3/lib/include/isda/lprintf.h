/*
 * ISDA CDS Standard Model
 *
 * Copyright (C) 2009 International Swaps and Derivatives Association, Inc.
 * Developed and supported in collaboration with Markit
 * 
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the ISDA CDS Standard Model Public License.
 */

#ifndef LPRINTF_H
#define LPRINTF_H

#include <stdarg.h>
#include "cfileio.h"
#include "cgeneral.h"


#ifdef __cplusplus
extern "C"
{
#endif


/*f
*******************************************************************************
** vprintf functionality using a TFile
*******************************************************************************
*/
int JpmcdsVfprintf(TFile *tFile, char *fmt, va_list arg);


#ifdef __cplusplus
}
#endif

#endif    /* LPRINTF_H */

