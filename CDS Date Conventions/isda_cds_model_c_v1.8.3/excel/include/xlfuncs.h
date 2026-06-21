/*
 * ISDA CDS Standard Model
 *
 * Copyright (C) 2009 International Swaps and Derivatives Association, Inc.
 * Developed and supported in collaboration with Markit
 * 
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the ISDA CDS Standard Model Public License.
 */

#ifndef XLFUNCS_H
#define XLFUNCS_H

#include "xlutils.h"
#include "xlapi.h"
#include "cgeneral.h"

/*f
***************************************************************************
** Standard Excel addin entry point.
***************************************************************************
*/
XLFUNC(short) xlAutoOpen(void);


/*f
***************************************************************************
** Standard Excel addin cleanup function.
***************************************************************************
*/
XLFUNC(short) xlAutoClose(void);


/*f
***************************************************************************
** Standard Excel addin function to release parameter memory.
***************************************************************************
*/
XLFUNC(void) xlAutoFree(XLOPER *x);


/*f
***************************************************************************
** Standard Excel function for addin loaded event.
***************************************************************************
*/
XLFUNC(int) xlAutoAdd(void);


/*f
***************************************************************************
** Standard Excel function for addin unloaded event.
***************************************************************************
*/
XLFUNC(int) xlAutoRemove(void);


/*f
***************************************************************************
** Standard Excel function to provide the library name to Excel.
***************************************************************************
*/
XLFUNC(LPXLOPER) xlAddInManagerInfo(LPXLOPER xAction);

#endif
