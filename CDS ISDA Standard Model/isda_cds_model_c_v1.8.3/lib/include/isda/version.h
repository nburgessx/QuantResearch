/*
 * ISDA CDS Standard Model
 *
 * Copyright (C) 2009 International Swaps and Derivatives Association, Inc.
 * Developed and supported in collaboration with Markit
 * 
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the ISDA CDS Standard Model Public License.
 */

#ifndef VERSION_H
#define VERSION_H

#include "cgeneral.h"


#ifdef __cplusplus
extern "C"
{
#endif

/*
***************************************************************************
** Return library version number.
***************************************************************************
*/
EXPORT int JpmcdsVersionString(char*);

#ifdef __cplusplus
}
#endif

#endif
