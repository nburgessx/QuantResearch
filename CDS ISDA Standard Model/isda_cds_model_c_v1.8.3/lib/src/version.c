/*
 * ISDA CDS Standard Model
 *
 * Copyright (C) 2009 International Swaps and Derivatives Association, Inc.
 * Developed and supported in collaboration with Markit
 * 
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the ISDA CDS Standard Model Public License.
 */

#include "version.h"

#include <stdio.h>
#include <string.h>

#define FIDR_LITERAL(x) FIDR_LITERAL_(x)
#define FIDR_LITERAL_(x) #x


/*
***************************************************************************
** Return library version number.
***************************************************************************
*/
EXPORT int JpmcdsVersionString (char *version)
{
    strcpy (version, "CDS version " FIDR_LITERAL(VERSION));
    return 0;
}
