/*
 * ISDA CDS Standard Model
 *
 * Copyright (C) 2009 International Swaps and Derivatives Association, Inc.
 * Developed and supported in collaboration with Markit
 * 
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the ISDA CDS Standard Model Public License.
 */

#ifndef CMEMORY_H
#define CMEMORY_H

#include "cgeneral.h"
#include <stdlib.h>

#ifdef __cplusplus
extern "C"
{
#endif


/*f
***************************************************************************
** Memory allocator.
**
** Returns NULL if no memory left.
***************************************************************************
*/
void* JpmcdsMallocSafe(size_t theSize);


/*f
***************************************************************************
** Memory deallocator.
***************************************************************************
*/
void JpmcdsFreeSafe(void *x);     /* (I) */


#ifdef __cplusplus
}
#endif

#endif    /* CMEMORY_H */
