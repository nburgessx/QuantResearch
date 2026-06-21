/*
 * ISDA CDS Standard Model
 *
 * Copyright (C) 2009 International Swaps and Derivatives Association, Inc.
 * Developed and supported in collaboration with Markit
 * 
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the ISDA CDS Standard Model Public License.
 */

#ifndef OBJECTS_H
#define OBJECTS_H

/*f
***************************************************************************
** Store data for addin object.  Returns handle (NULL if error).
***************************************************************************
*/
char* StoreObject(char* name, void* data);

/*f
***************************************************************************
** Retrieve data for addin object.  Returns object (NULL if error).
***************************************************************************
*/
void* RetrieveObject(char* handle);

/*f
***************************************************************************
** Clears all objects from cache.  
***************************************************************************
*/
void FreeObjects();

#endif
