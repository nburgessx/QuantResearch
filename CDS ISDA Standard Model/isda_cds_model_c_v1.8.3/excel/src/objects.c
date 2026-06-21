/*
 * ISDA CDS Standard Model
 *
 * Copyright (C) 2009 International Swaps and Derivatives Association, Inc.
 * Developed and supported in collaboration with Markit
 * 
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the ISDA CDS Standard Model Public License.
 */

#include <string.h>
#include "cgeneral.h"
#include "memory.h"
#include "macros.h"
#include "cerror.h"
#include "tcurve.h"

/*t
 * Holds the name, version and data for a spreadsheet object.
 */
typedef struct _TObject
{
    char *name;
    int   version; /* increments by 1 for every new store operation */
    void *data;
    void *next;
} TObject;


static TObject *cache = NULL;


#define SEPARATOR '#'


/*
***************************************************************************
** Find node for object name.  Returns cache node (NULL if not found).
***************************************************************************
*/
TObject* FindNode(char* name)
{
    TObject *node = cache;
    while (node)
    {
        if (stricmp(node->name, name) == 0)
            break;

        node = (TObject*) node->next;
    }
    
    return node;
}


/*
***************************************************************************
** Store data for addin object.  Returns handle (NULL if error).
***************************************************************************
*/
char* StoreObject(char* name, void* data)
{
    static char *routine = "StoreObject";
    int          status = FAILURE;
    char        *handle = NULL;
    TObject     *node;

    static char buffer[255];

    if (data == NULL)
    {
        JpmcdsErrMsg("%s: No data to store provided.\n", routine);
        goto done;
    }

    if (name == NULL)
    {
        JpmcdsErrMsg("%s: No object name provided.\n", routine);
        goto done;
    }

    if (strlen(name) > 200)
    {
        JpmcdsErrMsg("%s: Object name cannot exceed 200 characters.\n", routine);
        goto done;
    }

    node = FindNode(name);
    if (node == NULL)
    {
        node = NEW(TObject);
        if (node == NULL)
            goto done;

        node->name = strdup(name);
        node->version = 1;
        node->data = data;
        node->next = cache;
        cache = node;
    }
    else
    {
        node->version += 1;
        JpmcdsFreeTCurve((TCurve*)(node->data));
        node->data = data;
    }

    sprintf(buffer, "%s%c%d\0", node->name, SEPARATOR, node->version);
    handle = strdup(buffer);
    if (handle == NULL)
        goto done;

    strcpy(handle, buffer);
    status = SUCCESS;
    
done:
    if (status != SUCCESS)
    {
        JpmcdsErrMsg("%s: Failed!\n", routine);
        FREE(handle);
        handle = NULL;
    }

    return handle;
}


/*
***************************************************************************
** Retrieve data for addin object.  Returns object (NULL if error).
***************************************************************************
*/
void* RetrieveObject(char* handle)
{
    static char *routine = "RetrieveObject";
    int          status = FAILURE;
    void        *data = NULL;
    char        *name = NULL;
    TObject     *node;
    char        *sep;

    if (handle == NULL)
    {
        JpmcdsErrMsg("%s: No object handle specified.\n", routine);
        goto done;
    }

    name = strdup(handle);
    if (name == NULL)
        goto done;

    /* remove version number from lookup string */
    sep = strchr(name, SEPARATOR);
    if (sep != NULL)
        *sep = '\0';

    node = FindNode(name);
    if (node == NULL)
    {
        JpmcdsErrMsg("%s: No object called %s found.\n", routine, name);
        goto done;
    }

    data = node->data;
    status = SUCCESS;

done:
    if (status != SUCCESS)
    {
        JpmcdsErrMsg("%s: Failed!\n", routine);
        data = NULL;
    }

    FREE(name);
    return data;
}


/*
***************************************************************************
** Clears all objects from cache.  
***************************************************************************
*/
void FreeObjects()
{
    TObject *next = cache;
    TObject *node;

    while (next)
    {
        node = next;
        next = (TObject*) node->next;
        FREE(node->name);
        JpmcdsFreeTCurve((TCurve*)(node->data));
        FREE(node);
    }
}
