/*
 * ISDA CDS Standard Model
 *
 * Copyright (C) 2009 International Swaps and Derivatives Association, Inc.
 * Developed and supported in collaboration with Markit
 * 
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the ISDA CDS Standard Model Public License.
 */

#ifndef FUNCDEF_H
#define FUNCDEF_H

#include <setjmp.h>
#include "cgeneral.h"


#ifdef __cplusplus
extern "C"
{
#endif


/*t
 * TParamDef defines a single parameter of a wrapper function.
 */
typedef struct _TParamDef
{
    /*
    ** Name of the parameter. Used in function wizards, error messages etc.
    */
    char     *name;
    /*
    ** Description of add-in function.
    */
    char*      description;
} TParamDef;

/*t
 * TFuncDef defines an add-in function supported by the interface.
 */
typedef struct _TFuncDef
{
    /*
    ** Name of add-in function.
    */
    char*      name;
    /*
    ** Number of inputs.
    */
    size_t     numInputs;
    /*
    ** Definition of input parameters.
    */
    TParamDef *params;
    /*
    ** Description of add-in function.
    */
    char*      description;
} TFuncDef;


#ifdef __cplusplus
}
#endif

#endif    /* FUNCDEF_H */

