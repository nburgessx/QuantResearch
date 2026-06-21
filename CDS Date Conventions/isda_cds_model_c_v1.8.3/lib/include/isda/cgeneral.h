/*
 * ISDA CDS Standard Model
 *
 * Copyright (C) 2009 International Swaps and Derivatives Association, Inc.
 * Developed and supported in collaboration with Markit
 * 
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the ISDA CDS Standard Model Public License.
 */

#ifndef CGENERAL_H
#define CGENERAL_H

#include <stdlib.h>
#include <stdio.h>

#ifdef __cplusplus
extern "C"
{
#endif 

#ifdef NO_EXPORT
#undef EXPORT
#define EXPORT
#endif

#if defined WIN32 || defined WIN32   /* Windows */

#ifndef EXPORT
#define EXPORT  __declspec(dllexport)
#endif

#else                                   /* UNIX, etc */

#ifndef EXPORT
#define EXPORT                          /* EXPORT-> nothing */
#endif

#endif                                  /* if (WIN32, etc) */

typedef int TBoolean;

#ifndef TRUE
#define TRUE 1
#endif

#ifndef FALSE
#define FALSE 0
#endif

#ifndef SUCCESS
#define SUCCESS 0
#endif

#ifndef FAILURE
#define FAILURE -1
#endif

/* This special value means one of three things:  "use the default value",
   "this value does not exist", or "it doesn't matter what this value is".
   This symbol must be set to a negative number, although it need not be
   set to -1. */
#ifndef JPMCDS_DEFAULT
#define JPMCDS_DEFAULT -1
#endif

#ifdef __cplusplus
}
#endif

#endif    /* CGENERAL_H */
