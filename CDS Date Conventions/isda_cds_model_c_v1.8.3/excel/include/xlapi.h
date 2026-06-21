/*
 * ISDA CDS Standard Model
 *
 * Copyright (C) 2009 International Swaps and Derivatives Association, Inc.
 * Developed and supported in collaboration with Markit
 * 
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the ISDA CDS Standard Model Public License.
 */

#ifndef XLAPI_H 
#define XLAPI_H
                     
/* MS-Windows API */                     
#include "windows.h"

/* MS-Windows Excel API */   
#include "xlcall.h"    /* Must come after windows.h */
#include "framewrk.h"  /* Must come after xlcall.h */

/*m
** XLFUNC macro defines the type of function which can be called from
** Excel. For example, XLFUNC (LPXLOPER) expands to
**     __declspec(dllexport) LPXLOPER
** for Windows NT (32-bit)
*/
#define XLFUNC(retType) __declspec(dllexport) retType

/*t
** OPER structure - same shape as XLOPER but less types defined in the union.
**
** This structure was originally defined in an appendix to one of the Excel
** manuals before the release of the Excel development kit (which was for
** Excel4)
*/
typedef struct _oper
{
    union
    {
        double num;                     /* xltypeNum */
        LPSTR str;                      /* xltypeStr */
        WORD boolean;                   /* xltypeBool */
        WORD err;                       /* xltypeErr */
        struct 
        {
            struct _oper *lparray;
            WORD rows;
            WORD columns;
        } array;                        /* xltypeMulti */
    } val;
    WORD type;
} OPER, *LPOPER ;

#endif /* XLAPI_H */
