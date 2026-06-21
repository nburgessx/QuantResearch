/*
 * ISDA CDS Standard Model
 *
 * Copyright (C) 2009 International Swaps and Derivatives Association, Inc.
 * Developed and supported in collaboration with Markit
 * 
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the ISDA CDS Standard Model Public License.
 */

#include "xlfuncs.h"
#include <string.h>
#include "funcdef.h"
#include "objects.h"


extern TFuncDef *gtoFuncDefList[];
extern size_t gtoFuncDefCount;

#define JPMCDS_MAX_LIB_NAME 64
char gtoLibName[] = "CDS analytics";

#define FNPREFIX "CDS"
#define FUNCTION_CATEGORY "CDS"



/*
***************************************************************************
** Register addin function with Excel.
***************************************************************************
*/
static void RegisterFunction(LPXLOPER pxDLL, char *prefix, TFuncDef *funcDef);


/*
***************************************************************************
** Standard Excel addin entry point.
***************************************************************************
*/
XLFUNC(short) xlAutoOpen(void)
{
    int      xlret;
    char     buf[128];
    XLOPER   xDLL;
    LPXLOPER pxDLL = &xDLL;
    size_t   i;

    xlret = Excel(xlGetName, pxDLL, 0);
    if (xlret != xlretSuccess)
    {
        sprintf(buf, " Excel get DLL name function failed: %s", JpmcdsExcelFormatReturnStatus(xlret));
        Excel(xlcAlert,0,2,TempStr(buf),TempInt(2));
        return 1;
    }

    for (i = 0 ; i < gtoFuncDefCount ; i++)
        RegisterFunction(pxDLL, FNPREFIX, gtoFuncDefList[i]);

    Excel(xlFree,0,1,pxDLL);
    return 1;
}


/*
***************************************************************************
** Standard Excel addin cleanup function.
***************************************************************************
*/
XLFUNC(short) xlAutoClose(void)
{
    static char funcName[257]; /* Function name */
    OPER        xFuncName;
    size_t      i;

    for (i = 0; i < gtoFuncDefCount; i++)
    {
        /* Define the function name as a pascal string. */
        funcName[0] = strlen(FNPREFIX) + strlen(gtoFuncDefList[i]->name) + 1;
        strcpy (funcName+1, FNPREFIX);
        strcat (funcName, "_");
        strcat (funcName, gtoFuncDefList[i]->name);

        /* Create a string oper from the function name pascal string */
        xFuncName.type =    xltypeStr;
        xFuncName.val.str = funcName;

        /* Call Excel using xlfSetName to "unregister" the function */
        Excel(xlfSetName, 0, 1, &xFuncName);
    }

    FreeObjects();
    return 1;
}


/*
***************************************************************************
** Standard Excel addin function to release parameter memory.
***************************************************************************
*/
XLFUNC(void) xlAutoFree(XLOPER *x)
{
    JpmcdsFreeExcelOper(x);
}


/*
***************************************************************************
** Standard Excel function for addin loaded event.
***************************************************************************
*/
XLFUNC(int) xlAutoAdd(void)
{
    char buf[JPMCDS_MAX_LIB_NAME+32];
    strcpy(buf, " ");
    strcat(buf, gtoLibName);
    strcat(buf, " add-in has been loaded");
    
    /* Display a dialog box indicating that the XLL was successfully added */
    Excel(xlcAlert, 0, 2, TempStr(buf), TempInt(2));

    return 1;
}


/*
***************************************************************************
** Standard Excel function for addin unloaded event.
***************************************************************************
*/
XLFUNC(int) xlAutoRemove(void)
{
    char buf[JPMCDS_MAX_LIB_NAME+32];
    strcpy(buf, " ");
    strcat(buf, gtoLibName);
    strcat(buf, " add-in has been removed");
    
    Excel(xlcAlert, 0, 2, TempStr(buf), TempInt(2));
    
    return 1;
}


/*
***************************************************************************
** Standard Excel function to provide the library name to Excel.
***************************************************************************
*/
XLFUNC(LPXLOPER) xlAddInManagerInfo(LPXLOPER xAction)
{
    static XLOPER xInfo, xIntAction;
    static char xLibName[JPMCDS_MAX_LIB_NAME];
    LPXLOPER pxInfo;

    /*
    ** This code coerces the passed-in value to an integer. This is how the
    ** code determines what is being requested. If it receives a 1, 
    ** it returns a string representing the long name. If it receives 
    ** anything else, it returns a #VALUE! error.
    */
    Excel(xlCoerce, &xIntAction, 2, xAction, TempInt(xltypeInt));

    if(xIntAction.val.w == 1)
    {
        xInfo.xltype = xltypeStr;
        xLibName[0] = strlen(gtoLibName);
        strncpy (xLibName+1, gtoLibName, xLibName[0]);
        xInfo.val.str = xLibName;
    }
    else
    {
        xInfo.xltype = xltypeErr;
        xInfo.val.err = xlerrValue;
    }

    Excel(xlFree, 0, 1, (LPXLOPER) &xIntAction ,0);    
    pxInfo = &xInfo;

    return pxInfo;
}


/*
***************************************************************************
** Helper function for Excel registration
***************************************************************************
*/
static LPXLOPER ParameterDescription(TFuncDef* funcDef, size_t i)
{
    static char description[256];

    if (funcDef->numInputs > i)
    {
        TParamDef *paramDef = (TParamDef*) funcDef->params + i;

        if (paramDef->description != NULL)
        {
            description[0] = ' ';
            strncpy(description+1, paramDef->description, 255);
            description[255] = '\0';
            return TempStr(description);
        }
    }

    return TempStr(" ");
}


/*
***************************************************************************
** Function to register one function given by its add-in function 
** definition.
***************************************************************************
*/
static void RegisterFunction(LPXLOPER  pxDLL, char *prefix, TFuncDef *funcDef)
{
    int status = FAILURE;

    /*
    ** These are all pascal strings to be used with string type OPERs
    ** for calls to Excel. We will want to use standard C-string
    ** functionality while building the strings, so we need an extra
    ** character at the end. Thus 257 is the maximum length.
    */
    static char errMsg[257];
    static char argTypes[33];  /* Argument types - a string of P's */
    static char argNames[257]; /* Argument names */
    static char funcName[257]; /* Function name */
    static char cFuncName[257];/* C-function name */
    static char catName[257];  /* Category name */
    static char funcDesc[257]; /* function description */

    size_t    nameLen;
    size_t    numArgs;
    size_t    lenArgNames;
    TBoolean  namesAreTooLong;
    size_t    maxArgs = sizeof(argTypes)-2;
    size_t    maxLenArgNames = sizeof(argNames)-2;
    size_t    ipos;
    size_t    nargs;

    /* Define the category name as a pascal string. */
    catName[0] = strlen(FUNCTION_CATEGORY);
    strcpy (catName+1, FUNCTION_CATEGORY);

    /* Define the function name as a pascal string. */
    funcName[0] = strlen(prefix) + strlen(funcDef->name) + 1;
    strcpy (funcName+1, prefix);
    strcat (funcName, "_");
    strcat (funcName, funcDef->name);
    strncpy (cFuncName, funcName, sizeof(funcName));
    maxLenArgNames -= (funcName[0] + 2); /* it seems to be this way! */

    /*
    ***********************************************************************
    ** Define the list of argument types as a pascal string.
    **
    ** There is always one output, plus the standard inputs.
    **
    ** If the function is an object constructor then the base name is
    ** needed.
    ***********************************************************************
    */
    numArgs = 1 + funcDef->numInputs;
    if (numArgs > maxArgs)
    {
        sprintf (errMsg+1, "%s: Too many arguments (%d)\n", funcName+1, numArgs);
        goto done;
    }
    argTypes[0] = (char)numArgs;
    memset (argTypes+1, '\0', maxArgs+1);
    memset (argTypes+1, 'P', numArgs);

    /*
    ***********************************************************************
    ** Build up the list of argument names as a pascal string.
    **
    ** For convenience, there will be a comma at the beginning of the
    ** string in the short term!
    **
    ** If this exceeds the maximum length, then we will use a short format.
    ***********************************************************************
    */
    namesAreTooLong = FALSE;
    lenArgNames     = 0;
    argNames[1]     = '\0'; /* Start with an empty string */
    for (ipos = 0; ipos < funcDef->numInputs; ++ipos)
    {
        TParamDef *paramDef;

        paramDef = funcDef->params + ipos;
        nameLen  = strlen(paramDef->name);

        if ((lenArgNames + nameLen + 1) > maxLenArgNames)
        {
            namesAreTooLong = TRUE;
            break; /* No need to do any more... */
        }
        else
        {
            argNames[lenArgNames] = ',';
            ++lenArgNames;
            strcpy (argNames+lenArgNames, paramDef->name);
            lenArgNames += nameLen;
        }
    }

    if (namesAreTooLong)
    {
        sprintf (errMsg+1, "%s: Argument names are too long for Excel. Maximum length = %ld.", funcName+1, maxLenArgNames);
        goto done;
    }

    argNames[0] = (char)(lenArgNames);

    /* Define the function description as a pascal string. */
    if (funcDef->description != NULL)
    {
        funcDesc[0] = ' ';
        strncpy(funcDesc + 1, funcDef->description, 255);
        funcDesc[256] = '\0';
    }
    else
    {
        funcDesc[0] = ' ';
        funcDesc[1] = '\0';
    }

    /*
    ** Since we use TempStr, it doesn't matter what is in the first character
    ** of the strings, since TempStr replaces it with the string length.
    **
    ** As it happens, it should be the correct length already!
    */
    nargs = funcDef->numInputs;
    if (nargs > 20) nargs = 20;    /* Excel has limit of 30 LXOPERS - so cannot register remaining parameters descriptions */

    Excel (xlfRegister, 
           0,               /* not interested in any return details */
           10 + nargs,  /* number of input parameters in this registration call */
           pxDLL,
           TempStr(funcName),
           TempStr(argTypes),
           TempStr(cFuncName),
           TempStr(argNames),
           TempStr(" 1"),
           TempStr(catName),
           TempStr(" "),        /* hot key / short cut */
           TempStr(" "),        /* help id, eg. "alib.hlp!123" */
           TempStr(funcDesc),    /* function description */
           /* argument descriptions */
           ParameterDescription(funcDef, 0),
           ParameterDescription(funcDef, 1),
           ParameterDescription(funcDef, 2),
           ParameterDescription(funcDef, 3),
           ParameterDescription(funcDef, 4),
           ParameterDescription(funcDef, 5),
           ParameterDescription(funcDef, 6),
           ParameterDescription(funcDef, 7),
           ParameterDescription(funcDef, 8),
           ParameterDescription(funcDef, 9),
           ParameterDescription(funcDef, 10),
           ParameterDescription(funcDef, 11),
           ParameterDescription(funcDef, 12),
           ParameterDescription(funcDef, 13),
           ParameterDescription(funcDef, 14),
           ParameterDescription(funcDef, 15),
           ParameterDescription(funcDef, 16),
           ParameterDescription(funcDef, 17),
           ParameterDescription(funcDef, 18),
           ParameterDescription(funcDef, 19));

    status = SUCCESS;

done:

    if (status != SUCCESS)
    {
        Excel (xlcAlert, 0, 2, TempStr(errMsg), TempInt(2));
    }
}
