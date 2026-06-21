/*
 * ISDA CDS Standard Model
 *
 * Copyright (C) 2009 International Swaps and Derivatives Association, Inc.
 * Developed and supported in collaboration with Markit
 * 
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the ISDA CDS Standard Model Public License.
 */

#include "xlutils.h"          

#include <stdlib.h>
#include <string.h>

#include "cfileio.h"
#include "convert.h"
#include "macros.h"                     
#include "cerror.h"  
#include "xlapi.h"

#define EXCEL_ERROR_INPUT -2

typedef struct _EnumStruct
{
    BOOL bFuncWiz;
    short hwndXLMain;
} EnumStruct, FAR * LPEnumStruct;


/*
***************************************************************************
** EXCEL_DATE_OFFSET: Offset between excel date and TDate
**
** GetDateSystem:     Function which sets EXCEL_DATE_OFFSET - called when
**                    the offset is zero.
***************************************************************************
*/
static long EXCEL_DATE_OFFSET = 0;
static void GetDateSystem(void);

/*
***************************************************************************
** XSTR_SIZE:    Size of buffer which can always contain an Excel string.
**
** XSTR_TO_CSTR: Macro which converts an Excel string (which contains the
**               length in position zero and is not null-terminated) to a
**               C-string.
***************************************************************************
*/                               
#define XSTR_SIZE 257
#define XSTR_TO_CSTR(xstr,cstr) strncpy(cstr,xstr+1,(unsigned char)xstr[0]); cstr[(unsigned char)xstr[0]]='\0'

static void JpmcdsStringRemoveWhiteSpace(char *s)
{
    char *d;

    d = s;
    while (*s != '\0')
    {
        if (*s!=' ' && *s!='\t' && *s!='\n')
        {
            *d++ = *s;
        }
        s++;
    }
    *d = '\0';
}

/*
***************************************************************************
** Extract date from atomic OPER or return FAILURE.
***************************************************************************
*/
int JpmcdsExcelGetDate (OPER* oper, TDate *value)
{                             
    char buf[XSTR_SIZE]; /* For conversion from string to date */
    
    switch (oper->type)
    {
    case xltypeNum:
        *value = JpmcdsExcelDateToDate (oper->val.num);
        return SUCCESS;
    case xltypeBool:
        if (!oper->val.boolean)
        {
            *value = 0;
            return SUCCESS;
        }
        return FAILURE;

    case xltypeStr:
        XSTR_TO_CSTR(oper->val.str, buf);
        JpmcdsStringRemoveWhiteSpace(buf);
        if (buf[0]=='\0')
        {
            *value = 0;
            return SUCCESS;
        }
        if (JpmcdsStringToDate(buf, value) == SUCCESS)
        {
            return SUCCESS;
        }
        return FAILURE;

    case xltypeMissing:
    case xltypeNil:         /* TDRM 1/3/97 - for compatibility with old XL interface */
        *value = 0; 
        return SUCCESS;

    case xltypeErr:
        return EXCEL_ERROR_INPUT;

    default:
        return FAILURE;
    }
}


/*
***************************************************************************
** Extract double from atomic OPER or return FAILURE.
***************************************************************************
*/
int JpmcdsExcelGetDouble (OPER *oper, double *value)
{                  
    char buf[XSTR_SIZE]; /* For conversion from string to double */
    char *ep;            /* End pointer for strtod */

    switch (oper->type)
    {
    case xltypeNum:
        *value = oper->val.num;
        return SUCCESS;

    case xltypeBool:
        *value = oper->val.boolean;
        return SUCCESS;

    case xltypeStr:
        XSTR_TO_CSTR(oper->val.str, buf);
        *value = strtod (buf, &ep);
        if (*ep == '\0')
        { /* String successfully converted to double */
            return SUCCESS;
        }
        return FAILURE;

    case xltypeMissing:
    case xltypeNil:
        *value = 0.0;   
        return SUCCESS;

    case xltypeErr:
        return EXCEL_ERROR_INPUT;

    default:
        return FAILURE;
    }
}
    

/*
***************************************************************************
** Extract long from atomic OPER or return FAILURE.
***************************************************************************
*/
int JpmcdsExcelGetLong (OPER *oper, long *value)
{                                     
    char buf[XSTR_SIZE]; /* For conversion from string to long */
    char *ep;            /* End pointer for strtol */
 
    switch (oper->type)
    {
    case xltypeNum:
        *value = (long)oper->val.num;
        /* Excel uses doubles - make sure that is actually an integer */        
        if (ARE_ALMOST_EQUAL(oper->val.num, (double)(*value)))
        {
            return SUCCESS;
        }
        return FAILURE;            

    case xltypeStr:
        XSTR_TO_CSTR(oper->val.str, buf);
        *value = strtol (buf, &ep, 10);
        if (*ep == '\0')
        { /* String successfully converted to long */
            return SUCCESS;
        }
        return FAILURE; 
        
    case xltypeBool:
        *value = oper->val.boolean;
        return SUCCESS;        

    case xltypeMissing:
    case xltypeNil:
        *value = 0; 
        return SUCCESS;

    case xltypeErr:
        return EXCEL_ERROR_INPUT;

    default:
        return FAILURE;
    }
}


/*
***************************************************************************
** Extract string from atomic OPER or return FAILURE.
***************************************************************************
*/
int JpmcdsExcelGetString (OPER *oper, char **value)
{
    char buf[XSTR_SIZE]; /* For numeric output and Excel string */
     
    switch (oper->type)
    {
    case xltypeNum:
        if ((double)((long)(oper->val.num)) == oper->val.num)
        {
            /* i.e. it is a whole number */
            sprintf (buf, "%ld", (long)(oper->val.num));
        }
        else
        {
            /* i.e. it is a fractional number */
            sprintf (buf, "%g", oper->val.num);
        }
        break;

    case xltypeBool:
        strcpy(buf, oper->val.boolean ? "TRUE" : "FALSE");
        break;

    case xltypeStr:
        XSTR_TO_CSTR(oper->val.str, buf);
        break;

    case xltypeMissing:
    case xltypeNil:
        return SUCCESS;

    case xltypeErr:
        return EXCEL_ERROR_INPUT;

    default:
        return FAILURE;
    }
    
    *value = NEW_ARRAY(char, strlen(buf) + 1);
    if (*value == NULL)
    {
        return FAILURE;
    }
    strcpy(*value, buf);
    return SUCCESS;
}


/*
***************************************************************************
** Set an XLOPER to a string value.
***************************************************************************
*/
int JpmcdsExcelSetString(char *str, XLOPER *x)
{
    size_t len;

    if (str == NULL) 
        str = "";

    len = strlen(str);
    
    if (len > 255)
    {
        len = 255;
    }
    x->val.str = NEW_ARRAY(char, len+1);
    if (x->val.str == NULL)
    {
        return FAILURE;
    }
    x->xltype = xltypeStr;
    x->val.str[0] = len;
    if (len > 0)
    {
        strncpy (x->val.str+1, str, len);
    }
    return SUCCESS;
}    


/*
***************************************************************************
** Set an XLOPER to a date value.
***************************************************************************
*/
int JpmcdsExcelSetDate(TDate date, XLOPER *x)
{
    if (EXCEL_DATE_OFFSET == 0)
        GetDateSystem();

    x->xltype = xltypeNum;
    x->val.num = date - EXCEL_DATE_OFFSET;
    return SUCCESS;
}    


/*
***************************************************************************
** Convert from excel date to TDate
***************************************************************************
*/
TDate JpmcdsExcelDateToDate (double xlDate)
{         
    TDate tDate;
    if (EXCEL_DATE_OFFSET == 0)
    {
        GetDateSystem();
    }           
    if ((long)xlDate > 0)
    {
        tDate = (TDate)xlDate + EXCEL_DATE_OFFSET;
    }
    else
    {
        tDate = 0;
    }
            
    return (tDate);
}


/*
***************************************************************************
** The Excel date system must be re-calculated for each function call.
** This function clears out the previous value of the date system.
***************************************************************************
*/
void JpmcdsClearExcelDateSystem (void)
{
    EXCEL_DATE_OFFSET = 0;
}


/*
***************************************************************************
** Recalculated Excel date system.
***************************************************************************
*/
static void GetDateSystem(void)
{
    XLOPER x;
    
    EXCEL_DATE_OFFSET = 109205; /* Corresponds to 1/1/1901 */

/*
** We get the date for 1/1/1904.
** If date format is 1904, then the result will be zero.
** Otherwise the result will be 4 years worth.
*/    
    Excel (xlfDate, (LPXLOPER)&x, 3, TempNum(1904), TempNum(1), TempNum(1));
    if (x.xltype & xltypeNum)
    {
        if (x.val.num <= 1) 
        { /* 1904 date system in use */
            EXCEL_DATE_OFFSET += 1462;
        }
    }
}


/*
***************************************************************************
** Frees an Excel oper.
***************************************************************************
*/
void JpmcdsFreeExcelOper (XLOPER *x)
{       
    if (x == NULL)
        return;

    if (x->xltype & xltypeMulti)
    {   
        size_t size;   /* Size of the array */
        size_t pos;    /* Iterator for array */     
        LPXLOPER xpos; /* XLOPER at pos */

        size = x->val.array.rows * x->val.array.columns;

        for (pos = 0; pos < size; ++pos)
        {
            xpos = x->val.array.lparray + pos;
            if (xpos->xltype & xltypeStr)
            {
                FREE_ARRAY (xpos->val.str);
            }
        }
        
        FREE_ARRAY (x->val.array.lparray);
    }
    else if (x->xltype & xltypeStr)
    {
        FREE_ARRAY (x->val.str);
    }

    FREE(x);                        
    return;
}


/*
***************************************************************************
** Convert the output from a common routine into Excel format.
** This either results in an error output (if the output was NULL)
** or an output which will be freed subsequently by xlAutoFree.
***************************************************************************
*/
LPXLOPER ExcelOutput (void *output)
{
    static XLOPER xErr;
    LPXLOPER x;

    if (output == NULL)
    {
        x = ErrorOutput();
    }
    else
    {
        x = (LPXLOPER)output;
        x->xltype |= xlbitDLLFree;
    }

    return x;
}


/*
***************************************************************************
** Returns the error output.
***************************************************************************
*/
LPXLOPER ErrorOutput (void)
{
    static XLOPER xErr;

    xErr.xltype  = xltypeErr;
    xErr.val.err = xlerrNA;

    return &xErr;
}


/*
***************************************************************************
** Returns a string giving the return status from a call to the Excel 
** function.
***************************************************************************
*/
char* JpmcdsExcelFormatReturnStatus (int xlret)
{
    switch (xlret)
    {
    case xlretSuccess:
        return "SUCCESS";
    case xlretAbort:
        return "MACRO HALTED";
    case xlretInvXlfn:
        return "INVALID FUNCTION NUMBER";
    case xlretInvCount:
        return "INVALID NUMBER OF ARGUMENTS";
    case xlretInvXloper:
        return "INVALID OPER STRUCTURE";
    case xlretStackOvfl:
        return "STACK OVERFLOW";
    case xlretFailed:
        return "COMMAND FAILED";
    case xlretUncalced:
        return "UNCALCED CELL";
    default:
        return "UNKNOWN STATUS";
    }
}
