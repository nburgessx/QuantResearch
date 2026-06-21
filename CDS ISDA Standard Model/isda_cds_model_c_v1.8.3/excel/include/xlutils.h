/*
 * ISDA CDS Standard Model
 *
 * Copyright (C) 2009 International Swaps and Derivatives Association, Inc.
 * Developed and supported in collaboration with Markit
 * 
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the ISDA CDS Standard Model Public License.
 */

#ifndef XLUTILS_H
#define XLUTILS_H
                                                 
#include "cdate.h"                                                 
#include "cgeneral.h"
#include "xlapi.h"


/*f
***************************************************************************
** Extract date from atomic OPER or return FAILURE.
***************************************************************************
*/
int JpmcdsExcelGetDate (OPER* oper, TDate *value);                                   


/*f
***************************************************************************
** Extract double from atomic OPER or return FAILURE.
***************************************************************************
*/
int JpmcdsExcelGetDouble (OPER *oper, double *value);


/*f
***************************************************************************
** Extract long from atomic OPER or return FAILURE.
***************************************************************************
*/
int JpmcdsExcelGetLong (OPER *oper, long *value);


/*f
***************************************************************************
** Extract string from atomic OPER or return FAILURE.
***************************************************************************
*/
int JpmcdsExcelGetString (OPER *oper, char **value);


/*f
***************************************************************************
** Set an XLOPER to a string value.
***************************************************************************
*/
int JpmcdsExcelSetString (char *str, XLOPER *xloper);


/*f
***************************************************************************
** Set an XLOPER to a date value.
***************************************************************************
*/
int JpmcdsExcelSetDate(TDate date, XLOPER *x);


/*f
***************************************************************************
** Convert from excel date to TDate
***************************************************************************
*/
TDate JpmcdsExcelDateToDate (double xlDate);                                   


/*f
***************************************************************************
** The Excel date system must be re-calculated for each function call.
** This function clears out the previous value of the date system.
***************************************************************************
*/
void JpmcdsClearExcelDateSystem (void);                                   


/*f
***************************************************************************
** Frees an Excel oper.
***************************************************************************
*/
void JpmcdsFreeExcelOper (XLOPER *x);


/*f
***************************************************************************
** Convert the output from a common routine into Excel format.
** This either results in an error output (if the output was NULL) 
** or an output which will be freed subsequently by xlAutoFree.
***************************************************************************
*/
LPXLOPER ExcelOutput (void *output);


/*f
***************************************************************************
** Returns the error output.
***************************************************************************
*/
LPXLOPER ErrorOutput (void);


/*f
***************************************************************************
** Returns a string giving the return status from a call to the 
** Excel function.
***************************************************************************
*/
char* JpmcdsExcelFormatReturnStatus (int xlret);

#endif /* XLUTILS_H */

