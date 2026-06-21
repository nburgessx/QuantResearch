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
#include "objects.h"
#include "version.h"
#include "cerror.h"
#include "tcurve.h"
#include "cdsone.h"
#include "convert.h"
#include "date_sup.h"
#include "zerocurve.h"
#include "cds.h"
#include "cxzerocurve.h"
#include "defaulted.h"

#include "bastypes.h"
#include "cdate.h"
#include "ldate.h"
#include "cgeneral.h"
#include "dtlist.h"
#include "stub.h"
#include "cashflow.h"
#include "macros.h"


/*
***************************************************************************
** Read array of doubles from Excel.
***************************************************************************
*/
int ReadDoubleArray(TBoolean mandatory, OPER *oper, size_t n, double **values, char *routine, char *param)
{
    size_t i;
    
    if (mandatory && oper->type == xltypeMissing)
    {
        JpmcdsErrMsg("%s: %s parameter is required\n", routine, param);
        return FAILURE;
    }

    if (oper->type == xltypeMulti)
    {
        if (oper->val.array.rows != n)
        {
            JpmcdsErrMsg("%s: %s - incorrect number of elements, expected %d\n", routine, param, n);
            return FAILURE;
        }
        
        *values = NEW_ARRAY(double, n);
        if (*values == NULL)
            return FAILURE;
        
        for (i = 0 ; i < n ; i++)
        {
            if (JpmcdsExcelGetDouble(oper->val.array.lparray + i, *values + i) != SUCCESS)
                return FAILURE;
        }
    }
    else
    {
        if (n != 1)
        {
            JpmcdsErrMsg("%s: %s - incorrect number of elements, expected %d\n", routine, param, n);
            return FAILURE;
        }
        
        *values = NEW_ARRAY(double, n);
        if (*values == NULL)
            return FAILURE;
        
        if (JpmcdsExcelGetDouble(oper, *values) != SUCCESS)
            return FAILURE;
    }

    return SUCCESS;
}


/*
***************************************************************************
** Read date from Excel, 
** input either as an explicit dates or a date interval.
***************************************************************************
*/
int ReadDateOrInterval(OPER *oper, int i, TDate valueDate, TDate* date, char *routine, char *param)
{
    if (JpmcdsExcelGetDate(oper, date) != SUCCESS)
    {
        TDateInterval  ivl;
        char          *s;

        if (JpmcdsExcelGetString(oper, &s) != SUCCESS)
            return FAILURE;

        if (JpmcdsStringToDate(s, date) != SUCCESS)
        {
            if (JpmcdsStringToDateInterval(s, routine, &ivl) != SUCCESS)
            {
                JpmcdsErrMsg ("%s: %s - invalid interval for element[%d].\n", routine, param, i);
                FREE(s);
                return FAILURE;
            }
            
            FREE(s);
    
            if (JpmcdsDateFwdThenAdjust(valueDate, &ivl, JPMCDS_BAD_DAY_NONE, "None", date) != SUCCESS)
            {
                JpmcdsErrMsg ("%s: %s - invalid interval for element[%d].\n", routine, param, i);
                return FAILURE;
            }
        }
    }
    
    return SUCCESS;
}


/*
***************************************************************************
** Read array of dates from Excel, 
** input either as explicit dates or as date intervals.
***************************************************************************
*/
int ReadDateOrIntervalArray(TBoolean mandatory, OPER *oper, size_t n, TDate valueDate, TDate **values, char *routine, char *param)
{
    size_t i;

    if (mandatory && oper->type == xltypeMissing)
    {
        JpmcdsErrMsg("%s: %s parameter is required\n", routine, param);
        return FAILURE;
    }

    if (oper->type == xltypeMulti)
    {
        if (oper->val.array.rows != n)
        {
            JpmcdsErrMsg("%s: %s - incorrect number of elements, expected %d\n", routine, param, n);
            return FAILURE;
        }
    
        *values = NEW_ARRAY(TDate, n);
        if (*values == NULL)
            return FAILURE;
    
        for (i = 0 ; i < n ; i++)
        {
            if (ReadDateOrInterval(oper->val.array.lparray + i, i+1, valueDate, *values+i, routine, param) != SUCCESS)
                return FAILURE;
        }
    }
    else
    {
        if (n != 1)
        {
            JpmcdsErrMsg("%s: %s - incorrect number of elements, expected %d\n", routine, param, n);
            return FAILURE;
        }
        
        *values = NEW_ARRAY(TDate, n);
        if (*values == NULL)
            return FAILURE;
        
        if (ReadDateOrInterval(oper, 1, valueDate, *values, routine, param) != SUCCESS)
            return FAILURE;
    }

    return SUCCESS;
}


/* convenience macro for getting Excel array size */
#define ARRAY_SIZE(oper, n) \
    if (oper->type == xltypeMissing) \
        n = 0; \
    else if (oper->type == xltypeMulti) \
        n = oper->val.array.rows; \
    else \
        n = 1
    

/* convenience macro for getting parameter values, checking if mandatory parameters are defined */
#define PARAM(mandatory, paramtype, oper, name) \
    if (mandatory && (oper->type == xltypeMissing || oper->type == xltypeNil)) \
    { \
        JpmcdsErrMsg("%s: %s parameter is required\n", routine, #name); \
        goto done; \
    } \
    if (oper->type == xltypeMulti) \
    { \
        JpmcdsErrMsg("%s: %s scalar value is required\n", routine, #name); \
        goto done; \
    } \
    if (JpmcdsExcelGet##paramtype##(oper, &##name##) != SUCCESS) \
        goto done
    

/* convenience macro for setting default parameter values */
#define SET_DEFAULT(variable, value) \
    if (variable == NULL) variable = strdup(value) \


/* convenience macro for returning value to Excel, reporting any error */
#define EXCEL_RETURN \
    if (status != SUCCESS) \
    { \
        JpmcdsErrMsg("%s: Failed!\n", routine); \
        JpmcdsFreeExcelOper(a0); \
    } \
    \
    return (status == SUCCESS) ? ExcelOutput(a0) : ErrorOutput()


/*
***************************************************************************
** Wrapper for Excel addin Version function.
***************************************************************************
*/
XLFUNC(LPXLOPER) CDS_Version(void)
{
    static char *routine = "CDS_Version";
    int          status = FAILURE;
    LPXLOPER     a0 = NULL;
    char         version[256];

    a0 = NEW(XLOPER);
    if (a0 == NULL)
        goto done;

    if (JpmcdsVersionString(version) != SUCCESS)
        goto done;

    if (JpmcdsExcelSetString(version, a0) != SUCCESS)
        goto done;

    status = SUCCESS;

done:
    EXCEL_RETURN;
}


/*
***************************************************************************
** Wrapper for Excel addin ErrorLogStatus function.
***************************************************************************
*/
XLFUNC(LPXLOPER) CDS_ErrorLogStatus(void)
{
    static char *routine = "CDS_ErrorLogStatus";
    int          status = FAILURE;
    LPXLOPER     a0 = NULL;

    a0 = NEW(XLOPER);
    if (a0 == NULL)
        goto done;

    a0->xltype  = xltypeNum;
    a0->val.num = (double) JpmcdsErrMsgStatus();
    status = SUCCESS;

done:
    EXCEL_RETURN;
}


/*
***************************************************************************
** Wrapper for Excel addin ErrorLogContents function.
***************************************************************************
*/
XLFUNC(LPXLOPER) CDS_ErrorLogContents(void)
{
    int       status = FAILURE;
    LPXLOPER  a0 = NULL;
    char    **lines;
    int       rows = 0;
    int       i;
    
    a0 = NEW(XLOPER);
    if (a0 == NULL)
        goto done;

    lines = JpmcdsErrGetMsgRecord();
    if (lines == NULL)
    {
        status = JpmcdsExcelSetString("(no log contents)", a0);
        goto done;
    }

    while (lines[rows] != NULL)
        rows++;

    a0->xltype = xltypeMulti;
    a0->val.array.rows = rows;
    a0->val.array.columns = 1;
    a0->val.array.lparray = NEW_ARRAY(XLOPER, rows);
    if (a0->val.array.lparray == NULL)
        goto done;
            
    for (i = 0; i < rows; i++)
        (void) JpmcdsExcelSetString(lines[i], a0->val.array.lparray + i);

    status = SUCCESS;

done:
    if (status != SUCCESS)
        JpmcdsFreeExcelOper(a0);
    return (status == SUCCESS) ? ExcelOutput(a0) : ErrorOutput();
}


/*
***************************************************************************
** Wrapper for Excel addin ErrorLogFilename function.
***************************************************************************
*/
XLFUNC(LPXLOPER) CDS_ErrorLogFilename(void)
{
    static char *routine = "CDS_ErrorLogFilename";
    int          status = FAILURE;
    LPXLOPER     a0 = NULL;
    char*        filename;

    a0 = NEW(XLOPER);
    if (a0 == NULL)
        goto done;

    filename = JpmcdsErrMsgGetFileName();
    if (filename == NULL)
        goto done;

    if (JpmcdsExcelSetString(filename, a0) != SUCCESS)
        goto done;

    status = SUCCESS;

done:
    EXCEL_RETURN;
}


/*
***************************************************************************
** Wrapper for Excel addin SetErrorLogFilename function.
***************************************************************************
*/
XLFUNC(LPXLOPER) CDS_SetErrorLogFilename(OPER *a1,OPER *a2)
{
    static char *routine = "CDS_SetErrorLogFilename";
    int          status = FAILURE;
    LPXLOPER     a0 = NULL;
    char        *filename = NULL;
    long         append;

    a0 = NEW(XLOPER);
    if (a0 == NULL)
        goto done;

    PARAM(TRUE, String, a1, filename);
    PARAM(TRUE, Long  , a2, append);

    status = JpmcdsErrMsgFileName(filename, (TBoolean) (append));
    a0->xltype  = xltypeNum;
    a0->val.num = 1.0;
    status = SUCCESS;

done:
    FREE(filename);
    EXCEL_RETURN;
}


/*
***************************************************************************
** Wrapper for Excel addin SetErrorLogStatus function.
***************************************************************************
*/
XLFUNC(LPXLOPER) CDS_SetErrorLogStatus(OPER *a1)
{
    static char *routine = "CDS_SetErrorLogStatus";
    int          status = FAILURE;
    LPXLOPER     a0 = NULL;
    long         state;

    a0 = NEW(XLOPER);
    if (a0 == NULL)
        goto done;

    PARAM(TRUE, Long, a1, state);

    if (state != 0)
    {
        if (JpmcdsErrMsgEnableRecord(20, 128) != SUCCESS) /* ie. 20 lines, each of max length 128 */
            goto done;
        JpmcdsErrMsgOn();
    }
    else 
    {
        if (JpmcdsErrMsgDisableRecord() != SUCCESS)
            goto done;
        JpmcdsErrMsgOff();
    }

    a0->xltype  = xltypeNum;
    a0->val.num = (state ? 1.0 : 0.0);
    status = SUCCESS;

done:
    EXCEL_RETURN;   
}


/*
***************************************************************************
** Wrapper for Excel addin LoadHolidays function.
***************************************************************************
*/
XLFUNC(LPXLOPER) CDS_LoadHolidays(OPER *a1,OPER *a2)
{
    static char *routine = "CDS_LoadHolidays";
    int          status = FAILURE;
    LPXLOPER     a0 = NULL;
    char        *name = NULL;
    char        *filename = NULL;

    a0 = NEW(XLOPER);
    if (a0 == NULL)
        goto done;

    PARAM(TRUE, String, a1, name);
    PARAM(TRUE, String, a2, filename);

    status = JpmcdsHolidayLoadFromDisk(name, filename);

    a0->xltype  = xltypeNum;
    a0->val.num = 1.0;
    status = SUCCESS;

done:
    FREE(name);
    FREE(filename);
    EXCEL_RETURN;
}


/*
***************************************************************************
** Wrapper for Excel addin IRZeroCurveBuild function.
***************************************************************************
*/
XLFUNC(LPXLOPER) CDS_IRZeroCurveBuild(OPER *a1,OPER *a2,OPER *a3,OPER *a4,OPER *a5,
        OPER *a6,OPER *a7,OPER *a8,OPER *a9,OPER *a10,OPER *a11,OPER *a12)
{
    static char  *routine = "CDS_IRZeroCurveBuild";
    int           status = FAILURE;
    LPXLOPER      a0 = NULL;
    TDate         valueDate;
    char         *types = NULL;
    TDate        *dates = NULL;
    double       *rates = NULL;
    char         *mmDcc = NULL;
    char         *fixedIvl = NULL;
    char         *floatIvl = NULL;
    char         *fixedDcc = NULL;
    char         *floatDcc = NULL;
    char         *badDayConv = NULL;
    char         *holidays = NULL;
    TCurve       *zeroCurve = NULL;
    char         *name = NULL;
    char         *handle = NULL;

    int           i;
    int           nInstr;
    long          mmDCC;
    TDateInterval fixedIVL;
    TDateInterval floatIVL;
    long          fixedDCC;
    long          floatDCC;
    double        fixedFreq;
    double        floatFreq;

    int           dtSuccess;
    TDateInterval tmp;
    TDateAdjIntvl busday;
    TDate         baseDate;

    a0 = NEW(XLOPER);
    if (a0 == NULL)
        goto done;

    JpmcdsClearExcelDateSystem();
    ARRAY_SIZE(a2, nInstr);
    PARAM(TRUE , Date, a1, valueDate);

    /* instrument types */
    if (a2->type == xltypeMissing)
    {
        JpmcdsErrMsg("%s: Types parameter is required\n", routine);
        goto done;
    }

    types = NEW_ARRAY(char, nInstr);
    if (types == NULL)
        goto done;

    for (i = 0 ; i < nInstr ; i++)
    {
        char *type;
        if (JpmcdsExcelGetString(a2->val.array.lparray + i, &type) != SUCCESS)
            goto done;
        types[i] = type[0];
        FREE(type);
    }

    if (ReadDateOrIntervalArray(TRUE, a3, nInstr, valueDate, &dates, routine, "End Dates") != SUCCESS)
        goto done;

    if (ReadDoubleArray(TRUE, a4, nInstr, &rates, routine, "Rates") != SUCCESS)
        goto done;

    PARAM(FALSE, String, a5 , mmDcc);
    PARAM(TRUE , String, a6 , fixedIvl);
    PARAM(FALSE, String, a7 , floatIvl);
    PARAM(TRUE , String, a8 , fixedDcc);
    PARAM(FALSE, String, a9 , floatDcc);
    PARAM(FALSE, String, a10, badDayConv);
    PARAM(TRUE , String, a11, holidays);
    PARAM(TRUE , String, a12, name);

    if (mmDcc == NULL) mmDcc = "ACT/360";
    if (JpmcdsStringToDayCountConv(mmDcc, &mmDCC) != SUCCESS)
        goto done;

    if (JpmcdsStringToDateInterval(fixedIvl, routine, &fixedIVL) != SUCCESS)
        goto done;

    if (floatIvl == NULL)
        floatIvl = fixedIvl;

    if (JpmcdsStringToDateInterval(floatIvl, routine, &floatIVL) != SUCCESS)
        goto done;

    if (JpmcdsStringToDayCountConv(fixedDcc, &fixedDCC) != SUCCESS)
        goto done;

    if (floatDcc == NULL)
        floatDcc = fixedDcc;

    if (JpmcdsStringToDayCountConv(floatDcc, &floatDCC) != SUCCESS)
        goto done;

    SET_DEFAULT(badDayConv, "N");
    if (JpmcdsBadDayConvValid(routine, badDayConv[0]) != SUCCESS)
        goto done;

    if (JpmcdsDateIntervalToFreq(&fixedIVL, &fixedFreq) != SUCCESS)
        goto done;

    if (JpmcdsDateIntervalToFreq(&floatIVL, &floatFreq) != SUCCESS)
        goto done;

    /** adjust cash rates dates to business day */
    baseDate = valueDate;
    for(i = 0; i < nInstr; i++)
    {
        if(types[i] == 'M')
        {
            if(dates[i] - baseDate <= 3)
            {
                /* for business days */
                tmp.flag = 0;
                tmp.prd = dates[i] - baseDate;
                tmp.prd_typ = 'D';
                busday.holidayFile = holidays;
                busday.isBusDays = TRUE;
                busday.badDayConv = JPMCDS_BAD_DAY_FOLLOW;
                busday.interval = tmp;
                /* adjust to business day */
                dtSuccess = JpmcdsDtFwdAdj(baseDate, &busday, dates+i);
            }
            else if(dates[i] - baseDate <= 21)
            {
                /* for less than or equal to 3 weeks */
                /* adjust to business day */
                dtSuccess = JpmcdsBusinessDay(dates[i],'F', holidays, dates+i);
            }
            else
            {
                /* adjust to business day */
                dtSuccess = JpmcdsBusinessDay(dates[i],'M', holidays, dates+i);
            }
        }
    }

    zeroCurve = JpmcdsBuildIRZeroCurve(
            valueDate,
            types,
            dates,
            rates,
            nInstr,
            mmDCC,
            (long) fixedFreq,
            (long) floatFreq,
            fixedDCC,
            floatDCC,
            badDayConv[0],
            holidays);

    handle = StoreObject(name, zeroCurve);
    if (handle == NULL)
        goto done;

    if (JpmcdsExcelSetString(handle, a0) != SUCCESS)
        return NULL;

    status = SUCCESS;

done:
    FREE(types);
    FREE(dates);
    FREE(rates);
    FREE(mmDcc);
    FREE(fixedDcc);
    FREE(floatDcc);
    FREE(fixedIvl);
    FREE(floatIvl);
    FREE(badDayConv);
    FREE(holidays);
    FREE(name);
    FREE(handle);
    EXCEL_RETURN;
}


/*
***************************************************************************
** Wrapper for Excel addin IRZeroCurveMake function.
***************************************************************************
*/
XLFUNC(LPXLOPER) CDS_IRZeroCurveMake(OPER *a1,OPER *a2,OPER *a3,OPER *a4,OPER *a5,OPER *a6,OPER *a7)
{
    static char  *routine = "CDS_IRZeroCurveMake";
    int           status = FAILURE;
    LPXLOPER      a0 = NULL;
    TDate         baseDate;
    TDate        *dates = NULL;
    double       *rates = NULL;
    long          basis;
    char         *zcDcc = NULL;
    TCurve       *zeroCurve = NULL;
    char         *name = NULL;
    char         *handle = NULL;

    int           n;
    long          dcc;

    a0 = NEW(XLOPER);
    if (a0 == NULL)
        goto done;

    JpmcdsClearExcelDateSystem();
    ARRAY_SIZE(a2, n);
    PARAM(TRUE , Date, a1, baseDate);

    if (ReadDateOrIntervalArray(TRUE, a2, n, baseDate, &dates, routine, "Dates") != SUCCESS)
        goto done;

    if (ReadDoubleArray(TRUE, a3, n, &rates, routine, "Rates") != SUCCESS)
        goto done;

    PARAM(TRUE , Long  , a4 , basis);
    PARAM(FALSE, String, a5 , zcDcc);
    PARAM(TRUE , String, a6 , name);

    SET_DEFAULT(zcDcc, "ACT/365F");
    if (JpmcdsStringToDayCountConv(zcDcc, &dcc) != SUCCESS)
        goto done;

    zeroCurve = JpmcdsMakeTCurve(baseDate, dates, rates, n, basis, dcc);

    handle = StoreObject(name, zeroCurve);
    if (handle == NULL)
        goto done;

    if (JpmcdsExcelSetString(handle, a0) != SUCCESS)
        return NULL;

    status = SUCCESS;

done:
    FREE(dates);
    FREE(rates);
    FREE(zcDcc);
    FREE(name);
    FREE(handle);
    EXCEL_RETURN;
}


/*
***************************************************************************
** Wrapper for Excel addin CleanSpreadCurveBuild function.
***************************************************************************
*/
XLFUNC(LPXLOPER) CDS_CleanSpreadCurveBuild(OPER *a1,OPER *a2,OPER *a3,OPER *a4,
OPER *a5,OPER *a6,OPER *a7,OPER *a8,OPER *a9,OPER *a10,OPER *a11,OPER *a12,
OPER *a13,OPER *a14,OPER *a15,OPER *a16)
{
    static char   *routine = "CDS_CleanSpreadCurveBuild";
    int            status = FAILURE;
    LPXLOPER       a0 = NULL;
    TDate          today;               /* 1 */
    TDate          startDate;
    TDate          stepinDate;
    TDate          cashSettleDate;
    TDate         *dates = NULL;        /* 5 */
    double        *rates = NULL;
    long          *includes = NULL;
    long           payAccOnDefault;
    char          *couponInterval = NULL;
    char          *stubType = NULL;    /* 10 */
    char          *paymentDcc = NULL;
    char          *badDayConv = NULL;
    char          *holidays = NULL;
    double         recoveryRate;
    char          *name = NULL;
    TCurve        *discCurve = NULL;
    TCurve        *spreadCurve = NULL;
    char          *handle1 = NULL;
    char          *handle2 = NULL;

    int            i;
    int            n;
    TStubMethod    stub;
    long           dcc;
    TDateInterval  ivl;

    a0 = NEW(XLOPER);
    if (a0 == NULL)
        goto done;

    JpmcdsClearExcelDateSystem();
    ARRAY_SIZE(a5, n);
    PARAM(TRUE , Date  , a1 , today);
    PARAM(TRUE , Date  , a2 , startDate);
    PARAM(TRUE , Date  , a3 , stepinDate);
    PARAM(TRUE , Date  , a4 , cashSettleDate);

    if (ReadDateOrIntervalArray(TRUE, a5, n, startDate, &dates, routine, "End Dates") != SUCCESS)
        goto done;

    if (ReadDoubleArray(TRUE, a6, n, &rates, routine, "Rates") != SUCCESS)
        goto done;

    /* include flags - NB. missing entry = include everything */
    if (a7->type != xltypeMissing && a7->val.array.rows != 0)
    {
        includes = NEW_ARRAY(long, n);
        if (includes == NULL)
            goto done;
    
        if (a7->val.array.rows != n)
        {
            JpmcdsErrMsg("%s: Include Flags - incorrect number of elements, expected %d\n", routine, n);
            goto done;
        }
    
        for (i = 0 ; i < n ; i++)
        {
            if (JpmcdsExcelGetLong(a7->val.array.lparray + i, includes + i) != SUCCESS)
                goto done;      
        }
    }

    PARAM(TRUE , Long  , a8 , payAccOnDefault);
    PARAM(TRUE , String, a9 , couponInterval);
    PARAM(FALSE, String, a10,  stubType);
    PARAM(FALSE, String, a11, paymentDcc);
    PARAM(FALSE, String, a12, badDayConv);
    PARAM(TRUE , String, a13, holidays);
    PARAM(TRUE , String, a14, handle1);
    PARAM(TRUE , Double, a15, recoveryRate);
    PARAM(TRUE , String, a16, name);

    SET_DEFAULT(stubType, "f/s");
    if (JpmcdsStringToStubMethod(stubType, &stub) != SUCCESS)
        goto done;

    SET_DEFAULT(paymentDcc, "ACT/360");
    if (JpmcdsStringToDayCountConv(paymentDcc, &dcc) != SUCCESS)
        goto done;

    SET_DEFAULT(badDayConv, "N");
    if (JpmcdsBadDayConvValid(routine, badDayConv[0]) != SUCCESS)
        goto done;

    if (JpmcdsStringToDateInterval(couponInterval, routine, &ivl) != SUCCESS)
        goto done;

    discCurve = (TCurve*) RetrieveObject(handle1);
    if (discCurve == NULL)
        goto done;
    
     spreadCurve = JpmcdsCleanSpreadCurve(today,
                                          discCurve,
                                          startDate,
                                          stepinDate,
                                          cashSettleDate,
                                          n,
                                          dates,
                                          rates,
                                          includes,
                                          recoveryRate,
                                          (TBoolean)(payAccOnDefault),
                                          &ivl,
                                          dcc,
                                          &stub,
                                          badDayConv[0],
                                          holidays);

    handle2 = StoreObject(name, spreadCurve);
    if (handle2 == NULL)
        goto done;

    if (JpmcdsExcelSetString(handle2, a0) != SUCCESS)
        return NULL;

    status = SUCCESS;

done:
    FREE(dates);
    FREE(rates);
    FREE(includes);
    FREE(couponInterval);
    FREE(stubType);
    FREE(paymentDcc);
    FREE(badDayConv);
    FREE(holidays);
    FREE(handle1);
    FREE(handle2);
    FREE(name);
    EXCEL_RETURN;
}


/*
***************************************************************************
** Wrapper for Excel addin DiscountFactor function.
***************************************************************************
*/
XLFUNC(LPXLOPER) CDS_DiscountFactor(OPER *a1,OPER *a2)
{
    static char   *routine = "CDS_DiscountFactor";
    int            status = FAILURE;
    LPXLOPER       a0 = NULL;
    char          *handle = NULL;
    TDate          date;
    TCurve        *curve;
    double         result;

    a0 = NEW(XLOPER);
    if (a0 == NULL)
        goto done;

    JpmcdsClearExcelDateSystem();
    PARAM(TRUE , String, a1, handle);
    PARAM(TRUE , Date  , a2, date);

    curve = (TCurve*) RetrieveObject(handle);
    if (curve == NULL)
        goto done;
    
    result = JpmcdsZeroPrice(curve, date);

    a0->xltype  = xltypeNum;
    a0->val.num = result;
    status = SUCCESS;

done:
    FREE(handle);
    EXCEL_RETURN;
}


/*
***************************************************************************
** Wrapper for Excel addin DatesAndRates function.
***************************************************************************
*/
XLFUNC(LPXLOPER) CDS_DatesAndRates(OPER *a1)
{
    static char   *routine = "CDS_DatesAndRates";
    int            status = FAILURE;
    LPXLOPER       a0 = NULL;
    char          *handle = NULL;
    TCurve        *curve;
    int            i;

    a0 = NEW(XLOPER);
    if (a0 == NULL)
        goto done;

    JpmcdsClearExcelDateSystem();
    PARAM(TRUE, String, a1 , handle);

    curve = (TCurve*) RetrieveObject(handle);
    if (curve == NULL)
        goto done;    
    
    /* dates and rates to a 2 column Excel array */
    a0->xltype = xltypeMulti;
    a0->val.array.rows = curve->fNumItems;
    a0->val.array.columns = 2;
    a0->val.array.lparray = NEW_ARRAY(XLOPER, 2 * curve->fNumItems);
    if (a0->val.array.lparray == NULL)
        goto done;

    for (i = 0; i < curve->fNumItems; i++)
    {
        a0->val.array.lparray[2*i].xltype  = xltypeNum;
        JpmcdsExcelSetDate(curve->fArray[i].fDate, a0->val.array.lparray + 2*i);
        a0->val.array.lparray[2*i+1].xltype  = xltypeNum;
        a0->val.array.lparray[2*i+1].val.num  = curve->fArray[i].fRate;
    }

    status = SUCCESS;

done:
    FREE(handle);
    EXCEL_RETURN;
}


/*
***************************************************************************
** Wrapper for Excel addin UpfrontFlat function.
***************************************************************************
*/
XLFUNC(LPXLOPER) CDS_UpfrontFlat(OPER *a1,OPER *a2,OPER *a3,OPER *a4,OPER *a5,OPER *a6,OPER *a7,
                                 OPER *a8,OPER *a9,OPER *a10,OPER *a11,OPER *a12,OPER *a13,
                                 OPER *a14,OPER *a15,OPER *a16,OPER *a17)
{
    static char   *routine = "CDS_UpfrontFlat";
    int            status = FAILURE;
    LPXLOPER       a0 = NULL;
    TDate          today;
    TDate          valueDate;
    TDate          benchmarkStartDate;
    TDate          startDate;
    TDate          stepinDate;
    TDate          endDate;
    double         couponRate;
    long           payAccOnDefault;
    char          *couponInterval = NULL;
    char          *stubType = NULL;
    char          *paymentDcc = NULL;
    char          *badDayConv = NULL;
    char          *holidays = NULL;
    char          *handle = NULL;
    double         parSpread;
    double         recoveryRate;
    long           isPriceClean;
    double         result;

    TStubMethod    stub;
    long           dcc;
    TDateInterval  ivl;
    TCurve        *curve;

    a0 = NEW(XLOPER);
    if (a0 == NULL)
        goto done;

    JpmcdsClearExcelDateSystem();
    PARAM(TRUE , Date  , a1 , today);
    PARAM(TRUE , Date  , a2 , valueDate);
    PARAM(TRUE , Date  , a3 , benchmarkStartDate);
    PARAM(TRUE , Date  , a4 , stepinDate);
    PARAM(TRUE , Date  , a5 , startDate);
    PARAM(TRUE , Date  , a6 , endDate);
    PARAM(TRUE , Double, a7 , couponRate);
    PARAM(TRUE , Long  , a8 , payAccOnDefault);
    PARAM(TRUE , String, a9 , couponInterval);
    PARAM(FALSE, String, a10, stubType);
    PARAM(FALSE, String, a11, paymentDcc);
    PARAM(FALSE, String, a12, badDayConv);
    PARAM(TRUE , String, a13, holidays);
    PARAM(TRUE , String, a14, handle);
    PARAM(TRUE , Double, a15, parSpread);
    PARAM(TRUE , Double, a16, recoveryRate);
    PARAM(TRUE , Long  , a17, isPriceClean);

    SET_DEFAULT(stubType, "f/s");
    if (JpmcdsStringToStubMethod(stubType, &stub) != SUCCESS)
        goto done;

    SET_DEFAULT(paymentDcc, "ACT/360");
    if (JpmcdsStringToDayCountConv(paymentDcc, &dcc) != SUCCESS)
        goto done;

    SET_DEFAULT(badDayConv, "N");
    if (JpmcdsBadDayConvValid(routine, badDayConv[0]) != SUCCESS)
        goto done;

    if (JpmcdsStringToDateInterval(couponInterval, routine, &ivl) != SUCCESS)
        goto done;

    curve = (TCurve*) RetrieveObject(handle);
    if (curve == NULL)
        goto done;    

    if (JpmcdsCdsoneUpfrontCharge(today,
                                  valueDate,
                                  benchmarkStartDate,
                                  stepinDate,
                                  startDate,
                                  endDate,
                                  couponRate,
                                  (TBoolean)(payAccOnDefault),
                                  &ivl,
                                  &stub,
                                  dcc,
                                  badDayConv[0],
                                  holidays,
                                  curve,
                                  parSpread,
                                  recoveryRate,
                                  (TBoolean)(isPriceClean),
                                  &result) != SUCCESS) goto done;

    a0->xltype  = xltypeNum;
    a0->val.num = result;
    status = SUCCESS;

done:
    FREE(couponInterval);
    FREE(stubType);
    FREE(paymentDcc);
    FREE(badDayConv);
    FREE(holidays);
    FREE(handle);
    EXCEL_RETURN;
}


/*
***************************************************************************
** Wrapper for Excel addin ParSpreadFlat function.
***************************************************************************
*/
XLFUNC(LPXLOPER) CDS_ParSpreadFlat(OPER *a1,OPER *a2,OPER *a3,OPER *a4,OPER *a5,OPER *a6,OPER *a7,
                                   OPER *a8,OPER *a9,OPER *a10,OPER *a11,OPER *a12,OPER *a13,OPER *a14,OPER *a15,OPER *a16,
                                   OPER *a17)
{
    static char   *routine = "CDS_ParSpreadFlat";
    int            status = FAILURE;
    LPXLOPER       a0 = NULL;
    TDate          today;
    TDate          valueDate;
    TDate          benchmarkStartDate;
    TDate          stepinDate;
    TDate          startDate;
    TDate          endDate;
    double         couponRate;
    long           payAccOnDefault;
    char          *couponInterval = NULL;
    char          *stubType = NULL;
    char          *paymentDcc = NULL;
    char          *badDayConv = NULL;
    char          *holidays = NULL;
    char          *handle = NULL;
    double         upfrontCharge;
    double         recoveryRate;
    long           isPriceClean;
    double         result;

    TStubMethod    stub;
    long           dcc;
    TDateInterval  ivl;
    TCurve        *curve;

    a0 = NEW(XLOPER);
    if (a0 == NULL)
        goto done;

    JpmcdsClearExcelDateSystem();
    PARAM(TRUE , Date  , a1 , today);
    PARAM(TRUE , Date  , a2 , valueDate);
    PARAM(TRUE , Date  , a3 , benchmarkStartDate);
    PARAM(TRUE , Date  , a4 , stepinDate);
    PARAM(TRUE , Date  , a5 , startDate);
    PARAM(TRUE , Date  , a6 , endDate);
    PARAM(TRUE , Double, a7 , couponRate);
    PARAM(TRUE , Long  , a8 , payAccOnDefault);
    PARAM(TRUE , String, a9 , couponInterval);
    PARAM(FALSE, String, a10, stubType);
    PARAM(FALSE, String, a11, paymentDcc);
    PARAM(FALSE, String, a12, badDayConv);
    PARAM(TRUE , String, a13, holidays);
    PARAM(TRUE , String, a14, handle);
    PARAM(TRUE , Double, a15, upfrontCharge);
    PARAM(TRUE , Double, a16, recoveryRate);
    PARAM(TRUE , Long  , a17, isPriceClean);

    SET_DEFAULT(stubType, "f/s");
    if (JpmcdsStringToStubMethod(stubType, &stub) != SUCCESS)
        goto done;

    SET_DEFAULT(paymentDcc, "ACT/360");
    if (JpmcdsStringToDayCountConv(paymentDcc, &dcc) != SUCCESS)
        goto done;

    SET_DEFAULT(badDayConv, "N");
    if (JpmcdsBadDayConvValid(routine, badDayConv[0]) != SUCCESS)
        goto done;

    if (JpmcdsStringToDateInterval(couponInterval, routine, &ivl) != SUCCESS)
        goto done;

    curve = (TCurve*) RetrieveObject(handle);
    if (curve == NULL)
        goto done;    

    if (JpmcdsCdsoneSpread(today,
                           valueDate,
                           benchmarkStartDate,
                           stepinDate,
                           startDate,
                           endDate,
                           couponRate,
                           (TBoolean)(payAccOnDefault),
                           &ivl,
                           &stub,
                           dcc,
                           badDayConv[0],
                           holidays,
                           curve,
                           upfrontCharge,
                           recoveryRate,
                           (TBoolean)(isPriceClean),
                           &result) != SUCCESS) goto done;
    
    a0->xltype  = xltypeNum;
    a0->val.num = result;
    status = SUCCESS;

done:
    FREE(couponInterval);
    FREE(stubType);
    FREE(paymentDcc);
    FREE(badDayConv);
    FREE(holidays);
    FREE(handle);
    EXCEL_RETURN;
}


/*
***************************************************************************
** Wrapper for Excel addin CdsPrice function.
***************************************************************************
*/
XLFUNC(LPXLOPER) CDS_CdsPrice(OPER *a1,OPER *a2,OPER *a3,OPER *a4,OPER *a5,OPER *a6,OPER *a7,OPER *a8,
                              OPER *a9,OPER *a10,OPER *a11,OPER *a12,OPER *a13,OPER *a14,OPER *a15,
                              OPER *a16)
{
    static char  *routine = "CDS_CdsPrice";
    int           status = FAILURE;
    LPXLOPER      a0 = NULL;
    TDate         today;
    TDate         valueDate;
    TDate         stepinDate;
    TDate         startDate;
    TDate         endDate;
    double        couponRate;
    long          payAccOnDefault;
    char         *couponInterval = NULL;
    char         *stubType = NULL;
    char         *paymentDcc = NULL;
    char         *badDayConv = NULL;
    char         *holidays = NULL;
    char         *handle1 = NULL;
    char         *handle2 = NULL;
    double        recoveryRate;
    long          isPriceClean;
    double        result;

    TStubMethod    stub;
    long           dcc;
    TDateInterval  ivl;
    TCurve        *discCurve;
    TCurve        *spreadCurve;

    a0 = NEW(XLOPER);
    if (a0 == NULL)
        goto done;

    JpmcdsClearExcelDateSystem();
    PARAM(TRUE , Date  , a1 , today);
    PARAM(TRUE , Date  , a2 , valueDate);
    PARAM(TRUE , Date  , a3 , stepinDate);
    PARAM(TRUE , Date  , a4 , startDate);
    PARAM(TRUE , Date  , a5 , endDate);
    PARAM(TRUE , Double, a6 , couponRate);
    PARAM(TRUE , Long  , a7 , payAccOnDefault);
    PARAM(TRUE , String, a8 , couponInterval);
    PARAM(FALSE, String, a9 , stubType);
    PARAM(FALSE, String, a10,  paymentDcc);
    PARAM(FALSE, String, a11, badDayConv);
    PARAM(TRUE , String, a12, holidays);
    PARAM(TRUE , String, a13, handle1);
    PARAM(TRUE , String, a14, handle2);
    PARAM(TRUE , Double, a15, recoveryRate);
    PARAM(TRUE , Long  , a16, isPriceClean);
    
    SET_DEFAULT(stubType, "f/s");
    if (JpmcdsStringToStubMethod(stubType, &stub) != SUCCESS)
        goto done;

    SET_DEFAULT(paymentDcc, "ACT/360");
    if (JpmcdsStringToDayCountConv(paymentDcc, &dcc) != SUCCESS)
        goto done;

    SET_DEFAULT(badDayConv, "N");
    if (JpmcdsBadDayConvValid(routine, badDayConv[0]) != SUCCESS)
        goto done;

    if (JpmcdsStringToDateInterval(couponInterval, routine, &ivl) != SUCCESS)
        goto done;

    discCurve = (TCurve*) RetrieveObject(handle1);
    if (discCurve == NULL)
        goto done;    

    spreadCurve = (TCurve*) RetrieveObject(handle2);
    if (spreadCurve == NULL)
        goto done;    

    if (JpmcdsCdsPrice(today,
                       valueDate,
                       stepinDate,
                       startDate,
                       endDate,
                       couponRate,
                       (TBoolean)(payAccOnDefault),
                       &ivl,
                       &stub,
                       dcc,
                       badDayConv[0],
                       holidays,
                       discCurve,
                       spreadCurve,
                       recoveryRate,
                       (TBoolean)(isPriceClean),
                       &result) != SUCCESS) goto done;

    a0->xltype  = xltypeNum;
    a0->val.num = result;
    status = SUCCESS;

done:
    FREE(couponInterval);
    FREE(stubType);
    FREE(paymentDcc);
    FREE(badDayConv);
    FREE(holidays);
    FREE(handle1);
    FREE(handle2);
    EXCEL_RETURN;
}


/*
***************************************************************************
** Wrapper for Excel addin ParSpreads function.
***************************************************************************
*/
XLFUNC(LPXLOPER) CDS_ParSpreads(OPER *a1,OPER *a2,OPER *a3,OPER *a4,OPER *a5,OPER *a6,OPER *a7,
                                OPER *a8,OPER *a9,OPER *a10,OPER *a11,OPER *a12,OPER *a13)
{
    static char   *routine = "CDS_ParSpreads";
    int            status = FAILURE;
    LPXLOPER       a0 = NULL;
    TDate          today;
    TDate          stepinDate;
    TDate          startDate;
    TDate         *dates = NULL;
    long           payAccOnDefault;
    char          *couponInterval = NULL;
    char          *stubType = NULL;
    char          *paymentDcc = NULL;
    char          *badDayConv = NULL;
    char          *holidays = NULL;
    double         recoveryRate;
    double        *results = NULL;
    char          *handle1 = NULL;
    char          *handle2 = NULL;

    int            n, i;
    TStubMethod    stub;
    long           dcc;
    TDateInterval  ivl;
    TCurve        *discCurve;
    TCurve        *spreadCurve;

    a0 = NEW(XLOPER);
    if (a0 == NULL)
        goto done;

    JpmcdsClearExcelDateSystem();
    ARRAY_SIZE(a4, n);
    PARAM(TRUE , Date  , a1 , today);
    PARAM(TRUE , Date  , a2 , stepinDate);
    PARAM(TRUE , Date  , a3 , startDate);
    if (ReadDateOrIntervalArray(TRUE, a4, n, stepinDate, &dates, routine, "End Dates") != SUCCESS)
        goto done;
    PARAM(TRUE , Long  , a5 , payAccOnDefault);
    PARAM(TRUE , String, a6 , couponInterval);
    PARAM(FALSE, String, a7 , stubType);
    PARAM(FALSE, String, a8 , paymentDcc);
    PARAM(FALSE, String, a9, badDayConv);
    PARAM(TRUE , String, a10, holidays);
    PARAM(TRUE , String, a11, handle1);
    PARAM(TRUE , String, a12, handle2);
    PARAM(TRUE , Double, a13, recoveryRate);
    
    SET_DEFAULT(stubType, "f/s");
    if (JpmcdsStringToStubMethod(stubType, &stub) != SUCCESS)
        goto done;

    SET_DEFAULT(paymentDcc, "ACT/360");
    if (JpmcdsStringToDayCountConv(paymentDcc, &dcc) != SUCCESS)
        goto done;

    SET_DEFAULT(badDayConv, "N");
    if (JpmcdsBadDayConvValid(routine, badDayConv[0]) != SUCCESS)
        goto done;

    if (JpmcdsStringToDateInterval(couponInterval, routine, &ivl) != SUCCESS)
        goto done;

    discCurve = (TCurve*) RetrieveObject(handle1);
    if (discCurve == NULL)
        goto done;    

    spreadCurve = (TCurve*) RetrieveObject(handle2);
    if (spreadCurve == NULL)
        goto done;

    results = NEW_ARRAY(double, n);
    if(results == NULL)
        goto done;

    if (JpmcdsCdsParSpreads(today,
                            stepinDate,
                            startDate,
                            n,
                            dates,
                            (TBoolean)(payAccOnDefault),
                            &ivl,
                            &stub,
                            dcc,
                            badDayConv[0],
                            holidays,
                            discCurve,
                            spreadCurve,
                            recoveryRate,
                            results) != SUCCESS)
        goto done;

    /* results is an 1 column Excel array */
    a0->xltype = xltypeMulti;
    a0->val.array.rows = n;
    a0->val.array.columns = 1;
    a0->val.array.lparray = NEW_ARRAY(XLOPER, n);
    if (a0->val.array.lparray == NULL)
        goto done;

    for (i = 0; i < n; i++)
    {
        a0->val.array.lparray[i].xltype  = xltypeNum;
        a0->val.array.lparray[i].val.num  = results[i];
    }

    status = SUCCESS;

done:
    FREE(dates);
    FREE(stubType);
    FREE(couponInterval);
    FREE(paymentDcc);
    FREE(badDayConv);
    FREE(holidays);
    FREE(handle1);
    FREE(handle2);
    FREE(results);
    EXCEL_RETURN;
}


/*
***************************************************************************
** Wrapper for Excel addin FeeLegFlows function.
***************************************************************************
*/
XLFUNC(LPXLOPER) CDS_FeeLegFlows(OPER *a1,OPER *a2,OPER *a3,OPER *a4,OPER *a5,
                                 OPER *a6,OPER *a7,OPER *a8,OPER *a9)
{
    static char   *routine = "CDS_FeeLegFlows";
    int            status = FAILURE;
    LPXLOPER       a0 = NULL;
    TDate          startDate;
    TDate          endDate;
    double         rate;
    double         notional;
    char          *couponInterval = NULL;
    char          *stubType = NULL;
    char          *paymentDcc = NULL;
    char          *badDayConv = NULL;
    char          *holidays = NULL;
    TCashFlowList *cfl = NULL;

    int            i;
    TStubMethod    stub;
    long           dcc;
    TDateInterval  ivl;

    a0 = NEW(XLOPER);
    if (a0 == NULL)
        goto done;

    JpmcdsClearExcelDateSystem();
    PARAM(TRUE , Date  , a1 , startDate);
    if (ReadDateOrInterval(a2, TRUE, startDate, &endDate, routine, "End Date") != SUCCESS)
        goto done;
    PARAM(TRUE , Double, a3 , rate);
    PARAM(TRUE , Double, a4 , notional);
    PARAM(TRUE , String, a5 , couponInterval);
    PARAM(FALSE, String, a6 , stubType);
    PARAM(FALSE, String, a7 , paymentDcc);
    PARAM(FALSE, String, a8 , badDayConv);
    PARAM(TRUE , String, a9 , holidays);
    
    SET_DEFAULT(stubType, "f/s");
    if (JpmcdsStringToStubMethod(stubType, &stub) != SUCCESS)
        goto done;

    SET_DEFAULT(paymentDcc, "ACT/360");
    if (JpmcdsStringToDayCountConv(paymentDcc, &dcc) != SUCCESS)
        goto done;

    SET_DEFAULT(badDayConv, "N");
    if (JpmcdsBadDayConvValid(routine, badDayConv[0]) != SUCCESS)
        goto done;

    if (JpmcdsStringToDateInterval(couponInterval, routine, &ivl) != SUCCESS)
        goto done;

    cfl = JpmcdsCdsFeeLegFlows(startDate,
                               endDate,
                               &ivl,
                               &stub,
                               notional,
                               rate,
                               dcc,
                               badDayConv[0],
                               holidays);

    if (cfl == NULL)
        goto done;
    
    /* write cashflow list to a 2 column Excel array */
    a0->xltype = xltypeMulti;
    a0->val.array.rows = cfl->fNumItems;
    a0->val.array.columns = 2;
    a0->val.array.lparray = NEW_ARRAY(XLOPER, 2 * cfl->fNumItems);
    if (a0->val.array.lparray == NULL)
        goto done;
            
    for (i = 0; i < cfl->fNumItems; i++)
    {
        a0->val.array.lparray[2*i].xltype  = xltypeNum;
        JpmcdsExcelSetDate(cfl->fArray[i].fDate, a0->val.array.lparray + 2*i);
        a0->val.array.lparray[2*i+1].xltype  = xltypeNum;
        a0->val.array.lparray[2*i+1].val.num  = cfl->fArray[i].fAmount;
    }

    status = SUCCESS;

done:
    JpmcdsFreeCFL(cfl);
    FREE(couponInterval);
    FREE(stubType);
    FREE(paymentDcc);
    FREE(badDayConv);
    FREE(holidays);
    EXCEL_RETURN;
}

/*
***************************************************************************
** Wrapper for Excel addin DefaultedCDS function.
***************************************************************************
*/
XLFUNC(LPXLOPER) CDS_DefaultAccrual(OPER *a1, OPER *a2, OPER *a3, OPER *a4, 
        OPER *a5, OPER *a6, OPER *a7, OPER *a8, OPER *a9, OPER *a10, OPER *a11)
{
    
    static char   *routine = "CDS_DefaultAccrual";
    int            status = FAILURE;
    LPXLOPER       a0 = NULL;
    
    TDate          tradeDate;
    TDate          edd;
    TDate          startDate;
    TDate          endDate;
    double         couponRate;
    double         notional;
    char          *couponInterval = NULL;
    char          *stubType = NULL;
    char          *paymentDcc = NULL;
    char          *badDayConv = NULL;
    char          *holidays = NULL;
    
    double        defaultAccrual;
    double        accrualDays;

    TStubMethod    stub;
    long           dcc;
    TDateInterval  ivl;
   

    a0 = NEW(XLOPER);
    if (a0 == NULL)
        goto done;

    JpmcdsClearExcelDateSystem();
    
    PARAM(TRUE , Date  , a1 , tradeDate);
    PARAM(TRUE , Date  , a2 , edd);
    PARAM(TRUE , Date  , a3 , startDate);
    PARAM(TRUE , Date  , a4 , endDate);
    PARAM(TRUE , Double, a5 , notional);
    PARAM(TRUE , Double, a6 , couponRate);
    PARAM(TRUE , String, a7 , couponInterval);
    PARAM(FALSE, String, a8 , stubType);
    PARAM(FALSE, String, a9 , paymentDcc);
    PARAM(FALSE, String, a10, badDayConv);
    PARAM(TRUE , String, a11, holidays);

    SET_DEFAULT(stubType, "f/s");
    if (JpmcdsStringToStubMethod(stubType, &stub) != SUCCESS)
        goto done;
    
    SET_DEFAULT(paymentDcc, "ACT/360");
    if (JpmcdsStringToDayCountConv(paymentDcc, &dcc) != SUCCESS)
        goto done;
    
    if (JpmcdsStringToDateInterval(couponInterval, routine, &ivl) != SUCCESS)
        goto done;

    SET_DEFAULT(badDayConv, "N");
    if (JpmcdsBadDayConvValid(routine, badDayConv[0]) != SUCCESS)
        goto done;

    if (JpmcdsDefaultAccrual(
            tradeDate, 
            edd,
            startDate,
            endDate,
            &ivl,
            &stub,
            notional,
            couponRate,
            dcc,
            badDayConv[0],
            holidays,
            &accrualDays,
            &defaultAccrual
            ) != SUCCESS) goto done;
    
    a0->xltype = xltypeMulti;
    a0->val.array.rows = 2;
    a0->val.array.columns = 1;
    a0->val.array.lparray = NEW_ARRAY(XLOPER, 2);
    if (a0->val.array.lparray == NULL)
        goto done;
            
    a0->val.array.lparray[0].xltype  = xltypeNum;
    a0->val.array.lparray[0].val.num  = accrualDays;
    a0->val.array.lparray[1].xltype  = xltypeNum;
    a0->val.array.lparray[1].val.num  = defaultAccrual;

    status = SUCCESS;
done:
    FREE(couponInterval);
    FREE(stubType);
    FREE(paymentDcc);
    FREE(holidays);     
    
    EXCEL_RETURN;
}



