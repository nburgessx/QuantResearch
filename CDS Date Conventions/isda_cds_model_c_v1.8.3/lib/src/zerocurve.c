/*
 * ISDA CDS Standard Model
 *
 * Copyright (C) 2009 International Swaps and Derivatives Association, Inc.
 * Developed and supported in collaboration with Markit
 * 
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the ISDA CDS Standard Model Public License.
 */

#include "zerocurve.h"
#include "tcurve.h"
#include "macros.h"
#include "cerror.h"
#include "gtozc.h"
#include <ctype.h>

#define JpmcdsSWAPNAME   'S'
#define JpmcdsMONEYNAME  'M'


/*
***************************************************************************
** Build zero curve from money market, and swap instruments.
***************************************************************************
*/
EXPORT TCurve* JpmcdsBuildIRZeroCurve(
    TDate      valueDate,       /* (I) Value date                       */
    char      *instrNames,      /* (I) Array of 'M' or 'S'              */
    TDate     *dates,           /* (I) Array of swaps dates             */
    double    *rates,           /* (I) Array of swap rates              */
    long       nInstr,          /* (I) Number of benchmark instruments  */
    long       mmDCC,           /* (I) DCC of MM instruments            */
    long       fixedSwapFreq,   /* (I) Fixed leg freqency               */
    long       floatSwapFreq,   /* (I) Floating leg freqency            */
    long       fixedSwapDCC,    /* (I) DCC of fixed leg                 */
    long       floatSwapDCC,    /* (I) DCC of floating leg              */
    long       badDayConv,      /* (I) Bad day convention               */
    char      *holidayFile)     /* (I) Holiday file                     */
{
    static char routine[] = "BuildIRZeroCurve";
    int         status    = FAILURE;

    int     i;
    int     nCash = 0;
    int     nSwap = 0;
    char    fwdLength = '3';   /* not used */
    char    instr;

    TDate  *cashDates    = NULL;
    TDate  *swapDates    = NULL;
    double *cashRates    = NULL;
    double *swapRates    = NULL;
    TCurve *zcurveIni    = NULL;
    TCurve *zcurveCash   = NULL;
    TCurve *zcurveSwap   = NULL;

    /* Allocate enough spaces for cash and swap dates/rates */
    cashDates  = NEW_ARRAY(TDate,  nInstr);
    swapDates  = NEW_ARRAY(TDate,  nInstr);
    cashRates  = NEW_ARRAY(double, nInstr);
    swapRates  = NEW_ARRAY(double, nInstr);

    if (cashDates  == NULL || swapDates  == NULL || cashRates  == NULL || swapRates  == NULL)
        goto done;

    /* Sort out cash and swap separately */
    for(i = 0; i < nInstr; i++)
    {
        instr = toupper(instrNames[i]);
        if (instr != JpmcdsMONEYNAME && instr != JpmcdsSWAPNAME)
        {
            JpmcdsErrMsg("%s: unknown instrument type (%c)."
                      " Only (M)oney market or (S)wap is allowed.\n",
                      routine, instrNames[i]);
            goto done;
        }

        if (instr == JpmcdsMONEYNAME)
        {    /* MM Rate */
            cashDates[nCash]  = dates[i];
            cashRates[nCash]  = rates[i];
            nCash++;
        }
        else        /* Swap Rate */
        {
            swapDates[nSwap]  = dates[i];
            swapRates[nSwap]  = rates[i];
            nSwap++;
        }
    }

    /* Initialize the zero curve */
    zcurveIni = JpmcdsNewTCurve(valueDate, 0, (double) 1L, JPMCDS_ACT_365F);
    if (zcurveIni == NULL)
        goto done;

    /* Cash instruments */
    zcurveCash = JpmcdsZCCash(zcurveIni, cashDates, cashRates, nCash, mmDCC);
    if (zcurveCash == NULL)
        goto done;

    /* Swap instruments */
    zcurveSwap = JpmcdsZCSwaps(zcurveCash,
                            NULL,   /* discZC */
                            swapDates,
                            swapRates,
                            nSwap,
                            fixedSwapFreq,
                            floatSwapFreq,
                            fixedSwapDCC,
                            floatSwapDCC,
                            fwdLength,
                            badDayConv,
                            holidayFile);
    if (zcurveSwap == NULL)
        goto done;

    status = SUCCESS;

 done:
    FREE(cashDates);
    FREE(cashRates);
    FREE(swapDates);
    FREE(swapRates);
    JpmcdsFreeTCurve(zcurveIni);
    JpmcdsFreeTCurve(zcurveCash);
    if (status != SUCCESS)
    {
        JpmcdsFreeTCurve(zcurveSwap);
        zcurveSwap = NULL;
        JpmcdsErrMsgFailure(routine);
    }
    return (zcurveSwap);
}
