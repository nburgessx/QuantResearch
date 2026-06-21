/*
 * ISDA CDS Standard Model
 *
 * Copyright (C) 2009 International Swaps and Derivatives Association, Inc.
 * Developed and supported in collaboration with Markit
 * 
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the ISDA CDS Standard Model Public License.
 */

#include <math.h>
#include <ctype.h>
#include <string.h>
#include "cgeneral.h"
#include "bastypes.h"
#include "cerror.h"
#include "strutil.h"
#include "stub.h"
#include "convert.h"
#include "ldate.h"
#include "date_sup.h"
#include "macros.h"


/*
***************************************************************************
** Computes stub payment.
**
** Minus accrued only set for Bond option with Front Stub.
***************************************************************************
*/
int JpmcdsStubPayment
    (TDate   prevCouponDate,            /* (I) Needed for Actual/Actual */
     TDate   nextCouponDate,            /* (I) Needed for Actual/Actual */
     TDate   stubStart,                 /* (I) When accrual starts */
     TDate   stubEnd,                   /* (I) When accrual ends */
     double  rate,                      /* (I) Annual rate */
     long    couponDayCountConv,        /* (I)Day count convention */
     double *stubPayment)               /* (O) At accrue End date */
{
    static char routine[]="JpmcdsStubPayment";
    double stubYears;                   /* Number years in stub. */
    double stubFrac;                    /* Fraction of coupon  */
    double coupon, couponYearFrac;
    
    if (prevCouponDate >= nextCouponDate)
    {
        JpmcdsErrMsg("%s: prevCouponDate (%s) >= nextCouponDate (%s).\n", 
                  routine, JpmcdsFormatDate(prevCouponDate),
                  JpmcdsFormatDate(nextCouponDate));
        return FAILURE;
    }

    if (stubStart < prevCouponDate || stubStart >  nextCouponDate)
    {
        JpmcdsErrMsg
         ("%s: StubStart(%s) < prevCouponDate(%s) OR > nextCouponDate(%s).\n",
          routine, JpmcdsFormatDate(stubStart),
          JpmcdsFormatDate(prevCouponDate), 
          JpmcdsFormatDate(nextCouponDate));
        return FAILURE;
    }

    if (stubEnd > nextCouponDate || stubEnd < prevCouponDate)
    {
        JpmcdsErrMsg
        ("%s: StubEnd(%s) > nextCouponDate(%s) OR < prevCouponDate(%s).\n",
         routine, JpmcdsFormatDate(stubEnd),
         JpmcdsFormatDate(nextCouponDate), 
         JpmcdsFormatDate(prevCouponDate));
        return FAILURE ;
    }

    if (stubStart != prevCouponDate && stubEnd != nextCouponDate)
    {
        JpmcdsErrMsg("%s: StubStart(%s) <> prevCouponDate(%s) AND"
                  " StubEnd(%s) <> nextCouponDate(%s).\n",
                  routine, JpmcdsFormatDate(stubStart),
                  JpmcdsFormatDate(prevCouponDate), 
                  JpmcdsFormatDate(stubEnd), 
                  JpmcdsFormatDate(nextCouponDate));
        return FAILURE;
    }

    /* Compute coupon year fraction and coupon.
     */
    if (JpmcdsDayCountFraction(prevCouponDate, nextCouponDate,
                            couponDayCountConv, &couponYearFrac) == FAILURE)
        return FAILURE;
    coupon = rate * couponYearFrac;


    /* Compute stub year fraction and stub fraction. */
    if (JpmcdsDayCountFraction(stubStart, stubEnd, couponDayCountConv, &stubYears) == FAILURE)
        return FAILURE;

    stubFrac = stubYears / couponYearFrac;
    
    /* Compute stub */
    *stubPayment = rate * stubYears;

    return SUCCESS;
}


/*
***************************************************************************
** This function converts a string to a stub method.
***************************************************************************
*/
int JpmcdsStringToStubMethod
(char        *name,      /* (I) Stub method name */
 TStubMethod *stubMethod /* (O) Stub method returned */
)
{
    static char routine[] = "JpmcdsStringToStubMethod";
    static char invalidFormatMsg[]  = "%s: Invalid stub method format %s. "
        "Format is LOCATION/TYPE/[LENGTH].\n"
        "    TYPE and LENGTH are optional elements in this string.\n"
        "    LOCATION is F=front or B=back.\n"
        "    LENGTH is S=short (default), L=long.\n";

    int         status    = FAILURE;

    char       *buf = NULL; /* Deep copy of name - used during parsing */
    char       *locStr;
    char       *lengthStr;

    TBoolean    stubAtEnd = FALSE;
    TBoolean    longStub = FALSE;

    if (name == NULL || stubMethod == NULL)
    {
        JpmcdsErrMsg ("%s: NULL inputs.\n", routine);
        goto done;
    }

    buf = JpmcdsStringDuplicate (name);
    if (buf == NULL)
        goto done;

    if (JpmcdsStringParser (buf, "/", &lengthStr) != SUCCESS)
        goto done;

    locStr = buf;

    switch (toupper(*locStr))
    {
    case 'F':
        stubAtEnd = FALSE;
        break;
    case 'B':
        stubAtEnd = TRUE;
        break;
    default:
        JpmcdsErrMsg (invalidFormatMsg, routine, name);
        goto done;
    }

    if ( lengthStr == NULL )
    {
        longStub = FALSE;
    }
    else
    {
        switch (toupper(*lengthStr))
        {
        case 'S':
            longStub = FALSE;
            break;
        case 'L':
            longStub = TRUE;
            break;
        default:
            JpmcdsErrMsg(invalidFormatMsg, routine);
            goto done;
        }
    }

/*
** Populate the output.
*/
    stubMethod->stubAtEnd = stubAtEnd;
    stubMethod->longStub = longStub;

    status = SUCCESS;

 done: 
    FREE (buf);

    if (status != SUCCESS)
        JpmcdsErrMsgFailure (routine);

    return status;
}


/*
*********************************************************************
** Returns TRUE if a stub should be at the end.
**
** The decision is based on the following.
** 1. First the default stub position is checked. If this is a back
** stub then return TRUE - if a front stub return FALSE.
** 2. If the default stub position is auto then the function will
** return TRUE, unless there is a stub in whach case FALSE is
** returned.
*********************************************************************
*/
int JpmcdsIsEndStub
(TDate             startDate,       /* (I) Start date */
 TDate             maturityDate,    /* (I) Maturity date */
 TDateInterval    *ivl,             /* (I) Interval */
 TStubPos          stubPos,         /* (I) Stub position */
 TBoolean         *isEndStub)       /* (O) TRUE if end stub */
{
    static char  routine[] = "JpmcdsIsEndStub";
    int          status = FAILURE;
    int          numIvls;
    int          extraDays;

    switch(stubPos)
    {
    case JPMCDS_STUB_POS_DEFAULT_BACK:
        *isEndStub = TRUE;
        status = SUCCESS;
        break;
    case JPMCDS_STUB_POS_DEFAULT_FRONT:
        *isEndStub = FALSE;
        status = SUCCESS;
        break;
    case JPMCDS_STUB_POS_DEFAULT_AUTO:
        if (JpmcdsCountDates(startDate,
                          maturityDate,
                          ivl,
                          &numIvls,
                          &extraDays) != SUCCESS)
        {
            goto done;
        }

        if (extraDays > 0)
        {
            *isEndStub = FALSE;
        }
        else
        {
            *isEndStub = TRUE;
        }
        status = SUCCESS;
        break;
    }

done:
    if (status != SUCCESS)
        JpmcdsErrMsgFailure(routine);

    return status;
}


/*
*********************************************************************
** Creates a TStubPos from a long.
*********************************************************************
*/
int JpmcdsStubPosMake
(long      stubPosValue, /* (I) Number representing stub position */
 TStubPos *stubPos)      /* (O) Stub position */
{
    static char routine[] = "JpmcdsStubPosMake";
    int         status = FAILURE;

    switch(stubPosValue)
    {
    case JPMCDS_STUB_POS_DEFAULT_AUTO:
        *stubPos = JPMCDS_STUB_POS_DEFAULT_AUTO;
        break;
    case JPMCDS_STUB_POS_DEFAULT_FRONT:
        *stubPos = JPMCDS_STUB_POS_DEFAULT_FRONT;
        break;
    case JPMCDS_STUB_POS_DEFAULT_BACK:
        *stubPos = JPMCDS_STUB_POS_DEFAULT_BACK;
        break;
    default:
        JpmcdsErrMsg("%s: Unknown stub pos %ld.\n", routine, stubPosValue);
        goto done;
    }

    status = SUCCESS;

done:
    if (status != SUCCESS)
        JpmcdsErrMsgFailure(routine);

    return status;
}


/*
*********************************************************************
** Extracts TStubPos and bad day convention from TBadDayAndStubPos.
*********************************************************************
*/
int JpmcdsBadDayAndStubPosSplit
(TBadDayAndStubPos  badDayAndStubPos, /* (I) BDC and stub position */
 long              *badDayConv,       /* (O) Bad day convention */
 TStubPos          *stubPos)          /* (O) Stub position */
{
    static char routine[] = "JpmcdsBadDayAndStubPosSplit";
    int         status = FAILURE;

    long value = (long) badDayAndStubPos;
    long bdc = value & JPMCDS_BDC_MASK;
    long sp = (value & JPMCDS_STUB_POS_MASK)>>JPMCDS_STUB_POS_SHIFT;

    if (JpmcdsStubPosMake(sp, stubPos) != SUCCESS)
    {
        goto done;
    }
    *badDayConv = bdc;

    status = SUCCESS;

done:
    if (status != SUCCESS)
        JpmcdsErrMsgFailure(routine);

    return status;
}
