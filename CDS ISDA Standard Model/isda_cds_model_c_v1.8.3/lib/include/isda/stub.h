/*
 * ISDA CDS Standard Model
 *
 * Copyright (C) 2009 International Swaps and Derivatives Association, Inc.
 * Developed and supported in collaboration with Markit
 * 
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the ISDA CDS Standard Model Public License.
 */

#ifndef STUB_H
#define STUB_H

#include "yearfrac.h"

#ifdef __cplusplus
extern "C"
{
#endif


/* The following enum is used to pass the default stub position */
typedef enum
{
JPMCDS_STUB_POS_DEFAULT_AUTO  = 0,
JPMCDS_STUB_POS_DEFAULT_FRONT = 1,
JPMCDS_STUB_POS_DEFAULT_BACK  = 2
} TStubPos;


typedef long int TBadDayAndStubPos;
#define JPMCDS_BDC_MASK       0x000000FF
#define JPMCDS_STUB_POS_MASK  0x0000FF00
#define JPMCDS_STUB_POS_SHIFT 8

/* To construct "long stubPlacement" bitfield argument
 * (successor to "TBoolean stubAtEnd", encoding both position & length).
 */

/* Bit #0 */
#define JPMCDS_STUB_POSN_MASK    1
#define JPMCDS_STUB_POSN_FRONT   0
#define JPMCDS_STUB_POSN_BACK    1

/* Bit #1 */
#define JPMCDS_STUB_LENGTH_MASK  2
#define JPMCDS_STUB_LENGTH_SHORT 0
#define JPMCDS_STUB_LENGTH_LONG  2


/*t
 * TStubMethod combines all the types of stub information for swaps.
 */
typedef struct _TStubMethod
{
    TBoolean stubAtEnd;
    TBoolean longStub;
} TStubMethod;


/*f
***************************************************************************
** Computes stub payment.
**
** Minus accrued only set for Bond option with Front Stub.
***************************************************************************
*/
int JpmcdsStubPayment
   (TDate   prevCouponDate,      /* (I) Needed for Actual/Actual */
    TDate   nextCoupon,          /* (I) Needed for Actual/Actual*/
    TDate   stubStart,           /* (I) When stub starts */
    TDate   stubEnd,             /* (I) When stub ends */
    double  rate,                /* (I) Annual rate */
    long    couponDayCountConv,  /* (I) See JpmcdsDayCountFraction */
    double *stubPayment);        /* (O) At stub End date */


/*f
***************************************************************************
** This function converts a string to a stub method.
***************************************************************************
*/
int JpmcdsStringToStubMethod
(char        *name,      /* (I) Stub method name */
 TStubMethod *stubMethod /* (O) Stub method returned */
);


/*f
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
 TBoolean         *isEndStub);      /* (O) TRUE if end stub */


/*f
*********************************************************************
** Creates a TStubPos from a long.
*********************************************************************
*/
int JpmcdsStubPosMake
(long      stubPosValue, /* (I) Number representing stub position */
 TStubPos *stubPos);     /* (O) Stub position */


/*f
*********************************************************************
** Extracts TStubPos and bad day convention from TBadDayAndStubPos.
*********************************************************************
*/
int JpmcdsBadDayAndStubPosSplit
(TBadDayAndStubPos  badDayAndStubPos, /* (I) BDC and stub position */
 long              *badDayConv,       /* (O) Bad day convention */
 TStubPos          *stubPos);         /* (O) Stub position */


#ifdef __cplusplus
}
#endif

#endif    /* STUB_H */
