/* 
 * File:   defaulted.h
 * Author: john.yang
 *
 * Created on October 16, 2012, 5:31 PM
 */

#ifndef DEFAULTED_H
#define	DEFAULTED_H

#include "cdate.h"
#include "stub.h"

#ifdef	__cplusplus
extern "C" {
#endif

/*f
 */
EXPORT int JpmcdsDefaultAccrual(
        TDate           tradeDate, 
        TDate           edd, 
        TDate           startDate, 
        TDate           endDate, 
        TDateInterval  *couponInterval,
        TStubMethod    *stubType,
        double          notional,
        double          couponRate,
        long            paymentDcc,
        long            badDayConv,
        char           *calendar,
        double         *accrualDays,
        double         *defaultAccrual);

EXPORT int test(double *testvar);

#ifdef	__cplusplus
}
#endif

#endif	/* DEFAULTED_H */

