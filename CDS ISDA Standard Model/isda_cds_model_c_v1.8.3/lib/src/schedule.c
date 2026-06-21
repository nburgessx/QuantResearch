/*
 * ISDA CDS Standard Model
 *
 * Copyright (C) 2009 International Swaps and Derivatives Association, Inc.
 * Developed and supported in collaboration with Markit
 * 
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the ISDA CDS Standard Model Public License.
 */

#include <ctype.h>
#include <stdlib.h>
#include <string.h>

#include "cerror.h"
#include "cmemory.h"
#include "cashflow.h"
#include "datelist.h"
#include "date_sup.h"
#include "ldate.h"
#include "macros.h"
#include "schedule.h"


/*
***************************************************************************
** Checks arguments to JpmcdsNewCouponDatesSwap.
***************************************************************************
*/
static int CheckArgsSwapMakeSchedule(
    TDate           startDate,
    TDate           matDate,
    long            payBadDayConv,
    long            accBadDayConv,
    long            resetBadDayConv,
    int             resetOffsetDays,
    TDate           firstRollDate,
    TDate           lastRollDate,
    TDate           fullFirstCoupDate);


/*
***************************************************************************
** Creates a date schedule for one stream of a swap.
***************************************************************************
*/
TCouponDateList *   JpmcdsNewCouponDatesSwap(
    TDate           startDate,          /* (I) Effective date of stream     */
    TDate           matDate,            /* (I) Maturity date of stream      */
    TDateInterval   *couponInterval,    /* (I) Time between payments        */
    TBoolean        adjustLastAccDate,  /* (I) Maturity date adjusted       */
    TBoolean        inArrears,          /* (I) Arrears setting flag         */
    int             payOffsetDays,      /* (I) Days from acc end to payment */
    int             resetOffsetDays,    /* (I) Days from accrual to reset   */
    TBoolean        stubAtEnd,          /* (I) TRUE=back stub,FALSE=frontstub*/
    TBoolean        longStub,           /* (I) TRUE=long stub,FALSE=shortstub*/
    TDate           firstRollDate,      /* (I) Date of first cashflow       */
    TDate           lastRollDate,       /* (I) Acc start of final period    */
    TDate           fullFirstCoupDate,  /* (I) Acc start of first period    */
    long            payBadDayConv,      /* (I) Bad day convention for payment*/
    long            accBadDayConv,      /* (I) Bad day convention for accrual*/
    long            resetBadDayConv,    /* (I) Bad day convention for resets*/
    char            *holidayFile)       /* (I) See JpmcdsBusinessDay           */
{
   int              i;                    /* loops over the unadjusted dates*/
   int              j;                    /* */
   int              retCode = FAILURE;    /* return value                   */
   long             stubInfo;             /* holds stub location            */
   TDate            tempAccDate;          /* temporarily stores adjusted date*/
   TDateList        *unadjDates = NULL;   /* holds array of unadjusted dates*/
   TCouponDateList  *adjStreamSched = NULL;/* holds stream pending removal
                                               of empty accrual periods */
   TCouponDateList  *streamSched = NULL;  /* holds the final output stream  */

   static char routine[] = "JpmcdsNewCouponDatesSwap";

   if (CheckArgsSwapMakeSchedule(startDate,
                                 matDate,
                                 payBadDayConv,
                                 accBadDayConv,
                                 resetBadDayConv,
                                 resetOffsetDays,
                                 firstRollDate,
                                 lastRollDate,
                                 fullFirstCoupDate)==FAILURE)
   {
       JpmcdsErrMsg("%s: checking arguments uncovered an error\n", routine);
       goto error;
   }


   /* generate a list of unadjusted dates */

   if (JpmcdsGenerateUnadjDates(startDate,
                                matDate,
                                couponInterval,
                                stubAtEnd,
                                longStub,
                                firstRollDate,
                                lastRollDate,
                               &unadjDates,
                               &stubInfo)==FAILURE)
   {
       JpmcdsErrMsg("%s: failed to generate unadjusted coupon dates.\n", routine);
       goto error;
   }


   if (unadjDates == NULL)
   {
       JpmcdsErrMsg("%s: error generating unadjusted dates.\n", routine);
       goto error;
   }

   adjStreamSched = JpmcdsNewEmptyCouponDateList(unadjDates->fNumItems-1);
   if (adjStreamSched == NULL)
   {
       JpmcdsErrMsg("%s:Failed to generate new pre-adjusted stream schedule.\n", routine);
       goto error;
   }
   adjStreamSched->stubLocation = stubInfo;

   /* adjust accrual dates */

   for (i=0 ;i<unadjDates->fNumItems ; i++)
   {
       if (!(i==unadjDates->fNumItems-1) ||
             adjustLastAccDate == TRUE)
       {
           if (JpmcdsBusinessDay(
                    unadjDates->fArray[i],
                    accBadDayConv,
                    holidayFile,
                    &tempAccDate)==FAILURE)
           {
               JpmcdsErrMsg("%s: error when adjusting the accrual date.\n",
                          routine);
               goto error;
           }
       }
       else
       {
           tempAccDate = unadjDates->fArray[i];
       }
       if (i>0)
       {
           adjStreamSched->fArray[i-1].accrueEndDate = tempAccDate;
       }
       if (i<unadjDates->fNumItems-1)
       {
           adjStreamSched->fArray[i].accrueStartDate = tempAccDate;
       }
   }

   /* Prepare final output stream */
   streamSched = JpmcdsNewEmptyCouponDateList(unadjDates->fNumItems-1);
   if (streamSched == NULL)
   {
       JpmcdsErrMsg("%s: failed to generate new stream schedule.\n", routine);
       goto error;
   }
   streamSched->stubLocation = stubInfo;

   /*
    * Remove empty accrual periods and adjust stub location
    * as appropriate.
    * Thanks to Neill Penney for the solution.
    */
   for(i=0, j=0; i<adjStreamSched->fNumItems; i++ )
   {
      if ( adjStreamSched->fArray[i].accrueStartDate !=
           adjStreamSched->fArray[i].accrueEndDate )
      {
          streamSched->fArray[j] = adjStreamSched->fArray[i];
          j++;
      }
      else
      {
         if ( i==0 )
         {
            if ( adjStreamSched->stubLocation == JPMCDS_FRONT_STUB )
            {
               streamSched->stubLocation = JPMCDS_NO_STUB;
            }
            else
            if ( adjStreamSched->stubLocation == JPMCDS_TWO_STUB )
            {
               streamSched->stubLocation = JPMCDS_BACK_STUB;
            }
            else
            {
               /*
                * Removed empty accrual period from front
                * but no front stub. This is OK.
                */
            }
            streamSched->fNumItems--;
         }
         else
         if ( i==adjStreamSched->fNumItems-1 )
         {
            if ( adjStreamSched->stubLocation == JPMCDS_BACK_STUB )
            {
               streamSched->stubLocation = JPMCDS_NO_STUB;
            }
            else
            if ( adjStreamSched->stubLocation == JPMCDS_TWO_STUB )
            {
               streamSched->stubLocation = JPMCDS_FRONT_STUB;
            }
            else
            {
               /*
                * Removed empty accrual period from back
                * but no back stub. This is OK.
                */
            }
            streamSched->fNumItems--;
         }
         else
         {
            /*
             * Removed empty accrual period from middle
             * of schedule. This is OK.
             */
            streamSched->fNumItems--;
         }
      }
   }


   /* Confirm expected number of items */
   if ( j != streamSched->fNumItems )
   {
      JpmcdsErrMsg("%s: Error in number of dates in stream: \n"
                "Index=> %d, fNumItems=%d\n",
                routine, j, streamSched->fNumItems);
      PROGRAM_BUG();
      goto error;
   }

   /* calculate the payment dates */

   for (i=0 ;i<streamSched->fNumItems ;i++ )
   {
       if (JpmcdsBusinessDay(
              streamSched->fArray[i].accrueEndDate,
              payBadDayConv,
              holidayFile,
              &streamSched->fArray[i].payDate)==FAILURE)
       {
           JpmcdsErrMsg("%s: error when adjusting payment date(1).\n",
                      routine);
           goto error;
       }
   }

   if (payOffsetDays != 0)
   {
       for (i=0 ;i<streamSched->fNumItems ;i++ )
       {
           if (JpmcdsDateFromBusDaysOffset(
                   streamSched->fArray[i].payDate,
                   payOffsetDays,
                   holidayFile,
                   &streamSched->fArray[i].payDate)==FAILURE)
           {
                   JpmcdsErrMsg("%s: error when adjusting payment date(2).\n",
                              routine);
                   goto error;
           }
       }
   }


   /* calculate the reset dates */


   for (i=0 ;i<streamSched->fNumItems ; i++ )
   {
       if (!inArrears)
       {
           streamSched->fArray[i].resetDate = streamSched->fArray[i].accrueStartDate;
       }
       else
       {
          streamSched->fArray[i].resetDate = streamSched->fArray[i].accrueEndDate;
       }


       if (JpmcdsBusinessDay(
              streamSched->fArray[i].resetDate,
              resetBadDayConv,
              holidayFile,
              &streamSched->fArray[i].resetDate)==FAILURE)
       {
           JpmcdsErrMsg("%s: error when calculating reset date.\n",
                      routine);
           goto error;
       }

       if (resetOffsetDays != 0)
       {
           if (JpmcdsDateFromBusDaysOffset(
                  streamSched->fArray[i].resetDate,
                  -resetOffsetDays,
                  holidayFile,
                  &streamSched->fArray[i].resetDate)==FAILURE)
           {
               JpmcdsErrMsg("%s: error when calculating reset date.\n",
                          routine);
               goto error;
           }
       }
       /* According to STIRT Traders, ISDA rules are incorrect:
        * Payment can occur on reset date
      */
       while (streamSched->fArray[i].resetDate >
              streamSched->fArray[i].payDate)
       {
           if (JpmcdsBusinessDay(
                  --streamSched->fArray[i].resetDate,
                  JPMCDS_BAD_DAY_PREVIOUS,
                  holidayFile,
                  &streamSched->fArray[i].resetDate)==FAILURE)
           {
               JpmcdsErrMsg("%s: error when adjusting reset date.\n",
                          routine);
               goto error;
           }
       }
   }



   if (streamSched->stubLocation > JPMCDS_BACK_STUB &&
       fullFirstCoupDate > 0)
   {
       streamSched->fArray[0].accrueStartDate = fullFirstCoupDate;
   }

   retCode = SUCCESS;

   error:
   JpmcdsFreeCouponDateList(adjStreamSched);
       if (retCode == FAILURE)
       {
          JpmcdsErrMsg("%s: Failed.\n", routine);
          JpmcdsFreeCouponDateList(streamSched);
          streamSched = NULL;
       }
   JpmcdsFreeDateList(unadjDates);
   return streamSched;
}


/*
***************************************************************************
** Checks arguments to JpmcdsNewCouponDatesSwap.
***************************************************************************
*/
static int CheckArgsSwapMakeSchedule(
    TDate           startDate,
    TDate           matDate,
    long            payBadDayConv,
    long            accBadDayConv,
    long            resetBadDayConv,
    int             resetOffsetDays,
    TDate           firstRollDate,
    TDate           lastRollDate,
    TDate           fullFirstCoupDate)
{
   static char routine[] = "CheckArgsSwapMakeSchedule";

   if (startDate > matDate)
   {
       JpmcdsErrMsg("%s: maturity date (%s) is on or before start date (%s)\n",
                  routine,
                  JpmcdsFormatDate(matDate),
                  JpmcdsFormatDate(startDate));
       return FAILURE;
   }

   if (JpmcdsBadDayConvValid(routine, payBadDayConv)==FAILURE)
   {
       JpmcdsErrMsg("%s: failed to check pay bad day convention.\n", routine);
       return FAILURE;
   }
   if (JpmcdsBadDayConvValid(routine, accBadDayConv)==FAILURE)
   {
       JpmcdsErrMsg("%s: failed to check acc bad day convention.\n", routine);
       return FAILURE;
   }
   if (JpmcdsBadDayConvValid(routine, resetBadDayConv)==FAILURE)
   {
       JpmcdsErrMsg("%s: failed to check reset bad day convention.\n", routine);
       return FAILURE;
   }

   if (fullFirstCoupDate > startDate)
   {
       JpmcdsErrMsg("%s: full first coupon date (%s) is after"
                 "swap start date (%s).\n",
                 routine,
                 JpmcdsFormatDate(fullFirstCoupDate),
                 JpmcdsFormatDate(startDate));
       return FAILURE;
   }

   if (firstRollDate > 0)
   {
       if (lastRollDate > 0)
       {
           JpmcdsErrMsg("%s: firstRollDate(%s) and lastRollDate(%s) cannot"
                     "both be defined.\n",
                      routine,
                      JpmcdsFormatDate(firstRollDate),
                      JpmcdsFormatDate(lastRollDate));
           return FAILURE;
       }
       if (firstRollDate < startDate)
       {
           JpmcdsErrMsg("%s: firstRollDate(%s) is before start date (%s).\n",
                      routine,
                      JpmcdsFormatDate(firstRollDate),
                      JpmcdsFormatDate(startDate));
           return FAILURE;
       }
       if (firstRollDate > matDate)
       {
           JpmcdsErrMsg("%s: firstRollDate(%s) is after maturity (%s).\n",
                      routine,
                      JpmcdsFormatDate(firstRollDate),
                      JpmcdsFormatDate(matDate));
           return FAILURE;
       }
   }
   if (lastRollDate > 0)
   {
       if (lastRollDate < startDate)
       {
           JpmcdsErrMsg("%s: lastRollDate(%s) is before start date (%s).\n",
                      routine,
                      JpmcdsFormatDate(lastRollDate),
                      JpmcdsFormatDate(startDate));
           return FAILURE;
       }
       if (lastRollDate > matDate)
       {
           JpmcdsErrMsg("%s: lastRollDate(%s) is after maturity (%s).\n",
                      routine,
                      JpmcdsFormatDate(lastRollDate),
                      JpmcdsFormatDate(matDate));
           return FAILURE;
       }
   }

   if (resetOffsetDays < 0)
   {
       JpmcdsErrMsg("%s: offset for reset date cannot be negative.\n", routine);
       return FAILURE;
   }

   return SUCCESS;
}


/*
***************************************************************************
** Checks the details for a stream schedule.
***************************************************************************
*/
int JpmcdsCouponDateListCheck(TCouponDateList *schedule)
{
   int  i;

   static char routine[] = "JpmcdsCouponDateListCheck";

   if (schedule->fNumItems < 1)
   {
       JpmcdsErrMsg("%s: schedule must have at least one coupon period"
                 " but has only %d.\n",
                  routine,
                  schedule->fNumItems);
       return FAILURE;
   }
   if (schedule->stubLocation < JPMCDS_NO_STUB ||
       schedule->stubLocation > JPMCDS_TWO_STUB)
   {
       JpmcdsErrMsg("%s: incorrect number of stub periods (%ld).\n",
                  routine,
                  schedule->stubLocation);
       return FAILURE;

   }
   if (schedule->fArray[schedule->fNumItems-1].accrueEndDate <
       schedule->fArray[0].accrueStartDate)
   {
       JpmcdsErrMsg("%s: swap start date (%s) is after swap end date (%s).\n",
          routine,
          JpmcdsFormatDate(schedule->fArray[0].accrueStartDate),
          JpmcdsFormatDate(schedule->fArray[schedule->fNumItems-1].accrueEndDate));
       return FAILURE;
   }

   for (i=0; i<schedule->fNumItems; i++)
   {
      if (schedule->fArray[i].accrueEndDate <
           schedule->fArray[i].accrueStartDate)
       {
           JpmcdsErrMsg("%s: coupon #%d has accrue end date (%s) before"
                     " accrue start date (%s).\n",
                     routine,
                     i+1,
                     JpmcdsFormatDate(schedule->fArray[i].accrueEndDate),
                     JpmcdsFormatDate(schedule->fArray[i].accrueStartDate));
           return FAILURE;
       }
   }

   for (i=1 ;i<schedule->fNumItems ;i++ )
   {
       if (schedule->fArray[i-1].accrueStartDate >
           schedule->fArray[i].accrueStartDate)
       {
           JpmcdsErrMsg
               ("%s: acc start date #%d (%s) > acc start date #%d (%s).\n",
                routine,
                i,
                JpmcdsFormatDate(schedule->fArray[i-1].accrueStartDate),
                i+1,
                JpmcdsFormatDate(schedule->fArray[i].accrueStartDate));
           return FAILURE;
       }
       if (schedule->fArray[i-1].accrueEndDate >
           schedule->fArray[i].accrueEndDate)
       {
           JpmcdsErrMsg("%s: acc end date #%d (%s) > acc end date #%d (%s).\n",
                      routine,
                      i,
                      JpmcdsFormatDate(schedule->fArray[i-1].accrueEndDate),
                      i+1,
                      JpmcdsFormatDate(schedule->fArray[i].accrueEndDate));
           return FAILURE;
       }
       if (schedule->fArray[i-1].payDate >
           schedule->fArray[i].payDate)
       {
           JpmcdsErrMsg("%s: payment date #%d (%s) > payment date #%d (%s).\n",
                      routine,
                      i,
                      JpmcdsFormatDate(schedule->fArray[i-1].payDate),
                      i+1,
                      JpmcdsFormatDate(schedule->fArray[i].payDate));
           return FAILURE;
       }
       if (schedule->fArray[i-1].resetDate >
           schedule->fArray[i].resetDate)
       {
           JpmcdsErrMsg("%s: reset date #%d (%s) > reset date #%d (%s).\n",
                      routine,
                      i,
                      JpmcdsFormatDate(schedule->fArray[i-1].resetDate),
                      i+1,
                      JpmcdsFormatDate(schedule->fArray[i].resetDate));
           return FAILURE;
       }
   }

   return SUCCESS;
}


/*
***************************************************************************
** Creates a new empty TCouponDateList.
***************************************************************************
*/
TCouponDateList* JpmcdsNewEmptyCouponDateList(int numPeriods)
{
   TCouponDateList *streamSchedule = NULL;

   static char routine[] = "JpmcdsNewEmptyCouponDateList";

   if (numPeriods < 0)
   {
       JpmcdsErrMsg("%s: number of coupon periods (%d) must be >= 0 .\n",
                  routine,
                  numPeriods);
       return NULL;
   }

   streamSchedule = NEW(TCouponDateList);

   if (streamSchedule == NULL)
   {
       JpmcdsErrMsg("%s: out of memory when allocating TCouponDateList.\n", routine);
       return NULL;
   }

   streamSchedule->fNumItems = numPeriods;

   if (numPeriods>0)
   {
       streamSchedule->fArray = NEW_ARRAY(TCouponDates,numPeriods);
       if (streamSchedule->fArray == NULL)
       {
           JpmcdsErrMsg("%s: out of memory when allocating TCouponDateList.\n", routine);
           FREE(streamSchedule);
           return NULL;
       }
   }
   else
   {
       streamSchedule->fArray = NULL;
   }

   return streamSchedule;
}


/*
***************************************************************************
** Frees a TCouponDateList.
***************************************************************************
*/
void JpmcdsFreeCouponDateList(TCouponDateList *theSchedule)
{
   if (theSchedule != NULL)
   {
       FREE(theSchedule->fArray);
       FREE(theSchedule);
   }
}


/*
***************************************************************************
** Generate a list of unadjusted coupon dates for the stream.
***************************************************************************
*/
int JpmcdsGenerateUnadjDates(
    TDate           startDate,
    TDate           matDate,
    TDateInterval   *couponInterval,
    TBoolean        stubAtEnd,
    TBoolean        longStub,
    TDate           firstRollDate,
    TDate           lastRollDate,
    TDateList       **dateListPtr,
    long            *stubInfo)
{
   static char routine[] = "JpmcdsGenerateUnadjDates";
   int            status = FAILURE;   /* Until proven successful */
   int            i;                  /* loops over unadjusted dates        */
   int            index;              /* needed for reversing date order    */
   int            count;              /* used to calculate coupon dates     */
   int            numPeriods;         /* number of complete periods in stream*/
   int            dummyInt;           /* unused output from JpmcdsCountDates   */
   int            nDates;             /* number of unadjusted dates         */
   int            direction;          /* direction we are stepping in       */
   TDate          beginDate;          /* where we start stepping            */
   TDate          endDate;            /* where we finish stepping           */
   TDate          unadjustedDate;     /* holds date during calculations     */
   TDate          *dateArr = NULL;    /* stores the dates within routine    */
   TDateInterval  interval;           /* used with count to caculate dates  */
   TDateList      *dateList = NULL;   /* holds the output list of unadj dates*/


   *dateListPtr = NULL;
   *stubInfo = 0;

   /* figure out how big, roughly, the datelist is going to be */

   if (JpmcdsCountDates(startDate,
                    matDate,
                    couponInterval,
                    &numPeriods,
                    &dummyInt)==FAILURE)
       goto done;

   /* allocate memory to the temporary dates array */
   /* add potential stub to count of complete periods*/
   numPeriods = numPeriods + 3;

   dateArr = NEW_ARRAY(TDate,numPeriods);
   if (dateArr == NULL)
       goto done;


   /* do we have a first/last roll date? */
   if (firstRollDate > 0)
   {
       beginDate = firstRollDate;
       endDate   = matDate;
       direction = 1;

       nDates    = 0;
       dateArr[nDates++] = startDate;

       if (beginDate > startDate)
       {
           *stubInfo += JPMCDS_FRONT_STUB;
       }
   }
   else if (lastRollDate > 0)
   {
       beginDate = lastRollDate;
       endDate   = startDate;
       direction = -1;

       nDates    = 0;
       dateArr[nDates++] = matDate;

       if (matDate > beginDate)
       {
           *stubInfo += JPMCDS_BACK_STUB;
       }
   }
   else
   {
       if (stubAtEnd)
       {
           beginDate = startDate;
           endDate   = matDate;
           direction = 1;
       }
       else
       {
           beginDate = matDate;
           endDate   = startDate;
           direction = -1;
       }
       nDates = 0;
   }

   /* start looping */

   dateArr[nDates++] = beginDate;
   count = 1;
   do
   {
       if (JpmcdsMakeDateInterval(
                couponInterval->prd * direction * count,
                couponInterval->prd_typ,
                &interval)==FAILURE)
       {
           JpmcdsErrMsg("%s: Error generating coupon interval.\n", routine);
           goto done;
       }
       count++;
       if (JpmcdsDtFwdAny(beginDate,
                       &interval,
                       &unadjustedDate)==FAILURE)
       {
           JpmcdsErrMsg("%s: Error generating unadjusted dates.\n", routine);
           goto done;
       }
       if ((endDate-unadjustedDate) * direction > 0)
       {
           dateArr[nDates++] = unadjustedDate;
       }
       else
       {
           dateArr[nDates++] = endDate;
       }
   }
   while((endDate-unadjustedDate)*direction > 0 );

   /* check for stubs */

   if (endDate != unadjustedDate)
   {
       if (direction ==1)
       {
           *stubInfo += JPMCDS_BACK_STUB;
       }
       else
       {
           *stubInfo += JPMCDS_FRONT_STUB;
       }
   }

   /* if we are using a long stub, remove the appropriate date */
   /* Require 3 dates in order to perform this operation */

   if (firstRollDate == 0 && lastRollDate == 0)
   {
       if ((longStub                       &&
            couponInterval->prd_typ != 'I' &&
            couponInterval->prd_typ != 'E' &&
            couponInterval->prd_typ != 'J') &&
            *stubInfo > 0 &&
            nDates > 2 )
       {
           dateArr[nDates-2] = dateArr[nDates-1];
           --nDates;
       }
   }

   /* copy dates into the TDateList.
      Remember that they might be in reverse order */

   dateList = JpmcdsNewEmptyDateList(numPeriods);
   if (dateList == NULL)
       goto done;

   dateList->fNumItems = nDates;

   for (i=0 ;i<nDates ; ++i )
   {
       if (direction == 1)
       {
           index = i;
       }
       else
       {
           index = nDates-i-1;
       }
       dateList->fArray[i] = dateArr[index];
   }

   *dateListPtr = dateList;
   status = SUCCESS;

 done:
   if (status == FAILURE)
       JpmcdsErrMsg("%s: Failed.\n", routine);
   FREE(dateArr);

   return status;
}
