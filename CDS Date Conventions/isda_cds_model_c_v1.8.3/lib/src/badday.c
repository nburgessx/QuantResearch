/*
 * ISDA CDS Standard Model
 *
 * Copyright (C) 2009 International Swaps and Derivatives Association, Inc.
 * Developed and supported in collaboration with Markit
 * 
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the ISDA CDS Standard Model Public License.
 */

#include "cerror.h"
#include "convert.h"
#include "bastypes.h"
#include "date_sup.h"
#include "cmemory.h"
#include "badday.h"
#include "macros.h"
#include <stdlib.h>


/*
***************************************************************************
** This function converts a bad day to a good one (if it's in list).
***************************************************************************
*/
TDate JpmcdsBad2GoodBadDayList(
   TBadDayList *bdl,                  /* (I) bad day list */
   TDate        d)                     /* (I) bad date to convert */
{
   int i;

   if (bdl==NULL)                     /* allow NULL to ignore list */
   {
       return d;
   }

   for (i=0; i<bdl->count; i++)
   {
       if (bdl->badDay[i]==d) 
       {
           return bdl->goodDay[i];
       }
   }
   return d;
}


/*
***************************************************************************
** This function converts a good day to a bad one (if it's in list).
***************************************************************************
*/
TDate JpmcdsGood2BadBadDayList(
   TBadDayList *bdl,                  /* (I) bad day list */
   TDate        d)                     /* (I) bad date to convert */
{
   int i;

   if (bdl==NULL)                     /* allow NULL to ignore list */
   {
       return d;
   }

   for (i=0; i<bdl->count; i++)
   {
       if (bdl->goodDay[i]==d) 
       {
           return bdl->badDay[i];
       }
   }
   return d;
}


/*
***************************************************************************
** This function converts bad dates in a date list to good ones based on 
** a BadDayList.  DateList passed in is MODIFIED to contain "good" dates.
***************************************************************************
*/
void JpmcdsFixBadDLBadDayList(
   TBadDayList *bdl,                  /* (I) bad day list */
   TDateList   *dl)                    /* (M) date list to convert */
{
   int i;

   if (dl==NULL)
   {
       return;
   }

   for (i=0; i<dl->fNumItems; i++)
   {
      dl->fArray[i] = JpmcdsBad2GoodBadDayList(bdl, dl->fArray[i]);
   }
}

