/*
 * ISDA CDS Standard Model
 *
 * Copyright (C) 2009 International Swaps and Derivatives Association, Inc.
 * Developed and supported in collaboration with Markit
 * 
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the ISDA CDS Standard Model Public License.
 */

#include "buscache.h"
#include <ctype.h>
#include "cfileio.h"
#include "convert.h"
#include "cerror.h"
#include "dtlist.h"
#include "dateconv.h"
#include "macros.h"
#include "strutil.h"
#include <limits.h>

#ifdef LINUX
#define stricmp strcasecmp
#endif

/*---------------------------------------------------------------------------
 *                  MACROS used in this file.
 *---------------------------------------------------------------------------
 */

#define MAX_BUFFER                   (256)

/*
** Copy a string from source to target and upcase. The macro also
** does length-checking
*/
#define COPY_CAPITALIZE(pTarget,pSrc,len) \
if (pSrc != NULL && pTarget != NULL)\
{\
   char *  PSRC = (pSrc);\
   char *  PTARGET = (pTarget);\
   while ( (*PSRC) != 0  && (PTARGET - pTarget) <= (len) )\
   {\
          (*PTARGET) = (char)toupper((int)(*PSRC));\
          PSRC++;  PTARGET++;\
   }\
   (*PTARGET) = 0;\
}


/*t
 * Holiday cache information
 * */
typedef struct _THoliday
{
    char         *name;   /* Upper case name associated with this entry */
    THolidayList *hl;     /* Holiday list */
    void         *next;
} THoliday;


THoliday *cache = NULL;


/*
***************************************************************************
** Private function declarations
***************************************************************************
*/
static int verifyHolidayList
(THolidayList *hl /* (I/O) Holiday list */
);

/* Finds an entry associated with a name in the cache. */
static THoliday *holidayFind(char *name);

/* Delete holiday calendar from cache. */
static void holidayDelete(THoliday *hol);

/* Create holiday calendar. */
static THoliday* JpmcdsNewHoliday(THolidayList *hl, char *name);

/* Free memory used by holiday calendar. */
void JpmcdsFreeHoliday(THoliday *holiday);


/*
***************************************************************************
** This file is broken into two types of functions.
**
** 1. Holiday cache manipulation functions. These are functions that
**    manipulate the holiday cache but perform no significant analytics.
**
** 2. Holiday list manipulation functions. These are functions that
**    manipulate the holiday list structure but perform no significant
**    analytics.
***************************************************************************
*/


/*
***************************************************************************
** The function first checks if the input name is really an encoded
** pointer. If it is, the function decodes the pointer and returns
** it. Otherwise, it checks whether the name is "NONE" or "NO_WEEKENDS".
** If it is either of those, the function returns a pointer to an
** appropriate static data structure (for efficiency). If it is
** none of the above, then the function searches for the name in the
** cache. If it finds the name, it returns a pointer to the cached
** holiday list. Otherwise, it interprets the name as a file name
** and tries to load the holiday file from disk.
**
** You can use name = "NONE" or name = "NO_WEEKENDS" and the resulting
** holiday list will be defined appropriately.
**
** This behaviour exactly duplicates the behaviour when using a holiday
** name as an input to any analytics function.
**
** Returns NULL on failure, a valid THolidayList pointer on success.
***************************************************************************
*/
THolidayList* JpmcdsHolidayListFromCache
(char *name)    /* (I) Name associated with the holidays */
{
    static char   routine[] = "JpmcdsHolidayListFromCache";
    int           status    = FAILURE;

    THolidayList  *hl = NULL; /* To be returned */
    THoliday      *hol;

    /* Must have a name. */
    if (name == NULL)
    {
        JpmcdsErrMsg ("%s: NULL inputs.\n", routine);
        goto done;
    }

    /* Use holidayFind to find in the cache. Note that NONE and NO_WEEKENDS are always in the cache. */
    hol = holidayFind (name);
    if (hol != NULL)
    {
        hl = hol->hl;
    }
    else
    {
        /* Not in the cache, so use JpmcdsHolidayListRead to read from file. */
        hl = JpmcdsHolidayListRead (name);
        if (hl == NULL)
            goto done;

        /* Now add hl to the cache. This takes a shallow copy, but is
           guaranteed to be destructive if the operation fails. */
        if (JpmcdsHolidayListAddToCache (name, hl) != SUCCESS)
            goto done;
    }

    status = SUCCESS;

done:

    if (status != SUCCESS)
    {
        JpmcdsErrMsg ("%s: Failed.\n", routine);
        return NULL;
    }

    return hl;
}

/*
***************************************************************************
** Adds a holiday list to the holiday cache. If the entry already exists
** in the cache, then the old version will be deleted.
***************************************************************************
*/
int JpmcdsHolidayListAddToCache
(char         *name,   /* (I) Name to associate holidays with */
 THolidayList *hl      /* (I) Adds shallow copy */
)
{
    static char routine[] = "JpmcdsHolidayListAddToCache";
    int         status    = FAILURE;

    THoliday   *hol;
    THoliday   *oldHol;

    /* JpmcdsNewHoliday creates a holiday structure with the holiday list
       as a shallow copy, but a deep copy of the capitalized name.
       This shallow copy is for efficiency. We do this call first so
       that we can be sure that either hl is owned by hol or hl is
       deleted. */
    hol = JpmcdsNewHoliday(hl, name);
    if (hol == NULL)
        goto done;

    /* If necessary get rid of the old entry with the same name. 
       Use hol->name since this has been capitalized. */
    oldHol = holidayFind(hol->name);
    if (oldHol != NULL)
    {
        /*
        ** Cannot delete NONE or NO_WEEKENDS at this point. This is a sign
        ** that somebody is trying to overwrite these standard names.
        */
        if (hol->name[0] == 'N')
        {
            if (strcmp(hol->name, "NONE") == 0 || strcmp (hol->name, "NO_WEEKENDS") == 0)
            {
                JpmcdsErrMsg ("%s: Attempt to over-write standard holiday %s\n", routine, hol->name);
                goto done;
            }
        }

        holidayDelete(oldHol);
    }

    /* Insert into cache */
    hol->next = cache;
    cache = hol;
    hol = NULL; /* Now owned by cache */
    status = SUCCESS;

done:
    JpmcdsFreeHoliday (hol);

    if (status != SUCCESS)
        JpmcdsErrMsg ("%s: Failed.\n", routine);

    return status;
}


/*f
***************************************************************************
** Empty holiday cache
***************************************************************************
*/
void JpmcdsHolidayEmptyCache (void)
{
    THoliday *node;
    THoliday *next;

    if(cache != NULL)
    {
        node = cache;
        next = node->next;
        while (node)
        {
            JpmcdsFreeHoliday(node);
            node = next;
            if(node) next = node->next;
        }
        cache = NULL;
    }
}


/*
***************************************************************************
** Finds an entry associated with a name in the cache.
** Converts the name to upper case before searching in the hash table.
***************************************************************************
*/
static THoliday *holidayFind(char *name)
{
    THoliday *hol = NULL;

    /* ensure that built-in calendars are present in the cache */
    if (cache == NULL)
    {
        THoliday *h;

        h = JpmcdsNewHoliday(JpmcdsHolidayListNewGeneral(NULL, JPMCDS_WEEKEND_SATURDAY | JPMCDS_WEEKEND_SUNDAY), "NONE");
        if (h == NULL)
            return NULL;

        cache = h;
        
        h = JpmcdsNewHoliday(JpmcdsHolidayListNewGeneral(NULL, JPMCDS_WEEKEND_NO_WEEKENDS), "NO_WEEKENDS");
        if (h == NULL)
            return NULL;
        
        cache->next = h;
        h->next = NULL;
    }
    
    hol = cache;
    if (name != NULL)
    {
        while (hol)
        {
            if (stricmp(name, hol->name) == 0)
                return hol;
            
            hol = hol->next;
        }
    }

    /* not necessarily an error */
    return NULL;
}


/*
***************************************************************************
** Deletes an entry in the holiday cache.
***************************************************************************
*/
static void holidayDelete(THoliday *hol)
{
    THoliday *node;
    
    /* remove node from cache list */
    if (cache == hol)
        cache = hol->next;
    else
    {
        node = cache;
        while (node)
        {
            if (node->next == hol)
                node->next = hol->next;
            else
                node = node->next;
        }
    }

    /* Free memory. */
    JpmcdsFreeHoliday(hol);
}


/*
***************************************************************************
** Creates a new THoliday pointer. Takes shallow copy of the holiday list,
** takes a deep capitalized copy of the name.
***************************************************************************
*/
static THoliday* JpmcdsNewHoliday
(THolidayList *hl,   /* (I) Holiday list - takes shallow copy */
 char         *name  /* (I) Name - takes deep capitalized copy */
)
{
    static char routine[] = "JpmcdsNewHoliday";
    int         status    = FAILURE;

    THoliday   *hol = NULL;

    if (hl == NULL || name == NULL)
    {
        JpmcdsErrMsg ("%s: NULL inputs.\n", routine);
        goto done;
    }

    hol = NEW(THoliday);
    if (hol == NULL)
        goto done;

    hol->hl   = hl;
    hol->name = JpmcdsStringDuplicate(name);
    hl        = NULL; /* Now owned by hol */
    if (hol->name == NULL)
        goto done;

    hol->next = NULL;
    status = SUCCESS;

done:
    JpmcdsHolidayListDelete(hl);
    if (status != SUCCESS)
    {
        JpmcdsFreeHoliday (hol);
        hol = NULL;
        JpmcdsErrMsg ("%s: Failed.\n", routine);
    }

    return hol;
}

/*
***************************************************************************
** Deletes a THoliday pointer. Possibly needs to be exported since it will
** be used as a call-back function.
***************************************************************************
*/
void JpmcdsFreeHoliday
(THoliday *holiday)
{
    if (holiday != NULL)
    {
        JpmcdsHolidayListDelete(holiday->hl);
        FREE_ARRAY(holiday->name);
        FREE(holiday);
        holiday = NULL;
    }
}




/*
***************************************************************************
** 3. Holiday list manipulation functions. These are functions that
**    manipulate the holiday list structure but perform no significant
**    analytics.
***************************************************************************
*/

/*
***************************************************************************
** Creates a new holiday structure.
**
** The date list can be NULL, in which case the resulting date list in the
** holiday structure will be a date list with no dates, e.g.
**      hl->dateList->fNumItems = 0;
***************************************************************************
*/
THolidayList* JpmcdsHolidayListNewGeneral
(TDateList *dateList,      /* (I) Date list to use */
 long       weekends       /* (I) Weekends flag - use JPMCDS_WEEKEND_... */
)
{
    static char   routine[] = "JpmcdsHolidayListNewGeneral";
    int           status = FAILURE;

    THolidayList *hl = NULL;
    TDateList    *dl = NULL;

    /* get new date list */
    if (dateList == NULL)
    {
        dl = JpmcdsNewEmptyDateList (0);
    }
    else
    {
        dl = JpmcdsCopyDateList (dateList);
    }
    if (dl == NULL)
    {
        goto done;
    }

    /* get new holiday list */
    hl = NEW(THolidayList);
    if (hl == NULL)
    {
        goto done;
    }

    /* fill in holiday list */
    hl->dateList     = dl;
    hl->weekends     = weekends;
    dl               = NULL; /* Now owned by hl */

    if (verifyHolidayList (hl) != SUCCESS)
        goto done;

    status = SUCCESS;

done:

    JpmcdsFreeDateList (dl);

    if (status != SUCCESS)
    {
        JpmcdsHolidayListDelete (hl);
        hl = NULL;
        JpmcdsErrMsg ("%s: Failed.\n", routine);
    }

    return hl;
}


/*
***************************************************************************
** Deletes a holiday list structure.
***************************************************************************
*/
void JpmcdsHolidayListDelete
(THolidayList *hl)  /* (I) Holiday list to delete */
{
    if (hl != NULL)
    {
        JpmcdsFreeDateList (hl->dateList);
        FREE(hl);
    }
}


/*
***************************************************************************
** Reads a holiday file into memory as a holiday list.
**
** Structure of holiday file, ascii text file of lines like:
**     #            - commment (blank lines are also skipped)
**     19631026     - means 10/26/1963 is a holiday
**     # SATURDAY_NOT_ALWAYS_HOLIDAY  - sets "saturday isn't always a holiday"
**     # SUNDAY_NOT_ALWAYS_HOLIDAY    - sets "sunday isn't always a holiday"
**     # MONDAY_ALWAYS_HOLIDAY        - sets "monday as always a holiday"
**     # TUESDAY_ALWAYS_HOLIDAY       - sets "tuesday as always a holiday"
**     # WEDNESDAY_ALWAYS_HOLIDAY     - sets "wednesday as always a holiday"
**     # THURDSAY_ALWAYS_HOLIDAY      - sets "thursday as always a holiday"
**     # FRIDAY_ALWAYS_HOLIDAY        - sets "friday as always a holiday"
**
** Dates must be in increasing order.
***************************************************************************
*/
THolidayList* JpmcdsHolidayListRead
(char *fileName   /* (I) Name of file to read (may differ) */
)
{
    static char routine[] = "JpmcdsHolidayListRead";
    int         status = FAILURE;

    int        numHols;             /* count of holidays */
    int        idx;                 /* counts over holidays */
    char       buffer[MAX_BUFFER];  /* line of text from file */
    char       ucBuffer[MAX_BUFFER];/* upper case line of text from file */
    TFile     *fp = NULL;           /* file handle */
    THolidayList *hl = NULL;        /* to be returned */
    TDateList *dl = NULL;          /* date list read in */
    long       weekends;

    fp = JpmcdsFopen (fileName, JPMCDS_FREAD);
    if (fp == NULL)
    {
        JpmcdsErrMsg ("%s: Couldn't open file %s.\n", routine, fileName);
        goto done;
    }

    /* count number of holidays specified */
    numHols = 0;
    while (JpmcdsFgets (buffer, MAX_BUFFER, fp) == SUCCESS)
    {
        if (atol (buffer) > 16010101)
            numHols++;
    }

    if (JpmcdsFclose (fp) == FAILURE)  /* close file, reopen to read again */
    {
        fp = NULL;
        goto done;
    }
    fp = NULL;

    fp = JpmcdsFopen (fileName, JPMCDS_FREAD);
    if (fp == NULL)
    {
        JpmcdsErrMsg ("%s: Couldn't open file %s twice.\n", routine, fileName);
        goto done;
    }

    dl = JpmcdsNewEmptyDateList (numHols);

    if (dl == NULL)
        goto done;

    idx = 0;
    weekends = JPMCDS_WEEKEND_STANDARD; /* Includes SAT and SUN */

    while (JpmcdsFgets (buffer, MAX_BUFFER, fp) == SUCCESS)
    {
        if (buffer[0] == '#')
        {
            static char monString[] = "# MONDAY_ALWAYS_HOLIDAY";
            static char tueString[] = "# TUESDAY_ALWAYS_HOLIDAY";
            static char wedString[] = "# WEDNESDAY_ALWAYS_HOLIDAY";
            static char thuString[] = "# THURSDAY_ALWAYS_HOLIDAY";
            static char friString[] = "# FRIDAY_ALWAYS_HOLIDAY";
            static char satString[] = "# SATURDAY_NOT_ALWAYS_HOLIDAY";
            static char sunString[] = "# SUNDAY_NOT_ALWAYS_HOLIDAY";

            COPY_CAPITALIZE(ucBuffer, buffer, MAX_BUFFER-1);

            if (strncmp (buffer, monString, sizeof(monString)-1) == 0)
                weekends |= JPMCDS_WEEKEND_MONDAY;

            if (strncmp (buffer, tueString, sizeof(tueString)-1) == 0)
                weekends |= JPMCDS_WEEKEND_TUESDAY;

            if (strncmp (buffer, wedString, sizeof(wedString)-1) == 0)
                weekends |= JPMCDS_WEEKEND_WEDNESDAY;

            if (strncmp (buffer, thuString, sizeof(thuString)-1) == 0)
                weekends |= JPMCDS_WEEKEND_THURSDAY;

            if (strncmp (buffer, friString, sizeof(friString)-1) == 0)
                weekends |= JPMCDS_WEEKEND_FRIDAY;

            if (strncmp (buffer, satString, sizeof(satString)-1) == 0)
                weekends &= ~JPMCDS_WEEKEND_SATURDAY;

            if (strncmp (buffer, sunString, sizeof(sunString)-1) == 0)
                weekends &= ~JPMCDS_WEEKEND_SUNDAY;
        }
        else
        {
            long dateVal = atol (buffer);
            if (dateVal > 16010101)
            {
                TMonthDayYear   mdy;
                TDate           thisDate;

                if (idx >= numHols)
                {
                    JpmcdsErrMsg ("%s: More dates on second scan of file %s.\n", routine, fileName);
                    goto done;
                }

                mdy.year  = dateVal / 10000;
                mdy.month = (dateVal - mdy.year * 10000) / 100;
                mdy.day   = dateVal % 100;
                if (JpmcdsMDYToDate (&mdy, &thisDate) == FAILURE)
                {
                    JpmcdsErrMsg ("%s: invalid date: %s", routine, buffer);
                    goto done;
                }
                if (idx > 0 && dl->fArray[idx-1] >= thisDate)
                {
                    JpmcdsErrMsg ("%s: Date out of order: %s", routine, buffer);
                    goto done;
                }
                dl->fArray[idx] = thisDate;
                ++idx;
            }
        }
    }

    if (numHols < 1 && weekends == JPMCDS_WEEKEND_STANDARD)
    {
        JpmcdsErrMsg ("%s: No holiday information found in %s.\n", routine, fileName);
        JpmcdsErrMsg ("   Either week-end information or dates must be provided.\n");
        goto done;
    }

    if (idx != numHols)
    {
        JpmcdsErrMsg ("%s: Less dates on second scan of file %s\n", routine, fileName);
        goto done;
    }

    if (JpmcdsFclose (fp) == FAILURE)
    {
        fp = NULL;
        goto done;
    }
    fp = NULL;

    hl = JpmcdsHolidayListNewGeneral (dl, weekends);
    if (hl == NULL)
        goto done;

    status = SUCCESS;

done:

    JpmcdsFclose (fp);
    JpmcdsFreeDateList (dl);
    if (status != SUCCESS)
        JpmcdsErrMsg ("%s: Failed.\n", routine);

    return hl;
}


/*
***************************************************************************
** This function verifies a holiday list. This means the following actions
** are performed:
** 1. Validates that dates are in strictly increasing order.
** 2. Removes dates which are already defined by being a week-end.
**
** Currently does not perform stage 2 since this causes the buscache test
** driver to fail. We will consider re-introducing this feature after the
** release of version 9.4.
***************************************************************************
*/
static int verifyHolidayList
(THolidayList *hl /* (I/O) Holiday list */
)
{
#undef CHECK_FOR_WEEKENDS

    static char routine[] = "verifyHolidayList";
    int         status    = FAILURE;

    long        idx;
    long        count;
    TDate       lastDate = 0;

    count = 0;

    for (idx = 0; idx < hl->dateList->fNumItems; ++idx)
    {
        TDate thisDate = hl->dateList->fArray[idx];

        if (thisDate <= lastDate)
        {
            JpmcdsErrMsg ("%s: Dates are not in strictly increasing order.\n", routine);
            JpmcdsErrMsg ("    [%ld] = %s; LastDate = %s.\n",
                       idx, JpmcdsFormatDate(thisDate), JpmcdsFormatDate(lastDate));
            goto done;
        }

        lastDate = thisDate;

        /* remove weekends */
        if (JPMCDS_IS_WEEKDAY(thisDate, hl->weekends))
        {
            hl->dateList->fArray[count] = thisDate;
            ++count;
        }
    }

    hl->dateList->fNumItems = count;

    status = SUCCESS;

done:

    if (status != SUCCESS)
        JpmcdsErrMsg ("%s: Failed.\n", routine);

    return status;
}
