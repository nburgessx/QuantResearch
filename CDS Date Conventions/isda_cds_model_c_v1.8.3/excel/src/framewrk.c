/*
**  Microsoft Excel Developer's Kit
**
**  File:           SAMPLE\FRAMEWRK\FRAMEWRK.C
**  Description:    Framework library for Microsoft Excel
**  Platform:       Microsoft Windows
**
**  This library provides some basic functions
**  that help you write Excel DLLs. It includes
**  simple functions for managing memory with XLOPERs,
**  creating temporary XLOPERs, robustly calling
**  Excel4(), and printing debugging strings on
**  a terminal attached to COM1.
**
**  The main purpose of this library is to help
**  you to write cleaner C code for calling Excel.
**  For example, using the framework library you
**  can write
**
**      Excel(xlcDisplay, 0, 2, TempMissing(), TempBool(0));
**
**  instead of the more verbose
**
**      XLOPER xMissing, xBool;
**      xMissing.xltype = xltypeMissing;
**      xBool.xltype = xltypeBool;
**      xBool.val.boolean = 0;
**      Excel4(xlcDisplay, 0, 2, (LPXLOPER) &xMissing, (LPXLOPER) &xBool);
**
**
**  The library is non-reentrant: it assumes that
**  if there are multiple copies of Excel using the
**  DLL, they will not be preempted. This is
**  acceptable under Windows/DOS, but may present
**  some problems under Windows/NT. In particular,
**  the function Excel() frees all the temporary
**  memory, which means that this DLL must not be
**  reentered while temporary memory is in use.
**
**  Define DEBUG to use the debugging functions.
**
**  Source code is provided so that you may
**  enhance this library or optimize it for your
**  own application.
**
**  $Header$
*/


#include <windows.h>
#include "xlcall.h"
#include "framewrk.h"

#include <stdarg.h>

/*
** Globals (see the comment about reentrancy 
** earlier in this file)
*/

char vMemBlock[MEMORYSIZE]; /* Memory for temporary XLOPERs */
int vOffsetMemBlock=0;      /* Offset of next memory block to allocate */




/*
** GetTempMemory
**
** Allocates temporary memory. Temporary memory
** can only be freed in one chunk, by calling
** FreeAllTempMemory(). This is done by Excel().
**
** Arguments:
**
**      int cBytes      How many bytes to allocate
**
** Returns:
**
**      LPSTR           A pointer to the allocated memory,
**                      or 0 if more memory cannot be
**                      allocated. If this fails,
**                      check that you are initializing
**                      vOffsetMemBlock to 0, and check that
**                      MEMORYSIZE is big enough.
**
** Algorithm:
**
**      The memory allocation algorithm is extremely
**      simple: on each call, allocate the next cBytes
**      bytes of a static memory buffer. If the buffer
**      becomes too full, simply fail. To free memory,
**      simply reset the pointer (vOffsetMemBlock)
**      back to zero. This memory scheme is very fast
**      and is optimized for the assumption that the
**      only thing you are using temporary memory
**      for is to hold arguments while you call Excel().
**      We rely on the fact that you will free all the
**      temporary memory at the same time. We also
**      assume you will not need more memory than
**      the amount required to hold a few arguments
**      to Excel().
*/

LPSTR GetTempMemory(int cBytes)
{
    LPSTR lpMemory;

    if (vOffsetMemBlock + cBytes > MEMORYSIZE)
    {
        return 0;
    }
    else
    {
        lpMemory = (LPSTR) &vMemBlock + vOffsetMemBlock;
        vOffsetMemBlock += cBytes;

        /* Prevent odd pointers */
        if (vOffsetMemBlock & 1) vOffsetMemBlock++;
        return lpMemory;
    }
}


/*
** FreeAllTempMemory
**
** Frees all temporary memory that has been allocated.
**
** Arguments:
**
**      None.
**
** Return value:
**
**      None.
*/

void FreeAllTempMemory(void)
{
    vOffsetMemBlock = 0;
}


/*
** Excel
**
** A fancy wrapper for the Excel4() function. It also
** does the following:
**
**  (1) Checks that none of the LPXLOPER arguments are 0,
**      which would indicate that creating a temporary XLOPER
**      has failed. In this case, it doesn't call Excel
**      but it does print a debug message.
**  (2) If an error occurs while calling Excel,
**      print a useful debug message.
**  (3) When done, free all temporary memory.
**
**  #1 and #2 require DEBUG to be defined.
**
** Arguments (same as Excel4()):
**
**      int xlfn            Function number (xl...) to call
**      LPXLOPER pxResult   Pointer to a place to stuff the result,
**                          or 0 if you don't care about the result.
**      int count           Number of arguments
**      ...                 (all LPXLOPERs) - the arguments.
**
** Return value:
**
**      A return code (Some of the xlret... values, as defined
**      in XLCALL.H, OR'ed together).
**
** Note:
**
**      Be sure to cast all the arguments after the third
**      to LPXLOPERs. If you accidentally pass a near pointer
**      instead of a far pointer, you will probably crash Excel.
*/

int Excel(int xlfn, LPXLOPER pxResult, int count, ...)
{
    int xlret;
    LPXLOPER xArray[30]; /* Cannot have more than 30 arguments */
    int maxCount = sizeof(xArray) / sizeof(LPXLOPER);
    int pos;
    va_list ap;

/*
** In the published version of this function, this walking through the variable
** argument list was not done "properly". This is an attempt to do so...
*/
    va_start (ap, count);
    if (count > maxCount) count = maxCount;
    for (pos = 0; pos < count; ++pos)
    {
        xArray[pos] = va_arg(ap, LPXLOPER);
    }
    va_end (ap);

    xlret = Excel4v(xlfn,pxResult,count,xArray);

    FreeAllTempMemory();

    return xlret;
}


/*
** TempNum
**
** Creates a temporary numeric (IEEE floating point) XLOPER.
**
** Arguments:
**
**      double d        The value
**
** Returns:
**
**      LPXLOPER        The temporary XLOPER, or 0
**                      if GetTempMemory() failed.
**
*/
LPXLOPER TempNum(double d)
{
    LPXLOPER lpx;

    lpx = (LPXLOPER) GetTempMemory(sizeof(XLOPER));

    if (!lpx)
    {
        return 0;
    }

    lpx->xltype = xltypeNum;
    lpx->val.num = d;

    return lpx;
}


/*
** TempStr
**
** Creates a temporary string XLOPER.
**
** Arguments:
**
**      LPSTR lpstr     The string, as a null-terminated
**                      C string, with the first byte
**                      undefined. This function will
**                      count the bytes of the string
**                      and insert that count in the
**                      first byte of lpstr. Excel cannot
**                      handle strings longer than 255
**                      characters.
**
** Returns:
**
**      LPXLOPER        The temporary XLOPER, or 0
**                      if GetTempMemory() failed.
**
** Notes:
**
**      (1) This function has the side effect of inserting
**          the byte count as the first character of
**          the created string.
**
**      (2) For highest speed, with constant strings,
**          you may want to manually count the length of
**          the string before compiling, and then avoid
**          using this function.
**
**      (3) Behavior is undefined for non-null terminated
**          input or strings longer than 255 characters.
**
**      Modified as recommended in ALIB-95 for VC7
*/
LPXLOPER TempStr(LPSTR lpstr)
{
    LPXLOPER lpx;
    int chars;

    lpx = (LPXLOPER) GetTempMemory(sizeof(XLOPER));

    if (!lpx)
    {
        return 0;
    }

    chars = lstrlen(lpstr);
    if (chars > 255) chars = 255;

    lpx->val.str = (char*) GetTempMemory((sizeof(char) * chars) + 1);
    if (!lpx->val.str)
    {
        return 0;
    }
    strncpy(lpx->val.str, lpstr, chars);
    lpx->val.str[0] = (BYTE) chars;
    lpx->val.str[chars] = '\0';
    lpx->xltype = xltypeStr;

    return lpx;
}


/*
** TempBool
**
** Creates a temporary logical (true/false) XLOPER.
**
** Arguments:
**
**      int b           0 - for a FALSE XLOPER
**                      Anything else - for a TRUE XLOPER
**
** Returns:
**
**      LPXLOPER        The temporary XLOPER, or 0
**                      if GetTempMemory() failed.
**
*/
LPXLOPER TempBool(int b)
{
    LPXLOPER lpx;

    lpx = (LPXLOPER) GetTempMemory(sizeof(XLOPER));

    if (!lpx)
    {
        return 0;
    }

    lpx->xltype = xltypeBool;
    lpx->val.boolean = b?1:0;

    return lpx;
}


/*
** TempInt
**
** Creates a temporary integer XLOPER.
**
** Arguments:
**
**      short int i          The integer
**
** Returns:
**
**      LPXLOPER        The temporary XLOPER, or 0
**                      if GetTempMemory() failed.
**
*/
LPXLOPER TempInt(short int i)
{
    LPXLOPER lpx;

    lpx = (LPXLOPER) GetTempMemory(sizeof(XLOPER));

    if (!lpx)
    {
        return 0;
    }

    lpx->xltype = xltypeInt;
    lpx->val.w = i;

    return lpx;
}


/*
** TempErr
**
** Creates a temporary error XLOPER.
**
** Arguments:
**
**      WORD err        The error code. One of the xlerr...
**                      constants, as defined in XLCALL.H.
**                      See the Excel user manual for
**                      descriptions about the interpretation
**                      of various error codes.
**
** Returns:
**
**      LPXLOPER        The temporary XLOPER, or 0
**                      if GetTempMemory() failed.
**
*/
LPXLOPER TempErr(WORD err)
{
    LPXLOPER lpx;

    lpx = (LPXLOPER) GetTempMemory(sizeof(XLOPER));

    if (!lpx)
    {
        return 0;
    }

    lpx->xltype = xltypeErr;
    lpx->val.err = err;

    return lpx;
}



/*
** TempMissing
**
** This is used to simulate a missing argument when
** calling Excel(). It creates a temporary
** "missing" XLOPER.
**
** Arguments:
**
**      none.
**
** Returns:
**
**      LPXLOPER        The temporary XLOPER, or 0
**                      if GetTempMemory() failed.
**
*/
LPXLOPER TempMissing(void)
{
    LPXLOPER lpx;

    lpx = (LPXLOPER) GetTempMemory(sizeof(XLOPER));

    if (!lpx)
    {
        return 0;
    }

    lpx->xltype = xltypeMissing;

    return lpx;
}


/*
** InitFramework
**
** Initializes all the framework functions.
**
** Arguments:
**
**      None.
**
** Return value:
**
**      None.
*/
void InitFramework(void)
{
    vOffsetMemBlock = 0;
}


