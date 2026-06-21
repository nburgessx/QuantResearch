/*
** Microsoft framework functions for interfacing with Excel.
*/

#include "cgeneral.h"

/*
** Total amount of memory to allocate for all temporary XLOPERs
*/

#define MEMORYSIZE 10240


/* 
** Function prototypes
*/
LPSTR    GetTempMemory(int cBytes);
void     FreeAllTempMemory(void);
int      Excel(int xlfn, LPXLOPER pxResult, int count, ...);
LPXLOPER TempNum(double d);
LPXLOPER TempStr(LPSTR lpstr);
LPXLOPER TempBool(int b);
LPXLOPER TempInt(short int i);
LPXLOPER TempErr(WORD i);
LPXLOPER TempMissing(void);
void     InitFramework(void);
