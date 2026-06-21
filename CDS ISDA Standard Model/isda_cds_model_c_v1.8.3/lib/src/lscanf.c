/*
 * ISDA CDS Standard Model
 *
 * Copyright (C) 2009 International Swaps and Derivatives Association, Inc.
 * Developed and supported in collaboration with Markit
 * 
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the ISDA CDS Standard Model Public License.
 */

#include <stdio.h>
#include <stdarg.h>
#include <limits.h>
#include <ctype.h>
#include <math.h>
#include <errno.h>
#include <memory.h>
#include "cgeneral.h"
#include "cerror.h"
#include "cfileio.h"
#include "macros.h"


/*static void clearset(int);*/
static int testbit(int);

#define FFFLOAT      1
#define FFDOUBLE     2 
#define FFLONGDOUBLE    3

struct format {
    unsigned       suppress : 1;
    unsigned    haveWidth : 1;
    unsigned       fetchBase : 1;
    unsigned    negate : 1;
    unsigned       valid : 1;
    unsigned       unsigned_ : 1;
    unsigned    floating : 1;
    unsigned    dot : 1;
    unsigned       hSize : 1;
    unsigned       lSize : 1;
    unsigned       LSize : 1;
    int            fieldWidth;
};

#define JPMCDS_SIGDIGLEN 15
#define JPMCDS_SIGDIGLENOFFSET JPMCDS_SIGDIGLEN + 10

/*
 *  Has unusual behavior if the sig field is just JPMCDS_SIGDIGLEN.
 */
struct decrec {
    char        sgn;
    short       exp;
/*
 *  Added for float in part to work around asm code.
 */
    double         sigD;
    char        sig[JPMCDS_SIGDIGLENOFFSET];
};

static void ldtof(short, struct decrec *, void *);

#define TOP_TABLE 30
#define BOTTOM_TABLE -30

static double exptable[] = {
    1.0e-30,
    1.0e-29,
    1.0e-28,
    1.0e-27,
    1.0e-26,
    1.0e-25,
    1.0e-24,
    1.0e-23,
    1.0e-22,
    1.0e-21,
    1.0e-20,
    1.0e-19,
    1.0e-18,
    1.0e-17,
    1.0e-16,
    1.0e-15,
    1.0e-14,
    1.0e-13,
    1.0e-12,
    1.0e-11,
    1.0e-10,
    1.0e-9,
    1.0e-8,
    1.0e-7,
    1.0e-6,
    1.0e-5,
    1.0e-4,
    1.0e-3,
    1.0e-2,
    1.0e-1,
    1.0,
    1.0e1,
    1.0e2,
    1.0e3,
    1.0e4,
    1.0e5,
    1.0e6,
    1.0e7,
    1.0e8,
    1.0e9,
    1.0e10,
    1.0e11,
    1.0e12,
    1.0e13,
    1.0e14,
    1.0e15,
    1.0e16,
    1.0e17,
    1.0e18,
    1.0e19,
    1.0e20,
    1.0e21,
    1.0e22,
    1.0e23,
    1.0e24,
    1.0e25,
    1.0e26,
    1.0e27,
    1.0e28,
    1.0e29,
    1.0e30
};


/*
***************************************************************************
**  Local routine since acc version core dumps for no apparent reason.
**  Also this approach is faster.
***************************************************************************
*/
double JpmcdsLexp10(double dexp)
{
    int   exp;
    double   v = 1.0;

/*
 *  The looping may introduce a slight error value.
 *  For example, 10-15 might be 1.0000000000000000009.
 */
    exp = (int)dexp;
    if (exp <= TOP_TABLE) {
        if (exp >= BOTTOM_TABLE) {
            v = exptable[-1*BOTTOM_TABLE + exp];
        } else {
            v = pow((double)10,(double)dexp);
        }
    }  else {
            v = pow((double)10,(double)dexp);
    }
 
    return (v);
}



static int localGetc(TFile *tFile)
{
    int     c;

    if (tFile->hasLastChar == TRUE) 
    {
        c = tFile->lastChar;
        tFile->lastChar = 0;
        tFile->hasLastChar = FALSE;
        return (c);
    }

    if (tFile->type == TFILE_STRING) 
    {
        c = *(tFile->charPtr);
        (tFile->charPtr)++;
    } 
    else 
    {
        c = JpmcdsFgetc(tFile);
    }

    return (c);
}

static int localUngetc(TFile *tFile, int c)
{
    if (tFile->hasLastChar == TRUE) 
    {
        JpmcdsErrMsg("Unexpected ungetchar.\n");
    }

    tFile->lastChar = c;
    tFile->hasLastChar = TRUE;

    return 0;
}


static int localIsspace(int c)
{
    return isspace(c);
}


int JpmcdsLvfscanf(TFile *tFile, char *fmt, va_list arg)
{
    int nassigned = 0, nconverted = 0, nread = 0;
    int overflow, state;
    register short c, base = 0, digit;
    register long result = 0;
    register char *s = NULL;
    struct format F;
    struct decrec D;

    for (c = *fmt; c; c = *++fmt) {
        if (c != '%')
            goto match1;
        memset(&F, 0, sizeof(F));
        
            /*  check for assignment-suppression flag  */
        
        c = *++fmt;
        if (c == '*') {
            F.suppress = TRUE;
            c = *++fmt;
        }
        
            /*  decode field width  */
        
        if (isdigit(c)) {
            F.haveWidth = TRUE;
            do {
                F.fieldWidth = (10 * F.fieldWidth) + (c - '0');
                c = *++fmt;
            } while (isdigit(c));
            if (F.fieldWidth <= 0)
                goto done;
        }
        
            /*  determine appropriate conversion  */
        
conv: switch (c) {
                
                /*  'h' size modifier  */
                
            case 'h':
                F.hSize = TRUE;
                c = *++fmt;
                goto conv;
                
                /*  'l' size modifier  */
                
            case 'l':
                F.lSize = TRUE;
                c = *++fmt;
                goto conv;
                        
                /*  '?' base modifier  */
                
            case '?':
                F.fetchBase = TRUE;
                c = *++fmt;
                goto conv;
        
                /*  decimal (signed)  */
                
            case 'd':
                base = 10;
                goto conv_signed;
                
                /*  integer (signed)  */
                
            case 'i':
                base = 0;
                goto conv_signed;
                
                /*  octal (unsigned)  */
            
            case 'o':
                base = 8;
                goto conv_unsigned;
                
                /*  decimal (unsigned)  */
            
            case 'u':
                base = 10;
                goto conv_unsigned;
                
                /*  hexadecimal (unsigned)  */
            
            case 'p':
                F.lSize = TRUE;
                /* ... */
            case 'x':
            case 'X':
                base = 16;
                goto conv_unsigned;

                /*  floating point  */

            case 'e':
            case 'E':
            case 'f':
            case 'g':
            case 'G':
                F.floating = TRUE;
                state = -1;
                goto scan;

                /*  string  */
            
            case 's':
                do {
                    c = localGetc(tFile), ++nread;
                } while (localIsspace((char)c));
/*  No-op - only fixed scan set.
                clearset(1);
 */
                goto string;
                
                /*  character  */
            
            case 'c':
                if (!F.haveWidth)
                    F.fieldWidth = 1;
                if (!F.suppress)
                    s = va_arg(arg, char *);
                while (F.fieldWidth-- > 0) {
                    if ((c = localGetc(tFile)) == EOF)
                        goto done;
                    if (!F.suppress)
                        *s++ = (char)c;
                    ++nread;
                }
                if (!F.suppress)
                    ++nassigned;
                ++nconverted;
                continue;

                /*  store # bytes read so far  */
                
            case 'n':
                result = nread;
                if (!F.suppress)
                    --nassigned;
                goto assign;
                
                /*  others  */
            
            default:
                if (c != '%')
                    goto done;
            match1:
                if (localIsspace((char)c)) {
                    do {
                        c = localGetc(tFile), ++nread;
                    } while (localIsspace((char)c));
                    localUngetc(tFile, c), --nread;
                }
                else {
                    if ((c = localGetc(tFile)) != (char) *fmt) {
                        localUngetc(tFile, c);
                        goto done;
                    }
                    ++nread;
                }
                continue;
        }
        /*NOTREACHED*/

            /*  string scanner  */
            
string:
        if (!F.haveWidth)
            F.fieldWidth = INT_MAX;
        if (!F.suppress)
            s = va_arg(arg, char *);
        for (; c != EOF; c = localGetc(tFile), ++nread) {
            --F.fieldWidth;
            if (!testbit(c))
                break;
            F.valid = TRUE;
            if (!F.suppress)
                *s++ = (char)c;
            if (F.fieldWidth == 0)
                goto endstring;
        }
        localUngetc(tFile, c), --nread;
endstring:
        if (!F.valid)
            goto done;
        if (!F.suppress) {
            *s = 0;
            ++nassigned;
        }
        ++nconverted;
        continue;
        
            /*  integer conversions  */
    
conv_unsigned:
        F.unsigned_ = TRUE;
conv_signed:
        if (F.fetchBase)
            base = (short)va_arg(arg, int);
        state = 0;
        
            /*  numeric scanner  */
        
scan: result = 0;
        do {
            c = localGetc(tFile), ++nread;
        } while (localIsspace((char)c));
        if (!F.haveWidth)
            F.fieldWidth = INT_MAX;
        overflow = FALSE;
        for (; c != EOF; c = localGetc(tFile), ++nread) {
            --F.fieldWidth;
            switch (state) {
                     
                /*  (integer) start state  */
                     
                case 0:
                    state = 1;
                    if (c == '-') {
                        F.negate = TRUE;
                        break;
                    }
                    if (c == '+')
                        break;
                    /* ... */
                    
                /*  (integer) look for leading '0'  */
                        
                case 1:
                    state = 3;
                    if (c == '0') {
                        F.valid = TRUE;
                        if (F.fieldWidth) {
                            if (base == 0) {
                                base = 8;
                                state = 2;
                            }
                            else if (base == 16)
                                state = 2;
                        }
                        break;
                    }
                    goto state3;
                    
                /*  (integer) look for leading '0x'  */
                
                case 2:
                    state = 3;
                    if (c == 'x' || c == 'X') {
                        base = 16;
                        F.valid = FALSE;
                        break;
                    }
                    /* ... */
                
                /*  (integer) process each digit  */
                
                case 3:
                state3:
                    digit = c;
                    if (digit >= '0' && digit <= '9')
                        digit -= '0';
                    else if (digit >= 'A' && digit <= 'Z')
                        digit -= 'A' - 10;
                    else if (digit >= 'a' && digit <= 'z')
                        digit -= 'a' - 10;
                    else
                        goto pushback;
                    if (base == 0)
                        base = 10;
                    if (digit >= base)
                        goto pushback;
/*
 * Instead of assembly code, multiply times the base.
 */
                if (F.lSize) {
                    result *= base;
                    result += digit;
                } else if (F.hSize) {
                    result *= base;
                    result += digit;
                } else {
                    result *= base;
                    result += digit;
                }
                    
                    F.valid = TRUE;
                    break;

                /*  (floating) start state  */
                
                case -1:
                    state = -2;
                    D.exp = 0;
                    D.sig[0] = 0;
                    D.sigD = 0.0;
                    if (c == '-') {
                        D.sgn = 1;
                        break;
                    }
                    D.sgn = 0;
                    if (c == '+')
                        break;
                    /* ... */
                
                /*  (floating) process each digit  */
                
                case -2:
                    if (c >= '0' && c <= '9') {
                        F.valid = TRUE;
                        if (c != '0' || D.sig[0]) {
                            if (D.sig[0] >= sizeof D.sig - 1) {
                                if (!F.dot)
                                    ++D.exp;
                                break;
                            }
                            D.sig[(int)(++D.sig[0])] = (char)c;
                            D.sigD *= 10;
                            D.sigD += (int)c - 
                                 (int)'0';
                        }
                        if (F.dot)
                            --D.exp;
                    }
                    else if (c == '.' && !F.dot)
                        F.dot = TRUE;
                    else if ((c == 'e' || c == 'E') && F.valid) {
                        base = 10;
                        F.valid = FALSE;
                        state = 0;
                    }
                    else
                        goto pushback;
                    break;
            }
            
                /*  get next character  */
            
            if (F.fieldWidth == 0)
                goto endscan;
        }
        
            /*  push back last character read  */

pushback:   
        localUngetc(tFile, c), --nread;
endscan:
        if (!F.valid)
            goto done;
            
            /*  set sign of result  */
            
        if (F.negate && result) {
            result = -result;
            if (F.unsigned_ || result > 0)
                overflow = TRUE;
        }
        else if (!F.unsigned_ && result < 0)
            overflow = TRUE;
            
            /*  see if result fits  */

/*
 * Can't do this.  Just set result to D.exp.
 */
        if (F.floating) {
            result += D.exp;
        }
        if (F.hSize) {
            if (F.unsigned_) {
                if ((long)(unsigned short) result != result)
                    overflow = TRUE;
            }
            else {
                if ((short) result != result)
                    overflow = TRUE;
            }
        }
        else if (!F.lSize) {
            if (F.unsigned_) {
                if ((long)(unsigned int) result != result)
                    overflow = TRUE;
            }
            else {
                if ((int) result != result)
                    overflow = TRUE;
            }
        }
        
            /*  report overflow  */
            
        if (overflow) {
            if (F.unsigned_)
                result = 0;
            else if (F.hSize || F.floating)
                result = SHRT_MIN;
            else if (F.lSize)
                result = LONG_MIN;
            else
                result = INT_MIN;
            if (!F.negate)
                result = ~result;
            if (!F.floating)
                errno = ERANGE;
        }
        
            /*  assign result  */
    
assign:
        if (!F.suppress) {
            s = va_arg(arg, char *);
            if (F.floating) {
                D.exp = (short)result;
                /*
                 *  Only support double. 
                 */
                if (F.LSize)
                    ldtof(FFLONGDOUBLE, &D, s);
                else if (F.lSize)
                    ldtof(FFDOUBLE, &D, s);
                else
                    ldtof(FFFLOAT, &D, s);
            }
            else
            if (F.lSize)
                * (long *) s = result;
            else if (F.hSize)
                * (short *) s = (short)result;
            else
                * (int *) s = (int)result;
            ++nassigned;
        }
        ++nconverted;
    }

        /*  all done!  */
        
done:
    if (nconverted == 0 && c == EOF)
        return EOF;

    return nassigned;
}


/* ---------- scanset primitives ---------- */

/*
 *  testbit - see if character is in scanset
 *
 */

static int
testbit(int c)
{
/*  Only test against string characters. */
char  d;

    d = (char)c;

    if (d == '\t') 
        return  (0);
    else if (d == '\n')
        return  (0);
    else if (d == '\v') 
        return  (0);
    else if (d == '\f')
        return  (0);
    else if (d == '\r')
        return  (0);
    else if (d == ' ')
        return  (0);
    else if (d == '\0')
        return (0);

    return (1);
}


/* ---------- floating-point conversion ---------- */


/*
 *  Wrote a simple version of this routine.
 */
static void
ldtof(short code, struct decrec *d, void *p)
{
    double      t = (double)d->sigD;

    if (d->sgn == 1) {
        t *= -1;
    }

    switch(code) {
    case FFFLOAT:
        *(float *)p = (float)( JpmcdsLexp10((double)d->exp) * t);
        break;
    
    case FFDOUBLE:
        *(double *)p = JpmcdsLexp10((double)d->exp) * t;
        break;

    }
}
