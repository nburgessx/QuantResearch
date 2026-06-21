/*
 * ISDA CDS Standard Model
 *
 * Copyright (C) 2009 International Swaps and Derivatives Association, Inc.
 * Developed and supported in collaboration with Markit
 * 
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the ISDA CDS Standard Model Public License.
 */

#include "lprintf.h"
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <memory.h> 
#include <math.h>
#include "cgeneral.h"
#include "cfileio.h"
#include "cerror.h"
#include "macros.h"

#define BUFLEN 512

struct format {
    unsigned       leftJustify : 1;
    unsigned       forceSign : 1;
    unsigned       altForm : 1;
    unsigned       zeroPad : 1;
    unsigned       havePrecision : 1;
    unsigned       hSize : 1;
    unsigned       lSize : 1;
    unsigned       LSize : 1;
    char        sign;
    char        exponent;
    int            fieldWidth;
    int            precision;
};


/* default output precision for floating numbers */
#define JPMCDS_DEF_PRECISION 6



int JpmcdsLocalPutc(char c, TFile *tFile)
{
    if (tFile->type == TFILE_STRING) {
        *(tFile->charPtr) = c;
        (tFile->curSize)++;
        if (tFile->curSize == tFile->size) {
            return (-1);
        }
        (tFile->charPtr)++;
    } else {
        if (JpmcdsFputc(c, tFile) == FAILURE) {
            return (-1);
        }
    }

    return (c);
}

int JpmcdsLocalFwrite(char *ptr, int nitems, TFile *tFile)
{
    int   i;
    char  *tptr;

    if (tFile->type == TFILE_STRING) {
        tptr = ptr;
        for (i = 0; i < nitems; i++) {
            *(tFile->charPtr) = *tptr;
            (tFile->charPtr)++;
            tptr++;
            (tFile->curSize)++;
            if (tFile->curSize == tFile->size) {
                return (-1);
            }
        }
    } else {
        if (JpmcdsFwrite(ptr, nitems, tFile) == FAILURE) {
            return (-1);
        }
    }

    return (nitems);
}


/*
*******************************************************************************
** vprintf functionality using a TFile
*******************************************************************************
*/
int JpmcdsVfprintf(TFile *tFile, char *fmt, va_list arg)
{
    register int c, i, nwritten = 0;
    register unsigned long n;
    double x;
    register char *s;
    char buf[BUFLEN], *digits, *t;
    char buf2[BUFLEN];
    struct format F;
    
    for (c = *fmt; c; c = *++fmt) {
        if (c != '%')
            goto copy1;
        memset(&F, 0, sizeof(F));
        
        /*  decode flags  */
            
        for (;;) {
            c = *++fmt;
            if (c == '-')
                F.leftJustify = TRUE;
            else if (c == '+')
                F.forceSign = TRUE;
            else if (c == ' ')
                F.sign = ' ';
            else if (c == '#')
                F.altForm = TRUE;
            else if (c == '0')
                F.zeroPad = TRUE;
            else
                break;
        }
        
            /*  decode field width  */
            
        if (c == '*') {
            if ((F.fieldWidth = va_arg(arg, int)) < 0) {
                F.leftJustify = TRUE;
                F.fieldWidth = -F.fieldWidth;
            }
            c = *++fmt;
        }
        else {
            for (; c >= '0' && c <= '9'; c = *++fmt)
                F.fieldWidth = (10 * F.fieldWidth) + (c - '0');
        }
        
            /*  decode precision  */
            
        if (c == '.') {
            if ((c = *++fmt) == '*') {
                F.precision = va_arg(arg, int);
                c = *++fmt;
            }
            else {
                for (; c >= '0' && c <= '9'; c = *++fmt)
                    F.precision = (10 * F.precision) + (c - '0');
            }
            if (F.precision >= 0)
                F.havePrecision = TRUE;
        }
        
            /*  perform appropriate conversion  */
        
        s = &buf[BUFLEN];
        if (F.leftJustify)
            F.zeroPad = FALSE;
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
                        
                /*  decimal (signed)  */
                
            case 'd':
            case 'i':
                if (F.lSize)
                    n = va_arg(arg, long);
                else
                    n = va_arg(arg, int);
                if (F.hSize)
                    n = (short) n;
                if ((long) n < 0) {
                    n = -((long)n);
                    F.sign = '-';
                }
                else if (F.forceSign)
                    F.sign = '+';
                goto decimal;
                
                /*  decimal (unsigned)  */
                
            case 'u':
                if (F.lSize)
                    n = va_arg(arg, unsigned long);
                else
                    n = va_arg(arg, unsigned int);
                if (F.hSize)
                    n = (unsigned short) n;
                F.sign = 0;
                goto decimal;
                
                /*  decimal (common code)  */

            decimal:
                if (!F.havePrecision) {
                    if (F.zeroPad) {
                        F.precision = F.fieldWidth;
                        if (F.sign)
                            --F.precision;
                    }
                    if (F.precision < 1)
                        F.precision = 1;
                }
                for (i = 0; n; n /= 10, i++)
                    *--s = (char)(n % 10 + '0');
                for (; i < F.precision; i++)
                    *--s = '0';
                if (F.sign) {
                    *--s = F.sign;
                    i++;
                }
                break;
                
                /*  octal (unsigned)  */
                
            case 'o':
                if (F.lSize)
                    n = va_arg(arg, unsigned long);
                else
                    n = va_arg(arg, unsigned int);
                if (F.hSize)
                    n = (unsigned short) n;
                if (!F.havePrecision) {
                    if (F.zeroPad)
                        F.precision = F.fieldWidth;
                    if (F.precision < 1)
                        F.precision = 1;
                }
                for (i = 0; n; n /= 8, i++)
                    *--s = (char)(n % 8 + '0');
                if (F.altForm && i && *s != '0') {
                    *--s = '0';
                    i++;
                }
                for (; i < F.precision; i++)
                    *--s = '0';
                break;
                
                /*  hexadecimal (unsigned)  */
                
            case 'p':
                F.havePrecision = F.lSize = TRUE;
                F.precision = 8;
                /* ... */
            case 'X':
                digits = "0123456789ABCDEF";
                goto hexadecimal;
            case 'x':
                digits = "0123456789abcdef";
                /* ... */
            hexadecimal:
                if (F.lSize)
                    n = va_arg(arg, unsigned long);
                else
                    n = va_arg(arg, unsigned int);
                if (F.hSize)
                    n = (unsigned short) n;
                if (!F.havePrecision) {
                    if (F.zeroPad) {
                        F.precision = F.fieldWidth;
                        if (F.altForm)
                            F.precision -= 2;
                    }
                    if (F.precision < 1)
                        F.precision = 1;
                }
                for (i = 0; n; n /= 16, i++)
                    *--s = digits[n % 16];
                for (; i < F.precision; i++)
                    *--s = '0';
                if (F.altForm) {
                    *--s = (char)c;
                    *--s = '0';
                    i += 2;
                }
                break;
                
                /*  fixed-point  */

            case 'f':
                x = va_arg(arg, double);
                if (!F.havePrecision)
                    F.precision = JPMCDS_DEF_PRECISION;
                sprintf(buf2, "%.*f", F.precision, x);
                goto floating;

                /*  floating-point  */
                
            case 'e':
            case 'E':
                x = va_arg(arg, double);
                if (!F.havePrecision)
                    F.precision = 6;
                if (c == 'e')
                    sprintf(buf2, "%.*e", F.precision, x);
                else
                    sprintf(buf2, "%.*E", F.precision, x);
                goto floating;
                
                /*  "general" notation  */
                
            case 'g':
            case 'G':
                x = va_arg(arg, double);
                if (!F.havePrecision)
                    F.precision = 6;
                else if (F.precision == 0)
                    F.precision = 1;
                F.exponent = (char)(c - 2);
                if (c == 'g')
                    sprintf(buf2, "%.*g", F.precision, x);
                else
                    sprintf(buf2, "%.*G", F.precision, x);
                goto floating;
                
                    /*  floating (common code)  */

            floating:
                i = strlen(buf2);
                s -= i;
                strncpy(s, buf2, i);
                break;

                /*  character  */
                
            case 'c':
                *--s = (char)va_arg(arg, int);
                i = 1;
                break;
                
                /*  string  */
                
            case 's':
                s = va_arg(arg, char *);
                if (F.altForm) {
                    i = (unsigned char) *s++;
                    if (F.havePrecision && i > F.precision)
                        i = F.precision;
                }
                else {
                    if (!F.havePrecision)
                        i = strlen(s);
                    else if ((t = memchr(s, '\0', F.precision)) != 0)
                        i = t - s;
                    else
                        i = F.precision;
                }
                break;
                
                /*  store # bytes written so far  */
                
            case 'n':
                s = va_arg(arg, char *);
                if (F.hSize)
                    * (short *) s = (short)nwritten;
                else if (F.lSize)
                    * (long *) s = nwritten;
                else
                    * (int *) s = nwritten;
                continue;
            
                /*  oops - unknown conversion, abort  */
                
            default:
                if (c != '%')
                    goto done;
            copy1:
                if (JpmcdsLocalPutc((char)c, tFile) < 0)
                    return EOF;
                ++nwritten;
                continue;
        }
            
            /*  pad on the left  */
            
        if (i < F.fieldWidth && !F.leftJustify) {
            do {
                if (JpmcdsLocalPutc(' ', tFile) < 0)
                    return EOF;
                ++nwritten;
            } while (i < --F.fieldWidth);
        }
        
            /*  write the converted result  */

        if (i != 0)
        {
            if (JpmcdsLocalFwrite(s, i, tFile) != i)
                return EOF;
        }
        nwritten += i;
            
            /*  pad on the right  */
            
        for (; i < F.fieldWidth; i++) {
            if (JpmcdsLocalPutc(' ', tFile) < 0)
                return EOF;
            ++nwritten;
        }
    }
    
done:
    if (tFile->type == TFILE_STRING)
        JpmcdsLocalPutc('\0', tFile);

    return nwritten;
}
