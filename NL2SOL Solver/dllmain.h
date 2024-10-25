//
// This source code resides at www.jaeckel.org/LetsBeRational.7z .
//
// ======================================================================================
// Copyright © 2013-2024 Peter Jäckel.
//
// Permission to use, copy, modify, and distribute this software is freely granted,
// provided that this notice is preserved.
//
// WARRANTY DISCLAIMER
// The Software is provided "as is" without warranty of any kind, either express or implied,
// including without limitation any implied warranties of condition, uninterrupted use,
// merchantability, fitness for a particular purpose, or non-infringement.
// ======================================================================================
//
#ifndef   DLLMAIN_H
#define   DLLMAIN_H

#ifndef DLL_EXPORT
# if defined( WIN32 ) || defined( _WIN32 ) || defined( WIN64 ) || defined( _WIN64 )
#  define DLL_EXPORT __declspec(dllexport)
# elif defined(__GNUC__) // See http://gcc.gnu.org/wiki/Visibility
#  define DLL_EXPORT __attribute__((visibility("default")))
# else
#  define DLL_EXPORT
# endif
#endif

#if defined( RETURN_BUILDDATE_AS_XLOPER ) && !defined(NO_XL_API)
class XLOper;
extern "C" DLL_EXPORT XLOper* BuildDate();
#else
extern "C" DLL_EXPORT const char* BuildDate();
#endif

extern "C" DLL_EXPORT const char* CPUName();
extern "C" DLL_EXPORT const char* DLLName();
extern "C" DLL_EXPORT const char* DLLDirectory();
extern "C" DLL_EXPORT const char* CompilerVersion();
extern "C" DLL_EXPORT const char * BuildConfiguration();
extern "C" DLL_EXPORT int ExcelAPIMode();
extern "C" DLL_EXPORT int Revision();
extern "C" DLL_EXPORT int Bitness();

#include <string>

std::string to_string(const wchar_t* wcs, size_t length);
std::wstring to_wstring(const char* s, size_t length);

#ifdef NDEBUG
# define UNUSED_IN_RELEASE [[maybe_unused]]
# define DebugPrintf(FMT, ...)
#else
# define UNUSED_IN_RELEASE
void DebugPrintf(const char* format, ...);
#endif

// We advise against the use of the Microsoft-specific mapping functions (MultiByteToWideChar/WideCharToMultiByte):
//  - Code Page 1252 characters are then not mapped adequately under the XLOper12 API.
//  Note that UTF-8 characters are never mapped correctly under the XLOper4API.
//  The best compromise seems to be to use CP1252 in special strings such as the COPYRIGHT macro in version.h and not to use Microsoft-specific string conversion functions.
// #define USE_WINDOWS_SPECIFIC_STRING_CONVERSIONS_WHERE_APPLICABLE

// This sends DebugPrintf() output to the debugger's 'Output' window under Visual Studio (or to Windows' 'DebugView' monitor application).
#define USE_WINDOWS_SPECIFIC_DEBUG_OUTPUT_WHERE_APPLICABLE

#endif // DLLMAIN_H
