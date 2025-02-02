//
// This source code resides at www.jaeckel.org/LetsBeRational.7z .
//
// ======================================================================================
// Copyright © 2023-2024 Peter Jäckel.
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

#define NOMINMAX // to suppress MSVC's definition of min() and max()

#include "dllmain.h"

#if __has_include("version.h")
# if defined(_MSC_VER)
#  pragma warning(push)
#  pragma warning(disable : 4828) // Warning C4828: The file contains a character starting at offset 0xae that is illegal in the current source character set (codepage 65001).
# endif
# include "version.h"
# if defined(_MSC_VER)
#  pragma warning(pop) // #pragma warning(disable : 26812)
# endif
#endif

#if defined(NO_XL_API)
# define DECLARE_XL_FUNCTION( ... )
#else
# include "XLFunctions.h"
#endif

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <assert.h>
#include <vector>
#include <map>

#include <cstdio> // for vsnprintf

#ifdef COPYRIGHT
# define COPYRIGHT_STRING COPYRIGHT " "
#else
# define COPYRIGHT_STRING ""
#endif

#define  TOSTRING_AUX(x) #x
#define  TOSTRING(x)     TOSTRING_AUX(x)

#ifdef REVISION
# define REVISION_STRING " (revision " TOSTRING(REVISION) ") "
#else
# define REVISION_STRING " "
#endif

#if !defined(NDEBUG)

#include <stdarg.h>

# define VERBOSE_DEBUG_LOG

void DebugPrintf(const char* format, ...) {
# if (defined( WIN32 ) || defined( _WIN32 ) || defined( WIN64 ) || defined( _WIN64 )) && defined(USE_WINDOWS_SPECIFIC_DEBUG_OUTPUT_WHERE_APPLICABLE)
  const size_t buflen = 1024;
  char buffer[buflen]; // NULL termination.
  va_list arg;
  va_start(arg, format);
  vsnprintf(buffer, buflen, format, arg);
  va_end(arg);
  ::OutputDebugStringA(buffer);
# else
  va_list arg;
  va_start(arg, format);
  vfprintf(stderr, format, arg);
  va_end(arg);
  fflush(stderr);
# endif
}

#endif

#if !(defined(WIN32) || defined(_WIN32) || defined(WIN64) || defined(_WIN64)) || !defined(USE_WINDOWS_SPECIFIC_STRING_CONVERSIONS_WHERE_APPLICABLE)

# if defined(_MSC_VER)
#  pragma warning(disable : 4996) // Warning	warning C4996 : 'std::codecvt_utf8_utf16<wchar_t,1114111,(std::codecvt_mode)0>' : warning STL4017 : std::wbuffer_convert, std::wstring_convert, and the <codecvt> header ...
# endif

# include <codecvt>
# include <locale>

namespace {
  // See https://en.wikipedia.org/wiki/Windows-1252. The most common characters like ©äöüÄÖÜß etc. directly map to UTF16 via direct cast.
  // According to the information on Microsoft's and the Unicode Consortium's websites, positions 81, 8D, 8F, 90, and 9D are unused;
  // however, the Windows API MultiByteToWideChar maps these to the corresponding C1 control codes. The "best fit" mapping documents this behavior, too.
  const static std::map<char, wchar_t> cp1252_to_utf16_special_cases {
    /* € */ { 0x80, 0x20AC }, /* ‚ */ { 0x82, 0x201A }, /* ƒ */ { 0x83, 0x192 }, /* „ */ { 0x84, 0x201E }, /* … */ { 0x85, 0x2026 }, /* † */ { 0x86, 0x2020 }, /* ‡ */ { 0x87, 0x2021 }, /* ˆ */ { 0x88, 0x02C6 }, /* ‰ */ { 0x89, 0x2030 }, /* Š */ { 0x8A, 0x160 }, /* ‹ */ { 0x8B, 0x2039 }, /* Œ */ { 0x8C, 0x152 }, /* Ž */ { 0x8E, 0x017D },
    /* ‘ */ { 0x91, 0x2018 }, /* ’ */ { 0x92, 0x2019 }, /* “ */ { 0x93, 0x201C }, /* ” */ { 0x94, 0x201D }, /* • */ { 0x95, 0x2022 }, /* – */ { 0x96, 0x2013 }, /* — */ { 0x97, 0x2014 }, /* ˜ */ { 0x98, 0x02DC }, /* ™ */ { 0x99, 0x2122 }, /* š */ { 0x9A, 0x161 }, /* › */ { 0x9B, 0x203A }, /* œ */ { 0x9C, 0x153 }, /* ž */ { 0x9E, 0x017E }, /* Ÿ */ { 0x9F, 0x178 }
  };
  wchar_t cp1252_to_wchar(char c) {
    auto entry = cp1252_to_utf16_special_cases.find(c); // We cast all characters except those explicitly listed with a special mapping without compromise as a last resort.
    return entry == cp1252_to_utf16_special_cases.end() ? (wchar_t)(unsigned char)c : entry->second;
  }

}

#endif

std::string to_string(const wchar_t* wcs, size_t length) {
  std::string result;
  if (wcs) {
#if ( defined( WIN32 ) || defined( _WIN32 ) || defined( WIN64 ) || defined( _WIN64 ) ) && defined( USE_WINDOWS_SPECIFIC_STRING_CONVERSIONS_WHERE_APPLICABLE )
    int output_length = WideCharToMultiByte(CP_UTF8, 0, wcs, (int)length, 0, 0, 0, 0);
    std::vector<char> buffer(output_length);
    WideCharToMultiByte(CP_UTF8, 0, wcs, (int)length, &buffer[0], output_length, 0, 0);
# if defined(VERBOSE_DEBUG_LOG)
    DebugPrintf("\nConverted wide string '%ls' into string '%s'\nvia WideCharToMultiByte().\n\n", wcs, &buffer[0]);
# endif
    result = std::string(&buffer[0], output_length);
#else
    // See https://en.cppreference.com/w/cpp/locale/wstring_convert/to_bytes
    try {
      result = std::wstring_convert<std::codecvt_utf8_utf16<wchar_t>>().to_bytes(wcs, wcs + length);
    } catch (UNUSED_IN_RELEASE const std::exception& e) {
      DebugPrintf("\nException: %s.\n\n", e.what());
      std::vector<char> buffer(length);
      for (size_t i = 0; i < length; ++i)
        buffer[i] = (char)wcs[i]; // non-compromising downcast as a last resort.
      result = std::string(&buffer[0], length);
    }
# if defined(VERBOSE_DEBUG_LOG)
    DebugPrintf("\nConverted wide string '%ls' into string '%s'.\n\n", wcs, result.c_str());
# endif
#endif	
  }
  return result;
}

std::wstring to_wstring(const char* s, size_t length) {
  std::wstring result;
  if (s) {
#if ( defined( WIN32 ) || defined( _WIN32 ) || defined( WIN64 ) || defined( _WIN64 ) ) && defined( USE_WINDOWS_SPECIFIC_STRING_CONVERSIONS_WHERE_APPLICABLE )
    int output_length = MultiByteToWideChar(CP_UTF8, 0, s, (int)length, 0, 0);
    std::vector<wchar_t> buffer(output_length+1); // For NULL termination.
    MultiByteToWideChar(CP_UTF8, 0, s, (int)length, &buffer[0], output_length);
# if defined(VERBOSE_DEBUG_LOG)
    DebugPrintf("\nConverted string '%s' into wide string '%ls'\nvia WideCharToMultiByte().\n\n", s, &buffer[0]);
# endif
    result = std::wstring(&buffer[0], output_length);
#else
    // See https://en.cppreference.com/w/cpp/locale/wstring_convert/from_bytes
    try {
      result = std::wstring_convert<std::codecvt_utf8_utf16<wchar_t>>().from_bytes(s, s + length);
    } catch (UNUSED_IN_RELEASE const std::exception& e) {
      DebugPrintf("\nException: %s.\n\n", e.what());
      // The most likely cause for a conversion exception above is the presence of CP1252 characters (Windows-specific, sometimes erroneously labelled as 'ANSI').
      // We try a last resort recovery with fall-back to a direct cast.
      std::vector<wchar_t> buffer(length);
      for (size_t i = 0; i < length; ++i) 
        buffer[i] = cp1252_to_wchar(s[i]);
      result = std::wstring(&buffer[0], length);
    }
# if defined(VERBOSE_DEBUG_LOG)
    DebugPrintf("\nConverted string '%s' into wide string '%ls'.\n\n", s, result.c_str());
# endif
#endif	
  }
  return result;
}

#if defined( WIN32 ) || defined( _WIN32 ) || defined( WIN64 ) || defined( _WIN64 )
# define DL_OPEN( LIBNAME ) (HMODULE)::LoadLibraryA(LIBNAME)
# define DL_CLOSE( HANDLE ) (0!=::FreeLibrary((HMODULE)(HANDLE)))
# define DL_SYM( HANDLE , NAME ) (void*)::GetProcAddress((HMODULE)(HANDLE),(NAME))
# define sleep(length) Sleep(length*1000);
# include <processenv.h>
# include <shellapi.h>
# include <intrin.h>
#else
# include <unistd.h> // For sleep().
# include <dlfcn.h>
# ifdef RTLD_GROUP
#  define RTLD_FLAGS ( RTLD_LAZY | RTLD_GROUP )
# else
#  ifdef RTLD_DEEPBIND
#   define RTLD_FLAGS ( RTLD_LAZY | RTLD_DEEPBIND )
#  else
#   define RTLD_FLAGS ( RTLD_LAZY )
#  endif
# endif
# define DL_OPEN( LIBNAME ) (void*)::dlopen((LIBNAME), RTLD_FLAGS)
# define DL_CLOSE( HANDLE ) (0==::dlclose(HANDLE))
# define DL_SYM( HANDLE , NAME ) (void*)::dlsym((HANDLE),(NAME))
# if __has_include("filesystem")
#  include <filesystem>
#  define GET_ABSOLUTE_PATH_NAME std::filesystem::absolute
# else
#  include <limits.h>
   inline std::string GET_ABSOLUTE_PATH_NAME(const char*f){
     char buffer[PATH_MAX + 1] = { 0 };
     return realpath(f, buffer) ? buffer : f;
   }
# endif
#endif

#if defined(_MSC_VER)
# define wcsicmp _wcsicmp
#endif

namespace {

  std::string cpu_model_name() { // See https://learn.microsoft.com/de-de/cpp/intrinsics/cpuid-cpuidex?view=msvc-170
    int data[12+1] = { 0 }; // https://en.wikipedia.org/wiki/CPUID says "The string is specified in Intel/AMD documentation to be null-terminated, however this is not always the case [...] and software should not rely on it."
    for (int i = 0; i < 3; ++i)
#ifdef _WIN32
      __cpuidex((int*)&data[i * 4], 0x80000002 + i, 0);
#else
      asm volatile ("cpuid" : "=a" (data[i * 4 + 0]), "=b" (data[i * 4 + 1]), "=c" (data[i * 4 + 2]), "=d" (data[i * 4 + 3]) : "a" (0x80000002 + i), "c" (0));
#endif
    return (const char*)&data[0];
  }

#if !defined(NO_XL_API)

  static std::string toString(int x) { char buffer[32]; snprintf(buffer, 32, "%d", x); return buffer; }

  typedef int (WINAPI* PFN_EXCEL_CALLBACK)(int xlfn, void* pOperRes, int count, void** ppOpers);
  PFN_EXCEL_CALLBACK     pExcel4v = 0;

  typedef int (PASCAL* MDCALLBACK12PROC) (int xlfn, int coper, void** rgpxloper12, void* xloper12Res);
  MDCALLBACK12PROC pMdCallBack12 = 0;

  int RegisteredNumberOfFunctions = 0;

  bool HaveExcel() {
    return 0 != pExcel4v || 0 != pMdCallBack12;
  }

# if defined(_MSC_VER)
#  pragma warning(push)
#  pragma warning(disable : 26812) // Warning	C26812	The enum type XLOper::APIModeType is unscoped.
# endif
  int CallbackToExcelV(int xlfn, XLOper& operRes, int count, void** opers) {
    int result = -1;
    switch (XLOper::APIMode) {
      case XLOper::APIMode4:
        result = pExcel4v(xlfn, &operRes, count, opers);
        break;
      case XLOper::APIMode12:
        result = pMdCallBack12(xlfn, count, opers, &operRes);
        break;
      default:
        DebugPrintf("Excel callback interface not found.\n");
        exit(-1);
        break;
    }
    if (0 == result && xlfRegister == xlfn && xltypeErr != (int)operRes.type()) ++RegisteredNumberOfFunctions;
# if defined(VERBOSE_DEBUG_LOG)
    if (xlfRegister == xlfn) {
      if (0 == result && xltypeErr != (int)operRes.type())
        DebugPrintf("Registered function with Excel name '%s' and DLL name '%s'.\n", std::string(*(XLOper*)opers[3]).c_str(), ((XLOper*)opers[1])->operator std::string().c_str());
      else
        DebugPrintf("Failed to register function with Excel name '%s' and DLL name '%s'.\n", std::string(*(XLOper*)opers[3]).c_str(), ((XLOper*)opers[1])->operator std::string().c_str());
    }
# endif
    return result;
  }
# if defined(_MSC_VER)
#  pragma warning(pop) // #pragma warning(disable : 26812)
# endif

  // i=1; while [ $i -le 30 ] ; do  echo -n ", const XLOper &o$i=XLOper::Missing" ; i=$[ $i + 1 ] ; done
  int CallbackToExcel(int xlfn, XLOper& operRes, int count, const XLOper& o1 = XLOper::Missing, const XLOper& o2 = XLOper::Missing, const XLOper& o3 = XLOper::Missing, const XLOper& o4 = XLOper::Missing, const XLOper& o5 = XLOper::Missing, const XLOper& o6 = XLOper::Missing, const XLOper& o7 = XLOper::Missing, const XLOper& o8 = XLOper::Missing, const XLOper& o9 = XLOper::Missing, const XLOper& o10 = XLOper::Missing, const XLOper& o11 = XLOper::Missing, const XLOper& o12 = XLOper::Missing, const XLOper& o13 = XLOper::Missing, const XLOper& o14 = XLOper::Missing, const XLOper& o15 = XLOper::Missing, const XLOper& o16 = XLOper::Missing, const XLOper& o17 = XLOper::Missing, const XLOper& o18 = XLOper::Missing, const XLOper& o19 = XLOper::Missing, const XLOper& o20 = XLOper::Missing, const XLOper& o21 = XLOper::Missing, const XLOper& o22 = XLOper::Missing, const XLOper& o23 = XLOper::Missing, const XLOper& o24 = XLOper::Missing, const XLOper& o25 = XLOper::Missing, const XLOper& o26 = XLOper::Missing, const XLOper& o27 = XLOper::Missing, const XLOper& o28 = XLOper::Missing, const XLOper& o29 = XLOper::Missing, const XLOper& o30 = XLOper::Missing) {
    // i=1; while [ $i -le 30 ] ; do  echo -n "&o$i," ; i=$[ $i + 1 ] ; done
    const void* opers[] = { &o1,&o2,&o3,&o4,&o5,&o6,&o7,&o8,&o9,&o10,&o11,&o12,&o13,&o14,&o15,&o16,&o17,&o18,&o19,&o20,&o21,&o22,&o23,&o24,&o25,&o26,&o27,&o28,&o29,&o30 };
    return CallbackToExcelV(xlfn, operRes, count, (void**)opers);
  }

  std::string replace(std::string hay_stack, const char* needle, const std::string& needle_replacement) {
    size_t pos = 0, l = strlen(needle);
    while ((pos = hay_stack.find(needle, pos)) != std::string::npos) {
      hay_stack.replace(pos, l, needle_replacement);
      pos += needle_replacement.length();
    }
    return hay_stack;
  }

#endif // !defined(NO_XL_API)

  std::string get_directory(const std::string& file_name) {
    const size_t b = file_name.find_last_of("/\\", std::string::npos, 2);
    return std::string::npos == b ? std::string() : file_name.substr(0, b);
  }

#if defined( WIN32 ) || defined( _WIN32 ) || defined( WIN64 ) || defined( _WIN64 )
  HMODULE ThisDLLsModuleHandle;
#endif

#if defined(__GNUC__)
# pragma GCC diagnostic push
# pragma GCC diagnostic ignored "-Wdate-time"
#endif
  std::string ThisDLLsTimeStamp = __DATE__ " " __TIME__;
#if defined(__GNUC__)
# pragma GCC diagnostic pop
#endif

  // As returned from a call back to Excel's xlGetName function.
  std::string get_file_name() {
#if !defined(NO_XL_API)
    if (HaveExcel()) {
      XLOper xDLL, xlRet; // NIL initialised.
      CallbackToExcelV(xlGetName, xDLL, 0, 0); // This causes Excel to place memory into xDLL which we don't own - we need to delegate freeing it up back to Excel via the xlFree callback further down.
      std::string file_name(xDLL);
# if defined(VERBOSE_DEBUG_LOG)
      DebugPrintf("xlGetName returned %s.\n", file_name.c_str());
# endif
      const void* opers[] = { &xDLL };
      CallbackToExcelV(xlFree, xlRet, 1, (void**)opers);
      new (&xDLL) XLOper(); // Set to NIL without any cleaning up. This is to suppress the destructor of class XLOper attempting to release memory that is here actually managed by Excel.
      return file_name;
    }
#endif // !defined(NO_XL_API)
#if defined( WIN32 ) || defined( _WIN32 ) || defined( WIN64 ) || defined( _WIN64 )
    char buffer[_MAX_PATH * 2] = { 0 };
    ::GetModuleFileNameA(ThisDLLsModuleHandle, buffer, sizeof(buffer));
    return buffer;
#else
    Dl_info dlInfo;
    dladdr((void*)get_file_name, &dlInfo);
    return GET_ABSOLUTE_PATH_NAME(dlInfo.dli_fname);
#endif
  }

#ifdef XL_CATEGORY_NAME
# define THIS_XL_CATEGORY_NAME TOSTRING(XL_CATEGORY_NAME)
#else
# define THIS_XL_CATEGORY_NAME this_xll_base_name()
#endif

#ifdef XL_PREFIX
  const std::string& this_xl_prefix() { static std::string name = TOSTRING(XL_PREFIX); return name; }
# define THIS_XL_PREFIX this_xl_prefix()
#else
# define THIS_XL_PREFIX this_xl_category_name()
#endif

  // As returned from a call back to Excel's xlGetName function.
  const std::string& this_file_name() { static std::string name = get_file_name(); return name; }

  const char* skip_path(const char* file_name) { return (const char*)std::max((uintptr_t)file_name, std::max((uintptr_t)strrchr(file_name, '/'), (uintptr_t)strrchr(file_name, '\\')) + 1); }

  // The file name of the DLL, without path, but including the extension (i.e., ".xll").
  const std::string& this_xll_name() { static std::string name = skip_path(this_file_name().c_str()); return name; }

#if !defined(NO_XL_API)

# if defined(VERBOSE_DEBUG_LOG) || !defined(XL_CATEGORY_NAME)
  std::string strip_extension(const std::string& file_name) {
    const size_t e = file_name.find_last_of('.');
    return std::string::npos == e ? file_name : file_name.substr(0, e);
  }
  // The file name of the DLL, without path, and without any extension (i.e., ".xll").
  const std::string& this_xll_base_name() { static std::string name = strip_extension(this_xll_name()); return name; }
# endif

  // If XL_CATEGORY_NAME is not defined, the Excel function category is determined dynamically as the XLL's base file name (without path and without extension).
  const std::string& this_xl_category_name() { static std::string name = THIS_XL_CATEGORY_NAME; return name; }

#endif // !defined(NO_XL_API)

}

#if !defined(NO_XL_API)

void RegisterAllFunctions() {
# if defined(VERBOSE_DEBUG_LOG)
  DebugPrintf("RegisterAllFunctions().\n");
# endif
  RegisteredNumberOfFunctions = 0;
  const std::string& file_name = this_file_name(), & xll_name = this_xll_name(), & xl_category_name = this_xl_category_name(), & xl_prefix = THIS_XL_PREFIX;
  XLOper xDLL(file_name), xlRet/*NIL initialised.*/;
# if defined(VERBOSE_DEBUG_LOG)
  DebugPrintf("This dll name: %s .\nThis dll base name: %s.\n", xll_name.c_str(), this_xll_base_name().c_str());
# endif
  std::string message("Registering library '" + xl_category_name + "'" + REVISION_STRING + COPYRIGHT_STRING + file_name + " ...");
  const size_t max_arg_count = XLFunction::MaximumArgumentCount();
  CallbackToExcel(xlcMessage, xlRet, 2, true, message.c_str());
  const std::vector<XLFunction>& funcs = XLFunction::DeclaredFunctions();
  XLOper xl_function_type("1"), xl_category(xl_category_name), xl_shortcut_text(""), xl_help_topic("");
  for (size_t i = 0; i < funcs.size(); ++i) {
    XLOper dll_name(funcs[i].DLLName), data_codes(XLOper::replace_xloper_codes(funcs[i].DataCodes.c_str())), xl_name('.' == funcs[i].XLName[0] ? xl_prefix + funcs[i].XLName : funcs[i].XLName), arg_names(funcs[i].ArgumentNames), func_desc(replace(funcs[i].FunctionDescription, "@THIS_XLL_NAME@", xll_name));
    std::vector< XLOper*> opers{ &xDLL, & dll_name, & data_codes, & xl_name, & arg_names, & xl_function_type, & xl_category, & xl_shortcut_text, & xl_help_topic, & func_desc };
    const size_t n_base = opers.size(), arg_count = std::min(max_arg_count, std::min(funcs[i].DataCodeArgumentCount(), funcs[i].ArgumentDescriptions.size()));
    std::vector<XLOper> arg_descs(arg_count);
    opers.resize(n_base + arg_count);
    for (size_t j = 0; j < arg_count; ++j)
      new ((opers[n_base + j] = &arg_descs[j])) XLOper(funcs[i].ArgumentDescriptions[j]);
    // See https://learn.microsoft.com/en-us/office/client-developer/excel/xlfregister-form-1 .
    CallbackToExcelV(xlfRegister, xlRet, (int)opers.size(), (void**)&opers[0]);
  }
  message = "Registered " + toString(RegisteredNumberOfFunctions) + (1 == RegisteredNumberOfFunctions ? " function" : " functions") + " from library " + file_name + ".";
# if defined(VERBOSE_DEBUG_LOG)
  DebugPrintf("%s\n", message.c_str());
# endif
  // For some versions of Excel, the below message never appears - we only see the above "Registering library ...".
  CallbackToExcel(xlcMessage, xlRet, 2, true, message);
  sleep(2); // To give us a chance to notice the message 'Registered ...' at the bottom of the window.
  CallbackToExcel(xlcMessage, xlRet, 1, false);
}

#endif // !defined(NO_XL_API)

#if defined( WIN32 ) || defined( _WIN32 ) || defined( WIN64 ) || defined( _WIN64 )
std::string dll_time_stamp(void* dll_handle) {
  IMAGE_DOS_HEADER* dosHeader = (IMAGE_DOS_HEADER*)dll_handle;
  assert(dosHeader->e_magic == IMAGE_DOS_SIGNATURE);
  const BYTE* base_address = (BYTE*)(dosHeader);
  IMAGE_NT_HEADERS* ntHeaders = (IMAGE_NT_HEADERS*)(base_address + dosHeader->e_lfanew);
  assert(ntHeaders->Signature == IMAGE_NT_SIGNATURE);
  const time_t t = ntHeaders->FileHeader.TimeDateStamp;
  char str[26];
  ctime_s(str, sizeof str, &t);
  str[24] = 0; // Wipe the \n character.
  return str;
}
#endif

#if defined( RETURN_BUILDDATE_AS_XLOPER ) && !defined(NO_XL_API)
XLOper * BuildDate() { return toExcel(ThisDLLsTimeStamp); }
// The code '?' becomes 'P' in API version 4, and 'Q' in API version 12. The suffix '!$' means volatile and thread-safe.
// "." as the function name means it will be registered to appear in Excel with the same name as here in the code, prefixed with "<this_dll_name>.".
DECLARE_XL_FUNCTION(BuildDate, ".", "?!$", "", "returns the build date of this DLL (@THIS_XLL_NAME@).", {});
#else
const char* BuildDate() { return ThisDLLsTimeStamp.c_str(); }
// "." as the function name means it will be registered to appear in Excel with the same name as here in the code, prefixed with "<this_dll_name>.".
DECLARE_XL_FUNCTION(BuildDate, ".", "C!$", "", "returns the build date of this DLL (@THIS_XLL_NAME@).", {});
#endif

const char* CPUName() { static std::string cpu = cpu_model_name();  return cpu.c_str(); }
DECLARE_XL_FUNCTION(CPUName, { "", "." }, "C!$", "", "returns the CPU model name of the executing hardware.", {});

const char* DLLName() { return this_xll_name().c_str(); }
DECLARE_XL_FUNCTION(DLLName, ".", "C!$", "", "returns the name of this DLL (@THIS_XLL_NAME@).", {});

const char* DLLDirectory() { static std::string directory = get_directory(this_file_name());  return directory.c_str(); }
DECLARE_XL_FUNCTION(DLLDirectory, ".", "C!$", "", "returns the path to this DLL (@THIS_XLL_NAME@).", {});

const char* CompilerVersion() {
  const char* s =
#if   defined(_MSC_VER)
    "MSVC " TOSTRING(_MSC_VER);
#elif defined(__MINGW64__)
    "MinGW-w64 " TOSTRING(__MINGW64_VERSION_MAJOR) "." TOSTRING(__MINGW64_VERSION_MINOR) " (GCC " TOSTRING(__GNUC__) "." TOSTRING(__GNUC_MINOR__) "." TOSTRING(__GNUC_PATCHLEVEL__) ")";
#elif defined(__MINGW32__)
    "MinGW 32 " TOSTRING(__MINGW32_MAJOR_VERSION) "." TOSTRING(__MINGW32_MINOR_VERSION) " (GCC " TOSTRING(__GNUC__) "." TOSTRING(__GNUC_MINOR__) "." TOSTRING(__GNUC_PATCHLEVEL__) ")";
#elif defined(__GNUC__)
    "GCC " TOSTRING(__GNUC__) "." TOSTRING(__GNUC_MINOR__) "." TOSTRING(__GNUC_PATCHLEVEL__);
#else
    "<unknown>";
#endif
  return s;
}

// "." means it will be registered to appear in Excel with the same name as here in the code, prefixed with "<this_dll_name>.".
// The code '?' becomes 'P' in API version 4, and 'Q' in API version 12. The suffix '!$' means volatile and thread-safe.
DECLARE_XL_FUNCTION(CompilerVersion, ".", "C!$", "", "returns the compiler name and version number that this DLL (@THIS_XLL_NAME@) was built with.", {});

const char* BuildConfiguration() {
#ifdef NDEBUG
  return "RELEASE";
#else
  return "DEBUG";
#endif
}
DECLARE_XL_FUNCTION(BuildConfiguration, ".", "C!$", "", "returns the build configuration, i.e., DEBUG or RELEASE, of this DLL (@THIS_XLL_NAME@).", {})

#if !defined(NO_XL_API)
int ExcelAPIMode() { return (int)XLOper::APIMode; }
// "." means it will be registered to appear in Excel with the same name as here in the code, prefixed with "<this_dll_name>.".
// Data type 'J' is a 32 bit integer, see https://docs.microsoft.com/en-us/office/client-developer/excel/data-types-used-by-excel.
DECLARE_XL_FUNCTION(ExcelAPIMode, ".", "J!$", "", "returns the Excel API version in use.", {})
#endif

int Bitness() { return 8 * (int)sizeof(intptr_t); }

int Revision() { return REVISION; }
DECLARE_XL_FUNCTION(Revision, ".", "J!$", "", "returns the revision number of this DLL.", {})

// "." means it will be registered to appear in Excel with the same name as here in the code, prefixed with "<this_dll_name>.".
// Data type 'J' is a 32 bit integer, see https://docs.microsoft.com/en-us/office/client-developer/excel/data-types-used-by-excel.
DECLARE_XL_FUNCTION(Bitness, ".", "J!$", "", "returns the bitness of this DLL.", {})

#if defined( WIN32 ) || defined( _WIN32 ) || defined( WIN64 ) || defined( _WIN64 )

#if !defined(NO_XL_API)

extern "C" DLL_EXPORT double ExcelVersion() { return XLOper::ExcelVersion; }
DECLARE_XL_FUNCTION(ExcelVersion, { "", "." }, "B!$", "", "returns the Excel version number.", {});

namespace {
  const std::string& this_addin_info() {
    static std::string addin_info = this_xl_category_name() + REVISION_STRING + COPYRIGHT_STRING + "Excel add-in.";
    return addin_info;
  }
}

extern "C" DLL_EXPORT XLOper * AddinInfo() {
# if defined(VERBOSE_DEBUG_LOG)
  DebugPrintf("AddinInfo() returns \"%s\".\n", this_addin_info().c_str());
# endif
  return toExcel(this_addin_info());
}
// The code '?' becomes 'P' in API version 4, and 'Q' in API version 12. The suffix '!$' means volatile and thread-safe.
DECLARE_XL_FUNCTION(AddinInfo, ".", "?!$", "", "returns add-in information for this DLL (@THIS_XLL_NAME@).", {});

#endif // !defined(NO_XL_API)

extern "C" BOOL APIENTRY DllMain(HMODULE hModule, DWORD reason, LPVOID lpReserved) {
  (void)lpReserved;
  switch (reason) {
    case DLL_PROCESS_ATTACH:
    {
# if defined(VERBOSE_DEBUG_LOG)
      DebugPrintf("Process attaching.\n",0);
# endif
      ThisDLLsModuleHandle = hModule;
      ThisDLLsTimeStamp = dll_time_stamp(hModule);
# if defined(VERBOSE_DEBUG_LOG)
      DebugPrintf("DLL build date: %s.\n", ThisDLLsTimeStamp.c_str());
# endif
# if !defined(NO_XL_API)
      HMODULE module_handle = ::GetModuleHandle(0);
      if (0 != module_handle) {
        pMdCallBack12 = (MDCALLBACK12PROC)DL_SYM(module_handle, "MdCallBack12");
        if (0 != pMdCallBack12)
          XLOper::APIMode = XLOper::APIMode12;
      }
      ::SetErrorMode(SEM_FAILCRITICALERRORS); // avoid message box if library not found
      HMODULE xlcall32_handle = DL_OPEN("xlcall32.dll");
      ::SetErrorMode(0);
      if (0 != xlcall32_handle) {
        pExcel4v = (PFN_EXCEL_CALLBACK)DL_SYM(xlcall32_handle, "Excel4v");
        if (0 != pExcel4v) {
          int argc = 0;
          wchar_t** argv = CommandLineToArgvW(GetCommandLineW(), &argc);
          if (0 != argv) {
            for (int i = 1; i < argc; ++i) {
              const wchar_t* excel4_api_mode_flag = L"--APIMode=Excel4";
              if (0 == wcsicmp(argv[i], excel4_api_mode_flag) || 0 == wcsicmp(argv[i], excel4_api_mode_flag + 1) || 0 == wcsicmp(argv[i], excel4_api_mode_flag + 2)) {
#  ifdef _DEBUG
                DebugPrintf("Found command line flag '%ls'.\n", argv[i]);
                if (0 != pMdCallBack12)
                  DebugPrintf("Disabling API mode 'Excel12'\nFalling back to API mode 'Excel4'.\n");
                else
                  DebugPrintf("API mode 'Excel12' was not available.\nStaying in API mode 'Excel4'.\n");
#  endif
                pMdCallBack12 = 0;
                break;
              }
            }
            LocalFree(argv);
          }
          if (0 == pMdCallBack12)
            XLOper::APIMode = XLOper::APIMode4;
        }
      }
#  if defined(VERBOSE_DEBUG_LOG)
      if (0 != pExcel4v)      DebugPrintf("Found 'Excel4v'.\n");
      if (0 != pMdCallBack12) DebugPrintf("Found 'MdCallBack12'.\n");
      if (!HaveExcel()) DebugPrintf("Found no Excel callback interface.\n");
      DebugPrintf("Selected API version: %d.\n", (int)XLOper::APIMode);
#  endif
# endif // !defined(NO_XL_API)
      break;
    }
# if defined(VERBOSE_DEBUG_LOG)
    case DLL_THREAD_ATTACH:
      DebugPrintf("Thread attaching.\n");
      break;
    case DLL_THREAD_DETACH:
      DebugPrintf("Thread detaching.\n");
      break;
    case DLL_PROCESS_DETACH:
      DebugPrintf("Process detaching.\n");
      break;
# endif
  }
  return TRUE;
}

# if !defined(NO_XL_API)

// In order to avoid signatures such as '_xlAutoOpen@0', etc., to show up in the export table of the DLL in this configuration,
// we do not export the WINAPI functions via a compiler directive, but instead use a linker directive (see end of this file).
// Note that GCC, by default, also decorates exported WINAPI functions, e.g., 'xlAutoOpen@0'.

extern "C" int WINAPI xlAutoOpen() {
  // Set Excel version number for later use.
  XLOper xlVersionInfo, xlVersionNumber, xlTwo(2), xlRet, xlTypeNum(xltypeNum);
  const void* opers[2] = { &xlTwo, 0 };
  CallbackToExcelV(xlfGetWorkspace, xlVersionInfo, 1, (void**)opers);
  opers[0] = &xlVersionInfo;
  opers[1] = &xlTypeNum;
  CallbackToExcelV(xlCoerce, xlVersionNumber, 2, (void**)opers);
  CallbackToExcelV(xlFree, xlRet, 1, (void**)opers);
  new (&xlVersionInfo) XLOper(); // Set to NIL without any cleaning up. This is to suppress the destructor of class XLOper attempting to release memory that is here actually managed by Excel.
  XLOper::ExcelVersion = xlVersionNumber;

  // Set Excel window handle for later use.
  CallbackToExcelV(xlGetHwnd, xlRet, 0, 0);
  XLOper::ExcelWindowHandle = xlRet.to_int();

  // Store process id for later use.
  XLOper::ProcessID = GetCurrentProcessId();

  RegisterAllFunctions();
  return 1;
}

extern "C" int WINAPI xlAutoClose() {
  return 1;
}

extern "C" void WINAPI xlAutoFree(XLOper4 * p) {
#  if defined(VERBOSE_DEBUG_LOG)
  DebugPrintf("xlAutoFree(ptr=0x%p, type=0x%04X)\n", p, p ? p->xltype : (WORD)-1);
#  endif
  if (p && (p->xltype & xlbitDLLFree)) delete p;
}

extern "C" void WINAPI xlAutoFree12(XLOper12 * p) {
#  if defined(VERBOSE_DEBUG_LOG)
  DebugPrintf("xlAutoFree12(ptr=0x%p, type=0x%04X)\n", p, p ? p->xltype : -1);
#  endif
  if (p && (p->xltype & xlbitDLLFree)) delete p;
}

extern "C" XLOper4 * WINAPI xlAddInManagerInfo(XLOper4 * pxAction) {
#  if defined(VERBOSE_DEBUG_LOG)
  DebugPrintf("xlAddInManagerInfo(ptr=0x%p, type=0x%04X)\n", pxAction, pxAction ? pxAction->xltype : (WORD)-1);
#  endif
  // See https://learn.microsoft.com/en-us/office/client-developer/excel/xladdinmanagerinfo-xladdinmanagerinfo12.
  XLOper4 xIntAction, xDestType(xltypeInt);
  const void* opers[2] = { pxAction, &xDestType };
  pExcel4v(xlCoerce, &xIntAction, 2, (void**)opers);
  if (1 == xIntAction.val.w)
    return (new XLOper4(this_addin_info()))->toExcel();
  return 0;
}

extern "C" XLOper12 * WINAPI xlAddInManagerInfo12(XLOper12 * pxAction) {
#  if defined(VERBOSE_DEBUG_LOG)
  DebugPrintf("xlAddInManagerInfo12(ptr=0x%p, type=0x%04X)\n", pxAction, pxAction ? pxAction->xltype : -1);
#  endif
  // See https://learn.microsoft.com/en-us/office/client-developer/excel/xladdinmanagerinfo-xladdinmanagerinfo12.
  XLOper12 xIntAction, xDestType(xltypeInt);
  const void* opers[2] = { pxAction, &xDestType };
  pMdCallBack12(xlCoerce, 2, (void**)opers, &xIntAction);
  if (1 == xIntAction.val.w)
    return (new XLOper12(this_addin_info()))->toExcel();
  return 0;
}

#  if defined(_MSC_VER)
#   if defined( WIN64 ) || defined( _WIN64 )
#    pragma comment(linker, "/EXPORT:xlAutoOpen")
#    pragma comment(linker, "/EXPORT:xlAutoClose")
#    pragma comment(linker, "/EXPORT:xlAutoFree")
#    pragma comment(linker, "/EXPORT:xlAutoFree12")
#    pragma comment(linker, "/EXPORT:xlAddInManagerInfo")
#    pragma comment(linker, "/EXPORT:xlAddInManagerInfo12")
#   else
#    pragma comment(linker, "/EXPORT:xlAutoOpen=_xlAutoOpen@0")
#    pragma comment(linker, "/EXPORT:xlAutoClose=_xlAutoClose@0")
#    pragma comment(linker, "/EXPORT:xlAutoFree=_xlAutoFree@4")
#    pragma comment(linker, "/EXPORT:xlAutoFree12=_xlAutoFree12@4")
#    pragma comment(linker, "/EXPORT:xlAddInManagerInfo=_xlAddInManagerInfo@4")
#    pragma comment(linker, "/EXPORT:xlAddInManagerInfo12=_xlAddInManagerInfo12@4")
#   endif
#  else
asm(".section .drectve");
#   if defined( WIN64 ) || defined( _WIN64 )
asm(".ascii \"-export:xlAutoOpen,xlAutoClose,xlAutoFree,xlAutoFree12,xlAddInManagerInfo,xlAddInManagerInfo12\"");
#   else
asm(".ascii \"-export:xlAutoOpen=xlAutoOpen@0,xlAutoClose=xlAutoClose@0,xlAutoFree=xlAutoFree@4,xlAutoFree12=xlAutoFree12@4,xlAddInManagerInfo=xlAddInManagerInfo@4,xlAddInManagerInfo12=xlAddInManagerInfo12@4\"");
#   endif
#  endif

# endif // !defined(NO_XL_API)

#endif // defined( WIN32 ) || defined( _WIN32 ) || defined( WIN64 ) || defined( _WIN64 )
