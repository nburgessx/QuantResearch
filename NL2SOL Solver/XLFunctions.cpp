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

#include "XLFunctions.h"

#include <string.h>
#include <stdlib.h>
#include <limits.h>

static std::vector<XLFunction>& declared_functions() {
  static std::vector<XLFunction> functions;
  return functions;
}

const std::vector<XLFunction>& XLFunction::DeclaredFunctions() { return declared_functions(); }

// See https://learn.microsoft.com/en-us/office/client-developer/excel/data-types-used-by-excel
static const char* xl_data_type_codes = "ALBECFDGHIJMNOKPQRU?" /* The '?' is our own extension */;

size_t XLFunction::DataCodeArgumentCount() const {
  size_t n = 0;
  for (size_t i = 0; i < DataCodes.size(); ++i)
    if (0 != strchr(xl_data_type_codes, std::toupper(DataCodes[i])))
      ++n;
  return std::max(n, (size_t)1) - (size_t)1;
}

size_t XLFunction::MaximumArgumentCount() {
  switch (XLOper::APIMode) {
    case XLOper::APIMode4:  return 20;
    case XLOper::APIMode12: return 245;
    default: break;
  }
  return 0;
}

int XLFunction::Declare(const char* dll_name, const char* xl_name /* may be null */, const std::string& data_codes, const std::string& argument_names, const std::string& function_description, const std::vector<const char*>& argument_descriptions) {
  if (0 == xl_name || 0 == xl_name[0])
    xl_name = dll_name;
  std::vector<XLFunction>& functions = declared_functions();
  functions.push_back(XLFunction());
  XLFunction& f = functions[functions.size() - 1];
  f.DLLName = dll_name;
  f.XLName = ('.' == xl_name[0] && 0 == xl_name[1]) ? '.' + f.DLLName : xl_name;
  f.DataCodes = data_codes;
  f.ArgumentNames = argument_names;
  f.FunctionDescription = function_description;
  for (size_t i = 0; i < argument_descriptions.size() && 0 != argument_descriptions[i]; ++i)
    f.ArgumentDescriptions.push_back(argument_descriptions[i]);
  return (int)functions.size();
}

int XLFunction::Declare(const char* dll_name, const char* xl_name /* may be null */, const std::string& data_codes, const std::string& argument_names, const std::string& function_descripton,
                        const char* a1, const char* a2, const char* a3, const char* a4, const char* a5, const char* a6, const char* a7, const char* a8, const char* a9, const char* a10,
                        const char* a11, const char* a12, const char* a13, const char* a14, const char* a15, const char* a16, const char* a17, const char* a18, const char* a19, const char* a20) {
  return Declare(dll_name, xl_name, data_codes, argument_names, function_descripton, { a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15, a16, a17, a18, a19, a20 });
}

int XLFunction::Declare(const char* dll_name, const std::vector<const char*>& xl_names /* need not be null terminated */, const std::string& data_codes, const std::string& argument_names, const std::string& function_description, const std::vector<const char*>& argument_descriptions) {
  for (size_t i = 0; i < xl_names.size(); ++i)
    if (0 != xl_names[i])
      Declare(dll_name, xl_names[i], data_codes, argument_names, function_description, argument_descriptions);
  return 0;
}

#if defined( WIN32 ) || defined( _WIN32 ) || defined( WIN64 ) || defined( _WIN64 )

// See https://learn.microsoft.com/en-us/office/client-developer/excel/how-to-call-xll-functions-from-the-function-wizard-or-replace-dialog-boxes for an introduction.

typedef struct { int call_count; bool found; } xldlg_enum_data;

static BOOL CALLBACK xldlg_enum_thread_windows_proc(HWND hwnd, xldlg_enum_data* p) {
  ++p->call_count;
  const size_t BUFFSIZE = 63;
  char buffer[BUFFSIZE + 1] = { 0 };
  buffer[BUFFSIZE] = 0;
  GetClassName(hwnd, buffer, BUFFSIZE);
  if (_strnicmp(buffer, "bosa_sdm_xl", 11) == 0) {
    buffer[BUFFSIZE] = 0;
    GetWindowText(hwnd, buffer, BUFFSIZE);
    // We only want to detect function wizard windows, not "find and replace" dialogues.
    p->found = 0 == strstr(buffer, "Replace");
    return FALSE; // Tell Windows to stop iterating.
  }
  return TRUE; // Tell Windows to continue iterating.
}

#ifdef USE_FULL_ENUMWINDOWS_TO_IDENTIFY_FUNCTION_WIZARD

static unsigned GetWindowProcessID(HWND hWnd) {
  unsigned pid = 0;
  GetWindowThreadProcessId(hWnd, (LPDWORD)(&pid));
  return pid;
}

static BOOL CALLBACK xldlg_enum_proc(HWND hwnd, xldlg_enum_data* p) {
  if (XLOper::ExcelVersion < 15 ? (short)(intptr_t)GetParent(hwnd) == (short)XLOper::ExcelWindowHandle : GetWindowProcessID(hwnd) == XLOper::ProcessID)
    return xldlg_enum_thread_windows_proc(hwnd, p);
  ++p->call_count;
  return TRUE; // Tell Windows to continue iterating.
}

#endif

bool XLFunction::CalledFromFunctionWizard() {
  xldlg_enum_data data{ 0, false };
#ifdef USE_FULL_ENUMWINDOWS_TO_IDENTIFY_FUNCTION_WIZARD
  EnumWindows((WNDENUMPROC)xldlg_enum_proc, (LPARAM)&data);
#else
  // Note that, unlike as described at https://learn.microsoft.com/en-us/office/client-developer/excel/how-to-call-xll-functions-from-the-function-wizard-or-replace-dialog-boxes,
  // we use EnumThreadWindows() which, when invoked here from a regular worksheet function call actually does zero callbacks to xldlg_enum_thread_windows_proc().
  // This compares with, in my tests on a regular windows session, over 400 callbacks when using EnumWindows(). When invoked from under the function wizard,
  // I see ~30 callbacks, still fewer until we have a match than when using EnumWindows() [~40-50], though, when invoked from under the function wizard,
  // the number of calls is actually not material - we just have to make sure that we catch the fact that we are under the function wizard.
  //
  // Regarding the extra call to GetCurrentThreadId(), see, e.g., https://stackoverflow.com/questions/15007771/is-operation-of-getting-id-of-current-thread-time-expensive:
  // 
  //    GetCurrentThreadId() is quite cheap [effectively zero cost, ed.], and does not involve a system call. Basically it works out to:
  //
  //      int GetCurrentThreadId() {
  //        _asm {
  //          mov eax, fs: [18h]
  //          mov eax, [eax + 24h]
  //        }
  //      }
  EnumThreadWindows(GetCurrentThreadId(), (WNDENUMPROC)xldlg_enum_thread_windows_proc, (LPARAM)&data);
#endif
  return data.found;
}

#else

bool XLFunction::CalledFromFunctionWizard() { return false; }

#endif

void XLFunction::DisallowFromFunctionWizard() {
  if (CalledFromFunctionWizard())
    throw GENERIC_EXCEPTION_TYPE("Calculation suppressed in function wizard");
}
