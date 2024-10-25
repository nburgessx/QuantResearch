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

#ifndef XLFUNCTONS_H
#define XLFUNCTONS_H

#include "XLOper.h"
#include <vector>
#include <string>
#ifndef GENERIC_EXCEPTION_TYPE
#include <stdexcept>
#define GENERIC_EXCEPTION_TYPE std::runtime_error
#endif

class XLFunction {

public:

  std::string DLLName, XLName, DataCodes, ArgumentNames, FunctionDescription;
  std::vector<const char*> ArgumentDescriptions;

  static const std::vector<XLFunction>& DeclaredFunctions();

  size_t DataCodeArgumentCount() const;

  static size_t MaximumArgumentCount();

  //
  // Regarding the number of permissible arguments to Excel worksheet functions, taken from https://learn.microsoft.com/en-us/office/client-developer/excel/xlfregister-form-1 :
  //    In Excel 2003 and earlier, xlfRegister can take, at most, 30 arguments so that you can provide this help for the first 20 of your function arguments only.
  //    Starting in Excel 2007, xlfRegister can take up to 255 arguments so that you can provide this help for up to 245 function parameters.
  //
  // In other words, as of Excel 2007, we can support up to 245 arguments, though, at present, this code only supports up to 20.
  //

  // data_codes: see https://docs.microsoft.com/en-us/office/client-developer/excel/data-types-used-by-excel for the data type codes. The first character is the return type, then one character per argument.
  // argument_descriptions: one per argument, i.e., one fewer than the data code characters (the first one denotes the function's return type).
  static int Declare(const char* dll_name, const char* xl_name /* may be null */, const std::string& data_codes, const std::string& argument_names, const std::string& function_descripton, const std::vector<const char*>& argument_descriptions = {});

  // Legacy interface.
  // i=1; while [ $i -le 20 ] ; do  echo -n ", const char *a$i=0" ; i=$[ $i + 1 ] ; done
  // xl_name: if null, the dll_name is used. If it starts with period '.', then it is prefixed by the XLL name.
  // data_codes: see https://docs.microsoft.com/en-us/office/client-developer/excel/data-types-used-by-excel for the data type codes.
  static int Declare(const char* dll_name, const char* xl_name /* may be null */, const std::string& data_codes, const std::string& argument_names, const std::string& function_descripton,
                     const char* a1, const char* a2 = 0, const char* a3 = 0, const char* a4 = 0, const char* a5 = 0, const char* a6 = 0, const char* a7 = 0, const char* a8 = 0, const char* a9 = 0, const char* a10 = 0,
                     const char* a11 = 0, const char* a12 = 0, const char* a13 = 0, const char* a14 = 0, const char* a15 = 0, const char* a16 = 0, const char* a17 = 0, const char* a18 = 0, const char* a19 = 0, const char* a20 = 0);

  static int Declare(const char* dll_name, const std::vector<const char*>& xl_names /* allows synonyms */, const std::string& data_codes, const std::string& argument_names, const std::string& function_descripton, const std::vector<const char*>& argument_descriptions = {});

  static bool CalledFromFunctionWizard();

  static void DisallowFromFunctionWizard();

};

//
// Regarding parallel execution of Excel worksheet functions, taken from https://learn.microsoft.com/en-us/office/client-developer/excel/xlfregister-form-1 :
//
//    Registering worksheet functions as thread-safe
//    Starting in Excel 2007, Excel can perform multithreaded workbook recalculation. This means that it can assign different instances of a thread-safe function to concurrent threads for reevaluation.
//    Starting in Excel 2007, most of the built-in worksheet functions are thread-safe.Starting in Excel 2007, Excel also allows XLLs to register worksheet functions as thread-safe. To do this,
//    
//                 include a $ character after the last parameter code in pxTypeText.
//                 =================================================================
//
//      Note: Only worksheet functions can be declared as thread-safe. Excel does not consider a macro sheet equivalent function to be thread-safe,
//            so that you cannot append both #and $ characters to the pxTypeText argument.
//
//    If you have registered a function as thread-safe, you must ensure that it behaves in a thread-safe way, although Excel rejects any thread-unsafe calls via the C API.
//    For example, if a thread-safe function tries to call xlfGetCell, the call fails with the xlretNotThreadSafe error.

// Volatile functions and recalculation, also taken from https://learn.microsoft.com/en-us/office/client-developer/excel/xlfregister-form-1
//
//    On a worksheet, you can make a DLL function or code resource volatile, so that it recalculates every time the worksheet recalculates. To do this,
//
//                 add an exclamation mark (!) after the last argument code in the pxTypeText argument.
//                 ====================================================================================
//

// NOTE: any (worksheet) function name specified with a leading '.' will be registered preceded by this XLL's name.

#define DECLARE_XL_FUNCTION( FUNCTION, EXCEL_FUNCTION_NAME /* This can also be a simple list like { "", "." } */, DATA_CODES, ARGUMENT_NAMES, FUNCTION_DESCRIPTION, ... ) \
   [[maybe_unused]] static int xl_function_##FUNCTION = XLFunction::Declare( #FUNCTION, EXCEL_FUNCTION_NAME, DATA_CODES, ARGUMENT_NAMES, FUNCTION_DESCRIPTION, __VA_ARGS__ );

#endif // XLFUNCTONS_H
