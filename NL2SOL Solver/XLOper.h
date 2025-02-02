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

#ifndef XLOPER_H
#define XLOPER_H

#if defined( WIN32 ) || defined( _WIN32 ) || defined( WIN64 ) || defined( _WIN64 )
#  include <windows.h>
#else
#include <cstdint>
// https://learn.microsoft.com/en-us/windows/win32/winprog/windows-data-types
typedef unsigned char BYTE;
typedef uint16_t WORD;
typedef int32_t INT32;
typedef uint32_t DWORD;
typedef int BOOL;
typedef char CHAR;
typedef CHAR* LPSTR;
typedef wchar_t WCHAR;
typedef void* PVOID;
typedef void* LPVOID;
typedef PVOID HANDLE;
typedef HANDLE HWND;
typedef HANDLE HINSTANCE;
typedef HINSTANCE HMODULE;
typedef int32_t LONG;
typedef struct tagPOINT {
  LONG x;
  LONG y;
} POINT;
#define _cdecl
#define __stdcall
#define pascal
#define PASCAL
#define WINAPI __stdcall
#define APIENTRY WINAPI
#define CALLBACK __stdcall
#define VOID void
#endif

#include "XLCall.h"
#include <string>
#include <vector>

//
// Undocumented feature:
//
//   When Excel returns an xlOper of type 0x1000, which is strictly speaking illegal,
//   the returned xlOper contains a wide char string in Pascal encoding. Wide char, in
//   this context means that every character is two bytes long. ASCII characters are
//   encoded with the ASCII code first, followed by a zero (I am not sure if this depends
//   on the Endian convention of the hardware, though). So, the first two bytes, cast into
//   an unsigned short, tell you how many (two byte) characters follow. Excel surmounts
//   the 255 character limitation in this way.
//
//   Note that (xltypeWideCharStr&xltypeType) == 0x0000 !
//
#define xltypeWideCharStr        0x1000

// All bits that are mutually exclusive
#define xltypeType       0x0FFF

// Note that the XLType list is restricted to what we use and support.
enum class XLType {
  Num = xltypeNum,
  Str = xltypeStr,
  Bool = xltypeBool,
  Err = xltypeErr,
  Multi = xltypeMulti,
  Missing = xltypeMissing,
  Nil = xltypeNil,
  Int = xltypeInt
};

// Regarding XLType::Err, the integer member err.val determines the meaning, not all of which are actually 'errors':
enum class XLErr {
  Null = xlerrNull,
  Div0 = xlerrDiv0,
  Value = xlerrValue,
  Ref = xlerrRef,
  Name = xlerrName,
  Num = xlerrNum,          // "#NUM" was always meant to signal Not-A-Number, i.e., NaN
  NaN = xlerrNum,          // "#NUM" was always meant to signal Not-A-Number, i.e., NaN
  NA = xlerrNA,            //  <--  this means Not Applicable, which is not an error at all
  NotApplicable = xlerrNA  //  <--  this means Not Applicable, which is not an error at all
};

//
// A class that is bitwise identical to XLOPER.
//
class XLOper4 : public xloper {

public:

  XLOper4(const char*);
  XLOper4(const wchar_t* s);
  XLOper4(const std::string&);
  XLOper4(const std::wstring&);
  XLOper4(const std::vector<double>&);
  XLOper4(const std::vector<double>&, WORD rows, WORD columns);
  XLOper4(WORD rows, WORD columns, const XLOper4& value);
  XLOper4(WORD rows, WORD columns);
  XLOper4(double);
  XLOper4(bool);
  XLOper4(int);
  XLOper4(XLType);
  XLOper4(); // Nil
  XLOper4(const XLOper4& other);
  ~XLOper4();
  XLOper4* toExcel(); // Sets the xlbitDLLFree free in the xltype field and returns this.
  template <class T> static XLOper4* toExcel(const T& t) { return (new XLOper4(t))->toExcel(); }
  size_t size() const;
  size_t rows() const;
  size_t columns() const;
  XLOper4& operator[] (size_t k);
  const XLOper4& operator[] (size_t k) const;
  XLOper4& operator() (size_t k, size_t l);
  const XLOper4& operator() (size_t k, size_t l) const;
  XLOper4& operator=(const char*);
  XLOper4& operator=(const wchar_t* s);
  XLOper4& operator=(const std::string&);
  XLOper4& operator=(const std::wstring&);
  XLOper4& operator=(double);
  XLOper4& operator=(int);
};


//
// A class that is bitwise identical to XLOPER12.
//
class XLOper12 : public xloper12 {

public:

  XLOper12(const char*);
  XLOper12(const wchar_t* s);
  XLOper12(const std::string&);
  XLOper12(const std::wstring&);
  XLOper12(const std::vector<double>&);
  XLOper12(const std::vector<double>&, RW rows, COL columns);
  XLOper12(RW rows, COL columns, const XLOper12& value);
  XLOper12(RW rows, COL columns);
  XLOper12(double);
  XLOper12(bool);
  XLOper12(int);
  XLOper12(XLType);
  XLOper12(); // Nil
  XLOper12(std::nullptr_t); // Nil
  XLOper12(const XLOper12& other);
  ~XLOper12();
  XLOper12* toExcel(); // Sets the xlbitDLLFree free in the xltype field and returns this.
  XLOper12* Externalise(); // Same as to Excel().
  template <class T> static XLOper12* toExcel(const T& t) { return (new XLOper12(t))->toExcel(); }
  size_t size() const;
  size_t rows() const;
  size_t columns() const;
  XLOper12& operator[] (size_t k);
  const XLOper12& operator[] (size_t k) const;
  XLOper12& operator() (size_t k, size_t l);
  const XLOper12& operator() (size_t k, size_t l) const;
  XLOper12& operator=(const char*);
  XLOper12& operator=(const wchar_t* s);
  XLOper12& operator=(const std::string&);
  XLOper12& operator=(const std::wstring&);
  XLOper12& operator=(double);
  XLOper12& operator=(int);
};

#if defined(_MSC_VER)
# pragma warning(push)
# pragma warning(disable : 26812) // Warning	C26812	The enum type XLOper::APIModeType is unscoped. --- that's blatant nonsense here.
#endif
//
// A polymorphic class that either is an XLOper4 or XLOper12 depending on the value of XLOper::APIMode.
//
class XLOper {

public:

  enum APIModeType {
    APIMode4 = 4,
    APIMode12 = 12
  };

  static APIModeType APIMode;

  static double ExcelVersion;

  // Under APIMode4, this is only the low part of the window handle, but that's enough for the intended purpose (to identify if a function is being called from the function wizard).
  static int ExcelWindowHandle;

  static unsigned ProcessID;

  union {
    struct xloper    o4;
    struct xloper12 o12;
  };

  XLOper(const char*);
  XLOper(const wchar_t* s);
  XLOper(const std::string&);
  XLOper(const std::wstring&);
  XLOper(const std::vector<double>&);
  XLOper(const std::vector<double>&, size_t rows, size_t columns);
  XLOper(size_t rows, size_t columns, const XLOper& value);
  XLOper(size_t rows, size_t columns);
  XLOper(double);
  XLOper(bool);
  XLOper(int);
  XLOper(XLType);
  XLOper(); // Nil
  XLOper(const XLOper& other);
  ~XLOper();

  XLType type() const;
  XLErr err() const;

  const char* type_name() const;

  // Sets the xlbitDLLFree free in the xltype field and returns this.
  XLOper* toExcel();
  template <class T> static XLOper* toExcel(const T& t) { return (new XLOper(t))->toExcel(); }
  size_t size() const;
  size_t rows() const;
  size_t columns() const;
  XLOper& operator[] (size_t k);
  const XLOper& operator[] (size_t k) const;
  XLOper& operator() (size_t k, size_t l);
  const XLOper& operator() (size_t k, size_t l) const;

  XLOper& operator=(const char*);
  XLOper& operator=(const wchar_t* s);
  XLOper& operator=(const std::string&);
  XLOper& operator=(const std::wstring&);
  XLOper& operator=(double);
  XLOper& operator=(int);

  bool isUndefined() const; // Missing or Nil or Err && val.err == NotApplicable

  operator std::string() const;
  operator std::wstring() const;
  int to_int() const;
  bool to_bool() const;
  operator double() const;
  std::vector<double> to_vector() const;
  operator std::vector<double>() const;

  template <class T> T to() const;

  template <class T> T operator | (const T& t) const { return isUndefined() ? t : to<T>(); }

  // 'P' in API version 4, and 'Q' in API version 12.
  static char xloper_code();

  static std::string replace_xloper_codes(const char* codes); // '?' -> 'P' in API version 4, and 'Q' in API version 12.

  // This one is special. It appears as a missing xloper type for both Excel API versions.
  static const XLOper Missing;
  static XLOper NotApplicable();

};
#if defined(_MSC_VER)
# pragma warning(pop) // #pragma warning(disable : 26812)
#endif

template <class T> XLOper* toExcel(const T& t) { return XLOper::toExcel(t); }

#define CATCH_TO_EXCEL catch (const std::exception& e) { return toExcel(e.what()); } catch (const std::string& s) { return toExcel(s); } catch (const char* s) { return toExcel(s); } catch (...) { return toExcel("Unknown exception"); }

#endif // XLOPER_H
