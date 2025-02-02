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

#include "XLOper.h"
#include "dllmain.h"

#include <string.h>
#include <stdlib.h>
#include <limits.h>
#include <vector>
#include <assert.h>
#include <cmath>
#include <limits>

#include <cstdio> // for snprintf
#include <cwchar> // for swprintf

#if !defined(NDEBUG)
#define VERBOSE_DEBUG_LOG
#endif

#if defined(_MSC_VER)
# pragma warning(push)
# pragma warning(disable : 26812) // Warning	C26812	The enum type XLOper::APIModeType is unscoped. --- not true.
#endif
XLOper::APIModeType XLOper::APIMode = XLOper::APIMode12;
#if defined(_MSC_VER)
# pragma warning(pop) // #pragma warning(disable : 26812)
#endif

double XLOper::ExcelVersion = 0;

int XLOper::ExcelWindowHandle = 0;

unsigned XLOper::ProcessID = 0;

#if defined(VERBOSE_DEBUG_LOG)
void DebugPrint(const char* s) { DebugPrintf("%s", s); }
void DebugPrint(const wchar_t* s) { DebugPrintf("%ls", s); }
const char* string_type(const char*) { return "regular"; }
const char* string_type(const wchar_t*) { return "wide"; }
#endif

namespace {
  size_t length(const char* s) { return s ? strlen(s) : 0; }
  size_t length(const wchar_t* s) { return s ? wcslen(s) : 0; }
  template <typename S> struct pascal_string_intermediate {
    const S* m_s;
    const size_t m_l;
    template <typename T> static size_t max_val(size_t l) { return std::min(l, (size_t)(std::make_unsigned_t<T>) - 1); }
    template <typename T> operator T* () const {
      if (m_s) {
        const size_t l = max_val<T>(m_l);
        T* t = new T[l + 2]{ (T)(std::make_unsigned_t<T>)l, 0 }; // This ensures null-termination. Not required but often helpful.
        std::copy(m_s, m_s + l, t + 1);
#if defined(VERBOSE_DEBUG_LOG)
        DebugPrintf("Copied %s string '", string_type(m_s)); DebugPrint(m_s); DebugPrintf("' of length %d [assumed to be %d] to %s string '", (int)length(m_s), (int)m_l, string_type(t + 1)); DebugPrint(t + 1); DebugPrintf("' of length %d [assumed to be %d].\n", (int)length(t + 1), (int)l);
#endif
        return t;
      }
      return 0;
    }
  };
  template <typename S> pascal_string_intermediate<S> pascal_string(const S* s, size_t length) { return pascal_string_intermediate<S>{s, length}; }
  template <typename S> pascal_string_intermediate<S> pascal_string(const std::basic_string<S>& s) { return pascal_string_intermediate<S>{s.data(), s.length()}; }
  // Specialization for UTF16 -> UTF8.
  template <> template <> pascal_string_intermediate<wchar_t>::operator char* () const { return m_s ? (char*)pascal_string(::to_string(m_s, m_l)) : 0; }
  // Specialization for UTF8 -> UTF16.
  template <> template <> pascal_string_intermediate<char>::operator wchar_t* () const { return m_s ? (wchar_t*)pascal_string(::to_wstring(m_s, m_l)) : 0; }
}

XLOper4::XLOper4(const char* s) { if (0 != (val.str = pascal_string(s, length(s)))) xltype = (WORD)xltypeStr; }
XLOper4::XLOper4(const std::string& s) { if (0 != (val.str = pascal_string(s))) xltype = (WORD)xltypeStr; }
XLOper4::XLOper4(const wchar_t* s) { if (0 != (val.str = pascal_string(s, length(s)))) xltype = (WORD)xltypeStr; }
XLOper4::XLOper4(const std::wstring& s) { if (0 != (val.str = pascal_string(s))) xltype = (WORD)xltypeStr; }

XLOper4::XLOper4(const std::vector<double>& v) { new (this)XLOper4(v, 1, (WORD)v.size()); }

XLOper4::XLOper4(WORD rows, WORD columns) {
  val.array.rows = rows;
  val.array.columns = columns;
  val.array.lparray = new XLOper4[(size_t)rows * (size_t)columns];
  xltype = xltypeMulti;
}

XLOper4::XLOper4(WORD rows, WORD columns, const XLOper4& value) {
  val.array.rows = rows;
  val.array.columns = columns;
  val.array.lparray = new XLOper4[(size_t)rows * (size_t)columns];
  for (WORD i = 0; i < rows; ++i)
    for (WORD j = 0; j < columns; ++j)
      new (val.array.lparray + i * columns + j) XLOper4(value);
  xltype = xltypeMulti;
}

XLOper4::XLOper4(const std::vector<double>& v, WORD rows, WORD columns) {
  new (this)XLOper4(rows, columns);
  const INT32 m = (INT32)std::min((size_t)rows * (size_t)columns, (size_t)(v.size()));
  for (WORD i = 0, p = 0; i < rows && p < m; ++i)
    for (WORD j = 0; j < columns && p < m; ++j)
      new (val.array.lparray + i * columns + j) XLOper4(v[p++]);
}

XLOper4::XLOper4(int i) {
  xltype = xltypeInt;
  val.w = (short)i;
}

XLOper4::XLOper4(double x) {
  if (!std::isnan(x)) {
    xltype = xltypeNum;
    val.num = x;
  } else {
    xltype = xltypeErr;
    val.err = (WORD)XLErr::NaN;
  }
}

XLOper4::XLOper4(bool b) {
  xltype = xltypeBool;
  val.xbool = b;
}

XLOper4::XLOper4(XLType xlOperType) {
  xltype = (WORD)xlOperType;
}

XLOper4::XLOper4() {
  xltype = (WORD)XLType::Nil;
}

XLOper4::XLOper4(const XLOper4& other) {
  *this = other; // raw byte copy
  xltype &= ~(xlbitDLLFree | xlbitXLFree);
  switch (xltype) {
    case xltypeStr:
      if (val.str) val.str = pascal_string(val.str + 1, (size_t)val.str[0]);
      break;
    case xltypeMulti:
    {
      const size_t n = size();
      if (n) {
        XLOper4* src = (XLOper4*)val.array.lparray, * dst = (XLOper4*)(val.array.lparray = new XLOper4[n]);
        for (size_t i = 0; i < n; ++i) new (dst + i)XLOper4(src[i]);
      } else
        val.array.lparray = 0;
      break;
    }
    case xltypeBigData: xltype = (WORD)XLType::Nil; // not supported
      break;
  }
}

XLOper4& XLOper4::operator=(const char* s) { this->~XLOper4(); new(this)XLOper4(s); return *this; }
XLOper4& XLOper4::operator=(const wchar_t* s) { this->~XLOper4(); new(this)XLOper4(s); return *this; }
XLOper4& XLOper4::operator=(const std::string& s) { this->~XLOper4(); new(this)XLOper4(s); return *this; }
XLOper4& XLOper4::operator=(const std::wstring& s) { this->~XLOper4(); new(this)XLOper4(s); return *this; }
XLOper4& XLOper4::operator=(double s) { this->~XLOper4(); new(this)XLOper4(s); return *this; }
XLOper4& XLOper4::operator=(int s) { this->~XLOper4(); new(this)XLOper4(s); return *this; }

size_t XLOper4::size() const {
  return ((xltype & ~(xlbitDLLFree | xlbitXLFree)) == xltypeMulti) ? (size_t)val.array.rows * (size_t)val.array.columns : 1;
}

size_t XLOper4::rows() const {
  return ((xltype & ~(xlbitDLLFree | xlbitXLFree)) == xltypeMulti) ? (size_t)val.array.rows : 1;
}

size_t XLOper4::columns() const {
  return ((xltype & ~(xlbitDLLFree | xlbitXLFree)) == xltypeMulti) ? (size_t)val.array.columns : 1;
}

XLOper4& XLOper4::operator[] (size_t k) {
  assert(k < size());
  return ((xltype & ~(xlbitDLLFree | xlbitXLFree)) == xltypeMulti) ? ((XLOper4*)val.array.lparray)[k] : *this;
}

const XLOper4& XLOper4::operator[] (size_t k) const {
  assert(k < size());
  return ((xltype & ~(xlbitDLLFree | xlbitXLFree)) == xltypeMulti) ? ((XLOper4*)val.array.lparray)[k] : *this;
}

XLOper4& XLOper4::operator() (size_t k, size_t l) { return (*this)[k * val.array.columns + l]; }

const XLOper4& XLOper4::operator() (size_t k, size_t l) const { return (*this)[k * val.array.columns + l]; }

XLOper12::XLOper12(const char* s) { if (0 != (val.str = pascal_string(s, length(s)))) xltype = (WORD)xltypeStr; }
XLOper12::XLOper12(const std::string& s) { if (0 != (val.str = pascal_string(s))) xltype = (WORD)xltypeStr; }
XLOper12::XLOper12(const wchar_t* s) { if (0 != (val.str = pascal_string(s, length(s)))) xltype = (WORD)xltypeStr; }
XLOper12::XLOper12(const std::wstring& s) { if (0 != (val.str = pascal_string(s))) xltype = (WORD)xltypeStr; }

XLOper12::XLOper12(const std::vector<double>& v) { new (this)XLOper12(v, 1, (COL)v.size()); }

XLOper12::XLOper12(RW rows, COL columns) {
  val.array.rows = rows;
  val.array.columns = columns;
  val.array.lparray = new XLOper12[rows * columns];
  xltype = xltypeMulti;
}

XLOper12::XLOper12(RW rows, COL columns, const XLOper12& value) {
  val.array.rows = rows;
  val.array.columns = columns;
  val.array.lparray = new XLOper12[rows * columns];
  for (RW i = 0; i < rows; ++i)
    for (COL j = 0; j < columns; ++j)
      new (val.array.lparray + i * columns + j) XLOper12(value);
  xltype = xltypeMulti;
}

XLOper12::XLOper12(const std::vector<double>& v, RW rows, COL columns) {
  new (this)XLOper12(rows, columns);
  const INT32 m = std::min((INT32)rows * (INT32)columns, (INT32)(v.size()));
  for (RW i = 0, p = 0; i < rows && p < m; ++i)
    for (COL j = 0; j < columns && p < m; ++j)
      new (val.array.lparray + i * columns + j) XLOper12(v[p++]);
}

XLOper12::XLOper12(int i) {
  xltype = xltypeInt;
  val.w = i;
}

XLOper12::XLOper12(double x) {
  if (!std::isnan(x)) {
    xltype = xltypeNum;
    val.num = x;
  } else {
    xltype = xltypeErr;
    val.err = (int)XLErr::NaN;
  }
}

XLOper12::XLOper12(bool b) {
  xltype = xltypeBool;
  val.xbool = b;
}

XLOper12::XLOper12(XLType xlOperType) {
  xltype = (DWORD)xlOperType;
}

XLOper12::XLOper12() {
  xltype = (DWORD)XLType::Nil;
}

XLOper12::XLOper12(std::nullptr_t) { xltype = (DWORD)XLType::Nil; }

XLOper12::XLOper12(const XLOper12& other) {
  *this = other; // raw byte copy
  xltype &= ~(xlbitDLLFree | xlbitXLFree);
  switch (xltype) {
    case xltypeStr:
      if (val.str) val.str = pascal_string(val.str + 1, (size_t)val.str[0]);
      break;
    case xltypeMulti:
    {
      const size_t n = size();
      if (n) {
        XLOper12* src = (XLOper12*)other.val.array.lparray, * dst = (XLOper12*)(val.array.lparray = new XLOper12[n]);
        for (size_t i = 0; i < n; ++i) new (dst + i)XLOper12(src[i]);
      } else
        val.array.lparray = 0;
      break;
    }
    case xltypeBigData: xltype = (DWORD)XLType::Nil; // not supported
      break;
  }
}

XLOper12& XLOper12::operator=(const char* s) { this->~XLOper12(); new(this)XLOper12(s); return *this; }
XLOper12& XLOper12::operator=(const wchar_t* s) { this->~XLOper12(); new(this)XLOper12(s); return *this; }
XLOper12& XLOper12::operator=(const std::string& s) { this->~XLOper12(); new(this)XLOper12(s); return *this; }
XLOper12& XLOper12::operator=(const std::wstring& s) { this->~XLOper12(); new(this)XLOper12(s); return *this; }
XLOper12& XLOper12::operator=(double s) { this->~XLOper12(); new(this)XLOper12(s); return *this; }
XLOper12& XLOper12::operator=(int s) { this->~XLOper12(); new(this)XLOper12(s); return *this; }

size_t XLOper12::size() const {
  return ((xltype & ~(xlbitDLLFree | xlbitXLFree)) == xltypeMulti) ? (size_t)val.array.rows * (size_t)val.array.columns : 1;
}

size_t XLOper12::rows() const {
  return ((xltype & ~(xlbitDLLFree | xlbitXLFree)) == xltypeMulti) ? (size_t)val.array.rows : 1;
}

size_t XLOper12::columns() const {
  return ((xltype & ~(xlbitDLLFree | xlbitXLFree)) == xltypeMulti) ? (size_t)val.array.columns : 1;
}

XLOper12& XLOper12::operator[] (size_t k) {
  assert(k < size());
  return ((xltype & ~(xlbitDLLFree | xlbitXLFree)) == xltypeMulti) ? ((XLOper12*)val.array.lparray)[k] : *this;
}

const XLOper12& XLOper12::operator[] (size_t k) const {
  assert(k < size());
  return ((xltype & ~(xlbitDLLFree | xlbitXLFree)) == xltypeMulti) ? ((XLOper12*)val.array.lparray)[k] : *this;
}

XLOper12& XLOper12::operator() (size_t k, size_t l) { return (*this)[k * val.array.columns + l]; }

const XLOper12& XLOper12::operator() (size_t k, size_t l) const { return (*this)[k * val.array.columns + l]; }

namespace {

  const char* type_name(int type) {
    switch (type & ~(xlbitDLLFree | xlbitXLFree)) {
      case xltypeNum:     return "number";
      case xltypeStr:     return "string";
      case xltypeBool:    return "boolean";
      case xltypeErr:     return "#err";
      case xltypeMulti:   return "array";
      case xltypeMissing: return "missing";
      case xltypeNil:     return "nil";
      case xltypeInt:     return "integer";
      case xltypeBigData: return "binary";
      default: break;
    }
    return "<unknown>";
  }

  bool unchecked_to_bool(const XLOper& o) {
    switch (XLOper::APIMode) {
      case XLOper::APIMode4:  return o.o4.val.xbool != 0;
      case XLOper::APIMode12: return o.o12.val.xbool != 0;
      default: break;
    }
    return false;
  }

  int unchecked_to_int(const XLOper& o) {
    switch (XLOper::APIMode) {
      case XLOper::APIMode4:  return o.o4.val.w;
      case XLOper::APIMode12: return o.o12.val.w;
      default: break;
    }
    return -1;
  }

#if defined(VERBOSE_DEBUG_LOG)
  const char* xloper_version(const XLOper4&) { return "4"; }
  const char* xloper_version(const XLOper12&) { return "12"; }
#endif

  std::string to_string(double x) {
    const size_t buflen = 32;
    char buffer[buflen];
    snprintf(buffer, buflen, "%0.17g", x);
    return buffer;
  }

  std::wstring to_wstring(double x) {
    const size_t buflen = 32;
    wchar_t buffer[buflen];
    swprintf(buffer, buflen, L"%0.17g", x);
    return buffer;
  }

  std::string pascal_to_string(const char* pascal_string) { assert(0!=pascal_string); return std::string(pascal_string + 1, (size_t)(unsigned char)pascal_string[0]); }
  
  std::string pascal_to_string(const wchar_t* pascal_string) { assert(0 != pascal_string); return ::to_string(pascal_string + 1, (size_t)pascal_string[0]); }

  std::string unchecked_to_string(const XLOper& o) {
    switch (XLOper::APIMode) {
      case XLOper::APIMode4:  return pascal_to_string(o.o4.val.str);
      case XLOper::APIMode12: return pascal_to_string(o.o12.val.str);
      default: break;
    }
    return "unsupported API version (" + to_string((double)XLOper::APIMode) + ")";
  }

  std::wstring pascal_to_wstring(const char* pascal_string) { assert(0 != pascal_string); return ::to_wstring(pascal_string + 1, (size_t)(unsigned char)pascal_string[0]); }
  
  std::wstring pascal_to_wstring(const wchar_t* pascal_string) { assert(0 != pascal_string); return std::wstring(pascal_string + 1, (size_t)pascal_string[0]); }

  std::wstring unchecked_to_wstring(const XLOper& o) {
    switch (XLOper::APIMode) {
      case XLOper::APIMode4:  return pascal_to_wstring(o.o4.val.str);
      case XLOper::APIMode12: return pascal_to_wstring(o.o12.val.str);
      default: break;
    }
    return L"unsupported API version (" + to_wstring((double)XLOper::APIMode) + L")";
  }

  // Used for XLOper4 and XLOper12
  template <class X> void free_xloper_content(X* p) {
#if defined(VERBOSE_DEBUG_LOG)
    DebugPrintf("Freeing up XLOper%s (ptr=0x%p, type=0x%04X [%s]) of size %d.\n", xloper_version(*p), p, p->xltype, type_name(p->xltype), (int)sizeof(*p));
#endif
    // XLOpers that are populated by Excel via callbacks into it may contain secondary memory that we cannot release.
    // Excel signals this by setting the xlbitXLFree bit on such objects. See https://msdn.microsoft.com/en-us/library/office/bb687840.aspx.
    switch (p->xltype & ~(xlbitDLLFree | xlbitXLFree)) {
      case xltypeStr:
        if (p->val.str) { delete[] p->val.str; p->val.str = 0; }
        break;
      case xltypeMulti: // This is recursive
        if (p->val.array.lparray) { delete[](X*)p->val.array.lparray; p->val.array.lparray = 0; }
        break;
      default:
        break;
    }
    // The following statement is a safety precaution.
    p->xltype = xltypeNil;
  }

  inline const XLOper4*  pxloper4(const XLOper& o) { return (const XLOper4*)(&o.o4); }
  inline const XLOper12* pxloper12(const XLOper& o) { return (const XLOper12*)(&o.o12); }

  inline XLOper4*  pxloper4(XLOper* p) { return (XLOper4*)(&(p->o4)); }
  inline XLOper12* pxloper12(XLOper* p) { return (XLOper12*)(&(p->o12)); }

}

XLOper4::~XLOper4() { free_xloper_content(this); }

XLOper4* XLOper4::toExcel() {
#if defined(__GNUC__)
# pragma GCC diagnostic push
# pragma GCC diagnostic ignored "-Wnonnull-compare"
#endif
  if (0 != (void*)this)
#if defined(__GNUC__)
# pragma GCC diagnostic pop 
#endif
    this->xltype |= (WORD)xlbitDLLFree;
  return this;
}

XLOper12::~XLOper12() { free_xloper_content(this); }

XLOper12* XLOper12::toExcel() {
#if defined(__GNUC__)
# pragma GCC diagnostic push
# pragma GCC diagnostic ignored "-Wnonnull-compare"
#endif
  if (0 != (void*)this)
#if defined(__GNUC__)
# pragma GCC diagnostic pop 
#endif
    this->xltype |= xlbitDLLFree;
  return this;
}

XLOper12* XLOper12::Externalise() { return toExcel(); }

template <class A> void Construct(XLOper* p, const A& a) {
  switch (XLOper::APIMode) {
    case XLOper::APIMode4:  new (pxloper4(p))  XLOper4(a);    break;
    case XLOper::APIMode12: new (pxloper12(p)) XLOper12(a);   break;
    default:                memset((void*)p, -1, sizeof(*p)); break;
  }
}

#if defined(_MSC_VER)
// Warning C26495 Variable 'XLOper::<unnamed-tag>::o4' is uninitialized. Always initialize a member variable(type.6).
// Warning C26495 Variable 'XLOper::<unnamed-tag>::o12' is uninitialized. Always initialize a member variable(type.6).
# pragma warning(push)
# pragma warning(disable : 26495) // We are confident with what we are doing here.
#endif

XLOper::XLOper(const char* s) { Construct(this, s); }

XLOper::XLOper(const wchar_t* s) { Construct(this, s); }

XLOper::XLOper(const std::string& s) { Construct(this, s); }

XLOper::XLOper(const std::wstring& s) { Construct(this, s); }

XLOper::XLOper(const std::vector<double>& v) { Construct(this, v); }

XLOper::XLOper(int i) { Construct(this, i); }

XLOper::XLOper(double x) { Construct(this, x); }

XLOper::XLOper(bool b) { Construct(this, b); }

XLOper::XLOper(XLType t) { Construct(this, t); }

XLOper::XLOper() { Construct(this, XLType::Nil); }

XLOper::XLOper(size_t rows, size_t columns) {
  switch (APIMode) {
    case APIMode4:  new (pxloper4(this))  XLOper4((WORD)rows, (WORD)columns); break;
    case APIMode12: new (pxloper12(this)) XLOper12((RW)rows, (COL)columns);   break;
    default: break;
  }
}

XLOper::XLOper(size_t rows, size_t columns, const XLOper& value) {
  switch (APIMode) {
    case APIMode4:  new (pxloper4(this))  XLOper4((WORD)rows, (WORD)columns, *pxloper4(value)); break;
    case APIMode12: new (pxloper12(this)) XLOper12((RW)rows, (COL)columns, *pxloper12(value));  break;
    default: break;
  }
}

XLOper::XLOper(const std::vector<double>& v, size_t rows, size_t columns) {
  switch (APIMode) {
    case APIMode4:  new (pxloper4(this))  XLOper4(v, (WORD)rows, (WORD)columns); break;
    case APIMode12: new (pxloper12(this)) XLOper12(v, (RW)rows, (COL)columns);   break;
    default: break;
  }
}

XLOper::XLOper(const XLOper& o) {
  switch (APIMode) {
    case APIMode4:  new (pxloper4(this))  XLOper4(*pxloper4(o));   break;
    case APIMode12: new (pxloper12(this)) XLOper12(*pxloper12(o)); break;
    default: break;
  }
}
#if defined(_MSC_VER)
# pragma warning(pop) // #pragma warning(disable : 26495)
#endif

XLOper::~XLOper() {
  switch (APIMode) {
    case APIMode4:  pxloper4(this)->~XLOper4();   break;
    case APIMode12: pxloper12(this)->~XLOper12(); break;
    default: break;
  }
}

XLOper& XLOper::operator=(const char* s) { this->~XLOper(); new(this)XLOper(s); return *this; }
XLOper& XLOper::operator=(const wchar_t* s) { this->~XLOper(); new(this)XLOper(s); return *this; }
XLOper& XLOper::operator=(const std::string& s) { this->~XLOper(); new(this)XLOper(s); return *this; }
XLOper& XLOper::operator=(const std::wstring& s) { this->~XLOper(); new(this)XLOper(s); return *this; }
XLOper& XLOper::operator=(double s) { this->~XLOper(); new(this)XLOper(s); return *this; }
XLOper& XLOper::operator=(int s) { this->~XLOper(); new(this)XLOper(s); return *this; }

XLType XLOper::type() const {
  switch (APIMode) {
    case APIMode4:  return (XLType)o4.xltype;
    case APIMode12: return (XLType)o12.xltype;
    default: break;
  }
  return (XLType)-1;
}

const char* XLOper::type_name() const { return ::type_name((int)type()); }

XLErr XLOper::err() const {
  switch (APIMode) {
    case APIMode4:  return (XLErr)o4.val.err;
    case APIMode12: return (XLErr)o12.val.err;
    default: break;
  }
  return (XLErr)-1;
}

XLOper* XLOper::toExcel() {
#if defined(VERBOSE_DEBUG_LOG)
  DebugPrintf("toExcel() on XLOper of type=0x%04X (%s) of size %d with content '%ls'.\n", (unsigned int)type(), type_name(), (int)sizeof(*this), this->operator std::wstring().c_str());
#endif
#if defined(__GNUC__)
# pragma GCC diagnostic push
# pragma GCC diagnostic ignored "-Wnonnull-compare"
#endif
  if (0 != (void*)this)
#if defined(__GNUC__)
# pragma GCC diagnostic pop 
#endif
    switch (XLOper::APIMode) {
      case APIMode4:  o4.xltype |= (WORD)xlbitDLLFree; break;
      case APIMode12: o12.xltype |= xlbitDLLFree; break;
      default: break;
    }
  return this;
}

size_t XLOper::size() const {
  switch (APIMode) {
    case APIMode4:  return pxloper4(*this)->size();
    case APIMode12: return pxloper12(*this)->size();
    default: break;
  }
  return 0;
}

size_t XLOper::rows() const {
  switch (APIMode) {
    case APIMode4:  return pxloper4(*this)->rows();
    case APIMode12: return pxloper12(*this)->rows();
    default: break;
  }
  return 0;
}

size_t XLOper::columns() const {
  switch (APIMode) {
    case APIMode4:  return pxloper4(*this)->columns();
    case APIMode12: return pxloper12(*this)->columns();
    default: break;
  }
  return 0;
}

XLOper& XLOper::operator[] (size_t k) {
  switch (APIMode) {
    case APIMode4:  return *(XLOper*)&pxloper4(this)->operator[](k);
    case APIMode12: return *(XLOper*)&pxloper12(this)->operator[](k);
    default: break;
  }
  return *this;
}

const XLOper& XLOper::operator[] (size_t k) const {
  switch (APIMode) {
    case APIMode4:  return *(XLOper*)&pxloper4(*this)->operator[](k);
    case APIMode12: return *(XLOper*)&pxloper12(*this)->operator[](k);
    default: break;
  }
  return *this;
}

XLOper& XLOper::operator() (size_t k, size_t l) {
  switch (APIMode) {
    case APIMode4:  return *(XLOper*)&pxloper4(this)->operator()(k, l);
    case APIMode12: return *(XLOper*)&pxloper12(this)->operator()(k, l);
    default: break;
  }
  return *this;
}

const XLOper& XLOper::operator() (size_t k, size_t l) const {
  switch (APIMode) {
    case APIMode4:  return *(XLOper*)&pxloper4(*this)->operator()(k, l);
    case APIMode12: return *(XLOper*)&pxloper12(*this)->operator()(k, l);
    default: break;
  }
  return *this;
}

XLOper::operator std::string() const {
  int xltype = (int)type();
  switch (xltype & ~(xlbitDLLFree | xlbitXLFree)) {
    case xltypeNum:     return to_string(*(const double*)(this)); // field val.num is at offset 0 and of length 8 in API version 4 and 12 both for 32 bit and 64 bit.
    case xltypeStr:     return unchecked_to_string(*this);
    case xltypeInt:     return to_string((double)to_int());
    case xltypeBool:    return unchecked_to_bool(*this) ? "TRUE" : "FALSE";
    case xltypeErr:
    case xltypeMulti:
    case xltypeMissing:
    case xltypeNil:
    case xltypeBigData:
      return '<' + std::string(::type_name(xltype)) + '>';
  }
  return "unsupported XLOper type";
}

XLOper::operator std::wstring() const {
  int xltype = (int)type();
  switch (xltype & ~(xlbitDLLFree | xlbitXLFree)) {
    case xltypeNum:     return to_wstring(*(const double*)(this)); // field val.num is at offset 0 and of length 8 in API version 4 and 12 both for 32 bit and 64 bit.
    case xltypeStr:     return unchecked_to_wstring(*this);
    case xltypeInt:     return to_wstring((double)to_int());
    case xltypeBool:    return unchecked_to_bool(*this) ? L"TRUE" : L"FALSE";
    case xltypeErr:
    case xltypeMulti:
    case xltypeMissing:
    case xltypeNil:
    case xltypeBigData:
      const char* the_type = ::type_name(xltype);
      return L'<' + ::to_wstring(the_type, strlen(the_type)) + L'>';
  }
  return L"unsupported XLOper type";
}

int XLOper::to_int() const {
  int xltype = (int)type();
  switch (xltype & ~(xlbitDLLFree | xlbitXLFree)) {
    case xltypeNum:     return (int)*(const double*)(this); // field val.num is at offset 0 and of length 8 in API version 4 and 12 both for 32 bit and 64 bit.
    case xltypeInt:     return unchecked_to_int(*this);
    case xltypeBool:    return unchecked_to_bool(*this) ? 1 : 0;
    default:;
  }
  return -1;
}

bool XLOper::to_bool() const { return (bool)to_int(); }

XLOper::operator double() const {
  int xltype = (int)type();
  switch (xltype & ~(xlbitDLLFree | xlbitXLFree)) {
    case xltypeNum:     return *(const double*)(this); // field val.num is at offset 0 and of length 8 in API version 4 and 12 both for 32 bit and 64 bit.
    case xltypeInt:     return (double)unchecked_to_int(*this);
    case xltypeBool:    return unchecked_to_bool(*this) ? 1 : 0;
    default:;
  }
  return std::numeric_limits<double>::quiet_NaN();
}

std::vector<double> XLOper::to_vector() const { 
  int xltype = (int)type();
  switch (xltype & ~(xlbitDLLFree | xlbitXLFree)) {
    case xltypeNum:
    case xltypeInt:
    case xltypeBool:    return std::vector<double>(1, this->operator double());
    case xltypeMulti:
    {
      std::vector<double> v(size());
      for (size_t i = 0, n = size(); i < n; ++i)
        v[i] = (*this)[i];
      return v;
    }
    case xltypeMissing: return std::vector<double>();
    default:;
  }
  return std::vector<double>(1, std::numeric_limits<double>::quiet_NaN());
}

XLOper::operator std::vector<double>() const { return to_vector(); }

char XLOper::xloper_code() { // 'P' in API version 4, and 'Q' in API version 12.
  switch (XLOper::APIMode) {
    case XLOper::APIMode4:  return 'P';
    case XLOper::APIMode12: return 'Q';
    default: break;
  }
  return '?';
}

std::string XLOper::replace_xloper_codes(const char* codes) { // '?' -> 'P' in API version 4, and 'Q' in API version 12.
  int n = (int)strlen(codes);
  std::vector<char> s(n + 1);
  const char o = xloper_code();
  const bool permit_thread_safety_char = XLOper::APIMode4 != XLOper::APIMode;
  for (int i = 0, j = 0; i < n; ++i)
    if (permit_thread_safety_char || codes[i] != '$')
      s[j++] = '?' == codes[i] ? o : codes[i];
  return &s[0];
}

XLOper missing() {
  //
  // The XLOper4  type member 'xltype' is at offset 8 in 32-bit and at offset 16 in 64 bit, and of length 2 [all in bytes].
  // The XLOper12 type member 'xltype' is at offset 24 both in 32-bit and in 64 bit, and of length 4.
  // We can therefore set both of these fields to value 'xltypeMissing' since they don't overlap in memory, which works for both API versions.
  //
  XLOper missing;
  missing.o4.xltype = (WORD)XLType::Missing;
  missing.o12.xltype = (DWORD)XLType::Missing;
  return missing;
}

const XLOper XLOper::Missing = missing();

XLOper XLOper::NotApplicable() {
  // The XLOper4  type member 'xltype' is at offset 8 in 32-bit and at offset 16 in 64 bit, and of length 2 [all in bytes].
  // The XLOper12 type member 'xltype' is at offset 24 both in 32-bit and in 64 bit, and of length 4.
  // We can therefore set both of these fields since they don't overlap in memory, which works for both API versions.
  XLOper na;
  na.o4.xltype = (WORD)XLType::Err;
  na.o12.xltype = (DWORD)XLType::Err;
  switch (XLOper::APIMode) {
    case XLOper::APIMode4:  na.o4.val.err = (WORD)XLErr::NotApplicable; break;
    case XLOper::APIMode12: na.o12.val.err = (int)XLErr::NotApplicable; break;
    default: break;
  }
  return na;
}

bool XLOper::isUndefined() const {
  const int t = (int)type() & ~(xlbitDLLFree | xlbitXLFree);
  return  xltypeMissing == t || xltypeNil == t || (xltypeErr == t && XLErr::NotApplicable == err());
}

template <> int XLOper::to<int>() const { return to_int(); }
template <> bool XLOper::to<bool>() const { return to_bool(); }
template <> double XLOper::to<double>() const { return this->operator double(); }
template <> std::string XLOper::to<std::string>() const { return this->operator std::string(); }
