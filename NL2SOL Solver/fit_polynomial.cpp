//
// ======================================================================================
// Copyright © 2024 Peter Jäckel of OTC Analytics.
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

#include "nl2sol.h"

#include "dllmain.h"

#include <cmath>
#include <assert.h>
#include <float.h>
//#include <numeric>
//#include <functional>
#include <memory>
#include <stdexcept>

namespace {
  double evaluate_polynomial(double x, const double* p, size_t n) {
    double y = 0;
    if (n > 0) {
      y = p[--n];
      for (; n > 0;) {
        y *= x;
        y += p[--n];
      }
    }
    return y;
  }
  double evaluate_polynomial(double x, const std::vector<double>& coeffs) {
    return evaluate_polynomial(x, &coeffs[0], coeffs.size());
  }
  class PolynomialFit : public NL2SOL::Functor {
    const std::vector<double> m_abscissas, m_ordinates;
    const bool m_using_jacobian;
  public:
    PolynomialFit(const std::vector<double>& abscissas, const std::vector<double>& ordinates, bool use_jacobian) :
      m_abscissas(abscissas), m_ordinates(ordinates), m_using_jacobian(use_jacobian) {
      if (0 == abscissas.size())
        throw std::runtime_error("Please provide some data (0 points given).");
      if (ordinates.size() != abscissas.size())
        throw std::runtime_error("Please provide equal abscissa and ordinate values (" + std::to_string(abscissas.size()) + " abscissas given, " + std::to_string(ordinates.size()) + " ordinates given).");
    }
    size_t data_length() const { return m_abscissas.size(); }
    virtual bool operator()(const double* coeffs, size_t n_x, double* func_vec /* f=f(x) */, size_t n_f) {
      assert(m_abscissas.size() == n_f);
      assert(m_ordinates.size() == n_f);
      if (NL2SOL::any_not_is_number(coeffs, n_x)) return false;
      for (size_t i = 0; i < n_f; ++i)
        func_vec[i] = evaluate_polynomial(m_abscissas[i], coeffs, n_x) - m_ordinates[i];
      return !NL2SOL::any_not_is_number(func_vec, n_f);
    }
    virtual bool have_jacobian() const { return m_using_jacobian; }
    virtual bool jacobian(const double* coeffs, size_t n_x, size_t n_f,
                          double* jac /* jac[i+j*n_f] = ∂f[i]/∂x[j]. The length of jac[] must be at least n_f*n_x */) {
      assert(m_abscissas.size() == n_f);
      assert(m_ordinates.size() == n_f);
      if (n_x > 0) {
        if (NL2SOL::any_not_is_number(coeffs, n_x)) return false;
        for (size_t i = 0; i < n_f; ++i) {
          const double x_i = m_abscissas[i];
          double f_ij = jac[i] = 1;
          for (size_t j = 1; j < n_x; ++j)
            jac[i + j * n_f] = (f_ij *= x_i);
        }
      }
      return !NL2SOL::any_not_is_number(jac, n_f * n_x);
    }
  };
}

#include "XLFunctions.h"

extern "C" DLL_EXPORT XLOper * fit_polynomial(const XLOper & order, const XLOper & abscissas, const XLOper & ordinates, const XLOper & use_jacobian) {
  try {
    int m = order.to_int() + 1;
    if (m < 1)
      throw std::runtime_error("Please specify a non-negative polynomial order.");

    PolynomialFit poly(abscissas.to_vector(), ordinates.to_vector(), use_jacobian | false);

    std::vector<double> coeffs(m);
    NL2SOL nl2sol;
    nl2sol.tryToSolve(poly, coeffs, poly.data_length());

    std::unique_ptr<XLOper> pres(new XLOper(4, std::max(2, m), XLOper::NotApplicable()));
    XLOper& result = *pres;
    const auto na = XLOper::NotApplicable();
    for (int i = 0; i < m; ++i)
      result(0, i) = coeffs[i];
    result(1, 0) = "RMS:";
    result(1, 1) = nl2sol.RMS;
    result(2, 0) = "Functor evaluations:";
    result(2, 1) = nl2sol.FunctionEvaluations;
    result(3, 0) = poly.have_jacobian()? "Jacobian evaluations:": "Jacobian fin. diff. approximations:";
    result(3, 1) = nl2sol.JacobianEvaluations;
    return pres.release()->toExcel();

  }
  CATCH_TO_EXCEL
}

DECLARE_XL_FUNCTION(fit_polynomial, "", "?????$", "poly_order,abscissas,ordinates,use_jacobian",
                    "returns the coefficients of a polynomial fit to the given input data.",
                    {"a positive integer","the abscissa values of the data points","the ordinate values of the data points","[false] an optional boolean to indicate whether to use an analytical Jacobian"});

extern "C" DLL_EXPORT XLOper * polynomial(const XLOper & coeffs, const XLOper & abscissas) {
  try {
    const auto coefficients = coeffs.to_vector();
    const auto x = abscissas.to_vector();
    std::unique_ptr<XLOper> pres(new XLOper(abscissas));
    XLOper& result = *pres;
    for (size_t i = 0; i < x.size(); ++i)
      result[i] = evaluate_polynomial(x[i], coefficients);
    return pres.release()->toExcel();
  }
  CATCH_TO_EXCEL
}

DECLARE_XL_FUNCTION(polynomial, "", "???$", "coefficients,abscissas",
                    "returns the evaluation of the polynomial given by the coefficients at the given input abscissas.",
                    { "a vector of polynomial coefficients","the abscissa values where to evaluate the polynomial" });
