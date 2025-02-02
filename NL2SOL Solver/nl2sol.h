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

#ifndef NL2SOL_H
#define NL2SOL_H

#include <vector>
#include <stddef.h>

class NL2SOL {

   std::vector<double> m_workspace;
   std::vector<int> m_iv;

public:

   struct Defaults {
      static const double RelativeTolerance; // = √DBL_EPSILON, used both for x (the argument vector) and f (the function vector), as well as the finite differencing width.
      static const double AbsoluteTolerance; // = DBL_EPSILON/2
      static const size_t MaximumNumberOfIterations; // 10·DBL_DIG
   };

   // INPUT (but can be changed before a re-run)
   double ArgumentTolerance, FunctionTolerance, AbsoluteFunctionTolerance, FiniteDifferencingWidth; // for Jacobian computation
   size_t MaximumNumberOfIterations;

   // OUTPUT
   double AchievedAccuracy, RMS;
   int ReturnCode, Iterations, FunctionEvaluations, JacobianEvaluations;
   
#if defined( DEFINE_NL2SOL_RETURN_CODE_ENUM )
   enum ReturnCodeEnum {
       UNKNOWN = 0,
       ARGUMENT_CONVERGENCE = 3,
       RELATIVE_FUNCTION_CONVERGENCE = 4,
       ARGUMENT_AND_RELATIVE_FUNCTION_CONVERGENCE = 5,
       ABSOLUTE_FUNCTION_CONVERGENCE = 6,
       SINGULAR_CONVERGENCE = 7,
       FALSE_CONVERGENCE = 8,
       FUNCTION_EVALUATION_LIMIT_EXCEEDED = 9,
       ITERATION_LIMIT_EXCEEDED = 10,
       USER_INTERRUPT = 11,
       INVALID_STARTING_COORDINATES = 13,
       BAD_PARAMETERS = 14,
       UNABLE_TO_COMPUTE_JACOBIAN = 15,
       UNDERDETERMINED_PROBLEM = 16,
       RESTART_ATTEMPTED_WITH_CHANGED_DIMENSIONS = 17,
       IV_INITS_OUT_OF_RANGE = 18
   };
#endif

   class Functor {
   public:
      // The return value FALSE indicates that x[] was not inside the admissible domain, i.e., that constraints are violated.
      virtual bool operator()(const double *x, size_t n_x, double *f /* f=f(x) */, size_t n_f) = 0;
      // The Jacobian function may be left undefined in which case finite-differencing is used.
      virtual bool have_jacobian() const { return false; }
      // The return value FALSE indicates that x[] was not inside the admissible domain, i.e., that constraints are violated.
      // NOTE:
      //  The Jacobian function must store the output matrix J in the output vector J[] in column-major order (ForTran).
      //  For example, if n_f is the length of f[], and n_x is the length of x[], then J(k,l) is to be stored in J[k+l*n_f].
      //  The element J(k,l) represents  the derivative
      //                    ∂ f[k]
      //                    ------
      //                    ∂ x[l] .
      virtual bool jacobian(const double *x, size_t n_x, size_t n_f,
         double *J /* J[k+l*n_f] = ∂f[k]/∂x[l]. The length of J[] must be at least n_f*n_x */){ return false; }
      // The return value TRUE indicates a request to abort the current calculation.
      virtual bool stop(){ return false; }
   };

   //
   // Below is a quote from section 15.5, page 688 of "Numerical Recipes in C",
   // W. H. Press, S. A. Teukolsky, W. t. Vetterling, B. P. Flannery,
   // second edition, ISBN 0-521-43108-5, which might still reside at:
   //
   //    http://www.library.cornell.edu/nr/cbookcpdf.html
   //
   //   The Levenberg-Marquardt algorithm can be implemented as a model-trust
   //   region method for minimization (see x 9.7 and ref. [2]) applied to the special case
   //   of a least squares function. A code of this kind due to Moré [3] can be found in
   //   MINPACK [4] . Another algorithm for nonlinear least-squares keeps the second-derivative
   //   term we dropped in the Levenberg-Marquardt method whenever it would
   //   be better to do so. These methods are called “full Newton-type” methods and are
   //   reputed to be more robust than Levenberg-Marquardt, but more complex. One
   //   calculator is the code NL2SOL [5] .
   //
   //   CITED REFERENCES AND FURTHER READING:
   //   Bevington, P.R. 1969, Data Reduction and Error Analysis for the Physical Sciences(New York:
   //     McGraw-Hill), Chapter 11.
   //   Marquardt, D.W. 1963, Journal of the Society for Industrial and Applied Mathematics, vol. 11,
   //     pp. 431–441. [1]
   //   Jacobs, D.A.H. (ed.) 1977, The State of the Art in Numerical Analysis (London: Academic
   //     Press), Chapter III.2 (by J.E. Dennis).
   //   Dennis, J.E., and Schnabel, R.B. 1983, Numerical Methods for Unconstrained Optimization and
   //     Nonlinear Equations (Englewood Cliffs, NJ: Prentice-Hall). [2]
   //   Moré, J.J. 1977, in Numerical Analysis, Lecture Notes in Mathematics, vol. 630, G.A. Watson,
   //     ed. (Berlin: Springer-Verlag), pp. 105–116. [3]
   //   Moré, J.J., Garbow, B.S., and Hillstrom, K.E. 1980, User Guide for MINPACK-1, Argonne National
   //     Laboratory Report ANL-80-74. [4]
   //   Dennis, J.E., Gay, D.M, and Welsch, R.E. 1981, ACM Transactions on Mathematical Software,
   //     vol. 7, pp. 348–368; op. cit., pp. 369–383. [5].
   //
   //
   //

   NL2SOL(double argumentTolerance = Defaults::RelativeTolerance, // This is a *relative* convergence tolerance.
      double functionTolerance = Defaults::RelativeTolerance, // This is a *relative* convergence tolerance.
      double absoluteFunctionTolerance = Defaults::AbsoluteTolerance, // 0.5*functionTolerance*functionTolerance
      double finiteDifferencingWidth = Defaults::RelativeTolerance, // for Jacobian computation, ignored when jacobian function is given
      size_t maximumNumberOfIterations = Defaults::MaximumNumberOfIterations);

   //
   // x is the argument vector, n_x is its length
   // f=f(x) is the function vector, n_f is its length
   // Note : The function below uses half the sum of the squares of the function vector elements as its internal
   // error measure which is compared against the above 'functionTolerance'. Upon completion, m_iterations is the
   // number of iterations it took. A negative value for m_iterations indicates failure to converge.
   // The return value is the Euclidean norm of the function vector at the returned x[] coordinates.
   double tryToSolve(Functor& functor, double *x, size_t n_x, size_t n_f);
   double tryToSolve(Functor& functor, std::vector<double> &x, size_t n_f);

   static bool is_number(double x);
   static bool any_not_is_number(const double *x, size_t n_x);
   static bool any_not_is_number(const std::vector<double>&x);

};

#endif // NL2SOL_H
