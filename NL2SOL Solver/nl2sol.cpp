//
// ======================================================================================
// 2024 Peter Jäckel of OTC Analytics.
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

#include "nl2sol.h"
#include <cmath>
#include <string>
#include <cfloat>
#include <stdio.h>
#include <algorithm>
#include <stdexcept>
#include <limits>
#include <float.h>
#include <string.h>

//
// ForTran source code is provided in comments at the end of this file
//

#define  STRINGIFY(x)    #x
#define  TOSTRING(x)     STRINGIFY(x)
#define  __FILE__LINE__  __FILE__ ":" TOSTRING(__LINE__)

#ifdef _MSC_VER
# define snprintf(buffer,sizeOfBuffer,...) _snprintf_s((buffer),(sizeOfBuffer),_TRUNCATE,__VA_ARGS__)
namespace {
   inline bool isInf(double x){ return !_finite(x); }
   inline bool isNaN(double x){ return 0!=_isnan(x); }
}
# define finite(x) _finite(x)
# define isnan(x) _isnan(x)
# ifndef _DEBUG
#  ifndef NDEBUG
#   define NDEBUG
#  endif
# endif
#else
namespace {
   inline bool isInf(double x){ return std::isinf(x); }
   inline bool isNaN(double x){ return std::isnan(x); }
}
#endif

#ifdef NDEBUG
#  define ASSERT(condition)
#else
#  define ASSERT(condition) if (!(condition)) throw std::runtime_error(__FILE__LINE__": assertion failure. The violated condition is: "#condition)
#endif

#define REQUIRE(condition,description) if (!(condition)) throw std::runtime_error(description)
#define INSIST(condition) if (!(condition)) throw std::runtime_error("Runtime violation of condition: "#condition)

namespace {

  typedef int(*i_fp)(int *);
  typedef void(*S_fp)(int *, int *, double *, int *, double *, int *uiparm, double *urparm, i_fp);

  template <class T> void fill(std::vector< T > & v, const T & x) { std::fill(v.begin(), v.end(), x); }

  bool isFinite(double x) { return 0 != finite(x); }
  std::string toString(double x) {
    char buffer[128];
    snprintf(buffer, sizeof(buffer), "%.16g", x);
    return buffer;
  }
  int assess(double *, int *, int *, double *, double *, double *, double *, double *);
  double dotprd(int *, double *, double *);
  int linvrt(int *, double *, double *);
  int lsqrt(int, int *, double *, double *, int *);
  int ltsqar(int *, double *, double *);
  int nl2sol_defaults(int *iv, double *v);
  int parchk(int *, int *, int *, int *, double *);
  int dupdat(double *, int *, double *, int *, int *, int *, double *);
  double reldst(int *, double *, double *, double *);
  double rmdcon(int k);
  int rptmul(int, int *, double *, int *, int *, double *, double *, double *, double *);
  int slupdt(double *, double *, int *, double *, double *, double *, double *, double *, double *, double *);
  int slvmul(int *, double *, double *, double *);
  double v2norm(int *, double *);
  int vaxpy(int *, double *, double, double *, double *);
  int vcopy(int, double *, double *);
  int vscopy(int *, double *, double);
  int nl2itr(NL2SOL::Functor &functor, double *d_, int *iv, double *j, int *n, int *nn, int *p, double *r_, double *v, double *x);
  int qrfact(int *nm, int *m, int *n, double *qr, double *alpha, int *ipivot, int *ierr, int nopivk, double *sum);
  int qapply(int *nn, int *n, int *p, double *j, double *r_, int *ierr);
  int lmstep(double *d_, double *g, int *ierr, int *ipivot, int *ka, int *p, double *qtr, double *r_, double *step, double *v, double *w);
  int gqtstp(double *d_, double *dig, double *dihdi, int *ka, double *l, int *p, double *step, double *v, double *w);
  int litvmu(int *, double *, double *, double *);
  int livmul(int *, double *, double *, double *);
  int lmstep(double *, double *, int *, int *, int *, int *, double *, double *, double *, double *, double *);
  double lsvmin(int *, double *, double *, double *);
  int covclc(int *, double *, int *, double *, int *, int *, int *, double *, double *, double *);
  int qrfact(int *, int *, int *, double *, double *, int *, int *, int, double *);
  int qapply(int *, int *n, int *p, double *, double *, int *);
  int nl2sol(NL2SOL::Functor &functor, int *function_evaluations, int *jacobian_evaluations, int *n, int *p, double *x, int *iv, double *v);
  int nl2sno(NL2SOL::Functor &functor, int *function_evaluations, int *jacobian_evaluations, int *n, int *p, double *x, int *iv, double *v);

  int assess(double *d_, int *iv, int *p, double *step, double *stlstg, double *v, double *x, double *x0) {
    /* Initialized data */

    const double half = .5;
    const double one = 1.;
    const double two = 2.;
    const double zero = 0.;
    const int irc = 3;
    const int mlstgd = 4;
    const int model = 5;
    const int nfcall = 6;
    const int nfgcal = 7;
    const int radinc = 8;
    const int restor = 9;
    const int stage = 10;
    const int stglim = 11;
    const int switch_ = 12;
    const int toobig = 2;
    const int xirc = 13;
    const int afctol = 31;
    const int decfac = 22;
    const int dstnrm = 2;
    const int dst0 = 3;
    const int dstsav = 18;
    const int f = 10;
    const int fdif = 11;
    const int flstgd = 12;
    const int f0 = 13;
    const int gtslst = 14;
    const int gtstep = 4;
    const int incfac = 23;
    const int lmax0 = 35;
    const int nreduc = 6;
    const int plstgd = 15;
    const int preduc = 7;
    const int radfac = 16;
    const int rdfcmn = 24;
    const int rdfcmx = 25;
    const int reldx = 17;
    const int rfctol = 32;
    const int stppar = 5;
    const int tuner1 = 26;
    const int tuner2 = 27;
    const int tuner3 = 28;
    const int xctol = 33;
    const int xftol = 34;

    /* System generated locals */
    int i_1;
    double d_1, d_2;

    /* Local variables */
    double emax, xmax, rfac1;
    int i_;
    bool goodx;
    double reldx1;
    int nfc;
    double gts;


    //  ***  assess candidate step (nl2sol version 2.2)  *** 


    //  ***  purpose  *** 

    //        this subroutine is called by an unconstrained minimization 
    //             of decreases) so far this iteration. 
    // iv(restor) (out) set to 0 unless x and v(f) have been restored, in 
    //             which case assess sets iv(restor) = 1. 
    //  iv(stage) (i/o) count of the number of models tried so far in the 
    //             current iteration. 
    // iv(stglim) (in)  maximum number of models to consider. 
    // iv(switch) (out) set to 0 unless a new model is being tried and it 
    //             gives a smaller function value than the previous model, 
    //             in which case assess sets iv(switch) = 1. 
    // iv(toobig) (in)  is nonzero if step was too big (e.g. if it caused 
    //             overflow). 
    //   iv(xirc) (i/o) value that iv(irc) would have in the absence of 
    //             convergence, false convergence, and oversized steps. 

    //  ***  v values referenced  *** 

    // v(afctol) (in)  absolute function convergence tolerance.  if the 
    //             absolute value of the current function value v(f) is less 
    //             than v(afctol), then assess returns with iv(irc) = 10. 
    // v(decfac) (in)  factor by which to decrease radius when iv(toobig) is 
    //             nonzero. 
    // v(dstnrm) (in)  the 2-norm of d*step. 
    // v(dstsav) (i/o) value of v(dstnrm) on saved step. 
    //   v(dst0) (in)  the 2-norm of d times the newton step (when defined, 
    //             i.e., for v(nreduc) .ge. 0). 
    //      v(f) (i/o) on both input and output, v(f) is the objective func- 
    //             tion value at x.  if x is restored to a previous value, 
    //             then v(f) is restored to the corresponding value. 
    //   v(fdif) (out) the function reduction v(f0) - v(f) (for the output 
    //             value of v(f) if an earlier step gave a bigger function 
    //             decrease, and for the input value of v(f) otherwise). 
    // v(flstgd) (i/o) saved value of v(f). 
    //     v(f0) (in)  objective function value at start of iteration. 
    // v(gtslst) (i/o) value of v(gtstep) on saved step. 
    // v(gtstep) (in)  inner product between step and gradient. 
    // v(incfac) (in)  minimum factor by which to increase radius. 
    //  v(lmax0) (in)  maximum reasonable step size (and initial step bound). 
    //             what was predicted, if a return with iv(irc) = 7, 8, 9, 
    //             or 10 does not occur, if v(dstnrm) .gt. v(lmax0), and if 
    //             v(preduc) .le. v(rfctol) * abs(v(f0)), then assess re- 
    //             turns with iv(irc) = 11.  if so doing appears worthwhile, 
    //             then assess repeats this test with v(preduc) computed for 
    //             a step of length v(lmax0) (by a return with iv(irc) = 6). 
    // v(nreduc) (i/o)  function reduction predicted by quadratic model for 
    //             newton step.  if assess is called with iv(irc) = 6, i.e., 
    //             use in the singular convervence test, then v(nreduc) is 
    //             set to -v(preduc) before the latter is restored. 
    // v(plstgd) (i/o) value of v(preduc) on saved step. 
    // v(preduc) (i/o) function reduction predicted by quadratic model for 
    //             current step. 
    // v(radfac) (out) factor to be used in determining the new radius, 
    //             which should be v(radfac)*dst, where  dst  is either the 
    //             output value of v(dstnrm) or the 2-norm of 
    //             diag(newd)*step  for the output value of step and the 
    //             updated version, newd, of the scale vector d.  for 
    //             iv(irc) = 3, v(radfac) = 1.0 is returned. 
    // v(rdfcmn) (in)  minimum value for v(radfac) in terms of the input 
    //             value of v(dstnrm) -- suggested value = 0.1. 
    // v(rdfcmx) (in)  maximum value for v(radfac) -- suggested value = 4.0. 
    //  v(reldx) (out) scaled relative change in x caused by step, computed 
    //             by function  reldst  as 
    //                 max (d(i)*abs(x(i)-x0(i)), 1 .le. i .le. p) / 
    //                    max (d(i)*(abs(x(i))+abs(x0(i))), 1 .le. i .le. p). 
    //             puted using the output (possibly restored) values of x 
    //             and step.  otherwise it is computed using the input 
    //             values. 
    // iv(mxiter) 
    // v(rfctol) (in)  relative function convergence tolerance.  if the 
    //             actual function reduction is at most twice what was pre- 
    //             dicted and  v(nreduc) .le. v(rfctol)*abs(v(f0)),  then 
    //             assess returns with iv(irc) = 8 or 9.  see also v(lmax0). 
    // v(stppar) (in)  marquardt parameter -- 0 means full newton step. 
    // v(tuner1) (in)  tuning constant used to decide if the function 
    //             reduction was much less than expected.  suggested 
    //             value = 0.1. 
    // v(tuner2) (in)  tuning constant used to decide if the function 
    //             reduction was large enough to accept step.  suggested 
    //             value = 10**-4. 
    // v(tuner3) (in)  tuning constant used to decide if the radius 
    //             should be increased.  suggested value = 0.75. 
    //  v(xctol) (in)  x-convergence criterion.  if step is a newton step 
    //             (v(stppar) = 0) having v(reldx) .le. v(xctol) and giving 
    //             at most twice the predicted function decrease, then 
    //             assess returns iv(irc) = 7 or 9. 
    //  v(xftol) (in)  false convergence tolerance.  if step gave no or only 
    //             a small function decrease and v(reldx) .lt. v(xftol), 
    //             then assess returns with iv(irc) = 12. 

    // -------------------------------  notes  -------------------------------

    //  ***  application and usage restrictions  *** 
    //        this routine is called as part of the nl2sol (nonlinear 
    //     least-squares) package.  it may be used in any unconstrained 
    //     minimization solver that uses dogleg, goldfeld-quandt-trotter, 
    //     or levenberg-marquardt steps. 

    //  ***  algorithm notes  *** 

    //        see (1) for further discussion of the assessing and model 
    //     switching strategies.  while nl2sol considers only two models, 
    //     assess is designed to handle any number of models. 

    //  ***  usage notes  *** 

    //        on the first call of an iteration, only the i/o variables 
    //     step, x, iv(irc), iv(model), v(f), v(dstnrm), v(gtstep), and 
    //     v(preduc) need have been initialized.  between calls, no i/o 
    //     values execpt step, x, iv(model), v(f) and the stopping toler- 
    //     ances should be changed. 
    //        after a return for convergence or false convergence, one can 
    //     change the stopping tolerances and call assess again, in which 
    //     case the stopping tests will be repeated. 

    //  ***  references  *** 

    //     (1) dennis, j.e., jr., gay, d.m., and welsch, r.e. (1981), 
    //        an adaptive nonlinear least-squares algorithm, 
    //        acm trans. math. software, vol. 7, no. 3. 

    //     (2) powell, m.j.d. (1970)  a fortran subroutine for solving 
    //        systems of nonlinear algebraic equations, in numerical 
    //        methods for nonlinear algebraic equations, edited by 
    //        p. rabinowitz, gordon and breach, london. 

    //  ***  history  *** 

    //        john dennis designed much of this routine, starting with 
    //     ideas in (2). roy welsch suggested the model switching strategy. 
    //        david gay and stephen peters cast this subroutine into a more 
    //     portable form (winter 1977), and david gay cast it into its 
    //     present form (fall 1978). 

    //  ***  general  *** 

    //     this subroutine was written in connection with research 
    //     supported by the national science foundation under grants 
    //     mcs-7600324, dcr75-10143, 76-1431ss, mcs76-11989, and 
    //     mcs-7906671. 

    // ------------------------  external quantities  ------------------------

    //  ***  external functions and subroutines  *** 
    //
    // vcopy.... copies one vector to another. 

    // --------------------------  local variables  --------------------------

    /* Parameter adjustments */
    --iv;
    --x0;
    --x;
    --stlstg;
    --step;
    --d_;
    --v;

    /* Function Body */

    // +++++++++++++++++++++++++++++++  body  ++++++++++++++++++++++++++++++++


    nfc = iv[nfcall];
    iv[switch_] = 0;
    iv[restor] = 0;
    rfac1 = one;
    goodx = true;
    i_ = iv[irc];

    if (i_ >= 1 && i_ <= 12) {
      switch (static_cast<int>(i_)) {
      case 1:  goto L20;
      case 2:  goto L30;
      case 3:  goto L10;
      case 4:  goto L10;
      case 5:  goto L40;
      case 6:  goto L360;
      case 7:  goto L290;
      case 8:  goto L290;
      case 9:  goto L290;
      case 10:  goto L290;
      case 11:  goto L290;
      case 12:  goto L140;
      }
    }
    iv[irc] = 13;
    goto L999;

    //  ***  initialize for new iteration  *** 

  L10:
    iv[stage] = 1;
    iv[radinc] = 0;
    v[flstgd] = v[f0];
    if (iv[toobig] == 0) {
      goto L90;
    }
    iv[stage] = -1;
    iv[xirc] = i_;
    goto L60;

    //  ***  step was recomputed with new model or smaller radius  *** 
    //  ***  first decide which  *** 

  L20:
    if (iv[model] != iv[mlstgd]) {
      goto L30;
    }
    //        ***  old model retained, smaller radius tried  *** 
    //        ***  do not consider any more new models this iteration  *** 
    iv[stage] = iv[stglim];
    iv[radinc] = -1;
    goto L90;

    //  ***  a new model is being tried.  decide whether to keep it.  *** 

  L30:
    ++iv[stage];

    //     ***  now we add the possibiltiy that step was recomputed with  *** 
    //     ***  the same model, perhaps because of an oversized step.     *** 

  L40:
    if (iv[stage] > 0) {
      goto L50;
    }

    //        ***  step was recomputed because it was too big.  *** 

    if (iv[toobig] != 0) {
      goto L60;
    }

    //        ***  restore iv(stage) and pick up where we left off.  *** 

    iv[stage] = -iv[stage];
    i_ = iv[xirc];
    switch (static_cast<int>(i_)) {
    case 1:  goto L20;
    case 2:  goto L30;
    case 3:  goto L90;
    case 4:  goto L90;
    case 5:  goto L70;
    }

  L50:
    if (iv[toobig] == 0) {
      goto L70;
    }

    //  ***  handle oversize step  *** 

    if (iv[radinc] > 0) {
      goto L80;
    }
    iv[stage] = -iv[stage];
    iv[xirc] = iv[irc];

  L60:
    v[radfac] = v[decfac];
    --iv[radinc];
    iv[irc] = 5;
    goto L999;

  L70:
    if (v[f] < v[flstgd]) {
      goto L90;
    }

    //     *** the new step is a loser.  restore old model.  *** 

    if (iv[model] == iv[mlstgd]) {
      goto L80;
    }
    iv[model] = iv[mlstgd];
    iv[switch_] = 1;

    //     ***  restore step, etc. only if a previous step decreased v(f). 

  L80:
    if (v[flstgd] >= v[f0]) {
      goto L90;
    }
    iv[restor] = 1;
    v[f] = v[flstgd];
    v[preduc] = v[plstgd];
    v[gtstep] = v[gtslst];
    if (iv[switch_] == 0) {
      rfac1 = v[dstnrm] / v[dstsav];
    }
    v[dstnrm] = v[dstsav];
    nfc = iv[nfgcal];
    goodx = false;


    //  ***  compute relative change in x by current step  *** 

  L90:
    reldx1 = reldst(p, &d_[1], &x[1], &x0[1]);

    //  ***  restore x and step if necessary  *** 

    if (goodx) {
      goto L105;
    }
    i_1 = *p;
    for (i_ = 1; i_ <= i_1; ++i_) {
      step[i_] = stlstg[i_];
      x[i_] = x0[i_] + stlstg[i_];
    }

  L105:
    v[fdif] = v[f0] - v[f];
    if (v[fdif] > v[tuner2] * v[preduc]) {
      goto L120;
    }

    //        ***  no (or only a trivial) function decrease 
    //        ***  -- so try new model or smaller radius 

    v[reldx] = reldx1;
    if (v[f] < v[f0]) {
      goto L110;
    }
    iv[mlstgd] = iv[model];
    v[flstgd] = v[f];
    v[f] = v[f0];
    vcopy(*p, &x[1], &x0[1]);
    iv[restor] = 1;
    goto L115;
  L110:
    iv[nfgcal] = nfc;
  L115:
    iv[irc] = 1;
    if (iv[stage] < iv[stglim]) {
      goto L130;
    }
    iv[irc] = 5;
    --iv[radinc];
    goto L130;

    //  ***  nontrivial function decrease achieved  *** 

  L120:
    iv[nfgcal] = nfc;
    rfac1 = one;
    if (goodx) {
      v[reldx] = reldx1;
    }
    v[dstsav] = v[dstnrm];
    if (v[fdif] > v[preduc] * v[tuner1]) {
      goto L200;
    }

    //  ***  decrease was much less than predicted -- either change models 
    //  ***  or accept step with decreased radius. 

    if (iv[stage] >= iv[stglim]) {
      goto L125;
    }
    //        ***  consider switching models  *** 
    iv[irc] = 2;
    goto L130;

    //     ***  accept step with decreased radius  *** 

  L125:
    iv[irc] = 4;

    //  ***  set v(radfac) to fletcher*s decrease factor  *** 

  L130:
    iv[xirc] = iv[irc];
    emax = v[gtstep] + v[fdif];
    v[radfac] = half * rfac1;

    if (emax < v[gtstep]) {
      // Computing MAX 
      d_1 = v[rdfcmn], d_2 = half * v[gtstep] / emax;
      v[radfac] = rfac1 * std::max(d_1, d_2);
    }

    //  ***  do false convergence test  *** 

  L140:
    if (v[reldx] <= v[xftol]) {
      goto L160;
    }
    iv[irc] = iv[xirc];
    if (v[f] < v[f0]) {
      goto L230;
    }
    goto L300;

  L160:
    iv[irc] = 12;
    goto L310;

    //  ***  handle good function decrease  *** 

  L200:
    if (v[fdif] < -v[tuner3] * v[gtstep]) {
      goto L260;
    }

    //     ***  increasing radius looks worthwhile.  see if we just 
    //     ***  recomputed step with a decreased radius or restored step 
    //     ***  after recomputing it with a larger radius. 

    if (iv[radinc] < 0) {
      goto L260;
    }
    if (iv[restor] == 1) {
      goto L260;
    }

    //        ***  we did not.  try a longer step unless this was a newton 
    //        ***  step. 

    v[radfac] = v[rdfcmx];
    gts = v[gtstep];

    if (v[fdif] < (half / v[radfac] - one) * gts) {
      // Computing MAX 
      d_1 = v[incfac], d_2 = half * gts / (gts + v[fdif]);
      v[radfac] = std::max(d_1, d_2);
    }
    iv[irc] = 4;
    if (v[stppar] == zero) {
      goto L300;
    }
    //             ***  step was not a newton step.  recompute it with 
    //             ***  a larger radius. 
    iv[irc] = 5;
    ++iv[radinc];

    //  ***  save values corresponding to good step  *** 

  L230:
    v[flstgd] = v[f];
    iv[mlstgd] = iv[model];
    vcopy(*p, &stlstg[1], &step[1]);
    v[dstsav] = v[dstnrm];
    iv[nfgcal] = nfc;
    v[plstgd] = v[preduc];
    v[gtslst] = v[gtstep];
    goto L300;

    //  ***  accept step with radius unchanged  *** 

  L260:
    v[radfac] = one;
    iv[irc] = 3;
    goto L300;

    //  ***  come here for a restart after convergence  *** 

  L290:
    iv[irc] = iv[xirc];
    if (v[dstsav] >= zero) {
      goto L310;
    }
    iv[irc] = 12;
    goto L310;

    //  ***  perform convergence tests  *** 

  L300:
    iv[xirc] = iv[irc];
  L310:
    //
    //     if ((d_1 = v[f], abs(d_1)) < v[afctol]) {
    //
    if (fabs(v[f]) < v[afctol]) {
      iv[irc] = 10;
    }
    if (half * v[fdif] > v[preduc]) {
      goto L999;
    }
    //
    //     emax = v[rfctol] * (d_1 = v[f0], abs(d_1));
    //
    emax = v[rfctol] * fabs(v[f0]);
    if (v[dstnrm] > v[lmax0] && v[preduc] <= emax) {
      iv[irc] = 11;
    }
    if (v[dst0] < zero) {
      goto L320;
    }
    i_ = 0;

    if ((v[nreduc] > zero && v[nreduc] <= emax) || (v[nreduc] == zero && v[preduc] == zero)) {
      i_ = 2;
    }

    if (v[stppar] == zero && v[reldx] <= v[xctol] && goodx) {
      ++i_;
    }
    if (i_ > 0) {
      iv[irc] = i_ + 6;
    }

    //  ***  consider recomputing step of length v(lmax0) for singular 
    //  ***  convergence test. 

  L320:
    //
    //    if ((i_1 = iv[irc] - 3, abs(i_1)) > 2 && iv[irc] != 12) {
    //
    if (abs(iv[irc] - 3) > 2 && iv[irc] != 12) {
      goto L999;
    }
    if (v[dstnrm] > v[lmax0]) {
      goto L330;
    }
    if (v[preduc] >= emax) {
      goto L999;
    }
    if (v[dst0] <= zero) {
      goto L340;
    }
    if (half * v[dst0] <= v[lmax0]) {
      goto L999;
    }
    goto L340;
  L330:
    if (half * v[dstnrm] <= v[lmax0]) {
      goto L999;
    }
    xmax = v[lmax0] / v[dstnrm];
    if (xmax * (two - xmax) * v[preduc] >= emax) {
      goto L999;
    }
  L340:
    if (v[nreduc] < zero) {
      goto L370;
    }

    //  ***  recompute v(preduc) for use in singular convergence test  *** 

    v[gtslst] = v[gtstep];
    v[dstsav] = v[dstnrm];
    if (iv[irc] == 12) {
      v[dstsav] = -v[dstsav];
    }
    v[plstgd] = v[preduc];
    iv[irc] = 6;
    vcopy(*p, &stlstg[1], &step[1]);
    goto L999;

    //  ***  perform singular convergence test with recomputed v(preduc)  *** 

  L360:
    v[gtstep] = v[gtslst];
    //
    //    v[dstnrm] = (d_1 = v[dstsav], abs(d_1));
    //
    v[dstnrm] = fabs(v[dstsav]);
    vcopy(*p, &step[1], &stlstg[1]);
    iv[irc] = iv[xirc];
    if (v[dstsav] <= zero) {
      iv[irc] = 12;
    }
    v[nreduc] = -v[preduc];
    v[preduc] = v[plstgd];
  L370:
    //
    //    if (-v[nreduc] <= v[rfctol] * (d_1 = v[f0], abs(d_1))) {
    //
    if (-v[nreduc] <= v[rfctol] * fabs(v[f0])) {
      iv[irc] = 11;
    }

  L999:
    return 0;

  } /* assess */

  double dotprd(int *p, double *x, double *y) {
    /* Initialized data */

    const double one = 1.;
    double sqteta = 0.;
    const double zero = 0.;

    /* System generated locals */
    int i_1;
    double ret_val;

    /* Local variables */
    int i_;
    double t;

    //  ***  return the inner product of the p-vectors x and y.  *** 

    //  ***  rmdcon(2) returns a machine-dependent constant, sqteta, which 
    //  ***  is slightly larger than the smallest positive number that 
    //  ***  can be squared without underflowing. 

    /* Parameter adjustments */
    --y;
    --x;

    /* Function Body */

    ret_val = zero;
    if (*p <= 0) {
      goto L999;
    }
    if (sqteta == zero) {
      sqteta = rmdcon(2);
    }
    i_1 = *p;
    for (i_ = 1; i_ <= i_1; ++i_) {
      // Computing MAX 
      //
      //        d_3 = (d_1 = x[i_], abs(d_1)), d_4 = (d_2 = y[i_], abs(d_2));
      //        t = std::max(d_3,d_4);

      t = std::max(fabs(x[i_]), fabs(y[i_]));
      if (t > one) {
        goto L10;
      }
      if (t < sqteta) {
        goto L20;
      }
      t = x[i_] / sqteta * y[i_];
      if (fabs(t) < sqteta) {
        goto L20;
      }
    L10:
      ret_val += x[i_] * y[i_];
    L20:
      ;
    }

  L999:
    return ret_val;

  } /* dotprd */

  int linvrt(int *n, double *lin, double *l) {
    /* Initialized data */

    const double one = 1.;
    const double zero = 0.;

    /* System generated locals */
    int i_1, i_2, i_3;

    /* Local variables */
    int i_, k;
    double t;
    int j_0, j_1, k0, ii, jj, im1, np1;


    //  ***  compute  lin = l**-1,  both  n x n  lower triang. stored   *** 
    //  ***  compactly by rows.  lin and l may share the same storage.  *** 

    //  ***  parameters  *** 

    //     dimension l(n*(n+1)/2), lin(n*(n+1)/2) 

    //  ***  local variables  *** 


    /* Parameter adjustments */
    --l;
    --lin;

    /* Function Body */

    //  ***  body  *** 

    np1 = *n + 1;
    j_0 = *n * np1 / 2;
    i_1 = *n;
    for (ii = 1; ii <= i_1; ++ii) {
      i_ = np1 - ii;
      lin[j_0] = one / l[j_0];
      if (i_ <= 1) {
        goto L999;
      }
      j_1 = j_0;
      im1 = i_ - 1;
      i_2 = im1;
      for (jj = 1; jj <= i_2; ++jj) {
        t = zero;
        j_0 = j_1;
        k0 = j_1 - jj;
        i_3 = jj;
        for (k = 1; k <= i_3; ++k) {
          t -= l[k0] * lin[j_0];
          --j_0;
          k0 = k0 + k - i_;
        }
        lin[j_0] = t / l[k0];
      }
      --j_0;
    }
  L999:
    return 0;
  } /* linvrt */

  int lsqrt(int n1, int *n, double *l, double *a, int *irc) {
    /* Initialized data */

    const double zero = 0.;

    /* System generated locals */
    int i_1, i_2, i_3;

    /* Local variables */
    int i_, j, k;
    double t;
    int i0, j_0, ij, ik, jk;
    double td;
    int im1, jm1;


    //  ***  compute rows n1 through n of the cholesky factor  l  of 
    //  ***  a = l*(l**t),  where  l  and the lower triangle of  a  are both 
    //  ***  stored compactly by rows (and may occupy the same storage). 
    //  ***  irc = 0 means all went well.  irc = j means the leading 
    //  ***  principal  j x j  submatrix of  a  is not positive definite -- 
    //  ***  and  l(j*(j+1)/2)  contains the (nonpos.) reduced j-th diagonal. 

    //  ***  parameters  *** 

    //     dimension l(n*(n+1)/2), a(n*(n+1)/2) 

    //  ***  local variables  *** 


    //  ***  intrinsic functions  *** 

    /* Parameter adjustments */
    --a;
    --l;

    /* Function Body */

    //  ***  body  *** 

    i0 = n1 * (n1 - 1) / 2;
    i_1 = *n;
    for (i_ = n1; i_ <= i_1; ++i_) {
      td = zero;
      if (i_ == 1) {
        goto L40;
      }
      j_0 = 0;
      im1 = i_ - 1;
      i_2 = im1;
      for (j = 1; j <= i_2; ++j) {
        t = zero;
        if (j == 1) {
          goto L20;
        }
        jm1 = j - 1;
        i_3 = jm1;
        for (k = 1; k <= i_3; ++k) {
          ik = i0 + k;
          jk = j_0 + k;
          t += l[ik] * l[jk];
        }
      L20:
        ij = i0 + j;
        j_0 += j;
        t = (a[ij] - t) / l[j_0];
        l[ij] = t;
        td += t * t;
      }
    L40:
      i0 += i_;
      t = a[i0] - td;
      if (t <= zero) {
        goto L60;
      }
      l[i0] = sqrt(t);
    }

    *irc = 0;
    goto L999;

  L60:
    l[i0] = t;
    *irc = i_;

  L999:
    return 0;

  } /* lsqrt */

  int ltsqar(int *n, double *a, double *l) {
    /* System generated locals */
    int i_1, i_2, i_3;

    /* Local variables */
    int i_, j, k, m, i1, ii;
    double lj, lii;
    int iim1;


    //  ***  set a to lower triangle of (l**t) * l  *** 

    //  ***  l = n x n lower triang. matrix stored rowwise.  *** 
    //  ***  a is also stored rowwise and may share storage with l.  *** 

    //     dimension a(n*(n+1)/2), l(n*(n+1)/2) 


    /* Parameter adjustments */
    --l;
    --a;

    /* Function Body */
    ii = 0;
    i_1 = *n;
    for (i_ = 1; i_ <= i_1; ++i_) {
      i1 = ii + 1;
      ii += i_;
      m = 1;
      if (i_ == 1) {
        goto L30;
      }
      iim1 = ii - 1;
      i_2 = iim1;
      for (j = i1; j <= i_2; ++j) {
        lj = l[j];
        i_3 = j;
        for (k = i1; k <= i_3; ++k) {
          a[m] += lj * l[k];
          ++m;
        }
      }
    L30:
      lii = l[ii];
      i_2 = ii;
      for (j = i1; j <= i_2; ++j) {
        a[j] = lii * l[j];
      }
    }

    return 0;
  } /* ltsqar */

  int nl2sol_defaults(int *iv, double *v) {
    /* Initialized data */

    const double one = 1.;
    const double three = 3.;
    const int covprt = 14;
    const int covreq = 15;
    const int dtype = 16;
    const int inits = 25;
    const int mxfcal = 17;
    const int mxiter = 18;
    const int outlev = 19;
    const int parprt = 20;
    const int prunit = 21;
    const int solprt = 22;
    const int statpr = 23;
    const int x0prt = 24;
    const int afctol = 31;
    const int cosmin = 43;
    const int decfac = 22;
    const int delta0 = 44;
    const int dfac = 41;
    const int dinit = 38;
    const int dltfdc = 40;
    const int dltfdj = 36;
    const int d0init = 37;
    const int epslon = 19;
    const int fuzz = 45;
    const int incfac = 23;
    const int jtinit = 39;
    const int lmax0 = 35;
    const int phmnfc = 20;
    const int phmxfc = 21;
    const int rdfcmn = 24;
    const int rdfcmx = 25;
    const int rfctol = 32;
    const int rlimit = 42;
    const int tuner1 = 26;
    const int tuner2 = 27;
    const int tuner3 = 28;
    const int tuner4 = 29;
    const int tuner5 = 30;
    const int xctol = 33;
    const int xftol = 34;

    /* System generated locals */
    double d_1, d_2, d_3;

    /* Local variables */
    double machep;
    double mepcrt, sqteps;


    //  ***  supply nl2sol (version 2.2) default values to iv and v  *** 

    /* Parameter adjustments */
    --v;
    --iv;

    /* Function Body */
    // -----------------------------------------------------------------------

    iv[1] = 12;
    iv[covprt] = 1;
    iv[covreq] = 1;
    iv[dtype] = 1;
    iv[inits] = 0;
    iv[mxfcal] = 200;
    iv[mxiter] = 150;
    iv[outlev] = 1;
    iv[parprt] = 1;
    iv[prunit] = 6; // ForTran standard output uniti number, unused.
    iv[solprt] = 1;
    iv[statpr] = 1;
    iv[x0prt] = 1;

    machep = rmdcon(3);
    v[afctol] = 1e-20;
    if (machep > 1e-10) {
      // Computing 2nd power 
      d_1 = machep;
      v[afctol] = d_1 * d_1;
    }
    // Computing MAX 
    d_1 = 1e-6, d_2 = machep * 100.;
    v[cosmin] = std::max(d_1, d_2);
    v[decfac] = .5;
    sqteps = rmdcon(4);
    v[delta0] = sqteps;
    v[dfac] = .6;
    v[dinit] = 0.;
    d_1 = one / three;
    //    mepcrt = pow_dd(&machep, &d_1);
    mepcrt = pow(machep, d_1);
    v[dltfdc] = mepcrt;
    v[dltfdj] = sqteps;
    v[d0init] = 1.;
    v[epslon] = .1;
    v[fuzz] = 1.5;
    v[incfac] = 2.;
    v[jtinit] = 1e-6;
    v[lmax0] = 100.;
    v[phmnfc] = -.1;
    v[phmxfc] = .1;
    v[rdfcmn] = .1;
    v[rdfcmx] = 4.;
    // Computing MAX 
    // Computing 2nd power 
    d_3 = mepcrt;
    d_1 = 1e-10, d_2 = d_3 * d_3;
    v[rfctol] = std::max(d_1, d_2);
    v[rlimit] = rmdcon(5);
    v[tuner1] = .1;
    v[tuner2] = 1e-4;
    v[tuner3] = .75;
    v[tuner4] = .5;
    v[tuner5] = .75;
    v[xctol] = sqteps;
    v[xftol] = machep * 100.;

    return 0;
  } /* nl2sol_defaults */

  int parchk(int *iv, int *n, int *nn, int *p, double *v) {
    /* Initialized data */

    const int nvdflt = 27;
    const double zero = 0.;
    const int dtype = 16;
    const int dtype0 = 29;
    const int d0init = 37;
    const int epslon = 19;
    const int inits = 25;
    const int jtinit = 39;
    const int jtol0 = 86;
    const int jtol1 = 87;
    const int oldn = 45;
    const int oldnn = 46;
    const int oldp = 47;
    const int parprt = 20;
    const int parsv1 = 51;
    const int prunit = 21;
    double big = 0.;
    double tiny = 1.;

    double vm[27] = { .001,-.99,.001,.01,1.2,.01,1.2,0.,0.,.001,
      -1.,0.0,0.0,0.0,0.,0.,0.0,0.0,0.,-10.,0.,0.0,0.,1e10,0.0,0.0,1.01
    };
    double vx[27] = { .9,-.001,10.,.8,100.,.8,100.,.5,.5,1.,1.,0.0,
      0.0,.1,1.,1.,0.0,1.,0.0,0.0,0.0,1.,1.,0.0,1.,1.,100. };

    /* System generated locals */
    int i_1;

    /* Local variables */
    int i_, k, l, m;
    int jtolp;
    double machep, vk;
    int pu;
    int iv1;

    //  ***  check nl2sol (version 2.2) parameters, print changed values  *** 

    //     dimension iv(*), v(*) 

    // dfault -- supplies dfault parameter values. 
    // rmdcon -- returns machine-dependent constants. 
    // vcopy  -- copies one vector to another. 

    //  ***  local variables  *** 

    //     character*4 cngd(3), dflt(3), vn(2,27), which(3) 

    //  ***  iv and v subscripts  *** 

    /* Parameter adjustments */
    --v;
    --iv;

    /* Function Body */

    // .......................................................................

    if (iv[1] == 0) {
      nl2sol_defaults(&iv[1], &v[1]);
    }
    pu = iv[prunit];
    iv1 = iv[1];
    if (iv1 != 12) {
      goto L30;
    }
    if (*nn >= *n && *n >= *p && *p >= 1) {
      goto L20;
    }
    iv[1] = 16;
    goto L999;
  L20:
    k = iv[21];
    nl2sol_defaults(&iv[21], &v[33]);
    iv[21] = k;
    iv[dtype0] = iv[dtype + 20];
    iv[oldn] = *n;
    iv[oldnn] = *nn;
    iv[oldp] = *p;
    goto L80;
  L30:
    if (*n == iv[oldn] && *nn == iv[oldnn] && *p == iv[oldp]) {
      goto L50;
    }
    iv[1] = 17;
    goto L999;

  L50:
    if (iv1 <= 11 && iv1 >= 1) {
      goto L70;
    }
    iv[1] = 50;
    goto L999;

  L70:
  L80:
    if (big > tiny) {
      goto L90;
    }
    tiny = rmdcon(1);
    machep = rmdcon(3);
    big = rmdcon(6);
    vm[11] = machep;
    vx[11] = big;
    vm[12] = tiny;
    vx[12] = big;
    vm[13] = machep;
    vm[16] = tiny;
    vx[16] = big;
    vm[17] = machep;
    vx[18] = big;
    vx[19] = big;
    vx[20] = big;
    vm[21] = machep;
    vx[23] = rmdcon(5);
    vm[24] = machep;
    vm[25] = machep;
  L90:
    m = 0;
    if (iv[inits] >= 0 && iv[inits] <= 2) {
      goto L110;
    }
    m = 18;
  L110:
    k = epslon;
    i_1 = nvdflt;
    for (i_ = 1; i_ <= i_1; ++i_) {
      vk = v[k];
      if (vk >= vm[i_ - 1] && vk <= vx[i_ - 1]) {
        goto L130;
      }
      m = k;
    L130:
      ++k;
    }

    if (iv1 == 12 && v[jtinit] > zero) {
      goto L170;
    }

    //  ***  check jtol values  *** 

    jtolp = jtol0 + *p;
    i_1 = jtolp;
    for (i_ = jtol1; i_ <= i_1; ++i_) {
      if (v[i_] > zero) {
        goto L160;
      }
      k = i_ - jtol0;
      m = i_;
    L160:
      ;
    }

  L170:
    if (m == 0) {
      goto L180;
    }
    iv[1] = m;
    goto L999;

  L180:
    if (pu == 0 || iv[parprt] == 0) {
      goto L999;
    }
    if (iv1 != 12 || iv[inits] == 0) {
      goto L200;
    }
    m = 1;
  L200:
    if (iv[dtype] == iv[dtype0]) {
      goto L210;
    }
    m = 1;
  L210:
    k = epslon;
    l = parsv1;
    i_1 = nvdflt;
    for (i_ = 1; i_ <= i_1; ++i_) {
      if (v[k] == v[l]) {
        goto L230;
      }
      m = 1;
    L230:
      ++k;
      ++l;
    }
    iv[dtype0] = iv[dtype];
    vcopy(nvdflt, &v[parsv1], &v[epslon]);
    if (iv1 != 12) {
      goto L999;
    }
    if (v[jtinit] > zero) {
      goto L260;
    }
    jtolp = jtol0 + *p;
  L260:
    if (v[d0init] > zero) {
      goto L999;
    }
    k = jtol1 + *p;
    l = k + *p - 1;
  L999:
    return 0;
  } /* parchk */

  int dupdat(double *d_, int *iv, double *j, int *n, int *nn, int *p, double *v) {
    /* Initialized data */

    const int dfac = 41;
    const int dtype = 16;
    const int jtol0 = 86;
    const int niter = 31;
    const int s = 53;
    const double zero = 0.;

    /* System generated locals */
    int j_dim1, j_offset, i_1;
    double d_1, d_2;

    /* Local variables */
    int i_;
    double vdfac, t;
    int jtoli, d0, s1;
    double sii;


    //  ***  update scale vector d for nl2itr (nl2sol version 2.2)  *** 

    //  ***  parameter declarations  *** 

    /* Parameter adjustments */
    --iv;
    j_dim1 = *nn;
    j_offset = j_dim1 + 1;
    j -= j_offset;
    --d_;
    --v;

    /* Function Body */
    // -----------------------------------------------------------------------

    i_ = iv[dtype];
    if (i_ == 1) {
      goto L20;
    }
    if (iv[niter] > 0) {
      goto L999;
    }

  L20:
    vdfac = v[dfac];
    d0 = jtol0 + *p;
    s1 = iv[s] - 1;
    i_1 = *p;
    for (i_ = 1; i_ <= i_1; ++i_) {
      s1 += i_;
      sii = v[s1];
      t = v2norm(n, &j[i_ * j_dim1 + 1]);
      if (sii > zero) {
        t = sqrt(t * t + sii);
      }
      jtoli = jtol0 + i_;
      ++d0;
      if (t < v[jtoli]) {
        // Computing MAX 
        d_1 = v[d0], d_2 = v[jtoli];
        t = std::max(d_1, d_2);
      }
      // Computing MAX 
      d_1 = vdfac * d_[i_];
      d_[i_] = std::max(d_1, t);
    }

  L999:
    return 0;
  } /* dupdat */

  double reldst(int *p, double *d_, double *x, double *x0) {
    /* Initialized data */

    const double zero = 0.;

    /* System generated locals */
    int i_1;
    double ret_val;

    /* Local variables */
    double emax, xmax;
    int i_;
    double t;


    //  ***  compute and return relative difference between x and x0  *** 
    //  ***  nl2sol version 2.2  *** 

    /* Parameter adjustments */
    --x0;
    --x;
    --d_;

    /* Function Body */

    emax = zero;
    xmax = zero;
    i_1 = *p;
    for (i_ = 1; i_ <= i_1; ++i_) {
      //
      //        t = (d_1 = d_[i_] * (x[i_] - x0[i_]), abs(d_1));
      //
      t = fabs(d_[i_] * (x[i_] - x0[i_]));
      if (emax < t) {
        emax = t;
      }
      //
      //        t = d_[i_] * ((d_1 = x[i_], abs(d_1)) + (d_2 = x0[i_], abs(d_2)));
      //
      t = d_[i_] * (fabs(x[i_]) + fabs(x0[i_]));
      if (xmax < t) {
        xmax = t;
      }
    }
    ret_val = zero;
    if (xmax > zero) {
      ret_val = emax / xmax;
    }
    return ret_val;
  } /* reldst */

  double rmdcon(int k) {
    /* Initialized data */

    const double one001 = 1.001;
    const double pt999 = .999;
    const double machep = DBL_EPSILON, big = DBL_MAX, eta = DBL_MIN;

    /* System generated locals */
    double ret_val;

    //  ***  return machine dependent constants used by nl2sol  *** 

    // +++  comments below contain data statements for various machines.  +++ 
    // +++  to convert to another machine, place a c in column 1 of the   +++ 
    // +++  data statement line(s) that correspond to the current machine +++ 
    // +++  and remove the c from column 1 of the data statement line(s)  +++ 
    // +++  that correspond to the new machine.                           +++ 

    //  ***  the constant returned depends on k... 

    //  ***        k = 1... smallest pos. eta such that -eta exists. 
    //  ***        k = 2... square root of 1.001*eta. 
    //  ***        k = 3... unit roundoff = smallest pos. no. machep such 
    //  ***                 that 1 + machep .gt. 1 .and. 1 - machep .lt. 1. 
    //  ***        k = 4... square root of 0.999*machep. 
    //  ***        k = 5... square root of 0.999*big (see k = 6). 
    //  ***        k = 6... largest machine no. big such that -big exists. 

    // -------------------------------  body  --------------------------------

    switch (k) {
    case 1:  goto L10;
    case 2:  goto L20;
    case 3:  goto L30;
    case 4:  goto L40;
    case 5:  goto L50;
    case 6:  goto L60;
    }

  L10:
    ret_val = eta;
    goto L999;

  L20:
    ret_val = sqrt(one001 * eta);
    goto L999;

  L30:
    ret_val = machep;
    goto L999;

  L40:
    ret_val = sqrt(pt999 * machep);
    goto L999;

  L50:
    ret_val = sqrt(pt999 * big);
    goto L999;

  L60:
    ret_val = big;

  L999:
    return ret_val;
  } /* rmdcon */

  int rptmul(int func, int *ipivot, double *j, int *nn, int *p, double *rd, double *x, double *y, double *z_) {
    /* System generated locals */
    int j_dim1, j_offset, i_1, i_2;

    /* Local variables */
    int i_, k;
    double zk;
    int im1, km1;

    //  ***  func = 1... set  y = rmat * (perm**t) * x. 
    //  ***  func = 2... set  y = perm * (rmat**t) * rmat * (perm**t) * x. 
    //  ***  func = 3... set  y = perm * (rmat**t) x. 


    //  ***  perm = matrix whose i-th col. is the ipivot(i)-th unit vector. 
    //  ***  rmat is the upper triangular matrix whose strict upper triangle 
    //  ***       is stored in  j  and whose diagonal is stored in rd. 
    //  ***  z is a scratch vector. 
    //  ***  x and y may share storage. 


    //  ***  local variables  *** 


    //  ***  external function  *** 


    // -----------------------------------------------------------------------

    /* Parameter adjustments */
    --z_;
    --y;
    --x;
    --rd;
    j_dim1 = *nn;
    j_offset = j_dim1 + 1;
    j -= j_offset;
    --ipivot;

    /* Function Body */
    if (func > 2) {
      goto L50;
    }

    //  ***  first set  z = (perm**t) * x  *** 

    i_1 = *p;
    for (i_ = 1; i_ <= i_1; ++i_) {
      k = ipivot[i_];
      z_[i_] = x[k];
    }

    //  ***  now set  y = rmat * z  *** 

    y[1] = z_[1] * rd[1];
    if (*p <= 1) {
      goto L40;
    }
    i_1 = *p;
    for (k = 2; k <= i_1; ++k) {
      km1 = k - 1;
      zk = z_[k];
      i_2 = km1;
      for (i_ = 1; i_ <= i_2; ++i_) {
        y[i_] += j[i_ + k * j_dim1] * zk;
      }
      y[k] = zk * rd[k];
    }

  L40:
    if (func <= 1) {
      goto L999;
    }
    goto L70;

  L50:
    i_1 = *p;
    for (i_ = 1; i_ <= i_1; ++i_) {
      y[i_] = x[i_];
    }

    //  ***  set  z = (rmat**t) * y  *** 

  L70:
    z_[1] = y[1] * rd[1];
    if (*p == 1) {
      goto L90;
    }
    i_1 = *p;
    for (i_ = 2; i_ <= i_1; ++i_) {
      im1 = i_ - 1;
      z_[i_] = y[i_] * rd[i_] + dotprd(&im1, &j[i_ * j_dim1 + 1], &y[
        1]);
    }

    //  ***  now set  y = perm * z  *** 

  L90:
    i_1 = *p;
    for (i_ = 1; i_ <= i_1; ++i_) {
      k = ipivot[i_];
      y[k] = z_[i_];
    }

  L999:
    return 0;
  } /* rptmul */

  int slupdt(double *a, double *cosmin, int *p, double *size, double *step, double *u, double *w, double *wchmtd, double *wscale, double *y) {
    /* Initialized data */

    const double half = .5;
    const double one = 1.;
    const double zero = 0.;

    /* System generated locals */
    int i_1, i_2;

    /* Local variables */
    int i_, j, k;
    double t;
    double ui, wi, denmin;
    double sdotwm;

    //  ***  update symmetric  a  so that  a * step = y  *** 
    //  ***  (lower triangle of  a  stored rowwise       *** 

    /* Parameter adjustments */
    --a;
    --y;
    --wchmtd;
    --w;
    --u;
    --step;

    /* Function Body */

    // -----------------------------------------------------------------------

    sdotwm = dotprd(p, &step[1], &wchmtd[1]);
    denmin = *cosmin * v2norm(p, &step[1]) * v2norm(p, &wchmtd[1]);
    *wscale = one;
    if (denmin != zero) {
      // Computing MIN 
      //
      //        d_2 = one, d_3 = (d_1 = sdotwm / denmin, abs(d_1));
      //        *wscale = std::min(d_2,d_3);
      *wscale = std::min(one, fabs(sdotwm / denmin));
    }
    t = zero;
    if (sdotwm != zero) {
      t = *wscale / sdotwm;
    }
    i_1 = *p;
    for (i_ = 1; i_ <= i_1; ++i_) {
      w[i_] = t * wchmtd[i_];
    }
    slvmul(p, &u[1], &a[1], &step[1]);
    t = half * (*size * dotprd(p, &step[1], &u[1]) - dotprd(p, &step[1], &y[
      1]));
    i_1 = *p;
    for (i_ = 1; i_ <= i_1; ++i_) {
      u[i_] = t * w[i_] + y[i_] - *size * u[i_];
    }

    //  ***  set  a = a + u*(w**t) + w*(u**t)  *** 

    k = 1;
    i_1 = *p;
    for (i_ = 1; i_ <= i_1; ++i_) {
      ui = u[i_];
      wi = w[i_];
      i_2 = i_;
      for (j = 1; j <= i_2; ++j) {
        a[k] = *size * a[k] + ui * w[j] + wi * u[j];
        ++k;
      }
    }

    return 0;
  } /* slupdt */

  int slvmul(int *p, double *y, double *s, double *x) {
    /* System generated locals */
    int i_1, i_2;

    /* Local variables */
    int i_, j, k;
    double xi;
    int im1;


    //  ***  set  y = s * x,  s = p x p symmetric matrix.  *** 
    //  ***  lower triangle of  s  stored rowwise.         *** 

    // -----------------------------------------------------------------------

    /* Parameter adjustments */
    --x;
    --y;
    --s;

    /* Function Body */
    j = 1;
    i_1 = *p;
    for (i_ = 1; i_ <= i_1; ++i_) {
      y[i_] = dotprd(&i_, &s[j], &x[1]);
      j += i_;
    }

    if (*p <= 1) {
      goto L999;
    }
    j = 1;
    i_1 = *p;
    for (i_ = 2; i_ <= i_1; ++i_) {
      xi = x[i_];
      im1 = i_ - 1;
      ++j;
      i_2 = im1;
      for (k = 1; k <= i_2; ++k) {
        y[k] += s[j] * xi;
        ++j;
      }
    }

  L999:
    return 0;
  } /* slvmul */

  double v2norm(int *p, double *x) {
    /* Initialized data */

    const double one = 1.;
    const double zero = 0.;
    double sqteta = 0.;

    /* System generated locals */
    int i_1;
    double ret_val;

    /* Local variables */
    int i_, j;
    double r_, t, scale, xi;

    //  ***  return the 2-norm of the p-vector x, taking  *** 
    //  ***  care to avoid the most likely underflows.    *** 

    /* Parameter adjustments */
    --x;

    /* Function Body */

    if (*p > 0) {
      goto L10;
    }
    ret_val = zero;
    goto L999;
  L10:
    i_1 = *p;
    for (i_ = 1; i_ <= i_1; ++i_) {
      if (x[i_] != zero) {
        goto L30;
      }
    }
    ret_val = zero;
    goto L999;

  L30:
    //
    // The original code here was
    //
    //    scale = (d_1 = x[i_], abs(d_1));
    //
    // but this resulted in a call to the integer function of abs().
    scale = fabs(x[i_]);
    if (i_ < *p) {
      goto L40;
    }
    ret_val = scale;
    goto L999;
  L40:
    t = one;
    if (sqteta == zero) {
      sqteta = rmdcon(2);
    }

    //     ***  sqteta is (slightly larger than) the square root of the 
    //     ***  smallest positive floating point number on the machine. 
    //     ***  the tests involving sqteta are done to prevent underflows. 

    j = i_ + 1;
    i_1 = *p;
    for (i_ = j; i_ <= i_1; ++i_) {
      //
      //        xi = (d_1 = x[i_], abs(d_1));
      //
      xi = fabs(x[i_]);
      if (xi > scale) {
        goto L50;
      }
      r_ = xi / scale;
      if (r_ > sqteta) {
        t += r_ * r_;
      }
      goto L60;
    L50:
      r_ = scale / xi;
      if (r_ <= sqteta) {
        r_ = zero;
      }
      t = one + t * r_ * r_;
      scale = xi;
    L60:
      ;
    }

    ret_val = scale * sqrt(t);
  L999:
    return ret_val;
  } /* v2norm */

  int vaxpy(int *p, double *w, double a, double *x, double *y) {
    /* System generated locals */
    int i_1;

    /* Local variables */
    int i_;

    //  ***  set w = a*x + y  --  w, x, y = p-vectors, a = scalar  *** 

    /* Parameter adjustments */
    --y;
    --x;
    --w;

    /* Function Body */
    i_1 = *p;
    for (i_ = 1; i_ <= i_1; ++i_) {
      w[i_] = a * x[i_] + y[i_];
    }
    return 0;
  } /* vaxpy */

  int vcopy(int p, double *y, double *x) {
    /* System generated locals */
    int i_1;

    /* Local variables */
    int i_;


    //  ***  set y = x, where x and y are p-vectors  *** 

    /* Parameter adjustments */
    --x;
    --y;

    /* Function Body */
    i_1 = p;
    for (i_ = 1; i_ <= i_1; ++i_) {
      y[i_] = x[i_];
    }
    return 0;
  } /* vcopy */

  int vscopy(int *p, double *y, double s) {
    /* System generated locals */
    int i_1;

    /* Local variables */
    int i_;

    //  ***  set p-vector y to scalar s  *** 

    /* Parameter adjustments */
    --y;

    /* Function Body */
    i_1 = *p;
    for (i_ = 1; i_ <= i_1; ++i_) {
      y[i_] = s;
    }
    return 0;
  } /* vscopy */

void calcR(NL2SOL::Functor& functor, int *function_evaluations, int * n, int * p, double * x, int *withinbounds, double *residuals) {
  if(!functor(x,*p,residuals,*n)){
    *withinbounds = 0;   //   This indicates out of bounds.
    return;
  }
  ++(*function_evaluations);
}

void calcJ(NL2SOL::Functor& functor, int *jacobian_evaluations, int * n, int * p, double * x, int *withinbounds, double *jacobian) {
  INSIST(functor.have_jacobian());
  if (!functor.jacobian(x,*p,*n,jacobian)) {
    *withinbounds = 0;   //   This indicates out of bounds.
    return;
  }
  ++(*jacobian_evaluations);
}

int nl2sol(NL2SOL::Functor& functor, int* function_evaluations, int* jacobian_evaluations, int* n, int* p, double* x, int* iv, double* v) {
  /* Initialized data */
  const int d_ = 27;
  const int j = 33;
  const int r_ = 50;
  const int nfcall = 6;
  const int nfgcal = 7;
  const int toobig = 2;

  /* System generated locals */
  int i_1;

  /* Local variables */
  int d1, j_1, r1;
  int nf;
  bool strted;

  //  ***  minimize nonlinear sum of squares using analytic jacobian  *** 
  //  ***  (nl2sol version 2.2)  *** 

  //     dimension iv(60+p),  v(93 + n*p + 3*n + p*(3*p+33)/2) 
  //     dimension uiparm(*), urparm(*) 

  //  ***  purpose  *** 

  //     Given a p-vector x of parameters, calcr computes an n-vector 
  //     r = r(x) of residuals corresponding to x.  (r(x) probably arises 
  //     from a nonlinear model involving p parameters and n observations.) 
  //     this routine interacts with nl2itr to seek a parameter vector x 
  //     that minimizes the sum of the squares of (the components of) r(x), 
  //     i.e., that minimizes the sum-of-squares function 
  //     f(x) = (r(x)**t) * r(x) / 2.  r(x) is assumed to be a twice con- 
  //     tinuously differentiable function of x. 

  // --------------------------  parameter usage  --------------------------

  // n........ (input) the number of observations, i.e., the number of 
  //                  components in r(x).  n must be .ge. p. 
  // p........ (input) the number of parameters (components in x).  p must 
  //                  be positive. 
  // x........ (input/output).  on input, x is an initial guess at the 
  //                  desired parameter estimate.  on output, x contains 
  //                  the best parameter estimate found. 
  // calcr.... (input) a subroutine which, given x, computes r(x).  calcr 
  //                  must be declared external in the calling program. 
  //                  it is invoked by 
  //                  when calcr is called, nf is the invocation count 
  //                  for calcr.  it is included for possible use with 
  //                  calcj.  if x is out of bounds (e.g. if it would 
  //                  cause overflow in computing r(x)), then calcr should 
  //                  set nf to 0.  this will cause a shorter step to be 
  //                  attempted.  the other parameters are as described 
  //                  above and below.  calcr should not change n, p, or x. 
  // calcj.... (input) a subroutine which, given x, computes the jacobian 
  //                  matrix j of r at x, i.e., the n by p matrix whose 
  //                  (i,k) entry is the partial derivative of the i-th 
  //                  component of r with respect to x(k).  calcj must be 
  //                  declared external in the calling program.  it is 
  //                  invoked by 
  //                    call calcj(n,p,x,nf,j,uiparm,urparm,ufparm)
  //                  nf is the invocation count for calcr at the time 
  //                  r(x) was evaluated.  the x passed to calcj is 
  //                  usually the one passed to calcr on either its most 
  //                  recent invocation or the one prior to it.  if calcr 
  //                  saves intermediate results for use by calcj, then it 
  //                  is possible to tell from nf whether they are valid 
  //                  for the current x (or which copy is valid if two 
  //                  copies are kept).  if j cannot be computed at x, 
  //                  then calcj should set nf to 0.  in this case, nl2sol 
  //                  will return with iv(1) = 15.  the other parameters 
  //                  to calcj are as described above and below.  calcj 
  //                  should not change n, p, or x. 
  // iv....... (input/output) an int value array of length at least 
  //                  60 + p that helps control the nl2sol algorithm and 
  //                  that is used to store various intermediate quanti- 
  //                  ties.  of particular interest are the initialization/ 
  //                  printing and limit the number of iterations and func- 
  //                  tion evaluations.  see the section on iv input 
  //                  values below. 
  // v........ (input/output) a floating-point value array of length at 
  //                  least 93 + n*p + 3*n + p*(3*p+33)/2 that helps con- 
  //                  trol the nl2sol algorithm and that is used to store 
  //                  various intermediate quantities.  of particular in- 
  //                  terest are the entries in v that limit the length of 
  //                  the first step attempted (lmax0), specify conver- 
  //                  gence tolerances (afctol, rfctol, xctol, xftol), 
  //                  and help choose the step size used in computing the 
  //                  covariance matrix (delta0).  see the section on 
  //                  (selected) v input values below. 
  // uiparm... (input) user int parameter array passed without change 
  //                  to calcr and calcj. 
  // urparm... (input) user floating-point parameter array passed without 
  //                  change to calcr and calcj. 
  // ufparm... (input) user external subroutine or function passed without 
  //                  change to calcr and calcj. 

  //  ***  iv input values (from subroutine dfault)  *** 

  // iv(1)...  on input, iv(1) should have a value between 0 and 12...... 
  //             0 and 12 mean this is a fresh start.  0 means that 
  //             dfault(iv, v) is to be called to provide all default 
  //             values to iv and v.  12 (the value that dfault assigns to 
  //             iv(1)) means the caller has already called dfault(iv, v) 
  //             and has possibly changed some iv and/or v entries to non- 
  //             default values.  default = 12. 
  // iv(covprt)... iv(14) = 1 means print a covariance matrix at the solu- 
  //             tion.  (this matrix is computed just before a return with 
  //             iv(1) = 3, 4, 5, 6.) 
  //             iv(covprt) = 0 means skip this printing.  default = 1. 
  // iv(covreq)... iv(15) = nonzero means compute a covariance matrix 
  //             just before a return with iv(1) = 3, 4, 5, 6.  in 
  //             this case, an approximate covariance matrix is obtained 
  //             in one of several ways.  let k = abs(iv(covreq)) and let 
  //             scale = 2*f(x)/max(1,n-p),  where 2*f(x) is the residual 
  //             sum of squares.  if k = 1 or 2, then a finite-difference 
  //             hessian approximation h is obtained.  if h is positive 
  //             definite (or, for k = 3, if the jacobian matrix j at x 
  //             is nonsingular), then one of the following is computed... 
  //                  k = 1....  scale * h**-1 * (j**t * j) * h**-1. 
  //                  k = 2....  scale * h**-1. 
  //                  k = 3....  scale * (j**t * j)**-1. 
  //             (j**t is the transpose of j, while **-1 means inverse.) 
  //             ient values (calls on calcr and calcj) are used in com- 
  //             puting h (with step sizes determined using v(delta0) -- 
  //             see below), while if iv(covreq) is negative, then only 
  //             function values (calls on calcr) are used (with step 
  //             sizes determined using v(dltfdc) -- see below).  if 
  //             iv(covreq) = 0, then no attempt is made to compute a co- 
  //             variance matrix (unless iv(covprt) = 1, in which case 
  //             iv(covreq) = 1 is assumed).  see iv(covmat) below. 
  //             default = 1. 
  // iv(dtype).... iv(16) tells how the scale vector d (see ref. 1) should 
  //             be chosen.  iv(dtype) .ge. 1 means choose d as described 
  //             below with v(dfac).  iv(dtype) .le. 0 means the caller 
  //             has chosen d and has stored it in v starting at 
  //             v(94 + 2*n + p*(3*p + 31)/2).  default = 1. 
  // iv(inits).... iv(25) tells how the s matrix (see ref. 1) should be 
  //             initialized.  0 means initialize s to 0 (and start with 
  //             the gauss-newton model).  1 and 2 mean that the caller 
  //             has stored the lower triangle of the initial s rowwise in 
  //             v starting at v(87+2*p).  iv(inits) = 1 means start with 
  //             the gauss-newton model, while iv(inits) = 2 means start 
  //             with the augmented model (see ref. 1).  default = 0. 
  // iv(mxfcal)... iv(17) gives the maximum number of function evaluations 
  //             (calls on calcr, excluding those used to compute the co- 
  //             variance matrix) allowed.  if this number does not suf- 
  //             fice, then nl2sol returns with iv(1) = 9.  default = 200. 
  // iv(mxiter)... iv(18) gives the maximum number of iterations allowed. 
  //             it also indirectly limits the number of gradient evalua- 
  //             tions (calls on calcj, excluding those used to compute 
  //             the covariance matrix) to iv(mxiter) + 1.  if iv(mxiter) 
  //             iterations do not suffice, then nl2sol returns with 
  //             iv(1) = 10.  default = 150. 
  // iv(outlev)... iv(19) controls the number and length of iteration sum- 
  //             mary lines printed (by itsmry).  iv(outlev) = 0 means do 
  //             not print any summary lines.  otherwise, print a summary 
  //             line after each abs(iv(outlev)) iterations.  if iv(outlev) 
  //             is positive, then summary lines of length 117 (plus carri- 
  //             age control) are printed, including the following...  the 
  //             iteration and function evaluation counts, current func- 
  //             tion value (v(f) = half the sum of squares), relative 
  //             difference in function values achieved by the latest step 
  //             (i.e., reldf = (f0-v(f))/f0, where f0 is the function 
  //             value from the previous iteration), the relative function 
  //             reduction predicted for the step just taken (i.e., 
  //             preldf = v(preduc) / f0, where v(preduc) is described 
  //             below), the scaled relative change in x (see v(reldx) 
  //             below), the models used in the current iteration (g = 
  //             gauss-newton, s=augmented), the marquardt parameter 
  //             stppar used in computing the last step, the sizing factor 
  //             used in updating s, the 2-norm of the scale vector d 
  //             times the step just taken (see ref. 1), and npreldf, i.e., 
  //             v(nreduc)/f0, where v(nreduc) is described below -- if 
  //             npreldf is positive, then it is the relative function 
  //             reduction predicted for a newton step (one with 
  //             stppar = 0).  if npreldf is zero, either the gradient 
  //             vanishes (as does preldf) or else the augmented model 
  //             is being used and its hessian is indefinite (with preldf 
  //             positive).  if npreldf is negative, then it is the nega- 
  //             of the relative function reduction predicted for a step 
  //             computed with step bound v(lmax0) for use in testing for 
  //             singular convergence. 
  //             length 79 (or 55 is iv(covprt) = 0) are printed, includ- 
  //             ing only the first 6 items listed above (through reldx). 
  //             default = 1. 
  // iv(parprt)... iv(20) = 1 means print any nondefault v values on a 
  //             fresh start or any changed v values on a restart. 
  //             iv(parprt) = 0 means skip this printing.  default = 1. 
  // iv(prunit)... iv(21) is the output unit number on which all printing 
  //             is done.  iv(prunit) = 0 means suppress all printing. 
  //             (setting iv(prunit) to 0 is the only way to suppress the 
  //             one-line termination reason message printed by itsmry.) 
  //             default = standard output unit (unit 6 on most systems). 
  // iv(solprt)... iv(22) = 1 means print out the value of x returned (as 
  //             well as the corresponding gradient and scale vector d). 
  //             iv(solprt) = 0 means skip this printing.  default = 1. 
  // iv(statpr)... iv(23) = 1 means print summary statistics upon return- 
  //             ing.  these consist of the function value (half the sum 
  //             of squares) at x, v(reldx) (see below), the number of 
  //             function and gradient evaluations (calls on calcr and 
  //             calcj respectively, excluding any calls used to compute 
  //             the covariance), the relative function reductions predict- 
  //             ed for the last step taken and for a newton step (or per- 
  //             haps a step bounded by v(lmax0) -- see the descriptions 
  //             of preldf and npreldf under iv(outlev) above), and (if an 
  //             attempt was made to compute the covariance) the number of 
  //             calls on calcr and calcj used in trying to compute the 
  //             covariance.  iv(statpr) = 0 means skip this printing. 
  //             default = 1. 
  // iv(x0prt).... iv(24) = 1 means print the initial x and scale vector d 
  //             (on a fresh start only).  iv(x0prt) = 0 means skip this 
  //             printing.  default = 1. 

  //  ***  (selected) iv output values  *** 

  // iv(1)........ on output, iv(1) is a return code.... 
  //             3 = x-convergence.  the scaled relative difference be- 
  //                  tween the current parameter vector x and a locally 
  //                  optimal parameter vector is very likely at most 
  //                  v(xctol). 
  //             4 = relative function convergence.  the relative differ- 
  //                  ence between the current function value and its lo- 
  //                  cally optimal value is very likely at most v(rfctol). 
  //             5 = both x- and relative function convergence (i.e., the 
  //                  conditions for iv(1) = 3 and iv(1) = 4 both hold). 
  //             6 = absolute function convergence.  the current function 
  //                  value is at most v(afctol) in absolute value. 
  //             7 = singular convergence.  the hessian near the current 
  //                  iterate appears to be singular or nearly so, and a 
  //                  step of length at most v(lmax0) is unlikely to yield 
  //                  a relative function decrease of more than v(rfctol). 
  //             8 = false convergence.  the iterates appear to be converg- 
  //                  ing to a noncritical point.  this may mean that the 
  //                  convergence tolerances (v(afctol), v(rfctol), 
  //                  v(xctol)) are too small for the accuracy to which 
  //                  the function and gradient are being computed, that 
  //                  there is an error in computing the gradient, or that 
  //                  the function or gradient is discontinuous near x. 
  //             9 = function evaluation limit reached without other con- 
  //                  vergence (see iv(mxfcal)). 
  //            10 = iteration limit reached without other convergence 
  //                  (see iv(mxiter)). 
  //            11 = stopx returned .true. (external interrupt).  see the 
  //                  usage notes below. 
  //            13 = f(x) cannot be computed at the initial x. 
  //            14 = bad parameters passed to assess (which should not 
  //                  occur). 
  //            15 = the jacobian could not be computed at x (see calcj 
  //                  above). 
  //            16 = n or p (or parameter nn to nl2itr) out of range -- 
  //                  p .le. 0 or n .lt. p or nn .lt. n. 
  //            17 = restart attempted with n or p (or par. nn to nl2itr) 
  //                  changed. 
  //            18 = iv(inits) is out of range. 
  //            19...45 = v(iv(1)) is out of range. 
  //            50 = iv(1) was out of range. 
  //            87...(86+p) = jtol(iv(1)-86) (i.e., v(iv(1)) is not 
  //                  positive (see v(dfac) below). 
  // iv(covmat)... iv(26) tells whether a covariance matrix was computed. 
  //             the covariance matrix is stored rowwise in v starting at 
  //             v(iv(covmat)).  if iv(covmat) = 0, then no attempt was 
  //             made to compute the covariance.  if iv(covmat) = -1, 
  //             then the finite-difference hessian was indefinite.  and 
  //             and if iv(covmat) = -2, then a successful finite-differ- 
  //             encing step could not be found for some component of x 
  //             (i.e., calcr set nf to 0 for each of two trial steps). 
  //             note that iv(covmat) is reset to 0 after each successful 
  //             step, so if such a step is taken after a restart, then 
  //             the covariance matrix will be recomputed. 
  // iv(d)........ iv(27) is the starting subscript in v of the current 
  //             scale vector d. 
  // iv(g)........ iv(28) is the starting subscript in v of the current 
  //             least-squares gradient vector (j**t)*r. 
  // iv(nfcall)... iv(6) is the number of calls so far made on calcr (i.e., 
  //             function evaluations, including those used in computing 
  //             the covariance). 
  // iv(nfcov).... iv(40) is the number of calls made on calcr when 
  //             trying to compute covariance matrices. 
  // iv(ngcall)... iv(30) is the number of gradient evaluations (calls on 
  //             calcj) so far done (including those used for computing 
  //             the covariance). 
  // iv(ngcov).... iv(41) is the number of calls made on calcj when 
  //             trying to compute covariance matrices. 
  // iv(niter).... iv(31) is the number of iterations performed. 
  // iv(r)........ iv(50) is the starting subscript in v of the residual 
  //             vector r corresponding to x. 

  //  ***  (selected) v input values (from subroutine dfault)  *** 

  // v(afctol)... v(31) is the absolute function convergence tolerance. 
  //             the sum of squares) is less than v(afctol), and if nl2sol 
  //             does not return with iv(1) = 3, 4, or 5, then it returns 
  //             with iv(1) = 6.  default = std::max(10**-20, machep**2), where 
  //             machep is the unit roundoff. 
  // v(delta0)... v(44) is a factor used in choosing the finite-difference 
  //             step size used in computing the covariance matrix when 
  //             iv(covreq) = 1 or 2.  for component i, step size 
  //                  v(delta0) * std::max(abs(x(i)), 1/d(i)) * sign(x(i)) 
  //             is used, where d is the current scale vector (see ref. 1). 
  //             (if this step results in calcr setting nf to 0, then -0.5 
  //             times this step is also tried.)  default = machep**0.5, 
  //             where machep is the unit roundoff. 
  // v(dfac)..... v(41) and the d0 and jtol arrays (see v(d0init) and 
  //             v(jtinit)) are used in updating the scale vector d when 
  //             iv(dtype) .gt. 0.  (d is initialized according to 
  //             v(dinit).)  let d1(i) = 
  //               std::max(sqrt(jcnorm(i)**2 + std::max(s(i,i),0)), v(dfac)*d(i)), 
  //             where jcnorm(i) is the 2-norm of the i-th column of the 
  //             current jacobian matrix and s is the s matrix of ref. 1. 
  //             d1(i) .lt. jtol(i), in which case d(i) is set to 
  //                                std::max(d0(i), jtol(i)). 
  //             iteration as for iv(dtype) = 1 (after any initialization 
  //             due to v(dinit)) and is left unchanged thereafter. 
  //             default = 0.6. 
  // v(dinit).... v(38), if nonnegative, is the value to which the scale 
  //             vector d is initialized.  default = 0. 
  // v(dltfdc)... v(40) helps choose the step size used when computing the 
  //             covariance matrix when iv(covreq) = -1 or -2.  for 
  //             differences involving x(i), the step size first tried is 
  //                       v(dltfdc) * std::max(abs(x(i)), 1/d(i)), 
  //             where d is the current scale vector (see ref. 1).  (if 
  //             this step is too big the first time it is tried, i.e., if 
  //             calcr sets nf to 0, then -0.5 times this step is also 
  //             tried.)  default = machep**(1/3), where machep is the 
  //             unit roundoff. 
  // v(d0init)... v(37), if positive, is the value to which all components 
  //             of the d0 vector (see v(dfac)) are initialized.  if 
  //             v(dfac) = 0, then it is assumed that the caller has 
  //             stored d0 in v starting at v(p+87).  default = 1.0. 
  // v(jtinit)... v(39), if positive, is the value to which all components 
  //             of the jtol array (see v(dfac)) are initialized.  if 
  //             v(jtinit) = 0, then it is assumed that the caller has 
  //             stored jtol in v starting at v(87).  default = 10**-6. 
  // v(lmax0).... v(35) gives the maximum 2-norm allowed for d times the 
  //             very first step that nl2sol attempts.  it is also used 
  //             in testing for singular convergence -- if the function 
  //             reduction predicted for a step of length bounded by 
  //             v(lmax0) is at most v(rfctol) * abs(f0), where  f0  is 
  //             the function value at the start of the current iteration, 
  //             and if nl2sol does not return with iv(1) = 3, 4, 5, or 6, 
  //             then it returns with iv(1) = 7.    default = 100. 
  // v(rfctol)... v(32) is the relative function convergence tolerance. 
  //             reduction (see v(nreduc)) of at most v(rfctol)*abs(f0) at 
  //             the start of the current iteration, where  f0  is the 
  //             then current function value, and if the last step attempt- 
  //             ed achieved no more than twice the predicted function 
  //             decrease, then nl2sol returns with iv(1) = 4 (or 5). 
  //             default = std::max(10**-10, machep**(2/3)), where machep is 
  //             the unit roundoff. 
  // v(tuner1)... v(26) helps decide when to check for false convergence 
  //             and to consider switching models.  this is done if the 
  //             actual function decrease from the current step is no more 
  //             than v(tuner1) times its predicted value.  default = 0.1. 
  // v(xctol).... v(33) is the x-convergence tolerance.  if a newton step 
  //             (see v(nreduc)) is tried that has v(reldx) .le. v(xctol) 
  //             and if this step yields at most twice the predicted func- 
  //             tion decrease, then nl2sol returns with iv(1) = 3 (or 5). 
  //             (see the description of v(reldx) below.) 
  //             default = machep**0.5, where machep is the unit roundoff. 
  // v(xftol).... v(34) is the false convergence tolerance.  if a step is 
  //             tried that gives no more than v(tuner1) times the predict- 
  //             ed function decrease and that has v(reldx) .lt. v(xftol), 
  //             and if nl2sol does not return with iv(1) = 3, 4, 5, 6, or 
  //             7, then it returns with iv(1) = 8.  (see the description 
  //             of v(reldx) below.)  default = 100*machep, where 
  //             machep is the unit roundoff. 
  // v(*)........ dfault supplies to v a number of tuning constants, with 
  //             which it should ordinarily be unnecessary to tinker.  see 
  //             version 2.2 of the nl2sol usage summary (which is an 
  //             appendix to ref. 1). 

  //  ***  (selected) v output values  *** 

  // v(dgnorm)... v(1) is the 2-norm of (d**-1)*g, where g is the most re- 
  //             cently computed gradient and d is the corresponding scale 
  //             vector. 
  // v(dstnrm)... v(2) is the 2-norm of d*step, where step is the most re- 
  //             cently computed step and d is the current scale vector. 
  // v(f)........ v(10) is the current function value (half the sum of 
  //             squares).
  // v(f0)....... v(13) is the function value at the start of the current 
  //             iteration. 
  // v(nreduc)... v(6), if positive, is the maximum function reduction 
  //             possible according to the current model, i.e., the func- 
  //             tion reduction predicted for a newton step (i.e., 
  //             step = -h**-1 * g,  where  g = (j**t) * r  is the current 
  //             gradient and h is the current hessian approximation -- 
  //             h = (j**t)*j  for the gauss-newton model and 
  //             h = (j**t)*j + s  for the augmented model). 
  //                  v(nreduc) = zero means h is not positive definite. 
  //             the function reduction predicted for a step computed with 
  //             a step bound of v(lmax0) for use in testing for singular 
  //             convergence. 
  // v(preduc)... v(7) is the function reduction predicted (by the current 
  //             quadratic model) for the current step.  this (divided by 
  //             v(f0)) is used in testing for relative function 
  //             convergence. 
  // v(reldx).... v(17) is the scaled relative change in x caused by the 
  //             current step, computed as 
  //                  std::max(abs(d(i)*(x(i)-x0(i)), 1 .le. i .le. p) / 
  //                     std::max(d(i)*(abs(x(i))+abs(x0(i))), 1 .le. i .le. p), 
  //             where x = x0 + step. 

  // -------------------------------  notes  -------------------------------

  //  ***  algorithm notes  *** 

  //        see ref. 1 for a description of the algorithm used. 
  //        on problems which are naturally well scaled, better perform- 
  //     ance may be obtained by setting v(d0init) = 1.0 and iv(dtype) = 0, 
  //     which will cause the scale vector d to be set to all ones. 

  //  ***  usage notes  *** 

  //        after a return with iv(1) .le. 11, it is possible to restart, 
  //     i.e., to change some of the iv and v input values described above 
  //     and continue the algorithm from the point where it was interrupt- 
  //     ed.  iv(1) should not be changed, nor should any entries of iv 
  //     and v other than the input values (those supplied by dfault). 
  //        those who do not wish to write a calcj which computes the ja- 
  //     cobian matrix analytically should call nl2sno rather than nl2sol. 
  //     nl2sno uses finite differences to compute an approximate jacobian. 
  //        those who would prefer to provide r and j (the residual and 
  //     jacobian) by reverse communication rather than by writing subrou- 
  //     tines calcr and calcj may call on nl2itr directly.  see the com- 
  //     ments at the beginning of nl2itr. 
  //        those who use nl2sol interactively may wish to supply their 
  //     own stopx function, which should return .true. if the break key 
  //     has been pressed since stopx was last invoked.  this makes it pos- 
  //     sible to externally interrupt nl2sol (which will return with 
  //     iv(1) = 11 if stopx returns .true.). 
  //        storage for j is allocated at the end of v.  thus the caller 
  //     may make v longer than specified above and may allow calcj to use 
  //     elements of j beyond the first n*p as scratch storage. 

  //  ***  portability notes  *** 

  //        the nl2sol distribution tape contains both single- and double- 
  //     precision versions of the nl2sol source code, so it should be un- 
  //     necessary to change precisions. 
  //        only the functions imdcon and rmdcon contain machine-dependent 
  //     constants.  to change from one machine to another, it should 
  //     suffice to change the (few) relevant lines in these functions. 
  //        intrinsic functions are explicitly declared.  on certain com- 
  //     puters (e.g. univac), it may be necessary to comment out these 
  //     declarations.  so that this may be done automatically by a simple 
  //     program, such declarations are preceded by a comment having c/+ 
  //     in columns 1-3 and blanks in columns 4-72 and are followed by 
  //     a comment having c/ in columns 1 and 2 and blanks in columns 3-72. 
  //        the nl2sol source code is expressed in 1966 ansi standard 
  //     fortran.  it may be converted to fortran 77 by 
  //     commenting out all lines that fall between a line having c/6 in 
  //     columns 1-3 and a line having c/7 in columns 1-3 and by removing 
  //     (i.e., replacing by a blank) the c in column 1 of the lines that 
  //     follow the c/7 line and preceed a line having c/ in columns 1-2 
  //     and blanks in columns 3-72.  these changes convert some data 
  //     statements into parameter statements, convert some variables from 
  //     these variables use character strings delimited by primes instead 
  //     of hollerith constants.  (such variables and data statements 
  //     appear only in modules itsmry and parchk.  parameter statements 
  //     appear nearly everywhere.) 

  //  ***  references  *** 

  // 1.  dennis, j.e., gay, d.m., and welsch, r.e. (1981), an adaptive 
  //             nonlinear least-squares algorithm, acm trans. math. 
  //             software, vol. 7, no. 3. 

  //  ***  general  *** 

  //     coded by david m. gay (winter 1979 - winter 1980). 
  //     this subroutine was written in connection with research 
  //     supported by the national science foundation under grants 
  //     mcs-7600324, dcr75-10143, 76-1431ss, mcs76-11989, and 
  //     mcs-7906671. 

  // ----------------------------  declarations  ---------------------------

  // itsmry... prints iteration summary and info about initial and final x. 
  // nl2itr... reverse-communication routine that carries out nl2sol algorithm. 

  /* Parameter adjustments */
  --x;
  --iv;
  --v;

  /* Function Body */

  // +++++++++++++++++++++++++++++++  body  ++++++++++++++++++++++++++++++++

  d1 = (*n << 1) + 94 + *p * (*p * 3 + 31) / 2;
  iv[d_] = d1;
  r1 = d1 + *p;
  iv[r_] = r1;
  j_1 = r1 + *n;
  iv[j] = j_1;
  strted = true;
  if (iv[1] != 0 && iv[1] != 12) {
    goto L40;
  }
  strted = false;
  iv[nfcall] = 1;
  iv[nfgcal] = 1;

L10:
  nf = iv[nfcall];
  calcR(functor, function_evaluations, n, p, &x[1], &nf, &v[r1]);
  if (functor.stop()) {
    iv[1] = 11;
    goto L900;
  }
  if (strted) {
    goto L20;
  }
  if (nf > 0) {
    goto L30;
  }
  iv[1] = 13;
  goto L60;

L20:
  if (nf <= 0) {
    iv[toobig] = 1;
  }
  goto L40;

L30:
  calcJ(functor, jacobian_evaluations, n, p, &x[1], &iv[nfgcal], &v[j_1]);
  if (iv[nfgcal] == 0) {
    //
    // PJ, 2024-02-13:
    //
    // Execution can reach this point for well-specified problems when v[xftol], the false convergence threshold, is not
    // greater than DBL_EPSILON since the assessment of false convergence in assess() does effectively not allow for
    // that. However, since we often want to suppress spurious false convergence assessments, we typically have
    // v[xftol] == 0.0 which leads us here. This is because nl2itr/assess rely on false convergence to arise *before*
    // the situation that we are right on the edge of a constraint boundary, and thus will attempt to take a step and
    // directly compute the Jacobian after this step without first checking that the new guess is actually a valid
    // coordinate point.
    //
    {
      int ok = 1;
      calcR(functor, function_evaluations, n, p, &x[1], &ok, &v[r1]);
      if (ok) { // If we cannot compute the Jacobian but *can* compute the objective, then something is ill specified.
        goto L50;
      }
    }
    //
    // Restore the stored last good point and exit with FALSE_CONVERGENCE.
    //
    const int x0 = 60, rsave = 52;
    vcopy(*p, &x[1], &v[iv[x0]]);
    vcopy(*n, &v[r1], &v[iv[rsave]]);
    iv[1] = 8; // FALSE_CONVERGENCE
#ifdef _DEBUG
    {
      int ok = 1;
      calcR(functor, function_evaluations, n, p, &x[1], &ok, &v[r1]);
      ASSERT(ok);
    }
#endif
    goto L999;
  }
  strted = true;

L40:
  nl2itr(functor, &v[d1], &iv[1], &v[j_1], n, n, p, &v[r1], &v[1], &x[1]);
  if (functor.stop()) {
    iv[1] = 11;
    goto L900;
  }
  if ((i_1 = iv[1] - 2) < 0) {
    goto L10;
  } else if (i_1 == 0) {
    goto L30;
  } else {
    goto L999;
  }

L50:
  iv[1] = 15;
L60:
L900:
L999:
  return 0;
} /* nl2sol */

int nl2sno(NL2SOL::Functor& functor, int* function_evaluations, int* jacobian_evaluations, int* n, int* p, double* x, int* iv, double* v) {

  /* Initialized data */

  const double hfac = 1e3;
  const double negpt5 = -.5;
  const double one = 1.;
  const double zero = 0.;
  const int covprt = 14;
  const int covreq = 15;
  const int d_ = 27;
  const int dtype = 16;
  const int j = 33;
  const int nfcall = 6;
  const int nfgcal = 7;
  const int r_ = 50;
  const int toobig = 2;
  const int dltfdj = 36;
  const int dinit = 38;
  double hlim = 0.;

  /* System generated locals */
  int i_1, i_2;
  double d_1, d_2;

  /* Local variables */
  double h_;
  int i_, k, d1, j_1, r1;
  int dk, nf, rn;
  double xk;
  bool strted;
  int j1k;


  //  ***  like nl2sol, but without calcj -- minimize nonlinear sum of  *** 
  //  ***  squares using finite-difference jacobian approximations      *** 
  //  ***  (nl2sol version 2.2)  *** 

  //     dimension iv(60+p),  v(93 + n*p + 3*n + p*(3*p+33)/2) 

  // -----------------------------  discussion  ----------------------------

  //        the parameters for nl2sno are the same as those for nl2sol 
  //     (which see), except that calcj is omitted.  instead of calling 
  //     calcj to obtain the jacobian matrix of r at x, nl2sno computes 
  //     an approximation to it by finite (forward) differences -- see 
  //     v(dltfdj) below.  nl2sno uses function values only when comput- 
  //     the covariance matrix (rather than the functions and gradients 
  //     that nl2sol may use).  to do so, nl2sno sets iv(covreq) to -1 if 
  //     iv(covprt) = 1 with iv(covreq) = 0 and to minus its absolute 
  //     value otherwise.  thus v(delta0) is never referenced and only 
  //     v(dltfdc) matters -- see nl2sol for a description of v(dltfdc). 
  //        the number of extra calls on calcr used in computing the jaco- 
  //     bian approximation are not included in the function evaluation 
  //     count iv(nfcall) and are not otherwise reported. 

  // v(dltfdj)... v(36) helps choose the step size used when computing the
  //             finite-difference jacobian matrix.  for differences in- 
  //             volving x(i), the step size first tried is 
  //                       v(dltfdj) * std::max(abs(x(i)), 1/d(i)), 
  //             where d is the current scale vector (see ref. 1).  (if 
  //             this step is too big, i.e., if calcr sets nf to 0, then 
  //             smaller steps are tried until the step size is shrunk be- 
  //             low 1000 * machep, where machep is the unit roundoff. 
  //             default = machep**0.5. 

  //  ***  references  *** 

  // 1.  dennis, j.e., gay, d.m., and welsch, r.e. (1981), an adaptive 
  //             nonlinear least-squares algorithm, acm trans. math. 
  //             software, vol. 7, no. 3. 

  //  ***  general  *** 

  //     coded by david m. gay. 
  //     this subroutine was written in connection with research 
  //     supported by the national science foundation under grants 
  //     mcs-7600324, dcr75-10143, 76-1431ss, mcs76-11989, and 
  //     mcs-7906671. 

  // +++++++++++++++++++++++++++  declarations  ++++++++++++++++++++++++++++

  // dfault... supplies default parameter values. 
  // itsmry... prints iteration summary and info about initial and final x. 
  // nl2itr... reverse-communication routine that carries out nl2sol algo- 
  //             rithm. 
  // rmdcon... returns machine-dependent constants. 
  // vscopy... sets all elements of a vector to a scalar. 

  /* Parameter adjustments */
  --x;
  --iv;
  --v;

  /* Function Body */

  //  ***  iv subscript values  *** 

  // +++++++++++++++++++++++++++++++  body  ++++++++++++++++++++++++++++++++

  d1 = (*n << 1) + 94 + *p * (*p * 3 + 31) / 2;
  iv[d_] = d1;
  r1 = d1 + *p;
  iv[r_] = r1;
  j_1 = r1 + *n;
  iv[j] = j_1;
  rn = j_1 - 1;
  if (iv[1] == 0) {
    nl2sol_defaults(&iv[1], &v[1]);
  }
  //
  //    iv[covreq] = -(i_1 = iv[covreq], abs(i_1));
  //
  iv[covreq] = -abs(iv[covreq]);
  if (iv[covprt] != 0 && iv[covreq] == 0) {
    iv[covreq] = -1;
  }
  strted = true;
  if (iv[1] != 12) {
    goto L80;
  }
  strted = false;
  iv[nfcall] = 1;
  iv[nfgcal] = 1;
  //        ***  initialize scale vector d to ones for computing 
  //        ***  initial jacobian. 
  if (iv[dtype] > 0) {
    vscopy(p, &v[d1], one);
  }
  if (v[dinit] > zero) {
    vscopy(p, &v[d1], v[dinit]);
  }

L10:
  nf = iv[nfcall];
  calcR(functor, function_evaluations, n, p, &x[1], &nf, &v[r1]);
  if (functor.stop()) {
    iv[1] = 11;
    goto L900;
  }
  if (strted) {
    goto L20;
  }
  if (nf > 0) {
    goto L30;
  }
  iv[1] = 13;
  goto L90;

L20:
  if (nf <= 0) {
    iv[toobig] = 1;
  }
  goto L80;

  //  ***  compute finite-difference jacobian  *** 

L30:
  j1k = j_1;
  dk = d1;
  i_1 = *p;
  for (k = 1; k <= i_1; ++k) {
    xk = x[k];
    // Computing MAX 
    //        d_1 = abs(xk), d_2 = one / v[dk];
    //        h_ = v[dltfdj] * std::max(d_1,d_2);
    //        h_ = v[dltfdj] * std::max(d_1,d_2);
    d_1 = fabs(xk);
    d_2 = one / v[dk];
    if (d_2 > d_1)
      d_1 = d_2;
    h_ = v[dltfdj] * d_1;
    ++dk;
  L40:
    x[k] = xk + h_;
    nf = iv[nfgcal];
    calcR(functor, function_evaluations, n, p, &x[1], &nf, &v[j1k]);
    if (functor.stop()) {
      iv[1] = 11;
      goto L900;
    }
    if (nf > 0) {
      goto L50;
    }
    if (hlim == zero) {
      hlim = hfac * rmdcon(3);
    }
    //             ***  hlim = hfac times the unit roundoff  *** 
    h_ = negpt5 * h_;
    //
    //                if (abs(h_) >= hlim) {
    //
    if (fabs(h_) >= hlim) {
      goto L40;
    }

    //
    // PJ, 2024-02-13:
    //
    // Execution can reach this point for well-specified problems when v[xftol], the false convergence threshold, is not
    // greater than DBL_EPSILON since the assessment of false convergence in assess() does effectively not allow for
    // that. However, since we often want to suppress spurious false convergence assessments, we typically have
    // v[xftol] == 0.0 which leads us here. This is because nl2itr/assess rely on false convergence to arise *before*
    // the situation that we are right on the edge of a constraint boundary, and thus will attempt to take a step and
    // directly compute the Jacobian after this step without first checking that the new guess is actually a valid
    // coordinate point.
    //
    x[k] = xk; // Restore the base point of the finite-differencing perturbation.
    {
      int ok = 1;
      calcR(functor, function_evaluations, n, p, &x[1], &ok, &v[r1]);
      if (!ok) {
        //
        // Restore the stored last good point and exit with FALSE_CONVERGENCE.
        //
        const int x0 = 60, rsave = 52;
        vcopy(*p, &x[1], &v[iv[x0]]);
        vcopy(*n, &v[r1], &v[iv[rsave]]);
        iv[1] = 8; // FALSE_CONVERGENCE
#ifdef _DEBUG
        ok = 1;
        calcR(functor, function_evaluations, n, p, &x[1], &ok, &v[r1]);
        ASSERT(ok);
#endif
        goto L999;
      }
    }
    // If we cannot compute the Jacobian but *can* compute the objective, then something is ill specified.
    // In that case, we allow dropping through to UNABLE_TO_COMPUTE_JACOBIAN.
    iv[1] = 15; // UNABLE_TO_COMPUTE_JACOBIAN
    goto L90;

  L50:
    x[k] = xk;
    i_2 = rn;
    for (i_ = r1; i_ <= i_2; ++i_) {
      v[j1k] = (v[j1k] - v[i_]) / h_;
      ++j1k;
    }
  }
  ++(*jacobian_evaluations);

  strted = true;

L80:
  nl2itr(functor, &v[d1], &iv[1], &v[j_1], n, n, p, &v[r1], &v[1], &x[1]);
  if (functor.stop()) {
    iv[1] = 11;
    goto L900;
  }
  if ((i_1 = iv[1] - 2) < 0) {
    goto L10;
  } else if (i_1 == 0) {
    goto L30;
  } else {
    goto L999;
  }

L90:
L900:
L999:
  return 0;
} /* nl2sno */

int nl2itr(NL2SOL::Functor &functor, double *d_, int *iv, double *j, int *n, int *nn, int *p, double *r_, double *v, double *x) {
  /* Initialized data */

  const int jtol1 = 87;
  const int lmax0 = 35;
  const int nvsave = 9;
  const int phmxfc = 21;
  const int preduc = 7;
  const int radfac = 16;
  const int radius = 8;
  const int rad0 = 9;
  const int rlimit = 42;
  const int size = 47;
  const int stppar = 5;
  const int tuner4 = 29;
  const int tuner5 = 30;
  const int vsave1 = 78;
  const int wscale = 48;
  const double half = .5;
  const double negone = -1.;
  const double one = 1.;
  const double zero = 0.;
  const int cnvcod = 34;
  const int covmat = 26;
  const int covprt = 14;
  const int covreq = 15;
  const int dig = 43;
  const int dtype = 16;
  const int g = 28;
  const int h_ = 44;
  const int ierr = 32;
  const int inits = 25;
  const int ipivot = 61;
  const int ipiv0 = 60;
  const int irc = 3;
  const int kagqt = 35;
  const int kalm = 36;
  const int lky = 37;
  const int lmat = 58;
  const int mode = 38;
  const int model = 5;
  const int mxfcal = 17;
  const int mxiter = 18;
  const int nfcall = 6;
  const int nfgcal = 7;
  const int nfcov = 40;
  const int ngcov = 41;
  const int ngcall = 30;
  const int niter = 31;
  const int qtr = 49;
  const int radinc = 8;
  const int rd = 51;
  const int restor = 9;
  const int rsave = 52;
  const int s = 53;
  const int step = 55;
  const int stglim = 11;
  const int stlstg = 56;
  const int sused = 57;
  const int switch_ = 12;
  const int toobig = 2;
  const int w = 59;
  const int xirc = 13;
  const int x0 = 60;
  const int cosmin = 43;
  const int dgnorm = 1;
  const int dinit = 38;
  const int dstnrm = 2;
  const int d0init = 37;
  const int f = 10;
  const int fdif = 11;
  const int fuzz = 45;
  const int f0 = 13;
  const int gtstep = 4;
  const int incfac = 23;
  const int jtinit = 39;

  /* System generated locals */
  int j_dim1, j_offset, i_1, i_2;
  double d_1;

  /* Local variables */
  int pp1o2;
  double rdof1;
  int lmat1, temp1, temp2, ipiv1, step1;
  double e;
  int i_, k, l, m;
  double t;
  int ipivi, ipivk, g1, h0, h1, sstep;
  int s1;
  double t1;
  int w1;
  int rsave1, g01;
  int x01;
  int rd0, im1, rd1, km1, stpmod;
  int lstgst;
  double sttsst;
  int rdk, ipk, smh, dig1, lky1, qtr1;


  //  ***  carry out nl2sol (nonlinear least-squares) iterations  *** 
  //  ***  (nl2sol version 2.2)  *** 

  //  ***  parameter declarations  *** 

  //     dimension iv(60+p), v(93 + 2*n + p*(3*p+31)/2) 


  // --------------------------  parameter usage  --------------------------

  // d.... scale vector. 
  // iv... int value array. 
  // j.... n by p jacobian matrix (lead dimension nn). 
  // n.... number of observations (components in r). 
  // nn... lead dimension of j. 
  // p.... number of parameters (components in x). 
  // r.... residual vector. 
  // v.... floating-point value array. 
  // x.... parameter vector. 

  //  ***  discussion  *** 

  //        parameters iv, n, p, v, and x are the same as the correspond- 
  //     ing ones to nl2sol (which see), except that v can be shorter 
  //     (since the part of v that nl2sol uses for storing d, j, and r is 
  //     not needed).  moreover, compared with nl2sol, iv(1) may have the 
  //     two additional output values 1 and 2, which are explained below, 
  //     as is the use of iv(toobig) and iv(nfgcal).  the values iv(d), 
  //     iv(j), and iv(r), which are output values from nl2sol (and 
  //     nl2sno), are not referenced by nl2itr or the subroutines it calls. 
  //        on a fresh start, i.e., a call on nl2itr with iv(1) = 0 or 12, 
  //     nl2itr assumes that r = r(x), the residual at x, and j = j(x), 
  //     the corresponding jacobian matrix of r at x. 

  // iv(1) = 1 means the caller should set r to r(x), the residual at x, 
  //             and call nl2itr again, having changed none of the other 
  //             parameters.  an exception occurs if r cannot be evaluated 
  //             at x (e.g. if r would overflow), which may happen because 
  //             of an oversized step.  in this case the caller should set 
  //             iv(toobig) = iv(2) to 1, which will cause nl2itr to ig- 
  //             nore r and try a smaller step.  the parameter nf that 
  //             nl2sol passes to calcr (for possible use by calcj) is a 
  //             copy of iv(nfcall) = iv(6). 
  // iv(1) = 2 means the caller should set j to j(x), the jacobian matrix 
  //             of r at x, and call nl2itr again.  the caller may change 
  //             d at this time, but should not change any of the other 
  //             parameters.  the parameter nf that nl2sol passes to 
  //             calcj is iv(nfgcal) = iv(7).  if j cannot be evaluated 
  //             at x, then the caller may set iv(nfgcal) to 0, in which 
  //             case nl2itr will return with iv(1) = 15. 

  //  ***  general  *** 

  //     coded by david m. gay. 
  //     this subroutine was written in connection with research 
  //     supported by the national science foundation under grants 

  //     mcs-7600324, dcr75-10143, 76-1431ss, mcs76-11989, and 
  //     mcs-7906671. 
  //        (see nl2sol for references.) 

  // +++++++++++++++++++++++++++  declarations  ++++++++++++++++++++++++++++

  //  ***  local variables  *** 


  //     ***  constants  *** 

  //  ***  external functions and subroutines  *** 


  // assess... assesses candidate step. 
  // covclc... computes covariance matrix. 
  // dotprd... returns inner product of two vectors. 
  // dupdat... updates scale vector d. 
  // gqtstp... computes goldfeld-quandt-trotter step (augmented model). 
  // itsmry... prints iteration summary and info about initial and final x. 
  // lmstep... computes levenberg-marquardt step (gauss-newton model). 
  // parchk... checks validity of input iv and v values. 
  // qapply... applies orthogonal matrix q from qrfact to a vector. 
  // qrfact... computes qr decomposition of a matrix via householder trans. 
  // rptmul... multiplies vector by the r matrix (and/or its transpose) 
  //             stored by qrfact. 
  // slupdt... performs quasi-newton update on compactly stored lower tri- 
  //             angle of a symmetric matrix. 
  // vaxpy.... computes scalar times one vector plus another. 
  // vcopy.... copies one vector to another. 
  // vscopy... sets all elements of a vector to a scalar. 
  // v2norm... returns the 2-norm of a vector. 

  /* Parameter adjustments */
  --iv;
  --r_;
  --x;
  j_dim1 = *nn;
  j_offset = j_dim1 + 1;
  j -= j_offset;
  --d_;
  --v;

  /* Function Body */

  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  i_ = iv[1];
  if (i_ == 1) {
    goto L20;
  }
  if (i_ == 2) {
    goto L50;
  }

  //  ***  check validity of iv and v input values  *** 

  //     ***  note -- if iv(1) = 0, then parchk calls dfault(iv, v)  *** 
  parchk(&iv[1], n, nn, p, &v[1]);
  i_ = iv[1] - 2;
  if (i_ > 10) {
    goto L999;
  }
  switch (static_cast<int>(i_)) {
  case 1:  goto L350;
  case 2:  goto L350;
  case 3:  goto L350;
  case 4:  goto L350;
  case 5:  goto L350;
  case 6:  goto L350;
  case 7:  goto L195;
  case 8:  goto L160;
  case 9:  goto L195;
  case 10:  goto L10;
  }

  //  ***  initialization and storage allocation  *** 

L10:
  iv[niter] = 0;
  iv[nfcall] = 1;
  iv[ngcall] = 1;
  iv[nfgcal] = 1;
  iv[mode] = -1;
  iv[stglim] = 2;
  iv[toobig] = 0;
  iv[cnvcod] = 0;
  iv[covmat] = 0;
  iv[nfcov] = 0;
  iv[ngcov] = 0;
  iv[kalm] = -1;
  iv[radinc] = 0;
  iv[s] = jtol1 + (*p << 1);
  pp1o2 = *p * (*p + 1) / 2;
  iv[x0] = iv[s] + pp1o2;
  iv[step] = iv[x0] + *p;
  iv[stlstg] = iv[step] + *p;
  iv[dig] = iv[stlstg] + *p;
  iv[g] = iv[dig] + *p;
  iv[lky] = iv[g] + *p;
  iv[rd] = iv[lky] + *p;
  iv[rsave] = iv[rd] + *p;
  iv[qtr] = iv[rsave] + *n;
  iv[h_] = iv[qtr] + *n;
  iv[w] = iv[h_] + pp1o2;
  iv[lmat] = iv[w] + (*p << 2) + 7;
  //     +++ length of w = p*(p+9)/2 + 7.  lmat is contained in w. 
  if (v[dinit] >= zero) {
    vscopy(p, &d_[1], v[dinit]);
  }
  if (v[jtinit] > zero) {
    vscopy(p, &v[jtol1], v[jtinit]);
  }
  i_ = jtol1 + *p;
  if (v[d0init] > zero) {
    vscopy(p, &v[i_], v[d0init]);
  }
  v[rad0] = zero;
  v[stppar] = zero;
  v[radius] = v[lmax0] / (one + v[phmxfc]);

  //  ***  set initial model and s matrix  *** 

  iv[model] = 1;
  if (iv[inits] == 2) {
    iv[model] = 2;
  }
  s1 = iv[s];
  if (iv[inits] == 0) {
    vscopy(&pp1o2, &v[s1], zero);
  }

  //  ***  compute function value (half the sum of squares)  *** 

L20:
  t = v2norm(n, &r_[1]);
  if (t > v[rlimit]) {
    iv[toobig] = 1;
  }
  if (iv[toobig] != 0) {
    goto L30;
  }
  // Computing 2nd power 
  d_1 = t;
  v[f] = half * (d_1 * d_1);
L30:
  if ((i_1 = iv[mode]) < 0) {
    goto L40;
  } else if (i_1 == 0) {
    goto L350;
  } else {
    goto L730;
  }

L40:
  if (iv[toobig] == 0) {
    goto L60;
  }
  iv[1] = 13;
  goto L900;

  //  ***  make sure jacobian could be computed  *** 

L50:
  if (iv[nfgcal] != 0) {
    goto L60;
  }
  iv[1] = 15;
  goto L900;

  //  ***  compute gradient  *** 

L60:
  iv[kalm] = -1;
  g1 = iv[g];
  i_1 = *p;
  for (i_ = 1; i_ <= i_1; ++i_) {
    ASSERT(isFinite(r_[1]));
    ASSERT(isFinite(j[i_ * j_dim1 + 1]));
    v[g1] = dotprd(n, &r_[1], &j[i_ * j_dim1 + 1]);
    ASSERT(isFinite(v[g1]));
    ++g1;
  }
  if (iv[mode] > 0) {
    goto L710;
  }

  //  ***  update d and make copies of r for possible use later  *** 

  if (iv[dtype] > 0) {
    dupdat(&d_[1], &iv[1], &j[j_offset], n, nn, p, &v[1]);
  }
  rsave1 = iv[rsave];
  vcopy(*n, &v[rsave1], &r_[1]);
  qtr1 = iv[qtr];
  vcopy(*n, &v[qtr1], &r_[1]);

  //  ***  compute  d**-1 * gradient  *** 

  g1 = iv[g];
  dig1 = iv[dig];
  k = dig1;
  i_1 = *p;
  for (i_ = 1; i_ <= i_1; ++i_) {
    v[k] = v[g1] / d_[i_];
    ++k;
    ++g1;
  }
  v[dgnorm] = v2norm(p, &v[dig1]);

  if (iv[cnvcod] != 0) {
    goto L700;
  }
  if (iv[mode] == 0) {
    goto L570;
  }
  iv[mode] = 0;


  // -----------------------------  main loop  -----------------------------


  //  ***  print iteration summary, check iteration limit  *** 

L150:
L160:
  k = iv[niter];
  if (k < iv[mxiter]) {
    goto L170;
  }
  iv[1] = 10;
  goto L900;
L170:
  iv[niter] = k + 1;

  //  ***  update radius  *** 

  if (k == 0) {
    goto L185;
  }
  step1 = iv[step];
  i_1 = *p;
  for (i_ = 1; i_ <= i_1; ++i_) {
    v[step1] = d_[i_] * v[step1];
    ++step1;
  }
  step1 = iv[step];
  v[radius] = v[radfac] * v2norm(p, &v[step1]);

  //  ***  initialize for start of next iteration  *** 

L185:
  x01 = iv[x0];
  v[f0] = v[f];
  iv[kagqt] = -1;
  iv[irc] = 4;
  //
  //    iv[h_] = -(i_1 = iv[h_], abs(i_1));
  //
  iv[h_] = -abs(iv[h_]);
  iv[sused] = iv[model];

  //     ***  copy x to x0  *** 

  vcopy(*p, &v[x01], &x[1]);

  //  ***  check stopx and function evaluation limit  *** 

L190:
  if (!functor.stop()) {
    goto L200;
  }
  iv[1] = 11;
  goto L205;

  //     ***  come here when restarting after func. eval. limit or stopx. 

L195:
  if (v[f] >= v[f0]) {
    goto L200;
  }
  v[radfac] = one;
  k = iv[niter];
  goto L170;

L200:
  if (iv[nfcall] < iv[mxfcal] + iv[nfcov]) {
    goto L210;
  }
  iv[1] = 9;
L205:
  if (functor.stop() || v[f] >= v[f0]) {
    goto L900;
  }

  //        ***  in case of stopx or function evaluation limit with 
  //        ***  improved v(f), evaluate the gradient at x. 

  iv[cnvcod] = iv[1];
  goto L560;

  // . . . . . . . . . . . . .  compute candidate step  . . . . . . . . . . 

L210:
  step1 = iv[step];
  w1 = iv[w];
  if (iv[model] == 2) {
    goto L240;
  }

  //  ***  compute levenberg-marquardt step  *** 

  qtr1 = iv[qtr];
  if (iv[kalm] >= 0) {
    goto L215;
  }
  rd1 = iv[rd];
  if (-1 == iv[kalm]) {
    qrfact(nn, n, p, &j[j_offset], &v[rd1], &iv[ipivot], &iv[ierr], 0, &v[w1]);
  }
  qapply(nn, n, p, &j[j_offset], &v[qtr1], &iv[ierr]);
L215:
  h1 = iv[h_];
  if (h1 > 0) {
    goto L230;
  }

  //        ***  copy r matrix to h  *** 

  h1 = -h1;
  iv[h_] = h1;
  k = h1;
  rd1 = iv[rd];
  v[k] = v[rd1];
  if (*p == 1) {
    goto L230;
  }
  i_1 = *p;
  for (i_ = 2; i_ <= i_1; ++i_) {
    i_2 = i_ - 1;
    vcopy(i_2, &v[k + 1], &j[i_ * j_dim1 + 1]);
    k += i_;
    ++rd1;
    v[k] = v[rd1];
  }

L230:
  g1 = iv[g];
  lmstep(&d_[1], &v[g1], &iv[ierr], &iv[ipivot], &iv[kalm], p, &v[qtr1], &v[h1], &v[step1], &v[1], &v[w1]);
  goto L310;

  //  ***  compute goldfeld-quandt-trotter step (augmented model)  *** 

L240:
  if (iv[h_] > 0) {
    goto L300;
  }

  //     ***  set h to  d**-1 * ( (j**t)*j + s) ) * d**-1.  *** 

  h1 = -iv[h_];
  iv[h_] = h1;
  s1 = iv[s];
  if (-1 != iv[kalm]) {
    goto L270;
  }

  //        ***  j is in its original form  *** 

  i_1 = *p;
  for (i_ = 1; i_ <= i_1; ++i_) {
    t = one / d_[i_];
    i_2 = i_;
    for (k = 1; k <= i_2; ++k) {
      v[h1] = t * (dotprd(n, &j[i_ * j_dim1 + 1], &j[k * j_dim1 + 1])
        + v[s1]) / d_[k];
      ++h1;
      ++s1;
    }
  }
  goto L300;

  //  ***  lmstep has applied qrfact to j  *** 

L270:
  smh = s1 - h1;
  h0 = h1 - 1;
  ipiv1 = iv[ipivot];
  t1 = one / d_[ipiv1];
  rd0 = iv[rd] - 1;
  rdof1 = v[rd0 + 1];
  i_1 = *p;
  for (i_ = 1; i_ <= i_1; ++i_) {
    l = ipiv0 + i_;
    ipivi = iv[l];
    h1 = h0 + ipivi * (ipivi - 1) / 2;
    l = h1 + ipivi;
    m = l + smh;
    //             ***  v(l) = h(ipivot(i), ipivot(i))  *** 
    //             ***  v(m) = s(ipivot(i), ipivot(i))  *** 
    t = one / d_[ipivi];
    rdk = rd0 + i_;
    // Computing 2nd power 
    d_1 = v[rdk];
    e = d_1 * d_1;
    if (i_ > 1) {
      i_2 = i_ - 1;
      e += dotprd(&i_2, &j[i_ * j_dim1 + 1], &j[i_ * j_dim1 + 1]);
    }
    // Computing 2nd power 
    d_1 = t;
    v[l] = (e + v[m]) * (d_1 * d_1);
    if (i_ == 1) {
      goto L290;
    }
    l = h1 + ipiv1;
    if (ipivi < ipiv1) {
      l += (ipiv1 - ipivi) * (ipiv1 + ipivi - 3) / 2;
    }
    m = l + smh;
    //             ***  v(l) = h(ipivot(i), ipivot(1))  *** 
    //             ***  v(m) = s(ipivot(i), ipivot(1))  *** 
    v[l] = t * (rdof1 * j[i_ * j_dim1 + 1] + v[m]) * t1;
    if (i_ == 2) {
      goto L290;
    }
    im1 = i_ - 1;
    i_2 = im1;
    for (k = 2; k <= i_2; ++k) {
      ipk = ipiv0 + k;
      ipivk = iv[ipk];
      l = h1 + ipivk;
      if (ipivi < ipivk) {
        l += (ipivk - ipivi) * (ipivk + ipivi - 3) / 2;
      }
      m = l + smh;
      //                  ***  v(l) = h(ipivot(i), ipivot(k))  *** 
      //                  ***  v(m) = s(ipivot(i), ipivot(k))  *** 
      km1 = k - 1;
      rdk = rd0 + k;
      v[l] = t * (dotprd(&km1, &j[i_ * j_dim1 + 1], &j[k * j_dim1 + 1]
      ) + v[rdk] * j[k + i_ * j_dim1] + v[m]) / d_[ipivk];
    }
  L290:
    ;
  }

  //  ***  compute actual goldfeld-quandt-trotter step  *** 

L300:
  h1 = iv[h_];
  dig1 = iv[dig];
  lmat1 = iv[lmat];
  gqtstp(&d_[1], &v[dig1], &v[h1], &iv[kagqt], &v[lmat1], p, &v[step1], &v[1], &v[w1]);


  //  ***  compute r(x0 + step)  *** 

L310:
  if (iv[irc] == 6) {
    goto L350;
  }
  x01 = iv[x0];
  step1 = iv[step];
  vaxpy(p, &x[1], one, &v[step1], &v[x01]);
  ++iv[nfcall];
  iv[1] = 1;
  iv[toobig] = 0;
  goto L999;

  // . . . . . . . . . . . . .  assess candidate step  . . . . . . . . . . .

L350:
  step1 = iv[step];
  lstgst = iv[stlstg];
  x01 = iv[x0];
  assess(&d_[1], &iv[1], p, &v[step1], &v[lstgst], &v[1], &x[1], &v[x01]);

  //  ***  if necessary, switch models and/or restore r  *** 

  if (iv[switch_] == 0) {
    goto L360;
  }
  //
  //    iv[h_] = -(i_1 = iv[h_], abs(i_1));
  //
  iv[h_] = -abs(iv[h_]);
  iv[sused] += 2;
  vcopy(nvsave, &v[1], &v[vsave1]);
L360:
  if (iv[restor] == 0) {
    goto L390;
  }
  rsave1 = iv[rsave];
  vcopy(*n, &r_[1], &v[rsave1]);
L390:
  l = iv[irc] - 4;
  stpmod = iv[model];
  if (l > 0) {
    switch (static_cast<int>(l)) {
    case 1:  goto L410;
    case 2:  goto L440;
    case 3:  goto L450;
    case 4:  goto L450;
    case 5:  goto L450;
    case 6:  goto L450;
    case 7:  goto L450;
    case 8:  goto L450;
    case 9:  goto L640;
    case 10:  goto L570;
    }
  }

  //  ***  decide whether to change models  *** 

  e = v[preduc] - v[fdif];
  sstep = iv[lky];
  s1 = iv[s];
  slvmul(p, &v[sstep], &v[s1], &v[step1]);
  sttsst = half * dotprd(p, &v[step1], &v[sstep]);
  if (iv[model] == 1) {
    sttsst = -sttsst;
  }
  //
  //    if ((d_1 = e + sttsst, abs(d_1)) * v[fuzz] >= abs(e)) {
  //
  if (fabs(e + sttsst) * v[fuzz] >= fabs(e)) {
    goto L400;
  }

  //     ***  switch models  *** 

  iv[model] = 3 - iv[model];
  if (iv[model] == 1) {
    iv[kagqt] = -1;
  }
  if (iv[model] == 2 && iv[kalm] > 0) {
    iv[kalm] = 0;
  }
  if (-2 < l) {
    goto L480;
  }
  //
  //    iv[h_] = -(i_1 = iv[h_], abs(i_1));
  //
  iv[h_] = -abs(iv[h_]);
  iv[sused] += 2;
  vcopy(nvsave, &v[vsave1], &v[1]);
  goto L420;

L400:
  if (-3 < l) {
    goto L480;
  }

  //     ***  recompute step with decreased radius  *** 

  v[radius] = v[radfac] * v[dstnrm];
  goto L190;

  //  ***  recompute step, saving v values and r if necessary  *** 

L410:
  v[radius] = v[radfac] * v[dstnrm];
L420:
  if (v[f] >= v[f0]) {
    goto L190;
  }
  rsave1 = iv[rsave];
  vcopy(*n, &v[rsave1], &r_[1]);
  goto L190;

  //  ***  compute step of length v(lmax0) for singular convergence test 

L440:
  v[radius] = v[lmax0];
  goto L210;

  //  ***  convergence or false convergence  *** 

L450:
  iv[cnvcod] = l;
  if (v[f] >= v[f0]) {
    goto L700;
  }
  if (iv[xirc] == 14) {
    goto L700;
  }
  iv[xirc] = 14;

  // . . . . . . . . . . . .  process acceptable step  . . . . . . . . . . .

L480:
  iv[covmat] = 0;

  //  ***  set  lky = (j(x0)**t) * r(x)  *** 

  lky1 = iv[lky];
  if (iv[kalm] >= 0) {
    goto L500;
  }

  //     ***  jacobian has not been modified  *** 

  i_1 = *p;
  for (i_ = 1; i_ <= i_1; ++i_) {
    v[lky1] = dotprd(n, &j[i_ * j_dim1 + 1], &r_[1]);
    ++lky1;
  }
  goto L510;

  //  ***  qrfact has been applied to j.  store copy of r in qtr and  *** 
  //  ***  apply q to it.                                             *** 

L500:
  qtr1 = iv[qtr];
  vcopy(*n, &v[qtr1], &r_[1]);
  qapply(nn, n, p, &j[j_offset], &v[qtr1], &iv[ierr]);

  //  ***  multiply top p-vector in qtr by permuted upper triangle    *** 
  //  ***  stored by qrfact in j and rd.                              *** 

  rd1 = iv[rd];
  temp1 = iv[stlstg];
  rptmul(3, &iv[ipivot], &j[j_offset], nn, p, &v[rd1], &v[qtr1], &v[lky1], &v[temp1]);

  //  ***  see whether to set v(radfac) by gradient tests  *** 

L510:
  if (iv[irc] != 3) {
    goto L560;
  }
  step1 = iv[step];
  temp1 = iv[stlstg];
  temp2 = iv[x0];

  //     ***  set  temp1 = hessian * step  for use in gradient tests  *** 

  if (stpmod == 2) {
    goto L530;
  }

  //        ***  step computed using gauss-newton model  *** 
  //        ***  -- qrfact has been applied to j         *** 

  rd1 = iv[rd];
  rptmul(2, &iv[ipivot], &j[j_offset], nn, p, &v[rd1], &v[step1], &v[
    temp1], &v[temp2]);
  goto L560;

  //     ***  step computed using augmented model  *** 

L530:
  h1 = iv[h_];
  k = temp2;
  i_1 = *p;
  for (i_ = 1; i_ <= i_1; ++i_) {
    v[k] = d_[i_] * v[step1];
    ++k;
    ++step1;
  }
  slvmul(p, &v[temp1], &v[h1], &v[temp2]);
  i_1 = *p;
  for (i_ = 1; i_ <= i_1; ++i_) {
    v[temp1] = d_[i_] * v[temp1];
    ++temp1;
  }

  //  ***  save old gradient and compute new one  *** 

L560:
  ++iv[ngcall];
  g1 = iv[g];
  g01 = iv[w];
  vcopy(*p, &v[g01], &v[g1]);
  iv[1] = 2;
  goto L999;

  //  ***  initializations -- g0 = g - g0, etc.  *** 

L570:
  g01 = iv[w];
  g1 = iv[g];
  vaxpy(p, &v[g01], negone, &v[g01], &v[g1]);
  step1 = iv[step];
  temp1 = iv[stlstg];
  temp2 = iv[x0];
  if (iv[irc] != 3) {
    goto L600;
  }

  //  ***  set v(radfac) by gradient tests  *** 

  //     ***  set  temp1 = d**-1 * (hessian * step  +  (g(x0) - g(x)))  *** 

  k = temp1;
  l = g01;
  i_1 = *p;
  for (i_ = 1; i_ <= i_1; ++i_) {
    v[k] = (v[k] - v[l]) / d_[i_];
    ++k;
    ++l;
  }

  //        ***  do gradient tests  *** 

  if (v2norm(p, &v[temp1]) <= v[dgnorm] * v[tuner4]) {
    goto L590;
  }
  if (dotprd(p, &v[g1], &v[step1]) >= v[gtstep] * v[tuner5]) {
    goto L600;
  }
L590:
  v[radfac] = v[incfac];

  //  ***  finish computing lky = ((j(x) - j(x0))**t) * r  *** 

  //     ***  currently lky = (j(x0)**t) * r  *** 

L600:
  lky1 = iv[lky];
  vaxpy(p, &v[lky1], negone, &v[lky1], &v[g1]);

  //  ***  determine sizing factor v(size)  *** 

  //     ***  set temp1 = s * step  *** 
  s1 = iv[s];
  slvmul(p, &v[temp1], &v[s1], &v[step1]);

  //
  //    t1 = (d_1 = dotprd(p, &v[step1], &v[temp1]), abs(d_1));
  //
  t1 = fabs(dotprd(p, &v[step1], &v[temp1]));
  //
  //    t = (d_1 = dotprd(p, &v[step1], &v[lky1]), abs(d_1));
  //
  t = fabs(dotprd(p, &v[step1], &v[lky1]));
  v[size] = one;
  if (t < t1) {
    v[size] = t / t1;
  }

  //  ***  update s  *** 

  slupdt(&v[s1], &v[cosmin], p, &v[size], &v[step1], &v[temp1], &v[temp2], &v[g01], &v[wscale], &v[lky1]);
  iv[1] = 2;
  goto L150;

  // . . . . . . . . . . . . . .  misc. details  . . . . . . . . . . . . . .

  //  ***  bad parameters to assess  *** 

L640:
  iv[1] = 14;
  goto L900;

  //  ***  convergence obtained -- compute covariance matrix if desired *** 

L700:
  if (iv[covreq] == 0 && iv[covprt] == 0) {
    goto L760;
  }
  if (iv[covmat] != 0) {
    goto L760;
  }
  if (iv[cnvcod] >= 7) {
    goto L760;
  }
  iv[mode] = 0;
L710:
  covclc(&i_, &d_[1], &iv[1], &j[j_offset], n, nn, p, &r_[1], &v[1], &x[1]);
  switch (static_cast<int>(i_)) {
  case 1:  goto L720;
  case 2:  goto L720;
  case 3:  goto L740;
  case 4:  goto L750;
  }
L720:
  ++iv[nfcov];
  ++iv[nfcall];
  iv[restor] = i_;
  iv[1] = 1;
  goto L999;

L730:
  if (iv[restor] == 1 || iv[toobig] != 0) {
    goto L710;
  }
  iv[nfgcal] = iv[nfcall];
L740:
  ++iv[ngcov];
  ++iv[ngcall];
  iv[1] = 2;
  goto L999;

L750:
  iv[mode] = 0;
  if (iv[niter] == 0) {
    iv[mode] = -1;
  }

L760:
  iv[1] = iv[cnvcod];
  iv[cnvcod] = 0;

  //  ***  print summary of final iteration and other requested items  *** 

L900:
L999:
  return 0;

} /* nl2itr */

int covclc(int *covirc, double *d_, int *iv, double *j, int *n, int *nn, int *p, double *r_, double *v, double *x) {
  /* Initialized data */

  const double half = .5;
  const double negpt5 = -.5;
  const double one = 1.;
  const double two = 2.;
  const double zero = 0.;
  const int covmat = 26;
  const int covreq = 15;
  const int delta = 50;
  const int delta0 = 44;
  const int dltfdc = 40;
  const int f = 10;
  const int fx = 46;
  const int g = 28;
  const int h_ = 44;
  const int ierr = 32;
  const int ipivot = 61;
  const int ipiv0 = 60;
  const int kagqt = 35;
  const int kalm = 36;
  const int lmat = 58;
  const int mode = 38;
  const int nfgcal = 7;
  const int qtr = 49;
  const int rd = 51;
  const int rsave = 52;
  const int savei = 54;
  const int switch_ = 12;
  const int toobig = 2;
  const int w = 59;
  const int xmsave = 49;

  /* System generated locals */
  int j_dim1, j_offset, i_1, i_2, i_3;

  /* Local variables */
  int kind, mm1o2, pp1o2, stpi, stpm, i_, k, l, m;
  double t;
  bool havej;
  int ipivi, ipivk, g1;
  int w0, w1, gsave1, hc, gp, kl;
  double wk;
  int wl;
  int rd1;
  int ip1, mm1;
  double del;
  int hmi, irc, hpi, hpm, cov, stp0, qtr1;


  //  ***  compute covariance matrix for nl2itr (nl2sol version 2.2)  *** 

  //  ***  let k = iabs(iv(covreq).  for k .le. 2, a finite-difference 
  //  ***  hessian h is computed (using func. and grad. values if 
  //  ***  iv(covreq) is nonnegative, and using only func. values if 
  //  ***  iv(covreq) is negative).  for scale = 2*f(x) / std::max(1, n-p), 
  //  ***  where 2*f(x) is the residual sum of squares, covclc computes... 
  //  ***             k = 0 or 1...  scale * h**-1 * (j**t * j) * h**-1. 
  //  ***             k = 2...  scale * h**-1. 
  //  ***             k .ge. 3...  scale * (j**t * j)**-1. 

  //  ***  external subroutines  *** 

  // linvrt... invert lower triangular matrix. 
  // litvmu... apply inverse-transpose of compact lower triang. matrix. 
  // livmul... apply inverse of compact lower triang. matrix. 
  // lsqrt.... compute cholesky factor of (lower triang. of) a sym. matrix. 
  // ltsqar... given lower triang. matrix l, compute (l**t)*l. 
  // qrfact... compute qr decomposition of a matrix. 
  // vcopy.... copy one vector to another. 
  // vscopy... set all elements of a vector to a scalar. 

  //  ***  subscripts for iv and v  *** 

  /* Parameter adjustments */
  --iv;
  --r_;
  --x;
  j_dim1 = *nn;
  j_offset = j_dim1 + 1;
  j -= j_offset;
  --d_;
  --v;

  /* Function Body */

  // +++++++++++++++++++++++++++++++  body  ++++++++++++++++++++++++++++++++

  *covirc = 4;
  kind = iv[covreq];
  m = iv[mode];
  if (m > 0) {
    goto L10;
  }
  iv[kagqt] = -1;
  if (iv[kalm] > 0) {
    iv[kalm] = 0;
  }
  if (abs(kind) >= 3) {
    goto L300;
  }
  v[fx] = v[f];
  k = iv[rsave];
  vcopy(*n, &v[k], &r_[1]);
L10:
  if (m > *p) {
    goto L200;
  }
  if (kind < 0) {
    goto L100;
  }

  //  ***  compute finite-difference hessian using both function and 
  //  ***  gradient values. 

  gsave1 = iv[w] + *p;
  g1 = iv[g];
  if (m > 0) {
    goto L15;
  }
  //        ***  first call on covclc.  set gsave = g, take first step  *** 
  vcopy(*p, &v[gsave1], &v[g1]);
  iv[switch_] = iv[nfgcal];
  goto L80;

L15:
  del = v[delta];
  x[m] = v[xmsave];
  if (iv[toobig] == 0) {
    goto L30;
  }

  //     ***  handle oversize v(delta)  *** 

  if (del * x[m] > zero) {
    goto L20;
  }
  //             ***  we already tried shrinking v(delta), so quit  *** 
  iv[covmat] = -2;
  goto L190;

  //        ***  try shrinking v(delta)  *** 
L20:
  del = negpt5 * del;
  goto L90;

L30:
  cov = iv[lmat];
  gp = g1 + *p - 1;

  //  ***  set  g = (g - gsave)/del  *** 

  i_1 = gp;
  for (i_ = g1; i_ <= i_1; ++i_) {
    v[i_] = (v[i_] - v[gsave1]) / del;
    ++gsave1;
  }

  //  ***  add g as new col. to finite-diff. hessian matrix  *** 

  k = cov + m * (m - 1) / 2;
  l = k + m - 2;
  if (m == 1) {
    goto L60;
  }

  //  ***  set  h(i,m) = 0.5 * (h(i,m) + g(i))  for i = 1 to m-1  *** 

  i_1 = l;
  for (i_ = k; i_ <= i_1; ++i_) {
    v[i_] = half * (v[i_] + v[g1]);
    ++g1;
  }

  //  ***  add  h(i,m) = g(i)  for i = m to p  *** 

L60:
  ++l;
  i_1 = *p;
  for (i_ = m; i_ <= i_1; ++i_) {
    v[l] = v[g1];
    l += i_;
    ++g1;
  }

L80:
  ++m;
  iv[mode] = m;
  if (m > *p) {
    goto L190;
  }

  //  ***  choose next finite-difference step, return to get g there  *** 

  // Computing MAX 
  //
  //    d_2 = one / d_[m], d_3 = (d_1 = x[m], abs(d_1));
  //    del = v[delta0] * std::max(d_2,d_3);
  //
  del = v[delta0] * std::max(one / d_[m], fabs(x[m]));

  if (x[m] < zero) {
    del = -del;
  }
  v[xmsave] = x[m];
L90:
  x[m] += del;
  v[delta] = del;
  *covirc = 2;
  goto L999;

  //  ***  compute finite-difference hessian using function values only. 

L100:
  stp0 = iv[w] + *p - 1;
  mm1 = m - 1;
  mm1o2 = m * mm1 / 2;
  if (m > 0) {
    goto L105;
  }
  //        ***  first call on covclc.  *** 
  iv[savei] = 0;
  goto L180;

L105:
  i_ = iv[savei];
  // PJ: bug fix
  cov = iv[lmat];
  if (i_ > 0) {
    goto L160;
  }
  if (iv[toobig] == 0) {
    goto L120;
  }

  //     ***  handle oversize step  *** 

  stpm = stp0 + m;
  del = v[stpm];

  //
  // PJ:  BUG ALERT !  BUG ALERT !  BUG ALERT ! - this was wrong in the original version.
  //
  //
  // x() is dimensioned to have p elements !! The saved value is in v(xmsave) !
  //

  if (del * v[xmsave] > zero) {
    goto L110;
  }
  //             ***  we already tried shrinking the step, so quit  *** 
  iv[covmat] = -2;
  goto L999;

  //        ***  try shrinking the step  *** 
L110:
  del = negpt5 * del;

  // PJ:  BUG ALERT !  BUG ALERT !  BUG ALERT ! - this was wrong in the original version.
  //
  // x() is dimensioned to have p elements !! The saved value is in v(xmsave) !
  //
  //         x(m) = x(xmsave) + del 

  x[m] = v[xmsave] + del;
  v[stpm] = del;
  *covirc = 1;
  goto L999;

  //  ***  save f(x + stp(m)*e(m)) in h(p,m)  *** 

L120:
  pp1o2 = *p * (*p - 1) / 2;
  cov = iv[lmat];
  hpm = cov + pp1o2 + mm1;
  v[hpm] = v[f];

  //  ***  start computing row m of the finite-difference hessian h.  *** 

  hmi = cov + mm1o2;
  if (mm1 == 0) {
    goto L140;
  }
  hpi = cov + pp1o2;
  i_1 = mm1;
  for (i_ = 1; i_ <= i_1; ++i_) {
    v[hmi] = v[fx] - (v[f] + v[hpi]);
    ++hmi;
    ++hpi;
  }
L140:
  v[hmi] = v[f] - two * v[fx];

  //  ***  compute function values needed to complete row m of h.  *** 

  i_ = 1;

L150:
  iv[savei] = i_;
  stpi = stp0 + i_;
  v[delta] = x[i_];
  x[i_] += v[stpi];
  if (i_ == m) {
    x[i_] = v[xmsave] - v[stpi];
  }
  *covirc = 1;
  goto L999;

L160:
  x[i_] = v[delta];
  if (iv[toobig] == 0) {
    goto L170;
  }
  //        ***  punt in the event of an oversize step  *** 
  iv[covmat] = -2;
  goto L999;

  //  ***  finish computing h(m,i)  *** 

L170:
  stpi = stp0 + i_;
  hmi = cov + mm1o2 + i_ - 1;
  stpm = stp0 + m;
  v[hmi] = (v[hmi] + v[f]) / (v[stpi] * v[stpm]);
  ++i_;
  if (i_ <= m) {
    goto L150;
  }
  iv[savei] = 0;
  x[m] = v[xmsave];

L180:
  ++m;
  iv[mode] = m;
  if (m > *p) {
    goto L190;
  }

  //  ***  prepare to compute row m of the finite-difference hessian h. 
  //  ***  compute m-th step size stp(m), then return to obtain 
  //  ***  f(x + stp(m)*e(m)), where e(m) = m-th std. unit vector. 

  // Computing MAX 
  //
  //    d_2 = one / d_[m], d_3 = (d_1 = x[m], abs(d_1));
  //    del = v[dltfdc] * std::max(d_2,d_3);
  //
  del = v[dltfdc] * std::max(one / d_[m], fabs(x[m]));
  if (x[m] < zero) {
    del = -del;
  }
  v[xmsave] = x[m];
  x[m] += del;
  stpm = stp0 + m;
  v[stpm] = del;
  *covirc = 1;
  goto L999;

  //  ***  restore r, v(f), etc.  *** 

L190:
  k = iv[rsave];
  vcopy(*n, &r_[1], &v[k]);
  v[f] = v[fx];
  if (kind < 0) {
    goto L200;
  }
  iv[nfgcal] = iv[switch_];
  qtr1 = iv[qtr];
  vcopy(*n, &v[qtr1], &r_[1]);
  if (iv[covmat] < 0) {
    goto L999;
  }
  *covirc = 3;
  goto L999;

L200:
  cov = iv[lmat];

  //  ***  the complete finite-diff. hessian is now stored at v(cov).   *** 
  //  ***  use it to compute the requested covariance matrix.           *** 
  //
  //     ***  compute cholesky factor c of h = c*(c**t)  *** 
  //     ***  and store it at v(hc).  *** 

  hc = cov;
  if (abs(kind) == 2) {
    goto L210;
  }
  //
  //    hc = (i_1 = iv[h_], abs(i_1));
  //
  hc = abs(iv[h_]);
  iv[h_] = -hc;
L210:
  lsqrt(1, p, &v[hc], &v[cov], &irc);
  iv[covmat] = -1;
  if (irc != 0) {
    goto L999;
  }

  w1 = iv[w] + *p;
  if (abs(kind) > 1) {
    goto L350;
  }

  //  ***  covariance = scale * h**-1 * (j**t * j) * h**-1  *** 

  i_1 = *p * (*p + 1) / 2;
  vscopy(&i_1, &v[cov], zero);
  havej = iv[kalm] == -1;
  //     ***  havej = .true. means j is in its original form, while 
  //     ***  havej = .false. means qrfact has been applied to j. 

  m = *p;
  if (havej) {
    m = *n;
  }
  w0 = w1 - 1;
  rd1 = iv[rd];
  i_1 = m;
  for (i_ = 1; i_ <= i_1; ++i_) {
    if (havej) {
      goto L240;
    }

    //        ***  set w = ipivot * (row i of r matrix from qrfact).  *** 

    vscopy(p, &v[w1], zero);
    ipivi = ipiv0 + i_;
    l = w0 + iv[ipivi];
    v[l] = v[rd1];
    ++rd1;
    if (i_ == *p) {
      goto L260;
    }
    ip1 = i_ + 1;
    i_2 = *p;
    for (k = ip1; k <= i_2; ++k) {
      ipivk = ipiv0 + k;
      l = w0 + iv[ipivk];
      v[l] = j[i_ + k * j_dim1];
    }
    goto L260;

    //        ***  set w = (row i of j).  *** 

  L240:
    l = w0;
    i_2 = *p;
    for (k = 1; k <= i_2; ++k) {
      ++l;
      v[l] = j[i_ + k * j_dim1];
    }

    //        ***  set w = h**-1 * w.  *** 

  L260:
    livmul(p, &v[w1], &v[hc], &v[w1]);
    litvmu(p, &v[w1], &v[hc], &v[w1]);

    //        ***  add  w * w**t  to covariance matrix.  *** 

    kl = cov;
    i_2 = *p;
    for (k = 1; k <= i_2; ++k) {
      l = w0 + k;
      wk = v[l];
      i_3 = k;
      for (l = 1; l <= i_3; ++l) {
        wl = w0 + l;
        v[kl] += wk * v[wl];
        ++kl;
      }
    }
  }
  goto L380;

  //  ***  covariance = scale * (j**t * j)**-1.  *** 

L300:
  rd1 = iv[rd];
  if (iv[kalm] != -1) {
    goto L310;
  }

  //        ***  apply qrfact to j  *** 

  qtr1 = iv[qtr];
  vcopy(*n, &v[qtr1], &r_[1]);
  w1 = iv[w] + *p;

  qrfact(nn, n, p, &j[j_offset], &v[rd1], &iv[ipivot], &iv[ierr], 0, &v[w1]);
  iv[kalm] = -2;
L310:
  iv[covmat] = -1;
  if (iv[ierr] != 0) {
    goto L999;
  }
  cov = iv[lmat];
  //
  //    hc = (i_1 = iv[h_], abs(i_1));
  //
  hc = abs(iv[h_]);
  iv[h_] = -hc;

  //     ***  set hc = (r matrix from qrfact).  *** 

  l = hc;
  i_1 = *p;
  for (i_ = 1; i_ <= i_1; ++i_) {
    if (i_ > 1) {
      i_2 = i_ - 1;
      vcopy(i_2, &v[l], &j[i_ * j_dim1 + 1]);
    }
    l = l + i_ - 1;
    v[l] = v[rd1];
    ++l;
    ++rd1;
  }

  //  ***  the cholesky factor c of the unscaled inverse covariance matrix 
  //  ***  (or permutation thereof) is stored at v(hc). 

  //  ***  set c = c**-1. 

L350:
  linvrt(p, &v[hc], &v[hc]);

  //  ***  set c = c**t * c. 

  ltsqar(p, &v[hc], &v[hc]);

  if (hc == cov) {
    goto L380;
  }

  //     ***  c = permuted, unscaled covariance. 
  //     ***  set cov = ipivot * c * ipivot**t. 

  i_1 = *p;
  for (i_ = 1; i_ <= i_1; ++i_) {
    m = ipiv0 + i_;
    ipivi = iv[m];
    kl = cov - 1 + ipivi * (ipivi - 1) / 2;
    i_2 = i_;
    for (k = 1; k <= i_2; ++k) {
      m = ipiv0 + k;
      ipivk = iv[m];
      l = kl + ipivk;

      if (ipivk > ipivi) {
        l += (ipivk - ipivi) * (ipivk + ipivi - 3) / 2;
      }
      v[l] = v[hc];
      ++hc;
    }
  }

L380:
  iv[covmat] = cov;

  //  ***  apply scale factor = (resid. sum of squares) / std::max(1,n-p). 

  // Computing MAX 
  i_1 = 1, i_2 = *n - *p;
  t = v[f] / (half * (std::max(i_1, i_2)));
  k = cov - 1 + *p * (*p + 1) / 2;
  i_1 = k;
  for (i_ = cov; i_ <= i_1; ++i_) {
    v[i_] = t * v[i_];
  }

L999:
  return 0;
} /* covclc */

int gqtstp(double* d_, double* dig, double* dihdi, int* ka, double* l, int* p, double* step, double* v, double* w) {
  /* Initialized data */

  double alphak = 0.;
  double dst = 0.;
  double phi = 0.;
  double uk = 0.;
  const int dgnorm = 1;
  const int dstnrm = 2;
  const int dst0 = 3;
  const int epslon = 19;
  const int gtstep = 4;
  const int nreduc = 6;
  const int phmnfc = 20;
  const int phmxfc = 21;
  const int preduc = 7;
  const int radius = 8;
  const int rad0 = 9;
  const int stppar = 5;
  const double epsfac = 50.;
  const double four = 4.;
  const double half = .5;
  const double kappa = 2.;
  const double negone = -1.;
  const double one = 1.;
  const double six = 6.;
  const double three = 3.;
  const double two = 2.;
  const double zero = 0.;
  double dgxfac = 0.;
  double upperBoundEpsilon = DBL_EPSILON;
  /* System generated locals */
  int i_1, i_2;
  double d_1, d_2;

  /* Local variables */
  int diag, emin, emax;
  double root;
  int diag0;
  double epso6;
  int i_, j, k, q;
  double t;
  int x;
  double delta;
  int kalim, k1, q0;
  double t1;
  int x0;
  double lk, si, sk, wi, psifac;
  int dggdmx;
  double sw, oldphi;
  double phimin, phimax;
  int phipin;
  int dstsav, im1, lk0;
  int uk0;
  bool restrt;
  double twopsi, aki, akk, rad;
  int inc, irc;


  //  *** compute goldfeld-quandt-trotter step by more-hebden technique *** 
  //
  //  ***  (nl2sol version 2.2)  *** 

  //  ***  parameter declarations  *** 

  //     dimension dihdi(p*(p+1)/2), l(p*(p+1)/2), w(4*p+7) 

  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


  //  ***  purpose  *** 
  //     this routine is called as part of the nl2sol (nonlinear least- 
  //     squares) package (ref. 1), but it could be used in solving any 
  //     unconstrained minimization problem. 

  //  ***  algorithm notes  *** 

  //     The desired Goldfeld-Quandt-Trotter step (ref. 2, 3, 4) satisfies 
  //     (h + alpha*d**2)*step = -g  for some nonnegative alpha such that 
  //     h + alpha*d**2 is positive semidefinite.  alpha and step are 
  //     computed by a scheme analogous to the one described in ref. 5. 
  //     estimates of the smallest and largest eigenvalues of the hessian 
  //     are obtained from the gerschgorin circle theorem enhanced by a 
  //     simple form of the scaling described in ref. 6.  cases in which 
  //     h + alpha*d**2 is nearly (or exactly) singular are handled by 
  //     the technique discussed in ref. 2.  in these cases, a step of 
  //     (exact) length v(radius) is returned for which psi(step) exceeds 
  //     its optimal value by less than -v(epslon)*psi(step). 

  //  ***  functions and subroutines called  *** 

  // dotprd - returns inner product of two vectors. 
  // litvmu - applies inverse-transpose of compact lower triang. matrix. 
  // livmul - applies inverse of compact lower triang. matrix. 
  // lsqrt  - finds cholesky factor (of compactly stored lower triang.). 
  // lsvmin - returns approx. to min. sing. value of lower triang. matrix. 
  // rmdcon - returns machine-dependent constants. 
  // v2norm - returns 2-norm of a vector. 

  //  ***  references  *** 

  // 1.  dennis, j.e., gay, d.m., and welsch, r.e. (1981), an adaptive 
  //             nonlinear least-squares algorithm, acm trans. math. 
  //             software, vol. 7, no. 3. 
  // 2.  gay, d.m. (1981), computing optimal locally constrained steps, 
  //             siam j. sci. statist. computing, vol. 2, no. 2, pp. 
  //             186-197. 
  // 3.  goldfeld, s.m., quandt, r.e., and trotter, h.f. (1966), 
  //             maximization by quadratic hill-climbing, econometrica 34, 
  //             pp. 541-551. 
  // 4.  hebden, m.d. (1973), an algorithm for minimization using exact 
  //             second derivatives, report t.p. 515, theoretical physics 
  //             div., a.e.r.e. harwell, oxon., england. 
  // 5.  more, j.j. (1978), the levenberg-marquardt algorithm, implemen- 
  //             tation and theory, pp.105-116 of springer lecture notes 
  //             in mathematics no. 630, edited by g.a. watson, springer- 
  //             verlag, berlin and new york. 
  // 6.  varga, r.s. (1965), minimal gerschgorin sets, pacific j. math. 15, 
  //             pp. 719-729. 

  //  ***  general  *** 

  //     coded by david m. gay. 
  //     this subroutine was written in connection with research 
  //     supported by the national science foundation under grants 
  //     mcs-7600324, dcr75-10143, 76-1431ss, mcs76-11989, and 
  //     mcs-7906671. 

  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  /* Parameter adjustments */
  --dihdi;
  --l;
  --step;
  --dig;
  --d_;
  --v;
  --w;

  /* Function Body */

  //  ***  body  *** 

  //     ***  store largest abs. entry in (d**-1)*h*(d**-1) at w(dggdmx). 
  dggdmx = *p + 1;
  //     ***  store gerschgorin over- and underestimates of the largest 
  //     ***  and smallest eigenvalues of (d**-1)*h*(d**-1) at w(emax) 
  //     ***  and w(emin) respectively. 
  emax = dggdmx + 1;
  emin = emax + 1;
  //     ***  for use in recomputing step, the final values of lk, uk, dst, 
  //     ***  and the inverse derivative of more*s phi at 0 (for pos. def. 
  //     ***  h) are stored in w(lk0), w(uk0), w(dstsav), and w(phipin) 
  //     ***  respectively. 
  lk0 = emin + 1;
  phipin = lk0 + 1;
  uk0 = phipin + 1;
  dstsav = uk0 + 1;
  //     ***  store diag of (d**-1)*h*(d**-1) in w(diag),...,w(diag0+p). 
  diag0 = dstsav;
  diag = diag0 + 1;
  //     ***  store -d*step in w(q),...,w(q0+p). 
  q0 = diag0 + *p;
  q = q0 + 1;
  rad = v[radius];
  //     ***  phitol = max. error allowed in dst = v(dstnrm) = 2-norm of 
  //     ***  d*step. 
  phimax = v[phmxfc] * rad;
  phimin = v[phmnfc] * rad;
  //     ***  epso6 and psifac are used in checking for the special case 
  //     ***  of (nearly) singular h + alpha*d**2 (see ref. 2). 

  // Computing 2nd power 
  d_1 = rad;
  psifac = two * v[epslon] / (three * (four * (v[phmnfc] + one) * (kappa +
                                                                   one) + kappa + two) * (d_1 * d_1));
  //     ***  oldphi is used to detect limits of numerical accuracy.  if 
  //     ***  we recompute step and it does not change, then we accept it. 

  oldphi = zero;
  epso6 = v[epslon] / six;
  irc = 0;
  restrt = false;
  kalim = *ka + 50;

  //  ***  start or restart, depending on ka  *** 

  if (*ka >= 0) {
    goto L310;
  }

  //  ***  fresh start  *** 

  k = 0;
  uk = negone;
  *ka = 0;
  kalim = 50;

  //     ***  store diag(dihdi) in w(diag0+1),...,w(diag0+p)  *** 

  j = 0;
  i_1 = *p;
  for (i_ = 1; i_ <= i_1; ++i_) {
    j += i_;
    k1 = diag0 + i_;
    w[k1] = dihdi[j];
  }

  //     ***  determine w(dggdmx), the largest element of dihdi  *** 

  t1 = zero;
  j = *p * (*p + 1) / 2;
  i_1 = j;
  for (i_ = 1; i_ <= i_1; ++i_) {
    //
    //        t = (d_1 = dihdi[i_], abs(d_1));
    //
    t = fabs(dihdi[i_]);
    if (t1 < t) {
      t1 = t;
    }
  }
  w[dggdmx] = t1;

  //  ***  try alpha = 0  *** 

L40:
  lsqrt(1, p, &l[1], &dihdi[1], &irc);
  if (irc == 0) {
    goto L60;
  }
  //        ***  indef. h -- underestimate smallest eigenvalue, use this 
  //        ***  estimate to initialize lower bound lk on alpha. 
  j = irc * (irc + 1) / 2;
  t = l[j];
  l[j] = one;
  i_1 = irc;
  for (i_ = 1; i_ <= i_1; ++i_) {
    w[i_] = zero;
  }
  w[irc] = one;
  litvmu(&irc, &w[1], &l[1], &w[1]);
  t1 = v2norm(&irc, &w[1]);
  lk = -t / t1 / t1;
  v[dst0] = -lk;
  if (restrt) {
    goto L210;
  }
  v[nreduc] = zero;
  goto L70;

  //     ***  positive definite h -- compute unmodified newton step.  *** 
L60:
  lk = zero;
  livmul(p, &w[q], &l[1], &dig[1]);
  v[nreduc] = half * dotprd(p, &w[q], &w[q]);
  litvmu(p, &w[q], &l[1], &w[q]);
  dst = v2norm(p, &w[q]);
  v[dst0] = dst;
  phi = dst - rad;
  if (phi <= phimax) {
    goto L280;
  }
  if (restrt) {
    goto L210;
  }

  //  ***  prepare to compute gerschgorin estimates of largest (and 
  //  ***  smallest) eigenvalues.  *** 

L70:
  v[dgnorm] = v2norm(p, &dig[1]);
  if (v[dgnorm] == zero) {
    goto L450;
  }
  k = 0;
  i_1 = *p;
  for (i_ = 1; i_ <= i_1; ++i_) { //set w[j] = \Sigma_{j \neq i} |a_ij| for 1<=j<=*p, making use of the symmetry of A = (d**-1)*h*(d**-1)
    wi = zero;
    if (i_ == 1) {
      goto L90;
    }
    im1 = i_ - 1;
    i_2 = im1;
    for (j = 1; j <= i_2; ++j) {
      ++k;
      //
      //            t = (d_1 = dihdi[k], abs(d_1));
      //
      t = fabs(dihdi[k]);
      wi += t;
      w[j] += t;
    }
  L90:
    w[i_] = wi;
    ++k;
  }

  //  ***  (under-)estimate smallest eigenvalue of (d**-1)*h*(d**-1)  *** 

  k = 1;
  t1 = w[diag] - w[1];
  if (*p <= 1) {
    goto L120;
  }
  i_1 = *p;
  for (i_ = 2; i_ <= i_1; ++i_) { //get index k of the row minimising t = a_ii - \Sigma |a_ij|
    j = diag0 + i_;
    t = w[j] - w[i_];
    if (t >= t1) {
      goto L110;
    }
    t1 = t;
    k = i_;
  L110:
    ;
  }

L120:
  sk = w[k];
  j = diag0 + k;
  akk = w[j];               // i.e a_kk
  k1 = k * (k - 1) / 2 + 1; // index of first element of kth row
  inc = 1;
  t = zero;
  i_1 = *p;
  for (i_ = 1; i_ <= i_1; ++i_) {
    if (i_ == k) {
      goto L130;
    }
    //
    //        aki = (d_1 = dihdi[k1], abs(d_1));
    //
    aki = fabs(dihdi[k1]);
    si = w[i_];
    j = diag0 + i_;
    t1 = half * (akk - w[j] + si - aki);
    t1 += sqrt(t1 * t1 + sk * aki);
    if (t < t1) {
      t = t1;
    }
    if (i_ < k) {
      goto L140;
    }
  L130:
    inc = i_;
  L140:
    k1 += inc;
  }

  w[emin] = akk - t;
  uk = v[dgnorm] / rad - w[emin];
  // Due to round-off / subtractive cancellation, it is indeed possible that 'uk' is now negative. However, since 'uk' is an estimate for an upper bound
  // of alpha, and since alpha is intended to be a nonnegative number (see 'algorithm notes' above), it is necessary to ensure that uk cannot be negative.
  if (uk < 0)
    uk = 0;

  //  ***  compute gerschgorin (over-)estimate of largest eigenvalue  *** 

  k = 1;
  t1 = w[diag] + w[1];
  if (*p <= 1) {
    goto L170;
  }
  i_1 = *p;
  for (i_ = 2; i_ <= i_1; ++i_) {
    j = diag0 + i_;
    t = w[j] + w[i_];
    if (t <= t1) {
      goto L160;
    }
    t1 = t;
    k = i_;
  L160:
    ;
  }

L170:
  sk = w[k];
  j = diag0 + k;
  akk = w[j];
  k1 = k * (k - 1) / 2 + 1;
  inc = 1;
  t = zero;
  i_1 = *p;
  for (i_ = 1; i_ <= i_1; ++i_) {
    if (i_ == k) {
      goto L180;
    }
    //
    //        aki = (d_1 = dihdi[k1], abs(d_1));
    //
    aki = fabs(dihdi[k1]);
    si = w[i_];
    j = diag0 + i_;
    t1 = half * (w[j] + si - aki - akk);
    t1 += sqrt(t1 * t1 + sk * aki);
    if (t < t1) {
      t = t1;
    }
    if (i_ < k) {
      goto L190;
    }
  L180:
    inc = i_;
  L190:
    k1 += inc;
  }

  w[emax] = akk + t;
  // Computing MAX 
  //
  //    d_1 = lk, d_2 = v[dgnorm] / rad - w[emax];
  //    lk = std::max(d_1,d_2);
  lk = std::max(lk, v[dgnorm] / rad - w[emax]);
  // Due to round-off / subtractive cancellation, it is indeed possible that 'lk' is now negative, which is not admissible.
  if (lk < 0)
    lk = 0;

  //     ***  alphak = current value of alpha (see alg. notes above).  we 
  //     ***  use more*s scheme for initializing it. 
  //
  //    alphak = (d_1 = v[stppar], abs(d_1)) * v[rad0] / rad;
  //
  alphak = fabs(v[stppar]) * v[rad0] / rad;
  ASSERT(isFinite(alphak));

  if (irc != 0) {
    goto L210;
  }

  //  ***  compute l0 for positive definite h  *** 

  livmul(p, &w[1], &l[1], &w[q]);
  t = v2norm(p, &w[1]);
  w[phipin] = dst / t / t;
  // Computing MAX 
  //
  //    d_1 = lk, d_2 = phi * w[phipin];
  //    lk = std::max(d_1,d_2);
  lk = std::max(lk, phi * w[phipin]);
  // Due to round-off / subtractive cancellation, it is indeed possible that 'lk' is now negative, which is not admissible.
  if (lk < 0)
    lk = 0;

  //  ***  safeguard alphak and add alphak*i to (d**-1)*h*(d**-1)  *** 

L210:
  ++(*ka);
  if (-v[dst0] >= alphak || alphak < lk || alphak >= uk) { // if H is not positive definite, then -v[dst0]=lk. Otherwise v[dst0]=||d h**-1 g||
    // NOTE (pj): When lk and uk are within (almost) machine precision of each other, the square root of their ratio
    // is essentially just 1 plus noise. In that case, we use 1 instead.
    alphak = std::max(lk, uk * std::max(0.001, uk <= lk * (1 + 2 * DBL_EPSILON) ? 1 : sqrt(lk / uk)));
    ASSERT(isFinite(alphak));
    // NOTE (pj): The Goldfeldt-Quandt-Trotter step searches for an alpha such that (H + alpha · I) is a positive 
    // definite matrix with minimum alpha (it may also be replaced by a diagonal scale matrix). The algorithm works by
    // finding alpha between the upper bound (uk) and the lower bound (lk). In some cases, uk is, in fact, NOT a true
    // upper bound due to these bounds actually being only estimates of the true bounds. The algorithm fails
    // when the lower bound becomes equal to the upper bound. In that situation, we lift the upper bound. We gently
    // grow the growth factor, initially starting with twice the relative machine resolution DBL_EPSILON.
    if (uk <= lk * (1 + DBL_EPSILON)) {
      upperBoundEpsilon *= 2;
      if (upperBoundEpsilon >= 1.0)
        throw std::runtime_error("NL2SOL::gqtstp- the adjustment to the upper bound of the Goldfeld-Quandt-Trotter step to overcome numerical trunctation issues has become large enough to suggest that a different type of error has arisen.");
      uk = std::max(uk, lk) * (1. + upperBoundEpsilon);
    }
  }
  k = 0;
  i_1 = *p;
  for (i_ = 1; i_ <= i_1; ++i_) {
    k += i_;
    j = diag0 + i_;
    dihdi[k] = w[j] + alphak;
  }

  //  ***  try computing cholesky decomposition  *** 

  lsqrt(1, p, &l[1], &dihdi[1], &irc);
  if (irc == 0) {
    goto L250;
  }

  //  ***  (d**-1)*h*(d**-1) + alphak*i  is indefinite -- overestimate
  //  ***  smallest eigenvalue for use in updating lk  *** 

  j = irc * (irc + 1) / 2;
  t = l[j];
  l[j] = one;
  i_1 = irc;
  for (i_ = 1; i_ <= i_1; ++i_) {
    w[i_] = zero;
  }
  w[irc] = one;
  litvmu(&irc, &w[1], &l[1], &w[1]);
  t1 = v2norm(&irc, &w[1]);
  lk = alphak - t / t1 / t1;
  // Due to round-off / subtractive cancellation, it is indeed possible that 'lk' is now negative, which is not admissible.
  if (lk < 0)
    lk = 0;
  v[dst0] = -lk;

  goto L210;

  //  ***  alphak makes (d**-1)*h*(d**-1) positive definite. 
  //  ***  compute q = -d*step, check for convergence.  *** 

L250:
  livmul(p, &w[q], &l[1], &dig[1]);
  litvmu(p, &w[q], &l[1], &w[q]);
  dst = v2norm(p, &w[q]);
  phi = dst - rad;
  if (phi <= phimax && phi >= phimin) {
    goto L290;
  }
  if (phi == oldphi) {
    goto L290;
  }
  oldphi = phi;
  if (phi > zero) {
    goto L260;
  }
  //        ***  check for the special case of  h + alpha*d**2  (nearly) 
  //        ***  singular.  delta is .ge. the smallest eigenvalue of 
  //        ***  (d**-1)*h*(d**-1) + alphak*i. 
  if (v[dst0] > zero) {
    goto L260;
  }
  delta = alphak + v[dst0];
  twopsi = alphak * dst * dst + dotprd(p, &dig[1], &w[q]);
  if (delta < psifac * twopsi) {
    goto L270;
  }

  //  ***  unacceptable alphak -- update lk, uk, alphak  *** 

L260:
  if (*ka >= kalim) {
    goto L290;
  }
  livmul(p, &w[1], &l[1], &w[q]);
  t1 = v2norm(p, &w[1]);
  //     ***  the following dmin1 is necessary because of restarts  *** 
  if (phi < zero) {
    uk = std::min(uk, alphak);
    // Due to round-off / subtractive cancellation, it is indeed possible that 'uk' is now negative. However, since 'uk' is an estimate for an upper bound
    // of alpha, and since alpha is intended to be a nonnegative number (see 'algorithm notes' above), it is necessary to ensure that uk cannot be negative.
    if (uk < 0)
      uk = 0;
  }
  alphak += phi / t1 * (dst / t1) * (dst / rad);
  ASSERT(isFinite(alphak));
  lk = std::max(lk, alphak);
  // Due to round-off / subtractive cancellation, it is indeed possible that 'lk' is now negative, which is not admissible.
  if (lk < 0)
    lk = 0;
  upperBoundEpsilon = DBL_EPSILON;
  goto L210;

  //  ***  decide how to handle (nearly) singular h + alpha*d**2  *** 

  //     ***  if not yet available, obtain machine dependent value dgxfac. 

L270:
  if (dgxfac == zero) {
    dgxfac = epsfac * rmdcon(3);
  }

  //     ***  now decide.  *** 
  if (delta > dgxfac * w[dggdmx]) {
    goto L350;
  }
  //        ***  delta is so small we cannot handle the special case in 
  //        ***  the available arithmetic.  accept step as it is. 
  goto L290;

  //  ***  acceptable step on first try  *** 

L280:
  alphak = zero;
  ASSERT(isFinite(alphak));

  //  ***  successful step in general.  compute step = -(d**-1)*q  *** 

L290:
  i_1 = *p;
  for (i_ = 1; i_ <= i_1; ++i_) {
    j = q0 + i_;
    step[i_] = -w[j] / d_[i_];
  }
  v[gtstep] = -dotprd(p, &dig[1], &w[q]);
  v[preduc] = half * (fabs(alphak) * dst * dst - v[gtstep]);
  goto L430;


  //  ***  restart with new radius  *** 

L310:
  if (v[dst0] <= zero || v[dst0] - rad > phimax) {
    goto L330;
  }

  //     ***  prepare to return newton step  *** 

  restrt = true;
  ++(*ka);
  k = 0;
  i_1 = *p;
  for (i_ = 1; i_ <= i_1; ++i_) {
    k += i_;
    j = diag0 + i_;
    dihdi[k] = w[j];
  }
  uk = negone;
  goto L40;

L330:
  if (*ka == 0) {
    goto L60;
  }

  dst = w[dstsav];
  //
  //    alphak = (d_1 = v[stppar], abs(d_1));
  //
  alphak = fabs(v[stppar]);
  ASSERT(isFinite(alphak));
  phi = dst - rad;
  t = v[dgnorm] / rad;
  if (rad > v[rad0]) {
    goto L340;
  }

  //        ***  smaller radius  *** 
  uk = t - w[emin];
  // Due to round-off / subtractive cancellation, it is indeed possible that 'uk' is now negative. However, since 'uk' is an estimate for an upper bound
  // of alpha, and since alpha is intended to be a nonnegative number (see 'algorithm notes' above), it is necessary to ensure that uk cannot be negative.
  if (uk < 0)
    uk = 0;
  lk = zero;
  if (alphak > zero) {
    lk = w[lk0];
    // Due to round-off / subtractive cancellation, it is indeed possible that 'lk' is now negative, which is not admissible.
    if (lk < 0)
      lk = 0;
  }
  // Computing MAX 
  //
  //    d_1 = lk, d_2 = t - w[emax];
  //    lk = std::max(d_1,d_2);
  lk = std::max(lk, t - w[emax]);
  if (v[dst0] > zero) {
    // Computing MAX 
    //
    //        d_1 = lk, d_2 = (v[dst0] - rad) * w[phipin];
    //        lk = std::max(d_1,d_2);
    lk = std::max(lk, (v[dst0] - rad) * w[phipin]);
  }
  // Due to round-off / subtractive cancellation, it is indeed possible that 'lk' is now negative, which is not admissible.
  if (lk < 0)
    lk = 0;
  goto L260;

  //     ***  bigger radius  *** 
L340:
  uk = t - w[emin];
  if (alphak > zero) {
    // Computing MIN 
    d_1 = uk, d_2 = w[uk0];
    uk = std::min(d_1, d_2);
  }
  // Due to round-off / subtractive cancellation, it is indeed possible that 'uk' is now negative. However, since 'uk' is an estimate for an upper bound
  // of alpha, and since alpha is intended to be a nonnegative number (see 'algorithm notes' above), it is necessary to ensure that uk cannot be negative.
  if (uk < 0)
    uk = 0;
  // Computing MAX 
  d_1 = zero, d_2 = -v[dst0], d_1 = std::max(d_1, d_2), d_2 = t - w[emax];
  lk = std::max(d_1, d_2);
  if (v[dst0] > zero) {
    // Computing MAX 
    d_1 = lk, d_2 = (v[dst0] - rad) * w[phipin];
    lk = std::max(d_1, d_2);
  }
  // Due to round-off / subtractive cancellation, it is indeed possible that 'lk' is now negative, which is not admissible.
  if (lk < 0)
    lk = 0;
  goto L260;

  //  ***  handle (nearly) singular h + alpha*d**2  *** 

  //     ***  negate alphak to indicate special case  *** 
L350:
  alphak = -alphak;
  ASSERT(isFinite(alphak));
  //     ***  allocate storage for scratch vector x  *** 
  x0 = q0 + *p;
  x = x0 + 1;

  //  ***  use inverse power method with start from lsvmin to obtain 
  //  ***  approximate eigenvector corresponding to smallest eigenvalue 
  //  ***  of (d**-1)*h*(d**-1). 

  delta = kappa * delta;
  t = lsvmin(p, &l[1], &w[x], &w[1]);

  k = 0;
  //     ***  normalize w  *** 
L360:
  i_1 = *p;
  for (i_ = 1; i_ <= i_1; ++i_) {
    w[i_] = t * w[i_];
  }
  //     ***  complete current inv. power iter. -- replace w by (l**-t)*w. 
  litvmu(p, &w[1], &l[1], &w[1]);
  t1 = one / v2norm(p, &w[1]);
  t = t1 * t;
  if (t <= delta) {
    goto L390;
  }
  if (k > 30) {
    goto L290;
  }
  ++k;
  //     ***  start next inv. power iter. by storing normalized w in x. 
  i_1 = *p;
  for (i_ = 1; i_ <= i_1; ++i_) {
    j = x0 + i_;
    w[j] = t1 * w[i_];
  }
  //     ***  compute w = (l**-1)*x. 
  livmul(p, &w[1], &l[1], &w[x]);
  t = one / v2norm(p, &w[1]);
  goto L360;

L390:
  i_1 = *p;
  for (i_ = 1; i_ <= i_1; ++i_) {
    w[i_] = t1 * w[i_];
  }

  //  ***  now w is the desired approximate (unit) eigenvector and 
  //  ***  t*x = ((d**-1)*h*(d**-1) + alphak*i)*w. 

  sw = dotprd(p, &w[q], &w[1]);
  t1 = (rad + dst) * (rad - dst);
  root = sqrt(sw * sw + t1);
  if (sw < zero) {
    root = -root;
  }
  si = t1 / (sw + root);
  //     ***  accept current step if adding si*w would lead to a 
  //     ***  further relative reduction in psi of less than v(epslon)/3. 
  v[preduc] = half * twopsi;
  t1 = zero;
  t = si * (alphak * sw - half * si * (alphak + t * dotprd(p, &w[x], &w[1])));
  if (t < epso6 * twopsi) {
    goto L410;
  }
  v[preduc] += t;
  dst = rad;
  t1 = -si;
L410:
  i_1 = *p;
  for (i_ = 1; i_ <= i_1; ++i_) {
    j = q0 + i_;
    w[j] = t1 * w[i_] - w[j];
    step[i_] = w[j] / d_[i_];
  }
  v[gtstep] = dotprd(p, &dig[1], &w[q]);

  //  ***  save values for use in a possible restart  *** 

L430:
  v[dstnrm] = dst;
  v[stppar] = alphak;
  w[lk0] = lk;
  w[uk0] = uk;
  v[rad0] = rad;
  w[dstsav] = dst;

  //     ***  restore diagonal of dihdi  *** 

  j = 0;
  i_1 = *p;
  for (i_ = 1; i_ <= i_1; ++i_) {
    j += i_;
    k = diag0 + i_;
    dihdi[j] = w[k];
  }
  goto L999;

  //  ***  special case -- g = 0  *** 

L450:
  v[stppar] = zero;
  v[preduc] = zero;
  v[dstnrm] = zero;
  v[gtstep] = zero;
  i_1 = *p;
  for (i_ = 1; i_ <= i_1; ++i_) {
    step[i_] = zero;
  }

L999:
  return 0;

} /* gqtstp */

int litvmu(int *n, double *x, double *l, double *y) {
  /* Initialized data */

  const double zero = 0.;

  /* System generated locals */
  int i_1, i_2;

  /* Local variables */
  int i_, j, i0, ii, ij;
  double xi;
  int im1, np1;


  //  ***  solve  (l**t)*x = y,  where  l  is an  n x n  lower triangular 
  //  ***  matrix stored compactly by rows.  x and y may occupy the same 
  //  ***  storage.  *** 

  /* Parameter adjustments */
  --y;
  --x;
  --l;

  /* Function Body */

  i_1 = *n;
  for (i_ = 1; i_ <= i_1; ++i_) {
    ASSERT(isFinite(y[i_]));
    x[i_] = y[i_];
  }
  np1 = *n + 1;
  i0 = *n * (*n + 1) / 2;
  i_1 = *n;
  for (ii = 1; ii <= i_1; ++ii) {
    i_ = np1 - ii;
    xi = x[i_] / l[i0];
    ASSERT(isFinite(xi));
    x[i_] = xi;
    if (i_ <= 1) {
      goto L999;
    }
    i0 -= i_;
    if (xi == zero) {
      goto L30;
    }
    im1 = i_ - 1;
    i_2 = im1;
    for (j = 1; j <= i_2; ++j) {
      ij = i0 + j;
      ASSERT(isFinite(xi * l[ij]));
      x[j] -= xi * l[ij];
    }
  L30:
    ;
  }
L999:
  return 0;
} /* litvmu */

int livmul(int *n, double *x, double *l, double *y) {
  /* Initialized data */

  const double zero = 0.;

  /* System generated locals */
  int i_1, i_2;

  /* Local variables */
  int i_, j, k;
  double t;

  //  ***  solve  l*x = y, where  l  is an  n x n  lower triangular 
  //  ***  matrix stored compactly by rows.  x and y may occupy the same 
  //  ***  storage.  *** 

  /* Parameter adjustments */
  --y;
  --x;
  --l;

  /* Function Body */

  i_1 = *n;
  for (k = 1; k <= i_1; ++k) {
    if (y[k] != zero) {
      goto L20;
    }
    x[k] = zero;
  }
  goto L999;
L20:
  j = k * (k + 1) / 2;
  x[k] = y[k] / l[j];
  ASSERT(isFinite(x[k]));
  if (k >= *n) {
    goto L999;
  }
  ++k;
  i_1 = *n;
  for (i_ = k; i_ <= i_1; ++i_) {
    i_2 = i_ - 1;
    t = dotprd(&i_2, &l[j + 1], &x[1]);
    j += i_;
    x[i_] = (y[i_] - t) / l[j];
    ASSERT(isFinite(x[i_]));
  }
L999:
  return 0;
} /* livmul */

int lmstep(double* d_, double* g, int* ierr, int* ipivot, int* ka, int* p, double* qtr, double* r_, double* step, double* v, double* w) {
  /* Initialized data */

  double alphak = 0.;
  double psifac = 0.;
  const int dgnorm = 1;
  const int dstnrm = 2;
  const int dst0 = 3;
  const int epslon = 19;
  const int gtstep = 4;
  const int nreduc = 6;
  const int phmnfc = 20;
  const int phmxfc = 21;
  const int preduc = 7;
  const int radius = 8;
  const int rad0 = 9;
  const int stppar = 5;
  const double dfac = 256.;
  const double eight = 8.;
  const double half = .5;
  const double negone = -1.;
  const double one = 1.;
  const double p001 = .001;
  const double three = 3.;
  const double ttol = 2.5;
  const double zero = 0.;

  /* System generated locals */
  int i_1, i_2, i_3;
  double d_1, d_2;

  /* Local variables */
  int pp1o2, rmat;
  double dtol;
  int rmat0;
  double a, b;
  int i_, k, l;
  double t;
  int kalim;
  double d1, d2;
  int i1, j_1;
  double lk, si, sj, dfacsq, uk, wl, oldphi, phimin, phimax;
  int phipin;
  int dstsav;
  double sqrtak;
  int lk0;
  int ip1, uk0;
  double twopsi, adi, rad, phi;
  int res;
  double dst;
  int res0;


  //  ***  compute levenberg-marquardt step using more-hebden technique  ** 
  //  ***  nl2sol version 2.2.  *** 

  //  ***  parameter declarations  *** 

  //     dimension w(p*(p+5)/2 + 4) 

  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  //  ***  purpose  *** 

  //        given the r matrix from the qr decomposition of a jacobian 
  //             nonlinear least-squares algorithm, acm trans. math. 
  //             software, vol. 7, no. 3. 
  // 2.  gay, d.m. (1981), computing optimal locally constrained steps, 
  //             siam j. sci. statist. computing, vol. 2, no. 2, pp. 
  //             186-197. 
  // 3.  lawson, c.l., and hanson, r.j. (1974), solving least squares 
  //             problems, prentice-hall, englewood cliffs, n.j. 
  // 4.  more, j.j. (1978), the levenberg-marquardt algorithm, implemen- 
  //             tation and theory, pp.105-116 of springer lecture notes 
  //             in mathematics no. 630, edited by g.a. watson, springer- 
  //             verlag, berlin and new york. 

  //  ***  general  *** 

  //     coded by david m. gay. 
  //     this subroutine was written in connection with research 
  //     supported by the national science foundation under grants 
  //     mcs-7600324, dcr75-10143, 76-1431ss, mcs76-11989, and 
  //     mcs-7906671. 

  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  /* Parameter adjustments */
  --step;
  --qtr;
  --ipivot;
  --g;
  --d_;
  --r_;
  --v;
  --w;

  /* Function Body */

  //  ***  body  *** 

  //     ***  for use in recomputing step, the final values of lk and uk, 
  //     ***  the inverse derivative of more*s phi at 0 (for nonsing. j) 
  //     ***  and the value returned as v(dstnrm) are stored at w(lk0), 
  //     ***  w(uk0), w(phipin), and w(dstsav) respectively. 
  lk0 = *p + 1;
  phipin = lk0 + 1;
  uk0 = phipin + 1;
  dstsav = uk0 + 1;
  rmat0 = dstsav;
  //     ***  a copy of the r-matrix from the qr decomposition of j is 
  //     ***  stored in w starting at w(rmat), and a copy of the residual 
  //     ***  vector is stored in w starting at w(res).  the loops below 
  //     ***  that update the qr decomp. for a nonzero marquardt parameter 
  //     ***  work on these copies. 
  rmat = rmat0 + 1;
  pp1o2 = *p * (*p + 1) / 2;
  res0 = pp1o2 + rmat0;
  res = res0 + 1;
  rad = v[radius];
  if (rad > zero) {
    // Computing 2nd power 
    d_1 = rad;
    psifac = v[epslon] / ((eight * (v[phmnfc] + one) + three) * (d_1 * d_1));
  }
  phimax = v[phmxfc] * rad;
  phimin = v[phmnfc] * rad;
  //     ***  dtol, dfac, and dfacsq are used in rescaling the fast givens 
  //     ***  representation of the updated qr decomposition. 
  dtol = one / dfac;
  dfacsq = dfac * dfac;
  //     ***  oldphi is used to detect limits of numerical accuracy.  if 
  //     ***  we recompute step and it does not change, then we accept it. 
  oldphi = zero;
  lk = zero;
  uk = zero;
  kalim = *ka + 12;

  //  ***  start or restart, depending on ka  *** 

  if (*ka < 0) {
    goto L10;
  } else if (*ka == 0) {
    goto L20;
  } else {
    goto L370;
  }

  //  ***  fresh start -- compute v(nreduc)  *** 

L10:
  *ka = 0;
  kalim = 12;
  k = *p;
  if (*ierr != 0) {
    k = abs(*ierr) - 1;
  }
  v[nreduc] = half * dotprd(&k, &qtr[1], &qtr[1]);

  //  ***  set up to try initial gauss-newton step  *** 

L20:
  v[dst0] = negone;
  if (*ierr != 0) {
    goto L90;
  }

  //  ***  compute gauss-newton step  *** 

  //     ***  note -- the r-matrix is stored compactly by columns in 
  //     ***  r(1), r(2), r(3), ...  it is the transpose of a 
  //     ***  lower triangular matrix stored compactly by rows, and we 
  //     ***  treat it as such when using litvmu and livmul. 
  litvmu(p, &w[1], &r_[1], &qtr[1]);
  //     ***  temporarily store permuted -d*step in step. 
  i_1 = *p;
  for (i_ = 1; i_ <= i_1; ++i_) {
    j_1 = ipivot[i_];
    step[i_] = d_[j_1] * w[i_];
  }
  dst = v2norm(p, &step[1]);
  v[dst0] = dst;
  phi = dst - rad;
  if (phi <= phimax) {
    goto L410;
  }
  //     ***  if this is a restart, go to 110  *** 
  if (*ka > 0) {
    goto L110;
  }

  //  ***  gauss-newton step was unacceptable.  compute l0  *** 

  i_1 = *p;
  for (i_ = 1; i_ <= i_1; ++i_) {
    j_1 = ipivot[i_];
    step[i_] = d_[j_1] * (step[i_] / dst);
  }
  livmul(p, &step[1], &r_[1], &step[1]);
  t = one / v2norm(p, &step[1]);
  w[phipin] = t / dst * t;
  lk = phi * w[phipin];

  //  ***  compute u0  *** 

L90:
  i_1 = *p;
  for (i_ = 1; i_ <= i_1; ++i_) {
    w[i_] = g[i_] / d_[i_];
  }
  v[dgnorm] = v2norm(p, &w[1]);
  uk = v[dgnorm] / rad;
  if (uk <= zero) {
    goto L390;
  }

  //     ***  alphak will be used as the current marquardt parameter.  we 
  //     ***  use more*s scheme for initializing it. 
  //
  //    alphak = (d_1 = v[stppar], abs(d_1)) * v[rad0] / rad;
  //
  alphak = fabs(v[stppar]) * v[rad0] / rad;

  //  ***  top of loop -- increment ka, copy r to rmat, qtr to res  *** 

L110:
  ++(*ka);
  vcopy(pp1o2, &w[rmat], &r_[1]);
  vcopy(*p, &w[res], &qtr[1]);

  //  ***  safeguard alphak and initialize fast givens scale vector.  *** 

  if (alphak <= zero || alphak < lk || alphak >= uk) {
    // Computing MAX 
    d_1 = p001, d_2 = sqrt(lk / uk);
    alphak = uk * std::max(d_1, d_2);
    ASSERT(isFinite(alphak));
  }
  sqrtak = sqrt(alphak);
  ASSERT(isFinite(sqrtak));
  i_1 = *p;
  for (i_ = 1; i_ <= i_1; ++i_) {
    w[i_] = one;
  }

  //  ***  add alphak*d and update qr decomp. using fast givens trans.  *** 

  i_1 = *p;
  for (i_ = 1; i_ <= i_1; ++i_) {
    //        ***  generate, apply 1st givens trans. for row i of alphak*d
    //        ***  (use step to store temporary row)  *** 
    l = i_ * (i_ + 1) / 2 + rmat0;
    wl = w[l];
    d2 = one;
    d1 = w[i_];
    j_1 = ipivot[i_];
    adi = sqrtak * d_[j_1];
    ASSERT(isFinite(adi));
    if (adi >= fabs(wl)) {
      goto L150;
    }
  L130:
    a = adi / wl;
    ASSERT(isFinite(a));
    b = d2 * a / d1;
    ASSERT(isFinite(b));
    t = a * b + one;
    if (t > ttol) {
      goto L150;
    }
    w[i_] = d1 / t;
    d2 /= t;
    w[l] = t * wl;
    a = -a;
    i_2 = *p;
    for (j_1 = i_; j_1 <= i_2; ++j_1) {
      l += j_1;
      step[j_1] = a * w[l];
    }
    goto L170;

  L150:
    b = wl / adi;
    a = d1 * b / d2;
    t = a * b + one;
    if (t > ttol) {
      goto L130;
    }
    w[i_] = d2 / t;
    d2 = d1 / t;
    w[l] = t * adi;
    i_2 = *p;
    for (j_1 = i_; j_1 <= i_2; ++j_1) {
      l += j_1;
      wl = w[l];
      step[j_1] = -wl;
      w[l] = a * wl;
    }

  L170:
    if (i_ == *p) {
      goto L280;
    }

    //        ***  now use givens trans. to zero elements of temp. row  *

    ip1 = i_ + 1;
    i_2 = *p;
    for (i1 = ip1; i1 <= i_2; ++i1) {
      l = i1 * (i1 + 1) / 2 + rmat0;
      wl = w[l];
      ASSERT(isFinite(wl));
      si = step[i1 - 1];
      ASSERT(isFinite(si));
      d1 = w[i1];
      ASSERT(isFinite(d1));

      //             ***  rescale row i1 if necessary  *** 

      if (d1 >= dtol) {
        goto L190;
      }
      d1 *= dfacsq;
      wl /= dfac;
      k = l;
      i_3 = *p;
      for (j_1 = i1; j_1 <= i_3; ++j_1) {
        k += j_1;
        w[k] /= dfac;
      }

      //             ***  use givens trans. to zero next element of temp

    L190:
      if (fabs(si) > fabs(wl)) {
        goto L220;
      }
      if (si == zero) {
        goto L260;
      }
    L200:
      a = si / wl;
      ASSERT(isFinite(a));
      b = d2 * a / d1;
      ASSERT(isFinite(b));
      t = a * b + one;
      if (t > ttol) {
        goto L220;
      }
      w[l] = t * wl;
      w[i1] = d1 / t;
      d2 /= t;
      i_3 = *p;
      for (j_1 = i1; j_1 <= i_3; ++j_1) {
        l += j_1;
        wl = w[l];
        sj = step[j_1];
        w[l] = wl + b * sj;
        step[j_1] = sj - a * wl;
      }
      goto L240;

    L220:
      b = wl / si;
      a = d1 * b / d2;
      t = a * b + one;
      if (t > ttol) {
        goto L200;
      }
      w[i1] = d2 / t;
      d2 = d1 / t;
      w[l] = t * si;
      i_3 = *p;
      for (j_1 = i1; j_1 <= i_3; ++j_1) {
        l += j_1;
        wl = w[l];
        sj = step[j_1];
        w[l] = a * wl + sj;
        step[j_1] = b * sj - wl;
      }

      //             ***  rescale temp. row if necessary  *** 

    L240:
      if (d2 >= dtol) {
        goto L260;
      }
      d2 *= dfacsq;
      i_3 = *p;
      for (k = i1; k <= i_3; ++k) {
        step[k] /= dfac;
      }
    L260:
      ;
    }
  }

  //  ***  compute step  *** 

L280:
  litvmu(p, &w[res], &w[rmat], &w[res]);
  //     ***  recover step and store permuted -d*step at w(res)  *** 
  i_1 = *p;
  for (i_ = 1; i_ <= i_1; ++i_) {
    j_1 = ipivot[i_];
    k = res0 + i_;
    t = w[k];
    step[j_1] = -t;
    w[k] = t * d_[j_1];
  }
  dst = v2norm(p, &w[res]);
  phi = dst - rad;
  if (phi <= phimax && phi >= phimin) {
    goto L430;
  }
  if (oldphi == phi) {
    goto L430;
  }
  oldphi = phi;

  //  ***  check for (and handle) special case  *** 

  if (phi > zero) {
    goto L310;
  }
  if (*ka >= kalim) {
    goto L430;
  }
  twopsi = alphak * dst * dst - dotprd(p, &step[1], &g[1]);
  if (alphak >= twopsi * psifac) {
    goto L310;
  }
  v[stppar] = -alphak;
  goto L440;

  //  ***  unacceptable step -- update lk, uk, alphak, and try again  *** 

L300:
  if (phi < zero) {
    uk = std::min(uk, alphak);
  }
  goto L320;
L310:
  if (phi < zero) {
    uk = alphak;
  }
L320:
  i_1 = *p;
  for (i_ = 1; i_ <= i_1; ++i_) {
    j_1 = ipivot[i_];
    k = res0 + i_;
    step[i_] = d_[j_1] * (w[k] / dst);
  }
  livmul(p, &step[1], &w[rmat], &step[1]);
  i_1 = *p;
  for (i_ = 1; i_ <= i_1; ++i_) {
    step[i_] /= sqrt(w[i_]);
  }
  t = one / v2norm(p, &step[1]);
  alphak += t * phi * t / rad;
  ASSERT(isFinite(alphak));
  lk = std::max(lk, alphak);
  goto L110;

  //  ***  restart  *** 

L370:
  lk = w[lk0];
  uk = w[uk0];
  if (v[dst0] > zero && v[dst0] - rad <= phimax) {
    goto L20;
  }
  //
  //    alphak = (d_1 = v[stppar], abs(d_1));
  //
  alphak = fabs(v[stppar]);
  ASSERT(isFinite(alphak));
  dst = w[dstsav];
  phi = dst - rad;
  t = v[dgnorm] / rad;
  if (rad > v[rad0]) {
    goto L380;
  }

  //        ***  smaller radius  *** 
  uk = t;
  if (alphak <= zero) {
    lk = zero;
  }
  if (v[dst0] > zero) {
    // Computing MAX 
    d_1 = lk, d_2 = (v[dst0] - rad) * w[phipin];
    lk = std::max(d_1, d_2);
  }
  goto L300;

  //     ***  bigger radius  *** 
L380:
  if (alphak <= zero || uk > t) {
    uk = t;
  }
  lk = zero;
  if (v[dst0] > zero) {
    // Computing MAX 
    d_1 = lk, d_2 = (v[dst0] - rad) * w[phipin];
    lk = std::max(d_1, d_2);
  }
  goto L300;

  //  ***  special case -- rad .le. 0 or (g = 0 and j is singular)  *** 

L390:
  v[stppar] = zero;
  dst = zero;
  lk = zero;
  uk = zero;
  v[gtstep] = zero;
  v[preduc] = zero;
  i_1 = *p;
  for (i_ = 1; i_ <= i_1; ++i_) {
    step[i_] = zero;
  }
  goto L450;

  //  ***  acceptable gauss-newton step -- recover step from w  *** 

L410:
  alphak = zero;
  i_1 = *p;
  for (i_ = 1; i_ <= i_1; ++i_) {
    j_1 = ipivot[i_];
    step[j_1] = -w[i_];
  }

  //  ***  save values for use in a possible restart  *** 

L430:
  v[stppar] = alphak;
L440:
  v[gtstep] = dotprd(p, &step[1], &g[1]);
  v[preduc] = half * (alphak * dst * dst - v[gtstep]);
L450:
  v[dstnrm] = dst;
  w[dstsav] = dst;
  w[lk0] = lk;
  w[uk0] = uk;
  v[rad0] = rad;

  return 0;

} /* lmstep */

double lsvmin(int *p, double *l, double *x, double *y) {
  /* Initialized data */

  const double half = .5;
  const double one = 1.;
  const double r9973 = 9973.;
  const double zero = 0.;
  int ix = 2;

  /* System generated locals */
  int i_1, i_2;
  double ret_val;

  /* Local variables */
  double b;
  int i_, j;
  double t;
  int j_0;
  double splus, xplus;
  int ii, ji, jj, pplus1;
  int jm1;
  double sminus, xminus;
  int jjj;
  double psj;

  //  ***  estimate smallest sing. value of packed lower triang. matrix l 

  //  ***  parameter declarations  *** 

  //     dimension l(p*(p+1)/2) 

  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  //  ***  purpose  *** 

  //     this function returns a good over-estimate of the smallest 
  //     singular value of the packed lower triangular matrix l. 

  //  ***  parameter description  *** 

  //  p (in)  = the order of l.  l is a  p x p  lower triangular matrix. 
  //  l (in)  = array holding the elements of  l  in row order, i.e. 
  //             l(1,1), l(2,1), l(2,2), l(3,1), l(3,2), l(3,3), etc. 
  //  x (out) if lsvmin returns a positive value, then x is a normalized 
  //             approximate left singular vector corresponding to the 
  //             smallest singular value.  this approximation may be very 
  //             crude.  if lsvmin returns zero, then some components of x 
  //             are zero and the rest retain their input values. 
  //  y (out) if lsvmin returns a positive value, then y = (l**-1)*x is an 
  //             unnormalized approximate right singular vector correspond- 
  //             ing to the smallest singular value.  this approximation 
  //             may be crude.  if lsvmin returns zero, then y retains its 
  //             input value.  the caller may pass the same vector for x 
  //             and y (nonstandard fortran usage), in which case y over- 

  //  ***  application and usage restrictions  *** 

  //     there are no usage restrictions. 

  //  ***  algorithm notes  *** 

  //     the algorithm is based on (1), with the additional provision that 
  //     lsvmin = 0 is returned if the smallest diagonal element of l 
  //     (in magnitude) is not more than the unit roundoff times the 
  //     largest.  the algorithm uses a random number generator proposed 
  //     in (4), which passes the spectral test with flying colors -- see 
  //     (2) and (3). 

  //  ***  subroutines and functions called  *** 

  //        v2norm - function, returns the 2-norm of a vector. 

  //  ***  references  *** 

  //     (1) cline, a., moler, c., stewart, g., and wilkinson, j.h.(1977), 
  //         an estimate for the condition number of a matrix, report 
  //         tm-310, applied math. div., argonne national laboratory. 

  //     (2) hoaglin, d.c. (1976), theoretical properties of congruential 
  //         random-number generators --  an empirical view, 
  //         memorandum ns-340, dept. of statistics, harvard univ. 

  //     (3) knuth, d.e. (1969), the art of computer programming, vol. 2 
  //         (seminumerical algorithms), addison-wesley, reading, mass. 

  //     (4) smith, c.s. (1971), multiplicative pseudo-random number 
  //         generators with prime modulus, j. assoc. comput. mach. 18, 
  //         pp. 586-593. 

  //  ***  history  *** 

  //     designed and coded by david m. gay (winter 1977/summer 1978). 

  //  ***  general  *** 

  //     this subroutine was written in connection with research 
  //     supported by the national science foundation under grants 
  //     mcs-7600324, dcr75-10143, 76-1431ss, and mcs76-11989. 

  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  //  ***  local variables  *** 

  /* Parameter adjustments */
  --y;
  --x;
  --l;

  /* Function Body */

  //  ***  body  *** 

  //  ***  first check whether to return lsvmin = 0 and initialize x  *** 

  ii = 0;
  i_1 = *p;
  for (i_ = 1; i_ <= i_1; ++i_) {
    x[i_] = zero;
    ii += i_;
    if (l[ii] == zero) {
      goto L300;
    }
  }
  if (ix % 9973 == 0) {
    ix = 2;
  }
  pplus1 = *p + 1;

  //  ***  solve (l**t)*x = b, where the components of b have randomly 
  //  ***  chosen magnitudes in (.5,1) with signs chosen to make x large. 

  //     do j = p to 1 by -1... 
  i_1 = *p;
  for (jjj = 1; jjj <= i_1; ++jjj) {
    j = pplus1 - jjj;
    //       ***  determine x(j) in this iteration. note for i = 1,2,...,j
    //       ***  that x(i) holds the current partial sum for row i. 
    ix = ix * 3432 % 9973;
    b = half * (one + static_cast<double>(ix) / r9973);
    xplus = b - x[j];
    xminus = -b - x[j];
    splus = fabs(xplus);
    sminus = fabs(xminus);
    jm1 = j - 1;
    j_0 = j * jm1 / 2;
    jj = j_0 + j;
    xplus /= l[jj];
    xminus /= l[jj];
    if (jm1 == 0) {
      goto L30;
    }
    i_2 = jm1;
    for (i_ = 1; i_ <= i_2; ++i_) {
      ji = j_0 + i_;
      //
      //            splus += (d_1 = x[i_] + l[ji] * xplus, abs(d_1));
      //            sminus += (d_1 = x[i_] + l[ji] * xminus, abs(d_1));
      //
      splus += fabs(x[i_] + l[ji] * xplus);
      sminus += fabs(x[i_] + l[ji] * xminus);
    }
  L30:
    if (sminus > splus) {
      xplus = xminus;
    }
    x[j] = xplus;
    ASSERT(isFinite(x[j]));
    //       ***  update partial sums  *** 
    if (jm1 == 0) {
      goto L100;
    }
    i_2 = jm1;
    for (i_ = 1; i_ <= i_2; ++i_) {
      ji = j_0 + i_;
      x[i_] += l[ji] * xplus;
      ASSERT(isFinite(x[i_]));
    }
  L100:
    ;
  }

  //  ***  normalize x  *** 

  t = one / v2norm(p, &x[1]);
  i_1 = *p;
  for (i_ = 1; i_ <= i_1; ++i_) {
    x[i_] = t * x[i_];
    ASSERT(isFinite(x[i_]));
  }

  //  ***  solve l*y = x and return svmin = 1/twonorm(y)  *** 

  i_1 = *p;
  for (j = 1; j <= i_1; ++j) {
    psj = zero;
    jm1 = j - 1;
    j_0 = j * jm1 / 2;
    if (jm1 == 0) {
      goto L130;
    }
    i_2 = jm1;
    for (i_ = 1; i_ <= i_2; ++i_) {
      ji = j_0 + i_;
      psj += l[ji] * y[i_];
    }
  L130:
    jj = j_0 + j;
    y[j] = (x[j] - psj) / l[jj];
  }

  ret_val = one / v2norm(p, &y[1]);
  goto L999;

L300:
  ret_val = zero;
L999:
  return ret_val;
} /* lsvmin */

int qapply(int *nn, int *n, int *p, double *j, double *r_, int *ierr) {
  /* System generated locals */
  int j_dim1, j_offset, i_1, i_2;

  /* Local variables */
  int i_, k, l;
  double t;
  int nl1;

  //     *****parameters. 

  //     .................................................................. 
  //     .................................................................. 

  //     *****purpose. 
  //     this subroutine applies to r the orthogonal transformations 
  //     stored in j by qrfact 

  //     *****parameter description. 
  //     on input:

  //        nn is the row dimension of the matrix j as declared in 
  //             the calling program dimension statement 

  //        n is the number of rows of j and the size of the vector r 

  //        p is the number of columns of j and the size of sigma 

  //        j contains on and below its diagonal the column vectors 
  //             u which determine the householder transformations 
  //             ident - u*u.transpose 

  //        r is the right hand side vector to which the orthogonal 
  //             transformations will be applied 

  //        ierr if non-zero indicates that not all the transformations 
  //             were successfully determined and only the first 
  //             abs(ierr) - 1 transformations will be used 

  //     on output:

  //        r has been overwritten by its transformed image 

  //     *****application and usage restrictions. 
  //     none 

  //     *****algorithm notes. 
  //     the vectors u which determine the householder transformations 
  //     are normalized so that their 2-norm squared is 2.  the use of 
  //     these transformations here is in the spirit of (1). 

  //     *****subroutines and functions called. 

  //     dotprd - function, returns the inner product of vectors 

  //     *****references. 
  //     (1) businger, p. a., and golub, g. h. (1965), linear least squares 
  //        solutions by householder transformations, numer. math. 7, 
  //        pp. 269-276. 

  //     *****history. 
  //     designed by david m. gay, coded by stephen c. peters (winter 1977) 

  //     *****general. 

  //     this subroutine was written in connection with research 
  //     supported by the national science foundation under grants 
  //     mcs-7600324, dcr75-10143, 76-1431ss, and mcs76-11989. 

  //     .................................................................. 
  //     .................................................................. 

  /* Parameter adjustments */
  --r_;
  j_dim1 = *nn;
  j_offset = j_dim1 + 1;
  j -= j_offset;

  /* Function Body */
  k = *p;
  if (*ierr != 0) {
    k = abs(*ierr) - 1;
  }
  if (k == 0) {
    goto L999;
  }

  i_1 = k;
  for (l = 1; l <= i_1; ++l) {
    nl1 = *n - l + 1;
    t = -dotprd(&nl1, &j[l + l * j_dim1], &r_[l]);
    ASSERT(isFinite(t));
    i_2 = *n;
    for (i_ = l; i_ <= i_2; ++i_) {
      r_[i_] += t * j[i_ + l * j_dim1];
      ASSERT(isFinite(r_[i_]));
    }
  }
L999:
  return 0;
} /* qapply */

int qrfact(int* nm, int* m, int* n, double* qr, double* alpha, int* ipivot, int* ierr, int nopivk, double* sum) {
  /* Initialized data */

  const double one = 1.;
  const double p01 = .01;
  const double p99 = .99;
  const double zero = 0.;
  double rktol = 0.;
  double ufeta = 0.;

  /* System generated locals */
  int qr_dim1, qr_offset, i_1, i_2, i_3;
  double d_1;

  /* Local variables */
  double beta;
  int jbar;
  double temp, qrkk, sumj;
  int i_, j, k;
  double sigma;
  int minum, k1;
  double rktol1;
  double alphak;
  double qrkmax;
  int mk1;

  //  ***  compute the qr decomposition of the matrix stored in qr  *** 

  //     *****parameters. 
  //     *****local variables. 
  //     *****functions. 

  // dotprd... returns inner product of two vectors. 
  // rmdcon... returns machine-dependent constants. 
  // vaxpy... computes scalar times one vector plus another. 
  // vscopy... sets all elements of a vector to a scalar. 
  // v2norm... returns the 2-norm of a vector. 

  //     *****constants. 

  /* Parameter adjustments */
  --sum;
  --ipivot;
  --alpha;
  qr_dim1 = *nm;
  qr_offset = qr_dim1 + 1;
  qr -= qr_offset;

  /* Function Body */

  //     .................................................................. 
  //     .................................................................. 


  //     *****purpose. 

  //     this subroutine does a qr-decomposition on the m x n matrix qr, 
  //        with an optionally modified column pivoting, and returns the 
  //        upper triangular r-matrix, as well as the orthogonal vectors 
  //        used in the transformations. 

  //     *****parameter description. 
  //     on input. 

  //        nm must be set to the row dimension of the two dimensional 
  //             array parameters as declared in the calling program 
  //             dimension statement. 

  //        m must be set to the number of rows in the matrix. 

  //        n must be set to the number of columns in the matrix. 

  //        qr contains the real rectangular matrix to be decomposed. 

  //     nopivk is used to control pivotting.  columns 1 through 
  //        nopivk will remain fixed in position. 

  //        sum is used for temporary storage for the subroutine. 

  //     on output. 

  //        qr contains the non-diagonal elements of the r-matrix 
  //             in the strict upper triangle. the vectors u, which 
  //             define the householder transformations   i - u*u-transp, 
  //             are in the columns of the lower triangle. these vectors u 
  //             are scaled so that the square of their 2-norm is 2.0. 

  //        alpha contains the diagonal elements of the r-matrix. 

  //        ipivot reflects the column pivoting performed on the input 
  //             matrix to accomplish the decomposition. the j-th 
  //             element of ipivot gives the column of the original 
  //             matrix which was pivoted into column j during the 
  //             decomposition. 

  //        ierr is set to. 
  //             0 for normal return, 
  //             k if no non-zero pivot could be found for the k-th 
  //                  transformation, or 
  //             -k for an error exit on the k-th thansformation. 
  //             transformations are correct. 


  //     *****applications and usage restrictions. 
  //     this may be used when solving linear least-squares problems -- 
  //     see subroutine qr1 of rosepack.  it is called for this purpose 
  //     by llsqst in the nl2sol (nonlinear least-squares) package. 

  //     *****algorithm notes. 
  //     this version of qrfact tries to eliminate the occurrence of 
  //     underflows during the accumulation of inner products.  rktol1 
  //     is chosen below so as to insure that discarded terms have no 
  //     effect on the computed two-norms. 

  //     adapted from the algol routine solve (1). 

  //     *****references. 
  //     (1)     businger,p. and golub,g.h., linear least squares 
  //     solutions by housholder transformations, in wilkinson,j.h. 
  //     and reinsch,c.(eds.), handbook for automatic computation, 
  //     volume ii. linear algebra, springer-verlag, 111-118 (1971). 
  //     prepublished in numer.math. 7, 269-276 (1965). 

  //     *****history. 
  //     this amounts to the subroutine qr1 of rosepack with rktol1 used 
  //     in place of rktol below, with v2norm used to initialize (and 
  //     sometimes update) the sum array, and with calls on dotprd and 
  //     vaxpy in place of some loops. 

  //     *****general. 

  //     development of this program supported in part by 
  //     national science foundation grant gj-1154x3 and 
  //     national science foundation grant dcr75-08802 
  //     to national bureau of economic research, inc. 



  //     .................................................................. 
  //     .................................................................. 


  //     ..........  ufeta is the smallest positive floating point number 
  //        s.t. ufeta and -ufeta can both be represented. 

  //     ..........  rktol is the square root of the relative precision 
  //        of floating point arithmetic (machep). 
  //     *****body of program. 
  if (ufeta > zero) {
    goto L10;
  }
  ufeta = rmdcon(1);
  rktol = rmdcon(4);
L10:
  *ierr = 0;
  rktol1 = p01 * rktol;

  i_1 = *n;
  for (j = 1; j <= i_1; ++j) {
    sum[j] = v2norm(m, &qr[j * qr_dim1 + 1]);
    ASSERT(isFinite(sum[j]));
    ipivot[j] = j;
  }

  minum = std::min(*m, *n);

  i_1 = minum;
  for (k = 1; k <= i_1; ++k) {
    mk1 = *m - k + 1;
    //        ..........k-th householder transformation.......... 
    sigma = zero;
    jbar = 0;
    //        ..........find largest column sum.......... 
    if (k <= nopivk) {
      goto L50;
    }
    i_2 = *n;
    for (j = k; j <= i_2; ++j) {
      if (sigma >= sum[j]) {
        goto L30;
      }
      sigma = sum[j];
      jbar = j;
    L30:
      ;
    }

    if (jbar == 0) {
      goto L220;
    }
    if (jbar == k) {
      goto L50;
    }
    //        ..........column interchange.......... 
    i_ = ipivot[k];
    ipivot[k] = ipivot[jbar];
    ipivot[jbar] = i_;
    sum[jbar] = sum[k];
    sum[k] = sigma;

    i_2 = *m;
    for (i_ = 1; i_ <= i_2; ++i_) {
      sigma = qr[i_ + k * qr_dim1];
      qr[i_ + k * qr_dim1] = qr[i_ + jbar * qr_dim1];
      qr[i_ + jbar * qr_dim1] = sigma;
    }
    //        ..........end of column interchange.......... 
  L50:
    //        ..........  second inner product  .......... 
    qrkmax = zero;

    i_2 = *m;
    for (i_ = k; i_ <= i_2; ++i_) {
      //
      //            if ((d_1 = qr[i_ + k * qr_dim1], abs(d_1)) > qrkmax) {
      //                qrkmax = (d_2 = qr[i_ + k * qr_dim1], abs(d_2));
      if (fabs(qr[i_ + k * qr_dim1]) > qrkmax) {
        qrkmax = fabs(qr[i_ + k * qr_dim1]);
      }
    }

    if (qrkmax < ufeta) {
      goto L210;
    }
    alphak = v2norm(&mk1, &qr[k + k * qr_dim1]) / qrkmax;
    ASSERT(isFinite(alphak));
    // Computing 2nd power 
    d_1 = alphak;
    sigma = d_1 * d_1;

    //        ..........  end second inner product  .......... 
    qrkk = qr[k + k * qr_dim1];
    if (qrkk >= zero) {
      alphak = -alphak;
      ASSERT(isFinite(alphak));
    }
    alpha[k] = alphak * qrkmax;
    beta = qrkmax * sqrt(sigma - qrkk * alphak / qrkmax);
    qr[k + k * qr_dim1] = qrkk - alpha[k];
    i_2 = *m;
    for (i_ = k; i_ <= i_2; ++i_) {
      qr[i_ + k * qr_dim1] /= beta;
    }
    k1 = k + 1;
    if (k1 > *n) {
      goto L120;
    }

    i_2 = *n;
    for (j = k1; j <= i_2; ++j) {
      temp = -dotprd(&mk1, &qr[k + k * qr_dim1], &qr[k + j * qr_dim1]);

      //             ***  set qr(i,j) = qr(i,j) + temp*qr(i,k), i = k,..

      vaxpy(&mk1, &qr[k + j * qr_dim1], temp, &qr[k + k * qr_dim1], &qr[k + j * qr_dim1]);

      if (k1 > *m) {
        goto L110;
      }
      sumj = sum[j];
      if (sumj < ufeta) {
        goto L110;
      }
      //
      //            temp = (d_1 = qr[k + j * qr_dim1] / sumj, abs(d_1));
      //
      temp = fabs(qr[k + j * qr_dim1] / sumj);
      if (temp < rktol1) {
        goto L110;
      }
      if (temp >= p99) {
        goto L90;
      }
      // Computing 2nd power 
      d_1 = temp;
      sum[j] = sumj * sqrt(one - d_1 * d_1);
      goto L110;
    L90:
      i_3 = *m - k;
      sum[j] = v2norm(&i_3, &qr[k1 + j * qr_dim1]);
    L110:
      ;
    }
    //        ..........end of k-th householder transformation.......... 
  L120:
    ;
  }

  goto L999;
  //     ..........error exit on k-th transformation.......... 
L210:
  *ierr = -k;
  goto L230;
  //     ..........no non-zero acceptable pivot found.......... 
L220:
  *ierr = k;
L230:
  i_1 = *n;
  for (i_ = k; i_ <= i_1; ++i_) {
    alpha[i_] = zero;
    if (i_ > k) {
      i_2 = i_ - k;
      vscopy(&i_2, &qr[k + i_ * qr_dim1], zero);
    }
  }
  //     ..........return to caller.......... 
L999:
  return 0;
} /* qrfact */


}


const double NL2SOL::Defaults::RelativeTolerance = sqrt(DBL_EPSILON);
const double NL2SOL::Defaults::AbsoluteTolerance = 0.5*DBL_EPSILON;
const size_t NL2SOL::Defaults::MaximumNumberOfIterations = 10 * DBL_DIG;

bool NL2SOL::is_number(double x) { return !isNaN(x) && !isInf(x); }

bool NL2SOL::any_not_is_number(const std::vector<double>&x){ return any_not_is_number(&x[0],x.size()); }

bool NL2SOL::any_not_is_number(const double *x, size_t n_x) {
  for (size_t i = 0;i < n_x;++i)
    if (!is_number(x[i]))
      return true;
  return false;
}


NL2SOL::NL2SOL(double argumentTolerance, // This is a *relative* convergence tolerance.
               double functionTolerance, // This is a *relative* convergence tolerance.
               double absoluteFunctionTolerance, // 0.5*functionTolerance*functionTolerance
               double finiteDifferencingWidth,   // for Jacobian computation, ignored when jacobian function is given
               size_t maximumNumberOfIterations) : 
ArgumentTolerance(argumentTolerance), FunctionTolerance(functionTolerance), AbsoluteFunctionTolerance(absoluteFunctionTolerance),
FiniteDifferencingWidth(finiteDifferencingWidth) /* for Jacobian computation */, MaximumNumberOfIterations(maximumNumberOfIterations),
AchievedAccuracy(std::numeric_limits<double>::quiet_NaN()), RMS(std::numeric_limits<double>::quiet_NaN()),
ReturnCode(0), Iterations(0), FunctionEvaluations(0), JacobianEvaluations(0)
{
}

double NL2SOL::tryToSolve(Functor& functor, std::vector<double> &x, size_t n_f){
   return tryToSolve(functor, &x[0], x.size(), n_f);
}

double NL2SOL::tryToSolve(Functor& functor, double *x, size_t n_x, size_t n_f) {

  if (0 == n_x)
    throw std::runtime_error("cannot solve minimisation problem when the given number of parameters to solve for is zero. This is an ill-specified problem.");
  if (0 == n_f)
    throw std::runtime_error("cannot solve minimisation problem when the given number of objective function dimensions is zero. This is an ill-specified problem.");
  size_t lv = (105 + n_x*(n_f + 2 * n_x + 17) + 4 * n_f);
  size_t liv = (80 + n_x);
  m_workspace.resize(lv);
  m_iv.resize(liv);

  const double tolx = ArgumentTolerance, tolf = FunctionTolerance, atolf = AbsoluteFunctionTolerance;
  const int ntrial = (int)MaximumNumberOfIterations;

  const bool have_jacobian = functor.have_jacobian();

  ReturnCode = 0;
  Iterations = 0;
  FunctionEvaluations = 0;
  JacobianEvaluations = 0;
  AchievedAccuracy = RMS = std::numeric_limits<double>::quiet_NaN();

  fill(m_workspace, 0.0);
  fill(m_iv, 0);
  nl2sol_defaults(&m_iv[0], &m_workspace[0]);      //   Obtain default control values.

  int p = static_cast<int>(n_x);
  int n = static_cast<int>(n_f);

  if (!have_jacobian)
    m_workspace[35] = FiniteDifferencingWidth;               //   Finite differencing distance for the Jacobian. We leave it at the default value if we have an analytical Jacobian.
  m_workspace[36] = 1.;                                      //   Scaling of the errors.
  m_iv[13] = 0;                                              //   Don't calculate covariance matrix at the solution.
  m_iv[14] = 0;                                              //   Don't calculate covariance matrix at the solution.
  m_iv[15] = 0;
  m_workspace[37] = 1.;                                      //   Even more scaling. See nl2sol itself.
  m_iv[16] = ntrial*(p*p + 2);                               //   Maximum number of function calls.
  m_iv[17] = ntrial;                                         //   Maximum number of iterations.
  REQUIRE(atolf >= rmdcon(1), "NL2SOL::tryToSolve : the function absolute tolerance should be above " + toString(rmdcon(1)) + ", which is the smallest positive number that can be represented on this hardware.");
  m_workspace[30] = atolf;                                   //   This controls absolute tolerance.
  REQUIRE(tolf >= rmdcon(3), "NL2SOL::tryToSolve : the function relative tolerance should be above " + toString(rmdcon(3)) + ", which is the smallest number e such as 1+e!=e.");
  m_workspace[31] = tolf;                                    //   This is the actual relative tolerance.
  m_workspace[32] = tolx;                                    //   This is the relative convergence value (x-tolerance).
  m_workspace[33] = 0.;                                      //   This is the false convergence threshold. We need to avoid spurious false convergence reports.
  m_workspace[9] = std::numeric_limits<double>::quiet_NaN(); //   Preset the function value to something indicating that no evaluation has occurred yet.
                                                             //   Now call the ACM-TOMS-573 algorithm.
  if (have_jacobian)
    nl2sol(functor, &FunctionEvaluations, &JacobianEvaluations, &n, &p, x, &m_iv[0], &m_workspace[0]);
  else
    nl2sno(functor, &FunctionEvaluations, &JacobianEvaluations, &n, &p, x, &m_iv[0], &m_workspace[0]);
  // The "function value" is documented as "half the sum of squares".
  // We prefer the L₂-Norm.
  AchievedAccuracy = sqrt(2 * m_workspace[9]);
  RMS = AchievedAccuracy / sqrt((double)AchievedAccuracy);
  Iterations = m_iv[30];
  //
  // You'll find more detailed explanations of the error codes further down in NL2SOL::nl2sol().
  //
  ReturnCode = m_iv[0];
  switch (ReturnCode) {
  case 8:
    //             8 = false convergence.  the iterates appear to be converg- 
    //                  ing to a noncritical point.  this may mean that the 
    //                  convergence tolerances (v(afctol), v(rfctol), 
    //                  v(xctol)) are too small for the accuracy to which 
    //                  the function and gradient are being computed, that 
    //                  there is an error in computing the gradient, or that 
    //                  the function or gradient is discontinuous near x.
#ifdef SHOW_ACM_TOMS_573_WARNINGS
    fprintf(stderr, "Warning: The ACM-TOMS-573 algorithm (nl2sno) came to a false convergence halt.\n");
#endif
    break;
  case 9:
#ifdef SHOW_ACM_TOMS_573_WARNINGS
    fprintf(stderr, "Warning: The ACM-TOMS-573 algorithm (nl2sno) reached the function evaluation limit without other convergence.\n");
#endif
    break;
  case 10:
#ifdef SHOW_ACM_TOMS_573_WARNINGS
    fprintf(stderr, "Warning: The ACM-TOMS-573 algorithm (nl2sno) reached the iteration limit without other convergence.\n");
#endif
    ++Iterations;
    break;
  case 11:   // Interrupted by user or stopx() was flagged by the functor.
    Iterations = -Iterations;
    break;
  case 13:
    throw std::runtime_error("Invalid starting coordinates passed on to NL2SOL::tryToSolve()");
    break;
  case 15:
    throw std::runtime_error("Unable to calculate Jacobian in NL2SOL::tryToSolve()");
    break;
  case 16:
    if (p > n)
      throw std::runtime_error("The number of arguments must not be larger than the number of errors to minimise in NL2SOL::tryToSolve()");
    break;
  default:break;
  }
  if (ReturnCode > 11) {   //   A serious error from the NL2SOL algorithm ?
    std::string msg = "The ACM-TOMS-573 algorithm (nl2s";
    msg += (have_jacobian ? "ol" : "no");
    msg += ") returned error value ";
    msg += toString(ReturnCode);
    msg += " in NL2SOL::tryToSolve() which is almost certainly a programming error";
    throw std::runtime_error(msg);
  }
  return AchievedAccuracy;
}

struct FunctorWithFunctionPointers : public NL2SOL::Functor {
   typedef int (*Objective_fp)(const double *x, size_t n_x, double *f, size_t n_f);
   typedef int (*Jacobian_fp)(const double *x, size_t n_x, size_t n_f, double *jac);
   typedef int (*Stop_fp)();
   Objective_fp m_objective;
   Jacobian_fp  m_jacobian;
   Stop_fp m_stop;
   FunctorWithFunctionPointers(Objective_fp objective, Jacobian_fp jacobian, Stop_fp stop) :
      m_objective(objective), m_jacobian(jacobian), m_stop(stop){
   }
   virtual bool operator()(const double *x, size_t n_x, double *f /* f=f(x) */, size_t n_f){ return 0 != m_objective(x,n_x,f,n_f); }
   virtual bool have_jacobian() const { return 0 != m_jacobian; }
   virtual bool jacobian(const double *x, size_t n_x, size_t n_f, double *J) { return 0 != m_jacobian(x,n_x,n_f,J); }
   virtual bool stop(){ return 0 != m_stop && 0 != m_stop(); }
};

char* copy_string(const char* s) {
  char* d = 0;
  if (s) {
    size_t l = strlen(s);
    d = (char*)malloc(l + 1);
    if (0!=d) {
      memcpy(d, s, l);
      d[l] = 0;
    }
  }
  return d;
}

extern "C" DLL_EXPORT void c_free(void *s){
   free(s);
}

struct Data {
  // INPUT
  double ArgumentTolerance, FunctionTolerance, AbsoluteFunctionTolerance, FiniteDifferencingWidth; // for Jacobian computation
  size_t  MaximumNumberOfIterations;
  // OUTPUT
  double AchievedAccuracy; // L₂-Norm on exit.
  double RMS; // L₂-Norm/√n_f on exit, where n_f the number of function values in the last invocation to tryToSolve().
  int ReturnCode, Iterations, FunctionEvaluations, JacobianEvaluations;
  const char* ErrorMessage; // Requires manual free()
  //
  Data(double argumentTolerance = NL2SOL::Defaults::RelativeTolerance, // This is a *relative* convergence tolerance.
       double functionTolerance = NL2SOL::Defaults::RelativeTolerance, // This is a *relative* convergence tolerance.
       double absoluteFunctionTolerance = NL2SOL::Defaults::AbsoluteTolerance, // 0.5*functionTolerance*functionTolerance
       double finiteDifferencingWidth = NL2SOL::Defaults::RelativeTolerance, // for Jacobian computation, ignored when jacobian function is given
       size_t maximumNumberOfIterations = NL2SOL::Defaults::MaximumNumberOfIterations) :
    ArgumentTolerance(argumentTolerance), FunctionTolerance(functionTolerance), AbsoluteFunctionTolerance(absoluteFunctionTolerance),
    FiniteDifferencingWidth(finiteDifferencingWidth) /* for Jacobian computation */, MaximumNumberOfIterations(maximumNumberOfIterations),
    AchievedAccuracy(std::numeric_limits<double>::quiet_NaN()), RMS(std::numeric_limits<double>::quiet_NaN()), ReturnCode(0), Iterations(0),
    FunctionEvaluations(0), JacobianEvaluations(0), ErrorMessage(0) {}
};

void set_error_message(Data* data, const char *msg){
   if (data) data->ErrorMessage = copy_string(msg);
}

extern "C" DLL_EXPORT void nl2sol_default_data(Data *data){
   if (data) new (data) Data();
}

extern "C" DLL_EXPORT double nl2sol(int (*objective)(const double *x, size_t n_x, double *f, size_t n_f),
                              int (*jacobian)(const double *x, size_t n_x, size_t n_f, double *jac),
                              int (*stop)(), double *x, size_t n_x, size_t n_f, Data *data)
{
   NL2SOL solver;
   if (data) {
      data->ErrorMessage = 0;
      solver.ArgumentTolerance = data->ArgumentTolerance;
      solver.FunctionTolerance = data->FunctionTolerance;
      solver.AbsoluteFunctionTolerance = data->AbsoluteFunctionTolerance;
      solver.FiniteDifferencingWidth = data->FiniteDifferencingWidth;
      solver.MaximumNumberOfIterations = data->MaximumNumberOfIterations;
   }
   FunctorWithFunctionPointers functor(objective,jacobian,stop);
   try {
      solver.tryToSolve(functor,x,n_x,n_f);
   }
   catch (const std::exception& e) { set_error_message(data,e.what());        }
   catch (const std::string& e) {    set_error_message(data,e.c_str());       }
   catch (const char* e) {           set_error_message(data,e);               }
   catch (...) {                     set_error_message(data,"Unknown error"); }
   if (data) {
      data->AchievedAccuracy = solver.AchievedAccuracy;
      data->RMS = solver.RMS;
      data->ReturnCode = solver.ReturnCode;
      data->Iterations = solver.Iterations;
      data->FunctionEvaluations = solver.FunctionEvaluations;
      data->JacobianEvaluations = solver.JacobianEvaluations;
   }
   return solver.AchievedAccuracy;
}
