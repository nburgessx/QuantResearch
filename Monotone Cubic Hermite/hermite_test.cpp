# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

# include "hermite.hpp"

using namespace std;

int main ( );
void test01 ( );
void test02 ( );
void test03 ( );
void test04 ( );
void test05 ( );
void test06 ( );
void test07 ( );
void test08 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for HERMITE_TEST.
//
//  Discussion:
//
//    HERMITE_TEST tests the HERMITE library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 November 2011
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "HERMITE_TEST\n";
  cout << "  C++ version\n";
  cout << "  Test the HERMITE library.\n";

  test01 ( );
  test02 ( );
  test03 ( );
  test04 ( );
  test05 ( );
  test06 ( );
  test07 ( );
  test08 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "HERMITE_TEST\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void test01 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST01 uses f(x) = 1 + 2x + 3x^2 at x = 0, 1, 2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2011
//
//  Author:
//
//    John Burkardt
//
{
# define N 3

  int n = N;
  double x[N] = { 0.0, 1.0,  2.0 };
  double y[N] = { 1.0, 6.0, 17.0 };
  double yp[N] = { 2.0, 8.0, 14.0 };

  cout << "\n";
  cout << "TEST01\n";
  cout << "  HERMITE computes the Hermite interpolant to data.\n";
  cout << "  Here, f(x) = 1 + 2x + 3x^2.\n";

  hermite_demo ( n, x, y, yp );

  return;
# undef N
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 uses f(x) = 6 + 5x + 4x^2 + 3x^3 + 2x^4 + x^5 at x = 0, 1, 2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2011
//
//  Author:
//
//    John Burkardt
//
{
# define N 3

  int i;
  int n = N;
  double *x;
  double *y;
  double *yp;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  HERMITE computes the Hermite interpolant to data.\n";
  cout << "  Here, f(x) = 6 + 5x + 4x^2 + 3x^3 + 2x^4 + x^5.\n";

  x = new double[n];
  y = new double[n];
  yp = new double[n];

  for ( i = 0; i < n; i++ )
  {
    x[i] = ( double ) ( i );
    
    y[i] = 6.0 + x[i] * ( 
           5.0 + x[i] * ( 
           4.0 + x[i] * ( 
           3.0 + x[i] * ( 
           2.0 + x[i] ) ) ) );

    yp[i] = 5.0 + x[i] * ( 
            8.0 + x[i] * ( 
            9.0 + x[i] * ( 
            8.0 + x[i] *   
            5.0 ) ) );
  }
  hermite_demo ( n, x, y, yp );

  delete [] x;
  delete [] y;
  delete [] yp;

  return;
# undef N
}
//****************************************************************************80

void test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 uses f(x) = r1 + r2x + r3x^2 + r4x^3 + r5x^4 + r6x^5 at x = r7 r8 r9
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2011
//
//  Author:
//
//    John Burkardt
//
{
# define N 3

  double *c;
  int i;
  int n = N;
  int seed;
  double *x;
  double *y;
  double *yp;

  cout << "\n";
  cout << "TEST03\n";
  cout << "  HERMITE computes the Hermite interpolant to data.\n";
  cout << "  Here, f(x) is a fifth order polynomial with random\n";
  cout << "  coefficients, and the abscissas are random.\n";

  c = new double[2*n];
  x = new double[n];
  y = new double[n];
  yp = new double[n];

  seed = 123456789;

  r8vec_uniform_01 ( n, &seed, x );
  r8vec_print ( n, x, "  Random abscissas" );

  r8vec_uniform_01 ( 2 * n, &seed, c );
  r8vec_print ( 2 * n, c, "  Random polynomial coefficients." );

  for ( i = 0; i < n; i++ )
  {
    y[i] = c[0] + x[i] * ( 
           c[1] + x[i] * ( 
           c[2] + x[i] * ( 
           c[3] + x[i] * ( 
           c[4] + x[i] * ( 
           c[5] ) ) ) ) );

    yp[i] = c[1]       + x[i] * ( 
            c[2] * 2.0 + x[i] * ( 
            c[3] * 3.0 + x[i] * ( 
            c[4] * 4.0 + x[i] *   
            c[5] * 5.0 ) ) );
  }

  hermite_demo ( n, x, y, yp );

  delete [] c;
  delete [] x;
  delete [] y;
  delete [] yp;

  return;
# undef N
}
//****************************************************************************80

void test04 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST04 interpolates the Runge function using equally spaced data.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    31 October 2011
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  double max_dif;
  int n;
  int nd;
  int ndp;
  int ns;
  double *x;
  double *xd;
  double *xdp;
  double *xs;
  double xhi;
  double xlo;
  double xt;
  double *y;
  double *yd;
  double *ydp;
  double *yp;
  double *ys;
  double yt;

  cout << "\n";
  cout << "TEST04\n";
  cout << "  HERMITE computes the Hermite interpolant to data.\n";
  cout << "  Here, f(x) is the Runge function\n";
  cout << "  and the data is evaluated at equally spaced points.\n";
  cout << "  As N increases, the maximum error grows.\n";
  cout << "\n";
  cout << "     N     Max | F(X) - H(F(X)) |\n";
  cout << "\n";

  for ( n = 3; n <= 15; n = n + 2 )
  {
    y = new double[n];
    yp = new double[n];

    nd = 2 * n;
    xd = new double[nd];
    yd = new double[nd];

    ndp = 2 * n - 1;
    xdp = new double[ndp];
    ydp = new double[ndp];

    ns = 10 * ( n - 1 ) + 1;

    xlo = -5.0;
    xhi = +5.0;
    x = r8vec_linspace_new ( n, xlo, xhi );

    for ( i = 0; i < n; i++ )
    {
      y[i] = 1.0 / ( 1.0 + x[i] * x[i] );
      yp[i] = - 2.0 * x[i] / ( 1.0 + x[i] * x[i] ) / ( 1.0 + x[i] * x[i] );
    }

    hermite_interpolant ( n, x, y, yp, xd, yd, xdp, ydp );
//
//  Compare exact and interpolant at sample points.
//
    xs = r8vec_linspace_new ( ns, xlo, xhi );

    ys = dif_vals ( nd, xd, yd, ns, xs );

    max_dif = 0.0;
    for ( i = 0; i < ns; i++ )
    {
      xt = xs[i];
      yt = 1.0 / ( 1.0 + xt * xt );
      max_dif = r8_max ( max_dif, r8_abs ( ys[i] - yt ) );
    }

    cout << "  " << setw(4) << n
         << "  " << setw(14) << max_dif << "\n";

    delete [] x;
    delete [] xd;
    delete [] xdp;
    delete [] xs;
    delete [] y;
    delete [] yd;
    delete [] ydp;
    delete [] yp;
    delete [] ys;;
  }

  return;
}
//****************************************************************************80

void test05 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST05 interpolates the Runge function using Chebyshev spaced data.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    31 October 2011
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  double max_dif;
  int n;
  int nd;
  int ndp;
  int ns;
  double *x;
  double *xd;
  double *xdp;
  double *xs;
  double xhi;
  double xlo;
  double xt;
  double *y;
  double *yd;
  double *ydp;
  double *yp;
  double *ys;
  double yt;

  cout << "\n";
  cout << "TEST05\n";
  cout << "  HERMITE computes the Hermite interpolant to data.\n";
  cout << "  Here, f(x) is the Runge function\n";
  cout << "  and the data is evaluated at Chebyshev spaced points.\n";
  cout << "  As N increases, the maximum error decreases.\n";
  cout << "\n";
  cout << "     N     Max | F(X) - H(F(X)) |\n";
  cout << "\n";

  for ( n = 3; n <= 15; n = n + 2 )
  {
    y = new double[n];
    yp = new double[n];

    nd = 2 * n;
    xd = new double[nd];
    yd = new double[nd];

    ndp = 2 * n - 1;
    xdp = new double[ndp];
    ydp = new double[ndp];

    ns = 10 * ( n - 1 ) + 1;

    xlo = -5.0;
    xhi = +5.0;
    x = r8vec_chebyshev_new ( n, xlo, xhi );

    for ( i = 0; i < n; i++ )
    {
      y[i] = 1.0 / ( 1.0 + x[i] * x[i] );
      yp[i] = - 2.0 * x[i] / ( 1.0 + x[i] * x[i] ) / ( 1.0 + x[i] * x[i] );
    }

    hermite_interpolant ( n, x, y, yp, xd, yd, xdp, ydp );
//
//  Compare exact and interpolant at sample points.
//
    xs = r8vec_linspace_new ( ns, xlo, xhi );

    ys = dif_vals ( nd, xd, yd, ns, xs );

    max_dif = 0.0;
    for ( i = 0; i < ns; i++ )
    {
      xt = xs[i];
      yt = 1.0 / ( 1.0 + xt * xt );
      max_dif = r8_max ( max_dif, r8_abs ( ys[i] - yt ) );
    }

    cout << "  " << setw(4) << n
         << "  " << setw(14) << max_dif << "\n";

    delete [] x;
    delete [] xd;
    delete [] xdp;
    delete [] xs;
    delete [] y;
    delete [] yd;
    delete [] ydp;
    delete [] yp;
    delete [] ys;;
  }

  return;
}
//***************************************************************************80

void test06 ( )

//***************************************************************************80
//
//  Purpose:
//
//    TEST06 tests HERMITE_BASIS_0 and HERMITE_BASIS_1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 May 2011
//
//  Author:
//
//    John Burkardt
//
{
# define ND 2

  double f01;
  double f02;
  double f11;
  double f12;
  int i;
  int j;
  int nd = ND;
  double xd[ND] = { 0.0, 10.0 };
  double xv;
  double yd[ND];
  double yh;
  double ypd[ND];
  double yv;

  cout << "\n";
  cout << "TEST06:\n";
  cout << "  HERMITE_BASIS_0 and HERMITE_BASIS_1 evaluate the\n";
  cout << "  Hermite global polynomial basis functions\n";
  cout << "  of type 0: associated with function values, and\n";
  cout << "  of type 1: associated with derivative values.\n";
//
//  Let y = x^3 + x^2 + x + 1,
//  and compute the Hermite global polynomial interpolant based on two 
//  abscissas:
//
  for ( j = 0; j < nd; j++ )
  {
    yd[j] = pow ( xd[j], 3 ) + pow ( xd[j], 2 ) + xd[j] + 1.0;
    ypd[j] = 3.0 * pow ( xd[j], 2 ) + 2.0 * xd[j] + 1.0;
  }

  cout << "\n";
  cout << "  Interpolate y = x^3 + x^2 + x + 1.\n";
  cout << "\n";
  cout << "     XD         Y(XD)      Y'(XD)\n";
  cout << "\n";
  for ( j = 0; j < nd; j++ )
  {
    cout << "  " << setw(10) << xd[j]
         << "  " << setw(10) << yd[j]
         << "  " << setw(10) << ypd[j] << "\n";
  }

  cout << "\n";
  cout << "     XV         Y(XV)      H(XV)\n";
  cout << "\n";

  for ( j = 0; j <= 10; j++ )
  {
    xv = ( double ) ( j );

    yv = pow ( xv, 3 ) + pow ( xv, 2 ) + xv + 1.0;

    f01 = hermite_basis_0 ( 2, xd, 0, xv );
    f11 = hermite_basis_1 ( 2, xd, 0, xv );
    f02 = hermite_basis_0 ( 2, xd, 1, xv );
    f12 = hermite_basis_1 ( 2, xd, 1, xv );

    yh = yd[0] * f01 + ypd[0] * f11 + yd[1] * f02 + ypd[1] * f12;

    cout << "  " << setw(10) << xv
         << "  " << setw(10) << yv
         << "  " << setw(10) << yh << "\n";
  }
  return;
# undef ND
}
//****************************************************************************80

void test07 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST07 tests HERMITE_INTERPOLANT_RULE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 June 2011
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  int i;
  int k;
  int n;
  double q;
  double *w;
  double *x;

  cout << "\n";
  cout << "TEST07:\n";
  cout << "  HERMITE_INTERPOLANT_RULE\n";
  cout << "  is given a set of N abscissas for a Hermite interpolant\n";
  cout << "  and returns N pairs of quadrature weights\n";
  cout << "  for function and derivative values at the abscissas.\n";

  n = 3;
  a = 0.0;
  b = 10.0;
  x = r8vec_linspace_new ( n, a, b );
  w = hermite_interpolant_rule ( n, a, b, x );

  cout << "\n";
  cout << "     I       X               W(F(X))        W(F'(X))\n";
  cout << "\n";
  k = 0;
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(4) << i
         << "  " << setw(14) << x[i]
         << "  " << setw(14) << w[k]
         << "  " << setw(14) << w[k+1] << "\n";
    k = k + 2;
  }

  cout << "\n";
  cout << "  Use the quadrature rule over interval " << a << " to " << b << "\n";
  cout << "\n";

  q = 0.0;
  k = 0;
  for ( i = 0; i < n; i++ )
  {
    q = q + w[k] * 1 + w[k+1] * 0.0;
    k = k + 2;
  }
  cout << "  Estimate integral of 1 = " << q << "\n";

  q = 0.0;
  k = 0;
  for ( i = 0; i < n; i++ )
  {
    q = q + w[k] * x[i] + w[k+1] * 1.0;
    k = k + 2;
  }
  cout << "  Estimate integral of X = " << q << "\n";

  q = 0.0;
  k = 0;
  for ( i = 0; i < n; i++ )
  {
    q = q + w[k] * x[i] * x[i] + w[k+1] * 2.0 * x[i];
    k = k + 2;
  }
  cout << "  Estimate integral of X^2 = " << q << "\n";

  q = 0.0;
  k = 0;
  for ( i = 0; i < n; i++ )
  {
    q = q + w[k] / ( 1.0 + x[i] * x[i] ) 
          - w[k+1] * 2.0 * x[i] / pow ( 1.0 + x[i] * x[i], 2 );
    k = k + 2;
  }
  cout << "  Estimate integral of 1/(1+x^2) = " << q << "\n";

  delete [] w;
  delete [] x;

  n = 3;
  a = 0.0;
  b = 1.0;
  x = r8vec_linspace_new ( n, a, b );
  w = hermite_interpolant_rule ( n, a, b, x );

  cout << "\n";
  cout << "     I       X               W(F(X))        W(F'(X))\n";
  cout << "\n";
  k = 0;
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(4) << i
         << "  " << setw(14) << x[i]
         << "  " << setw(14) << w[k]
         << "  " << setw(14) << w[k+1] << "\n";
    k = k + 2;
  }

  cout << "\n";
  cout << "  Use the quadrature rule over interval " << a << " to " << b << "\n";
  cout << "\n";

  q = 0.0;
  k = 0;
  for ( i = 0; i < n; i++ )
  {
    q = q + w[k] * 1 + w[k+1] * 0.0;
    k = k + 2;
  }
  cout << "  Estimate integral of 1 = " << q << "\n";

  q = 0.0;
  k = 0;
  for ( i = 0; i < n; i++ )
  {
    q = q + w[k] * x[i] + w[k+1] * 1.0;
    k = k + 2;
  }
  cout << "  Estimate integral of X = " << q << "\n";

  q = 0.0;
  k = 0;
  for ( i = 0; i < n; i++ )
  {
    q = q + w[k] * x[i] * x[i] + w[k+1] * 2.0 * x[i];
    k = k + 2;
  }
  cout << "  Estimate integral of X^2 = " << q << "\n";

  q = 0.0;
  k = 0;
  for ( i = 0; i < n; i++ )
  {
    q = q + w[k] / ( 1.0 + x[i] * x[i] ) 
          - w[k+1] * 2.0 * x[i] / pow ( 1.0 + x[i] * x[i], 2 );
    k = k + 2;
  }
  cout << "  Estimate integral of 1/(1+x^2) = " << q << "\n";

  delete [] w;
  delete [] x;

  n = 11;
  a = 0.0;
  b = 10.0;
  x = r8vec_linspace_new ( n, a, b );
  w = hermite_interpolant_rule ( n, a, b, x );

  cout << "\n";
  cout << "     I       X               W(F(X))        W(F'(X))\n";
  cout << "\n";
  k = 0;
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(4) << i
         << "  " << setw(14) << x[i]
         << "  " << setw(14) << w[k]
         << "  " << setw(14) << w[k+1] << "\n";
    k = k + 2;
  }

  cout << "\n";
  cout << "  Use the quadrature rule over interval " << a << " to " << b << "\n";
  cout << "\n";

  q = 0.0;
  k = 0;
  for ( i = 0; i < n; i++ )
  {
    q = q + w[k] * 1 + w[k+1] * 0.0;
    k = k + 2;
  }
  cout << "  Estimate integral of 1 = " << q << "\n";

  q = 0.0;
  k = 0;
  for ( i = 0; i < n; i++ )
  {
    q = q + w[k] * x[i] + w[k+1] * 1.0;
    k = k + 2;
  }
  cout << "  Estimate integral of X = " << q << "\n";

  q = 0.0;
  k = 0;
  for ( i = 0; i < n; i++ )
  {
    q = q + w[k] * x[i] * x[i] + w[k+1] * 2.0 * x[i];
    k = k + 2;
  }
  cout << "  Estimate integral of X^2 = " << q << "\n";

  q = 0.0;
  k = 0;
  for ( i = 0; i < n; i++ )
  {
    q = q + w[k] / ( 1.0 + x[i] * x[i] ) 
          - w[k+1] * 2.0 * x[i] / pow ( 1.0 + x[i] * x[i], 2 );
    k = k + 2;
  }
  cout << "  Estimate integral of 1/(1+x^2) = " << q << "\n";

  delete [] w;
  delete [] x;

  n = 11;
  a = 0.0;
  b = 1.0;
  x = r8vec_linspace_new ( n, a, b );
  w = hermite_interpolant_rule ( n, a, b, x );

  cout << "\n";
  cout << "     I       X               W(F(X))        W(F'(X))\n";
  cout << "\n";
  k = 0;
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(4) << i
         << "  " << setw(14) << x[i]
         << "  " << setw(14) << w[k]
         << "  " << setw(14) << w[k+1] << "\n";
    k = k + 2;
  }

  cout << "\n";
  cout << "  Use the quadrature rule over interval " << a << " to " << b << "\n";
  cout << "\n";

  q = 0.0;
  k = 0;
  for ( i = 0; i < n; i++ )
  {
    q = q + w[k] * 1 + w[k+1] * 0.0;
    k = k + 2;
  }
  cout << "  Estimate integral of 1 = " << q << "\n";

  q = 0.0;
  k = 0;
  for ( i = 0; i < n; i++ )
  {
    q = q + w[k] * x[i] + w[k+1] * 1.0;
    k = k + 2;
  }
  cout << "  Estimate integral of X = " << q << "\n";

  q = 0.0;
  k = 0;
  for ( i = 0; i < n; i++ )
  {
    q = q + w[k] * x[i] * x[i] + w[k+1] * 2.0 * x[i];
    k = k + 2;
  }
  cout << "  Estimate integral of X^2 = " << q << "\n";

  q = 0.0;
  k = 0;
  for ( i = 0; i < n; i++ )
  {
    q = q + w[k] / ( 1.0 + x[i] * x[i] ) 
          - w[k+1] * 2.0 * x[i] / pow ( 1.0 + x[i] * x[i], 2 );
    k = k + 2;
  }
  cout << "  Estimate integral of 1/(1+x^2) = " << q << "\n";

  delete [] w;
  delete [] x;

  n = 11;
  a = 0.0;
  b = 1.0;
  x = r8vec_chebyshev_new ( n, a, b );
  w = hermite_interpolant_rule ( n, a, b, x );

  cout << "\n";
  cout << "     I       X               W(F(X))        W(F'(X))\n";
  cout << "\n";
  k = 0;
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(4) << i
         << "  " << setw(14) << x[i]
         << "  " << setw(14) << w[k]
         << "  " << setw(14) << w[k+1] << "\n";
    k = k + 2;
  }

  cout << "\n";
  cout << "  Use the quadrature rule over interval " << a << " to " << b << "\n";
  cout << "\n";

  q = 0.0;
  k = 0;
  for ( i = 0; i < n; i++ )
  {
    q = q + w[k] * 1 + w[k+1] * 0.0;
    k = k + 2;
  }
  cout << "  Estimate integral of 1 = " << q << "\n";

  q = 0.0;
  k = 0;
  for ( i = 0; i < n; i++ )
  {
    q = q + w[k] * x[i] + w[k+1] * 1.0;
    k = k + 2;
  }
  cout << "  Estimate integral of X = " << q << "\n";

  q = 0.0;
  k = 0;
  for ( i = 0; i < n; i++ )
  {
    q = q + w[k] * x[i] * x[i] + w[k+1] * 2.0 * x[i];
    k = k + 2;
  }
  cout << "  Estimate integral of X^2 = " << q << "\n";

  q = 0.0;
  k = 0;
  for ( i = 0; i < n; i++ )
  {
    q = q + w[k] * sin ( x[i] ) + w[k+1] * cos ( x[i] );
    k = k + 2;
  }
  cout << "  Estimate integral of SIN(X) = " << q << "\n";

  delete [] w;
  delete [] x;

  return;
}
//****************************************************************************80

void test08 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST08 tabulates the interpolant and its derivative. 
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 November 2011
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int n;
  int nd;
  int ndp;
  int ns;
  double *x;
  double *xd;
  double *xdp;
  double *xs;
  double *y;
  double *yd;
  double *ydp;
  double *yp;
  double *ys;
  double *ysp;

  cout << "\n";
  cout << "TEST08\n";
  cout << "  HERMITE_INTERPOLANT sets up the Hermite interpolant.\n";
  cout << "  HERMITE_INTERPOLANT_VALUE evaluates it.\n";
  cout << "  Consider data for y=sin(x) at x=0,1,2,3,4.\n";

  n = 5;
  y = new double[n];
  yp = new double[n];

  nd = 2 * n;
  xd = new double[nd];
  yd = new double[nd];

  ndp = 2 * n - 1;
  xdp = new double[ndp];
  ydp = new double[ndp];

  x = r8vec_linspace_new ( n, 0.0, 4.0 );
  for ( i = 0; i < n; i++ )
  {
    y[i] = sin ( x[i] );
    yp[i] = cos ( x[i] );
  }

  hermite_interpolant ( n, x, y, yp, xd, yd, xdp, ydp );
/*
  Now sample the interpolant at NS points, which include data values.
*/
  ns = 4 * ( n - 1 ) + 1;
  ys = new double[ns];
  ysp = new double[ns];

  xs = r8vec_linspace_new ( ns, 0.0, 4.0 );

  hermite_interpolant_value ( nd, xd, yd, xdp, ydp, ns, xs, ys, ysp );

  cout << "\n";
  cout << "  In the following table, there should be perfect\n";
  cout << "  agreement between F and H, and F' and H'\n";
  cout << "  at the data points X = 0, 1, 2, 3, and 4.\n";
  cout << "\n";
  cout << "  In between, H and H' approximate F and F'.\n";
  cout << "\n";
  cout << "     I       X(I)          F(X(I))         H(X(I)) ";
  cout << "        F'(X(I))        H'(X(I))\n";
  cout << "\n";
  for ( i = 0; i < ns; i++ )
  {
    cout << "  " << setw(4) << i
         << "  " << setw(14) << xs[i]
         << "  " << setw(14) << sin ( xs[i] )
         << "  " << setw(14) << ys[i]
         << "  " << setw(14) << cos ( xs[i] )
         << "  " << setw(14) << ysp[i] << "\n";
  }

  delete [] x;
  delete [] xd;
  delete [] xdp;
  delete [] xs;
  delete [] y;
  delete [] yd;
  delete [] ydp;
  delete [] yp;
  delete [] ys;
  delete [] ysp;

  return;
}