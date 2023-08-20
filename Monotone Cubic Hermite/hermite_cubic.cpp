# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>

using namespace std;

# include "hermite_cubic.hpp"

//****************************************************************************80

double hermite_cubic_integral ( double x1, double f1, double d1, double x2,
  double f2, double d2 )

//****************************************************************************80
//
//  Purpose:
//
//    HERMITE_CUBIC_INTEGRAL returns the integral of a Hermite cubic polynomial.
//
//  Discussion:
//
//    The integral is taken over the definition interval [X1,X2].
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 February 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Fred Fritsch, Ralph Carlson,
//    Monotone Piecewise Cubic Interpolation,
//    SIAM Journal on Numerical Analysis,
//    Volume 17, Number 2, April 1980, pages 238-246.
//
//  Parameters:
//
//    Input, double X1, F1, D1, the left endpoint, function value
//    and derivative.
//
//    Input, double X2, F2, D2, the right endpoint, function value
//    and derivative.
//
//    Output, double HERMITE_CUBIC_INTEGRAL, the integral of the Hermite
//    cubic polynomial over the interval X1 <= X <= X2.
//
{
  double h;
  double q;

  h = x2 - x1;

  q = 0.5 * h * ( f1 + f2 + h * ( d1 - d2 ) / 6.0 );

  return q;
}
//****************************************************************************80

double hermite_cubic_integrate ( double x1, double f1, double d1, double x2,
  double f2, double d2, double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    HERMITE_CUBIC_INTEGRATE integrates a Hermite cubic polynomial from A to B.
//
//  Discussion:
//
//    A and B may be scalars, or one may be a vector, or both
//    may be vectors of the same size.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 February 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Fred Fritsch, Ralph Carlson,
//    Monotone Piecewise Cubic Interpolation,
//    SIAM Journal on Numerical Analysis,
//    Volume 17, Number 2, April 1980, pages 238-246.
//
//  Parameters:
//
//    Input, double X1, F1, D1, the left endpoint, function value
//    and derivative.
//
//    Input, double X2, F2, D2, the right endpoint, function value
//    and derivative.
//
//    Input, double A, B, the left and right endpoints of the interval
//    of integration.
//
//    Output, double HERMITE_CUBIC_INTEGRATE, the integral of the
//    Hermite cubic polynomial over the interval A <= X <= B.
//
{
  double dterm;
  double fterm;
  double h;
  double phia1;
  double phia2;
  double phib1;
  double phib2;
  double psia1;
  double psia2;
  double psib1;
  double psib2;
  double q;
  double ta1;
  double ta2;
  double tb1;
  double tb2;
  double ua1;
  double ua2;
  double ub1;
  double ub2;

  h = x2 - x1;

  ta1 = ( a - x1 ) / h;
  ta2 = ( x2 - a ) / h;
  tb1 = ( b - x1 ) / h;
  tb2 = ( x2 - b ) / h;

  ua1 = ta1 * ta1 * ta1;
  phia1 = ua1 * ( 2.0 - ta1 );
  psia1 = ua1 * ( 3.0 * ta1 - 4.0 );

  ua2 = ta2 * ta2 * ta2;
  phia2 =  ua2 * ( 2.0 - ta2 );
  psia2 = -ua2 * ( 3.0 * ta2 - 4.0 );

  ub1 = tb1 * tb1 * tb1;
  phib1 = ub1 * ( 2.0 - tb1 );
  psib1 = ub1 * ( 3.0 * tb1 - 4.0 );

  ub2 = tb2 * tb2 * tb2;
  phib2 =  ub2 * ( 2.0 - tb2 );
  psib2 = -ub2 * ( 3.0 * tb2 - 4.0 );

  fterm =   f1 * ( phia2 - phib2 ) + f2 * ( phib1 - phia1 );
  dterm = ( d1 * ( psia2 - psib2 ) + d2 * ( psib1 - psia1 ) ) * ( h / 6.0 );

  q = 0.5 * h * ( fterm + dterm );

  return q;
}
//****************************************************************************80

double *hermite_cubic_lagrange_integral ( double x1, double x2 )

//****************************************************************************80
//
//  Purpose:
//
//    HERMITE_CUBIC_LAGRANGE_INTEGRAL: Hermite cubic Lagrange integrals.
//
//  Discussion:
//
//    The Hermite cubic polynomial P(X) for interval (X1,X2) and data
//    (F1,D1,F2,D2) satisfies:
//
//      P(X1) = F1,
//      P'(X1) = D1,
//      P(X2) = F2,
//      P'(X2) = D2.
//
//    We can determine four Lagrange polynomials L1(X) through L4(X) so that
//
//      P(X) = F1 * L1(X) + D1 * L2(X) + F2 * L3(X) + D2 * L4(X).
//
//    This function returns the integrals of these four polynomials over
//    the domain of definition [X1,X2].
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 February 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Fred Fritsch, Ralph Carlson,
//    Monotone Piecewise Cubic Interpolation,
//    SIAM Journal on Numerical Analysis,
//    Volume 17, Number 2, April 1980, pages 238-246.
//
//  Parameters:
//
//    Input, double X1, X2, the endpoints.
//
//    Output, double HERMITE_CUBIC_LAGRANGE_INTEGRAL[4], the integrals of the
//    Hermite cubic Lagrange polynomials from X1 to X2.
//
{
  double h;
  double *q;

  q = new double[4];

  h = x2 - x1;

  q[0] =   h     / 2.0;
  q[1] =   h * h / 12.0;
  q[2] =   h     / 2.0;
  q[3] = - h * h / 12.0;

  return q;
}
//****************************************************************************80

double *hermite_cubic_lagrange_integrate ( double x1, double x2, double a,
  double b )

//****************************************************************************80
//
//  Purpose:
//
//    HERMITE_CUBIC_LAGRANGE_INTEGRATE integrates Hermite cubic Lagrange polynomials.
//
//  Discussion:
//
//    A and B may be scalars, or one may be a vector, or both
//    may be vectors of the same size.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 February 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Fred Fritsch, Ralph Carlson,
//    Monotone Piecewise Cubic Interpolation,
//    SIAM Journal on Numerical Analysis,
//    Volume 17, Number 2, April 1980, pages 238-246.
//
//  Parameters:
//
//    Input, double X1, X2, the endpoints of the interval of definition.
//
//    Input, double A, B, the left and right endpoints of the interval
//    of integration.
//
//    Output, double HERMITE_CUBIC_LAGRANGE_INTEGRATE[4], the integrals of the
//    Hermite cubic Lagrange polynomials over the interval A <= X <= B.
//
{
  double h;
  double phia1;
  double phia2;
  double phib1;
  double phib2;
  double psia1;
  double psia2;
  double psib1;
  double psib2;
  double *q;
  double ta1;
  double ta2;
  double tb1;
  double tb2;
  double ua1;
  double ua2;
  double ub1;
  double ub2;

  h = x2 - x1;
  ta1 = ( a - x1 ) / h;
  ta2 = ( x2 - a ) / h;
  tb1 = ( b - x1 ) / h;
  tb2 = ( x2 - b ) / h;

  ua1 = ta1 * ta1 * ta1;
  phia1 = ua1 * ( 2.0 - ta1 );
  psia1 = ua1 * ( 3.0 * ta1 - 4.0 );

  ua2 = ta2 * ta2 * ta2;
  phia2 =  ua2 * ( 2.0 - ta2 );
  psia2 = -ua2 * ( 3.0 * ta2 - 4.0 );

  ub1 = tb1 * tb1 * tb1;
  phib1 = ub1 * ( 2.0 - tb1 );
  psib1 = ub1 * ( 3.0 * tb1 - 4.0 );

  ub2 = tb2 * tb2 * tb2;
  phib2 =  ub2 * ( 2.0 - tb2 );
  psib2 = -ub2 * ( 3.0 * tb2 - 4.0 );

  q = new double[4];

  q[0] = 0.5 * h * ( phia2 - phib2 );
  q[1] = 0.5 * h * ( psia2 - psib2 ) * ( h / 6.0 );
  q[2] = 0.5 * h * ( phib1 - phia1 );
  q[3] = 0.5 * h * ( psib1 - psia1 ) * ( h / 6.0 );

  return q;
}
//****************************************************************************80

void hermite_cubic_lagrange_value ( double x1, double x2, int n, double x[],
  double f[], double d[], double s[], double t[] )

//****************************************************************************80
//
//  Purpose:
//
//    HERMITE_CUBIC_LAGRANGE_VALUE evaluates the Hermite cubic Lagrange polynomials.
//
//  Discussion:
//
//    The Hermite cubic polynomial P(X) for interval (X1,X2) and data
//    (F1,D1,F2,D2) satisfies:
//
//      P(X1) = F1,
//      P'(X1) = D1,
//      P(X2) = F2,
//      P'(X2) = D2.
//
//    We can determine four Lagrange polynomials L1(X) through L4(X) so that
//
//      P(X) = F1 * L1(X) + D1 * L2(X) + F2 * L3(X) + D2 * L4(X).
//
//    This function returns the values and derivatives of these four
//    polynomials.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 February 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Fred Fritsch, Ralph Carlson,
//    Monotone Piecewise Cubic Interpolation,
//    SIAM Journal on Numerical Analysis,
//    Volume 17, Number 2, April 1980, pages 238-246.
//
//  Parameters:
//
//    Input, double X1, X2, the endpoints.
//
//    Input, int N, the number of sample points.
//
//    Input, double X[N], the sample points.
//
//    Output, double F[4*N], D[4*N], S[4*N], T[4*N], the value
//    and first three derivatives of the Hermite cubic Lagrange polynomials
//    at X.
//
{
  double dx;
  double h;
  int j;

  h = x2 - x1;

  for ( j = 0; j < n; j++ )
  {
    dx = x[j] - x1;
//
//  F1.
//
    f[0+j*4] = 1.0 + ( ( dx * dx ) / ( h * h )     ) * ( - 3.0 + ( dx / h ) *  2.0 );
    d[0+j*4] =       ( dx          / ( h * h )     ) * ( - 6.0 + ( dx / h ) *  6.0 );
    s[0+j*4] =       ( 1.0         / ( h * h )     ) * ( - 6.0 + ( dx / h ) * 12.0 );
    t[0+j*4] =       ( 1.0         / ( h * h * h ) )                        * 12.0;
//
//  D1
//
    f[1+j*4] = dx  + ( ( dx * dx ) / h         ) * ( - 2.0 + ( dx / h )       );
    d[1+j*4] = 1.0 + ( dx          / h         ) * ( - 4.0 + ( dx / h ) * 3.0 );
    s[1+j*4] =       ( 1.0         / h         ) * ( - 4.0 + ( dx / h ) * 6.0 );
    t[1+j*4] =       ( 1.0         / ( h * h ) )                        * 6.0;
//
//  F2
//
    f[2+j*4] = ( ( dx * dx ) / ( h * h )     ) * ( 3.0 -  2.0 * ( dx / h ) );
    d[2+j*4] = ( dx          / ( h * h )     ) * ( 6.0 -  6.0 * ( dx / h ) );
    s[2+j*4] = ( 1.0         / ( h * h )     ) * ( 6.0 - 12.0 * ( dx / h ) );
    t[2+j*4] = ( 1.0         / ( h * h * h ) ) * (     - 12.0              );
//
//  D2
//
    f[3+j*4] = ( ( dx * dx ) / h ) * ( - 1.0 + ( dx / h )       );
    d[3+j*4] = ( dx          / h ) * ( - 2.0 + ( dx / h ) * 3.0 );
    s[3+j*4] = ( 1.0         / h ) * ( - 2.0 + ( dx / h ) * 6.0 );
    t[3+j*4] = ( 1.0         / h )                        * 6.0;
  }
  return;
}
//****************************************************************************80

double hermite_cubic_spline_integral ( int nn, double xn[], double fn[],
  double dn[] )

//****************************************************************************80
//
//  Purpose:
//
//    HERMITE_CUBIC_SPLINE_INTEGRAL: Hermite cubic spline integral.
//
//  Discussion:
//
//    The integral is taken over the definition interval ( X[0], X[NN-1] ).
//
//    Note that if the intervals are equal in size, then the derivative
//    information in DN has no effect on the integral value,
//    except for the first and last entries.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 February 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Fred Fritsch, Ralph Carlson,
//    Monotone Piecewise Cubic Interpolation,
//    SIAM Journal on Numerical Analysis,
//    Volume 17, Number 2, April 1980, pages 238-246.
//
//  Parameters:
//
//    Input, int NN, the number of data points.
//
//    Input, double XN[NN], the coordinates of the data points.
//    The entries in XN must be in strictly ascending order.
//
//    Input, double FN[NN], the function values.
//
//    Input, double DN[NN], the derivative values.
//
//    Output, double HERMITE_CUBIC_SPLINE_INTEGRAL, the integral of the
//    Hermite cubic spline over the interval X[0] <= X <= X[NN-1].
//
{
  int i;
  double q;

  q = 0.0;

  for ( i = 0; i < nn - 1; i++ )
  {
    q = q +
    0.5 * ( xn[i+1] - xn[i] ) * ( fn[i] + fn[i+1]
        + ( xn[i+1] - xn[i] ) * ( dn[i] - dn[i+1] ) / 6.0 );
  }
  return q;
}
//****************************************************************************80

double *hermite_cubic_spline_integrate ( int nn, double xn[], double fn[],
  double dn[], int n, double a[], double b[] )

//****************************************************************************80
//
//  Purpose:
//
//    HERMITE_CUBIC_SPLINE_INTEGRATE integrate Hermite cubic spline over [A,B].
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 February 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Fred Fritsch, Ralph Carlson,
//    Monotone Piecewise Cubic Interpolation,
//    SIAM Journal on Numerical Analysis,
//    Volume 17, Number 2, April 1980, pages 238-246.
//
//  Parameters:
//
//    Input, int NN, the number of data points.
//
//    Input, double XN[NN], the coordinates of the data points.
//    The entries in XN must be in strictly ascending order.
//
//    Input, double FN[NN], the function values.
//
//    Input, double DN[NN], the derivative values.
//
//    Input, int N, the number of integration intervals.
//
//    Input, double A[N], B[N], the integration endpoints.
//
//    Output, double HERMITE_CUBIC_SPLINE_INTEGRATE[N], the integral
//    over the interval [A,B].
//
{
  double aa;
  double bb;
  int i;
  int ii;
  int j;
  int k;
  double *q;
  double s;

  q = new double[n];

  i = n / 2;
  j = n / 2;

  for ( ii = 0; ii < n; ii++ )
  {
    q[ii] = 0.0;

    if ( a[ii] <= b[ii] )
    {
      aa = a[ii];
      bb = b[ii];
      s = + 1.0;
    }
    else
    {
      aa = b[ii];
      bb = a[ii];
      s = - 1.0;
    }

    r8vec_bracket3 ( nn, xn, aa, &i );
    r8vec_bracket3 ( nn, xn, bb, &j );
//
//  Evaluate the polynomial with the appropriate data.
//
    if ( i == j )
    {
      q[ii] = hermite_cubic_integrate ( xn[i], fn[i], dn[i],
        xn[i+1], fn[i+1], dn[i+1], aa, bb );
    }
    else
    {
      q[ii] = hermite_cubic_integrate ( xn[i], fn[i], dn[i],
        xn[i+1], fn[i+1], dn[i+1], aa, xn[i+1] );

      for ( k = i + 1; k < j; k++ )
      {
        q[ii] = q[ii] +  hermite_cubic_integral ( xn[k], fn[k], dn[k],
          xn[k+1], fn[k+1], dn[k+1] );
      }
      q[ii] = q[ii] + hermite_cubic_integrate ( xn[j], fn[j], dn[j],
        xn[j+1], fn[j+1], dn[j+1], xn[j], bb );
    }
    q[ii] = s * q[ii];
  }
  return q;
}
//****************************************************************************80

double *hermite_cubic_spline_quad_rule ( int nn, double xn[] )

//****************************************************************************80
//
//  Purpose:
//
//    HERMITE_CUBIC_SPLINE_QUAD_RULE: Hermite cubic spline quadrature rule.
//
//  Discussion:
//
//    The integral is taken over the definition interval ( X[0], X[NN-1] ).
//
//    Note that if the intervals are equal in size, then the derivative
//    information in DN has no effect on the integral value,
//    except for the first and last entries.
//
//    The quadrature rule is
//
//      Integral ( XN[0] <= X <= XN[NN-1] ) F(X) dX is approximately
//
//      Sum ( 0 <= I <= NN-1 ) W[0,I] * F(X[I]) + W[1,I] * F'(X[I])
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 March 2011
//
//  Author:
//
//    John Burkardt.
//
//  Reference:
//
//    Fred Fritsch, Ralph Carlson,
//    Monotone Piecewise Cubic Interpolation,
//    SIAM Journal on Numerical Analysis,
//    Volume 17, Number 2, April 1980, pages 238-246.
//
//  Parameters:
//
//    Input, int NN, the number of data points.
//
//    Input, double XN[NN], the coordinates of the data points.
//    The entries in XN must be in strictly ascending order.
//
//    Output, double W[2*NN], the quadrature weights for F(1:NN)
//    and DN(1:NN).
//
{
  int j;
  double *w;

  w = new double[2*nn];

  w[0+0*2] = 0.5 * ( xn[1] - xn[0] );
  for ( j = 1; j < nn - 1; j++ )
  {
    w[0+j*2] = 0.5 * ( xn[j+1] - xn[j-1] );
  }
  w[0+(nn-1)*2] = 0.5 * ( xn[nn-1] - xn[nn-2]   );

  w[1+0*2] = pow ( xn[1] - xn[0], 2 ) / 12.0;
  for ( j = 1; j < nn - 1; j++ )
  {
    w[1+j*2] = ( xn[j+1] - xn[j-1] )
      * ( xn[j+1] - 2.0 * xn[j] + xn[j-1] ) / 12.0;
  }
  w[1+(nn-1)*2] = - pow ( xn[nn-2] - xn[nn-1], 2 ) / 12.0;

  return w;
}
//****************************************************************************80

void hermite_cubic_spline_value ( int nn, double xn[], double fn[],
  double dn[], int n, double x[], double f[], double d[], double s[],
  double t[] )

//****************************************************************************80
//
//  Purpose:
//
//    HERMITE_CUBIC_SPLINE_VALUE evaluates a Hermite cubic spline.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 February 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Fred Fritsch, Ralph Carlson,
//    Monotone Piecewise Cubic Interpolation,
//    SIAM Journal on Numerical Analysis,
//    Volume 17, Number 2, April 1980, pages 238-246.
//
//  Parameters:
//
//    Input, int NN, the number of data points.
//
//    Input, double XN[NN], the coordinates of the data points.
//    The entries in XN must be in strictly ascending order.
//
//    Input, double FN[NN], the function values.
//
//    Input, double DN[NN], the derivative values.
//
//    Input, int N, the number of sample points.
//
//    Input, double X[N], the coordinates of the sample points.
//
//    Output, double F[N], the function value at the sample points.
//
//    Output, double D[N], the derivative value at the sample points.
//
//    Output, double S[N], the second derivative value at the
//    sample points.
//
//    Output, double T[N], the third derivative value at the
//    sample points.
//
{
  int i;
  int left;

  left = n / 2;

  for ( i = 0; i < n; i++ )
  {
    r8vec_bracket3 ( nn, xn, x[i], &left );

    hermite_cubic_value ( xn[left], fn[left], dn[left], xn[left+1],
      fn[left+1], dn[left+1], 1, x+i, f+i, d+i, s+i, t+i );
  }
  return;
}
//****************************************************************************80

void hermite_cubic_to_power_cubic ( double x1, double f1, double d1, double x2,
  double f2, double d2, double *c0, double *c1, double *c2, double *c3 )

//****************************************************************************80
//
//  Purpose:
//
//    HERMITE_CUBIC_TO_POWER_CUBIC converts a Hermite cubic to power form.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 February 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Fred Fritsch, Ralph Carlson,
//    Monotone Piecewise Cubic Interpolation,
//    SIAM Journal on Numerical Analysis,
//    Volume 17, Number 2, April 1980, pages 238-246.
//
//  Parameters:
//
//    Input, double X1, F1, D1, the left endpoint, function value
//    and derivative.
//
//    Input, double X2, F2, D2, the right endpoint, function value
//    and derivative.
//
//    Output, double *C0, *C1, *C2, *C3, the power form of the polynomial.
//
{
  double df;
  double h;

  h =    x2 - x1;
  df = ( f2 - f1 ) / h;
//
//  Polynomial in terms of X - X1:
//
  *c0 = f1;
  *c1 = d1;
  *c2 = - ( 2.0 * d1 - 3.0 * df + d2 ) / h;
  *c3 =   (       d1 - 2.0 * df + d2 ) / h / h;
//
//  Shift polynomial to X.
//
  *c2 = *c2 - x1 * *c3;
  *c1 = *c1 - x1 * *c2;
  *c0 = *c0 - x1 * *c1;
  *c2 = *c2 - x1 * *c3;
  *c1 = *c1 - x1 * *c2;
  *c2 = *c2 - x1 * *c3;

  return;
}
//****************************************************************************80

void hermite_cubic_value ( double x1, double f1, double d1, double x2,
  double f2, double d2, int n, double x[], double f[], double d[],
  double s[], double t[] )

//****************************************************************************80
//
//  Purpose:
//
//    HERMITE_CUBIC_VALUE evaluates a Hermite cubic polynomial.
//
//  Discussion:
//
//    The input arguments can be vectors.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 February 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Fred Fritsch, Ralph Carlson,
//    Monotone Piecewise Cubic Interpolation,
//    SIAM Journal on Numerical Analysis,
//    Volume 17, Number 2, April 1980, pages 238-246.
//
//  Parameters:
//
//    Input, double X1, F1, D1, the left endpoint, function value
//    and derivative.
//
//    Input, double X2, F2, D2, the right endpoint, function value
//    and derivative.
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], the points at which the Hermite cubic
//    is to be evaluated.
//
//    Output, double F[N], D[N], S[N], T[N], the value and first
//    three derivatives of the Hermite cubic at X.
//
{
  double c2;
  double c3;
  double df;
  double h;
  int i;

  h =    x2 - x1;
  df = ( f2 - f1 ) / h;

  c2 = - ( 2.0 * d1 - 3.0 * df + d2 ) / h;
  c3 =   (       d1 - 2.0 * df + d2 ) / h / h;

  for ( i = 0; i < n; i++ )
  {
    f[i] =       f1 + ( x[i] - x1 ) * ( d1
                    + ( x[i] - x1 ) * ( c2
                    + ( x[i] - x1 ) *   c3 ) );
    d[i] =       d1 + ( x[i] - x1 ) * ( 2.0 * c2
                    + ( x[i] - x1 ) * 3.0 * c3 );
    s[i] = 2.0 * c2 + ( x[i] - x1 ) * 6.0 * c3;
    t[i] = 6.0 * c3;
  }
  return;
}
//****************************************************************************80

void power_cubic_to_hermite_cubic ( double c0, double c1, double c2, double c3,
  double x1, double x2, double *f1, double *d1, double *f2, double *d2 )

//****************************************************************************80
//
//  Purpose:
//
//    POWER_CUBIC_TO_HERMITE_CUBIC converts a power cubic to Hermite form.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 February 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Fred Fritsch, Ralph Carlson,
//    Monotone Piecewise Cubic Interpolation,
//    SIAM Journal on Numerical Analysis,
//    Volume 17, Number 2, April 1980, pages 238-246.
//
//  Parameters:
//
//    Input, double C0, C1, C2, C3, the power form of the
//    polynomial.
//
//    Input, double X1, X2, the left and right endpoints of
//    the Hermite form.
//
//    Output, double *F1, *D1, the function and derivative values at X1.
//
//    Output, double *F2, *D2, the function and derivative values at X2.
//
{
  *f1 = c0 + x1 * ( c1 + x1 * (       c2 + x1       * c3 ) );
  *d1 =             c1 + x1 * ( 2.0 * c2 + x1 * 3.0 * c3 );

  *f2 = c0 + x2 * ( c1 + x2 * (       c2 + x2       * c3 ) );
  *d2 =             c1 + x2 * ( 2.0 * c2 + x2 * 3.0 * c3 );

  return;
}
//****************************************************************************80

double r8_abs ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_ABS returns the absolute value of an R8.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 November 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the quantity whose absolute value is desired.
//
//    Output, double R8_ABS, the absolute value of X.
//
{
  double value;

  if ( 0.0 <= x )
  {
    value = + x;
  }
  else
  {
    value = - x;
  }
  return value;
}
//****************************************************************************80

double r8_uniform_01 ( int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8_UNIFORM_01 returns a unit pseudorandom R8.
//
//  Discussion:
//
//    This routine implements the recursion
//
//      seed = ( 16807 * seed ) mod ( 2^31 - 1 )
//      u = seed / ( 2^31 - 1 )
//
//    The integer arithmetic never requires more than 32 bits,
//    including a sign bit.
//
//    If the initial seed is 12345, then the first three computations are
//
//      Input     Output      R8_UNIFORM_01
//      SEED      SEED
//
//         12345   207482415  0.096616
//     207482415  1790989824  0.833995
//    1790989824  2035175616  0.947702
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Paul Bratley, Bennett Fox, Linus Schrage,
//    A Guide to Simulation,
//    Second Edition,
//    Springer, 1987,
//    ISBN: 0387964673,
//    LC: QA76.9.C65.B73.
//
//    Bennett Fox,
//    Algorithm 647:
//    Implementation and Relative Efficiency of Quasirandom
//    Sequence Generators,
//    ACM Transactions on Mathematical Software,
//    Volume 12, Number 4, December 1986, pages 362-376.
//
//    Pierre L'Ecuyer,
//    Random Number Generation,
//    in Handbook of Simulation,
//    edited by Jerry Banks,
//    Wiley, 1998,
//    ISBN: 0471134031,
//    LC: T57.62.H37.
//
//    Peter Lewis, Allen Goodman, James Miller,
//    A Pseudo-Random Number Generator for the System/360,
//    IBM Systems Journal,
//    Volume 8, Number 2, 1969, pages 136-143.
//
//  Parameters:
//
//    Input/output, int *SEED, the "seed" value.  Normally, this
//    value should not be 0.  On output, SEED has been updated.
//
//    Output, double R8_UNIFORM_01, a new pseudorandom variate,
//    strictly between 0 and 1.
//
{
  int i4_huge = 2147483647;
  int k;
  double r;

  if ( *seed == 0 )
  {
    std::cerr << "\n";
    std::cerr << "R8_UNIFORM_01 - Fatal error!\n";
    std::cerr << "  Input value of SEED = 0.\n";
    std::exit ( 1 );
  }

  k = *seed / 127773;

  *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

  if ( *seed < 0 )
  {
    *seed = *seed + i4_huge;
  }
//
//  Although SEED can be represented exactly as a 32 bit integer,
//  it generally cannot be represented exactly as a 32 bit real number.
//
  r = ( double ) ( *seed ) * 4.656612875E-10;

  return r;
}
//****************************************************************************80

void r8vec_bracket3 ( int n, double t[], double tval, int *left )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_BRACKET3 finds the interval containing or nearest a given value.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//    The routine always returns the index LEFT of the sorted array
//    T with the property that either
//    *  T is contained in the interval [ T[LEFT], T[LEFT+1] ], or
//    *  T < T[LEFT] = T[0], or
//    *  T > T[LEFT+1] = T[N-1].
//
//    The routine is useful for interpolation problems, where
//    the abscissa must be located within an interval of data
//    abscissas for interpolation, or the "nearest" interval
//    to the (extreme) abscissa must be found so that extrapolation
//    can be carried out.
//
//    This version of the function has been revised so that the value of
//    LEFT that is returned uses the 0-based indexing natural to C++.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 April 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, length of the input array.
//
//    Input, double T[N], an array that has been sorted into ascending order.
//
//    Input, double TVAL, a value to be bracketed by entries of T.
//
//    Input/output, int *LEFT.
//    On input, if 0 <= LEFT <= N-2, LEFT is taken as a suggestion for the
//    interval [ T[LEFT-1] T[LEFT] ] in which TVAL lies.  This interval
//    is searched first, followed by the appropriate interval to the left
//    or right.  After that, a binary search is used.
//    On output, LEFT is set so that the interval [ T[LEFT], T[LEFT+1] ]
//    is the closest to TVAL; it either contains TVAL, or else TVAL
//    lies outside the interval [ T[0], T[N-1] ].
//
{
  int high;
  int low;
  int mid;
//
//  Check the input data.
//
  if ( n < 2 )
  {
    std::cerr << "\n";
    std::cerr << "R8VEC_BRACKET3 - Fatal error!\n";
    std::cerr << "  N must be at least 2.\n";
    std::exit ( 1 );
  }
//
//  If *LEFT is not between 0 and N-2, set it to the middle value.
//
  if ( *left < 0 || n - 2 < *left )
  {
    *left = ( n - 1 ) / 2;
  }
//
//  CASE 1: TVAL < T[*LEFT]:
//  Search for TVAL in (T[I],T[I+1]), for I = 0 to *LEFT-1.
//
  if ( tval < t[*left] )
  {
    if ( *left == 0 )
    {
      return;
    }
    else if ( *left == 1 )
    {
      *left = 0;
      return;
    }
    else if ( t[*left-1] <= tval )
    {
      *left = *left - 1;
      return;
    }
    else if ( tval <= t[1] )
    {
      *left = 0;
      return;
    }
//
//  ...Binary search for TVAL in (T[I],T[I+1]), for I = 1 to *LEFT-2.
//
    low = 1;
    high = *left - 2;

    for ( ; ; )
    {
      if ( low == high )
      {
        *left = low;
        return;
      }

      mid = ( low + high + 1 ) / 2;

      if ( t[mid] <= tval )
      {
        low = mid;
      }
      else
      {
        high = mid - 1;
      }
    }
  }
//
//  CASE 2: T[*LEFT+1] < TVAL:
//  Search for TVAL in (T[I],T[I+1]) for intervals I = *LEFT+1 to N-2.
//
  else if ( t[*left+1] < tval )
  {
    if ( *left == n - 2 )
    {
      return;
    }
    else if ( *left == n - 3 )
    {
      *left = *left + 1;
      return;
    }
    else if ( tval <= t[*left+2] )
    {
      *left = *left + 1;
      return;
    }
    else if ( t[n-2] <= tval )
    {
      *left = n - 2;
      return;
    }
//
//  ...Binary search for TVAL in (T[I],T[I+1]) for intervals I = *LEFT+2 to N-3.
//
    low = *left + 2;
    high = n - 3;

    for ( ; ; )
    {

      if ( low == high )
      {
        *left = low;
        return;
      }

      mid = ( low + high + 1 ) / 2;

      if ( t[mid] <= tval )
      {
        low = mid;
      }
      else
      {
        high = mid - 1;
      }
    }
  }
//
//  CASE 3: T[*LEFT] <= TVAL <= T[*LEFT+1]:
//  T is just where the user said it might be.
//
  else
  {
  }

  return;
}
//****************************************************************************80

double *r8vec_even_new ( int n, double alo, double ahi )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_EVEN_NEW returns an R8VEC evenly spaced between ALO and AHI.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 May 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of values.
//
//    Input, double ALO, AHI, the low and high values.
//
//    Output, double R8VEC_EVEN_NEW[N], N evenly spaced values.
//    Normally, A[0] = ALO and A[N-1] = AHI.
//    However, if N = 1, then A[0] = 0.5*(ALO+AHI).
//
{
  double *a;
  int i;

  a = new double[n];

  if ( n == 1 )
  {
    a[0] = 0.5 * ( alo + ahi );
  }
  else
  {
    for ( i = 0; i < n; i++ )
    {
      a[i] = ( ( double ) ( n - i - 1 ) * alo
             + ( double ) (     i     ) * ahi )
             / ( double ) ( n     - 1 );
    }
  }

  return a;
}
//****************************************************************************80

double *r8vec_uniform_01_new ( int n, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_UNIFORM_01_NEW returns a new unit pseudorandom R8VEC.
//
//  Discussion:
//
//    This routine implements the recursion
//
//      seed = ( 16807 * seed ) mod ( 2^31 - 1 )
//      u = seed / ( 2^31 - 1 )
//
//    The integer arithmetic never requires more than 32 bits,
//    including a sign bit.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Paul Bratley, Bennett Fox, Linus Schrage,
//    A Guide to Simulation,
//    Second Edition,
//    Springer, 1987,
//    ISBN: 0387964673,
//    LC: QA76.9.C65.B73.
//
//    Bennett Fox,
//    Algorithm 647:
//    Implementation and Relative Efficiency of Quasirandom
//    Sequence Generators,
//    ACM Transactions on Mathematical Software,
//    Volume 12, Number 4, December 1986, pages 362-376.
//
//    Pierre L'Ecuyer,
//    Random Number Generation,
//    in Handbook of Simulation,
//    edited by Jerry Banks,
//    Wiley, 1998,
//    ISBN: 0471134031,
//    LC: T57.62.H37.
//
//    Peter Lewis, Allen Goodman, James Miller,
//    A Pseudo-Random Number Generator for the System/360,
//    IBM Systems Journal,
//    Volume 8, Number 2, 1969, pages 136-143.
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double R8VEC_UNIFORM_01_NEW[N], the vector of pseudorandom values.
//
{
  int i;
  int i4_huge = 2147483647;
  int k;
  double *r;

  if ( *seed == 0 )
  {
    cerr << "\n";
    cerr << "R8VEC_UNIFORM_01_NEW - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  r = new double[n];

  for ( i = 0; i < n; i++ )
  {
    k = *seed / 127773;

    *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

    if ( *seed < 0 )
    {
      *seed = *seed + i4_huge;
    }

    r[i] = ( double ) ( *seed ) * 4.656612875E-10;
  }

  return r;
}
//****************************************************************************80

void timestamp ( )

//****************************************************************************80
//
//  Purpose:
//
//    TIMESTAMP prints the current YMDHMS date as a time stamp.
//
//  Example:
//
//    31 May 2001 09:45:54 AM
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 July 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    None
//
{
# define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct std::tm *tm_ptr;
  std::time_t now;

  now = std::time ( NULL );
  tm_ptr = std::localtime ( &now );

  std::strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm_ptr );

  std::cout << time_buffer << "\n";

  return;
# undef TIME_SIZE
}