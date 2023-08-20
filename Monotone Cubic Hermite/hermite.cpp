# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <cstring>
# include <ctime>

using namespace std;

# include "hermite.hpp"

//****************************************************************************80

void dif_deriv ( int nd, double xd[], double yd[], int *ndp, double xdp[], 
  double ydp[] )

//****************************************************************************80
//
//  Purpose:
//
//    DIF_DERIV computes the derivative of a polynomial in divided difference form.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 June 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Carl deBoor,
//    A Practical Guide to Splines,
//    Springer, 2001,
//    ISBN: 0387953663,
//    LC: QA1.A647.v27.
//
//  Parameters:
//
//    Input, int ND, the size of the input table.
//
//    Input, double XD[ND], the abscissas for the divided
//    difference table.
//
//    Input, double YD[ND], the divided difference table.
//
//    Output, int *NDP, the size of the output table, which is ND-1.
//
//    Input, double XDP[NP], the abscissas for the divided
//    difference table for the derivative.
//
//    Output, double YDP[NDP], the divided difference
//    table for the derivative.
//
{
  int i;
  double *xd_temp;
  double *yd_temp;
//
//  Using a temporary copy of the difference table, shift the
//  abscissas to zero.
//
  xd_temp = new double[nd];
  yd_temp = new double[nd];

  for ( i = 0; i < nd; i++ )
  {
    xd_temp[i] = xd[i];
  }
  for ( i = 0; i < nd; i++ )
  {
    yd_temp[i] = yd[i];
  }

  dif_shift_zero ( nd, xd_temp, yd_temp );
//
//  Construct the derivative.
//
  *ndp = nd - 1;

  for ( i = 0; i < *ndp; i++ )
  {
    xdp[i] = 0.0;
  }

  for ( i = 0; i < *ndp; i++ )
  {
    ydp[i] = ( double ) ( i + 1 ) * yd_temp[i+1];
  }

  delete [] xd_temp;
  delete [] yd_temp;

  return;
}
//****************************************************************************80

void dif_shift_x ( int nd, double xd[], double yd[], double xv )

//****************************************************************************80
//
//  Purpose:
//
//    DIF_SHIFT_X replaces one abscissa of a divided difference table with a new one.
//
//  Discussion:
//
//    This routine shifts the representation of a divided difference polynomial by
//    dropping the last X value in XD, and adding a new X value to the
//    beginning of the Xd array, suitably modifying the coefficients stored
//    in YD.
//
//    The representation of the polynomial is changed, but the polynomial itself
//    should be identical.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 June 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Carl deBoor,
//    A Practical Guide to Splines,
//    Springer, 2001,
//    ISBN: 0387953663,
//    LC: QA1.A647.v27.
//
//  Parameters:
//
//    Input, int ND, the number of divided difference coefficients, and
//    the number of entries in XD.
//
//    Input/output, double XD[ND], the X values used in the representation of
//    the divided difference polynomial.  After a call to this routine, the 
//    last entry of XD has been dropped, the other
//    entries have shifted up one index, and XV has been inserted at the
//    beginning of the array.
//
//    Input/output, double YD[ND], the divided difference coefficients
//    corresponding to the XD array.  On output, this array has been
//    adjusted.
//
//    Input, double XV, a new X value which is to be used in the representation
//    of the polynomial.  On output, XD[0] equals XV and the representation
//    of the polynomial has been suitably changed.
//    Note that XV does not have to be distinct from any of the original XD
//    values.
//
{
  int i;
//
//  Recompute the divided difference coefficients.
//
  for ( i = nd - 2; 0 <= i; i-- )
  {
    yd[i] = yd[i] + ( xv - xd[i] ) * yd[i+1];
  }
//
//  Shift the X values up one position and insert XV.
//
  for ( i = nd - 1; 0 < i; i-- )
  {
    xd[i] = xd[i-1];
  }

  xd[0] = xv;

  return;
}
//****************************************************************************80

void dif_shift_zero ( int nd, double xd[], double yd[] )

//****************************************************************************80
//
//  Purpose:
//
//    DIF_SHIFT_ZERO shifts a divided difference table so that all abscissas are zero.
//
//  Discussion:
//
//    When the abscissas are changed, the coefficients naturally
//    must also be changed.
//
//    The resulting pair (XD, YD) still represents the
//    same polynomial, but the entries in YD are now the
//    standard polynomial coefficients.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 June 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Carl deBoor,
//    A Practical Guide to Splines,
//    Springer, 2001,
//    ISBN: 0387953663,
//    LC: QA1.A647.v27.
//
//  Parameters:
//
//    Input, int ND, the length of the XD and YD arrays.
//
//    Input/output, double XD[ND], the X values that correspond to the
//    divided difference table.  On output, XD contains only zeroes.
//
//    Input/output, double YD[ND], the divided difference table
//    for the polynomial.  On output, YD is also
//    the coefficient array for the standard representation
//    of the polynomial.
//
{
  int i;
  double xv;

  xv = 0.0;

  for ( i = 1; i <= nd; i++ )
  {
    dif_shift_x ( nd, xd, yd, xv );
  }

  return;
}
//****************************************************************************80

void dif_to_r8poly ( int nd, double xd[], double yd[], double c[] )

//****************************************************************************80
//
//  Purpose:
//
//    DIF_TO_R8POLY converts a divided difference table to a standard polynomial.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 February 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Carl deBoor,
//    A Practical Guide to Splines,
//    Springer, 2001,
//    ISBN: 0387953663,
//    LC: QA1.A647.v27.
//
//  Parameters:
//
//    Input, int ND, the number of coefficients, and abscissas.
//
//    Input, double XD[ND], the X values used in the divided difference
//    representation of the polynomial.
//
//    Input, double YD[ND], the divided difference table.
//
//    Output, double C[ND], the standard form polyomial coefficients.
//    C[0] is the constant term, and C[ND-1] is the coefficient
//    of X^(ND-1).
//
{
  int i;
  int j;

  for ( i = 0; i < nd; i++ )
  {
    c[i] = yd[i];
  }
//
//  Recompute the divided difference coefficients.
//
  for ( j = 1; j <= nd - 1; j++ )
  {
    for ( i = 1; i <= nd - j; i++ )
    {
      c[nd-i-1] = c[nd-i-1] - xd[nd-i-j] * c[nd-i];
    }
  }

  return;
}
//****************************************************************************80

double *dif_vals ( int nd, double xd[], double yd[], int nv, double xv[] )

//****************************************************************************80
//
//  Purpose:
//
//    DIF_VALS evaluates a divided difference polynomial at a set of points.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Carl deBoor,
//    A Practical Guide to Splines,
//    Springer, 2001,
//    ISBN: 0387953663,
//    LC: QA1.A647.v27.
//
//  Parameters:
//
//    Input, int ND, the order of the difference table.
//
//    Input, double XD[ND], the X values of the difference table.
//
//    Input, double YD[ND], the divided differences.
//
//    Input, int NV, the number of evaluation points.
//
//    Input, double XV[NV], the evaluation points.
//
//    Output, double DIF_VALS[NV], the value of the divided difference
//    polynomial at the evaluation points.
//
{
  int i;
  int j;
  double *yv;

  yv = new double[nv];

  for ( j = 0; j < nv; j++ )
  {
    yv[j] = yd[nd-1];
    for ( i = 2; i <= nd; i++ )
    {
      yv[j] = yd[nd-i] + ( xv[j] - xd[nd-i] ) * yv[j];
    }
  }
  return yv;
}
//****************************************************************************80

double hermite_basis_0 ( int n, double x[], int i, double xv )

//****************************************************************************80
//
//  Purpose:
//
//    HERMITE_BASIS_0 evaluates a zero-order Hermite interpolation basis function.
//
//  Discussion:
//
//    Given ND points XD, with values YD and derivative values YPD, the
//    Hermite interpolant can be written as:
//
//      H(X) = sum ( 1 <= I <= ND ) YD(I)  * H0(I;X)
//           + sum ( 1 <= I <= ND ) YPD(I) * H1(I;X)
//
//    where H0(I;X) is the I-th zero order Hermite interpolation basis function,
//    and H1(I;X) is the I-th first order Hermite interpolation basis function.
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
//  Parameters:
//
//    Input, int N, the number of abscissas.
//
//    Input, double X[N], the abscissas.
//
//    Input, int I, the index of the first-order basis function.
//    Indices are 0-based
//
//    Input, double XV, the evaluation point.
//
//    Output, double HERMITE_BASIS_0, the value of the function.
//
{
  double *factor;
  int j;
  double li;
  double lp;
  double lpp;
  double value;

  if ( i < 0 || n - 1 < i )
  {
    cerr << "\n";
    cerr << "HERMITE_BASIS_0 - Fatal error!\n";
    cerr << "  I < 0 or N - 1 < I.\n";
    exit ( 1 );
  }

  factor = new double[n];
//
//  L(X) = product ( X - X(1:N) )
//
//  L'(X(I)).
//
  for ( j = 0; j < n; j++ )
  {
    factor[j] = x[i] - x[j];
  }
  factor[i] = 1.0;

  lp = r8vec_product ( n, factor );
//
//  LI(X) = L(X) / ( X - X(I) ) / L'(X(I))
//
  for ( j = 0; j < n; j++ )
  {
    factor[j] = xv - x[j];
  }
  factor[i] = 1.0;

  li = r8vec_product ( n, factor ) / lp;
//
//  L''(X(I)).
//
  lpp = 0.0;
  for ( j = 0; j < n; j++ )
  {
    factor[j] = x[i] - x[j];
  }
  factor[i] = 1.0;

  for ( j = 0; j < n; j++ )
  {
    if ( j != i )
    {
      factor[j] = 1.0;
      lpp = lpp + 2.0 * r8vec_product ( n, factor );
      factor[j] = x[i] - x[j];
    }
  }

  value = ( 1.0 - ( xv - x[i] ) * lpp / lp ) * li * li;

  delete [] factor;

  return value;
}
//****************************************************************************80

double hermite_basis_1 ( int n, double x[], int i, double xv )

//****************************************************************************80
//
//  Purpose:
//
//    HERMITE_BASIS_1 evaluates a first-order Hermite interpolation basis function.
//
//  Discussion:
//
//    Given ND points XD, with values YD and derivative values YPD, the
//    Hermite interpolant can be written as:
//
//      H(X) = sum ( 1 <= I <= ND ) YD(I)  * H0(I;X)
//           + sum ( 1 <= I <= ND ) YPD(I) * H1(I;X)
//
//    where H0(I;X) is the I-th zero order Hermite interpolation basis function,
//    and H1(I;X) is the I-th first order Hermite interpolation basis function.
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
//  Parameters:
//
//    Input, int N, the number of abscissas.
//
//    Input, double X[N], the abscissas.
//
//    Input, int I, the index of the first-order basis function.
//    Indices are 0-based
//
//    Input, double XV, the evaluation point.
//
//    Output, double VALUE, the value of the function.
//
{
  double bot;
  double *factor;
  int j;
  double top;
  double value;

  if ( i < 0 || n - 1 < i )
  {
    cout << "\n";
    cout << "HERMITE_BASIS_1 - Fatal error!\n";
    cout << "  I < 0 or N - 1 < I.\n";
    exit ( 1 );
  }

  factor = new double[n];

  for ( j = 0; j < n; j++ )
  {
    factor[j] = xv - x[j];
  }
  factor[i] = 1.0;
  top = r8vec_product ( n, factor );

  for ( j = 0; j < n; j++ )
  {
    factor[j] = x[i] - x[j];
  }
  factor[i] = 1.0;
  bot = r8vec_product ( n, factor );

  value = ( xv - x[i] ) * ( top / bot ) * ( top / bot );

  delete [] factor;

  return value;
}
//****************************************************************************80

void hermite_demo ( int n, double x[], double y[], double yp[] )

//****************************************************************************80
//
//  Purpose:
//
//    HERMITE_DEMO computes and prints Hermite interpolant information for data.
//
//  Discussion:
//
//    Given a set of Hermite data, this routine calls HERMITE_INTERPOLANT to 
//    determine and print the divided difference table, and then DIF_TO_R8POLY to 
//    determine and print the coefficients of the polynomial in standard form.
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
//  Parameters:
//
//    Input, int N, the number of data points.
//
//    Input, double X[N], the abscissas.
//
//    Input, double Y[N], YP[N], the function and derivative
//    values at the abscissas.
//
{
  double *cd;
  int i;
  int nd;
  int ndp;
  int nv;
  double *xd;
  double *xdp;
  double *xv;
  double *yd;
  double *ydp;
  double *yv;
  double *yvp;

  cout << "\n";
  cout << "HERMITE_DEMO\n";
  cout << "  Compute coefficients CD of the Hermite polynomial\n";
  cout << "  interpolant to given data (x,y,yp).\n";

  cout << "\n";
  cout << "  Data:\n";
  cout << "              X           Y           Y'\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(4) << i
         << "  " << setw(10) << x[i]
         << "  " << setw(10) << y[i]
         << "  " << setw(10) << yp[i] << "\n";
  }

  nd = 2 * n;
  xd = new double[nd];
  yd = new double[nd];

  ndp = 2 * n - 1;
  xdp = new double[ndp];
  ydp = new double[ndp];

  hermite_interpolant ( n, x, y, yp, xd, yd, xdp, ydp );

  cout << "\n";
  cout << "  Difference table:\n";
  cout << "              XD          YD\n";
  cout << "\n";

  for ( i = 0; i < nd; i++ )
  {
    cout << "  " << setw(4) << i
         << "  " << setw(10) << xd[i]
         << "  " << setw(10) << yd[i] << "\n";
  }

  cout << "\n";
  cout << "  Difference table:\n";
  cout << "              XDP          YDP\n";
  cout << "\n";
  for ( i = 0; i < nd-1; i++ )
  {
    cout << "  " << setw(4) << i
         << "  " << setw(10) << xdp[i]
         << "  " << setw(10) << ydp[i] << "\n";
  }
  cd = new double[2*n];

  dif_to_r8poly ( nd, xd, yd, cd );

  r8poly_print ( nd - 1, cd, "  Hermite interpolating polynomial:" );
//
//  Verify interpolation claim!
//
  nv = n;
  xv = new double[nv];
  yv = new double[nv];
  yvp = new double[nv];

  for ( i = 0; i < nv; i++ )
  {
    xv[i] = x[i];
  }
  hermite_interpolant_value ( nd, xd, yd, xdp, ydp, nv, xv, yv, yvp );

  cout << "\n";
  cout << "  Data Versus Interpolant:\n";
  cout << "              X           Y           H           YP          HP\n";
  cout << "\n";
  for ( i = 0; i < nv; i++ )
  {
    cout << "  " << setw(4) << i
         << "  " << setw(10) << xv[i]
         << "  " << setw(10) << y[i]
         << "  " << setw(10) << yv[i]
         << "  " << setw(10) << yp[i]
         << "  " << setw(10) << yvp[i] << "\n";
  }

  delete [] cd;
  delete [] xd;
  delete [] xdp;
  delete [] xv;
  delete [] yd;
  delete [] ydp;
  delete [] yv;
  delete [] yvp;

  return;
}
//****************************************************************************80

void hermite_interpolant ( int n, double x[], double y[], double yp[], 
  double xd[], double yd[], double xdp[], double ydp[] )

//****************************************************************************80
//
//  Purpose:
//
//    HERMITE_INTERPOLANT sets up a divided difference table from Hermite data.
//
//  Discussion:
//
//    The polynomial represented by the divided difference table can be
//    evaluated by calling DIF_VALS.
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
//  Reference:
//
//    Carl deBoor,
//    A Practical Guide to Splines,
//    Springer, 2001,
//    ISBN: 0387953663,
//    LC: QA1.A647.v27.
//
//  Parameters:
//
//    Input, int N, of items of data 
//    ( X(I), Y(I), YP(I) ).
//
//    Input, double X[N], the abscissas.
//    These values must be distinct.
//
//    Input, double Y[N], YP[N], the function and derivative values.
//
//    Output, double XD[2*N], YD[2*N], the divided difference table
//    for the interpolant value.
//
//    Output, double XDP[2*N-1], YDP[2*N-1], the divided difference 
//    table for the interpolant derivative.
//
{
  int i;
  int j;
  int nd;
  int ndp;
//
//  Copy the data.
//
  nd = 2 * n;

  for ( i = 0; i < n; i++ )
  {
    xd[0+i*2] = x[i];
    xd[1+i*2] = x[i];   
  }
//
//  Carry out the first step of differencing.
//
  yd[0] = y[0];
  for ( i = 1; i < n; i++ )
  {
    yd[0+2*i] = ( y[i] - y[i-1] ) / ( x[i] - x[i-1] );
  }
  for ( i = 0; i < n; i++ )
  {
    yd[1+2*i] = yp[i];
  }
//
//  Carry out the remaining steps in the usual way.
//
  for ( i = 2; i < nd; i++ )
  {
    for ( j = nd - 1; i <= j; j-- )
    {
      yd[j] = ( yd[j] - yd[j-1] ) / ( xd[j] - xd[j-i] );
    }
  }
//
//  Compute the difference table for the derivative.
//
  dif_deriv ( nd, xd, yd, &ndp, xdp, ydp );

  return;
}
//****************************************************************************80

double *hermite_interpolant_rule ( int n, double a, double b, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    HERMITE_INTERPOLANT_RULE: quadrature rule for a Hermite interpolant.
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
//  Parameters:
//
//    Input, int N, the number of abscissas.
//
//    Input, double A, B, the integration limits.
//
//    Input, double X[N], the abscissas.
//
//    Output, double HERMITE_INTERPOLANT_RULE[2*N], the quadrature 
//    coefficients, given as pairs for function and derivative values 
//    at each abscissa.
//
{
  double a_value;
  double b_value;
  double *c;
  int i;
  int k;
  int nd;
  int ndp;
  double *w;
  double *xd;
  double *xdp;
  double *y;
  double *yd;
  double *ydp;
  double *yp;

  y = new double[n];
  yp = new double[n];

  nd = 2 * n;
  c = new double[nd];
  w = new double[nd];
  xd = new double[nd];
  yd = new double[nd];

  ndp = 2 * n - 1;
  xdp = new double[ndp];
  ydp = new double[ndp];

  for ( i = 0; i < n; i++ )
  {
    y[i] = 0.0;
    yp[i] = 0.0;
  }

  k = 0;

  for ( i = 0; i < n; i++ )
  {
    y[i] = 1.0;
    hermite_interpolant ( n, x, y, yp, xd, yd, xdp, ydp );
    dif_to_r8poly ( nd, xd, yd, c );
    a_value = r8poly_ant_val ( n, c, a );
    b_value = r8poly_ant_val ( n, c, b );
    w[k] = b_value - a_value;
    y[i] = 0.0;
    k = k + 1;

    yp[i] = 1.0;
    hermite_interpolant ( n, x, y, yp, xd, yd, xdp, ydp );
    dif_to_r8poly ( nd, xd, yd, c );
    a_value = r8poly_ant_val ( n, c, a );
    b_value = r8poly_ant_val ( n, c, b );
    w[k] = b_value - a_value;
    yp[i] = 0.0;
    k = k + 1;
  }

  delete [] c;
  delete [] xd;
  delete [] xdp;
  delete [] y;
  delete [] yd;
  delete [] ydp;
  delete [] yp;

  return w;
}
//****************************************************************************80

void hermite_interpolant_value ( int nd, double xd[], double yd[], double xdp[], 
  double ydp[], int nv, double xv[], double yv[], double yvp[] )

//****************************************************************************80
//
//  Purpose:
//
//    HERMITE_INTERPOLANT_VALUE evaluates the Hermite interpolant polynomial.
//
//  Discussion:
//
//    In fact, this function will evaluate an arbitrary polynomial that is
//    represented by a difference table.
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
//  Reference:
//
//    Carl deBoor,
//    A Practical Guide to Splines,
//    Springer, 2001,
//    ISBN: 0387953663,
//    LC: QA1.A647.v27.
//
//  Parameters:
//
//    Input, int ND, the order of the difference table.
//
//    Input, double XD[ND], YD[ND], the difference table for the
//    interpolant value.
//
//    Input, double XDP[ND-1], YDP[ND-1], the difference table for
//    the interpolant derivative.
//
//    Input, int NV, the number of evaluation points.
//
//    Input, double XV[NV], the evaluation points.
//
//    Output, double YV[NV], YVP[NV], the value of the interpolant and
//    its derivative at the evaluation points.
//
{
  int i;
  int j;
  int ndp;

  ndp = nd - 1;

  for ( j = 0; j < nv; j++ )
  {
    yv[j] = yd[nd-1];
    for ( i = nd - 2; 0 <= i; i-- )
    {
      yv[j] = yd[i] + ( xv[j] - xd[i] ) * yv[j];
    }

    yvp[j] = ydp[ndp-1];
    for ( i = ndp - 2; 0 <= i; i-- )
    {
      yvp[j] = ydp[i] + ( xv[j] - xdp[i] ) * yvp[j];
    }
  }
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

double r8_max ( double x, double y )

//****************************************************************************80
//
//  Purpose:
//
//    R8_MAX returns the maximum of two R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, Y, the quantities to compare.
//
//    Output, double R8_MAX, the maximum of X and Y.
//
{
  double value;

  if ( y < x )
  {
    value = x;
  }
  else
  {
    value = y;
  }
  return value;
}
//****************************************************************************80

double r8poly_ant_val ( int n, double poly_cof[], double xval )

//****************************************************************************80
//
//  Purpose:
//
//    R8POLY_ANT_VAL evaluates the antiderivative of an R8POLY in standard form.
//
//  Discussion:
//
//    The constant term of the antiderivative is taken to be zero.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the polynomial.
//
//    Input, double POLY_COF[N], the polynomial coefficients.  POLY_COF[0]
//    is the constant term, and POLY_COF[N-1] is the coefficient of X**(N-1).
//
//    Input, double XVAL, the point where the antiderivative is to be
//    evaluated.
//
//    Output, double R8POLY_ANT_VAL, the value of the antiderivative of the polynomial
//    at XVAL.
//
{
  int i;
  double value;

  value = 0.0;

  for ( i = n - 1; 0 <= i; i-- )
  {
    value = ( value + poly_cof[i] / ( double ) ( i + 1 ) ) * xval;
  }

  return value;
}
//****************************************************************************80

int r8poly_degree ( int na, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8POLY_DEGREE returns the degree of a polynomial.
//
//  Discussion:
//
//    The degree of a polynomial is the index of the highest power
//    of X with a nonzero coefficient.
//
//    The degree of a constant polynomial is 0.  The degree of the
//    zero polynomial is debatable, but this routine returns the
//    degree as 0.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 September 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NA, the dimension of A.
//
//    Input, double A[NA+1], the coefficients of the polynomials.
//
//    Output, int R8POLY_DEGREE, the degree of A.
//
{
  int degree;

  degree = na;

  while ( 0 < degree )
  {
    if ( a[degree] != 0.0 )
    {
      return degree;
    }
    degree = degree - 1;
  }

  return degree;
}
//****************************************************************************80

void r8poly_print ( int n, double a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8POLY_PRINT prints out a polynomial.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 September 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the dimension of A.
//
//    Input, double A[N+1], the polynomial coefficients.
//    A(0) is the constant term and
//    A(N) is the coefficient of X**N.
//
//    Input, string TITLE, a title.
//
{
  int i;
  double mag;
  int n2;
  char plus_minus;

  cout << "\n";
  cout << title << "\n";
  cout << "\n";

  n2 = r8poly_degree ( n, a );

  if ( n2 <= 0 )
  {
    cout << "  p(x) = 0\n";
    return;
  }

  if ( a[n2] < 0.0 )
  {
    plus_minus = '-';
  }
  else
  {
    plus_minus = ' ';
  }

  mag = r8_abs ( a[n2] );

  if ( 2 <= n2 )
  {
    cout << "  p(x) = " << plus_minus
         << setw(14) << mag << " * x ^ " << n2 << "\n";
  }
  else if ( n2 == 1 )
  {
    cout << "  p(x) = " << plus_minus
         << setw(14) << mag << " * x\n";
  }
  else if ( n2 == 0 )
  {
    cout << "  p(x) = " << plus_minus
         << setw(14) << mag << "\n";
  }

  for ( i = n2-1; 0 <= i; i-- )
  {
    if ( a[i] < 0.0 )
    {
      plus_minus = '-';
    }
    else
    {
      plus_minus = '+';
    }

    mag = r8_abs ( a[i] );

    if ( mag != 0.0 )
    {
      if ( 2 <= i )
      {
        cout << "         " << plus_minus
             << setw(14) << mag << " * x ^ " << i << "\n";
      }
      else if ( i == 1 )
      {
        cout << "         " << plus_minus
             << setw(14) << mag << " * x\n";
      }
      else if ( i == 0 )
      {
        cout << "         " << plus_minus
             << setw(14) << mag << "\n";
      }
    }
  }

  return;
}
//****************************************************************************80

double *r8vec_chebyshev_new ( int n, double a_first, double a_last )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_CHEBYSHEV_NEW creates a vector of Chebyshev spaced values.
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
//    08 June 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input, double A_FIRST, A_LAST, the first and last entries.
//
//    Output, double R8VEC_CHEBYSHEV_NEW[N], a vector of Chebyshev spaced data.
//
{
  double *a;
  double c;
  int i;
  double pi = 3.141592653589793;
  double theta;

  a = new double[n];

  if ( n == 1 )
  {
    a[0] = ( a_first + a_last ) / 2.0;
  }
  else
  {
    for ( i = 0; i < n; i++ )
    {
      theta = ( double ) ( n - i - 1 ) * pi / ( double ) ( n - 1 );

      c = cos ( theta );

      if ( ( n % 2 ) == 1 )
      {
        if ( 2 * i + 1 == n )
        {
          c = 0.0;
        }
      }

      a[i] = ( ( 1.0 - c ) * a_first  
             + ( 1.0 + c ) * a_last ) 
             /   2.0;

    }
  }

  return a;
}
//****************************************************************************80

double *r8vec_linspace_new ( int n, double a_first, double a_last )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_LINSPACE_NEW creates a vector of linearly spaced values.
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
//    29 March 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input, double A_FIRST, A_LAST, the first and last entries.
//
//    Output, double R8VEC_LINSPACE_NEW[N], a vector of linearly spaced data.
//
{
  double *a;
  int i;

  a = new double[n];

  if ( n == 1 )
  {
    a[0] = ( a_first + a_last ) / 2.0;
  }
  else
  {
    for ( i = 0; i < n; i++ )
    {
      a[i] = ( ( double ) ( n - 1 - i ) * a_first 
             + ( double ) (         i ) * a_last ) 
             / ( double ) ( n - 1     );
    }
  }
  return a;
}
//****************************************************************************80

void r8vec_print ( int n, double a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_PRINT prints an R8VEC.
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
//    16 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of components of the vector.
//
//    Input, double A[N], the vector to be printed.
//
//    Input, string TITLE, a title.
//
{
  int i;

  cout << "\n";
  cout << title << "\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(8)  << i
         << ": " << setw(14) << a[i]  << "\n";
  }

  return;
}
//****************************************************************************80

double r8vec_product ( int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_PRODUCT returns the product of the entries of an R8VEC.
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
//    17 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input, double A[N], the vector.
//
//    Output, double R8VEC_PRODUCT, the product of the vector.
//
{
  int i;
  double product;

  product = 1.0;
  for ( i = 0; i < n; i++ )
  {
    product = product * a[i];
  }

  return product;
}
//****************************************************************************80

void r8vec_uniform_01 ( int n, int *seed, double r[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_UNIFORM_01 returns a unit pseudorandom R8VEC.
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
//    Output, double R[N], the vector of pseudorandom values.
//
{
  int i;
  static int i4_huge = 2147483647;
  int k;

  if ( *seed == 0 )
  {
    cerr << "\n";
    cerr << "R8VEC_UNIFORM_01 - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

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

  return;
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