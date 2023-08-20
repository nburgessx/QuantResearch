# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

using namespace std;

# include "hermite_cubic.hpp"

int main ( );
void test01 ( );
void test02 ( );
void test03 ( );
void test04 ( );
void test05 ( );
void test06 ( );
void test07 ( );
void test08 ( );
void test09 ( );
void test10 ( );
void test11 ( );
void test12 ( );
void test13 ( );
void test14 ( );
void test15 ( );
double cubic_antiderivative ( double x );
double cubic_integrate ( double a, double b );
void cubic_value ( double x, double *f, double *d, double *s, double *t );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for HERMITE_CUBIC_TEST.
//
//  Discussion:
//
//    HERMITE_CUBIC_TEST tests the HERMITE_CUBIC library.
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
{
  timestamp ( );
  std::cout << "\n";
  std::cout << "HERMITE_CUBIC_TEST\n";
  std::cout << "  C++ version\n";
  std::cout << "  Test the HERMITE_CUBIC library.\n";

  test01 ( );
  test02 ( );
  test03 ( );
  test04 ( );
  test05 ( );
  test06 ( );
  test07 ( );
  test08 ( );
  test09 ( );
  test10 ( );
  test11 ( );
  test12 ( );
  test13 ( );
  test14 ( );
  test15 ( );
//
//  Terminate.
//
  std::cout << "\n";
  std::cout << "HERMITE_CUBIC_TEST\n";
  std::cout << "  Normal end of execution.\n";
  std::cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void test01 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST01 tests HERMITE_CUBIC_VALUE.
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
{
  double d[1];
  double d1;
  double d2;
  double f[1];
  double f1;
  double f2;
  int i;
  int j;
  int n = 1;
  double s[1];
  double t[1];
  double x[1];
  int x_interval;
  double x1;
  double x2;

  std::cout << "\n";
  std::cout << "TEST01:\n";
  std::cout << "  HERMITE_CUBIC_VALUE evaluates a Hermite cubic polynomial.\n";
  std::cout << "  Try out four sets of data:\n";
  std::cout << "  (F1,D1,F2,D2) = (1,0,0,0), (0,1,0,0), (0,0,1,0), (0,0,0,1)\n";
  std::cout << "  on [0,1] and [1.0,-2.0] (interval reversed)\n";

  for ( x_interval = 1; x_interval <= 2; x_interval++ )
  {
    if ( x_interval == 1 )
    {
      x1 = 0.0;
      x2 = 1.0;
    }
    else
    {
      x1 = 1.0;
      x2 = -2.0;
    }

    for ( i = 1; i <= 4; i++ )
    {
      f1 = 0.0;
      d1 = 0.0;
      f2 = 0.0;
      d2 = 0.0;

      if ( i == 1 )
      {
        f1 = 1.0;
      }
      else if ( i == 2 )
      {
        d1 = 1.0;
      }
      else if ( i == 3 )
      {
        f2 = 1.0;
      }
      else if ( i == 4 )
      {
        d2 = 1.0;
      }

      std::cout << "\n";
      std::cout << "    J      X           F           D\n";
      std::cout << "\n";

      for ( j = -3; j <= 12; j++ )
      {
        x[0] = ( ( double ) ( 10 - j ) * x1
               + ( double ) (      j ) * x2 )
               / ( double ) ( 10     );

        hermite_cubic_value ( x1, f1, d1, x2, f2, d2, n, x, f, d, s, t );

        if ( j == 0 )
        {
          std::cout << "*Data"
               << "  " << setw(10) << x1
               << "  " << setw(10) << f1
               << "  " << setw(10) << d1 << "\n";
        }
        std::cout << "  " << setw(3) << j
             << "  " << setw(10) << x[0]
             << "  " << setw(10) << f[0]
             << "  " << setw(10) << d[0] << "\n";
        if ( j == 10 )
        {
          std::cout << "*Data"
               << "  " << setw(10) << x2
               << "  " << setw(10) << f2
               << "  " << setw(10) << d2 << "\n";
        }
      }
    }
  }
  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests HERMITE_CUBIC_VALUE.
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
{
  double d[1];
  double dc;
  double d1;
  double d2;
  double f[1];
  double fc;
  double f1;
  double f2;
  int j;
  int n = 1;
  double s[1];
  double s1;
  double s2;
  double sc;
  double t[1];
  double t1;
  double t2;
  double tc;
  double x[1];
  int x_interval;
  double x1;
  double x2;

  std::cout << "\n";
  std::cout << "TEST02:\n";
  std::cout << "  HERMITE_CUBIC_VALUE evaluates a Hermite cubic polynomial.\n";
  std::cout << "  Try out data from a cubic function:\n";
  std::cout << "  on [0,10] and [-1.0,1.0] and [0.5,0.75]\n";

  for ( x_interval = 1; x_interval <= 3; x_interval++ )
  {
    if ( x_interval == 1 )
    {
      x1 = 0.0;
      x2 = 10.0;
    }
    else if ( x_interval == 2 )
    {
      x1 = -1.0;
      x2 = +1.0;
    }
    else if ( x_interval == 3 )
    {
      x1 = 0.5;
      x2 = 0.75;
    }

    cubic_value ( x1, &f1, &d1, &s1, &t1 );
    cubic_value ( x2, &f2, &d2, &s2, &t2 );

    std::cout << "\n";
    std::cout << "    J      X           F           D           S           T\n";
    std::cout << "\n";

    for ( j = -3; j <= 12; j++ )
    {
      x[0] = ( ( double ) ( 10 - j ) * x1
             + ( double ) (      j ) * x2 )
             / ( double ) ( 10 );

      hermite_cubic_value ( x1, f1, d1, x2, f2, d2, n, x, f, d, s, t );
      cubic_value ( x[0], &fc, &dc, &sc, &tc );

      if ( j == 0 )
      {
        std::cout << "*Data"
             << "  " << setw(10) << x1
             << "  " << setw(10) << f1
             << "  " << setw(10) << d1 << "\n";
      }
      std::cout << "Exact"
           << "  " << setw(10) << x[0]
           << "  " << setw(10) << fc
           << "  " << setw(10) << dc
           << "  " << setw(10) << sc
           << "  " << setw(10) << tc << "\n";
      std::cout << "  " << setw(3) << j
           << "  " << setw(10) << x[0]
           << "  " << setw(10) << f[0]
           << "  " << setw(10) << d[0]
           << "  " << setw(10) << s[0]
           << "  " << setw(10) << t[0] << "\n";
      if ( j == 10 )
      {
        std::cout << "*Data"
             << "  " << setw(10) << x2
             << "  " << setw(10) << f2
             << "  " << setw(10) << d2 << "\n";
      }
    }
  }
  return;
}
//****************************************************************************80

void test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 tests HERMITE_CUBIC_INTEGRATE.
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
{
  double a;
  double b;
  double d1;
  double d2;
  double f1;
  double f2;
  int j;
  double q_computed;
  double q_exact;
  double s1;
  double s2;
  double t1;
  double t2;
  int x_interval;
  double x1;
  double x2;

  std::cout << "\n";
  std::cout << "TEST03:\n";
  std::cout << "  HERMITE_CUBIC_INTEGRATE integrates a Hermite cubic\n";
  std::cout << "  polynomial from A to B.\n";

  for ( x_interval = 1; x_interval <= 3; x_interval++ )
  {
    if ( x_interval == 1 )
    {
      x1 = 0.0;
      x2 = 10.0;
    }
    else if ( x_interval == 2 )
    {
      x1 = -1.0;
      x2 = +1.0;
    }
    else if ( x_interval == 3 )
    {
      x1 = 0.5;
      x2 = 0.75;
    }

    cubic_value ( x1, &f1, &d1, &s1, &t1 );
    cubic_value ( x2, &f2, &d2, &s2, &t2 );

    std::cout << "\n";
    std::cout << "                                     Exact           Computed\n";
    std::cout << "    J          A           B         Integral        Integral\n";
    std::cout << "\n";

    a = x1 - 1.0;

    for ( j = - 3; j <= 12; j++ )
    {
      b = ( ( double ) ( 10 - j ) * x1
          + ( double ) (      j ) * x2 )
          / ( double ) ( 10     );

      q_exact = cubic_integrate ( a, b );

      q_computed = hermite_cubic_integrate ( x1, f1, d1, x2, f2, d2, a, b );

      std::cout << "  " << setw(3) << j
           << "  " << setw(10) << a
           << "  " << setw(10) << b
           << "  " << setw(14) << q_exact
           << "  " << setw(14) << q_computed << "\n";
    }
  }
  return;
}
//****************************************************************************80

void test04 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST04 tests HERMITE_CUBIC_SPLINE_VALUE.
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
{
  double *d;
  double *dn;
  double *f;
  double *fn;
  int i;
  int n = 51;
  int nn = 11;
  double *s;
  double *t;
  double u;
  double v;
  double x1;
  double x2;
  double *x;
  double *xn;

  std::cout << "\n";
  std::cout << "TEST04:\n";
  std::cout << "  HERMITE_CUBIC_SPLINE_VALUE evaluates a Hermite cubic spline.\n";

  x1 = 0.0;
  x2 = 10.0;

  xn = r8vec_even_new ( nn, x1, x2 );
  fn = new double[nn];
  dn = new double[nn];

  for ( i = 0; i < nn; i++ )
  {
    fn[i] = std::sin ( xn[i] );
    dn[i] = std::cos ( xn[i] );
  }

  x = r8vec_even_new ( n, x1, x2 );
  f = new double[n];
  d = new double[n];
  s = new double[n];
  t = new double[n];

  hermite_cubic_spline_value ( nn, xn, fn, dn, n, x, f, d, s, t );

  std::cout << "\n";
  std::cout << "     I      X       F computed     F exact      Error\n";
  std::cout << "\n";

  for ( i = 0; i < n; i++ )
  {
    u = std::sin ( x[i] );
    v = r8_abs ( f[i] - u );
    std::cout << "  " << setw(4) << i
         << "  " << setw(10) << x[i]
         << "  " << setw(10) << f[i]
         << "  " << setw(10) << u
         << "  " << setw(14) << v << "\n";
  }

  std::cout << "\n";
  std::cout << "     I      X       D computed     D exact      Error\n";
  std::cout << "\n";

  for ( i = 0; i < n; i++ )
  {
    u = std::cos ( x[i] );
    v = r8_abs ( d[i] - u );
    std::cout << "  " << setw(4) << i
         << "  " << setw(10) << x[i]
         << "  " << setw(10) << d[i]
         << "  " << setw(10) << u
         << "  " << setw(14) << v << "\n";
  }

  delete [] d;
  delete [] dn;
  delete [] f;
  delete [] fn;
  delete [] s;
  delete [] t;
  delete [] x;
  delete [] xn;

  return;
}
//****************************************************************************80

void test05 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST05 tests HERMITE_CUBIC_TO_POWER_CUBIC
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
{
  double c0;
  double c1;
  double c2;
  double c3;
  double d[1];
  double d1;
  double d1r;
  double d2;
  double d2r;
  double f[1];
  double f1;
  double f1r;
  double f2;
  double f2r;
  double fp;
  int j;
  int n = 1;
  double s[1];
  double s1;
  double s2;
  double t[1];
  double t1;
  double t2;
  double x[1];
  double x1;
  double x2;

  std::cout << "\n";
  std::cout << "TEST05:\n";
  std::cout << "  HERMITE_CUBIC_TO_POWER_CUBIC converts the Hermite data\n";
  std::cout << "  to the coefficients of the power form of the polynomial\n";
  std::cout << "  POWER_CUBIC_TO_HERMITE_CUBIC converts the power form\n";
  std::cout << "  to Hermite form\n";

  x1 = -1.0;
  x2 = +1.0;

  cubic_value ( x1, &f1, &d1, &s1, &t1 );
  cubic_value ( x2, &f2, &d2, &s2, &t2 );

  std::cout << "\n";
  std::cout << "  Hermite data:\n";
  std::cout << "\n";
  std::cout << "  X1, F1, D1:" << setw(10) << x1
       << "  " << setw(10) << f1
       << "  " << setw(10) << d1 << "\n";
  std::cout << "  X2, F2, D2:" << setw(10) << x2
       << "  " << setw(10) << f2
       << "  " << setw(10) << d2 << "\n";

  hermite_cubic_to_power_cubic ( x1, f1, d1, x2, f2, d2, &c0, &c1, &c2, &c3 );

  std::cout << "\n";
  std::cout << "  Power form:\n";
  std::cout << "  p(x) = " << c0 << " + " << c1 << " * x + "
       << c2 << " * x^2 + " << c3 << " * x^3\n";
  std::cout << "\n";
  std::cout << "      X       F (Hermite)  F (power)\n";
  std::cout << "\n";

  for ( j = -3; j <= 12; j++ )
  {
    x[0] = ( ( double ) ( 10 - j ) * x1
           + ( double ) (      j ) * x2 )
           / ( double ) ( 10     );

    hermite_cubic_value ( x1, f1, d1, x2, f2, d2, n, x, f+0, d+0, s+0, t+0 );

    fp = c0 + x[0] * ( c1 + x[0] * ( c2 + x[0] * c3 ) );

    std::cout << "  " << setw(10) << x[0]
         << "  " << setw(10) << f[0]
         << "  " << setw(10) << fp << "\n";
  }

  power_cubic_to_hermite_cubic ( c0, c1, c2, c3, x1, x2, &f1r, &d1r,
    &f2r, &d2r );

  std::cout << "\n";
  std::cout << "  Use POWER_CUBIC_TO_HERMITE_CUBIC to recover the\n";
  std::cout << "  original Hermite data:\n";
  std::cout << "\n";
  std::cout << "         Original   Recovered\n";
  std::cout << "\n";
  std::cout << "  F1:  " << "  " << setw(10) << f1
                    << "  " << setw(10) << f1r << "\n";
  std::cout << "  D1:  " << "  " << setw(10) << d1
                    << "  " << setw(10) << d1r << "\n";
  std::cout << "  F2:  " << "  " << setw(10) << f2
                    << "  " << setw(10) << f2r << "\n";
  std::cout << "  D2:  " << "  " << setw(10) << d2
                    << "  " << setw(10) << d2r << "\n";
  return;
}
//****************************************************************************80

void test06 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST06 tests HERMITE_CUBIC_INTEGRATE using vectors.
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
{
  double a;
  double b;
  double d1;
  double d2;
  double f1;
  double f2;
  int i;
  double q_computed;
  double q_exact;
  double s1;
  double s2;
  double t1;
  double t2;
  double x1;
  double x2;

  std::cout << "\n";
  std::cout << "TEST06:\n";
  std::cout << "  HERMITE_CUBIC_INTEGRATE integrates a Hermite cubic\n";
  std::cout << "  polynomial from A to B.\n";
  std::cout << "  Use A, B vectors for the calculation.\n";

  x1 = 0.0;
  x2 = 10.0;

  cubic_value ( x1, &f1, &d1, &s1, &t1 );
  cubic_value ( x2, &f2, &d2, &s2, &t2 );

  std::cout << "\n";
  std::cout << "                                 Exact       Computed\n";
  std::cout << "    J      A           B         Integral    Integral\n";
  std::cout << "\n";

  for ( i = -3; i <= 12; i++ )
  {
    a = x1 - 1.0;
    b = ( ( double ) ( 10 - i ) * x1
        + ( double ) (      i ) * x2 )
        / ( double ) ( 10     );

    q_exact = cubic_integrate ( a, b );

    q_computed = hermite_cubic_integrate ( x1, f1, d1, x2, f2, d2, a, b );

    std::cout << "  " << setw(3) << i
         << "  " << setw(10) << a
         << "  " << setw(10) << b
         << "  " << setw(14) << q_exact
         << "  " << setw(14) << q_computed << "\n";
  }
  return;
}
//****************************************************************************80

void test07 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST07 tests HERMITE_CUBIC_INTEGRAL.
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
{
  double d1;
  double d2;
  double f1;
  double f2;
  double q_computed;
  double q_exact;
  double s1;
  double s2;
  double t1;
  double t2;
  int x_interval;
  double x1;
  double x2;

  std::cout << "\n";
  std::cout << "TEST07:\n";
  std::cout << "  HERMITE_CUBIC_INTEGRAL integrates a Hermite cubic\n";
  std::cout << "  polynomial over the definition interval [X1,X2].\n";
  std::cout << "\n";
  std::cout << "                            Exact       Computed\n";
  std::cout << "     X1          X2         Integral    Integral\n";
  std::cout << "\n";

  for ( x_interval = 1; x_interval <= 3; x_interval++ )
  {
    if ( x_interval == 1 )
    {
      x1 = 0.0;
      x2 = 10.0;
    }
    else if ( x_interval == 2 )
    {
      x1 = -1.0;
      x2 = +1.0;
    }
    else if ( x_interval == 3 )
    {
      x1 = 0.5;
      x2 = 0.75;
    }

    cubic_value ( x1, &f1, &d1, &s1, &t1 );
    cubic_value ( x2, &f2, &d2, &s2, &t2 );

    q_exact = cubic_integrate ( x1, x2 );

    q_computed = hermite_cubic_integral ( x1, f1, d1, x2, f2, d2 );

    std::cout << "  " << setw(10) << x1
         << "  " << setw(10) << x2
         << "  " << setw(14) << q_exact
         << "  " << setw(14) << q_computed << "\n";
  }
  return;
}
//****************************************************************************80

void test08 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST08 tests HERMITE_CUBIC_SPLINE_INTEGRAL.
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
{
  double a;
  double b;
  double *dn;
  double *fn;
  int i;
  int nn = 11;
  double pi = 3.141592653589793;
  double q_computed;
  double q_exact;
  int test;
  double *xn;

  std::cout << "\n";
  std::cout << "TEST08:\n";
  std::cout << "  HERMITE_CUBIC_SPLINE_INTEGRAL integrates a Hermite\n";
  std::cout << "  cubic spline over the definition interval [X1,XNN].\n";
  std::cout << "\n";
  std::cout << "                            Exact       Computed\n";
  std::cout << "     X1          XNN        Integral    Integral\n";
  std::cout << "\n";

  fn = new double[nn];
  dn = new double[nn];

  for ( test = 1; test <= 3; test++ )
  {
    if ( test == 1 )
    {
      a = 0.0;
      b = 1.0;

      xn = r8vec_even_new ( nn, a, b );
      for ( i = 0; i < nn; i++ )
      {
        fn[i] = xn[i] * ( 4.0 * xn[i] - 1.0 ) * ( xn[i] - 1.0 );
        dn[i] = 1.0 + xn[i] * ( - 10.0 + xn[i] * 12.0 );
      }
      q_exact =
        ( xn[nn-1] * xn[nn-1] * ( 0.5 + xn[nn-1] * ( - ( 5.0 / 3.0 ) + xn[nn-1] ) ) )
      - ( xn[0]  * xn[0]  * ( 0.5 + xn[0]  * ( - ( 5.0 / 3.0 ) + xn[0]  ) ) );
    }
//
//  Use variable spacing.
//
    else if ( test == 2 )
    {
      a = 0.0;
      b = 1.0;

      xn = r8vec_even_new ( nn, a, b );
      for ( i = 0; i < nn; i++ )
      {
        xn[i] = std::sqrt ( xn[i] );
        fn[i] = xn[i] * ( 4.0 * xn[i] - 1.0 ) * ( xn[i] - 1.0 );
        dn[i] = 1.0 + xn[i] * ( - 10.0 + xn[i] * 12.0 );
      }
      q_exact =
        ( xn[nn-1] * xn[nn-1] * ( 0.5 + xn[nn-1] * ( - ( 5.0 / 3.0 ) + xn[nn-1] ) ) )
      - ( xn[0]  * xn[0]  * ( 0.5 + xn[0]  * ( - ( 5.0 / 3.0 ) + xn[0]  ) ) );
    }
//
//  Try a non-cubic.
//
    else if ( test == 3 )
    {
      a = 0.0;
      b = pi;

      xn = r8vec_even_new ( nn, a, b );
      for ( i = 0; i < nn; i++ )
      {
        fn[i] = std::sin ( xn[i] );
        dn[i] = std::cos ( xn[i] );
      }
      q_exact = - std::cos ( xn[nn-1] ) + std::cos ( xn[0] );
    }

    q_computed = hermite_cubic_spline_integral ( nn, xn, fn, dn );

    std::cout << "  " << setw(10) << xn[0]
         << "  " << setw(10) << xn[nn-1]
         << "  " << setw(14) << q_exact
         << "  " << setw(14) << q_computed << "\n";

    delete [] xn;
  }
  delete [] dn;
  delete [] fn;

  return;
}
//****************************************************************************80

void test09 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST09 tests HERMITE_CUBIC_SPLINE_INTEGRATE.
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
{
  double *a;
  double *b;
  double *dn;
  double *fn;
  int i;
  int n = 25;
  int nn = 11;
  double *q;
  double q_exact;
  double *sn;
  double *tn;
  double x1;
  double x2;
  double *xn;

  std::cout << "\n";
  std::cout << "TEST09:\n";
  std::cout << "  HERMITE_CUBIC_SPLINE_INTEGRATE integrates a Hermite\n";
  std::cout << "  cubic spline from A to B.\n";
//
//  Define the cubic spline.
//
  x1 = 0.0;
  x2 = 10.0;

  xn = r8vec_even_new ( nn, x1, x2 );
  fn = new double[nn];
  dn = new double[nn];
  sn = new double[nn];
  tn = new double[nn];

  for ( i = 0; i < nn; i++ )
  {
    cubic_value ( xn[i], fn+i, dn+i, sn+i, tn+i );
  }

  a = new double[n];
  for ( i = 0; i < n; i++ )
  {
    a[i] = 2.5;
  }
  b = r8vec_even_new ( n, x1 - 1.0, x2 + 1.0 );

  q = hermite_cubic_spline_integrate ( nn, xn, fn, dn, n, a, b );

  std::cout << "\n";
  std::cout << "                                 Exact       Computed\n";
  std::cout << "    I      A           B         Integral    Integral\n";
  std::cout << "\n";

  for  ( i = 0; i < n; i++ )
  {
    q_exact = cubic_integrate ( a[i], b[i] );

    std::cout << "  " << setw(3) << i
         << "  " << setw(10) << a[i]
         << "  " << setw(10) << b[i]
         << "  " << setw(10) << q_exact
         << "  " << setw(10) << q[i] << "\n";
  }

  delete [] a;
  delete [] b;
  delete [] dn;
  delete [] fn;
  delete [] q;
  delete [] sn;
  delete [] tn;
  delete [] xn;

  return;
}
//****************************************************************************80

void test10 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST10 tests HERMITE_CUBIC_SPLINE_INTEGRAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 February 2011
//
//  Author:
//
//    John Burkardt
//
{
  string comment;
  double *dn;
  double *fn;
  int i;
  int nn = 11;
  double pi = 3.141592653589793;
  double q_computed;
  double q_exact;
  int seed;
  int test;
  double *xn;

  seed = 123456789;

  std::cout << "\n";
  std::cout << "TEST10:\n";
  std::cout << "  HERMITE_CUBIC_SPLINE_INTEGRAL integrates a Hermite\n";
  std::cout << "  cubic spline over the definition interval [X1,XNN].\n";
  std::cout << "\n";
  std::cout << "  If the subintervals are equally spaced, the derivative\n";
  std::cout << "  information has no effect on the result, except for\n";
  std::cout << "  the first and last values, DN(1) and DN(NN).\n";
  std::cout << "\n";
  std::cout << "                            Exact       Computed\n";
  std::cout << "     X1          XNN        Integral    Integral  Comment\n";
  std::cout << "\n";

  fn = new double[nn];
  dn = new double[nn];

  for ( test = 1; test <= 5; test++ )
  {
//
//  Equal spacing.
//
    if ( test == 1 )
    {
      xn = r8vec_even_new ( nn, 0.0, pi );
      for ( i = 0; i < nn; i++ )
      {
        fn[i] = std::sin ( xn[i] );
        dn[i] = std::cos ( xn[i] );
      }
      q_exact = - std::cos ( xn[nn-1] ) + std::cos ( xn[0] );
      comment = "Equal spacing, correct DN";
    }
//
//  Equal spacing, reset DN(2:NN-1) to random numbers.
//
    else if ( test == 2 )
    {
      xn = r8vec_even_new ( nn, 0.0, pi );
      for ( i = 0; i < nn; i++ )
      {
        fn[i] = std::sin ( xn[i] );
        if ( i == 0 || i == nn - 1 )
        {
          dn[i] = std::cos ( xn[i] );
        }
        else
        {
          dn[i] = 1000.0 * r8_uniform_01 ( &seed );
        }
      }
      q_exact = - std::cos ( xn[nn-1] ) + std::cos ( xn[0] );
      comment = "Equal spacing, DN(2:N-1) random";
    }
//
//  Equal spacing, now reset all of DN to random numbers.
//
    else if ( test == 3 )
    {
      xn = r8vec_even_new ( nn, 0.0, pi );
      for ( i = 0; i < nn; i++ )
      {
        fn[i] = std::sin ( xn[i] );
        dn[i] = 1000.0 * r8_uniform_01 ( &seed );
      }
      q_exact = - std::cos ( xn[nn-1] ) + std::cos ( xn[0] );
      comment = "Equal spacing, DN(1:N) random";
    }
//
//  Variable spacing, correct data.
//
    else if ( test == 4 )
    {
      xn = r8vec_even_new ( nn, 0.0, pi * pi );
      for ( i = 0; i < nn; i++ )
      {
        xn[i] = std::sqrt ( xn[i] );
        fn[i] = std::sin ( xn[i] );
        dn[i] = std::cos ( xn[i] );
      }
      q_exact = - std::cos ( xn[nn-1] ) + std::cos ( xn[0] );
      comment = "Variable spacing, correct DN";
    }
//
//  Variable spacing, change one entry in DN.
//
    else if ( test == 5 )
    {
      xn = r8vec_even_new ( nn, 0.0, pi * pi );
      for ( i = 0; i < nn; i++ )
      {
        xn[i] = std::sqrt ( xn[i] );
        fn[i] = std::sin ( xn[i] );
        dn[i] = std::cos ( xn[i] );
      }
      dn[ ( nn - 1 ) / 2 ] = 1000.0 * r8_uniform_01 ( &seed );
      q_exact = - std::cos ( xn[nn-1] ) + std::cos ( xn[0] );
      comment = "Variable spacing, a single internal DN randomized.";
    }

    q_computed = hermite_cubic_spline_integral ( nn, xn, fn, dn );

    std::cout << "  " << setw(10) << xn[0]
         << "  " << setw(10) << xn[nn-1]
         << "  " << setw(14) << q_exact
         << "  " << setw(14) << q_computed
         << "  " << comment << "\n";

    delete [] xn;
  }

  delete [] dn;
  delete [] fn;

  return;
}
//****************************************************************************80

void test11 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST11 tests HERMITE_CUBIC_LAGRANGE_VALUE.
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
{
  double *d;
  double *f;
  int j;
  int n = 11;
  double *s;
  double *t;
  double *x;
  double x1;
  double x2;

  std::cout << "\n";
  std::cout << "TEST11:\n";
  std::cout << "  HERMITE_CUBIC_LAGRANGE_VALUE evaluates the four\n";
  std::cout << "  Lagrange basis functions associated with F1, D1,\n";
  std::cout << "  F2 and D2 such that\n";
  std::cout << "\n";
  std::cout << "  P(X) = F1 * LF1(X) + D1 * LD1(X)\n";
  std::cout << "       + F2 * LF2(X) + D2 * LD2(X).\n";
  std::cout << "\n";
  std::cout << "  The first, second and third derivatives of these four\n";
  std::cout << "  Lagrange basis functions are also computed.\n";

  x1 = 1.0;
  x2 = 2.0;
  x = r8vec_even_new ( n, 0.0, 2.5 );

  f = new double[4*n];
  d = new double[4*n];
  s = new double[4*n];
  t = new double[4*n];

  hermite_cubic_lagrange_value ( x1, x2, n, x, f, d, s, t );

  std::cout << "\n";
  std::cout << "  The Lagrange basis functions:\n";
  std::cout << "\n";
  std::cout << "     I        X           LF1         LD1         LF2         LD2\n";
  std::cout << "\n";
  for ( j = 0; j < n; j++ )
  {
    std::cout << "  " << setw(4) << j
         << "  " << setw(10) << x[j]
         << "  " << setw(10) << f[0+j*4]
         << "  " << setw(10) << f[1+j*4]
         << "  " << setw(10) << f[2+j*4]
         << "  " << setw(10) << f[3+j*4] << "\n";
  }

  std::cout << "\n";
  std::cout << "  The derivative of the Lagrange basis functions:\n";
  std::cout << "\n";
  std::cout << "     I        X           LF1         LD1         LF2         LD2\n";
  std::cout << "\n";
  for ( j = 0; j < n; j++ )
  {
    std::cout << "  " << setw(4) << j
         << "  " << setw(10) << x[j]
         << "  " << setw(10) << d[0+j*4]
         << "  " << setw(10) << d[1+j*4]
         << "  " << setw(10) << d[2+j*4]
         << "  " << setw(10) << d[3+j*4] << "\n";
  }

  delete [] d;
  delete [] f;
  delete [] s;
  delete [] t;
  delete [] x;

  return;
}
//****************************************************************************80

void test12 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST12 tests HERMITE_CUBIC_LAGRANGE_INTEGRAL.
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
{
  int i;
  double *q;
  double x1;
  double x2;

  std::cout << "\n";
  std::cout << "TEST12:\n";
  std::cout << "  HERMITE_CUBIC_LAGRANGE_INTEGRAL returns the integrals\n";
  std::cout << "  of the four Lagrange basis functions associated\n";
  std::cout << "  with F1, D1, F2 and D2 such that\n";
  std::cout << "\n";
  std::cout << "  P(X) = F1 * LF1(X) + D1 * LD1(X)\n";
  std::cout << "       + F2 * LF2(X) + D2 * LD2(X).\n";
  std::cout << "\n";
  std::cout << "  The Lagrange basis function integrals:\n";
  std::cout << "\n";
  std::cout << "        X1          X2          LF1         LD1         LF2         LD2\n";
  std::cout << "\n";

  x2 = 1.0;

  for ( i = -6; i <= 2; i++ )
  {
    x1 = ( double ) ( i );
    q = hermite_cubic_lagrange_integral ( x1, x2 );
    std::cout << "  " << setw(10) << x1
         << "  " << setw(10) << x2
         << "  " << setw(10) << q[0]
         << "  " << setw(10) << q[1]
         << "  " << setw(10) << q[2]
         << "  " << setw(10) << q[3] << "\n";
    delete [] q;
  }
  return;
}
//****************************************************************************80

void test13 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST13 tests HERMITE_CUBIC_LAGRANGE_INTEGRATE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 February 2011
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  double d1;
  double d2;
  double f1;
  double f2;
  int j;
  double p[4];
  double *q;
  double x1;
  double x2;

  std::cout << "\n";
  std::cout << "TEST13:\n";
  std::cout << "  HERMITE_CUBIC_LAGRANGE_INTEGRATE integrates a Hermite cubic\n";
  std::cout << "  Lagrange polynomial from A to B.\n";
  std::cout << "\n";
  std::cout << "  Compute each result TWICE:\n";
  std::cout << "  First row computed using HERMITE_CUBIC_INTEGRATE.\n";
  std::cout << "  Second row computed using HERMITE_CUBIC_LAGRANGE_INTEGRATE.\n";

  x1 = 0.0;
  x2 = 10.0;

  std::cout << "\n";
  std::cout << "        A           B           LF1         LD1         LF2         LD2\n";
  std::cout << "\n";

  a = x1 - 1.0;

  for ( j = -3; j <= 12; j++ )
  {
    b = ( ( double ) ( 10 - j ) * x1
        + ( double ) (      j ) * x2 )
        / ( double ) ( 10     );

    f1 = 1.0;
    d1 = 0.0;
    f2 = 0.0;
    d2 = 0.0;
    p[0] = hermite_cubic_integrate ( x1, f1, d1, x2, f2, d2, a, b );

    f1 = 0.0;
    d1 = 1.0;
    f2 = 0.0;
    d2 = 0.0;
    p[1] = hermite_cubic_integrate ( x1, f1, d1, x2, f2, d2, a, b );

    f1 = 0.0;
    d1 = 0.0;
    f2 = 1.0;
    d2 = 0.0;
    p[2] = hermite_cubic_integrate ( x1, f1, d1, x2, f2, d2, a, b );

    f1 = 0.0;
    d1 = 0.0;
    f2 = 0.0;
    d2 = 1.0;
    p[3] = hermite_cubic_integrate ( x1, f1, d1, x2, f2, d2, a, b );

    q = hermite_cubic_lagrange_integrate ( x1, x2, a, b );

    std::cout << "  " << setw(10) << a
         << "  " << setw(10) << b
         << "  " << setw(10) << p[0]
         << "  " << setw(10) << p[1]
         << "  " << setw(10) << p[2]
         << "  " << setw(10) << p[3] << "\n";
    std::cout << "  " << "          "
         << "  " << "          "
         << "  " << setw(10) << q[0]
         << "  " << setw(10) << q[1]
         << "  " << setw(10) << q[2]
         << "  " << setw(10) << q[3] << "\n";

    delete [] q;
  }

  return;
}
//****************************************************************************80

void test14 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST14 tests HERMITE_CUBIC_SPLINE_QUAD_RULE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 February 2011
//
//  Author:
//
//    John Burkardt
//
{
# define N 11

  double dn[N];
  double fn[N];
  int i;
  int j;
  int k;
  int l;
  int n = N;
  double q;
  double *r;
  int seed;
  double *w;
  double x[N];

  cout << "\n";
  cout << "TEST14:\n";
  cout << "  HERMITE_CUBIC_SPLINE_QUAD_RULE returns a quadrature rule\n";
  cout << "  for Hermite cubic splines.\n";

  seed = 123456789;

  for ( k = 1; k <= 2; k++ )
  {
    cout << "\n";
    if ( k == 1 )
    {
      cout << "  Case 1: Random spacing\n";
      r = r8vec_uniform_01_new ( n, &seed );

      x[0] = r[0];
      for ( i = 1; i < n; i++ )
      {
        x[i] = x[i-1] + r[i];
      }
      delete [] r;
    }
    else if ( k == 2 )
    {
      cout << "  Case 2: Uniform spacing\n";
      cout << "  F(2:N-1) have equal weight.\n";
      cout << "  D(2:N-1) have zero weight.\n";
      for ( i = 0; i < n; i++ )
      {
        x[i] = ( double ) ( 10 + i ) / 20.0;
      }
    }

    w = hermite_cubic_spline_quad_rule ( n, x );

    cout << "\n";
    cout << "   I   J        X         W                Q\n";
    cout << "\n";

    for ( i = 0; i <= 1; i++ )
    {
      for ( j = 0; j < n; j++ )
      {
        for ( l = 0; l < n; l++ )
        {
          fn[l] = 0.0;
          dn[l] = 0.0;
        }
        if ( i == 0 )
        {
          fn[j] = 1.0;
        }
        else
        {
          dn[j] = 1.0;
        }

        q = hermite_cubic_spline_integral ( n, x, fn, dn );

        cout << "  " << setw(2) << i
             << "  " << setw(2) << j
             << "  " << setw(10) << x[j]
             << "  " << setw(14) << q
             << "  " << setw(14) << w[i+j*2] << "\n";
      }
    }
    delete [] w;
  }
  return;
# undef N
}
//****************************************************************************80

void test15 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST15 tests HERMITE_CUBIC_SPLINE_QUAD_RULE.
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
{
# define N 11

  double dn;
  double fn;
  int j;
  int n = N;
  double q;
  double q_exact;
  double *r;
  double s;
  int seed;
  double t;
  double *w;
  double x[N];

  cout << "\n";
  cout << "TEST15:\n";
  cout << "  HERMITE_CUBIC_SPLINE_QUAD_RULE returns a quadrature rule\n";
  cout << "  for Hermite cubic splines.\n";

  seed = 123456789;

  r = r8vec_uniform_01_new ( n, &seed );

  x[0] = r[0];
  for ( j = 1; j < n; j++ )
  {
    x[j] = x[j-1] + r[j];
  }
  delete [] r;

  cout << "\n";
  cout << "  Random spacing\n";
  cout << "  Number of points N = " << n << "\n ";
  cout << "  Interval = [" << x[0] << ", " << x[n-1] << "]\n";

  w = hermite_cubic_spline_quad_rule ( n, x );

  q = 0.0;

  for ( j = 0; j < n; j++ )
  {
    cubic_value ( x[j], &fn, &dn, &s, &t );
    q = q + w[0+j*2] * fn + w[1+j*2] * dn;
  }

  q_exact = cubic_integrate ( x[0], x[n-1] );

  cout << "\n";
  cout << "  Q         = " << q << "\n";
  cout << "  Q (exact) = " << q_exact << "\n";

  delete [] w;

  return;
# undef N
}
//****************************************************************************80

double cubic_antiderivative ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    CUBIC_ANTIDERIVATIVE evaluates the antiderivative function of a cubic.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 January 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the argument.
//
//    Output, double CUBIC_ANTIDERIVATIVE, the value.
//
{
  double value;

  value = x * x * ( 5.0 + x * ( - 7.0 / 3.0 + x * 1.0 / 4.0 ) );

  return value;
}
//****************************************************************************80

double cubic_integrate ( double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    CUBIC_INTEGRATE integrates the cubic from A to B.
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
//  Parameters:
//
//    Input, double A, B, the integration interval.
//
//    Output, double Q, the integral from A to B.
//
{
  double q;

  q = cubic_antiderivative ( b ) - cubic_antiderivative ( a );

  return q;
}
//****************************************************************************80

void cubic_value ( double x, double *f, double *d, double *s, double *t )

//****************************************************************************80
//
//  Purpose:
//
//    CUBIC_VALUE evaluates a cubic function.
//
//  Discussion:
//
//    f(x) =   x^3 -  7 x^2 + 10 x
//    d(x) = 3 x^2 - 14 x   + 10
//    s(x) = 6 x   - 14
//    t(x) = 6
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
//  Parameters:
//
//    Input, double X, the argument.
//
//    Output, double F, D, S, T, the value and first three
//    derivatives of the cubic function.
//
{
  *f = 0.0 + x * ( 10.0 + x * ( -  7.0 + x * 1.0 ) );
  *d =             10.0 + x * ( - 14.0 + x * 3.0 );
  *s =                          - 14.0 + x * 6.0;
  *t =                                       6.0;

  return;
}