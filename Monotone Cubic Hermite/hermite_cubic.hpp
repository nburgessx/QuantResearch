double hermite_cubic_integral ( double x1, double f1, double d1, double x2, double f2, double d2 );

double hermite_cubic_integrate ( double x1, double f1, double d1, double x2, double f2, double d2, double a, double b );

double *hermite_cubic_lagrange_integral ( double x1, double x2 );

double *hermite_cubic_lagrange_integrate ( double x1, double x2, double a, double b );

void hermite_cubic_lagrange_value ( double x1, double x2, int n, double x[], double f[], double d[], double s[], double t[] );

double hermite_cubic_spline_integral ( int nn, double xn[], double fn[], double dn[] );

double *hermite_cubic_spline_integrate ( int nn, double xn[], double fn[], double dn[], int n, double a[], double b[] );

double *hermite_cubic_spline_quad_rule ( int nn, double xn[] );

void hermite_cubic_spline_value ( int nn, double xn[], double fn[], double dn[], int n, double x[], double f[], double d[], double s[], double t[] );

void hermite_cubic_to_power_cubic ( double x1, double f1, double d1, double x2, double f2, double d2, double *c0, double *c1, double *c2, double *c3 );

void hermite_cubic_value ( double x1, double f1, double d1, double x2, double f2, double d2, int n, double x[], double f[], double d[], double s[], double t[] );

void power_cubic_to_hermite_cubic ( double c0, double c1, double c2, double c3, double x1, double x2, double *f1, double *d1, double *f2, double *d2 );

double r8_abs ( double x );

double r8_uniform_01 ( int *seed );

void r8vec_bracket3 ( int n, double t[], double tval, int *left );

double *r8vec_even_new ( int n, double alo, double ahi );

double *r8vec_uniform_01_new ( int n, int *seed ); 

void timestamp ( );