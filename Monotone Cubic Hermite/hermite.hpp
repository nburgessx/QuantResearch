void dif_deriv ( int nd, double xd[], double yd[], int *ndp, double xdp[], 
  double ydp[] );
void dif_shift_x ( int nd, double xd[], double yd[], double xv );
void dif_shift_zero ( int nd, double xd[], double yd[] );
void dif_to_r8poly ( int nd, double xd[], double yd[], double c[] );
double *dif_vals ( int nd, double xd[], double yd[], int nv, double xv[] );
double hermite_basis_0 ( int n, double x[], int i, double xv );
double hermite_basis_1 ( int n, double x[], int i, double xv );
void hermite_demo ( int n, double x[], double y[], double yp[] );
void hermite_interpolant ( int n, double x[], double y[], double yp[], 
  double xd[], double yd[], double xdp[], double ydp[] );
double *hermite_interpolant_rule ( int n, double a, double b, double x[] );
void hermite_interpolant_value ( int nd, double xd[], double yd[], double xdp[], 
  double ydp[], int nv, double xv[], double yv[], double yvp[] );
double r8_abs ( double x );
double r8_max ( double x, double y );
double r8poly_ant_val ( int n, double poly_cof[], double xval );
int r8poly_degree ( int na, double a[] );
void r8poly_print ( int n, double a[], std::string title );
double *r8vec_chebyshev_new ( int n, double a_first, double a_last );
double *r8vec_linspace_new ( int n, double a_first, double a_last );
void r8vec_print ( int n, double a[], std::string title );
double r8vec_product ( int n, double a[] );
void r8vec_uniform_01 ( int n, int *seed, double r[] );
void timestamp ( );