KSpline.h 
=========

Documentation
=============
Source: https://github.com/ttk592/spline/
Documentation: https://kluge.in-chemnitz.de/opensource/spline/

Examples
========
Splines can be built as follows:

// Our function values for y = f(x)
const std::vector<double> xValues{1.0, 2.0, 3.0};
const std::vector<double> yValues{10.0, 20.0, 30.0};

// We have a flag to set if we want the spline to be monotonic
const bool make_montonic = false;

// Spline Type - Can be cubic spline or cubic hermite spline as per below comment
kluge::KSpline::spline_type spline_type = kluge::KSpline::cspline;
//kluge::KSpline::spline_type spline_type = kluge::KSpline::cspline_hermite;

// Set-Up Interpolator
kluge::KSpline myCubicSpline(xValues, yValues, spline_type, make_montonic);

// Example Getting an Interpolated Value at point 1.5
const double interp = myCubicSpline(1.5);

// Example Getting a First Order Derivative at point 1.5
const double deriv = myCubicSpline.deriv(1, 1.5);

// Example Integrating from 1.0 to 2.0
const double integral = myCubicSpline.integral(1.0, 2.0);