// Program to demonstrate tangent and adjoint differentiation for a simple function f(x)
#include <iostream>
using namespace std;

// Forward 'Tangent' Method
// Risk variables are denoted 'dot' with '_d' suffix
void tangent(const double & x)
{
    // f(x) = 2x^2 + 3x
    // df/dx = 4x + 3
    
    double x_d = 1.0;           // init
    
    double a = x*x;             // 1.
    double a_d = 2*x*x_d;
    
    double b = x;               // 2.
    double b_d = x_d;
    
    double f = 2*a + 3*b;       // 3.
    double f_d = 2*a_d + 3*b_d;
    
    // Tangent Results
    cout << "Tangent Results" << endl;
    cout << "f(" << x << ")=" << f << endl;
    cout << "df/dx = " << f_d << endl;
    cout << endl;
}

// Backward 'Adjoint' Method
// Risk variables are denoted 'bar' with '_b' suffix
void adjoint(const double & x)
{
    // f(x) = 2x^2 + 3x
    // df/dx = 4x + 3
    
    // Forward Sweep
    double a = x*x;             // 1.
    double b = x;               // 2.
    double f = 2*a + 3*b;       // 3.
    
    // Back Propogation
    double f_b = 1.0;           // init
    double a_b = 2*f_b;         // 3.
    double b_b = 3*f_b;         // 3.
    double x_b = b_b;           // 2.
    x_b += 2*x*a_b;             // 1.

    // Adjoint Results
    cout << "Adjoint Results" << endl;
    cout << "f(" << x << ")=" << f << endl;
    cout << "df/dx = x_bar = " << x_b << endl;
    cout << endl;
}

void showFunction()
{
    cout << "Function" << endl;
    cout << "f(x) = 2x^2 + 3x" << endl;
    cout << "df/dx = 4x+3" << endl;
    cout << endl;
}

int main()
{
    showFunction();
    tangent(2);
    adjoint(2);
}



