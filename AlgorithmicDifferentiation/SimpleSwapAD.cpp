// We price a simple swap and demonstrate how to calculate the risk using algorithmic differentation (AD)
// Here we demo how to compute swap DV01 using tangent and adjoint mode as part of the pricing process.
#include <cmath>
#include <iostream>
#include <iomanip>
using namespace std;

void price_swap()
{
    // Swap Data: 1 year swap fixed vs float
    double phi = 1.0;   // ReceiveFixed=1.0 or PayFixed=-1.0
    double n = 1000000; // Notional
    double r = 0.02;    // Fixed Rate 2.0%
    double tau = 1.0;   // Coupon Year Fraction
    double t = 1.0;     // Swap Maturity
    
    // Market Data
    double f = 0.01;    // Forward Rate 1.0%
    double z = 0.02;    // Zero Rate 2.0% used by Discount Factor = exp(-z.t)
    
    // Swap Price
    double df           = exp(-z*t); // 1.
    double pv_fixed     = phi*n*r*tau*df; // 2.
    double pv_float     = -phi*n*f*tau*df; // 3.
    double pv_swap      = pv_fixed+pv_float; // 4.
    
    // Analytical Risk
    // pv01 = forward rate risk - change in swap pv for 1 bps change in forwards
    // df01 = discount factor risk - change in swap pv for change in discount factors (from 1bp in zero rate change)
    // dv01 = forward + df risk
    double pv01_swap    = -phi*n*tau*df*0.0001;             // dFloatPV/dForwardRate x 1 bps change in forward rate
    double df01_swap    = (phi*n*r*tau - phi*n*f*tau)       // dSwapPrice/dDiscountFactor
                        * (exp(-(z+0.0001)*t)-exp(-z*t));   // x change in df for 1 bps change in zero rate
    double dv01_swap    = pv01_swap + df01_swap;
    
    // Risk Shift Sizes
    double shift_size_f = 0.0001; // 1 bps shift 
    double shift_size_z = 0.0001; // 1 bps shift
    double original_df  = exp(-z*t);
    double shifted_df   = exp(-(z+shift_size_z)*t);
    double shift_size_df = shifted_df-original_df; // Change in df for a 1 bps change in zero rate

    // Tangent Mode: Forward Sweep
    // Tangent risks are denoted "dot". This mode works forwards and shifts the inputs
    // Consequently we only get one risk output at a time, bit like numerical bumping
    //
    // Example:     Function    y = 2x^2 
    //              Tangent     y_dot = 4x.x_dot 
    //              Given       x_dot = 1           // Risk Shift_size
    //              Result      dy/dx = 4x
    //
    double z_dot        = 1 * shift_size_z; // init shift size - shift applied to inputs in tangent mode
    double f_dot        = 1 * shift_size_f; // init shift size - shift applied to inputs in tangent mode
    double df_dot       = -t*exp(-z*t)*z_dot; // 1.
    double pv_fixed_dot = phi*n*r*tau*df_dot; // 2.
    double pv_float_dot = -phi*n*tau*df*f_dot -phi*n*f*tau*df_dot; // 3.
    double pv_swap_dot  = pv_fixed_dot+pv_float_dot; // 4.
    
    // Adjoint Mode: Backwards Sweep
    // Adjoint risks are denoted "bar". This mode works in reverse and shifts the outputs
    // Consequently we get a full break-down of all risks in one go
    //
    // Example:     Function    y = 2x^2
    //              Tangent     x_bar = 4x.y_bar
    //              Given       y_bar = 1           // Risk Shift_size
    //              Result      dy/dx = 4x
    //
    double pv_swap_bar  = 1.0; // init to one as shift-sizes are applied to outputs in adjoint mode
    double pv_fixed_bar = pv_swap_bar; // 4.
    double pv_float_bar = pv_swap_bar; // 4.
    double f_bar        = -phi*n*tau*df*pv_float_bar    *shift_size_f;  // 3.   Note: f_bar output  - apply required shift size
    double df_bar       = -phi*n*f*tau*pv_float_bar     *shift_size_df; // 3.   Note: df_bar output - apply required shift size
    df_bar              += phi*n*r*tau*pv_fixed_bar     *shift_size_df; // 2.   Note: df_bar output - apply required shift size
    double z_bar        = -t*exp(-z*t)*df_bar; // 1.
    
    // Results
    cout << "Swap Price" << endl;
    cout << "pv_swap:       " << std::fixed << std::setprecision(2) << pv_swap << endl << endl;
    
    cout << "Analytical Swap Risk" << endl;
    cout << "forward risk:  " << std::fixed << std::setprecision(2) << pv01_swap << " (pv01)" << endl;
    cout << "discount risk: " << std::fixed << std::setprecision(2) << df01_swap << endl;
    cout << "dv01:          " << std::fixed << std::setprecision(2) << dv01_swap << endl << endl;
    
    cout << "Tangent Mode Risk" << endl;
    cout << "pv_swap_dot:   " << std::fixed << std::setprecision(2) << pv_swap_dot << endl;
    cout << "dv01:          " << std::fixed << std::setprecision(2) << pv_swap_dot << endl << endl;
    
    cout << "Adjoint Mode Risk" << endl;
    cout << "f_bar:         " << std::fixed << std::setprecision(2) << f_bar << " (forward risk)" << endl;
    cout << "df_bar:        " << std::fixed << std::setprecision(2) << df_bar << " (discount risk)" << endl;
    cout << "dv01:          " << std::fixed << std::setprecision(2) << f_bar+df_bar << endl;
}

int main()
{
    // Price swap has algorithmic differentiation (AD) risk embedded
    // Therefore when pricing we get the exact risk in real-time for close to free
    price_swap();
    return 0;
}

