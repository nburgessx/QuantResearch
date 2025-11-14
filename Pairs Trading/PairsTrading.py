# -*- coding: utf-8 -*-
"""
Created on Fri Nov 14 00:00:53 2025

@author: Nicholas Burgess
"""

import numpy as np
import yfinance as yf
from statsmodels.tsa.stattools import adfuller, coint
from statsmodels.api import OLS, add_constant

# ADF test function: returns True if stationary, False otherwise
def adf_test(series, name="Series", significance=0.05):
    res = adfuller(series, autolag='AIC')
    p_value = res[1]
    is_stationary = p_value < significance
    print(f"\nADF Test for {name}: Statistic={res[0]:.4f}, p-value={p_value:.4f}")
    if is_stationary:
        print(f"Result: Reject H0 → {name} is stationary")
    else:
        print(f"Result: Fail to reject H0 → {name} is non-stationary")
    return is_stationary

# Half-life function: only meaningful for stationary series
def half_life(series):
    x = series.values
    dx = x[1:] - x[:-1]
    x_lag = x[:-1]
    X = add_constant(x_lag)
    est = OLS(dx, X).fit()
    phi_hat = 1 + est.params[1]
    if abs(phi_hat) < 1:
        hl = -np.log(2)/np.log(abs(phi_hat))
        print(f"Estimated phi: {phi_hat:.4f}, Half-life: {hl:.2f} periods")
        return hl
    else:
        print(f"Estimated phi: {phi_hat:.4f}, half-life undefined (non-stationary)")
        return np.nan

# Engle-Granger co-integration test with meaningful output
def engle_granger_test(series1, series2, name1="Series1", name2="Series2", significance=0.05):
    coint_res = coint(series1, series2)
    t_stat, p_value = coint_res[0], coint_res[1]
    print(f"\nEngle-Granger Test for {name1} & {name2}: t-stat={t_stat:.4f}, p-value={p_value:.4f}")
    if p_value < significance:
        print(f"Result: Reject H0 → {name1} and {name2} are co-integrated (stationary spread). Candidate for pairs trading.")
    else:
        print(f"Result: Fail to reject H0 → {name1} and {name2} are not co-integrated. Not suitable for pairs trading.")
    return p_value < significance

if __name__ == "__main__":
    
    # Fetch example data
    data_xom = yf.download("XOM", start="2024-11-12", end="2025-11-12", auto_adjust=True)["Close"].dropna()
    data_cvx = yf.download("CVX", start="2024-11-12", end="2025-11-12", auto_adjust=True)["Close"].dropna()

    # Test if XOM is stationary and compute half-life only if stationary
    if adf_test(data_xom, "XOM"):
        half_life(data_xom)
    
    # Test if CVX is stationary and compute half-life only if stationary
    if adf_test(data_cvx, "CVX"):
        half_life(data_cvx)
        
    # Engle-Granger co-integration test
    engle_granger_test(data_xom, data_cvx, "XOM", "CVX")
