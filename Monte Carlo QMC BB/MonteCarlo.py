# -*- coding: utf-8 -*-
"""
Created on Sun Sep 17 15:36:36 2023
Monte Carlo Pricing
@author: Nicholas Burgess
"""

import math
import numpy as np

# Black-Scholes Geometric Brownian Motion
# s(t) = s(0) * exp( ( r - 0.5 * sigma**2 ) * dt) + sigma*sqrt(dt)*x )

# Stock Price Process
s0 = 100        
k = 105
r = 0.05
t = 1.0
sigma = 0.2

# Monte Carlo Settings
nSimulations = 100000
nTimeSteps = 100
dt = t / nTimeSteps

# Discount Factor
df = math.exp(-r*t)
print(f"Discount Factor: {df}")

# Random Numbers (Mersenne Twister)
np.random.seed(100)
rn = np.random.standard_normal((nTimeSteps+1, nSimulations))
print(f"Expected Value(Random Variates): {rn.mean():.6f}")
print(f"Variance(Random Variates): {rn.var():.6f}")

# Stock Price Process
s = np.zeros_like(rn)
s[0] = s0
for i in range(1, nTimeSteps+1):
    s[i] = s[i-1] * np.exp((r-0.5*sigma**2)*dt + sigma*math.sqrt(dt)*rn[i])

# Plot Monte Carlo Simulations
import matplotlib.pyplot as plt
plt.style.use("seaborn-v0_8")
plt.rcParams["font.family"] = "serif"
plt.rcParams["savefig.dpi"] = 300

plt.figure(figsize=(10,6))
plt.plot(s[:,:20]);

# Histogram to Distribution of Stock Price at Maturity, t

nStdDevs = 2 # Set Plot Boundary for +/- Number of Standard Deviations

s_t = s[-1]
plt.figure(figsize=(10,6))
plt.hist(s_t, bins=35, color="b", label="frequency")
plt.axvline(s_t.mean(), color="r", label="mean")
plt.axvline(s_t.mean() + nStdDevs*s_t.std(), color="g", label=f"+{nStdDevs} Standard Deviation(s)")
plt.axvline(s_t.mean() - nStdDevs*s_t.std(), color="g", label=f"-{nStdDevs} Standard Deviation(s)")
plt.legend(loc=0);


