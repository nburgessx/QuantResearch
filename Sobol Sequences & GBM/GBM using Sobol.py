# -*- coding: utf-8 -*-
"""
Created on Sat Sep  9 10:50:54 2023

@author: nburgessx
"""

import sobol_seq
import numpy as np
from scipy.stats import norm

from matplotlib import pyplot


def generate_sobol_standard_normal_variates(dimension, nPoints, skip=1):
    ''' 
    skip Parameter
    Indicates the number of initial points to skip, we usually skip the first value of zero
    '''
    sobols = sobol_seq.i4_sobol_generate(dimension, nPoints, skip)
    normals = norm.ppf(sobols)
    return normals


def generate_gbm_paths(s0, r, vol, t, n_paths=8192, dt_days=3, sbol_skip=0):
    
    # Set the Time Interval Periodicity, here daily
    dt = dt_days / 365.25

    # Build Time Points
    time_points = []
    t_ = 0
    while t_ + dt < t:
        t_ += dt
        time_points.append(t_)

    time_points.append(t)
    time_points = np.array(time_points)

    n_time_points = len(time_points)

    # Generate Random Numbers before Evaluating Paths
    rand = generate_sobol_standard_normal_variates(n_time_points, n_paths, skip=sbol_skip)

    # Initialize Paths
    paths = np.zeros((n_paths, n_time_points))
    paths[:,0] = s0
    
    # Update GBM Paths
    for t_i in range(1, n_time_points):
        dt_ = time_points[t_i] - time_points[t_i-1]
        paths[:, t_i] = paths[:, t_i-1] * np.exp((r - 0.5*vol**2) * dt_ + np.sqrt(dt_) * vol * rand[:, t_i])

    # Return time points and paths
    return time_points, paths

# Monte Carlo Paths for Sobol Sequences must be a power of 2 e.g. 2^10 = 1,024
n_paths = 1024
time_points, paths = generate_gbm_paths(100, 0, 0.2, 0.3, n_paths=n_paths)

# Plot Paths as Subplots
fig = pyplot.figure()
thisPath = fig.add_subplot(1,1,1)

# Plot Paths and Apply Chart Colour and Transparency Settings
for i in range(n_paths):
    thisPath.plot(time_points, paths[i, :], c=(1,0,0,0.05))

# Show Chart
pyplot.show()