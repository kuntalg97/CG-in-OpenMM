# Kuntal Ghosh
# Code to extrapolate data

import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline

filename = 'u_bumper.dat'
n_lines = len(open(filename).readlines())
n_total = 2001

# given values
xi, yi = np.loadtxt(filename, unpack=True)

# positions to inter/extrapolate
x = np.linspace(0.0, 25.022, n_total)

# spline order: 1 linear, 2 quadratic, 3 cubic ... 
order = 1 # quadratic splines work best for extrapolating pair potentials and linear splines for forces

# do inter/extrapolation
s = InterpolatedUnivariateSpline(xi, yi, k=order)
y = s(x)

out_file = open ('u_extrap.dat', 'w')

for i in range (1,n_total+1):
    out_file.write ("%16.6f%16.6f\n"%(x[i-1],y[i-1]))
