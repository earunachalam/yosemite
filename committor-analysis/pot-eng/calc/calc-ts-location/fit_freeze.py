#!/usr/bin/python3

# fit_freeze.py
# fit potential energy as a function of time to a logistic function. This is applied to a thermo output file from LAMMPS (for a trajectory during which ice forms from liquid water) which has been cleaned so that the first column is the time (number of 10fs timesteps elapsed) and the second is the potential energy of the system.


from matplotlib import rc
import matplotlib.pyplot as plt
import numpy as np
import os
from scipy.optimize import curve_fit
import sys


rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

def func(x, s, L, k, x0):
    return s + L/(1. + np.exp(-k*(x-x0)))

# data file with col1 = time, col2 = U
data = np.loadtxt('/tmp/cleaned.dat')
t = data[:,0]/1.e5  # now in ns
U = data[:,1]       # still in kcal/mol

try:
    popt, pcov = curve_fit(func, t, U, p0=[-23500, 1850, -8, 1.0])
    
    s = popt[0]
    L = popt[1]
    k = popt[2]
    t0 = popt[3]
    def x_at_f_of_L(f):
        return t0 - (1./k)*np.log(1./f - 1.)

    print(s, L, k, t0, x_at_f_of_L(0.99), func(x_at_f_of_L(0.99), *popt), x_at_f_of_L(0.01), func(x_at_f_of_L(0.01), *popt))

    # # plot
    # plt.plot(t,U)
    # plt.plot(t, func(t,*popt), 'k-', label='fit')

    # plt.scatter([x_at_f_of_L(0.99)], [func(x_at_f_of_L(0.99), *popt)], s=20, c='r', zorder=4)
    # plt.scatter([x_at_f_of_L(0.01)], [func(x_at_f_of_L(0.01), *popt)], s=20, c='r', zorder=5)

    # plt.xlabel(r'$t$/ns')
    # plt.ylabel(r'$U$/[kcal/mol]')

    # plt.show()
except:
    # do nothing
    ignore = 0
