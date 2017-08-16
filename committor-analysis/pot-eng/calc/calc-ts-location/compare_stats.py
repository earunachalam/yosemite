#!/usr/bin/python3

from matplotlib import rc
import matplotlib.pyplot as plt
import numpy as np
import os
from scipy.optimize import curve_fit
import sys

rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

median_tau_125 = 0
for init_loc in range(16,256,16):
    
    currstats = np.loadtxt("../../pass2/dat/" + str(init_loc) + "/allreps.dat")
    curr99 = currstats[:,4]
    curr01 = currstats[:,6]
    currtau = curr01 - curr99
    
    if init_loc == 125:
        median_tau_125 = np.median(currtau)
        print("Median tau for 125 is %.2f ns" % median_tau_125)

    print(np.sum(curr99 <= median_tau_125))
    print("%.2f of trajectories initiated from configuration at timestep %d end up in the ice basin." % (float(np.sum(curr99 <= median_tau_125))/len(curr99), init_loc))
    
    plt.scatter([init_loc for i in currtau], currtau)

plt.xlabel("Starting timestep (from original trajectory)")
plt.ylabel("Time to start of descent")
plt.ylim([-100, 100])
plt.show()
