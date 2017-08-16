#!/bin/bash

# find_log_lammps.sh
# find each of the log.lammps files (from each replicate of each starting configuration) and pass their paths to another script which extracts the potential energy as a function of time

find ../../pass2 -name allreps.dat -exec rm {} +
find ../../pass2 -name log.lammps -exec ./fit_freeze.sh {} \;
