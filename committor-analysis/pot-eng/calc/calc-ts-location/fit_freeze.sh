#!/bin/bash

# fit_freeze.sh
# extract potential energy as a function of time from log.lammps file
# use this 'cleaned' file to fit U9t) to a logistic function (using the Python program of the same name) and dump the fit parameters from each replicate for a given starting configuration into a single file

echo $1
PYDIR=$(pwd)
FILEDIR=$(dirname $1)
echo FILEDIR=${FILEDIR}
(
	cd ${FILEDIR}
	awk 'NF==8' log.lammps > /tmp/tmp.dat
	sed -i '1,4d;$d' /tmp/tmp.dat;
	awk -F' ' '{print $1 " " $2}' /tmp/tmp.dat > /tmp/cleaned.dat
	cp /tmp/cleaned.dat tmp.dat
	#python3 ${PYDIR}/fit_freeze.py >> ../allreps.dat
)
