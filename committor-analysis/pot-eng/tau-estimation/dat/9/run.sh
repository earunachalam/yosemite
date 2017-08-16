#!/bin/bash

#$ -pe orte 8
#$ -cwd
#$ -q all.q
#$ -j N
#$ -S /bin/bash
#$ -N 290-220K-hold_c9

mpirun -np ${NSLOTS} lmp_mpi -in in.npt -screen none
