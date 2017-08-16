#!/bin/bash

#$ -pe orte 8
#$ -cwd
#$ -q all.q
#$ -j N
#$ -S /bin/bash
#$ -N template_taskname

mpirun -np ${NSLOTS} lmp_mpi -in in.npt -screen none
