#!/bin/bash

#$ -pe orte 8
#$ -cwd
#$ -q all.q
#$ -j N
#$ -S /bin/bash
#$ -N nctrj_f100000

CMD="s/template_seed/$RANDOM/g"
perl -i -pe ${CMD} in.npt

mpirun -np ${NSLOTS} lmp_mpi -in in.npt -screen none
