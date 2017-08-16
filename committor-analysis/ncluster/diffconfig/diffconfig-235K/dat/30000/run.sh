#!/bin/bash

#$ -pe orte 64
#$ -cwd
#$ -q ibnet
#$ -j N
#$ -S /bin/bash
#$ -N nctrj_f30000

CMD="s/template_seed/$RANDOM/g"
perl -i -pe ${CMD} in.npt

mpirun -np ${NSLOTS} lmp_mpi -in in.npt -screen none