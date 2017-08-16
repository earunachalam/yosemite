#!/bin/bash

#$ -pe orte 64
#$ -q ibnet
#$ -cwd
#$ -j N
#$ -S /bin/bash
#$ -N long-run-290K

rm long-run-290K.*

CMD="s/template_seed/$RANDOM/g"
cp in.npt.template in.npt
perl -i -pe ${CMD} in.npt
mpirun -np ${NSLOTS} lmp_mpi -in in.npt -screen none
