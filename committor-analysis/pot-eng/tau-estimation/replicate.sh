#!/bin/bash --login

# replicate
# rerun LAMMPS simulation N_REP times with identical starting configuration but different velocities

N_REP=10

for isim in `seq 1 ${N_REP}`
do
	ISIMDIR=dat/${isim}
	mkdir -p ${ISIMDIR}
	cp {mW-C.sw,init.data} ${ISIMDIR}/
	cp run.sh.template ${ISIMDIR}/run.sh
	cp in.npt.template ${ISIMDIR}/in.npt
	
	(
		cd ${ISIMDIR}

		# set seed for new velocities
		CMD=s/template_seed/$RANDOM/g
		perl -i -pe ${CMD} in.npt

		# set task name for submission
		CMD=s/template_taskname/lowT_c${isim}/g
		perl -i -pe ${CMD} run.sh

		# submit
		qsub run.sh
	)
done
