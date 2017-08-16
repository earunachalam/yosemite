#!/bin/bash --login

# extract_spaced_all.sh
# extract_spaced for all trajectories

# name of trajec
TRAJNAME=traj.xtc

# name to assign conf.pdb for original coordinates (at start of seed trajectory) in each directory
LOCALCONF=init_conf.pdb

# parameters for each extract_spaced, in order of directory name
# N_STEPS T_START T_END N_REP
params=( (`seq 1 10`) (`seq 2 11`) )
echo $params
exit

# get absolute path to conf.pdb file
ORIGCONF=$(readlink -f conf.pdb)

DIRS=$(find .. -type d)
for DIR in ${DIRS}
do
	# remove ./from filename
	CDIR=${DIR///}
	CDIR=${CDIR//.}
	if [[ ${CDIR} =~ ^[0-9]+$ ]]
	then
		echo Starting ${CDIR}...
		cd ${DIR}

		#cp ${ORIGCONF} ${LOCALCONF}
		
		echo extract_spaced traj.xtc ${LOCALCONF} 10 

		exit
	fi
done
