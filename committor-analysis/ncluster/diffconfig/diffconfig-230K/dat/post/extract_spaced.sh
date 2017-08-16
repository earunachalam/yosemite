#!/bin/bash --login

# extract_spaced.sh
# chop up trajectory into separate init files, then initiate sets of simulations from each of these coordinates
# takes as input arguments:

# 1: the relative path to the trajectory file
# 2: the relative path to the atomic configuration (pdb) file
# 3: the number of points from which trajectories are initiated
# 4: first timestep from which we should initiate trajectories
# 5: last timestep
# 6: number of trajectories to initiate from each starting set of coordinates

TRAJECTORY=$1
CONFPDB=$2
DAT=$3
VMD_EXEC=/tmp/extract_spaced.vmd

N_STEPS=$4		# number of points from which trajectories are initiated
T_START=$5		# first timestep from which we should initiate trajectories
T_END=$6		# last timestep

N_REP=$7		# number of trajectories to initiate from each starting set of coordinates

# do not edit below this line

T_LENGTH=$(bc <<< "${T_END}-${T_START}")
T_INTERVAL=$(bc <<< "${T_LENGTH}/${N_STEPS}")

# create separate files to initiate trajectories


echo "Creating..."

# delete previously created files
rm -rf ${DAT}
mkdir ${DAT}

# write vmd conversion script
echo package require topotools > ${VMD_EXEC}
echo mol new {${CONFPDB}} type {pdb} first 0 last -1 step 1 waitfor all >> ${VMD_EXEC}
echo mol addfile {${TRAJECTORY}} type {xtc} first 0 last -1 step 1 waitfor all >> ${VMD_EXEC}

for tstep in `seq ${STARTIDX} ${N_STEPS}`
do
	T_CURRENT=$(bc <<< "${tstep}*${T_INTERVAL}+${T_START}")
	echo animate write lammpstrj init${T_CURRENT}.lammpstrj beg ${T_CURRENT} end ${T_CURRENT} skip 1 0 >> ${VMD_EXEC}
done

echo quit >> ${VMD_EXEC}

# run vmd script
vmd -dispdev text -e ${VMD_EXEC}

# delete the temporary vmd script
rm ${VMD_EXEC}

# move chunked dump files and other required files to appropriate data directory
for tstep in `seq ${STARTIDX} ${N_STEPS}`
do

	T_CURRENT=$(bc <<< "${tstep}*${T_INTERVAL}+${T_START}")
	DIR_CURRENT=${DAT}/${T_CURRENT}
	mkdir -p ${DIR_CURRENT}

	# put required files in run directory
	mv init${T_CURRENT}.lammpstrj ${DIR_CURRENT}/
	cp {mW-C.sw,init.data} ${DIR_CURRENT}/
	cp in.npt.template ${DIR_CURRENT}/in.npt
	cp run.sh.template ${DIR_CURRENT}/run.sh

	# set init file (depends on frame from which you're starting)
	CMD=s/template_iframe/${T_CURRENT}/g
	(cd ${DIR_CURRENT} && perl -i -pe ${CMD} in.npt)
	
	# set taskname (depends on frame from which you're starting)
	CMD=s/template_taskname/nctrj_f${T_CURRENT}/g
	(cd ${DIR_CURRENT} && perl -i -pe ${CMD} run.sh)
	
done
