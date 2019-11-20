#!/bin/bash

prefix=$1
N=$2
NUM_REPLICAS=$3
NUM_NODES=$((NUM_REPLICAS/16))
if ((NUM_NODES==0))
then
	NUM_NODES=1
fi

for ((i=1 ; i<=N ; i++))
do
	echo $i
	if [[ -n "${jobids}" ]]
	then
		echo "waiting for jobids=$jobids"
		jobids="--dependency=afterok:${jobids}"
	fi
	jobname="exo_${prefix}_${i}"
	jobfile="${prefix}job_${i}.sh"
	sed -e "s|PREFIX|${prefix}|g" -e "s|STEP|${i}|g" -e "s|NUM_REPLICAS|${NUM_REPLICAS}|g" job.sh > ${jobfile}
	echo submitting...
	jobids=$(sbatch --job-name ${jobname} ${jobids} -N ${NUM_NODES} -n ${NUM_REPLICAS} ${jobfile} | awk '{print $NF}')
	if [[ -n "${jobids}" ]]
	then
		echo "waiting for jobids=$jobids"
		jobids="--dependency=afterok:${jobids}"
	fi

	seedname="exo_seed_${prefix}_${i}"
	seedsfile="${prefix}seeds_${i}.sh"
	sed -e "s|PREFIX|${prefix}|g" -e "s|STEP|${i}|g" -e "s|NUM_REPLICAS|${NUM_REPLICAS}|g" get_seeds.sh > ${seedsfile}
	echo submitting...
	jobids=$(sbatch --job-name ${seedname} ${jobids} ${seedsfile} | awk '{print $NF}')
done
