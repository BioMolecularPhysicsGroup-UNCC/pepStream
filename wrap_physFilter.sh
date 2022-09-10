#!/bin/bash
#SBATCH --job-name       PhysicalFilter_pepStream
#SBATCH --time           14:00:00
#SBATCH --ntasks         1
#SBATCH --cpus-per-task  4
#SBATCH --mem            10G
#SBATCH --error          PhysicalFilter_error.err
#SBATCH --output         PhysicalFilterRun.log
##SBATCH --array="0-"${1}"%25"

#JNUM= ls *.fasta | wc 1
#JNUM=$1
#echo $JNUM
fileNames=( "$@" )


file="${fileNames[$SLURM_ARRAY_TASK_ID]}"
echo ${file:14:-5}
cd ${file:14:-5}
./run_iTasser_zDock.sh "${fileNames[$SLURM_ARRAY_TASK_ID]}"
cd ..
