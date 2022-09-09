#!/bin/bash

##SBATCH --job-name       PhysicalFilter_pepStream
##SBATCH --time           14:00:00
##SBATCH --ntasks         1
##SBATCH --cpus-per-task  4
##SBATCH --mem            10G
##SBATCH --error          PhysicalFilter_error.err
##SBATCH --output         PhysicalFilterRun.log
##SBATCH --partition     orion
#JNUM= ls *.fasta | wc 1
##SBATCH --array 0-"${JNUM}"%25

mkdir Filters
cp ./structures/PDB/contigous_seq_cMoRFs_*.fasta ./Filters
cd ./Filters
#cd ./Fragments
# runs first round of docking: coarse grain full strucutre
#cp ../run_iTasser_zDock.sh .
fileNames=()
cp ../run_iTasser_zDock.sh .
JNUM=0
for file in *.fasta;
do 
	echo "Current file is $file ..."
	tempName=${file:14:-5}
	fileNames+=($tempName)
	mkdir $tempName
	#mv $file ./$tempName
	cd $tempName
	mkdir iTasser
	mkdir zDock
	cp ../$file .
	cp -r ../../zDockModels/ .
	cp -r ../../zDockScripts/ .
#	cp ../run_iTasser_zDock.sh .
	chmod 775 *
#	sbatch run_iTasser_zDock.sh
	cd ..

	JNUM=$[ $JNUM + 1 ]
#	if [ $JNUM -gt 100 ]; then
#	break
#	fi
done

echo $JNUM
cp ../wrap_physFilter.sh .
###sbatch 
#./run_iTasser_zDock.sh "${fileNames[$SLURM_ARRAY_TASK_ID]}"

#sbatch wrap_physFilter.sh "$JNUM" "${fileNames[@]}"
sbatch --array "0-"${JNUM}"%25" wrap_physFilter.sh "${fileNames[@]}"
echo $JNUM
