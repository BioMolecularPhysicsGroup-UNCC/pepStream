#!/bin/bash
#SBATCH --job-name       AnalyticalFilter_pepStream
#SBATCH --time           24:00:00
#SBATCH --ntasks         1
#SBATCH --cpus-per-task  4
#SBATCH --mem            4G
#SBATCH --error          AnalyticalFilter_error.err
#SBATCH --output         AnalyticalFilterRun.log
cd $SLURM_SUBMIT_DIR

module load R
module load parallel/20210222 i-tasser/5.1
#module load jdk/9.0.4
module load blast

#cd iTasserSpecific
#cp ../*fasta .

# store fasta file as a var, subsection for the cMoRF name
file=$(ls ../*.fasta)
systemName=${file:14:-6}
# ------------------------------
# it's not worth the wait to re-run iTasser i guess. Too parametric.
# Might be that i-Tasser/AF2 run is sufficent ...
# ------------------------------
#mv *fasta seq.fasta
# make variable input for homological cut off,idcut, want two runs at 0.3 and 0.8
#cutOff=0.8
#/apps/pkg/i-tasser/5.1/I-TASSERmod/runI-TASSER.pl \
#-runstyle gnuparallel -pkgdir /apps/pkg/i-tasser/5.1 \
#-libdir /users/jpatte75/ITLIB \
#-LBS false -EC true -GO false \
#-homoflag real -idcut $cutOff \
#-nmodel 10 -ntemp 40 \
#-seqname $systemName \
#-datadir "$SLURM_SUBMIT_DIR/iTasserSpecific" \
#-java_home /usr -light true \
#-outdir "$SLURM_SUBMIT_DIR/SpecificiTasser" > output_${systemName}_${SLURM_JOB_ID}.log


#wait
#cd ..
# change name of old models first!
#rename model first_model *pdb
#cp ./iTasserSpecific/model*pdb .
rm -r zDockSpecific
mkdir zDockSpecific
cd zDockSpecific
#cp ../*pdb .

#rm first_first*
# copy the specific docking models
cp ../zDockModels/specific_models/*_m.pdb .
cp ../zDockScripts/* .
cp ../zDock/*_m.pdb .

# not needed if models already exist as _m
# for loop through models
#for model in *model*.pdb
#do
#mark_sur $model ${model:0:-4}_m.pdb
#done

# fix into for file in ?/dirWithSpecifcs/*
#for VARIABLE in 1 2 3 4 5 6 7 8 9 10

targetNames=( ../zDockModels/specific_models/*_m.pdb)
for TARGET in ${targetNames[@]##*/}
do
#the inner loop should be fine
	for model in model*_m.pdb
	do
	~/bin/zdock -R $TARGET -L $model -N 200 -o ${model:0:-3}_${TARGET}.out
#	~/bin/zdock -R DELLA_Model${VARIABLE}_m_Loop.pdb -L $model -N 200 -o ${model:0:-3}_DELLA_Model${VARIABLE}_Loop.out
#	~/bin/zdock -R DELLA_Model${VARIABLE}_m_MoRF1.pdb -L $model -N 200 -o ${model:0:-3}_DELLA_Model${VARIABLE}_MoRF1.out
#	~/bin/zdock -R DELLA_Model${VARIABLE}_m_MoRF12.pdb -L $model -N 200 -o ${model:0:-3}_DELLA_Model${VARIABLE}_MoRF12.out
#	~/bin/zdock -R DELLA_Model${VARIABLE}_m_MoRF2.pdb -L $model -N 200 -o ${model:0:-3}_DELLA_Model${VARIABLE}_MoRF2.out
#	~/bin/zdock -R DELLA_Model${VARIABLE}_m_MoRF3.pdb -L $model -N 200 -o ${model:0:-3}_DELLA_Model${VARIABLE}_MoRF3.out
	done
done


echo "Work up"
#grep someting out of summary shell
chmod +x summarize_zdock_spec.sh
./summarize_zdock_spec.sh
Rscript bootstrapping_zdockScores_spec.R > specific_outputfile.txt
cp specific_outputfile.txt ../
cd ..
#rm -r *iTasser*

