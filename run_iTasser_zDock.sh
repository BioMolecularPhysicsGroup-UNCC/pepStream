#!/bin/bash
##SBATCH --job-name       PhysicalFilter_pepStream
##SBATCH --time           14:00:00
##SBATCH --ntasks         1
##SBATCH --cpus-per-task  8
##SBATCH --mem            32G
##SBATCH --error          PhysicalFilter_error.err
##SBATCH --output         PhysicalFilterRun.log
##SBATCH --partition     orion
cd $SLURM_SUBMIT_DIR
cd $1

echo $1


module load R
module load parallel/20210222 i-tasser/5.1
#module load jdk/9.0.4
module load blast

cd iTasser
cp ../*fasta .
mv *fasta seq.fasta
wait

# store fasta file as a var, subsection for the cMoRF name
file=$(ls ../*.fasta)
#systemName=${file:14:-5}
# make variable input for homological cut off,idcut, want two runs at 0.3 and 0.8
cutOff=0.3

/apps/pkg/i-tasser/5.1/I-TASSERmod/runI-TASSER.pl \
-runstyle gnuparallel -pkgdir /apps/pkg/i-tasser/5.1 \
-libdir /users/jpatte75/ITLIB \
-LBS false -EC true -GO false \
-homoflag real -idcut $cutOff \
-nmodel 6 -ntemp 40 \
-seqname $1 \
-datadir "$SLURM_SUBMIT_DIR"/$1/iTasser \
-java_home /usr -light true \
-outdir "$SLURM_SUBMIT_DIR"/$1/iTasser > output_${systemName}_${SLURM_JOB_ID}.log

### Alphafold2

#module load alphafold/2.1.1

#FASTA_FILE="seq.fasta"
#MAX_DATE="2022-05-14"

#cd $SLURM_SUBMIT_DIR
#echo -e "\nRunning: run_alphafold.sh \\ \n  -d ${AF2_DATA} \\ \n  -o ${SLURM_SUBMIT_DIR}/out \\ \n  -f ${SLURM_SUBMIT_DIR}/${FASTA_FILE} \\ \n  -t ${MAX_DATE} \\ \n  -n $SLURM_CPUS_ON_NODE \\ \n  -a $CUDA_VISIBLE_DEVICES\n\n"

#run_alphafold.sh \
#  -d ${AF2_DATA} \
#  -o ${SLURM_SUBMIT_DIR}/out \
#  -f ${SLURM_SUBMIT_DIR}/${FASTA_FILE} \
#  -t ${MAX_DATE} \
#  -n $SLURM_CPUS_ON_NODE \
#  -a $CUDA_VISIBLE_DEVICES


########################################## make models
wait
cd ..
cp ./iTasser/model*pdb .
# alphafold structs??


cd zDock
cp ../iTasser/model*pdb .
cp ../zDockModels/*_m.pdb .
#cp ../../../DELLA_models_docking/* .
cp ../../../zDockScripts/* .
#cp ../zDockScipts/* .
#cp ../runZdock.slurm .
cp ../../bootstrapping_zdockScores.R .
#sbatch -W runZdock.slurm 
mark_sur model1.pdb model1_m.pdb
mark_sur model2.pdb model2_m.pdb
mark_sur model3.pdb model3_m.pdb
mark_sur model4.pdb model4_m.pdb
mark_sur model5.pdb model5_m.pdb
mark_sur model6.pdb model6_m.pdb
mark_sur model7.pdb model7_m.pdb
mark_sur model8.pdb model8_m.pdb
mark_sur model9.pdb model9_m.pdb
mark_sur model10.pdb model10_m.pdb

# replace with names like in run_specific
for VARIABLE in 1 2 3 4 5 6 7 8 9 10
do
zdock -R RIN4_*model${VARIABLE}_m.pdb -L model1_m.pdb -N 200 -o model1_RIN4_Model${VARIABLE}.out
zdock -R RIN4_*model${VARIABLE}_m.pdb -L model2_m.pdb -N 200 -o model2_RIN4_Model${VARIABLE}.out
zdock -R RIN4_*model${VARIABLE}_m.pdb -L model3_m.pdb -N 200 -o model3_RIN4_Model${VARIABLE}.out
zdock -R RIN4_*model${VARIABLE}_m.pdb -L model4_m.pdb -N 200 -o model4_RIN4_Model${VARIABLE}.out
zdock -R RIN4_*model${VARIABLE}_m.pdb -L model5_m.pdb -N 200 -o model5_RIN4_Model${VARIABLE}.out
zdock -R RIN4_*model${VARIABLE}_m.pdb -L model6_m.pdb -N 200 -o model6_RIN4_Model${VARIABLE}.out
zdock -R RIN4_*model${VARIABLE}_m.pdb -L model7_m.pdb -N 200 -o model7_RIN4_Model${VARIABLE}.out
zdock -R RIN4_*model${VARIABLE}_m.pdb -L model8_m.pdb -N 200 -o model8_RIN4_Model${VARIABLE}.out
zdock -R RIN4_*model${VARIABLE}_m.pdb -L model9_m.pdb -N 200 -o model9_RIN4_Model${VARIABLE}.out
zdock -R RIN4_*model${VARIABLE}_m.pdb -L model10_m.pdb -N 200 -o model10_RIN4_Model${VARIABLE}.out
done


wait
#grep someting out of summary shell
./summarize_zdock.sh
Rscript bootstrapping_zdockScores.R > outputfile.txt
cp outputfile.txt ../
cd ..
#rm -r zDock
#rm -r iTasser

