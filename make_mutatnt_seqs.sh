#!/bin/bash
#SBATCH --job-name       DivSeqs
#SBATCH --time           14:00:00      # Allow 100 CPU hrs / GB of blastn query seq
#SBATCH --ntasks         1
#SBATCH --cpus-per-task  18            # 1 whole node.
##SBATCH --mem            100G          # 1 whole node. Allow for whole database.
#SBATCH --error          error_DivSeqs.err
#SBATCH --output         DivSeqs.log
cd $SLURM_SUBMIT_DIR

module load emboss
module load anaconda3
for VARIABLE in 0 1 2 3 4 5 6
        do

cd Footprint_${VARIABLE}/
cp ../../PSS2Seqs.py .
water -sprotein footprint_*.fasta -bsequence homologs.fasta -gapextend 0.5 -gapopen 10 -aformat a2m -outfile my_aln_seqs
python3 PSS2Seqs.py
cd ..
	done
python3 mutantFootprints_2.py
