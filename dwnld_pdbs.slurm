#!/bin/bash

#SBATCH --job-name=download_pdbs
#SBATCH --partition=DTN
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --output=a_dwnld.out
#SBATCH --error=a_dwnld.err
#SBATCH --time=4:00:00

cd $SLURM_SUBMIT_DIR

#module load sra-tools/2.11.0
#SRR=SRR2121685
#OPTS="--gzip --skip-technical  --readids --read-filter pass --dumpbase --split-3"
#echo -e "\nRunning: prefetch $SRR ..."
#prefetch $SRR
#echo -e "\nRunning: fastq-dump $OPTS --outdir fastq $SRR/$SRR.sra ..."
#fastq-dump $OPTS --outdir fastq $SRR/$SRR.sra

#./batch_download.sh -f /dev/stdin <<<'1HTZ' -p

grep -oP '(?<=pdb\|)(?s).*(?=\|)' test_homologs.fasta > Pdb_files.txt
#sed 's/"//g' < PDB_files > PDB_files.txt
sed -i '$!s/$/,/' Pdb_files.txt
tr --delete '\n' < Pdb_files.txt > PDB_files.txt
rm Pdb_files.txt
./batch_download.sh -f PDB_files.txt -p
wait
mkdir structures
mkdir ./structures/PDB
gunzip *.pdb*
mv *.pdb ./structures/PDB
