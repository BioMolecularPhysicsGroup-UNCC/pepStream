#!/bin/bash
#SBATCH --job-name	 pBLAST_and_Frag
#SBATCH --time           14:00:00      # Allow 100 CPU hrs / GB of blastn query seq
#SBATCH --hint           nomultithread # Unless trying 72 threads
#SBATCH --ntasks         1
#SBATCH --cpus-per-task  18            # 1 whole node.
##SBATCH --mem            100G          # 1 whole node. Allow for whole database.
#SBATCH --error		 error_BLAST.err
#SBATCH --output	 shortblast.log
##SBATCH --partition	 bigmem		#mahuika02
cd $SLURM_SUBMIT_DIR

module load blast/2.9.0+
module load anaconda3

#module load BLASTDB/2019-01
cp ../fragSO.py .
cp ../getFragSO_v2.py .
cat mutated_footprints*.0.fasta > mutated_footprints.fasta

FORMAT="6 qseqid qstart qend qseq sseqid sgi sacc slen sstart send sseq staxids sscinames stitle length evalue bitscore nident mismatch"
BLASTOPTS="-evalue 100 -word_size 3 "
#-max_target_seqs 1000"
BLASTAPP=blastp
DB=$BLASTDB/pdbaa/pdbaa  #swissprot #pdbaa
FILE=mutated_footprints.fasta
# Keep the database in RAM if searching multiple 
# query sequences against a database of over 10GB.
#cp $BLASTDB/{$DB,taxdb}* $TMPDIR/ 
#export BLASTDB=$TMPDIR

# Single node multithreaded BLAST.
#srun $BLASTAPP $BLASTOPTS -db $DB -query footprintInput.fasta -outfmt "$FORMAT" \
#    -out results.out -num_threads $SLURM_CPUS_PER_TASK

#srun $BLASTAPP $BLASTOPTS -db $DB -query $FILE -outfmt "$FORMAT" \
#    -out cMoRF_results.out -num_threads $SLURM_CPUS_PER_TASK
#srun $BLASTAPP $BLASTOPTS -db $DB -query $FILE -outfmt '6 sseqid sseq' | awk 'BEGIN{FS="\t"; OFS="\n"}{gsub(/-/, "", $2); print ">"$1,$2}' \
#    > cMoRF_homologs.fasta

#srun blastn -db nr -word_size 6 -query footprintInput.fasta -outfmt 8 -out results.out -perc_identity 100 -max_target_seqs 100 
python3 getFragSO_v2.py
