#!/bin/bash
#SBATCH --job-name	 Short_BLAST_Run
#SBATCH --time           02:30:00      # Allow 100 CPU hrs / GB of blastn query seq
#SBATCH --hint           nomultithread # Unless trying 72 threads
#SBATCH --ntasks         1
#SBATCH --cpus-per-task  18            # 1 whole node.
##SBATCH --mem            100G          # 1 whole node. Allow for whole database.
#SBATCH --error		 error_BLAST.err
#SBATCH --output	 shortblast.log
##SBATCH --partition	 bigmem		#mahuika02
cd $SLURM_SUBMIT_DIR

module load blast
#module load BLASTDB/2019-01

FORMAT="6 qseqid qstart qend qseq sseqid sgi sacc slen sstart send sseq staxids sscinames stitle length evalue bitscore nident mismatch"
BLASTOPTS="-evalue 1000 -word_size 3 -max_target_seqs 1000"
BLASTAPP=psiblast
DB=swissprot #pdbaa

# Keep the database in RAM if searching multiple
# query sequences against a database of over 10GB.
#cp $BLASTDB/{$DB,taxdb}* $TMPDIR/
#export BLASTDB=$TMPDIR
DB=$BLASTDB/swissprot/swissprot

for VARIABLE in 0 1 2 3 4 5 6
	do
	
	cd Footprint_${VARIABLE}/

	FILE=footprint_*.fasta 

# Single node multithreaded BLAST.
#srun $BLASTAPP $BLASTOPTS -db $DB -query footprintInput.fasta -outfmt "$FORMAT" \
#    -out results.out -num_threads $SLURM_CPUS_PER_TASK

	srun $BLASTAPP $BLASTOPTS -db $DB -out_pssm my_protein.ckp -query $FILE -outfmt "$FORMAT" \
    -out results.out -out_ascii_pssm ascii_mtx_file -num_iterations 3 -num_threads $SLURM_CPUS_PER_TASK
	srun $BLASTAPP $BLASTOPTS -db $DB -query $FILE -outfmt '6 sseqid sseq' | awk 'BEGIN{FS="\t"; OFS="\n"}{gsub(/-/, "", $2); print ">"$1,$2}' \
    > homologs.fasta
	cd ..
done

#srun blastn -db nr -word_size 6 -query footprintInput.fasta -outfmt 8 -out results.out -perc_identity 100 -max_target_seqs 100 
