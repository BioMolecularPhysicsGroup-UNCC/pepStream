#!/usr/bin/env python

from Bio import *
from Bio import SeqIO
for current_seq in SeqIO.parse("full_output_DELLA_MoRF2.fasta", "fasta"):
    if str(current_seq.seq) != '':
        print( ">", current_seq.id, "\n", str(current_seq.seq))


# Combine a bunch of fasta files for whatever....
# apparently this is what i do
cat *.fasta > full_output_DELLA_MoRF1.fasta
sed 's/ *$//g' in.fasta > out.fasta
awk 'BEGIN {RS = ">" ; FS = "\n" ; ORS = ""} $2 {print ">"$0}' input.fas > output.fas

