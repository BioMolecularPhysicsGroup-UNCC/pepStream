import os
import sys
import random
import numpy as np
from Bio import AlignIO
from Bio import SeqIO
from Bio import Seq
#from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from optparse import OptionParser

def windowSeq(seq, n):
    j = 1
    seqs = [""]
    while len(seqs[-1]) > n or seqs[0] == "":
        partition = 2 * j
        partitionLen = len(seq)/partition
        if j == 1:
            seqs = []
        for i in range(partition):
            seqs.append( seq[int(i*partitionLen) : int(i*partitionLen+partitionLen)] )
        j = j + 1
    return seqs

# --------------------------------------------------        Parse the options
usage = "USAGE: python WindowSequence.py --f1 sequenceOfInterest --c1 window  [--o1 output1] \n"
parser = OptionParser(usage=usage)

#options
parser.add_option("--f1",help="sequence to window over", dest="f1")
parser.add_option("--c1",help="window size", dest="c1")
parser.add_option("--jobid",help="Job name for the output folder", dest="jobid")

(options, args) = parser.parse_args()
window = 7
out_folder = ""
if (options.f1 and options.c1):
	if (options.c1):
		window = float(options.c1)
	else:
		print( "No window supplied, using default of 7")
	
	if (options.jobid):
		out_folder = "out_"+options.jobid
	else:
		out_folder = "out_seq_"+str(random.randint(0,999999))
		print( "No job_id supplied, creating a folder with a random name: "+out_folder)
else:
	print( "Not enough input arguments supplied")
	print( usage)
	sys.exit()        
# --------------------------------------------------   	
# runtime : run WindowSequence.py --f1 2zsh.DELLA.fasta.txt --c1 4
cwd = os.getcwd()
if os.path.isdir(out_folder) == False:
    os.makedirs(out_folder)
           
# file = open(options.f1 , 'r')
# seqOfInterest = file.read().strip() # load your sequence to array
# file.close()

seqOfInterest = SeqIO.parse(options.f1, 'fasta')
seqOfInterest = list(seqOfInterest)

results = windowSeq(str(seqOfInterest[0].seq), window)
print(results)

# Come back to this for logical outprinting 23/8/19
# for i in results:
#     f = open(cwd+"\\"+out_folder+"\\parameters.txt",'w')
#     f.write("The results in this folder were obtained using the following parameters:\n")
#     f.close()

#split these into dir with logicals for the split
#ex WHOLEFRAGMENT == 111111111111
#ex WHOLE == 1111100000000