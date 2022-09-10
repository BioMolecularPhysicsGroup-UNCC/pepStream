import os
import re
import sys
import random
import Bio.PDB
import numpy as np
import collections
import pandas as pd
from Bio import AlignIO
from Bio import SeqIO
from Bio import Seq
#from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from optparse import OptionParser
#Author: JP
#Date Last Edit:19.9.19
#Disc:  Takes a fata sequence, binarily breaks it into user defined chunks iteratively
#       printing to file new fasta and the logical split in new directories
# Assume --

# .
def windowSeq(seq, n):
    j = 1
    seqs = [""]
    logicals = []
    footprintLen = len(seq)
    logicalFootprint = np.ones(footprintLen)
    logicals.append(logicalFootprint)
    while len(seqs[-1]) > n or seqs[0] == "":#poor method to start/stop the parsing
        partition = 2 * j
        partitionLen = len(seq)/partition
        if j == 1:
            seqs = [seq]
        for i in range(partition):
            seqs.append( seq[int(i*partitionLen) : int(i*partitionLen+partitionLen)] )#simply take slices from original seq
            logicalFootprint = np.zeros(footprintLen)
            logicalFootprint[int(i*partitionLen) : int(i*partitionLen+partitionLen)] = 1
            logicals.append(logicalFootprint)
        j = j + 1
    return seqs, logicals

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
if (options.f1):
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
# runtime : run WindowFootprint.py --f1 SeqOfInterest2.fa --c1 5
cwd = os.getcwd()
if os.path.isdir(out_folder) == False:
    os.makedirs(out_folder)
cwd = "./" + out_folder
seqOfInterest = SeqIO.parse(options.f1, 'fasta')
seqOfInterest = list(seqOfInterest)
results, splitPattern = windowSeq(str(seqOfInterest[0].seq), window)
print(results)
print(splitPattern)
os.chdir(cwd)
for i in range(len(results)):
#        
    # dirName = "Footprint_" + "".join(map(str,splitPattern[i].astype(int)))
    dirName = "Footprint_" + str(i)
    os.makedirs(dirName)
    cwdTemp = "./" + dirName #changed for cluster use, generalize?
    os.chdir(cwdTemp)
    recordFootprint = SeqRecord(results[i],'footprint_%i' % (i+1), '', "".join(map(str,splitPattern[i].astype(int))) )
    recordFootprint.seq = Seq(re.sub('[^ARNDCQEGHILKMFPSTWYV]',"",str(recordFootprint.seq).upper()))
    Bio.SeqIO.write(recordFootprint,'footprint_%i.fasta' % (i+1), "fasta")
    os.chdir("..")

# Come back to this for logical outprinting 23/8/19
# for i in results:
#     f = open(cwd+"\\"+out_folder+"\\parameters.txt",'w')
#     f.write("The results in this folder were obtained using the following parameters:\n")
#     f.close()

# Bio.SeqIO.write(new_frags,"mutated_sequences.fasta", "fasta")

#split these into dir with logicals for the split
#ex WHOLEFRAGMENT == 111111111111
#ex WHOLE == 1111100000000
