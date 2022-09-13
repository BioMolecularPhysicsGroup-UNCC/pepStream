# ======================================================================================= #
#                                                 ###                                     #
#    File:    PSS2Seqs.py                         ###    BioMolecular Physics Group       #
#    Author:  John Patterson                      ###   University of North Carolina      #
#    Created: 04-01-2019                          ###           at Charlotte              #
#                                                 ###                                     #
# ======================================================================================= #
# DESCRIPTION:
#    Takes as input a position-specific scoring matrix (PSSM) that is output from
#    psi-blast, along with an alignment (i.e. EMBOSSS, clustal, et cetera). As output,
#    a fasta file the unique sequences with substituted AA
# --------------------------------------------------------------------------------------- #
# DEPENDENCIES:
#    Bio 1.4.0
#    numpy 1.22.0
# --------------------------------------------------------------------------------------- #
# INPUT:
#    sequenceOfInterest.fasta = String containing N characters, where N is the number
#                               of residues/AAs. Each character in the sequence of
#                               interest follows the 1-letter residue naming convention.
# --------------------------------------------------------------------------------------- #
# OUTPUT:
#    footprint_i = Directory for the i-th set of fragments. This directory is
#                  saved to local disk at the out_folder path defined below.
# --------------------------------------------------------------------------------------- #
# USAGE:
#    Bash: python WindowSequence.py --f1 sequenceOfInterest --c1 window [--o1 output1]
# --------------------------------------------------------------------------------------- #
# BASH OPTIONS:
#    f1 = Input sequence.
#    c1 = Number of desired sequence partitions/fragments.
#    o1 = Desired output file name.
# ======================================================================================= #
import os
import re
import sys
import random
import Bio.PDB
import numpy as np
import collections
import pandas as pd
from optparse import OptionParser
from Bio import AlignIO
from Bio import SeqIO
from Bio import Seq
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

#Author: JP
#Date Last Edit:16.9.19
#Disc:  Takes a pssm from psi-blast and and alignment (clustal or emboss)
#       returns in fasta file the unique sequences with substituted AA
# Assume --

def pssScore(alnString, origString, pssMatrix):
    probs = []
    origString = list(origString)
    alnString = list(alnString)
    for k in range(len(alnString)):
        if alnString[k] == '-': # take note: - is from clustal _ is emboss: actually maybe not
        # alnString.replace(alnString[k],origString[k])
            alnString[k] = origString[k]
            pass
            # print( alnString)
            # print( alnString[k])
            # print(k)
            # print( pssMatrix[k][alnString[k]] )
            probs.append( pssMatrix[k][alnString[k]] /100 )
        
        # print(k)
        if pssMatrix[k] == None:
            pssMatrix[k] = {}
        if alnString[k] not in pssMatrix[k].keys():
        # if pssMatrix[k][alnString[k]] == None:
            pssMatrix[k][alnString[k]] = 1
        probs.append( pssMatrix[k][alnString[k]] /100 )
        
    
    # print(probs)
    # print(- np.log(probs))
    # print(np.prod( - np.log(probs) )) #using prod of -log(p) will yield 0s
    score = np.sum( - np.log(probs) )
    return ''.join(alnString), score
    
# --------------------------------------------------        Parse the options
# usage = "USAGE: python WindowSequence.py --f1 sequenceOfInterest --c1 window  [--o1 output1] \n"
# parser = OptionParser(usage=usage)
# 
# #options
# parser.add_option("--f1",help="sequence to window over", dest="f1")
# parser.add_option("--c1",help="window size", dest="c1")
# parser.add_option("--jobid",help="Job name for the output folder", dest="jobid")
# 
# (options, args) = parser.parse_args()
# window = 7
# out_folder = ""
# if (options.f1 and options.c1):
# 	if (options.c1):
# 		window = float(options.c1)
# 	else:
# 		print( "No window supplied, using default of 7")
# 	
# 	if (options.jobid):
# 		out_folder = "out_"+options.jobid
# 	else:
# 		out_folder = "out_PSSGseq_"+str(random.randint(0,999999))
# 		print( "No job_id supplied, creating a folder with a random name: "+out_folder)
# else:
# 	print( "Not enough input arguments supplied")
# 	print( usage)
# 	sys.exit()        
# --------------------------------------------------   	
# runtime : run WindowSequence.py --f1 SeqOfInterest.txt --c1 4
cwd = os.getcwd()
# if os.path.isdir(out_folder) == False:
#     os.makedirs(out_folder)

pssMatrix = pd.read_table("ascii_mtx_file",delim_whitespace=True, skiprows=3, skipfooter=5,header=None,index_col=0,engine='python')
column_numbers = [x for x in range(pssMatrix.shape[1])]  # list of columns' integer indices
transitionMatrix = pssMatrix.iloc[:,np.arange(21,41)] # just the frequency probabilities from the pssm
AAletters = 'ARNDCQEGHILKMFPSTWYV' #specifically the letter list used by psi blast
inputSequence = pssMatrix.iloc[:,0].astype(str).str.cat()
transitionProbs = dict.fromkeys(np.arange(0,len(pssMatrix.iloc[:,0])) ) 

# dict of dicts where the outter keys are position from original input seq
#   and the inside keys are the transition probabilities from the pssm
for i in np.arange(0,len(transitionMatrix)): 
    for j in np.arange(0,len(AAletters)):
        if transitionMatrix.iloc[i,:].iloc[j] != 0: #capture the change frequency 
            if transitionProbs[i] == None:
                transitionProbs[i]  = {}
            transitionProbs[i][AAletters[j]] = transitionMatrix.iloc[i,:].iloc[j]
            
            
# pads aln file from emboss
records = SeqIO.parse("my_aln_seqs", 'fasta')
records = list(records) # make a copy, otherwise our generator
                        # is exhausted after calculating maxlen
# maxlen = max(len(record.seq) for record in records)
maxlen = len(records[0].seq)
# pad sequences so that they all have the same length
for record in records:
    if len(record.seq) < maxlen:
        sequence = str(record.seq).rjust(maxlen, '-')
        record.seq = Seq(sequence)
    if len(record.seq) > maxlen:
        print(record.seq)
        sequence = str(record.seq)[-maxlen-1:-1]
        record.seq = Seq(sequence)
assert all(len(record.seq) == maxlen for record in records)

# write to temporary file and do alignment
output_file = '{}_padded'.format(os.path.splitext("my_aln_seqs")[0])
with open(output_file, 'w') as f:
    SeqIO.write(records, f, 'fasta')
               
# read and score every match seq using transition matrix dict. 'transitionProbs'           
align = AlignIO.read("my_aln_seqs_padded", "fasta")
scoreList = {} 
# make scores a dict. where keys are scores and values are align/mut sequences 
#need to make these unique!!!!!!!
for i in np.arange(0,len(align) ):
    wordsOf = str(align[int(i)].seq)
    mutSeq, tempScore = pssScore(wordsOf, inputSequence, transitionProbs)
    if scoreList.values() == mutSeq:
        print("repeate sequence generated")
    scoreList[mutSeq] = [tempScore]
# oscoreList = collections.OrderedDict(sorted(scoreList.items()))
oscoreList = sorted(scoreList, key=scoreList.get)
print(oscoreList)
new_frags=[]
i=0
for key in oscoreList:
    record = SeqRecord(Seq(key),'footprint_%i' % (i+1), '', str( scoreList[key]))
    new_frags.append(record)
    i=i+1
Bio.SeqIO.write(new_frags,"mutated_sequences.fasta", "fasta")

# for recombining the sequences into full length, find homology pssm seqs, then do all combos of fragments that will be resubmitted to pdb. 

