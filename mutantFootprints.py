import os
import re
import sys
import math
import random
import itertools
import Bio.PDB
import numpy as np
import collections
import pandas as pd
from Bio import AlignIO
from Bio import SeqIO
from Bio import Seq
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from optparse import OptionParser
#Author: JP
#Date Last Edit:19.9.19
#Disc:  Takes output from WindowFootprint.py and PSS2Seqs.py to create a fasta file
#       with new mutants variants from the original footprint
# Assume -- needs the logicals to recombine inside the fasta file in disc. section

# .
# --------------------------------------------------        Parse the options
# Unsure how to implement at the current moment. Could add specific 'range' of mutant outputs to use
# usage = "USAGE: python WindowSequence.py --f1 sequenceOfInterest --c1 window  [--o1 output1] \n"
# parser = OptionParser(usage=usage)
# 

def combineParts(fragments, parts):
    newFootprints = []
    previousFootprints = []
    count = 0
    for i in range(parts):
        if i+1 == 1:
            for sequence in fragments[i+1]:
                previousFootprints.append( sequence  )
                count = count + 1
        else:
            count = 0
            for j in range(len( previousFootprints )):
                # for k in range(len( fragments[i+1] )):
                for sequence in fragments[i+1]:
                    print(str(previousFootprints[j].seq))
                    print(str( sequence.seq ))
                    # newFootprints.append( SeqRecord( Seq( str(previousFootprints[j]) + str(fragments[i+1][k].seq) ),'footprint_%i' % (count), '', '')  )
                    newFootprints.append( SeqRecord( Seq( str(previousFootprints[j].seq) + str( sequence.seq ) ),'footprint_%i' % (count), '', '')  )
                    count = count + 1
            print(count)
            count = 0 
            previousFootprints = newFootprints
            newFootprints = []
        # print(previousFootprints)
    return previousFootprints


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
# 		out_folder = "out_seq_"+str(random.randint(0,999999))
# 		print( "No job_id supplied, creating a folder with a random name: "+out_folder)
# else:
# 	print( "Not enough input arguments supplied")
# 	print( usage)
# 	sys.exit()        
# --------------------------------------------------   	
# runtime : run mutantFootprints.py --f1 SeqOfInterest2.fa --c1 6
cwd = os.getcwd()
# if os.path.isdir(out_folder) == False:
#     os.makedirs(out_folder)
# cwd = cwd + "\\" + out_folder
dirList = os.listdir()
# logicals = []
# for dirName in dirList:
#     print( str.split(dirName, 'Footprint_')[1] )
#     logicals.append(str.split(dirName, 'Footprint_')[1])
# start making the new dict of seq for pdbaa blastp 
newMutantFootprints = {}
fragments = {}
part = 0
setionCount = 1
collectFragments = True
createFootprints = False
for i in range(len(dirList)):
    cwdTemp = cwd + '\\' + 'Footprint_%i' % (i+1)
    os.chdir(cwdTemp)
    seqOfInterest = SeqIO.parse('mutated_sequences.fasta', 'fasta')
    seqOfInterest = list(seqOfInterest)
    
    if createFootprints == True:
        part = part + 1 
        fragments[part] = seqOfInterest
        if i+1 == setionCount * 2 + 1 :
                # seqOfInterest = list(itertools.permutations(fragments, axis=1) )
                #itertools.combinations(iterable, r)
                seqOfInterest = combineParts(fragments,part)
                setionCount = setionCount * 2 + 1
                collectFragments = True
                fragments = {}
                part = 0
        else:
            continue
    if collectFragments == True:
        for footprint in seqOfInterest:
            newMutantFootprints[str(footprint.seq)] = footprint.id
        collectFragments = False
        createFootprints = True
# end for
#write new footprints out
os.chdir(cwd)
new_frags=[]
j = 0
for key in newMutantFootprints:
    record = SeqRecord(Seq(key),'footprint_%i' % (j+1), '', '')
    new_frags.append(record)
    j=j+1
Bio.SeqIO.write(new_frags,"mutated_footprints.fasta", "fasta")
    
    



# seqOfInterest = SeqIO.parse(options.f1, 'fasta')
# seqOfInterest = list(seqOfInterest)
# results, splitPattern = windowSeq(str(seqOfInterest[0].seq), window)
# print(results)
# print(splitPattern)
# os.chdir( cwd)
# for i in range(len(results)):
#     dirName = "Footprint_" + "".join(map(str,splitPattern[i].astype(int)))
#     os.makedirs(dirName)
#     cwdTemp = cwd + "\\" + dirName
#     os.chdir(cwdTemp)
#     recordFootprint = SeqRecord(results[i],'footprint_%i' % (i+1), '', "".join(map(str,splitPattern[0].astype(int))) )
#     recordFootprint.seq = Seq(re.sub('[^ARNDCQEGHILKMFPSTWYV]',"",str(recordFootprint.seq).upper()))
#     Bio.SeqIO.write(recordFootprint,'footprint_%i.fasta' % (i+1), "fasta")
#     os.chdir(cwd)
    