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

def combineParts(fragments, energies, parts):
    newFootprints = []
    previousFootprints = []
    previousEnergies = {}
    newEnergies = {}
    count = 0
    
    for i in range(parts):
        print("Combined stuff")
        if i+1 == 1:
            for sequence in fragments[i+1]:
                previousFootprints.append( sequence  )
                previousEnergies[ str(sequence.seq) ] = energies[i+1][str(sequence.seq)]
                count = count + 1 # important, to move to else
        else:
            count = 0
            for ii in range(len( previousFootprints )):
                # for k in range(len( fragments[i+1] )):
                for sequence in fragments[i+1]:
                    #print(str( sequence.seq ))
                    newFootprints.append( SeqRecord( Seq( str(previousFootprints[ii].seq) + str( sequence.seq ) ),'footprint_%i' % (count), '', '')  )
                    newEnergies[ newFootprints[-1].seq ] = previousEnergies[str(previousFootprints[ii].seq)] +  energies[i+1][sequence.seq]
                    count = count + 1
            print(count)    
            count = 0 
            previousEnergies = newEnergies
            previousFootprints = newFootprints
            newFootprints = []
        # print(previousFootprints)
    return previousFootprints, newEnergies


# #options
# parser.add_option("--f1",help="sequence to window over", dest="f1")
# parser.add_option("--c1",help="window size", dest="c1")
# parser.add_option("--jobid",help="Job name for the output folder", dest="jobid")
# 
# (options, args) = parser.parse_args()
# window = 7
# out_folder = ""
# if (options.f1 and options.c1):
#     if (options.c1):
#         window = float(options.c1)
#     else:
#         print( "No window supplied, using default of 7")
#     
#     if (options.jobid):
#         out_folder = "out_"+options.jobid
#     else:
#         out_folder = "out_seq_"+str(random.randint(0,999999))
#         print( "No job_id supplied, creating a folder with a random name: "+out_folder)
# else:
#     print( "Not enough input arguments supplied")
#     print( usage)
#     sys.exit()        
# --------------------------------------------------       
# runtime : run mutantFootprints.py --f1 SeqOfInterest2.fa --c1 6
cwd = os.getcwd()
# if os.path.isdir(out_folder) == False:
#     os.makedirs(out_folder)
# cwd = cwd + "\\" + out_folder
dirList = next(os.walk('.'))[1] #os.listdir()
# logicals = []
# for dirName in dirList:
#     print( str.split(dirName, 'Footprint_')[1] )
#     logicals.append(str.split(dirName, 'Footprint_')[1])
# start making the new dict of seq for pdbaa blastp 
newMutantFootprints = {}
newMutantEnergy = {}
fragments = {}
energies = {}
listOfEnergies = {}
footprintNum = 0
part = 0
setionCount = 1
collectFragments = True
createFootprints = False
for i in range(len(dirList)-1):
    cwdTemp = cwd + '/' + 'Footprint_%i' % (i)
    os.chdir(cwdTemp)
    seqOfInterest = SeqIO.parse('mutated_sequences.fasta', 'fasta')
    seqOfInterest = list(seqOfInterest)
    for footprint in seqOfInterest:
            listOfEnergies[str(footprint.seq)] = float(footprint.description.split()[1][1:-2])
            
    if createFootprints == True:
        part = part + 1 
        fragments[part] = seqOfInterest
        energies[part] = listOfEnergies
        if i+1 == setionCount * 2 +1:
                # seqOfInterest = list(itertools.permutations(fragments, axis=1) )
                #itertools.combinations(iterable, r)
                seqOfInterest, listOfEnergies = combineParts(fragments,energies,part)
                setionCount = setionCount * 2+1
                collectFragments = True
                fragments = {}
                part = 0
        else:
            continue
            
    if collectFragments == True:
        for footprint in seqOfInterest:
            if footprint.seq not in newMutantFootprints:
                footprintNum = footprintNum + 1
                newMutantFootprints[str(footprint.seq)] = 'footprint_%i' % int(footprintNum)
                newMutantEnergy[str(footprint.seq)] = listOfEnergies[str(footprint.seq)]
        collectFragments = False
        createFootprints = True
# end for
#write new footprints out
os.chdir(cwd)
new_frags=[]
footprints_df = pd.DataFrame( columns=["Name","Energy","Sequence"] )
j = 0
for key in newMutantFootprints:
    record = SeqRecord( Seq(key),'footprint_%i' % (j+1),'', str(newMutantEnergy[key]) )
    footprints_df = footprints_df.append({"Name":'footprint_%i' % (j+1),"Energy":newMutantEnergy[key] ,"Sequence":Seq(key)},ignore_index = True)
    new_frags.append(record)
    j=j+1
    if j%50000 == 0 :
        print("Print out "+str(j))
        Bio.SeqIO.write(new_frags,"mutated_footprints"+str(j/50000)+".fasta", "fasta")
        new_frags.clear()
        #new_frags=[]
        footprints_df.to_csv(r'./mutated_footprints'+str(j/50000)+'.csv',header=None, index=None, sep='\t', mode='a')
        footprints_df.iloc[0:0]
        footprints_df = pd.DataFrame( columns=["Name","Energy","Sequence"] )

#df = footprints_df.sort_values('Energy',ascending=True)
#footprints_df.to_csv(r'./mutated_footprints.csv',header=None, index=None, sep='\t', mode='a')



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
    

