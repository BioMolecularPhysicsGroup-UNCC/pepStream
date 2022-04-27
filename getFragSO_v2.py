try:
    import os
    import re
    import sys
    import glob
    import Bio.PDB
    import random
    import numpy as np
    import math
    import PyQt5
    import itertools
    import collections
    import pandas as pd
    # import matplotlib
    # matplotlib.use('TKAgg')
    #import matplotlib.pyplot as plt
    from itertools import groupby
    from operator import itemgetter
    
    from Bio.PDB import *
    from Bio.PDB.PDBParser import PDBParser
    from Bio.PDB.Polypeptide import PPBuilder
    from Bio.PDB import PDBList
    from Bio import AlignIO
    from Bio import SeqIO
    from Bio import Seq
    from Bio.Alphabet import IUPAC
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.PDB.Polypeptide import three_to_one
    from Bio.PDB import PDBIO
    
    from optparse import OptionParser
    
except ImportError:
        import sys
        sys.stderr.write("\nERROR: Something did not import correctly\n")
        sys.exit()

#Author: JP
#Date Last Edit:14.12.20
#Disc:  Parses unique pdb: IDs, Chains, Seq start & end from a blastp output
#Update v2 : This now also gets the significant other fragments from the parsed 
#           BLASTp output. 
# Assume -- result.out is a BLASTp output formatted specifically (see manual)
#           
#
# Runtime: simply call this script in the directory with the result.out from pBLAST
# 

cwd = os.getcwd()
os.chdir(cwd)
# seqOfInterest = SeqIO.parse('homologs.fasta', 'fasta')
# seqOfInterest = list(seqOfInterest)
# pdbList = {}
# freqList = {}
footprintStartEnd = {} 
        
fullBLASTPout = pd.read_table("results.out",header=None)

#Pulls organism name and prints to file
organismList = {}
for i in range(len( fullBLASTPout )):
    cMoRFinfo = fullBLASTPout.iloc[i]
        
    if cMoRFinfo[6][0:4] in footprintStartEnd:
            #i could look for different segments of the same crystal at this point
            #this does ensure that any organism getting into list is from non-redundant crystal
        continue
    else:
        footprintStartEnd[cMoRFinfo[6][0:4]] = [cMoRFinfo[8],cMoRFinfo[9],cMoRFinfo[6][-1],cMoRFinfo[12]]
        
    if str(cMoRFinfo[12]) in organismList:
        organismList[ cMoRFinfo[12] ] =  int(organismList[ cMoRFinfo[12] ]) + 1
    else:
        organismList[ cMoRFinfo[12] ] = 1 
        
    
# Print tabel of organisms to file
df = pd.DataFrame.from_dict(organismList, orient="index")
df.to_csv("organism_list.csv")

df = pd.DataFrame.from_dict(footprintStartEnd, orient="index")
df.to_csv("footprint_list.csv")

if "strucutres" not in os.listdir():
    os.makedirs("strucutres")
cwd = cwd + '\\' + "strucutres"
os.chdir(cwd)
# Retreive pdb files into local dir.
pdbl = PDBList()
PDBlist2 = list( footprintStartEnd.keys() )
# this turns on/off the downloading of pdb files !!!!!!!!
for i in PDBlist2:
    # urllib.request.urlretrieve('http://files.rcsb.org/download/101M.pdb', '101m.pdb')
    pdbl.retrieve_pdb_file(i,pdir='PDB',file_format="pdb")
# if "PDB" not in os.listdir():
#     os.makedirs("PDB")
cwd = cwd + '\\' + 'PDB'
os.chdir(cwd)
os.system("cp ../../fragSO.py .")
for key in footprintStartEnd.keys():
    startRes = str( footprintStartEnd[key][0] )
    endRes = str( footprintStartEnd[key][1] )
    chain = str( footprintStartEnd[key][2] )
    pdbName = 'pdb'+str(key).lower()+'.ent'
    outputName = str(key) + '_output'
    
    #os.system('python fragSO.py --f1 '+ pdbName + ' --s1 '+ startRes +' --s2 '+ endRes +' --c1 '+ chain +' --c 5.5 --i 0 --jobid '+ outputName)
    

# plot a histogram of unique names and their frequency 
# Bring some raw data.
# frequencies = organismList.values()
# # In my original code I create a series and run on that, 
# # so for consistency I create a series from the list.
# freq_series = pd.Series(list(frequencies))
# logicalsToPlot = freq_series > 5
# x_labels = pd.Series(list(organismList.keys()))
# for j in range(len( logicalsToPlot )):
#     if logicalsToPlot[j] == False:
#         freq_series[j].pop()
#         x_labels[j].pop()
# Plot the figure.
# plt.figure(figsize=(12, 6))
# # plt.tight_layout()
# plt.subplots(constrained_layout=True)
# ax = freq_series[logicalsToPlot ].plot(kind='bar')
# ax.set_title('Organism Frequency')
# ax.set_xlabel('Organism')
# ax.set_ylabel('Frequency')
# # ax.set_xticklabels(x_labels[logicalsToPlot ])
# ax.legend(freq_series[logicalsToPlot ],x_labels[logicalsToPlot ], bbox_to_anchor=(0, 1), loc="lower center", mode="expand")
# # ax.get_xaxis().set_visible(False)
# rects = ax.patches
# plt.show()   

# Make some labels.
# labels = ["label%d" % i for i in range(len(rects))]

# for rect, label in zip(rects, labels):
#     height = rect.get_height()
#     ax.text(rect.get_x() + rect.get_width() / 2, height + 5, label,
#             ha='center', va='bottom')     
# plot a histogram of the pdbIDs and their frequency
# plot freezes fucking computer

# plt.bar(freqList.keys(), freqList.values())
# plt.xticks(rotation=90)
# plt.show()






