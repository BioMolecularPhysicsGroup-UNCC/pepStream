import os
import re
import sys
import Bio.PDB
import random
import numpy as np
import pandas as pd
from optparse import OptionParser
from Bio import AlignIO
from Bio.Alphabet import IUPAC


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

# pssMatrix = pd.read_csv( "ascii_mtx_file" , sep='\t',skipfooter=5,skiprows=1,lineterminator='\n')
pssMatrix = pd.read_table("ascii_mtx_file",delim_whitespace=True, skiprows=2, skipfooter=5,header=None,index_col=0)
# pssMatrix = pd.read_table("ascii_mtx_file",skiprows=1, skipfooter=5,sep=' ')
# pssMatrix = pd.read_fwf("ascii_mtx_file",skiprows=1, skipfooter=5)
# pssMatrix = PSSMParser( "ascii_mtx_file" )

column_numbers = [x for x in range(pssMatrix.shape[1])]  # list of columns' integer indices
# [column_numbers.remove(x) for x in np.arange(0,18) ] #removing columns via array
# transitionMatrix = pssMatrix.iloc[:, column_numbers] #return all columns except the 0th column
transitionMatrix = pssMatrix.iloc[:,np.arange(21,41)] # just the frequency probabilities from the pssm
AAletters = 'ARNDCQEGHILKMFPSTWYV' #specifically the letter list used by psi blast
inputSequence = pssMatrix.iloc[:,0].astype(str).str.cat()
transitionProbs = dict.fromkeys(np.arange(0,len(pssMatrix.iloc[:,0])) ) 

#i need a dict of dicts where the primary keys are the original seq....
for i in np.arange(0,len(transitionMatrix)): 
    for j in np.arange(0,len(AAletters)):
        if transitionMatrix.iloc[i,:].iloc[j] != 0: #capture the change frequency 
            if transitionProbs[i] == None:
                transitionProbs[i]  = {}
            transitionProbs[i][AAletters[j]] = transitionMatrix.iloc[i,:].iloc[j]
            
# read and score every match seq using transition matrix dict.           
align = AlignIO.read("my-out-seqs.fa", "fasta")







