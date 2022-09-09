import os
import sys
import numpy as np
import random
from Bio import *
from os import path
import pandas as pd
from Bio import SeqIO
from optparse import OptionParser
# script to work up numerical output from physcial filter in pepStream
usage = "USAGE: python format_output.py --f1 NameOfFile --cmorfs NameOfFile --outid NameOfOutput \n"
parser = OptionParser(usage=usage)
parser.add_option("--f1",help="File for input containing dock scores", dest="f1")
parser.add_option("--outid",help="Name for the output folder", dest="outid")
parser.add_option("--cmorfs",help="fasta file containing cMoRF seqs", dest="cmorfs")
# ------------
df = pd.DataFrame( columns=["Name","Mean","Std. Dev.","Sequence"] ) # columns=list({"Name","Mean","Std. Dev."})
(options, args) = parser.parse_args()
rawFile = ""
out_folder = ""
if (options.f1):
	if path.exists(options.f1):
		rawFile = options.f1
	else:
		print("File not found, check spelling")

	if (options.outid):
		out_folder = "out_"+options.outid
	else:
		out_folder = "out_"+str(random.randint(0,999))+"_sortReport"
		print( "No job_id supplied, creating a folder with a random name: "+out_folder)
else:
	print( "Not a valid file,file name is " + options.f1 + " , path logic " + str(path.exists(options.f1)) )
	sys.exit()

tempName = ""
tempMean = 0
tempSD = 0
index=0
seqDict={}
# --------------------
for current_seq in SeqIO.parse(options.cmorfs, "fasta"):
	seqDict[current_seq.id[-13:]]= current_seq.seq
for line in open(rawFile, "r"):
	if line.startswith("==>"):
		if line != '':
			tempName = str(line[6:19])
	if line.startswith("[1]"):
		tempList = line.split()
		if tempList[1] != '"':
			try:
				tempMean = float(tempList[1])
				tempSD = float(tempList[2])
				tempSeq = seqDict[tempName]
# pd.concat([series1, series2])
				df = df.append({"Name":tempName,"Mean":tempMean,"Std. Dev.":tempSD,"Sequence":tempSeq},ignore_index = True)
			except:
				continue

df = df.sort_values('Mean',ascending=False)
df.to_csv(r'./'+out_folder,header=None, index=None, sep='\t', mode='a')

# topTenPercent_cMoRFs.txt
length=len(df)
tenPercent=int(np.ceil(length*0.1))
#print(tenPercent)
#print(df.head(tenPercent)["Name"])
#print(df["Name"].head(tenPercent))
df.head(tenPercent)["Name"].to_csv(r'topTenPercent_cMoRFs.txt',header=None, index=None, sep='\t', mode='a')

# theBottomTenPercent_cMoRFs.txt
#print(df["Name"].tail(tenPercent))
df.tail(tenPercent)["Name"].to_csv(r'theBottomTenPercent_cMoRFs.txt',header=None, index=None, sep='\t', mode='a')

print('fin')


