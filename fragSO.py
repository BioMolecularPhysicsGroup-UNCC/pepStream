try:
    import os
    import re
    import sys
    import glob
    import Bio.PDB
    import random
    import numpy as np
    
    from itertools import groupby
    from operator import itemgetter
    
    from Bio.PDB import *
    from Bio.PDB.PDBParser import PDBParser
    from Bio.PDB.Polypeptide import PPBuilder
    from Bio.PDB.Polypeptide import three_to_one
    from Bio.PDB import PDBIO
    from Bio import SeqIO
    
    from optparse import OptionParser
    
except ImportError:
        import sys
        sys.stderr.write("\nERROR: This script requires that Biopython (http://biopython.org) is installed.\n\n")
        sys.exit()
        
#Author: JP
#Date Last Edit:12.10.19 v2
#Disc:  Enter sequence from a text file (f2) and finds sequence "signficant other"
#       Returns pdb of SO frag. & text file of number;residue 1 letter; and longest possible contiguous SO frag.
# Assume pdb accession and sequence of match is given

# add a option to allow for residues next to query to remove.
# 
def pdb2fasta(inputFile, chain_id, outputFile):
    #author: Jordan Willis ; from Rosetta modelling software

    def usage():
        print("Pulls fasta file of a specific chain up to a specified output file")
        print("python get_fasta_from_pdb_by_chain.py <pdb> <chain> <output>")
        
    # def get_three_letter(pdb, chain_id):
    #     fasta_three = []
    #     for model in pdb.get_list():
    #         for chain in model.get_list():
    #             if chain.get_id() == chain_id:
    #                 for resi in chain.get_list():
    #                     fasta_three.append(resi.get_resname())
    #     return fasta_three
    
    def get_one_letter(list_of_three):
        fasta_one=[]
        for x in list_of_three:
            if x == "S56" or x == "TES":
                    x="N"
            else:
                    x=three_to_one(x)
            fasta_one.append(x)
        return fasta_one
    
    def output_to_file(one_letter, file):
        file_handle = open(outputFile, 'a')
        file_handle.write(">" + file.split('.')[0] + "\n")
        file_handle.write("".join(one_letter) + "\n")
        file_handle.close()

    print(inputFile)
    print(outputFile)
    print(chain_id)
    
    parser = PDBParser(PERMISSIVE=1)

    handle = open(inputFile)
    pdb = parser.get_structure("random",handle)
    
    three_letter = []
    for model in pdb.get_list():
        print(model)
        for chain in model.get_list():
            print(chain)
            print(chain.get_id())
            if chain.get_id() == chain_id:
                for resi in chain.get_list():
                    three_letter.append(resi.get_resname())
    #three_letter = get_three_letter(struct, chain)
    print(three_letter)
    one_letter = get_one_letter(three_letter)
    print(one_letter)
    output_to_file(one_letter,outputFile)

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

def get_atom_list(structure):
	output = dict()
	for chain in structure:
		#if chain.id in chains:
            for residue in chain.get_residues():
                hetflag, resseq, icode=residue.get_id()
                the_id = (chain.id+"_"+str(resseq)+"_"+icode).strip()
                for atom in residue.get_unpacked_list():
                    if hetflag==' ':
                        if the_id in output:
                            output[the_id].append(atom)
                        else:
                            output[the_id] = [atom]
	return output

def get_atom_sublist(structure, chains, seqFrag):
	i = 0
	output = dict()
	for chain in structure:
		if chain.id in chains:
			for residue in chain.get_residues():
			        i=i+1
			        if i in seqFrag:
        				hetflag, resseq, icode=residue.get_id()
        				the_id = (chain.id+"_"+str(resseq)+"_"+icode).strip()
        				for atom in residue.get_unpacked_list():
           					if hetflag==' ':
          					     if the_id in output:
         			                        output[the_id].append(atom)
          					     else:
         		                                 output[the_id] = [atom]
	return output

def is_contact(res_1,other_atoms,cutoff):
	for atom in res_1:
		ns = NeighborSearch(other_atoms)
		center = atom.get_coord()
		neighbors = ns.search(center, cutoff) # 5.0 for distance in angstrom
		residue_list = Selection.unfold_entities(neighbors, 'R') # R for residues
		if len(residue_list)>0:
			return True
	return False 	

def get_contacts(struc,all_atoms,verbose,cutoff):
	progress = 0 
	contacts = []
	for residue in struc:
		progress+=1
		if len(verbose)>0:
			print( verbose,progress,"out of",len(struc))
		atom_list = struc[residue]
		outcome = is_contact(atom_list,all_atoms,cutoff)
		if outcome:
		    if residue != struc:
		      contacts.append(residue)
	return contacts				

#Filter out all the atoms from the chain,residue map given by residue_map
def get_all_atoms(residue_map):
	all_atoms_out = []
	for residue in residue_map:
		for atom in residue_map[residue]:
			all_atoms_out.append(atom)
			#Set the b-factor to zero for coloring by contacts
			atom.set_bfactor(0.0)
	return all_atoms_out
	
#Save the structures with B-factor field replaced with contact (100) and interface neighborhood (50)
def save_contacts(structure,chains,out_file):
	
	#Save only those chains that we are supposed to
        Select = Bio.PDB.Select
        class ConstrSelect(Select):
            def accept_chain(self, chain):
                #print dir(residue)
                if chain.id in chains:
                    return 1
                else:
                    return 0
   	
        w = PDBIO()
        w.set_structure(structure)
        randint = random.randint(0,9999999)
        w.save("TMP"+str(randint)+".pdb",ConstrSelect())
        #Remove the HETATM and TER lines
        f_tmp = open("TMP"+str(randint)+".pdb", 'r')
        f_out = open(out_file, 'w')
        for line in f_tmp.readlines():
            if line[0:3]!="TER" and line[0:6]!="HETATM":
                f_out.write(line)
        f_tmp.close()
        f_out.close()	
        os.remove("TMP"+str(randint)+".pdb")

#Save the residues which are contacts or neighborhood interface in a space-delimited file
def save_residues(filename,interface,contacts):
	f = open(filename,'w')	
	for elem in interface:
		splitted = elem.split("_")
		#print("Splitted is")
		#print( splitted)

		resname = str((splitted[1]+splitted[2]).strip())
		chain = splitted[0]
		#contact or neighbor of interface?
		coninf = "I" #Interface neighbor
		if elem in contacts:
			coninf = "C" #Contact
		f.write(chain+" "+resname+" "+resname+" "+coninf+"\n")
	f.close()

#only save the residues of the file that are in thelist
def save_constrained(filename_in,filename_out,thelist):
	f = open(filename_in,'r')
	f_out = open(filename_out,'w')
	for line in f.readlines():
		if "ATOM" in line:
			line = line.strip()
			resname=line[23:28].strip()
			icode = "" 
			if (not is_number(resname[len(resname)-1])):
				icode = resname[len(resname)-1]
				resname = resname[0:(len(resname)-1)]
				
			resname = line[21]+"_"+resname+"_"+icode
			if resname in thelist:
				f_out.write(line+"\n")
	f_out.close()
	
def save_contigous_constrained(filename_in,filename_out,thelist):
    contactList = {}
    for elem in thelist:
        splitted = elem.split("_")

        if splitted[2]	== "":
            if splitted[0] not in contactList:
                contactList[splitted[0]] = [splitted[1] ]
                print(contactList[splitted[0]])
                continue
            # tempList = contactList[splitted[0]]
            # print(tempList)
            contactList[splitted[0]].append( splitted[1] )
    
        else:
            if splitted[2] not in contactList:
                contactList[splitted[2]] = [splitted[1]] 
                continue
            # tempList = contactList[splitted[2]]
            contactList[splitted[2]].append( splitted[1] )
    
    for key in contactList.keys():
        contactList[key] = list(map(int, contactList[key]))
        contactList[key].sort()
    print(contactList)	 
    contigousList = {}
    numberOfPeptides = 0
    for key in contactList.keys():
        for k, g in groupby( enumerate(list(map(int, contactList[key]))), lambda ix : ix[0] - ix[1] ) :
            tempList = list( map(itemgetter(1), g) )
            print( tempList )
            print( len( tempList ) )
            if len( tempList ) > 9:#this is controll for peptide length! Make an option!!
                contigousList[numberOfPeptides] = []
                for residue in tempList:
                    contigousList[numberOfPeptides].append( key+"_"+str(residue)+"_"+"" )
                numberOfPeptides = numberOfPeptides + 1
        
    # for group in consecutive_groups(iterable):
    #     print(list(group))    
    print(contigousList)    
    for key_2 in contigousList.keys():
        f = open(filename_in,'r')
        f_out = open(filename_out+"_"+str(key_2)+".pdb",'w')
        for line in f.readlines():
            if "ATOM" in line:
                line = line.strip()
                chain=str(line[21].strip())
                resname=line[23:28].strip()
                #print(resname)
                icode = "" 
                if (not is_number(resname[len(resname)-1])):
                    icode = resname[len(resname)-1]
                    resname = resname[0:(len(resname)-1)]
    
                resname = line[21]+"_"+resname+"_"+icode
                if resname in contigousList[key_2]:
                    f_out.write(line+"\n")
        f_out.close()
        print("The chain in this is " + chain)
        if len( contigousList.keys() ) > 0:
            fileInputName=filename_out+"_"+str(key_2)+".pdb"
            fileOutputName="contigous_seq_cMoRFs_"+str(key_2)+"_"+pdbID+".fasta"
            pdb2fasta(fileInputName,chain, fileOutputName)	
       	    #os.system('python get_fasta_from_pdb.py' + ' ' + fileInputName + ' ' + chain + ' ' + fileOutputName)
			
	
# --------------------------------------------------        Parse the options
usage = "USAGE: python fragSO.py --f1 FirstPDB --s1 SeqStart --s2 SeqEnd --c1 FirstPDBChains  [--c ContactCutoff] [--i InterfaceCutoff] [--o1 output1] [--o2 output2] \n"
parser = OptionParser(usage=usage)

#options
parser.add_option("--f1",help="First molecule pdb", dest="f1")
parser.add_option("--s1",help="Start of sequence to querery contacts", dest="s1")
parser.add_option("--s2",help="End of sequence to querery contacts", dest="s2")
parser.add_option("--c1",help="Molecule chains of interest", dest="c1")
parser.add_option("--c",help="Distance cutoff between contacting residues", dest="cutoff")
parser.add_option("--i",help="Interface distance cutoff between contacting residues", dest="interfacial")
parser.add_option("--jobid",help="Job name for the output folder", dest="jobid")

(options, args) = parser.parse_args()
cutoff = 4.5
i_cutoff = 10.0
out_folder = ""
if (options.f1 and options.f1 and options.c1):
	if (options.cutoff):
		cutoff = float(options.cutoff)
	else:
		print( "No contact cutoff supplied, using default of 4.5A")
	if (options.interfacial):
		i_cutoff = float(options.interfacial)
	else:
		print( "No interfacial cutoff supplied, using default of 10.0A")
	if (options.jobid):
		out_folder = "out_"+options.jobid
	else:
		out_folder = "out_"+str(random.randint(0,999999))
		print( "No job_id supplied, creating a folder with a random name: "+out_folder)
else:
	print( "Not enough input arguments supplied")
	print( usage)
	sys.exit() 

# -----------------------------------------------------                Runtime 
# Actual input from cmd: run fragSO.py --f1 1A2Y.pdb --s1 SeqStart --s2 SeqEnd --c1 B --c 5.5 --i 0 --jobid trial_output

# get input data, set output data
cwd = os.getcwd()
if os.path.isdir(out_folder) == False:
    os.makedirs(out_folder)
pdbID = str(options.f1[3:7])
struct_1 = PDBParser().get_structure('first_one', options.f1) # load your molecule
chains_1 = options.c1 
chainList = list()
ppb=PPBuilder() # generate polypeptides, for sequence 
for pp in ppb.build_peptides(struct_1):
    chainList = pp.get_sequence() 

# may not be needed; might turn it into a read in for teh seq # of the fragement

# file = open(options.f2 , 'r')
# seqOfInterest = file.read().strip() # load your sequence to array
# file.close()

# Check type of strucutre: Yes XRD, No CryoEM, No NMR?
methodType=struct_1.header['structure_method'];
print(methodType)
if any([x in methodType for x in ["NMR","nmr"]]):
    print(struct_1.header["name"])
    print("NMR strucutre, will take first frame only")
    
    
if "X-RAY" in methodType:
    print(struct_1.header["name"])
    
    
if "ELECTRON MICROSCOPY" in methodType:
    print("Cryo EM not handled by design")
    exit()

#Load the structure - 
atoms_1 = Selection.unfold_entities(struct_1[0], 'C') # C for chains, in pdb input
chain_ids = []	
for chain in struct_1[0]:
    chain_ids.append(chain.get_id())	


#get the mapping from chain,residue id to the atom lists
#this is where I get the matching fragment for finding its SO
startRes = int( options.s1)
endRes = int( options.s2 )
input_1 = get_atom_sublist(atoms_1,chains_1,np.arange(startRes,endRes,1))
input_2 = get_atom_list(atoms_1)

for residue in input_1: del input_2[residue]
# sanity check
if len(input_1) == 0:
	print( "Are you sure that file "+options.f1+" has chain(s) "+options.c1+"? Quitting.")
	quit()
	


#get the full atom lists for neighbor searchs
all_atoms_1 = get_all_atoms(input_1)
all_atoms_2 = get_all_atoms(input_2)

#run neighbor search on both instances - not optimal but good enough for most imaginable applications.

contacts_1 = get_contacts(input_2,all_atoms_1,"",cutoff)

contact_map_1 = []
#Color molecule 1
for atom in all_atoms_1:
    atom.set_bfactor(10.0)
        
for residue in contacts_1:
	for atom in input_2[residue]:
		atom.set_bfactor(100.0)
		contact_map_1.append(atom)
# colour the contacts	  
for residue in contacts_1:
	for atom in input_2[residue]:
		atom.set_bfactor(100.0)
		contact_map_1.append(atom)

# #Color molecule 2
# contact_map_2 = []
# for residue in contacts_2:
# 	for atom in input_2[residue]:
# 		atom.set_bfactor(100.0)
# 		contact_map_2.append(atom)

#Get interfacial residues
# I'm not even sure I need this. Kinda wanna drop the whole thing...
#run neighbor search on both instances - not optimal but good enough for most imaginable applications.
interface_1 = contacts_1
# interface_2 = contacts_2
if (i_cutoff>0):
	interface_1 = get_contacts(input_1,contact_map_1,"First molecule, interfacial residue ",i_cutoff)
	# interface_2 = get_contacts(input_2,contact_map_2,"Second molecule, interfacial residue ",i_cutoff)

#Color molecule 1
# for residue in interface_1:
# 	if residue in contacts_1:
# 		continue
# 		
# 	for atom in input_1[residue]:
# 		atom.set_bfactor(50.0)

#Save the contact-colored structures
save_contacts(struct_1[0],chain_ids,cwd+"\\"+out_folder+"\\molecule_1.pdb")

#save the interface residues
save_residues(cwd+"\\"+out_folder+"\\molecule_1.txt",interface_1,contacts_1)

f = open(cwd+"\\"+out_folder+"\\parameters.txt",'w')
f.write("The results in this folder were obtained using the following parameters:\n")
f.write("Contact distance: "+str(cutoff)+"A\n")
f.write("Interface distance: "+str(i_cutoff)+"A\n")
f.close()

#Constrain the pdb files to the interface residues only
save_constrained(cwd+"\\"+out_folder+"\\molecule_1.pdb",cwd+"\\"+out_folder+"\\molecule_1_constrained.pdb",interface_1)
save_contigous_constrained(cwd+"\\"+out_folder+"\\molecule_1.pdb",cwd+"\\"+out_folder+"\\molecule_contigous",interface_1)
print( "Results saved in "+cwd+"\\"+out_folder)

