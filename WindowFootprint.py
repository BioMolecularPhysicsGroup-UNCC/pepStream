# ======================================================================================= #
#                                                 ###                                     #
#    File:    WindowFootprint.py                  ###    BioMolecular Physics Group       #
#    Author:  John Patterson                      ###   University of North Carolina      #
#    Created: 03-01-2019                          ###           at Charlotte              #
#                                                 ###                                     #
# ======================================================================================= #
# DESCRIPTION:
#    Takes a fasta sequence file, partitions it into user-defined segments then
#    iteratively prints to file a new fasta and the logical split in new directories.
#    WindowFootprint.py is intended to be called from bash; as a result, the inputs
#    and output structure is described under this context.
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
import os                            # Import miscellaneous operating system interfaces
import random                        # Import random module for rng
import re                            # Import regular expression operations
import sys                           # Import system-specific parameters and functions
from optparse import OptionParser    # Import option parser module
import numpy as np                   # Import numpy module
from Bio import SeqIO                # Import required BIO modules
from Bio.Seq import Seq              # Import callable function from Seq module
from Bio.SeqRecord import SeqRecord  # Import callable function from SeqRecord module
fsep = os.sep                        # Get operating system file separator
# ======================================================================================= #
# ============================= Sequence Partition Function ============================= #
def windowSeq(seq, n):  # Function: windowSeq(inputSequence,threshold for # of fragments)
    j = 1                                      # Initialize counter
    seqs = [""]                                # Initialize empty string
    logicals = []  # Initialize list to populate with binary representation of seqs
    footprintLen = len(seq)                    # Get length of input sequence
    logicalFootprint = np.ones(footprintLen)   # Initialize vector of ones
    logicals.append(logicalFootprint)          # Add current footprint to list
    while len(seqs[-1]) > n or seqs[0] == "":  # Poor method to start/stop the parsing (TG fix)
        partition = 2 * j                      # Working partition size decreases as j increases
        partitionLen = len(seq) / partition    # Partition length is the ratio: seq/partition
        if j == 1:                             # When counter == 1
            seqs = [seq]                       # ^--> Set sequence to initial sequence of interest
        for k in range(partition):             # Loop over number of desired partitions
            seqs.append(seq[int(k * partitionLen):int(k * partitionLen + partitionLen)])      # Add i-th fragment
            logicalFootprint = np.zeros(footprintLen)  # Initialize current working binarized fragment
            logicalFootprint[int(k * partitionLen):int(k * partitionLen + partitionLen)] = 1  # Get i-th fragment
            logicals.append(logicalFootprint)          # Add current working binarized fragment to list
        j = j + 1          # Update counter
    return seqs, logicals  # Output list of sequence fragments (seqs) and binarized fragments (logicals)
# ======================================================================================= #
# --------------------------------- Parse Input Options --------------------------------- #
usage = "USAGE: python WindowSequence.py --f1 sequenceOfInterest --c1 window  [--o1 output1] \n"  # Usage for parser
parser = OptionParser(usage=usage)  # Initialize parser configuration
parser.add_option("--f1", help="sequence to window over", dest="f1")               # Option 1
parser.add_option("--c1", help="window size", dest="c1")                           # Option 2
parser.add_option("--jobid", help="Job name for the output folder", dest="jobid")  # Option 3
(options, args) = parser.parse_args()        # Execute parsing of inputs
window = 7                                   # Set default number of partitions/fragments.
out_folder = ""                              # Initialize empty output folder string
if options.f1:                               # If user provided input 1
    if options.c1:                           # If user provided input 2
        window = float(options.c1)           # ^-- Set option c1
    else:
        print("No window supplied, using default of 7")
    if options.jobid:                        # If output file name is provided by user
        out_folder = "out_" + options.jobid  # ^--> Create output folder string
    else:  # Otherwise, generate random nomenclature for output folder
        out_folder = "out_seq_" + str(random.randint(0, 999999))  # Set output file name randomly
        print("No job_id supplied, creating a folder with a random name: " + out_folder)
else:
    print("Not enough input arguments supplied")
    print(usage)
    sys.exit()  # Exit process
# =========================================================================================#
# ----------------------------- Execute Sequence Partitioning -----------------------------#
cwd = os.getcwd()                                 # Get operating system current working directory (cwd)
if not os.path.isdir(out_folder):                 # If output folder does not exist
    os.makedirs(out_folder)                       # ^--> Make output directory using random name
cwd = "./" + out_folder                           # Get current working path to output directory
seqOfInterest = SeqIO.parse(options.f1, 'fasta')  # Parse input sequence into .fasta format
seqOfInterest = list(seqOfInterest)               # Convert input sequence to list
results, splitPattern = windowSeq(str(seqOfInterest[0].seq), window)  # Execute sequence partition function
print(results)
print(splitPattern)
os.chdir(cwd)                            # Change directory on current operating system
for i in range(len(results)):            # Loop over number of footprints
    dirName = "Footprint_" + str(i + 1)  # Current footprint output directory name
    os.makedirs(dirName)                 # Create current footprint output directory
    cwdTemp = "./" + dirName             # Changed for cluster use, generalize? (TG fix)
    os.chdir(cwdTemp)                    # Change directory on current operating system
    recordFootprint = SeqRecord(results[i], 'footprint_%i' % (i + 1), '',
                                "".join(map(str, splitPattern[i].astype(int))))
    recordFootprint.seq = Seq(re.sub('[^ARNDCQEGHILKMFPSTWYV]', "", str(recordFootprint.seq).upper()))
    SeqIO.write(recordFootprint, 'footprint_%i.fasta' % (i + 1), "fasta")  # Write current footprint
    os.chdir("..")  # Go back one directory level
# =========================================================================================#