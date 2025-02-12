#!/usr/bin/python

import argparse
import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

################################################################################
# This script removes stop codon positions from all isolates in fasta aln
################################################################################

class FullPaths(argparse.Action):
    """Expand user- and relative-paths"""
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest,
            os.path.abspath(os.path.expanduser(values)))

def listdir_fullpath(d):
    return [os.path.join(d, f) for f in os.listdir(d)]

def is_dir(dirname):
    """Checks if a path is a directory"""
    if not os.path.isdir(dirname):
        msg = "{0} is not a directory".format(dirname)
        raise argparse.ArgumentTypeError(msg)
    else:
        return dirname

def is_file(filename):
    """Checks if a file exists"""
    if not os.path.isfile(filename):
        msg = "{0} is not a file".format(filename)
        raise argparse.ArgumentTypeError(msg)
    else:
        return filename

def get_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='Remove stop codons from fasta alignment')
    parser.add_argument("aln", help="original_alignment",
        action=FullPaths,
        type=is_file)
    return parser.parse_args()


args = get_args()

codon_stop_array = ["tag", "TAG", "tga", "TGA", "taa", "TAA", "uga", "UGA", "uaa", "UAA", "uag", "UAG"]

codonIndices = []

inFile = args.aln
for record in SeqIO.parse(inFile, "fasta"):
    tempRecordSeq = list(record.seq)
    for index in range(0, len(record.seq), 3):
        codon = record.seq[index:index+3]
        if codon in codon_stop_array:
            codonIndices.append([index,index+3])

# remove duplicate indices
noDupes = []
[noDupes.append(x) for x in codonIndices if x not in noDupes]

# sorts in descending order to remove last to first in alignment
# so indices don't get messed up
desc = sorted(noDupes,reverse=True)
print(desc)

toParse = SeqIO.parse(inFile, "fasta")

### this is def the part that takes the longest ........
"""
for codon in desc:
    print(codon)
    startIdx = codon[0]
    endIdx = codon[1]
    newRecords = []
    for record in toParse:
        tempRecordSeq = list(record.seq)
        del tempRecordSeq[startIdx:endIdx]
        record.seq = Seq("".join(tempRecordSeq))
        newRecords.append(record)
    SeqIO.write(newRecords,"temp.fasta","fasta")
    toParse = SeqIO.parse("temp.fasta","fasta")
"""
#My version
newRecords = []

for record in toParse:
    for codon in desc:
        if codon[1]<=len(record.seq):
            record.seq=record.seq[0:codon[0]]+record.seq[codon[1]:]
    newRecords.append(record)

SeqIO.write(newRecords,"noStopCodons.fasta","fasta")


