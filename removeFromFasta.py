#!/usr/bin/env python

import sys
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC,Gapped

###############
# This script removes seqences from a FASTA alignment
# takes txt file listing isolate id's to remove from aln
##################


output_handle = open("DCCs_and_ERR009.fasta", "w")
sequences = []

seq_to_remove = []
with open(sys.argv[1], 'r') as infile:
    for i, line in enumerate(infile):
        line = line.strip().split('\t')
        seqID = line[0]
        seq_to_remove.append(seqID)

for seq_record in SeqIO.parse("ESX_concatenated_alignment.fasta", "fasta"):
    print(seq_record.id)
    if seq_record.id in seq_to_remove:
        print("removed "+seq_record.id)
    else:
        print("we got one!")
        sequences.append(seq_record)

SeqIO.write(sequences, output_handle, "fasta")
output_handle.close()

