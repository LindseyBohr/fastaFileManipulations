#!/usr/bin/env python

import sys
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC,Gapped

###############
# This script separates seqences from a multiFASTA alignment
##################


multiAln = open(sys.argv[1],'r')

for seq_record in SeqIO.parse(multiAln, "fasta"):
    output_handle = open(seq_record.id+".fasta","w")
    SeqIO.write(seq_record,output_handle,"fasta")
    output_handle.close()



