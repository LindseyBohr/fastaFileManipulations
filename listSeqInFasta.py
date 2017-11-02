#!/usr/bin/env python

import sys
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC,Gapped
from operator import itemgetter

###############
# This script lists sequence ids in a fasta alignment
# should eventually update for a more universal use
# right now have to enter fasta name in file
##################

if len(sys.argv) !=2:
    print("Usage: ListSeqinFasta.py <multiFASTAalignment>")
    sys.exit(0)

for seq_record in SeqIO.parse(sys.argv[1], "fasta"):
    print(seq_record.id)

