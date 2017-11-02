#!/usr/bin/env python

import sys
import os
from Bio import SeqIO
from operator import itemgetter

if len(sys.argv) !=2:
    print("Usage: alignmentLength.py <multiFASTAalignment>")
    sys.exit(0)

for rec in SeqIO.parse(sys.argv[1], 'fasta'):
    name = rec.id
    seq = rec.seq
    seqLen = len(rec)
    print name, seqLen

