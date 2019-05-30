#!/usr/bin/env python

import sys
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC,Gapped

###############
# this script separates a multi fasta into alns of defined clades
##################

in_fasta = sys.argv[1]
groups = sys.argv[2]

c = []
jm = []
ds6 = []
other = []

c_out = open("clade_c.fasta", "w")
jm_out = open("clade_jm.fasta","w")
ds6_out = open("clade_ds6.fasta","w")
other_out = open("clade_other.fasta","w")


with open(groups, 'r') as infile:
    for i, line in enumerate(infile):
        line = line.strip().split('\t')
        if i > 0:
            seqID = line[0]
            bf_pheno = line[1]
            clade = line[2]
            if clade == "c":
                c.append(seqID)
            elif clade == "ds6":
                ds6.append(seqID)
            elif clade == "jm":
                jm.append(seqID)
            elif clade == "other": 
                other.append(seqID)
            else:
                pass

print("number of c in tree is "+str(len(c)))
print("number of jm in tree is "+str(len(jm)))
print("number of ds6 in tree is "+str(len(ds6)))
print("number of other in tree is "+str(len(other)))

c_seq = []
jm_seq = []
ds6_seq = []
other_seq = []

for seq_record in SeqIO.parse(in_fasta, "fasta"):
    print(seq_record.id)
    if seq_record.id in c:
        c_seq.append(seq_record)
    elif seq_record.id in jm:
        jm_seq.append(seq_record)
    elif seq_record.id in ds6:
        ds6_seq.append(seq_record)
    elif seq_record.id in other:
        other_seq.append(seq_record)

SeqIO.write(c_seq, c_out, "fasta")
SeqIO.write(jm_seq, jm_out, "fasta")
SeqIO.write(ds6_seq, ds6_out, "fasta")
SeqIO.write(other_seq, other_out, "fasta")

c_out.close()
jm_out.close()
ds6_out.close()
other_out.close()


