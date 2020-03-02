#!/usr/bin/env python

import sys
import os
from Bio import SeqIO
import argparse
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC,Gapped
import numpy as np

###############
# This script randomly subsamples an aln
##################

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
    parser = argparse.ArgumentParser(description='subsample fasta')
    parser.add_argument("input_fasta", help="OG fasta aln",
        action=FullPaths,
        type=is_file)
    return parser.parse_args()

def remove_isolates(toKeep_list,input_fasta):
    """this removes isolates from the fasta and saves them as individual fastas"""
    sequences_to_keep = []
    tot_sequences = []
    for seq_record in SeqIO.parse(input_fasta, "fasta"):
        tot_sequences.append(seq_record)
    print(toKeep_list)
    for i,seq_record in enumerate(tot_sequences):
        print i
        if i not in toKeep_list:
            print("removed "+seq_record.id)
        else:
            sequences_to_keep.append(seq_record)
    return sequences_to_keep

def print_out_aln(keepList,out_handle):
    output_handle = open(out_handle,"w")
    SeqIO.write(keepList, output_handle, "fasta")
    output_handle.close()

args = get_args()

Start = 0
Stop = 44 # 45 starting isolates
limit = 27 # subsampling to 27 isolates

randomKeepList = np.random.choice(Stop,limit,replace = False)
seq_to_keep_list = remove_isolates(randomKeepList, args.input_fasta)
print_out_aln(seq_to_keep_list,"subsampled.fasta")


