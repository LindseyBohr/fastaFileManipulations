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
    parser = argparse.ArgumentParser(description='remove isoaltes from aln, and export as separate fastas')
    parser.add_argument("input_fasta", help="OG fasta aln",
        action=FullPaths,
        type=is_file)
    parser.add_argument("sequences_to_remove",
        help=".txt file with seq to remove and separate",
        action=FullPaths, type=is_file),
    parser.add_argument("out_handle",
        help="output file name")
    return parser.parse_args()

def get_seq_to_remove(in_seq):
    seq_to_remove = []
    with open(in_seq, 'r') as infile:
        for i, line in enumerate(infile):
            line = line.strip().split('\t')
            seqID = line[0]
            seq_to_remove.append(seqID)
    return seq_to_remove

def remove_isolates(toRemove_list,input_fasta):
    """this removes isolates from the fasta and saves them as individual fastas"""
    sequences_to_keep = []
    for seq_record in SeqIO.parse(input_fasta, "fasta"):
        print(seq_record.id)
        if seq_record.id in toRemove_list:
            print("removed "+seq_record.id)
            sep_out_handle = open(seq_record.id+".fasta","w")
            SeqIO.write(seq_record,sep_out_handle,"fasta")
            sep_out_handle.close()
        else:
            sequences_to_keep.append(seq_record)
    return sequences_to_keep

def print_out_aln(keepList,out_handle):
    output_handle = open(out_handle,"w")
    SeqIO.write(keepList, output_handle, "fasta")
    output_handle.close()

args = get_args()
seq_to_remove_list = get_seq_to_remove(args.sequences_to_remove)
seq_to_keep_list = remove_isolates(seq_to_remove_list,args.input_fasta)
print_out_aln(seq_to_keep_list,args.out_handle)


