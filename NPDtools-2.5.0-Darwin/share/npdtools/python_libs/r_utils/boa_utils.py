#!/usr/bin/env python

from common import warning
from fasta_utils import write_fasta

'''
Usage:
python boa2fasta.py boa_annotated.txt [converted_output.fasta]

Convert BOA result into a fasta file 
(only peptide sequences of bacteriocins are extracted).
If an output filepath is not specified, 
the output file will be created near the input file (with .boa.fasta extension)  
'''

NAME_COLUMN = 'bacteriocin_name'
SEQ_COLUMN = 'sequence_of_bacteriocin'


def boa2fasta(boa_fpath, converted_fpath=None):
    if converted_fpath is None:
        converted_fpath = boa_fpath + '.boa.fasta'
    with open(boa_fpath, 'r') as fr:
        header = fr.readline().strip().split()
        if NAME_COLUMN not in header or SEQ_COLUMN not in header:
            warning("\t\tUnexpected format of BOA file (%s): missing columns %s/%s! Skipping this file.."
                  % (boa_fpath, NAME_COLUMN, SEQ_COLUMN))
            return None
        name_cln_idx = header.index(NAME_COLUMN)
        seq_cln_idx = header.index(SEQ_COLUMN)
        fasta_entries = []
        for line in fr:
            row = line.strip().split()
            name = row[name_cln_idx]
            seq = row[seq_cln_idx]
            fasta_entries.append((name, seq))
        write_fasta(fasta_entries, converted_fpath)
    return converted_fpath



