#!/usr/bin/env python

import sys
import config
import random
import re


def read_fasta(fpath=None):
    """
        Returns list of FASTA entries (in tuples: name, seq)
    """
    first = True
    seq = []
    name = ''

    if fpath:
        fasta_file = open(fpath)
    else:
        fasta_file = sys.stdin
    for raw_line in fasta_file:
        if raw_line.find('\r') != -1:
            lines = raw_line.split('\r')
        else:
            lines = [raw_line]
        for line in lines:
            if not line:
                continue
            if line[0] == '>':
                if not first:
                    yield name, "".join(seq)

                first = False
                name = line.strip()[1:]
                seq = []
            else:
                seq.append(line.strip())
    if name or seq:
        yield name, "".join(seq)
    if fpath:
        fasta_file.close()


def write_fasta_entry(outfile, name, seq, width):
    outfile.write('>%s\n' % name)
    if width is None:
        outfile.write(seq + '\n')
    else:
        for i in range(0, len(seq), width):
            outfile.write(seq[i:i + width] + '\n')


def write_fasta(fasta_entries, fpath=None, width=60):
    if fpath is None:
        fpath = sys.stdout
    with open(fpath, 'w') as outfile:
        for name, seq in fasta_entries:
            write_fasta_entry(outfile, name, seq, width)


def shuffle_fasta(fpath, out_fpath):
    def __shuffle(seq):
        return ''.join(random.sample(seq, len(seq)))

    shuffled_fasta = []
    for name, seq in read_fasta(fpath):
        shuffled_fasta.append((name, __shuffle(seq)))
    write_fasta(shuffled_fasta, out_fpath, width=None)


def is_nucleotide_fasta(fpath):
    MIN_LENGTH = 10
    MIN_ENTRIES = 5

    not_nucleotide_regexp = re.compile(r'[^acgtnACGTN]')
    for idx, (_, sample_seq) in enumerate(read_fasta(fpath)):
        if len(sample_seq) >= MIN_LENGTH:
            return not bool(not_nucleotide_regexp.search(sample_seq))
        if bool(not_nucleotide_regexp.search(sample_seq)):
            return False
        if idx + 1 >= MIN_ENTRIES:
            break
    return None


# for replacing fasta_formatter from FASTX package
def main():
    write_fasta(read_fasta(), width=None)


if __name__ == '__main__':
    main()