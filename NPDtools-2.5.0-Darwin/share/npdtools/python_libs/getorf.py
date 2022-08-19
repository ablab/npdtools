#!/usr/bin/python

# Substitution for 'getorf' utility from EMBOSS package
# common usage in RiPPquest project:  getorf input.fasta -o output.fasta.orf.format -w

# usage of EMBOSS getorf:
# getorf input.fasta input.fasta.orf -find 0  # between stop codons
# cat input.fasta.orf | $BASEDIR/fasta_utils.py | tr -d 'X' | sed 's/] (REVERSE SENSE)/_1/g' | sed 's/]/_0/g' | sed 's/ - /_/g' | sed 's/ \[/_/g' > input.fasta.orf.format


import sys
import getopt
from fasta_utils import read_fasta, write_fasta
from transeq import error, info, translate_file, translate_dna_single, rev_comp
from os.path import abspath

short_options = "o:w"
long_options = "output= wide".split()
MIN_ORF_SIZE = 30  # in nucleotides
MAX_ORF_SIZE = 1000000  # in nucleotides

class ORF(object):
    def __init__(self, basename, seq, start, reverse=False):
        self.name = basename
        self.seq = seq
        self.start = start
        self.reverse = reverse
        if self.reverse:
            self.end = start - 3 * len(self.seq) + 1
        else:
            self.end = start + 3 * len(self.seq) - 1

    def get_fasta_entry(self, idx):
        return self.name + '_%d_%d_%d_%d' % (idx, self.start, self.end, self.reverse), self.seq

    def is_good(self):
        return MIN_ORF_SIZE // 3 <= len(self.seq) <= MAX_ORF_SIZE // 3


def usage():
    sys.stderr.write('Usage: python ' + sys.argv[0] + ' input.fasta -o out.fasta [--wide]\n')
    # --wide for printing each FASTA entry in one line (without \n after each 60 symbols)


def get_orfs(name, translated_seq, start_pos, reverse=False):
    orfs = []
    position = start_pos
    for seq in translated_seq.split('*'):
        seq = seq.rstrip('X')
        first_x = True
        for split_seq in seq.split('X'):
            if split_seq:
                orfs.append(ORF(name, split_seq, position, reverse))
                position += 3 * len(split_seq) * (-1 if reverse else 1)
                first_x = True
            else:
                position += 3 * (-1 if reverse else 1)
                if first_x:
                    position += 3 * (-1 if reverse else 1)  # for 'X' hidden by split()
                    first_x = False
        position += 3 * (-1 if reverse else 1)  # for '*' hidden by 'split()
    return orfs


def filter_and_convert_orfs(orfs, idx=1):
    filtered = []
    for orf in orfs:
        if orf.is_good():
            filtered.append(orf.get_fasta_entry(idx))
            idx += 1
    return filtered


def main():
    if len(sys.argv) < 2:
        usage()
        sys.exit(0)

    try:
        options, input_fpath = getopt.gnu_getopt(sys.argv[1:], short_options, long_options)
    except getopt.GetoptError:
        _, exc_value, _ = sys.exc_info()
        sys.stderr.write(str(exc_value) + '\n')
        sys.stderr.write('\n')
        usage()
        sys.exit(1)

    if len(input_fpath) != 1:
        error('Incorrect number of input files (%d), should be exactly one!' % len(input_fpath))
    input_fpath = input_fpath[0]
    output_fpath = None
    wide = False
    for opt, arg in options:
        if opt in ('-o', "--output"):
            output_fpath = abspath(arg)
        elif opt in ('-w', "--wide"):
            wide = True

    all_orfs_fasta_entries = []
    for (name, seq) in read_fasta(input_fpath):
        # forward strain
        orfs = []
        for frame in range(1, 4):
            translated_seq = translate_dna_single(seq, frame)
            orfs += get_orfs(name, translated_seq, frame)
        orfs.sort(key=lambda x: x.end)
        all_orfs_fasta_entries += filter_and_convert_orfs(orfs)

        # reverse strain
        rev_seq = rev_comp(seq)
        rev_orfs = []
        for frame in range(1, 4):
            translated_seq = translate_dna_single(rev_seq, frame)
            rev_orfs += get_orfs(name, translated_seq, len(rev_seq) - frame + 1, reverse=True)
        rev_orfs.sort(key=lambda x: x.end, reverse=True)
        all_orfs_fasta_entries += filter_and_convert_orfs(rev_orfs, idx=len(all_orfs_fasta_entries) + 1)

    if wide:
        write_fasta(all_orfs_fasta_entries, output_fpath, width=None)
    else:
        write_fasta(all_orfs_fasta_entries, output_fpath)


if __name__ == '__main__':
    main()