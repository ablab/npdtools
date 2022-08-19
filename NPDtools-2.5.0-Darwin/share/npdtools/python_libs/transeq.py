#!/usr/bin/python

# Substitution for 'transeq' utility from EMBOSS package
# common usage in RiPPquest project:  transeq input.fasta output.fasta.6_frame -frame 6

# based on http://pythonforbiologists.com/index.php/applied-python-for-biologists/applied-python-4/

import sys
import getopt
from fasta_utils import read_fasta, write_fasta_entry
from os.path import abspath

short_options = "o:f:w"
long_options = "output= frame= wide".split()
allowed_frames = ['1', '2', '3', '-1', '-2', '-3', '6']
dna_letters = ['A', 'C', 'G', 'T']
default_fasta_width = 60  # use --wide to write each sequence in one line

gencode = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*', 'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W'}


def comp(letter):
    return {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}[letter.upper()]


def rev_comp(seq):
    c = dict(zip('ATCGNatcgn', 'TAGCNtagcn'))
    return ''.join(c.get(nucleotide, '') for nucleotide in reversed(seq))


# a function to translate a single codon
def translate_codon(codon):
    assert len(codon) <= 3
    pure_codon = codon.upper().replace('N', '')
    if len(pure_codon) == 3:
        return gencode.get(codon.upper(), 'X')
    if len(pure_codon) < 2:
        return 'X'
    if len(codon) == 2:
        codon += 'N'
    candidate = gencode.get(codon.upper().replace('N', dna_letters[0]), 'X')
    for l in dna_letters[1:]:
        if candidate != gencode.get(codon.upper().replace('N', l), 'X'):
            return 'X'
    return candidate


# a function to split a sequence into codons
def split_into_codons(dna, frame):
    codons = []
    i = frame - 4
    for i in range(frame - 1, len(dna)-2, 3):
        codon = dna[i:i+3]
        codons.append(codon)
    if dna[i+3:]:
        codons.append(dna[i+3:])
    return codons


# a function to translate a dna sequence in a single frame
def translate_dna_single(dna, frame=1):
    codons = split_into_codons(dna, frame)
    amino_acids = ''
    for codon in codons:
        amino_acids = amino_acids + translate_codon(codon)
    return amino_acids


def translate_file(fpath, frame, out_fpath=None, wide=True):
    translated = []
    if frame == '6':
        frames = map(int, allowed_frames[:-1])
        need_revers = True
    else:
        frames = [int(frame)]
        if int(frame) < 0:
            need_revers = True
        else:
            need_revers = False

    write_to_file = True
    if out_fpath is None:
        out_fpath = sys.stdout
        write_to_file = False
    if wide:
        width = None
    else:
        width = default_fasta_width
    with open(out_fpath, 'w') as f:
        for (name, seq) in read_fasta(fpath):
            if need_revers:
                reversed_seq = rev_comp(seq)

            for cur_frame in frames:
                if cur_frame > 0:
                    translated_name = name + '_' + str(cur_frame)
                    translated_seq = translate_dna_single(seq, cur_frame)
                else:
                    translated_name = name + '_' + str(3 - cur_frame)
                    translated_seq = translate_dna_single(reversed_seq, (len(seq) + 4 + cur_frame) % 3 + 1)
                    # complicated frame recomputation formula for
                    # "The order of the frames follows the Staden convention:
                    # Frame -1 is the reverse-complement of the sequence having the same codon phase as frame 1.
                    # Frame -2 is the same phase as frame 2. Frame -3 is the same phase as frame 3."
                    # http://www.ebi.ac.uk/Tools/st/emboss_transeq/help/index.html
                if write_to_file:
                    write_fasta_entry(f, translated_name, translated_seq, width)
                else:
                    translated.append((translated_name, translated_seq))
    if write_to_file:
        return None
    return translated


def usage():
    sys.stderr.write('Usage: python ' + sys.argv[0] + ' input.fasta -o out.fasta --frame 6 [--wide]\n')
    # --wide for printing each FASTA entry in one line (without \n after each 60 symbols)


def error(msg):
    sys.stderr.write('ERROR: ' + msg + '\n')
    sys.exit(1)


def info(msg):
    sys.stdout.write('INFO: ' + msg + '\n')


def main():
    if len(sys.argv) < 2:
        usage()
        sys.exit(0)

    try:
        options, input_fpath = getopt.gnu_getopt(sys.argv[1:], short_options, long_options)
    except getopt.GetoptError:
        _, exc_value, _ = sys.exc_info()
        sys.stderr.write(str(exc_value) + '\n\n')
        usage()
        sys.exit(1)

    if len(input_fpath) != 1:
        error('Incorret number of input files (%d), should be exactly one!' % len(input_fpath))
    input_fpath = input_fpath[0]
    output_fpath = None
    frame = 6
    wide = False
    for opt, arg in options:
        if opt in ('-o', "--output"):
            output_fpath = abspath(arg)
        elif opt in ('-f', "--frame"):
            if arg in allowed_frames:
                frame = arg
            else:
                error('Incorrect value for FRAME (%s), should be one of %s' % (arg, ", ".join(allowed_frames)))
        elif opt in ('-w', "--wide"):
            wide = True

    translate_file(input_fpath, frame, output_fpath, wide)
    info('Translated nucleic acid sequences')


if __name__ == '__main__':
    main()
