import os
from r_utils.pipeline import get_prod_class_mod_fpath, get_prod_classes_to_process, get_genome_analysis_fpath
from file_utils import remove_if_exists, is_complete_log
import config
from common import info, warning, sys_call
from os.path import isfile, basename
from fasta_utils import write_fasta, read_fasta

rippquest_main_sh_fpath = None


def process_regular_fasta_file(cfg, sequence_fpath, prod_class):
    # stage 1: genome analysis
    # sh rippquest_genome_analysis.sh $seq_path $compound_type $mod_file
    classes_to_process = get_prod_classes_to_process(prod_class)

    ga_fpaths = []
    for prod_class in classes_to_process:
        prod_class_mod_fpath = get_prod_class_mod_fpath(prod_class)
        ga_fpath = get_genome_analysis_fpath(sequence_fpath, prod_class)
        log_fpath = sequence_fpath + '_' + prod_class + '_' + config.LOGGER_METAMINER_NAME + '.log'
        if cfg.pipeline.reuse and is_complete_log(log_fpath):
            info('    Reusing existing Genome Analysis results: ' + ga_fpath)
            if not isfile(ga_fpath):
                warning('    Genome Analysis failed on {sequence_fpath} ({prod_class}) -- '
                        '(based on reused results)!'.format(**locals()))
                continue
        else:
            remove_if_exists(ga_fpath)
            cur_output = sequence_fpath + '_' + prod_class + '_' + config.LOGGER_METAMINER_NAME + '.log'
            command = 'sh ' + rippquest_main_sh_fpath + ' "{sequence_fpath}" {prod_class} ' \
                                                        '"{prod_class_mod_fpath}" ' \
                                                        '> "{cur_output}"'.format(**locals())
            if not sys_call(command, ga_fpath, indent='    '):
                warning('    Genome Analysis failed on {sequence_fpath} ({prod_class}), '
                        'see the log for details ({cur_output})!'.format(**locals()))
                remove_if_exists(ga_fpath)
                continue
        ga_fpaths.append(ga_fpath)
    return ga_fpaths


def process_aminoacid_fasta_file(cfg, sequence_fpath, prod_class, actual_in_fpath=None):
    classes_to_process = get_prod_classes_to_process(prod_class)
    ga_fpaths = []

    prod_class = classes_to_process[0]
    ga_fpath = get_genome_analysis_fpath(sequence_fpath, prod_class)
    if cfg.pipeline.reuse and isfile(ga_fpath):
        info('    Reusing existing Genome Analysis results for ' + sequence_fpath)
        ga_fpaths.append(ga_fpath)
    else:
        remove_if_exists(ga_fpath)
        fasta_entries = []
        in_sequence_fpath = actual_in_fpath if actual_in_fpath is not None else sequence_fpath
        for name, seq in read_fasta(in_sequence_fpath):
            fasta_entries.append((name.replace('_', '.') + '_1_allenzyme', seq))
        if not fasta_entries:
            warning('    Genome Analysis failed on {sequence_fpath}: no sequences in this file!'.format(**locals()))
            return []
        write_fasta(fasta_entries, ga_fpath, width=None)
        ga_fpaths.append(ga_fpath)

    # lightweight process, do not need to reuse
    if len(classes_to_process) > 1:
        for prod_class in classes_to_process[1:]:
            symlink_ga_fpath = get_genome_analysis_fpath(sequence_fpath, prod_class)
            remove_if_exists(symlink_ga_fpath)
            os.symlink(basename(ga_fpath), symlink_ga_fpath)
            ga_fpaths.append(symlink_ga_fpath)
    return ga_fpaths