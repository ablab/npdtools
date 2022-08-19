#!/usr/bin/env python

import sys
import os
from os.path import join, basename, isdir
import logging
from datetime import datetime
from site import addsitedir
from optparse import OptionParser

import npdtools_init
npdtools_init.init()

addsitedir(npdtools_init.python_modules_dir)

import common
import config
import log_utils
from log_utils import info, error, warning
from common import set_threads, check_file, clean, _clean_matches, parse_params_xml, \
    extract_unique_peptides, write_dict_to_tsv, get_command_line, print_stats
from fname_utils import update_fnames, generate_internal_fpath_mapping, _get_mapping_fpath
from file_utils import verify_file
from spectra_utils import get_spectra_fpaths
from config import prod_classes, platform_name
from spec_nets_utils import check_spec_network_files, parse_and_report_propagation
from r_utils import pipeline
from r_utils import ga_antismash
from r_utils import ga_fasta
from r_utils.pipeline import process_one_spectra_file, process_one_sequence_file, \
    save_results, save_params, get_sequence_fpaths
from r_utils.correspondence import get_spec_seq_correspondence, save_spec_seq_correspondence
from ordered_set import OrderedSet

ga_antismash.hmmsearch_fpath = join(npdtools_init.external_tools_dir, 'hmmer', platform_name,
                                      'hmmer', 'binaries', 'hmmsearch')
ga_antismash.hmmsearch_models_dir = join(npdtools_init.configs_dir, 'pfam_anal', 'ripp_hmm')
ga_fasta.rippquest_main_sh_fpath = join(npdtools_init.shell_scripts_dir, 'rippquest_genome_analysis.sh')
pipeline.rippquest_binary_fpath = join(npdtools_init.bin_dir, 'rippquest_ms')
pipeline.external_tools_dir = npdtools_init.external_tools_dir

log_utils.log = logging.getLogger(config.LOGGER_METAMINER_NAME)


class Config(object):
    def __init__(self, opts, params):
        # key dirs
        self.output_dirpath = opts.output_dirpath
        self.work_dirpath = opts.work_dirpath if opts.work_dirpath else join(self.output_dirpath, 'work')
        self.configs_dir = None  # it is set in "check_and_set"

        # other
        self.params = params  # specific MetaMiner parameters, e.g. mode, isotope
        self.pipeline = opts  # debug, reuse, etc
        self.log_fpath = join(self.work_dirpath, config.LOGGER_METAMINER_NAME + '.log')
        self.threads = opts.threads
        self.fdr = opts.fdr

        # special cases
        if 'fdr' in self.params:
            self.fdr = self.params['fdr']
        else:
            self.params['fdr'] = self.fdr
        self.params['debug'] = self.pipeline.debug

        # special pipeline modes
        setattr(self.pipeline, 'spec_nets_mode', opts.spec_network is not None)

    def check_and_set(self):
        if not isdir(self.output_dirpath):
            os.makedirs(self.output_dirpath)
        if not isdir(self.work_dirpath):
            os.makedirs(self.work_dirpath)
        self.configs_dir = common.copy_configs(npdtools_init.configs_dir, self.work_dirpath)
        pipeline.configs_dir = self.configs_dir


def get_used_sequence_fpaths(spec_seq_corr):
    used_sequence_fpaths = OrderedSet()
    for seq_set in spec_seq_corr.values():
        used_sequence_fpaths = used_sequence_fpaths | seq_set
    return list(used_sequence_fpaths)


def report_unused_spec_seq_files(spec_seq_corr, spectra_fpaths, sequence_fpaths):
    used_sequence_fpaths = get_used_sequence_fpaths(spec_seq_corr)
    # filtering spectra
    for spectra_fpath in spectra_fpaths:
        if len(spec_seq_corr[spectra_fpath]) == 0:
            warning("Skipping {spectra_fpath} because it does not have corresponding sequence files".format(**locals()))
            del spec_seq_corr[spectra_fpath]
    # filtering sequences
    if len(used_sequence_fpaths) != len(sequence_fpaths):
        for sequence_fpath in sequence_fpaths:
            if sequence_fpath not in used_sequence_fpaths:
                warning("Skipping {sequence_fpath} because it does not have corresponding spectra files".format(**locals()))


def parse_arguments():
    version = common.get_version()

    parser = OptionParser(description='MetaMiner: A Peptidogenomics Approach for the Discovery of Ribosomally '
                                      'Synthesized and Post-translationally Modified Peptides',
                          usage="%prog [options] -s sequence_files spectra_files", version=version)

    parser.add_option('-o', '--output', dest='output_dirpath', help='Output dir (required)')
    parser.add_option('-s', '--sequence', dest='sequence_paths', action='append',
                      help='Paths to sequence files (nucleotide/amino acid sequences in FASTA format, '
                           'or antiSMASH output in .gbk format, or BOA output in .txt format). '
                           'At least one file is required unless correspondence file (-C) with RefSeq IDs is specified.')

    parser.add_option('-t', '--threads', type='int', help='Threads number')

    # advanced settings TODO: put in a separate section if possible
    parser.add_option('-e', '--max-charge', help='Max possible charge [default: %default]', default='2')

    parser.add_option('-m', '--mode', help='Running mode (LL, HL, HH, custom) [default: %default]',
                      action='store', choices=['LL','HL','HH', 'custom'], default='custom')
    # parser.add_option('--ppm', help='Thresholds are given in ppm (defaults are 20 for "high" and '
    #                                 '500 for "low"). If not specified, thresholds are given in Da '
    #                                 '(defaults are 0.5 for "high" and 0.02 for "low")',
    #                   action="store_true", default=False)
    parser.add_option('--pm_thresh', help='Precursor ion threshold for "custom" Running mode (in Da) [default: %default]', default='0.02')
    parser.add_option('--product_ion_thresh', help='Product ion threshold for "custom" Running mode (in Da) [default: %default]', default='0.02')
    default_prod_class = 'lantibiotic'
    all_prod_class_options = prod_classes + ['all']
    parser.add_option('-c', '--class', dest='prod_class', default=default_prod_class, action='store',
                      choices=all_prod_class_options,
                      help='Product class. Valid choices are: %s [default: %s]' %
                           (', '.join(all_prod_class_options), default_prod_class))

    parser.add_option('--blind', help='Search with blind modifications (slow)', action="store_true", default=False)
    parser.add_option('--fdr', help='Compute FDR (doubles computation time)', action="store_true", default=False)

    parser.add_option('-C', '--correspondence', help='Path to sequence-spectra correspondence file '
                                                     '(tab-separated file with two columns indicating '
                                                     'spectra and sequence file basenames). Sequence column may include '
                                                     'RefSeq IDs prefixed with "#RefSeq:". In the latter case, '
                                                     'the corresponding references are automatically downloaded from NCBI. '
                                                     'If not specified the all-vs-all analysis is performed.')
    parser.add_option('--spec-network', help='Folder with Spectral Network output. Should contain params.xml, '
                                             'clusterinfo/<...>.clusterinfo and networkedges_selfloop/<...>.pairsinfo')
    parser.add_option('--keep-ga-files', help='Keep genome analysis files, e.g. '
                                              'for reusing them later with other spectra and --reuse option',
                      action="store_true", default=False)

    parser.add_option('--debug', help='Debug mode', action="store_true", default=False)
    parser.add_option('--reuse', help='Reuse existing files', action="store_true", default=False)

    if config.DEVELOPER_VERSION:
        parser.add_option('-w', '--work', dest='work_dirpath', help='Work dir')
        parser.add_option('-p', '--params', help='params.xml with all settings')
        parser.add_option('-P', '--pass-to-rippquest-ms',
                          help='Options string for passing to rippquest_ms binary "as is"')

    if len(sys.argv) < 2:
        parser.print_help()
        sys.exit(0)

    (opts, spectra) = parser.parse_args()
    if not opts.output_dirpath:
        parser.error("Output dir is required!")
    if opts.spec_network:
        check_spec_network_files(opts.spec_network)

    if not config.DEVELOPER_VERSION:
        setattr(opts, "work_dirpath", None)
        setattr(opts, "params", None)
        setattr(opts, "pass_to_rippquest_ms", "")

    return opts, spectra


def get_params_and_mapping(opts):
    if opts.params is not None:
        verify_file(opts.params, description="params file")
        params, file_mapping = parse_params_xml(opts.params)
    else:
        params = opts.__dict__  #TODO: take a subset of all options
        file_mapping = None
    return params, file_mapping


def metaminer(cfg, spec_seq_corr):
    _clean_matches(cfg)
    info('Command line: ' + get_command_line())
    info('Starting MetaMiner in %d threads' % cfg.threads)
    info('  Parameters are saved to ' + save_params(cfg))

    used_seq_fpaths = get_used_sequence_fpaths(spec_seq_corr)
    info('  Spectra-Sequence correspondence is saved to ' + save_spec_seq_correspondence(cfg, spec_seq_corr))

    if not len(used_seq_fpaths) or not len(spec_seq_corr):
        error('No input spectra or sequence files after performing Spectra-Sequence correspondence!')

    internal_fpath_mapping = generate_internal_fpath_mapping(cfg, used_seq_fpaths + list(spec_seq_corr.keys()),
                                                             num_seq=len(used_seq_fpaths))
    mapping_fpath = _get_mapping_fpath(cfg)
    write_dict_to_tsv(mapping_fpath, ['Input file path', 'Internal filename'], internal_fpath_mapping,
                      filter_func=basename)
    info('  Mapping of input file paths to internal filenames is saved to ' + mapping_fpath)

    num_PSMs = 0
    num_sig_PSMs = 0
    num_seq_failed = 0
    num_spec_failed = 0
    num_decoy = 0
    num_sig_decoy = 0
    # first stage: processing sequences
    info('*** Running Genome Analysis ***')
    chunk_size = cfg.threads
    for chunk_id, sequence_fpaths_chunk in enumerate([used_seq_fpaths[i:i+chunk_size] for i in range(0, len(used_seq_fpaths), chunk_size)]):
        n_jobs = min(len(sequence_fpaths_chunk), cfg.threads)
        if n_jobs > 1:
            parallel_args = [(cfg, sequence_fpath, cfg.params['prod_class'])
                             for sequence_fpath in sequence_fpaths_chunk]
            results = common.run_parallel(process_one_sequence_file, parallel_args, n_jobs=n_jobs)
            num_seq_failed += results.count(False)
        else:
            for sequence_fpath in sequence_fpaths_chunk:
                if not process_one_sequence_file(cfg, sequence_fpath, cfg.params['prod_class']):
                    num_seq_failed += 1
        info('-- Processed %d out of %d' % (min(len(used_seq_fpaths), (chunk_id + 1) * chunk_size), len(used_seq_fpaths)))

    if num_seq_failed == len(used_seq_fpaths):
        error('All sequence files failed to complete Genome Analysis! Skipping the rest.')

    # second stage: processing spectra
    info('*** Running Spectra Analysis ***')
    chunk_size = cfg.threads
    used_spec_fpaths = list(spec_seq_corr.keys())
    num_spectra_files = len(used_spec_fpaths)
    PSMs_for_spec_nets = []
    for chunk_id, spectra_fpaths_chunk in enumerate([used_spec_fpaths[i:i+chunk_size] for i in range(0, len(used_spec_fpaths), chunk_size)]):
        chunk_PSMs = []
        chunk_significant_PSMs = []
        n_jobs = min(len(spectra_fpaths_chunk), cfg.threads)
        if n_jobs > 1:
            parallel_args = [(cfg, spectra_fpath, spec_seq_corr[spectra_fpath])
                             for spectra_fpath in spectra_fpaths_chunk]
            results = common.run_parallel(process_one_spectra_file, parallel_args, n_jobs=n_jobs)
            PSMs_list, significant_PSMs_list, is_ok_list = \
                [x[0][0] for x in results], [x[0][1] for x in results], [x[1] for x in results]
            chunk_PSMs += [item for sublist in PSMs_list for item in sublist]
            chunk_significant_PSMs += [item for sublist in significant_PSMs_list for item in sublist]
            num_spec_failed += is_ok_list.count(False)
        else:
            for spectra_fpath in spectra_fpaths_chunk:
                (PSMs_list, significant_PSMs_list), is_ok = process_one_spectra_file(cfg,
                                                                spectra_fpath, spec_seq_corr[spectra_fpath])
                chunk_PSMs += PSMs_list
                chunk_significant_PSMs += significant_PSMs_list
                if not is_ok:
                    num_spec_failed += 1
        save_results(cfg, chunk_PSMs, chunk_significant_PSMs, first_time=(chunk_id == 0))
        if cfg.fdr:
            chunk_decoy, chunk_sig_decoy = save_results(cfg, chunk_PSMs, chunk_significant_PSMs,
                                                        first_time=(chunk_id == 0), decoy=True)
            num_decoy += chunk_decoy
            num_sig_decoy += chunk_sig_decoy
        num_PSMs += len(chunk_PSMs)
        num_sig_PSMs += len(chunk_significant_PSMs)
        if cfg.pipeline.spec_nets_mode:  # TODO: think which PSMs are needed for SpecNets (all/significant/etc)
            PSMs_for_spec_nets += chunk_significant_PSMs
        info('-- Processed %d out of %d' % (min(len(used_spec_fpaths), (chunk_id + 1) * chunk_size), num_spectra_files))
    stats_summary = {'num_spectra_files': num_spectra_files, 'num_sequence_files': len(used_seq_fpaths),
                     'num_failed_processing': num_spec_failed,
                     'num_PSMs': num_PSMs, 'num_sig_PSMs': num_sig_PSMs,
                     'num_decoy': num_decoy, 'num_sig_decoy': num_sig_decoy}
    return stats_summary, PSMs_for_spec_nets


def main():
    log_utils.log.setLevel(logging.DEBUG)
    log_utils.setup_log_handler(logging.StreamHandler(sys.stdout))

    opts, spectra = parse_arguments()
    params, file_mapping = get_params_and_mapping(opts)
    cfg = Config(opts, params)
    cfg.check_and_set()
    set_threads(cfg)

    log_handler = logging.FileHandler(cfg.log_fpath, mode='w')
    log_utils.setup_log_handler(log_handler)
    config.global_config = cfg

    spectra_fpaths = get_spectra_fpaths(spectra)
    if not spectra_fpaths:
        error('No input spectra files were found!')

    sequence_fpaths = get_sequence_fpaths(opts.sequence_paths)
    if not sequence_fpaths:
        if opts.correspondence and os.path.isfile(opts.correspondence):
            warning('No input sequence files were found! Will check whether there are downloable RefSeq IDs in the correpondence file')
        else:
            error('No input sequence files were found! Please specify them with -s option!')

    spec_seq_corr = get_spec_seq_correspondence(cfg, opts.correspondence, spectra_fpaths, sequence_fpaths, file_mapping)
    report_unused_spec_seq_files(spec_seq_corr, spectra_fpaths, sequence_fpaths)

    start_time = datetime.now()
    stats_summary, PSMs_for_spec_nets = metaminer(cfg, spec_seq_corr)
    if cfg.pipeline.spec_nets_mode and len(PSMs_for_spec_nets):
        global_output_dirpath = cfg.output_dirpath
        cfg.output_dirpath = join(cfg.output_dirpath, 'spec_nets')  # TODO: move "spec_nets" somewhere as a constant
        if not isdir(cfg.output_dirpath):
            os.makedirs(cfg.output_dirpath)
        parse_and_report_propagation(cfg, spectra_fpaths, opts.spec_network, PSMs_for_spec_nets)
        cfg.output_dirpath = global_output_dirpath
    finish_time = datetime.now()

    if file_mapping is not None:
        update_fnames(cfg, file_mapping)
    extract_unique_peptides(cfg, append_fdr=cfg.fdr, peptide_column='FragmentSeq')
    print_stats(cfg, stats_summary, finish_time - start_time, config.rippquest_p_value_threshold)
    log_utils.log.removeHandler(log_handler)


if __name__ == "__main__":
    try:
        main()
    except Exception:
        _, exc_value, _ = sys.exc_info()
        log_utils.exception(exc_value)
        error('Exception caught!')
    finally:
        if config.global_config is not None:
            clean(config.global_config)

