#!/usr/bin/env python

#####################################################################################
# Copyright (c) 2015-2017 Saint Petersburg State University, St. Petersburg, Russia
# Copyright (c) 2015-2017 University of California San Diego, La Jolla, CA, USA
# All Rights Reserved
# See file LICENSE for details.
#####################################################################################

import sys
import os
import traceback
import logging
from os.path import join, abspath, isdir, isfile
from datetime import datetime
from site import addsitedir
from optparse import OptionParser

import npdtools_init
npdtools_init.init()

addsitedir(npdtools_init.python_modules_dir)

from spectra_utils import get_spectra_fpaths
import common
import config
import log_utils
from log_utils import info, error
from common import set_threads, clean, _clean_matches, parse_params_xml, extract_unique_peptides, get_command_line, \
    print_stats
from fname_utils import update_fnames
from file_utils import verify_file, verify_dir
from d_utils.pipeline import append_smiles, save_params, process_single_spectra, save_full_results

from d_utils import pipeline
pipeline.dereplicate_fpath = join(npdtools_init.bin_dir, 'dereplicate')
pipeline.bin_dir = npdtools_init.bin_dir
pipeline.external_tools_dir = npdtools_init.external_tools_dir
import db_preprocessing

log_utils.log = logging.getLogger(config.LOGGER_DEREPLICATOR_NAME)


class Config(object):
    def __init__(self, opts, params):
        # key dirs
        self.output_dirpath = opts.output_dirpath
        self.work_dirpath = opts.work_dirpath if opts.work_dirpath else join(self.output_dirpath, 'work')
        self.configs_dir = None  # it is set in "check_and_set"

        # db stuff
        self.db_path = opts.db_path if opts.db_path else npdtools_init.configs_dir
        self.library_info = opts.library_info if opts.library_info else \
            join(self.db_path, 'library.info')
        self.library_smiles = opts.smiles if opts.smiles else \
            (self.library_info.replace('.info', '.smiles') if '.info' in self.library_info else '')
        self.db_entries = None

        # other
        self.params = params  # specific Dereplicator parameters, e.g. mode, isotope
        self.pipeline = opts  # debug, reuse, etc
        self.log_fpath = join(self.work_dirpath, config.LOGGER_DEREPLICATOR_NAME + '.log')
        self.threads = opts.threads
        self.fdr = opts.fdr

        # special cases
        if 'fdr' in self.params:
            self.fdr = self.params['fdr']
        else:
            self.params['fdr'] = self.fdr
        self.params['debug'] = self.pipeline.debug
        if 'pass_to_binary' in self.params and self.params['pass_to_binary']:
            setattr(self.pipeline, 'pass_to_dereplicate', self.params['pass_to_binary'])
        if 'min_score' in self.params and self.params["min_score"]:
            try:
                config.der_plus_score_threshold = float(str(self.params["min_score"]))
            except ValueError:
                pass

        # special settings for Dereplicator+
        if 'fragmentation_mode' not in self.params:
            self.params['fragmentation_mode'] = self.pipeline.fragmentation_mode
        setattr(self.pipeline, 'vis_type', 'none')
        #self.params['mode'] = 'custom'
        self.params['dereplicator+'] = True
        # absolute path of self.fragmentation_file depends on self.configs_dir and is set in "check_and_set"
        self.fragmentation_file = join('configs', 'common', 'fragment_list_%s.txt' % self.params['fragmentation_mode'])
        config.USE_SCORE_ONLY = True
        self.params['min_score'] = config.DER_PLUS_ALL_MIN_SCORE_RATE * config.der_plus_score_threshold
        # FT stuff
        self.preprocessed_ft = opts.preprocessed_ft if opts.preprocessed_ft else \
            (self.library_info + '.' + self.params['fragmentation_mode'] + '.ft' + db_preprocessing.DB_SUFFIX)
        self.preprocessed_ftd = opts.preprocessed_ftd if opts.preprocessed_ftd else \
            (self.library_info + '.' + self.params['fragmentation_mode'] + '.ftd' + db_preprocessing.DB_SUFFIX)


    def check_and_set(self):
        # verifying mandatory things
        verify_file(self.library_info, description='library info')
        verify_dir(self.db_path, description='database')
        # creating dirs
        if not isdir(self.output_dirpath):
            os.makedirs(self.output_dirpath)
        if not isdir(self.work_dirpath):
            os.makedirs(self.work_dirpath)
        self.configs_dir = common.copy_configs(npdtools_init.configs_dir, self.work_dirpath)
        self.fragmentation_file = join(self.configs_dir, self.fragmentation_file)
        verify_file(self.fragmentation_file, description='fragmentation file')
        #
        self.__read_db_entries()

    def __read_db_entries(self):
        self.db_entries = []
        with open(self.library_info) as f:
            for line in f:
                entry = line.split()
                if len(entry) > 2:
                    fpath = entry[0]
                    self.db_entries.append(abspath(join(self.db_path, fpath)))


def parse_arguments():
    version = common.get_version()

    parser = OptionParser(description='Dereplicator+: identification of metabolites '
                                      'through database search of mass spectra',
                          usage="%prog [options] spectra_files", version=version)
    parser.add_option('-o', '--output', dest='output_dirpath', help='Output dir (required)')
    parser.add_option('-d', '--db-path', help='Path to a metabolite database (required)')
    parser.add_option('-l', '--library-info', help='Path to the database description file (default: library.info inside the database dir)')
    parser.add_option('-s', '--smiles', help='Path to the database list of SMILES (default: library.smiles inside the database dir)')

    parser.add_option('-t', '--threads', type='int', help='Threads number [default: %default]', default=1)

    # Dereplicator+ specific options
    parser.add_option('--preprocess', help='Perform database preprocessing before dereplication. '
                                           'The preprocessing takes non-trivial time but may give a substantial speedup '
                                           'if you plan to dereplicate many spectra. The preprocessed database files may be reused '
                                           'in consecutive runs against the same database and using the same fragmentation mode.',
                      action="store_true", default=False)
    parser.add_option('--preprocessed_ft', help='Path to the preprocessed fragmentation trees '
                                                '(default: <library_info>.<fragmentation_mode>.ft.bin); '
                                                'you can preprocess your database using --preprocess')
    parser.add_option('--preprocessed_ftd', help='Path to the preprocessed fragmentation tree decoys '
                                                 '[used for FDR computation only, so --fdr should be specified] '
                                                '(default: <library_info>.<fragmentation_mode>.ftd.bin); '
                                                'you can preprocess your database using --preprocess')
    parser.add_option('--fragmentation_mode', help='Mode of fragmentation (e.g. "general_6_1_6" or "general_3_1_3") '
                                                   '[default: %default]', default='general_6_1_6')

    # advanced settings TODO: put in a separate section if possible
    parser.add_option('-e', '--max-charge', help='Max considered charge [default: %default]', default='2')

    parser.add_option('-m', '--mode', help='Running mode (LL, HL, HH, custom) [default: %default]',
                      action='store', choices=['LL', 'HL', 'HH', 'custom'], default='custom')
    parser.add_option('--pm_thresh', help='Precursor ion threshold in Da [default: %default]', default='0.005')
    parser.add_option('--product_ion_thresh', help='Product ion threshold in Da [default: %default]', default='0.01')
    parser.add_option('--ppm', help='Thresholds are given in ppm', action="store_true", default=False)

    parser.add_option('--fdr', help='Compute FDR (time consuming)', action="store_true", default=False)
    parser.add_option('--fdr-limit', help='Maximum allowed FDR in percents for significant matches (in 0.0-100.0 range) '
                                          '[default: %default]', default=1)
    parser.add_option('--min-score', help='Minimum score for significant matches (integer number) '
                                          '[default: %default]', default=config.DEFAULT_DER_PLUS_SIG_MIN_SCORE)

    parser.add_option('--debug', help='Debug mode', action="store_true", default=False)
    parser.add_option('--reuse', help='Reuse existing files', action="store_true", default=False)

    if config.DEVELOPER_VERSION:
        parser.add_option('-w', '--work', dest='work_dirpath',
                          help='Specific dir for storing intermediate results and logs [default: <out_dir>/work]')
        parser.add_option('-p', '--params', help='params.xml with all settings')
        parser.add_option('-P', '--pass-to-dereplicate',
                          help='Options string for passing to dereplicate binary "as is"',
                          default='')

    if len(sys.argv) < 2:
        parser.print_help()
        sys.exit(0)

    (opts, spectra) = parser.parse_args()
    if not opts.output_dirpath:
        parser.error("Output dir is required!")
    if not opts.db_path:
        parser.error("Database is required!")

    for fpath in [opts.library_info, opts.smiles, opts.preprocessed_ft, opts.preprocessed_ftd]:
        if fpath and not isfile(fpath):
            parser.error("Filepath %s was specified but such file does not exist!" % fpath)

    if not config.DEVELOPER_VERSION:
        setattr(opts, "work_dirpath", None)
        setattr(opts, "params", None)
        setattr(opts, "pass_to_dereplicate", "")

    return opts, spectra


def get_params_and_mapping(opts):
    if opts.params is not None:
        verify_file(opts.params, description="params file")
        params, file_mapping = parse_params_xml(opts.params)
    else:
        params = opts.__dict__  #TODO: take a subset of all options
        file_mapping = None
    return params, file_mapping


def preprocess_db_if_needed(cfg):
    info('Preprocessing database')
    preproc_dir = join(cfg.output_dirpath, 'db_preproc')
    from db_preprocessing import preprocess_db
    is_needed = False
    if not cfg.preprocessed_ft or not isfile(cfg.preprocessed_ft):
        cfg.preprocessed_ft = preprocess_db('ft', preproc_dir, cfg.db_path,
                                            cfg.library_info, cfg.configs_dir, cfg.pipeline.__dict__)
        is_needed = True
    if cfg.fdr and (not cfg.preprocessed_ftd or not isfile(cfg.preprocessed_ftd)):
        cfg.preprocessed_ftd = preprocess_db('ftd', preproc_dir, cfg.db_path,
                                             cfg.library_info, cfg.configs_dir, cfg.pipeline.__dict__)
        is_needed = True
    if is_needed:
        info('Preprocessing finished and its results can be reused later (see %s)' % preproc_dir)
    else:
        info('  Preprocessing was not performed (not needed)')


def dereplicate(cfg, spectra_fpaths):
    _clean_matches(cfg)
    info('Command line: ' + get_command_line())
    info('Starting Dereplicator+ in %d thread' % cfg.threads + ('s' if cfg.threads > 1 else ''))
    info('  Parameters are saved to ' + save_params(cfg, dereplicator_plus_mode=True))
    info('  Total number of spectra: %d ' % len(spectra_fpaths))

    if cfg.pipeline.preprocess:
        preprocess_db_if_needed(cfg)

    all_PSMs = []
    n_jobs = min(len(spectra_fpaths), cfg.threads)
    if n_jobs > 1:
        parallel_args = [(cfg, spectra_fpath, idx)
            for idx, spectra_fpath in enumerate(spectra_fpaths)]
        results = common.run_parallel(process_single_spectra, parallel_args, n_jobs=n_jobs)
        prepared_fpaths, all_stats = \
            [x[0] for x in results], [x[2] for x in results]
    else:
        prepared_fpaths = []
        all_stats = []
        for idx, spectra_fpath in enumerate(spectra_fpaths):
            prepared_fpath, _, stats = process_single_spectra(cfg, spectra_fpath, idx)
            prepared_fpaths.append(prepared_fpath)
            all_stats.append(stats)

    # stats
    num_failed_conversion = prepared_fpaths.count(None)
    num_failed_processing = all_stats.count(None) - num_failed_conversion
    num_spectra_files = len(spectra_fpaths)
    num_scans_processed = 0
    num_PSMs, num_sig_PSMs = save_full_results(cfg, prepared_fpaths, headers_to_exclude=['VisualizationID', 'P-Value'])
    if cfg.fdr:
        num_decoy, num_sig_decoy = save_full_results(cfg, prepared_fpaths, decoy=True,
                                                     headers_to_exclude=['VisualizationID', 'P-Value'])
    else:
        num_decoy, num_sig_decoy = 0, 0

    for stats in all_stats:
        if stats is None:
            continue
        num_scans_processed += stats[0]

    stats_summary = {'num_spectra_files': num_spectra_files, 'num_spectra_scans': num_scans_processed,
                     'num_failed_conversion': num_failed_conversion, 'num_failed_processing': num_failed_processing,
                     'num_PSMs': num_PSMs + num_decoy, 'num_sig_PSMs': num_sig_PSMs + num_sig_decoy,
                     'num_decoy': num_decoy, 'num_sig_decoy': num_sig_decoy}
    return stats_summary, all_PSMs


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
        error('No input spectra were found!')

    start_time = datetime.now()
    stats_summary, all_PSMs = dereplicate(cfg, spectra_fpaths)
    finish_time = datetime.now()

    append_smiles(cfg)
    if file_mapping is not None:
        update_fnames(cfg, file_mapping)
    extract_unique_peptides(cfg, append_fdr=cfg.fdr,
                            fdr_limit=cfg.params['fdr_limit'] if 'fdr_limit' in cfg.params else None)
    print_stats(cfg, stats_summary, finish_time - start_time, p_value_thresh=None, score_thresh=config.der_plus_score_threshold)
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

