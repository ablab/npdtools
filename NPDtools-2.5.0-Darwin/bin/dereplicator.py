#!/usr/bin/env python

#####################################################################################
# Copyright (c) 2015-2017 Saint Petersburg State University, St. Petersburg, Russia
# Copyright (c) 2015-2017 University of California San Diego, La Jolla, CA, USA
# All Rights Reserved
# See file LICENSE for details.
#####################################################################################

import sys
import os
import logging
from os.path import join, abspath, isdir
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
from spec_nets_utils import check_spec_network_files, parse_and_propagate
from d_utils.pipeline import append_smiles, save_params, process_single_spectra, save_full_results
from d_utils.aux import extract_best

from d_utils import pipeline
pipeline.dereplicate_fpath = join(npdtools_init.bin_dir, 'dereplicate')
pipeline.bin_dir = npdtools_init.bin_dir
pipeline.external_tools_dir = npdtools_init.external_tools_dir

log_utils.log = logging.getLogger(config.LOGGER_DEREPLICATOR_NAME)


class Config(object):
    def __init__(self, opts, params):
        # key dirs
        self.output_dirpath = opts.output_dirpath
        self.work_dirpath = opts.work_dirpath if opts.work_dirpath else join(self.output_dirpath, 'work')
        self.vis_dirpath = opts.vis_dirpath if hasattr(opts, 'vis_dirpath') and opts.vis_dirpath else join(self.output_dirpath, 'visualizations')
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

        # special pipeline modes
        setattr(self.pipeline, 'spec_nets_mode', opts.spec_network is not None)
        setattr(self.pipeline, 'varquest_mode', 'varquest' in self.params and self.params['varquest'])
        setattr(self.pipeline, 'variable_mode', self.pipeline.varquest_mode or self.pipeline.spec_nets_mode)

        # special cases
        if 'fdr' in self.params:
            self.fdr = self.params['fdr']
        else:
            self.params['fdr'] = self.fdr
        self.params['debug'] = self.pipeline.debug
        if 'mode' not in self.params:
            self.params['mode'] = self.pipeline.mode
        if 'accurate_p_value' in self.params and self.params["accurate_p_value"]:
            config.dereplicator_p_value_threshold = config.ACCURATE_P_VALUE_THRESHOLD
        if self.pipeline.varquest_mode:
            self.params['theor_spectra'] = join(self.work_dirpath, "ts.dat")
        if 'pass_to_binary' in self.params and self.params['pass_to_binary']:
            setattr(self.pipeline, 'pass_to_dereplicate', self.params['pass_to_binary'])
        if 'p_value_limit' in self.params and self.params["p_value_limit"]:
            try:
                config.dereplicator_p_value_threshold = float(str(self.params["p_value_limit"]))
            except ValueError:
                pass

    def check_and_set(self):
        # verifying mandatory things
        verify_file(self.library_info, description='library info')
        verify_dir(self.db_path, description='database')
        # creating dirs
        if not isdir(self.output_dirpath):
            os.makedirs(self.output_dirpath)
        if not isdir(self.work_dirpath):
            os.makedirs(self.work_dirpath)
        if self.pipeline.vis_type != 'none' and not isdir(self.vis_dirpath):
            os.makedirs(self.vis_dirpath)
        self.configs_dir = common.copy_configs(npdtools_init.configs_dir, self.work_dirpath)
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
    if os.path.basename(sys.argv[0]) == "varquest.py":
        varquest_mode = True
    else:
        varquest_mode = False
    version = common.get_version()

    if varquest_mode:
        parser = OptionParser(description='VarQuest: modification-tolerant identification of novel variants of '
                                          'peptidic antibiotics and other natural products',
                              usage="%prog [options] spectra_files", version=version)
    else:
        parser = OptionParser(description='Dereplicator: in silico identification of peptidic natural products '
                                          'through database search of mass spectra',
                              usage="%prog [options] spectra_files", version=version)
    parser.add_option('-o', '--output', dest='output_dirpath', help='Output dir (required)')
    parser.add_option('-d', '--db-path', help='Path to a PNP database (required)')
    parser.add_option('-l', '--library-info', help='Path to the database description file (default: library.info inside the database dir)')
    parser.add_option('-s', '--smiles', help='Path to the database list of SMILES (default: library.smiles inside the database dir)')

    parser.add_option('-t', '--threads', type='int', help='Threads number')

    # advanced settings TODO: put in a separate section if possible
    parser.add_option('-e', '--max-charge', help='Max possible charge [default: %default]', default='2')

    parser.add_option('-m', '--mode', help='Running mode (LL, HL, HH, custom) [default: %default]',
                      action='store', choices=['LL', 'HL', 'HH', 'custom'], default='custom')

    parser.add_option('--ppm', help='Thresholds are given in ppm (defaults are 20 for "high" and '
                                    '500 for "low"). If not specified, thresholds are given in Da '
                                    '(defaults are 0.5 for "high" and 0.02 for "low")',
                      action="store_true", default=False)
    parser.add_option('--pm_thresh', help='Precursor ion threshold for "custom" Running mode [default: %default]', default='0.02')
    parser.add_option('--product_ion_thresh', help='Product ion threshold for "custom" Running mode [default: %default]', default='0.02')

    if not varquest_mode:
        parser.add_option('-i', '--isotope', help='Maximum isotopic shift accepted [default: %default]',
                          action='store', choices=['0', '1', '2'], default='0')

    parser.add_option('--fdr', help='Compute FDR (time consuming)', action="store_true", default=False)
    parser.add_option('--fdr-limit', help='Maximum allowed FDR in percents for significant matches (in 0.0-100.0 range) '
                                          '[default: %default]', default=1)
    parser.add_option('--p-value-limit', help='Minimum allowed P-value for significant matches (in 0.0-1.0 range) '
                                              '[default: %default]', default=config.dereplicator_p_value_threshold)
    parser.add_option('--nps', help='Use NPS scoring and significance estimation model (slower than default method)',
                      action="store_true", default=False)

    # Variable identification settings
    if not varquest_mode:
        parser.add_option('--varquest', help='Use VarQuest algorithm for variable PNP identification',
                          action="store_true", default=False)
        parser.add_option('--max-mod',
                          help='Maximum allowed modification size (in Da; VarQuest only) [default: %default]',
                          default='300')
    else:
        parser.add_option('--max-mod', help='Maximum allowed modification size (in Da) [default: %default]',
                          default='300')

    parser.add_option('--debug', help='Debug mode (verbose, keep intermediate files)', action="store_true",
                      default=False)
    parser.add_option('--reuse', help='Reuse existing files', action="store_true", default=False)

    # options mostly needed for GNPS release
    if config.DEVELOPER_VERSION:
        parser.add_option('-w', '--work', dest='work_dirpath',
                          help='Specific dir for storing intermediate results and logs [default: <out_dir>/work]')
        parser.add_option('-z', '--vis-dir', dest='vis_dirpath',
                          help='Specific dir for storing visualization results [default: <out_dir>/visualizations]')
        parser.add_option('-p', '--params', help='params.xml with all settings')
        parser.add_option('-V', '--vis-type', dest='vis_type', help='Type of visualization: '
                                                                    'none -- disable; html -- plain html; '
                                                                    'portable_html -- html with all aux scripts encapsulated; '
                                                                    'json -- json only (suitable for GNPS); txt -- plain text (structure only); '
                                                                    'extended_txt -- plain text (structure, matched peaks, etc). '
                                                                    '[default: %default]',
                          choices=['none', 'html', 'portable_html', 'json', 'txt', 'extended_txt'], default='none')
        parser.add_option('--extract-best', help='Extract spectra with best hits into MGF file', action="store_true",
                          default=False)
        parser.add_option('--adduct_Na', help='Consider +Na adduct', action="store_true", default=False)
        parser.add_option('--adduct_K', help='Consider +K adduct', action="store_true", default=False)
        parser.add_option('--accurate_p_value',
                          help='Use accurate algorithm for P-value calculation (MS-DPR, slower than default method)',
                          action="store_true", default=False)
        parser.add_option('-P', '--pass-to-dereplicate', help='Options string for passing to dereplicate binary "as is"')
        parser.add_option('--spec-network', help='Folder with Spectral Network output. Should contain params.xml, '
                                                 'clusterinfo/<...>.clusterinfo and networkedges_selfloop/<...>.pairsinfo')

    if len(sys.argv) < 2:
        parser.print_help()
        sys.exit(0)

    (opts, spectra) = parser.parse_args()
    if not opts.output_dirpath:
        parser.error("Output dir is required!")
    if not opts.db_path:
        parser.error("PNP database is required!")

    if varquest_mode:
        setattr(opts, "varquest", True)
        setattr(opts, "spec_network", None)
        setattr(opts, "isotope", 0)

    ### Setting currently disabled attributes in Release mode
    if not config.DEVELOPER_VERSION:
        setattr(opts, "work_dirpath", None)
        setattr(opts, "vis_dirpath", None)
        setattr(opts, "params", None)
        setattr(opts, "vis_type", "none")
        setattr(opts, "accurate_p_value", True)
        setattr(opts, "adduct_Na", False)
        setattr(opts, "adduct_K", False)
        setattr(opts, "pass_to_dereplicate", "--min_score_for_pvalue 7 --min_score 5 --min_shared_peaks 5 --min_num_comp 0")
        setattr(opts, "extract_best", False)
        setattr(opts, "spec_network", None)

    if opts.spec_network and opts.varquest:
        parser.error("You cannot specify SpecNets and VarQuest simultaneously!")
    if opts.spec_network:
        check_spec_network_files(opts.spec_network)
    return opts, spectra


def get_params_and_mapping(opts):
    if opts.params is not None:
        verify_file(opts.params, description="params file")
        params, file_mapping = parse_params_xml(opts.params)
    else:
        params = opts.__dict__  #TODO: take a subset of all options
        file_mapping = None
    return params, file_mapping


def dereplicate(cfg, spectra_fpaths, return_all_PSMs=False, msg='', save_in_files=True):
    _clean_matches(cfg)
    info('Command line: ' + get_command_line())
    if msg:
        info(msg)
    info('Starting Dereplicator in %d thread' % cfg.threads + ('s' if cfg.threads > 1 else ''))
    info('  Parameters are saved to ' + save_params(cfg))
    info('  Total number of spectra: %d ' % len(spectra_fpaths))

    all_PSMs = []
    n_jobs = min(len(spectra_fpaths), cfg.threads)
    if n_jobs > 1:
        parallel_args = [(cfg, spectra_fpath, idx, return_all_PSMs)
                         for idx, spectra_fpath in enumerate(spectra_fpaths)]
        results = common.run_parallel(process_single_spectra, parallel_args, n_jobs=n_jobs)
        prepared_fpaths, PSMs_list, all_stats = \
            [x[0] for x in results], [x[1] for x in results], [x[2] for x in results]
        if return_all_PSMs:
            all_PSMs += [item for sublist in PSMs_list for item in sublist]
    else:
        prepared_fpaths = []
        all_stats = []
        for idx, spectra_fpath in enumerate(spectra_fpaths):
            prepared_fpath, PSMs, stats = process_single_spectra(cfg, spectra_fpath, idx, return_all_PSMs)
            prepared_fpaths.append(prepared_fpath)
            all_PSMs += PSMs
            all_stats.append(stats)

    # stats
    num_failed_conversion = prepared_fpaths.count(None)
    num_failed_processing = all_stats.count(None) - num_failed_conversion
    num_spectra_files = len(spectra_fpaths)
    num_scans_processed = 0
    num_visualizations = 0
    num_PSMs, num_sig_PSMs = 0, 0
    num_decoy, num_sig_decoy = 0, 0
    if save_in_files:
        num_PSMs, num_sig_PSMs = save_full_results(cfg, prepared_fpaths)
        if cfg.fdr:
            num_decoy, num_sig_decoy = save_full_results(cfg, prepared_fpaths, decoy=True)

    for stats in all_stats:
        if stats is None:
            continue
        num_scans_processed += stats[0]
        num_visualizations += stats[1]

    stats_summary = {'num_spectra_files': num_spectra_files, 'num_spectra_scans': num_scans_processed,
                     'num_failed_conversion': num_failed_conversion, 'num_failed_processing': num_failed_processing,
                     'num_PSMs': num_PSMs + num_decoy, 'num_sig_PSMs': num_sig_PSMs + num_sig_decoy,
                     'num_visualizations': num_visualizations, 'num_decoy': num_decoy, 'num_sig_decoy': num_sig_decoy}
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

    global_work_dir = cfg.work_dirpath
    if cfg.pipeline.spec_nets_mode:
        cfg.work_dirpath = join(cfg.work_dirpath, 'stage_1_files')
        if not isdir(cfg.work_dirpath):
            os.makedirs(cfg.work_dirpath)

    start_time = datetime.now()
    stats_summary, all_PSMs = dereplicate(cfg, spectra_fpaths,
                                          return_all_PSMs=cfg.pipeline.spec_nets_mode,
                                          msg="\n\tStage 1: initial run of Dereplicator to get annotations "
                                              "for SpecNet nodes" if cfg.pipeline.spec_nets_mode else "",
                                          save_in_files=not cfg.pipeline.spec_nets_mode)
    finish_time = datetime.now()

    if cfg.pipeline.spec_nets_mode:
        cfg.work_dirpath = global_work_dir
        parse_and_propagate(cfg, spectra_fpaths, opts.spec_network, all_PSMs)
        cfg.params['spec_net_propagation'] = True
        stats_summary, _ = dereplicate(cfg, spectra_fpaths,
                                       msg="\n\tStage 3: running Dereplicator with list of candidate peptides")
        finish_time = datetime.now()

    append_smiles(cfg)
    if file_mapping is not None:
        update_fnames(cfg, file_mapping)
    extract_unique_peptides(cfg, append_fdr=cfg.fdr,
                            fdr_limit=cfg.params['fdr_limit'] if 'fdr_limit' in cfg.params else None)
    if hasattr(cfg.pipeline, 'extract_best') and cfg.pipeline.extract_best:
        extract_best(cfg, npdtools_init.external_tools_dir)
    print_stats(cfg, stats_summary, finish_time - start_time, config.dereplicator_p_value_threshold)
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

