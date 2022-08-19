#!/usr/bin/env python

import sys
import os
import subprocess
from os.path import isdir, isfile, join, basename
import logging
from site import addsitedir
from optparse import OptionParser

import npdtools_init
npdtools_init.init()

addsitedir(npdtools_init.python_modules_dir)
import config
import common
import log_utils
from log_utils import info, error, warning

bin_dir = npdtools_init.bin_dir
configs_dir = npdtools_init.configs_dir
# FOR Fragmentation trees
print_structure = join(bin_dir, 'print_structure')
ft_cacher = join(bin_dir, 'cache')

# HARD LIMITS/CONSTANTS
MAX_MASS = None  # in Da or None
MAX_MASS_FT = 3000  # in Da
MIN_NUM_BONDS = 0
MAX_TIME = 60*60  # in seconds; set via `ulimit -t MAX_TIME`
MAX_MEM = 10**7  # in kbytes; set via `ulimit -v MAX_MEM` (max memory size and virtual memory)
SUFFIX = '.cereal'
DECOY_SUFFIX = SUFFIX + '.decoy'
DB_SUFFIX = '.bin'
LIB_INFO_DEFAULT_NAME = 'library.info'

log_utils.log = logging.getLogger(config.LOGGER_DEREPLICATOR_NAME)


def get_mol_fpaths(dirpath, library_info_fpath, min_num_bonds=None, max_mass=None):
    mol_fpaths = []
    with open(library_info_fpath) as f:
        for line in f:
            info = line.strip().split(' ')
            mol_fpath = join(dirpath, info[0])
            mass = float(info[2])
            num_bonds = int(info[3])
            if (min_num_bonds is None or num_bonds >= min_num_bonds) and \
               (max_mass is None or mass <= max_mass):
                mol_fpaths.append(mol_fpath)
    return mol_fpaths


def create_ft_for_mol(mol_fpath, output_dirpath, fragmentation_rule_fpath, configs_dir, decoy=False, debug=False):
    serialized_fpath = join(output_dirpath, basename(mol_fpath) + (DECOY_SUFFIX if decoy else SUFFIX))
    ulimit = 'ulimit -t %d; ulimit -v %d; ' % (MAX_TIME, MAX_MEM)

    cmd = [print_structure, mol_fpath, '--fragment', fragmentation_rule_fpath,
           '--export_fragmentation_summary', serialized_fpath,
           '--fragmentation_tree', '--break_double', '--multiple_cut', '--configs_dir', configs_dir]
    if decoy:
        cmd.append('--decoy')

    cmd = ulimit + ' '.join(cmd)
    if debug:
        info("Running: " + cmd)
    ret_code = subprocess.call(cmd, shell=True, stderr=open('/dev/null', 'w'))
    if ret_code != 0:
        warning('Failed to process ' + basename(mol_fpath))
        return False
    return True


def create_ft_for_mols(mol_fpaths, output_dirpath, fragmentation_rule_fpath, configs_dir, threads=1, decoy=False, debug=False):
    n_jobs = min(len(mol_fpaths), threads)

    info("  Processing %d MOL files in %d threads" % (len(mol_fpaths), n_jobs))

    if n_jobs > 1:
        parallel_args = [(mol_fpath, output_dirpath, fragmentation_rule_fpath, configs_dir, decoy, debug)
                         for mol_fpath in mol_fpaths]
        results = common.run_parallel(create_ft_for_mol, parallel_args, n_jobs=n_jobs)
    else:
        results = []
        for mol_fpath in mol_fpaths:
            result = create_ft_for_mol(mol_fpath, output_dirpath, fragmentation_rule_fpath, configs_dir, decoy=decoy, debug=debug)
            results.append(result)
    info("  Correctly processed %d out of %d, saved under %s" % (results.count(True), len(mol_fpaths), output_dirpath))


def create_combined_db(output_dirpath, library_info, result_fpath, decoy=False, debug=False):
    cmd = [ft_cacher, '-i', library_info, '-c', output_dirpath, '-o', result_fpath,
           '--suffix', DECOY_SUFFIX if decoy else SUFFIX]
    if debug:
        cmd.append('--verbose')

    cmd = ' '.join(cmd)
    if debug:
        info("Running: " + cmd)
    ret_code = subprocess.call(cmd, shell=True, stderr=open('/dev/null', 'w'))
    if ret_code != 0:
        error('Failed to created the preprocessed DB in ' + result_fpath)
        return False
    return True


def clean_up(output_dirpath, new_library_info, decoy=False):
    with open(new_library_info) as f:
        for line in f:
            mol_file = line.split(' ')[0]
            to_remove = join(output_dirpath, basename(mol_file) + (DECOY_SUFFIX if decoy else SUFFIX))
            if isfile(to_remove):
                os.remove(to_remove)


def preprocess_db(mode, work_dirpath, db_path, library_info, configs_dir, opts):
    if 'max_mem' in opts:
        global MAX_MEM
        MAX_MEM = opts['max_mem']
    threads = int(opts['threads']) if 'threads' in opts else 1
    debug = bool(opts['debug']) if 'debug' in opts else False

    if not isdir(work_dirpath):
        os.makedirs(work_dirpath)
    if debug:
        log_fpath = join(work_dirpath, 'preprocessing.log')
        log_handler = logging.FileHandler(log_fpath, mode='w')
        log_utils.setup_log_handler(log_handler)

    if mode.startswith('ft'):
        global MAX_MASS
        MAX_MASS = MAX_MASS_FT

    mol_fpaths = get_mol_fpaths(db_path, library_info, min_num_bonds=MIN_NUM_BONDS, max_mass=MAX_MASS)
    if mode.startswith('ft'):
        assert 'fragmentation_mode' in opts, "Parameter fragmentation_mode is not specified for FT/FTD mode!"
        fragmentation_rule_fpath = join(configs_dir, 'configs', 'common', 'fragment_list_' +
                                        opts['fragmentation_mode'] + '.txt')
        final_result_fpath = join(work_dirpath, basename(library_info) + '.' + opts['fragmentation_mode'] + '.'
                                                                             + mode + DB_SUFFIX)
        if debug:
            info("Preprocessed DB will be saved in " + final_result_fpath)
        create_ft_for_mols(mol_fpaths, work_dirpath, fragmentation_rule_fpath, configs_dir, decoy=(mode == 'ftd'),
                           threads=threads, debug=debug)
        create_combined_db(work_dirpath, library_info, final_result_fpath, decoy=(mode == 'ftd'),
                           debug=debug)
        if not debug:
            clean_up(work_dirpath, library_info, decoy=(mode == 'ftd'))
        return final_result_fpath
    return None


def main():
    log_utils.log.setLevel(logging.DEBUG)
    log_utils.setup_log_handler(logging.StreamHandler(sys.stdout))

    parser = OptionParser(description='')
    parser.add_option('-o', '--output', dest='output_dirpath', help='Output dir (required)')
    parser.add_option('--db-path', help='Path to a PNP database (required)')
    parser.add_option('-l', '--library-info',
                      help='Path to the database description file (default: library.info inside the database dir)')
    parser.add_option('-m', '--mode', help='Type of preprocessing: '
                                           'ft -- fragmentation trees; '
                                           'ftd -- fragmentation trees (decoys); '
                                           'vq -- VarQuest [default: %default]',
                      choices=['ft', 'ftd', 'vq'], default='ft')
    parser.add_option('--fragmentation_mode', help='Mode of fragmentation (only for FT preprocessing, '
                                                   'e.g. "general_6_1_6" or "general_3_1_3") '
                                                   '[default: %default]', default='general_6_1_6')
    parser.add_option('-t', '--threads', type='int', help='Threads number [default: %default]', default=2)
    parser.add_option('--max-mem', type='int', help='Max allowed RAM (in kbytes), skip preprocessing of a structure '
                                                    'if it exceeds the limit [default: %default]', default=MAX_MEM)
    parser.add_option('-d', '--debug', help='Debug mode (more verbose and does not remove intermediate files)',
                      action="store_true", default=False)

    if len(sys.argv) < 2:
        parser.print_help()
        sys.exit(0)

    (opts, _) = parser.parse_args()
    if not opts.output_dirpath:
        parser.error("Output dir is required!")
    if not opts.db_path or not isdir(opts.db_path):
        parser.error("PNP database dir is not specified or does not exist!")
    if not opts.library_info:
        opts.library_info = join(opts.output_dirpath, LIB_INFO_DEFAULT_NAME)
    if not opts.library_info or not isfile(opts.library_info):
        parser.error("Library info file is not specified or does not exist!")

    assert opts.mode.startswith('ft'), "Only FT/FTD is currently implemented!"

    preprocess_db(opts.mode, opts.output_dirpath, opts.db_path, opts.library_info,
                  npdtools_init.configs_dir, opts.__dict__)


if __name__ == '__main__':
    main()