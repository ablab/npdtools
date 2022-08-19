import sys
import os
import glob
import shutil
import shlex
import subprocess
import xml.etree.ElementTree
import config
from distutils import dir_util
from log_utils import warning, info, error
from os.path import join, isfile, isdir, basename
from file_utils import check_file
from fname_utils import _get_matches_fpath, _get_summary_fpath, _get_best_spectra_fpath
import npdtools_init


def set_threads(cfg):
    try:
        from joblib import Parallel, delayed
    except ImportError:
        warning('Failed to import joblib, will do everything in 1 thread')
        cfg.threads = 1
        return
    if 'threads' in cfg.params and cfg.params['threads'] and str(cfg.params['threads']).isdigit():  # highest priority
        cfg.threads = int(cfg.params['threads'])
        return
    if 'variable_brute_force' in cfg.params and cfg.params['variable_brute_force']:  # special case for Dereplicator
        cfg.threads = 1
        return
    if cfg.threads:  # user defined number of threads
        return
    try:
        import multiprocessing
        cfg.threads = max(1, multiprocessing.cpu_count() // 2)
    except:
        warning('Failed to determine the number of CPUs, setting number of threads to default (%d)' % config.DEFAULT_THREADS)
        cfg.threads = config.DEFAULT_THREADS


def get_version():
    version_fpath = os.path.join(npdtools_init.docs_dir, 'VERSION.txt')
    if not os.path.isfile(version_fpath):
        return 'unknown'
    with open(version_fpath) as f:
        return f.readline().strip()


def sys_call(command, output_fpath=None, indent='  ', silent=False):
    info(indent + 'Running: ' + command, silent=silent)
    ret_code = os.system(command)  # TODO: switch to subprocess
    if ret_code == 0 and (output_fpath is None or check_file(output_fpath)):
        return True
    return False


def sys_call_with_output(command, indent='  ', silent=False):
    info(indent + 'Running: ' + command, silent=silent)
    p = subprocess.Popen(shlex.split(command), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()
    if p.returncode != 0:
        return None
    return out


# based on the similar function in https://github.com/ablab/quast/blob/master/quast_libs/qutils.py
def run_parallel(_fn, fn_args, n_jobs):
    parallel_args = {'n_jobs': n_jobs}
    try:
        from joblib import Parallel, delayed
        try:
            # starting from joblib 0.10 the default backend has changed to 'loky' which causes NPDtools crashes,
            # since it uses default values in config.py module. So, we explicitly require 'multiprocessing' here.
            # Note that Parallel doesn't have 'require' argument in joblib 0.9 and earlier.
            new_style_parallel_args = {'backend': 'multiprocessing'}
            Parallel(**new_style_parallel_args)
            parallel_args.update(new_style_parallel_args)
        except TypeError:
            pass

        results = Parallel(**parallel_args)(delayed(_fn)(*args) for args in fn_args)
        return results
    except ImportError:
        error('Critical problem: failed to import joblib while intended to work in more than 1 thread!')


def clean(cfg, prefix="", exception=None, intermediate=False):
    if cfg.pipeline.debug:
        return
    if exception is None and hasattr(cfg, "clean_exception"):
        exception = cfg.clean_exception
    if exception and type(exception) is not list:
        exception = [exception]
    for path in glob.glob(join(cfg.work_dirpath, prefix + "*")):
        if path.endswith("log") or path.endswith(".tsv"):
            continue
        # special case: '.edit' is a suffix of RiPPquest genome analysis files
        if path.endswith(".edit") and (intermediate or
                                       'keep_ga_files' in cfg.pipeline.__dict__ and cfg.pipeline.keep_ga_files):
            continue
        if path.startswith(cfg.configs_dir):
            continue
        if exception:
            skip_this = False
            for exc in exception:
                if exc and path.endswith(exc):
                    skip_this = True
                    break
            if skip_this:
                continue
        if isfile(path):
            os.remove(path)
        elif isdir(path):
            shutil.rmtree(path)


def _clean_matches(cfg):
    for type in ['all', 'sig']:
        for decoy in [True, False]:
            report_fpath = _get_matches_fpath(cfg, type, decoy)
            if isfile(report_fpath):
                os.remove(report_fpath)


def parse_params_xml(fpath):
    with open(fpath) as f:
        content = f.read()
    params = dict()
    file_mapping = dict()
    for e in xml.etree.ElementTree.fromstring(content).findall('parameter'):
        if e.attrib['name'] == 'upload_file_mapping':
            fname, real_fname = e.text.split('|')[0:2]
            file_mapping[basename(fname)] = real_fname
            continue
        if e.text == 'on':
            value = True
        elif e.text == 'off':
            value = False
        else:
            value = e.text
        params[e.attrib['name']] = value
    return params, file_mapping


def get_param_command(params, name, option=None):
    if name in params and params[name] is not False:
        value = '' if type(params[name]) == bool else params[name]
        if option is None:
            option = '--' + name
        return ' {option} {value}'.format(**locals())
    else:
        return ''


def get_param(params, name, boolean=False):
    if name in params:
        if boolean:
            return 'Yes' if params[name] else 'No'
        return str(params[name])
    return 'No' if boolean else 'Not set'


def __get_mode_threshs(mode, params):
    if mode == "custom":
        pm_thresh = float(params['pm_thresh'])
        product_ion_thresh = float(params['product_ion_thresh'])
    else:
        pm_thresh = 0.5 if mode.startswith('L') else 0.02
        product_ion_thresh = 0.5 if mode.endswith('L') else 0.02
    if 'ppm' in params and params['ppm']:
        pm_thresh /= 1000000
        product_ion_thresh /= 1000000
        if mode != "custom":  # adjust default values
            pm_thresh *= 1000
            product_ion_thresh *= 1000
    return pm_thresh, product_ion_thresh


def get_mode_command(params):
    mode = params['mode']
    pm_thresh, product_ion_thresh = __get_mode_threshs(mode, params)
    mode_str = ' --product_ion_thresh {product_ion_thresh}'.format(**locals())
    mode_str += ' --pm_thresh {pm_thresh}'.format(**locals())
    if 'ppm' in params and params['ppm']:
        mode_str += ' --ppm'
    return mode_str


def get_mode(params):
    mode = params['mode']
    if mode == 'custom':
        return mode, __get_mode_threshs(mode, params)
    human_readable_names = {'L': 'low', 'H': 'high'}
    return human_readable_names[mode[0]] + '-' + human_readable_names[mode[1]], __get_mode_threshs(mode, params)


def parse_extra_info(extra_info):
    pairs = extra_info.split('|')
    try:
        keys = [pair[:pair.find(':')] for pair in pairs]
        values = [pair[pair.find(':') + 1:] for pair in pairs]
    except IndexError:
        warning('Failed to parse extra_info!')
        keys = []
        values = []
    return '\t' + '\t'.join(keys), '\t' + '\t'.join(values)


def __fdr_str(targets, decoys):
        return '%.2f' % (decoys * 100.0 / (targets + decoys))


def extract_unique_peptides(cfg, append_fdr=True, fdr_limit=None, peptide_column='Name'):
    sig_matches_fpath = _get_matches_fpath(cfg, type='sig')
    sig_decoy_matches_fpath = _get_matches_fpath(cfg, type='sig', decoy=True)
    sig_unique_matches_fpath = _get_matches_fpath(cfg, type='sig', unique=True)
    try:
        fdr_limit = float(str(fdr_limit))  # in case it is in string format
    except ValueError:
        fdr_limit = None

    class Hit(object):
        def __init__(self, match_id, peptide, quality, charge, target=True):
            self.match_id = match_id  # make sense only for targets
            self.peptide = peptide
            self.quality = quality
            self.charge = charge
            self.target = target

    hits = []
    matches = []
    header_columns = None
    max_charge = int(cfg.params["max_charge"])
    if config.USE_SCORE_ONLY:
        quality_column = 'Score'
    else:
        quality_column = 'P-Value'
    for fpath, target in [(sig_matches_fpath, True), (sig_decoy_matches_fpath, False)]:
        if check_file(fpath):
            with open(fpath) as f:
                for i, line in enumerate(f):
                    line = line.rstrip('\n')
                    if not line:
                        continue
                    l = line.split('\t')
                    if i == 0:
                        header = l
                        if target:
                            header_columns = header
                        continue
                    hits.append(Hit(match_id=len(matches),
                                    peptide=l[header.index(peptide_column)],
                                    quality=float(l[header.index(quality_column)]) * (-1. if config.USE_SCORE_ONLY else 1.),
                                    charge=int(l[header.index('Charge')]), target=target))
                    if hits[-1].charge > max_charge:
                        max_charge = hits[-1].charge
                    if target:
                        matches.append(l)

    hits = sorted(hits, key=lambda x: (x.quality, x.charge, x.target))
    unique_peptides = set()
    unique_decoy_peptides = set()

    unique_target_per_charge = [0] * max_charge
    unique_decoy_per_charge = [0] * max_charge
    all_target_per_charge = [0] * max_charge
    all_decoy_per_charge = [0] * max_charge

    processed_all_matches = []
    processed_unique_matches = []

    for cur_hit in hits:
        cur_charge = cur_hit.charge - 1
        if cur_hit.target:
            cur_match = matches[cur_hit.match_id]
            all_target_per_charge[cur_charge] += 1
            if append_fdr:
                cur_match.append(__fdr_str(all_target_per_charge[cur_charge],
                                           all_decoy_per_charge[cur_charge]))
            processed_all_matches.append(cur_match)
            if cur_hit.peptide not in unique_peptides:
                unique_target_per_charge[cur_charge] += 1
                unique_peptides.add(cur_hit.peptide)
                cur_unique_match = list(cur_match)
                if append_fdr:
                    cur_unique_match[-1] = __fdr_str(unique_target_per_charge[cur_charge],
                                           unique_decoy_per_charge[cur_charge])
                processed_unique_matches.append(cur_unique_match)
        else:
            all_decoy_per_charge[cur_charge] += 1
            if cur_hit.peptide not in unique_decoy_peptides:
                unique_decoy_per_charge[cur_charge] += 1
                unique_decoy_peptides.add(cur_hit.peptide)

    if append_fdr:
        header_columns.append('FDR')
    with open(sig_matches_fpath, 'w') as f:
        f.write('\t'.join(header_columns) + '\n')
        for match in processed_all_matches:
            if fdr_limit is not None and append_fdr and float(match[-1]) >= fdr_limit:  # skipping entries with high FDR
                continue
            f.write('\t'.join(match) + '\n')
    with open(sig_unique_matches_fpath, 'w') as f:
        f.write('\t'.join(header_columns) + '\n')
        for match in processed_unique_matches:
            if fdr_limit is not None and append_fdr and float(match[-1]) >= fdr_limit:  # skipping entries with high FDR
                continue
            f.write('\t'.join(match) + '\n')


def write_dict_to_tsv(fpath, header_columns, data, filter_func=lambda x: x):
    with open(fpath, 'w') as f:
        f.write('\t'.join(header_columns) + '\n')
        for k, v in data.items():
            f.write(k + '\t' + filter_func(v) + '\n')


def copy_configs(configs_src_base, work_dir):
    configs_dst = join(work_dir, 'configs')
    if not isdir(configs_dst):
        os.makedirs(configs_dst)
    for subdir in config.config_dirs:
        dir_util.copy_tree(join(configs_src_base, subdir), join(configs_dst, subdir),
                           preserve_times=False, preserve_mode=False)
    return configs_dst


def get_command_line(args=None):
    if args is None:
        args = sys.argv
    cmd = ''
    for arg in args:
        if ' ' in arg:
            cmd += '"' + arg + '" '
        else:
            cmd += arg + ' '
    return cmd


def print_stats(cfg, stats_summary, elapsed=None, p_value_thresh=config.DEFAULT_P_VALUE_THRESHOLD,
                score_thresh=config.DEFAULT_DER_PLUS_SIG_MIN_SCORE):
    '''
    mandatory fields:
        num_spectra_files, num_failed_processing,
        num_PSMs, num_sig_PSMs
        num_decoy, num_sig_decoy
    optional fields:
        num_sequence_files (RiPPquest only)
        num_failed_conversion
        num_visualizations
    '''

    num_spectra_files = stats_summary['num_spectra_files']
    num_failed_processing = stats_summary['num_failed_processing']
    num_PSMs = stats_summary['num_PSMs']
    num_sig_PSMs = stats_summary['num_sig_PSMs']
    num_decoy = stats_summary['num_decoy']
    num_sig_decoy = stats_summary['num_sig_decoy']
    if 'num_sequence_files' in stats_summary:
        num_sequence_files = stats_summary['num_sequence_files']
    else:
        num_sequence_files = None
    if 'num_failed_conversion' in stats_summary:
        num_failed_conversion = stats_summary['num_failed_conversion']
    else:
        num_failed_conversion = None
    if 'num_visualizations' in stats_summary:
        num_visualizations = stats_summary['num_visualizations']
    else:
        num_visualizations = None
    if 'num_spectra_scans' in stats_summary:
        num_spectra_scans = stats_summary['num_spectra_scans']
    else:
        num_spectra_scans = None

    summary_fpath = _get_summary_fpath(cfg)
    if not config.USE_SCORE_ONLY:
        sig_thresh_msg = 'p-value < %.0e' % p_value_thresh
    else:
        sig_thresh_msg = 'score >= %.1f' % score_thresh
    if cfg.fdr:
        decoy_PSM_str = ' [decoy: %d]' % stats_summary['num_decoy']
        sig_decoy_PSM_str = ' [decoy: %d]' % stats_summary['num_sig_decoy']
    else:
        decoy_PSM_str = ''
        sig_decoy_PSM_str = ''
    msg = '\n\n'
    msg += '  Summary (saved to {summary_fpath}):\n'.format(**locals())
    if num_sequence_files is not None:
        msg += '    Number of sequence files: {num_sequence_files}\n'.format(**locals())
    msg += '    Number of spectrum files: {num_spectra_files}\n'.format(**locals())
    if num_spectra_scans is not None:
        msg += '    Number of spectrum scans: {num_spectra_scans}\n'.format(**locals())
    if num_failed_conversion is not None:
        msg += '      out of them: {num_failed_conversion} failed conversion to mgf\n'.format(**locals())
    msg += '      out of them: {num_failed_processing} failed processing\n'.format(**locals())
    msg += '    Found PSMs: {num_PSMs}{decoy_PSM_str}\n' \
           '      out of them: {num_sig_PSMs} significant ({sig_thresh_msg}){sig_decoy_PSM_str}\n'.format(**locals())
    if num_visualizations is not None:
        msg += '      out of them: {num_visualizations} visualizations generated\n'.format(**locals())
    msg += '\n  Matches are in %s and %s\n' % (_get_matches_fpath(cfg), basename(_get_matches_fpath(cfg, 'sig')))
    if hasattr(cfg.pipeline, 'vis_type') and cfg.pipeline.vis_type != 'none':
        msg += '  Visualizations are in ' + cfg.vis_dirpath + '\n'
    if hasattr(cfg.pipeline, 'extract_best') and cfg.pipeline.extract_best:
        msg += '  Best spectra are in ' + _get_best_spectra_fpath(cfg) + '\n'
    msg += '  Log is in ' + cfg.log_fpath + '\n'
    if cfg.pipeline.debug:
        msg += '  Intermediate files are in ' + cfg.work_dirpath + '\n'
    if elapsed is not None:
        msg += '\n  Elapsed time: ' + str(elapsed) + '\n'
    info(msg)

    with open(summary_fpath, 'w') as f:
        f.write('Statistic\tValue\n')
        if num_sequence_files is not None:
            f.write('Number of sequence files\t{num_sequence_files}\n'.format(**locals()))
        f.write('Number of spectrum files\t{num_spectra_files}\n'.format(**locals()))
        if num_spectra_scans is not None:
            f.write('Number of spectrum scans\t{num_spectra_scans}\n'.format(**locals()))
        if num_failed_conversion:
            f.write('Failed to convert to MGF\t{num_failed_conversion}\n'.format(**locals()))
        if num_failed_processing:
            f.write('Failed to process\t{num_failed_processing}\n'.format(**locals()))
        f.write('Number of all PSMs\t%d\n' % (num_PSMs - num_decoy))
        f.write('Number of significant PSMs (%s)\t%d\n' %
                (sig_thresh_msg, num_sig_PSMs - num_sig_decoy))
        #TODO: number of unique peptides
        #TODO: FDR at multiple levels
        total_sig = num_sig_PSMs  # this variable includes both decoy and target
        f.write('FDR %% for %s\t%s\n' %
                (sig_thresh_msg, ('%.1f' % (100.0 * num_sig_decoy / total_sig)) if total_sig else 'N/A'))
        if elapsed is not None:
            f.write('Elapsed time\t{elapsed}\n'.format(**locals()))
