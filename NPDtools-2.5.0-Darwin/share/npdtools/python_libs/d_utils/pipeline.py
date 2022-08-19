import os
import config
import common
from os.path import join, abspath, basename, isfile, islink, isdir, getsize, splitext
from log_utils import info, error, warning
from fname_utils import _get_summary_fpath, _get_best_spectra_fpath, _get_fpath_by_idx, _get_prefix_by_fpath, \
    _get_matches_fpath, update_fnames, _get_params_fpath, _get_local_matches_fpath
from common import sys_call, clean, get_param_command, get_mode_command, parse_extra_info
from spectra_utils import spectra_need_convertion, convert_to_mgf, spectra_final_ext
from file_utils import count_lines, is_complete_log
from d_utils.visual_utils import generate_visualisations

dereplicate_fpath = None
bin_dir = None
external_tools_dir = None


class PSM(object):
    header_line = "SpecFile\tLocalSpecIdx\tScan\tLocalPeptideIdx\tName\tScore\tP-Value\tPeptideMass\tSpectrumMass\tRetention\tAdduct\tCharge"

    def __init__(self, spectrum_fpath, line):
        self.spectrum_fpath = spectrum_fpath
        #CPP OUTPUT:  nlp_idx compound_name compound_idx nlp_mass spectrum_idx spec_mz rt charge score pvalue
        #EXAMPLE: antimarin_45013 [Val7]-Surfactin_C13ai_monomethyl_ester 1007.65 0 1008.69 -1 1 0 0.030884 0 -41 4.32898e-08
        self.nlp_idx, self.compound_name, self.compound_idx, self.nlp_mass, self.spectrum_idx, self.scan_id, self.spec_mz, self.rt, self.adduct, \
            self.charge, self.isotop, self.mass_error, self.decoy, self.score, self.p_value = line.split()
        self.p_value = "%.1e" % float(self.p_value)
        self.decoy = (self.decoy == '1')

    def __str__(self):
        # see header_line
        return "\t".join([self.spectrum_fpath, self.spectrum_idx, self.scan_id, self.compound_idx, self.compound_name,
                          self.score, self.p_value, self.nlp_mass, self.spec_mz, self.rt, self.adduct, self.charge])


def __get_dereplication_result_fpath_by_prefix(cfg, prefix):
    return join(cfg.work_dirpath, prefix + '_' + config.LOGGER_DEREPLICATOR_NAME + '.log')


def parse_output(spectrum_fpath, output_fpath):
    PSMs = []
    significant_PSMs = []
    num_scans = 0
    with open(output_fpath) as f:
        for line in f:
            if line.startswith("#"):
                if 'spectra loaded' in line:
                    try:
                        num_scans = int(line.split()[1])
                    except ValueError:
                        warning('Failed to parse number of loaded spectra (scans)')
                continue
            if line.startswith('.'):
                warning('Skipping "%s" from the Dereplicator output (starts with "." -- does not seem to be a PSM)!' % line.strip())
                continue
            if len(line.split()) < len(PSM.header_line.split('\t')) - 1:  # 1 PSM field is not read from the line
                warning('Skipping "%s" from the Dereplicator output (too few columns -- does not seem to be a PSM)!' % line.strip())
                continue
            PSMs.append(PSM(spectrum_fpath, line))
            if (config.USE_SCORE_ONLY and float(PSMs[-1].score) >= config.der_plus_score_threshold) or \
               (not config.USE_SCORE_ONLY and float(PSMs[-1].p_value) < config.dereplicator_p_value_threshold):
                significant_PSMs.append(PSMs[-1])
    return PSMs, significant_PSMs, num_scans


def save_full_results(cfg, prepared_fpaths, decoy=False, headers_to_exclude=[]):
    if 'extra_info' in cfg.params and cfg.params['extra_info'].strip():
        extra_headers, _ = parse_extra_info(cfg.params['extra_info'])
    else:
        extra_headers, _ = '', ''

    def __get_line_after_exclusion(line, idx_to_exclude, separator='\t'):
        if not idx_to_exclude:
            return line
        return separator.join([item for idx, item in enumerate(line.strip().split(separator))
                               if idx not in idx_to_exclude]) + '\n'

    # all matches
    num_all = 0
    with open(_get_matches_fpath(cfg, decoy=decoy), 'w') as all_matches:
        header_line = PSM.header_line + extra_headers + '\n'
        col_idx_to_exclude = [header_line.split('\t').index(col_name_to_exclude)
                              for col_name_to_exclude in headers_to_exclude
                              if col_name_to_exclude in header_line.split('\t')]
        all_matches.write(__get_line_after_exclusion(header_line, col_idx_to_exclude))
        for fpath in prepared_fpaths:
            if fpath is None or not isfile(_get_local_matches_fpath(fpath, decoy=decoy)):
                continue
            with open(_get_local_matches_fpath(fpath, decoy=decoy)) as cur_matches:
                for line in cur_matches:
                    all_matches.write(__get_line_after_exclusion(line, col_idx_to_exclude))
                    num_all += 1
    # significant matches
    num_sig = 0
    with open(_get_matches_fpath(cfg, 'sig', decoy=decoy), 'w') as sig_matches:
        header_line = 'VisualizationID\t' + PSM.header_line + extra_headers + '\n'
        col_idx_to_exclude = [header_line.split('\t').index(col_name_to_exclude)
                              for col_name_to_exclude in headers_to_exclude
                              if col_name_to_exclude in header_line.split('\t')]
        sig_matches.write(__get_line_after_exclusion(header_line, col_idx_to_exclude))
        for fpath in prepared_fpaths:
            if fpath is None or not isfile(_get_local_matches_fpath(fpath, 'sig', decoy=decoy)):
                continue
            with open(_get_local_matches_fpath(fpath, 'sig', decoy=decoy)) as cur_matches:
                for line in cur_matches:
                    sig_matches.write(__get_line_after_exclusion(line, col_idx_to_exclude))
                    num_sig += 1
    return num_all, num_sig


def save_single_spectra_results(cfg, base_fpath, PSMs, significant_PSMs_and_VisIDs):
    if 'extra_info' in cfg.params and cfg.params['extra_info'].strip():
        _, extra_values = parse_extra_info(cfg.params['extra_info'])
    else:
        _, extra_values = '', ''

    # all matches
    with open(_get_local_matches_fpath(base_fpath), 'w') as target:
        with open(_get_local_matches_fpath(base_fpath, decoy=True), 'w') as decoy:
            for psm in PSMs:
                if not psm.decoy:
                    target.write(str(psm) + extra_values + "\n")
                else:
                    decoy.write(str(psm) + extra_values + "\n")

    # significant matches
    with open(_get_local_matches_fpath(base_fpath, type='sig'), 'w') as target:
        with open(_get_local_matches_fpath(base_fpath, type='sig', decoy=True), 'w') as decoy:
            for i, (psm, vis_id) in enumerate(significant_PSMs_and_VisIDs):
                if not psm.decoy:
                    target.write(str(vis_id) + "\t" + str(psm) + extra_values + "\n")
                else:
                    decoy.write(str(vis_id) + "\t" + str(psm) + extra_values + "\n")


def preprocess_one(cfg, fpath, idx):
    prepared_fpath = _get_fpath_by_idx(cfg, idx, spectra_final_ext(fpath))
    prefix = _get_prefix_by_fpath(prepared_fpath)
    info('Starting: ' + fpath + ' --> ' + prefix)
    if cfg.pipeline.reuse and cfg.pipeline.vis_type == 'none':  # we need spectra for visualizing
        dereplication_result_fpath = __get_dereplication_result_fpath_by_prefix(cfg, prefix)
        if is_complete_log(dereplication_result_fpath):
            info('  Skipping conversion of input spectra (completed dereplication results exist)')
            return prepared_fpath

    if cfg.pipeline.reuse and isfile(prepared_fpath) and getsize(prepared_fpath) > 0:
        info('  Reusing existing spectra file ({prepared_fpath})'.format(**locals()))
    else:
        if isfile(prepared_fpath) or islink(prepared_fpath):
            os.remove(prepared_fpath)
        if spectra_need_convertion(fpath):
            if not convert_to_mgf(external_tools_dir, fpath, prepared_fpath):
                warning('Conversion of input spectra for %s (%s) failed, skipping this file' %
                        (prefix, basename(fpath)))
                return None
        else:
            os.symlink(abspath(fpath), prepared_fpath)
    return prepared_fpath


def process_one(cfg, fpath, prepared_fpath):
    PSMs, stats = [], None
    if prepared_fpath is not None:
        prefix = _get_prefix_by_fpath(prepared_fpath)
        cur_output = __get_dereplication_result_fpath_by_prefix(cfg, prefix)
        if cfg.pipeline.reuse and is_complete_log(cur_output):
            info('  Reusing completed dereplication results: ' + cur_output)
        else:
            command = dereplicate_fpath + ' "{prepared_fpath}" "{cfg.db_path}" "{cfg.library_info}" ' \
                                          '--configs_dir "{cfg.configs_dir}"'.format(**locals())
            command += get_mode_command(cfg.params)
            command += get_param_command(cfg.params, 'nps')
            command += get_param_command(cfg.params, 'max_charge')
            command += get_param_command(cfg.params, 'min_num_comp')
            command += get_param_command(cfg.params, 'accurate_p_value', '--x-metropolis')
            command += get_param_command(cfg.params, 'fdr', '--decoy_type 2 --decoy_fraction 1')
            command += get_param_command(cfg.params, 'debug')
            command += get_param_command(cfg.params, 'min_score')
            if hasattr(cfg.pipeline, 'variable_mode') and cfg.pipeline.variable_mode:
                command += get_param_command(cfg.params, 'max_mod', '--variable_filter')
                command += get_param_command(cfg.params, 'spec_net_propagation', '--variable_mode 3')
                command += get_param_command(cfg.params, 'varquest', '--variable_mode 2')
                command += get_param_command(cfg.params, 'min_shared_peaks')
                command += get_param_command(cfg.params, 'theor_spectra')
            else:
                command += get_param_command(cfg.params, 'isotope')
                command += get_param_command(cfg.params, 'adduct_Na', '--add_adduct Na')
                command += get_param_command(cfg.params, 'adduct_K', '--add_adduct K')
            if 'dereplicator+' in cfg.params:
                command += ' --fragmentation_tree --no_pvalue --min_num_comp 0' \
                           ' --fragment "{cfg.fragmentation_file}"'.format(**locals())
                if cfg.preprocessed_ft and isfile(cfg.preprocessed_ft):
                    command += ' --parse_precomputed_fragmentation_trees ' + cfg.preprocessed_ft
                if cfg.preprocessed_ftd and isfile(cfg.preprocessed_ftd):
                    command += ' --parse_precomputed_fragmentation_tree_decoys ' + cfg.preprocessed_ftd

            if cfg.pipeline.pass_to_dereplicate:
                command += ' ' + cfg.pipeline.pass_to_dereplicate

            command += ' > "{cur_output}"'.format(**locals())
            if not sys_call(command, cur_output):
                warning('Run on %s (%s) finished abnormally, skipping its results' % (prefix, basename(fpath)))
                cur_output = None
        if cur_output is not None:
            PSMs, significant_PSMs, num_scans = parse_output(fpath, cur_output)
            vis_ids = generate_visualisations(bin_dir, cfg, prepared_fpath, significant_PSMs)
            vis_count = len(vis_ids) - vis_ids.count(None)
            save_single_spectra_results(cfg, prepared_fpath, PSMs, zip(significant_PSMs, vis_ids))
            info('Run on %s (%s) finished correctly (%d scans processed). PSMs found: %d (significant: %d); '
                 'visualizations generated: %d' % \
                 (prefix, basename(fpath), num_scans, len(PSMs), len(significant_PSMs), vis_count))
            stats = (num_scans, vis_count)
        clean(cfg, prefix, exception='.matches')
    return PSMs, stats


def process_single_spectra(cfg, fpath, idx, return_PSMs=False):
    prepared_fpath = preprocess_one(cfg, fpath, idx)
    PSMs, stats = process_one(cfg, fpath, prepared_fpath)
    return prepared_fpath, PSMs if return_PSMs else [], stats


def append_smiles(cfg):
    if not isfile(cfg.library_smiles):
        return
    if count_lines(cfg.library_smiles) != len(cfg.db_entries):
        return

    smiles = list()
    with open(cfg.library_smiles) as f:
        for line in f:
            smiles.append(line.rstrip())

    sig_fpath = _get_matches_fpath(cfg, 'sig')
    with open(sig_fpath) as f:
        header_entries = f.readline().rstrip('\n').split('\t')
        compound_idx_column = header_entries.index('LocalPeptideIdx')
        content = f.readlines()
    with open(sig_fpath, 'w') as f:
        header_entries.append('SMILES')
        f.write('\t'.join(header_entries) + '\n')
        for line in content:
            entries = line.rstrip('\n').split('\t')
            compound_idx = int(entries[compound_idx_column])
            entries.append(smiles[compound_idx])
            f.write('\t'.join(entries) + '\n')


def save_params(cfg, dereplicator_plus_mode=False):
    params_fpath = _get_params_fpath(cfg)
    with open(params_fpath, 'w') as f:
        f.write('Parameter\tValue\n')

        mode, (pm_thresh, product_ion_thresh) = common.get_mode(cfg.params)
        #f.write('Mode\t{mode}\n'.format(**locals()))  # it is not informative because of the next two items
        f.write('Precursor Ion Mass Tolerance\t{pm_thresh} Da\n'.format(**locals()))
        f.write('Fragment Ion Mass Tolerance\t{product_ion_thresh} Da\n'.format(**locals()))

        f.write('Max charge\t' + common.get_param(cfg.params, 'max_charge') + '\n')
        if dereplicator_plus_mode:
            f.write('Min Score to Consider a PSM\t' + common.get_param(cfg.params, 'min_score') + '\n')
            f.write('Fragmentation Mode\t' + common.get_param(cfg.params, 'fragmentation_mode') + '\n')
        else:
            f.write('Max isotopic shift\t' + common.get_param(cfg.params, 'isotope') + '\n')

            adducts = []
            for adduct in ['Na', 'K']:
                if common.get_param(cfg.params, 'adduct_' + adduct) == 'True':
                    adducts.append(adduct)
            f.write('Adducts\t' + (', '.join(adducts) if adducts else 'None') + '\n')

            f.write('Minimum Number of Bonds\t' + common.get_param(cfg.params, 'min_num_comp') + '\n')
            f.write('Accurate P-values\t' + common.get_param(cfg.params, 'accurate_p_value', boolean=True) + '\n')
            f.write('Search Analogs\t' + common.get_param(cfg.params, 'varquest', boolean=True) + '\n')

    return params_fpath