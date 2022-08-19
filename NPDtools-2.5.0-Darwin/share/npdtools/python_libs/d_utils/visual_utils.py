#!/usr/bin/env python

import re
import json
import sys
from common import warning, info, sys_call_with_output, sys_call, get_mode_command
from config import MASS_PROTON
from os.path import join, dirname, splitext, basename, realpath, isdir, isfile
from fname_utils import _get_prefix_by_fpath
from distutils import dir_util
import string
import random
import os


def id_generator(size=10, chars=string.ascii_uppercase + string.digits):
    return ''.join(random.choice(chars) for _ in range(size))


class ScoreAnnotation(object):
    def __init__(self, anno_line):
        #Format: mass num_break (idx type) [num_break times] err peak_ind mz charge num_comp (0/1) [num_comp times] ion_type
        #Example: 99.0684 2 1 amide_1 3 amide_1 0.202702 56 909.388 1 9 1 0 0 0 0 0 0 0 0 y

        anno_entries = anno_line.strip().split(' ')
        self.mass = anno_entries[0]
        self.num_break = int(anno_entries[1])
        last_break_idx = 1 + self.num_break * 2
        break_iter = iter(anno_entries[2:last_break_idx + 1])
        self.breaks = zip(break_iter, break_iter)
        self.err, self.peak_ind = anno_entries[last_break_idx + 1:last_break_idx + 3]
        self.mz = anno_entries[last_break_idx + 3]
        self.charge = anno_entries[last_break_idx + 4]
        self.num_comp = int(anno_entries[last_break_idx + 5])
        last_comp_idx = last_break_idx + 5 + self.num_comp
        self.comps = anno_entries[last_break_idx + 6:last_comp_idx + 1]
        self.ion_type = anno_entries[last_comp_idx + 1]

    def get_dict(self):
        selected_component_flag = '1' if self.ion_type == 'b' else '0'
        components = [idx for idx, flag in enumerate(self.comps) if flag == selected_component_flag]
        emp_mass = float(self.mz) - MASS_PROTON

        result = {'peakIdx': int(self.peak_ind),
                  'theorMass': float(self.err) + emp_mass,
                  'charge': int(self.charge),
                  'components': components}
        return result


def parse_spectrum_annotation(ann_file):
    peaks = []
    annotation = []
    modification = {}
    struct_rule_fragmented_graph = []
    rfg_len = 0
    with open(ann_file) as f:
        is_spectrum = False
        for line in f:
            if not is_spectrum and not peaks:  # spectrum content is not started yet
                if line.startswith('BEGIN IONS'):
                    is_spectrum = True
            elif is_spectrum:  # in the middle of spectrum processing
                if line.startswith('END IONS'):
                    is_spectrum = False
                    continue
                if line[0].isdigit():
                    peaks.append(map(float, [line.split()[0], line.split()[1]]))
            else:  # spectrum already processed
                if line.startswith('number'):
                    struct_rule_fragmented_graph.append(line)
                    rfg_len += int(line.split()[-1]) + 1
                elif struct_rule_fragmented_graph and len(struct_rule_fragmented_graph) < rfg_len:
                    struct_rule_fragmented_graph.append(line)
                elif line.startswith('mod'):
                    value = line.split(':')[1].strip()
                    if 'position' in line:
                        modification['position'] = int(value)
                    elif 'node' in line:
                        modification['componentFormula'] = value.split()[0]
                        modification['componentMass'] = float(value.split()[1])
                    elif 'mass shift' in line:
                        modification['massShift'] = float(value)
                elif line[0].isdigit():
                    annotation.append(ScoreAnnotation(line).get_dict())
    if 'position' in modification and modification['position'] < 0:
        modification = {}
    return peaks, annotation, modification, struct_rule_fragmented_graph


def get_psm_details(bin_dir, cfg, spectra_fpath, PSM, vis_id):
    psm_identifier = "" + PSM.spectrum_fpath + " (P-Value: " + PSM.p_value + ")"
    info("    Creating visualization for " + psm_identifier, silent=not cfg.pipeline.debug)
    if int(PSM.compound_idx) >= len(cfg.db_entries):
        warning("Something is wrong with PSM compound idx (%d): it is larger than DB size (%d), skipping"
                % (int(PSM.compound_idx), len(cfg.db_entries)))
        return
    matched_mol_fpath = cfg.db_entries[int(PSM.compound_idx)]
    matched_mol_aux = join(dirname(matched_mol_fpath), 'aux', splitext(basename(matched_mol_fpath))[0] + '.json')
    if not isfile(matched_mol_aux) and not cfg.pipeline.vis_type.endswith('txt'):
        warning("Aux JSON is not found for %s, skipping" % basename(matched_mol_fpath))
        return

    anno_fpath = join(cfg.work_dirpath, vis_id + '.ann')
    psm_details = None
    try:
        # create ann
        print_score_fpath = join(bin_dir, 'print_score')
        command = print_score_fpath + ' "{spectra_fpath}" "{matched_mol_fpath}" ' \
                                      '--charge {PSM.charge} --no_pvalue --print_matches --print_spectrum ' \
                                      '--scan_num {PSM.spectrum_idx} --configs_dir {cfg.configs_dir}'.format(**locals())
        command += get_mode_command(cfg.params).replace("pm_thresh", "blind_search_pm_thresh")
        if cfg.pipeline.vis_type.endswith('txt'):
            command += ' --blind_rule_fragmented_graph' if cfg.pipeline.variable_mode else ' --print_rule_fragmented_graph'
        if cfg.pipeline.variable_mode:
            command += ' --blind_search'
        command += ' > "{anno_fpath}"'.format(**locals())
        if not sys_call(command, anno_fpath, indent='      ', silent=not cfg.pipeline.debug):
            warning("Annotation file was not created for " + psm_identifier + ", skipping")
            return

        spectrum, annotation, modification, structure = parse_spectrum_annotation(anno_fpath)
        psm_details = {'spectrum': spectrum,
                       'annotation': annotation}
        if modification:
            psm_details['modification'] = modification
        if structure:
            psm_details['structure'] = structure

        if not cfg.pipeline.vis_type.endswith('txt'):
            with open(matched_mol_aux) as f:
                aux_data = json.load(f)
            if 'mol' not in aux_data or 'atomComponents' not in aux_data or 'fragmentedBonds' not in aux_data:
                warning("Aux JSON for %s does not include 'mol', 'atomComponents', or 'fragmentedBonds', skipping"
                        % basename(matched_mol_fpath))
                return
            psm_details.update(aux_data)

        ## TODO: in case of non standard fragmentation (Non Peptidic)
        ## TODO: we need to recalculate atomComponents and fragmentedBonds
        # print_structure_fpath = join(bin_dir, 'print_structure')
        # if 'fragmentedBonds' not in aux_data:
        #     command = print_structure_fpath + ' "{matched_mol_fpath}" --print_rule_fragmented_bonds ' \
        #                                       '--remove_h ' \
        #                                       '--configs_dir {cfg.configs_dir}'.format(**locals())
        #     peptide_bonds_output = sys_call_with_output(command, indent='      ', silent=not cfg.pipeline.debug)
        #     if peptide_bonds_output is None:
        #         raise Exception("print_structure failed to process " + matched_mol_fpath)
        #     peptide_bonds = peptide_bonds_output.rstrip().split('\n')
        #     aux_data['peptideBonds'] = map(str, peptide_bonds)
    except Exception:
        exc_type, exc_value, _ = sys.exc_info()
        warning("Something went wrong on creating visualization for " + psm_identifier + ", skipping")
        warning("Exception: " + str(exc_type) + "; Details: " + str(exc_value))
    finally:
        return psm_details


def embed_html_aux_files(target_html_content, html_aux_dirpath):
    js_line_tmpl = '<script type="text/javascript" src="{file_rel_path}"></script>'
    js_l_tag = '<script type="text/javascript" name="{name}">'
    js_r_tag = '    </script>'

    css_line_tmpl = '<link rel="stylesheet" type="text/css" href="{file_rel_path}" />'
    css_l_tag = '<style type="text/css" rel="stylesheet" name="{name}">'
    css_r_tag = '    </style>'

    css_dirpath = join(html_aux_dirpath, 'css')
    js_dirpath = join(html_aux_dirpath, 'js')
    css_files = [join(css_dirpath, f) for f in os.listdir(css_dirpath) if f.endswith('css')]
    js_files = [join(js_dirpath, f) for f in os.listdir(js_dirpath) if f.endswith('js')]

    for line_tmpl, files, l_tag, r_tag in [
            (js_line_tmpl, js_files, js_l_tag, js_r_tag),
            (css_line_tmpl, css_files, css_l_tag, css_r_tag)]:
        for fpath in files:
            if not os.path.exists(fpath):
                continue

            rel_fpath = os.path.relpath(fpath, dirname(html_aux_dirpath))
            line = line_tmpl.format(file_rel_path=rel_fpath)
            l_tag_formatted = l_tag.format(name=rel_fpath)
            with open(fpath) as f:
                contents = f.read()
                contents = '\n'.join(' ' * 8 + l for l in contents.split('\n'))
                try:
                    target_html_content = ''.join(i for i in target_html_content if ord(i) < 128)
                    contents = ''.join(i for i in contents if ord(i) < 128)
                    target_html_content = target_html_content.replace(line, l_tag_formatted + '\n' + contents + '\n' + r_tag)
                except Exception:
                    continue
    return target_html_content


def save_in_text_format(psm_details, verbose=False):
    def __format_modification():
        if 'modification' not in psm_details:
            return ''
        result = '# MODIFICATION INFO:\n'
        for k in sorted(psm_details['modification'].keys(), key=len):
            result += str(k) + ': ' + str(psm_details['modification'][k]) + '\n'
        return result

    def __format_matched_peaks():
        result = '# MATCHED PEAKS INFO:\n'
        result += '#' + ' '.join(['m/z', 'intensity', 'theorCharge', 'theorMass']) + '\n'
        for matched_peak in psm_details['annotation']:
            peakIdx = matched_peak['peakIdx']
            theorCharge = matched_peak['charge']
            theorMass = matched_peak['theorMass']
            peak = psm_details['spectrum'][peakIdx]
            peakMZ, peakIntensity = peak[0], peak[1]
            result += ' '.join(map(str, [peakMZ, peakIntensity, theorCharge, theorMass])) + '\n'
        return result

    def __format_structure():
        if 'structure' not in psm_details:
            return ''
        if 'modification' in psm_details:
            result = '# MODIFIED/MUTATED STRUCTURE:\n'
            mod_pos = psm_details['modification']['position']
            mod_mass = psm_details['modification']['massShift']
            # appedning chemical formula in the modified residue with the mass shift
            mod_line_items = psm_details['structure'][mod_pos + 1].split(' ')
            mod_line_items[1] += ('+' if mod_mass >= 0.0 else '') + str(mod_mass)
            psm_details['structure'][mod_pos + 1] = ' '.join(mod_line_items)
        else:
            result = '# STRUCTURE:\n'
        if not verbose:
            result = ''
        result += ''.join(psm_details['structure'])
        return result

    full_result = ''
    if verbose:
        full_result += __format_modification()
        full_result += __format_matched_peaks()
    full_result += __format_structure()
    return full_result


def generate_visualisations(bin_dir, cfg, spectra_fpath, significant_PSMs):
    def __process_str(entry, max_len=100):
        entry = entry.replace('"', '').replace("'", '')  # quotes filter
        if len(entry) > max_len:
            entry = entry[:max_len - 3] + '...'
        return entry

    if cfg.pipeline.vis_type == 'none' or not significant_PSMs:
        return [None] * len(significant_PSMs)

    prefix = _get_prefix_by_fpath(spectra_fpath)
    html_basedir = dirname(realpath(__file__))
    html_aux_src = join(html_basedir, 'html_aux_files')
    html_aux_dst = join(cfg.vis_dirpath, basename(html_aux_src))

    vis_ids = []
    for i, PSM in enumerate(significant_PSMs):
        if PSM.decoy:
            vis_ids.append(None)
        else:
            vis_id = prefix + '_PSM_' + str(i) + ('_' + id_generator() if cfg.pipeline.vis_type == 'json' else '')
            psm_details = get_psm_details(bin_dir, cfg, spectra_fpath, PSM, vis_id)
            if not psm_details:
                vis_ids.append(None)
                continue
            vis_ids.append(vis_id)
            if cfg.pipeline.vis_type == 'json':
                json_fpath = join(cfg.vis_dirpath, vis_id + '.json')
                with open(json_fpath, 'w') as f:
                    f.write(json.dumps(psm_details))
                info("    JSON saved to " + json_fpath, silent=not cfg.pipeline.debug)
            elif cfg.pipeline.vis_type.endswith('html'):
                # add aux info about the spectrum and the compound
                psm_details['spectrumInfo'] = {'filename': __process_str(basename(PSM.spectrum_fpath)),
                                               'scan': int(PSM.scan_id),
                                               'charge': int(PSM.charge), 'mz': float(PSM.spec_mz),
                                               'retention': float(PSM.rt) if float(PSM.rt) > 0 else None}
                psm_details['compoundInfo'] = {'name': __process_str(PSM.compound_name),
                                               'mass': float(PSM.nlp_mass), 'adduct': PSM.adduct}
                #TODO: add chemical formula to compound info ('formula': 'C2H5OH')

                html_fpath = join(cfg.vis_dirpath, vis_id + '.html')
                html_template = join(html_basedir, 'template.html')
                with open(html_template) as template_file:
                    html_content = template_file.read()
                html_content = html_content.replace('{{ JSON }}', json.dumps(psm_details))
                html_content = html_content.replace('\\n', '\\\\n')
                if cfg.pipeline.vis_type == 'portable_html':
                    html_content = embed_html_aux_files(html_content, html_aux_src)
                elif not isdir(html_aux_dst):  # regular 'html' -- html_aux_files should be present in the output dir
                    dir_util.copy_tree(html_aux_src, html_aux_dst, preserve_times=False, preserve_mode=False)
                with open(html_fpath, 'w') as output_file:
                    output_file.write(html_content)
            elif cfg.pipeline.vis_type.endswith('txt'):
                txt_fpath = join(cfg.vis_dirpath, vis_id + '.txt')
                with open(txt_fpath, 'w') as f:
                    f.write(save_in_text_format(psm_details, verbose=(cfg.pipeline.vis_type == 'extended_txt')))
    return vis_ids