############################################################
#  For processing antiSMASH ORFs
#  All of them are combined in <output_dir>/*.final.gbk
############################################################


from os.path import isfile, join, isdir
import shutil
from r_utils.pipeline import get_genome_analysis_fpath, get_prod_classes_to_process, get_prod_class_mod_fpath
from log_utils import info, warning
from fasta_utils import write_fasta
from file_utils import is_complete_log
from common import sys_call
import config
from config import neighbours_distance
from collections import OrderedDict
import re
import os


hmmsearch_fpath = None
hmmsearch_models_dir = None


def get_all_CDS(gbk_fpath):
    ## NOTE: Unfortunately, we can't use:
    #  from Bio import SeqIO
    #  genome=SeqIO.read(gbk_fpath,'genbank')
    #  ...
    ## because BioPython is not available at ccms-login-temp.
    ## However we need rather simple logic here

    # GenBank format keywords
    features_start_keyword = 'FEATURES'
    features_end_keyword = 'ORIGIN'
    feature_margin = '     '
    cds_feature_name = 'CDS'  # we are interested in CDS only
    field_prefix = '/'
    # interesting fields
    sec_met_field = 'sec_met'
    sec_met_kind = 'biosynthetic'
    locus_tag_field = 'locus_tag'
    translation_field = 'translation'

    class CDS(object):
        def __init__(self, cds):
            self.start, self.end = self._parse_coordinates(cds[0])
            self.locus_tag = None
            self.translation = None
            self.is_proper_sec_met = False

            cur_field = ''
            for line in cds:
                line = line.strip()
                if line.startswith(field_prefix):
                    self._assign_attributes(cur_field)
                    cur_field = line
                else:
                    cur_field += line
            self._assign_attributes(cur_field)

        def _parse_field(self, line):
            field_id, field_value = None, None
            if line.startswith(field_prefix):
                field_id = line.split('=')[0][len(field_prefix):]
                field_value = line[len(field_prefix) + len(field_id) + 1:].strip('"')
            return field_id, field_value

        def _parse_coordinates(self, header_line):
            start, end = map(int, re.findall(r'\d+', header_line))
            return start, end

        def _assign_attributes(self, cur_field):
                cur_id, cur_value = self._parse_field(cur_field)
                if cur_id == locus_tag_field:
                    self.locus_tag = cur_value
                elif cur_id == translation_field:
                    self.translation = cur_value
                elif cur_id == sec_met_field and cur_value == 'Kind: ' + sec_met_kind:
                    self.is_proper_sec_met = True

    def _is_proper_locus(locus_cds):
        for cds in locus_cds:
            if cds.is_proper_sec_met:
                return True
        return False

    all_cds = []
    locus_cds = []
    inside_features = False
    cur_feature = []
    cur_feature_is_cds = False
    with open(gbk_fpath) as f:
        for line in f:
            if inside_features:
                # end of features list or a new feature is started
                if line.startswith(features_end_keyword) or line[len(feature_margin) + 1].isalpha():
                    if line.startswith(features_end_keyword):
                        inside_features = False
                        if _is_proper_locus(locus_cds):
                            all_cds.append(locus_cds)
                        locus_cds = []
                    if cur_feature_is_cds:
                        locus_cds.append(CDS(cur_feature))
                    cur_feature = []
                    cur_feature_is_cds = (line.split()[0] == cds_feature_name)
                cur_feature.append(line)
            elif line.startswith(features_start_keyword):
                inside_features = True
    return all_cds


def write_sec_metabolites_with_neighbours(all_cds, fpath):
    def _get_neighbour_seqs(target_cds, locus_cds):
        start = target_cds.start - neighbours_distance
        end = target_cds.end + neighbours_distance
        neighbour_seqs = []
        for cds in locus_cds:
            if start <= cds.start <= end or start <= cds.end <= end:
                neighbour_seqs.append(cds.translation)
        return neighbour_seqs

    sec_metabolites_itself = OrderedDict()
    sec_metabolites_with_neighbours = []
    for locus_cds in all_cds:
        for cds in locus_cds:
            if cds.is_proper_sec_met:
                name = cds.locus_tag
                sec_metabolites_itself[name] = cds.translation
                seq = '*'.join(_get_neighbour_seqs(cds, locus_cds))
                sec_metabolites_with_neighbours.append((name, seq))
    write_fasta(sec_metabolites_with_neighbours, fpath, width=None)
    return sec_metabolites_itself


def get_enzymes(mod_fpath, type='enabler'):  # types: enabler, single, double
    type_to_enzyme_column = {'enabler': 1, 'single': 6, 'double': 9}
    assert type in type_to_enzyme_column, "incorrect enzyme type!"

    enzymes = []
    with open(mod_fpath) as f:
        for line in f:
            if line.startswith(type):
                enzymes.append(line.split()[type_to_enzyme_column[type]])
                if enzymes[-1] == 'all':  # special case: do not need "all"
                    enzymes.pop()
    return enzymes


def run_hmmsearch(output_dir, sec_met_fpath, enzyme, silent=True):
    log_fpath = join(output_dir, enzyme + '.log')
    result_fpath = join(output_dir, enzyme + '.txt')
    model_fpath = join(hmmsearch_models_dir, enzyme + '.hmm')
    if isfile(result_fpath):
        return result_fpath
    cmd = hmmsearch_fpath + ' --max -o {log_fpath} --domtblout {result_fpath} ' \
                            '{model_fpath} {sec_met_fpath}'.format(**locals())
    if sys_call(cmd, result_fpath, indent='        ', silent=silent):
        return result_fpath
    else:
        return None


def get_ids_with_enzyme(hmmsearch_output_fpath):
    e_value_column_idx = 6
    e_value_threshold = 1e-05
    seq_ids = []
    with open(hmmsearch_output_fpath) as f:
        for line in f:
            if line.startswith('#'):
                continue
            try:
                if float(line.split()[e_value_column_idx]) < e_value_threshold:
                    seq_ids.append(line.split()[0])
            except ValueError:
                pass
    return seq_ids


def get_corresponding_enzymes(cfg, hmmsearch_output_dir, sec_met_fpath, prod_class):
    mod_fpath = get_prod_class_mod_fpath(prod_class)
    sec_met_name_to_enzymes = {}
    enablers = get_enzymes(mod_fpath, 'enabler')
    mutations = get_enzymes(mod_fpath, 'single') + get_enzymes(mod_fpath, 'double')
    for enzyme in enablers + mutations:
        hmmsearch_output = run_hmmsearch(hmmsearch_output_dir, sec_met_fpath, enzyme, silent=not cfg.pipeline.debug)
        if hmmsearch_output:
            sec_met_ids = get_ids_with_enzyme(hmmsearch_output)
            for id in sec_met_ids:
                if id not in sec_met_name_to_enzymes and enzyme in enablers:
                    sec_met_name_to_enzymes[id] = set()
                if id in sec_met_name_to_enzymes and enzyme in mutations:
                    sec_met_name_to_enzymes[id].add(enzyme)
        else:
            warning('hmmsearch failed for {enzyme} enzyme on {sec_met_fpath}'.format(**locals()))
    return sec_met_name_to_enzymes


def write_sec_metabolites_genome_analysis_results(ga_fpath, sec_met_name_to_seq, sec_met_name_to_enzymes):
    def _split_and_remove_X_if_needed(seq):
        split = seq.split('XXX')
        x_free = [subseq.replace('X', '') for subseq in split]
        return [subseq for subseq in x_free if subseq]

    fasta_entries = []
    for name, seq in sec_met_name_to_seq.items():
        if name in sec_met_name_to_enzymes:
            enzymes = sec_met_name_to_enzymes[name]
            enzymes_suffix = '_%d' % len(enzymes) + ('_' if enzymes else '') + '_'.join(enzymes)
            correct_seqs = _split_and_remove_X_if_needed(seq)
            for idx, cor_seq in enumerate(correct_seqs):
                part_suffix = '.part%d' % idx if idx else ''
                fasta_entries.append((name.replace('_', '.') + part_suffix + enzymes_suffix, cor_seq))
    if fasta_entries:
        write_fasta(fasta_entries, ga_fpath, width=None)
        return ga_fpath
    return None


def sort_sec_metabolites_into_production_classes(cfg, all_cds, sequence_fpath, prod_class):
    sec_met_fpath = sequence_fpath + '.sec_met_with_neighbours.fasta'
    sec_met_name_to_seq = write_sec_metabolites_with_neighbours(all_cds, sec_met_fpath)

    hmmsearch_output_dir = sequence_fpath + '_ripp_hmm_search'
    if isdir(hmmsearch_output_dir) and not cfg.pipeline.reuse:
        shutil.rmtree(hmmsearch_output_dir)
    if not isdir(hmmsearch_output_dir):
        os.makedirs(hmmsearch_output_dir)
    ga_fpaths = []
    tt = get_prod_classes_to_process(prod_class)
    for prod_class in get_prod_classes_to_process(prod_class):
        ga_fpath = get_genome_analysis_fpath(sequence_fpath, prod_class)
        if cfg.pipeline.reuse and isfile(ga_fpath):
            info('      Reusing existing Genome Analysis results: ' + ga_fpath)
        else:
            info('      Extracting biosynthetic secondary metabolites to {ga_fpath} (based on hmmsearch results)'.format(**locals()))
            sec_met_name_to_enzymes = get_corresponding_enzymes(cfg, hmmsearch_output_dir, sec_met_fpath, prod_class)
            if not sec_met_name_to_enzymes:
                warning('      No {prod_class} metabolites found in {sequence_fpath}!'.format(**locals()))
                continue
            ga_fpath = write_sec_metabolites_genome_analysis_results(ga_fpath, sec_met_name_to_seq, sec_met_name_to_enzymes)
            if not ga_fpath:  # this should not happen actually, previous if should consider this case too
                warning('      Failed to write {prod_class} metabolites found in {sequence_fpath}!'.format(**locals()))
                continue
        ga_fpaths.append(ga_fpath)
    return ga_fpaths


def process_antismash_file(cfg, sequence_fpath, prod_class):
    ga_fpaths = []

    log_fpath = sequence_fpath + '_antismash_' + config.LOGGER_METAMINER_NAME + '.log'
    if cfg.pipeline.reuse and is_complete_log(log_fpath):
        info('    Reusing existing Genome Analysis results for ' + sequence_fpath)
        for pc in get_prod_classes_to_process(prod_class):
            ga_fpath = get_genome_analysis_fpath(sequence_fpath, pc)
            if not isfile(ga_fpath):
                warning('      No {pc} metabolites found in {sequence_fpath} '
                        '(based on reused results)!'.format(**locals()))
                continue
            ga_fpaths.append(ga_fpath)
        return ga_fpaths

    info('    Extracting biosynthetic secondary metabolites from {sequence_fpath} '
         '(logging to {log_fpath})'.format(**locals()))

    all_cds = get_all_CDS(sequence_fpath)
    if not all_cds:
        warning('    Genome Analysis failed on {sequence_fpath} (no CDS were found)!')
    else:
        ga_fpaths = sort_sec_metabolites_into_production_classes(cfg, all_cds, sequence_fpath, prod_class)

    with open(log_fpath, 'w') as f:
        f.write('#DONE\n')
    return ga_fpaths