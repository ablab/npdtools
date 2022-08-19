from os.path import join, abspath, basename, isfile, isdir, islink
from common import info, warning, check_file, sys_call, get_mode_command, get_param_command, clean
from config import prod_classes, excluded_prod_class
from fname_utils import _get_matches_fpath, _get_params_fpath, _get_internal_name, \
    _get_internal_fpath, _get_decoy_fpath
from file_utils import remove_if_exists, is_complete_log
from fasta_utils import shuffle_fasta, is_nucleotide_fasta
from spectra_utils import spectra_need_convertion, convert_to_mgf
import os
import config
import common

rippquest_binary_fpath = None
configs_dir = None
external_tools_dir = None


class PSM_ORF(object):
    header_line = "SpecFile\tScan\tSeqFile\tClass\tFragmentSeq\tModifiedSeq\tScore\tP-Value\tPeptideMass\tSpectrumMass\tRetention\tCharge\n"

    def __init__(self, spectrum_fpath, sequence_fpath, prod_class, line, decoy=False):
        self.spectrum_fpath = spectrum_fpath
        self.sequence_fpath = sequence_fpath
        self.prod_class = prod_class

        #CPP OUTPUT:  orf_idx fragment_idx orf_name fragment_seq nlp_mass scan_id spec_mz
        #             charge rt mass_error score pvalue modified_seq
        #EXAMPLE:     290 24 all_0_2224960_2254944_291_10740_10904_0_0
        #             SSVVIVNPAAFAIATTSAPSALSCSCSEGSA 2897.38 0 2840.25
        #             1 -1 0.00001 2 0.00642013 S-17SVVIVNPAAFAIATT-18S-18APS-18ALS+16CS+16CS-18EGSA
        self.orf_idx, self.fragment_idx, self.orf_name, self.fragment_seq, self.nlp_mass, \
        self.scan_id, self.spec_mz, self.charge, self.rt, self.mass_error, self.score, self.p_value, \
        self.modified_seq = line.split()
        self.p_value = "%.1e" % float(self.p_value)
        self.decoy = decoy

    def __str__(self):
        # see header_line above
        return "\t".join([self.spectrum_fpath, self.scan_id, self.sequence_fpath, self.prod_class,
                          self.fragment_seq, self.modified_seq,
                          self.score, self.p_value, self.nlp_mass, self.spec_mz, self.rt, self.charge])


def get_prod_class_mod_fpath(prod_class):
    return join(configs_dir, 'configs', 'rippquest', "mod_%s.txt" % prod_class)


def get_genome_analysis_fpath(sequence_fpath, prod_class):
    return sequence_fpath + '.format.' + prod_class + '.intervals.fasta.orf.format.edit'


def get_prod_classes_to_process(prod_class):
    if prod_class == 'all':
        return [pc for pc in prod_classes if pc not in excluded_prod_class]
    return [prod_class]


def process_genome_analysis_file(cfg, spec_fpath, ga_fpath, prod_class, output_fpath):
    # ./rippquest_ms $spec_file $seq_path/all.fa.format.${compound_type}.intervals.fasta.orf.format.edit
    #                $experiment_name $mod_file -max_charge $charge $extra_flag
    experiment_name = cfg.params['mode'].lower()
    prod_class_mod_fpath = get_prod_class_mod_fpath(prod_class)
    command = rippquest_binary_fpath + ' "{spec_fpath}" "{ga_fpath}" "{experiment_name}"' \
                                       ' "{prod_class_mod_fpath}" '.format(**locals()) + \
                                       ' --configs_dir ' + configs_dir
    command += get_mode_command(cfg.params)
    command += get_param_command(cfg.params, 'max_charge')
    command += get_param_command(cfg.pipeline.__dict__, 'antismash', '--external_orf')
    command += get_param_command(cfg.pipeline.__dict__, 'boa', '--external_orf')
    command += get_param_command(cfg.params, 'blind', '--single_blind')
    command += get_param_command(cfg.params, 'debug')
    if cfg.pipeline.pass_to_rippquest_ms:
        command += ' ' + cfg.pipeline.pass_to_rippquest_ms

    command += ' > "{output_fpath}"'.format(**locals())
    if not sys_call(command, output_fpath, indent='    '):
        warning('    RiPPquest-MS failed on {ga_fpath}!'.format(**locals()))
        return None
    return output_fpath


def process_one_sequence_file(cfg, fpath, prod_class):
    def __get_sequence_type(fpath):
        ext = os.path.splitext(fpath)[1].lower()
        raw_sequence_type = None
        for sequence_type, extensions in config.sequence_extensions_dict.items():
            if ext in extensions:
                raw_sequence_type = sequence_type
                break
        if raw_sequence_type is None:
            warning('   Failed to determine sequence file type for %s, assuming nucleotide fasta!' % fpath)
            return 'nucleotide_fasta'
        if raw_sequence_type == 'fasta':
            is_nucl = is_nucleotide_fasta(fpath)
            if is_nucl is None:
                warning('   Failed to determine whether FASTA file contains nucleotides or amino acids for %s, assuming nucleotide fasta!' % fpath)
                return 'nucleotide_fasta'
            if is_nucl:
                return 'nucleotide_fasta'
            else:
                return 'aminoacid_fasta'
        return raw_sequence_type

    sequence_fpath = _get_internal_fpath(fpath)
    info('Starting: ' + fpath + ' --> ' + _get_internal_name(sequence_fpath, is_internal_fpath=True))
    if isfile(sequence_fpath) and not abspath(sequence_fpath) == abspath(fpath):
        os.remove(sequence_fpath)
    if not isfile(sequence_fpath):
        os.symlink(abspath(fpath), sequence_fpath)

    sequence_type = __get_sequence_type(sequence_fpath)
    info('\t\t %s would be processed as %s file' % (_get_internal_name(sequence_fpath, is_internal_fpath=True), sequence_type))
    from r_utils import ga_fasta, ga_antismash, boa_utils
    if sequence_type == 'nucleotide_fasta':
        ga_fpaths = ga_fasta.process_regular_fasta_file(cfg, sequence_fpath, prod_class)
    elif sequence_type == 'aminoacid_fasta':
        ga_fpaths = ga_fasta.process_aminoacid_fasta_file(cfg, sequence_fpath, prod_class)
    elif sequence_type == 'antismash':
        ga_fpaths = ga_antismash.process_antismash_file(cfg, sequence_fpath, prod_class)
    elif sequence_type == 'boa':
        converted_sequence_fpath = boa_utils.boa2fasta(sequence_fpath)
        if converted_sequence_fpath is None:
            ga_fpaths = []
        else:
            ga_fpaths = ga_fasta.process_aminoacid_fasta_file(cfg, sequence_fpath, prod_class,
                                                              actual_in_fpath=converted_sequence_fpath)
    else:
        assert False, "incorrect sequence_type!"

    clean(cfg, basename(sequence_fpath), intermediate=True)  # we needed to keep Genome Analysis files for a while
    if not ga_fpaths:
        return False
    for ga_fpath in ga_fpaths:
        if cfg.fdr:
            decoy_fpath = _get_decoy_fpath(ga_fpath)
            if cfg.pipeline.reuse and check_file(decoy_fpath):
                info('    Reusing existing Decoy Genome Analysis file: ' + decoy_fpath)
            else:
                remove_if_exists(decoy_fpath)
                info('    Creating Decoy Genome Analysis file: ' + decoy_fpath)
                shuffle_fasta(ga_fpath, decoy_fpath)
    return True


def process_one_spectra_file(cfg, fpath, sequence_fpaths):
    info('Starting: ' + fpath + ' --> ' + _get_internal_name(fpath))
    empty_output = ([], []), False

    prepared_fpath = _get_internal_fpath(fpath)
    if cfg.pipeline.reuse and check_file(prepared_fpath):
        info('    Reusing existing spectra file ({prepared_fpath})'.format(**locals()))
    else:
        if isfile(prepared_fpath):
            os.remove(prepared_fpath)
        if spectra_need_convertion(fpath):
            if not convert_to_mgf(external_tools_dir, fpath, prepared_fpath):
                warning('Conversion of input spectra for %s failed, skipping this file' %
                        (basename(fpath)))
                return empty_output
        else:
            os.symlink(abspath(fpath), prepared_fpath)

    PSMs = []
    sig_PSMs = []
    at_least_one_is_ok = False
    for prod_class in get_prod_classes_to_process(cfg.params['prod_class']):
        genome_analysis_output = list()  # list of pairs (seq_fpath, ga_fpath, is_decoy)
        for seq_fpath in sequence_fpaths:
            ga_fpath = get_genome_analysis_fpath(_get_internal_fpath(seq_fpath), prod_class)
            if not isfile(ga_fpath) and not islink(ga_fpath):
                continue
            genome_analysis_output.append((seq_fpath, ga_fpath, False))
            if cfg.fdr and check_file(_get_decoy_fpath(ga_fpath)):
                genome_analysis_output.append((seq_fpath, _get_decoy_fpath(ga_fpath), True))

        if not genome_analysis_output:
            warning('    No Genome Analysis files for {fpath} ({prod_class})! Skipping..'.format(**locals()))
        else:
            for seq_fpath, ga_fpath, decoy in genome_analysis_output:
                log_fpath = join(cfg.work_dirpath, _get_internal_name(fpath) + '_' +
                                  _get_internal_name(seq_fpath) + '_' + prod_class + '_' + config.LOGGER_METAMINER_NAME + '.log')
                if decoy:
                    log_fpath = _get_decoy_fpath(log_fpath)
                if cfg.pipeline.reuse and is_complete_log(log_fpath):
                    info('    Reusing existing RiPPquest-MS results: ' + log_fpath)
                else:
                    log_fpath = process_genome_analysis_file(cfg, prepared_fpath, ga_fpath, prod_class, log_fpath)
                if log_fpath is not None:
                    cur_PSMs, cur_sig_PSMs = parse_output(fpath, seq_fpath, prod_class, log_fpath, decoy=decoy)
                    PSMs += cur_PSMs
                    sig_PSMs += cur_sig_PSMs
                    at_least_one_is_ok = True
    clean(cfg, prepared_fpath, intermediate=True)  # we needed to keep Genome Analysis files for a while
    return (PSMs, sig_PSMs), at_least_one_is_ok


def parse_output(spectrum_fpath, sequence_fpath, prod_class, output_fpath, decoy=False):
    PSMs = []
    significant_PSMs = []
    with open(output_fpath) as f:
        for line in f:
            if line.startswith("#"):
                continue
            if line.startswith('.'):
                warning('Skipping "%s" from the RiPPquest-MS output (starts with "." -- does not seem to be a PSM)!' % line.strip())
                continue
            if len(line.split()) < len(PSM_ORF.header_line.split('\t')) - 4:  # 4 PSM fields are not read from the line
                warning('Skipping "%s" from the RiPPquest-MS output (too few columns -- does not seem to be a PSM)!' % line.strip())
                continue
            PSMs.append(PSM_ORF(spectrum_fpath, sequence_fpath, prod_class, line, decoy=decoy))
            if float(PSMs[-1].p_value) < config.rippquest_p_value_threshold:
                significant_PSMs.append(PSMs[-1])
    return PSMs, significant_PSMs


def save_results(cfg, PSMs, significant_PSMs, first_time=False, decoy=False):
    num_all = 0
    with open(_get_matches_fpath(cfg, decoy=decoy), 'w' if first_time else 'a') as all:
        if first_time:
            all.write(PSM_ORF.header_line)
        for psm in PSMs:
            if psm.decoy == decoy:
                all.write(str(psm) + "\n")
                num_all += 1
    # significant matches
    num_sig = 0
    with open(_get_matches_fpath(cfg, 'sig', decoy=decoy), 'w' if first_time else 'a') as sig:
        if first_time:
            sig.write(PSM_ORF.header_line)
        for psm in significant_PSMs:
            if psm.decoy == decoy:
                sig.write(str(psm) + "\n")
                num_sig += 1
    return num_all, num_sig


def save_params(cfg):
    params_fpath = _get_params_fpath(cfg)
    with open(params_fpath, 'w') as f:
        f.write('Parameter\tValue\n')

        mode, (pm_thresh, product_ion_thresh) = common.get_mode(cfg.params)
        f.write('Mode\t{mode}\n'.format(**locals()))
        f.write('Parent Mass Tolerance\t{pm_thresh} Da\n'.format(**locals()))
        f.write('Product Ion Tolerance\t{product_ion_thresh} Da\n'.format(**locals()))

        f.write('Product Class\t' + common.get_param(cfg.params, 'prod_class') + '\n')
        f.write('Max charge\t' + common.get_param(cfg.params, 'max_charge') + '\n')
    return params_fpath


def get_sequence_fpaths(all_fpaths):
    def __is_seq_file(fpath):
        if not isfile(fpath):
            return False
        ext = os.path.splitext(fpath)[1].lower()
        if ext in config.allowed_sequence_extensions:
            return True
        return False

    sequence_fpaths = []
    if all_fpaths:
        for entry in all_fpaths:
            if isdir(entry):
                sequence_fpaths += [join(path, fpath) for (path, dirs, files) in os.walk(entry)
                                    for fpath in files if __is_seq_file(join(path, fpath))]
            elif __is_seq_file(entry):
                sequence_fpaths.append(entry)
    sequence_fpaths = list(set(sequence_fpaths))  # taking only unique sequences
    sequence_fpaths.sort()
    return sequence_fpaths