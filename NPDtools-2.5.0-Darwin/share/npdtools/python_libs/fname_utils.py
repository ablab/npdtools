import os
from os.path import basename, splitext, abspath, join
from log_utils import error

try:
    from collections import OrderedDict
except ImportError:
    from ordered_dict import OrderedDict

internal_fpath_mapping = OrderedDict()


def _get_internal_fpath(fpath):
    try:
        return internal_fpath_mapping[abspath(fpath)]
    except KeyError:
        error('Internal error -- internal mapping for {fpath} is not found'.format(**locals()))


def _get_internal_name(fpath, is_internal_fpath=False):
    if is_internal_fpath:
        internal_fpath = fpath
    else:
        internal_fpath = _get_internal_fpath(fpath)
    return splitext(basename(internal_fpath))[0]


def generate_internal_fpath_mapping(cfg, fpaths, num_seq=0):
    from spectra_utils import spectra_final_ext
    global internal_fpath_mapping
    for idx, fpath in enumerate(fpaths):
        if idx < num_seq:
            prefix = 'sequence-' + str(idx).zfill(5)
            ext = _get_ext(fpath)
        else:
            prefix = 'spectra-' + str(idx - num_seq).zfill(5)
            ext = spectra_final_ext(fpath)
        internal_fpath_mapping[abspath(fpath)] = join(cfg.work_dirpath, prefix + ext)
    return internal_fpath_mapping


def update_fnames(cfg, file_mapping):
    for report_fname in [_get_matches_fpath(cfg), _get_matches_fpath(cfg, 'sig')]:
        with open(report_fname) as f:
            header = f.readline()
            spec_fname_index = header.split('\t').index('SpecFile')
            if 'SeqFile' in header:
                seq_fname_index = header.split('\t').index('SeqFile')
            else:
                seq_fname_index = None
            content = f.readlines()
        with open(report_fname, 'w') as f:
            f.write(header)
            for line in content:
                entries = line.split('\t')
                if basename(entries[spec_fname_index]) in file_mapping:
                    entries[spec_fname_index] = file_mapping[basename(entries[spec_fname_index])]
                if seq_fname_index is not None and basename(entries[seq_fname_index]) in file_mapping:
                    entries[seq_fname_index] = file_mapping[basename(entries[seq_fname_index])]
                f.write('\t'.join(entries))


def _get_fpath_by_idx(cfg, idx, ext='.mgf'):
    prefix = 'spectra-' + str(idx).zfill(5)
    return join(cfg.work_dirpath, prefix + ext)


def _get_prefix_by_fpath(fpath):
    return splitext(basename(fpath))[0]


def _get_matches_fpath(cfg, type='all', decoy=False, unique=False):
    return join(cfg.output_dirpath, '%s_%s%smatches.tsv' %
                                    ('all' if type == 'all' else 'significant',
                                     'unique_' if unique else '',
                                     'decoy_' if decoy else ''))


def _get_local_matches_fpath(fpath, type='all', decoy=False):
    return fpath + '.%s.%smatches' % ('all' if type == 'all' else 'sig',
                                      'decoy.' if decoy else '')


def _get_best_spectra_fpath(cfg):
    return join(cfg.output_dirpath, 'best_spectra.mgf')


def _get_params_fpath(cfg):
    return join(cfg.work_dirpath, 'params.tsv')


def _get_corresp_fpath(cfg):
    return join(cfg.work_dirpath, 'corresp.tsv')


def _get_summary_fpath(cfg):
    return join(cfg.output_dirpath, 'summary.tsv')


def _get_mapping_fpath(cfg):
    return join(cfg.work_dirpath, 'mapping.tsv')


def _get_decoy_fpath(sequence_fpath):
    ext = splitext(sequence_fpath)[1]
    return sequence_fpath[:len(sequence_fpath) - len(ext)] + '.decoy' + ext


def _get_ext(fpath):
    return splitext(fpath)[1]


def _set_ext(fname, ext):
    if splitext(fname)[1] != ext:
        return fname + ext
    return fname


def compare_basename_wo_extension(path1, path2, extensions):
    def __simplify_path(path):
        base, ext = os.path.splitext(os.path.basename(path))
        if ext in extensions:
            return base
        return base + ext

    return __simplify_path(path1) == __simplify_path(path2)