import config
from log_utils import info, error, warning
from fname_utils import _get_matches_fpath, _get_best_spectra_fpath
from spectra_utils import get_scan_by_id


def extract_best(cfg, external_tools_dir):
    info('Extracting spectra with best hits (p-value < %.1e, FDR is None or < %.2f)' %
         (config.BS_MAX_P_VALUE, config.BS_MAX_FDR))
    best_hits = []
    with open(_get_matches_fpath(cfg, type='sig')) as f:
        header = f.readline().strip().split('\t')
        spec_fpath_idx = header.index('SpecFile')
        scan_idx = header.index('Scan')
        p_value_idx = header.index('P-Value')
        fdr_idx = header.index('FDR') if 'FDR' in header else None
        for line in f:
            values = line.strip().split('\t')
            if float(values[p_value_idx]) < config.BS_MAX_P_VALUE and \
                    (fdr_idx is None or float(values[fdr_idx]) < config.BS_MAX_FDR):
                best_hits.append((values[spec_fpath_idx], values[scan_idx]))
    with open(_get_best_spectra_fpath(cfg), 'w') as f:
        new_scan_id = 1
        for (spec_fpath, scan_id) in best_hits:
            scan = get_scan_by_id(spec_fpath, scan_id, external_tools_dir, cfg.work_dirpath)
            if scan is None:
                warning('Failed to extract scan #%s from %s' % (scan_id, spec_fpath))
            else:
                for scan_line in scan:
                    if scan_line.startswith('SCAN'):
                        f.write('SCANS=%d\n' % new_scan_id)
                    elif scan_line.startswith('TITLE') and 'SCAN' in scan_line:
                        f.write('TITLE=SCAN=%d\n' % new_scan_id)
                    else:
                        f.write(scan_line)
                new_scan_id += 1