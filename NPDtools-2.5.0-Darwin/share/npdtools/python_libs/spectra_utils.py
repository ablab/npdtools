#!/usr/bin/env python

import time
import os
from os.path import join, isfile, splitext, getsize, dirname, basename, isdir
import subprocess
import shlex
from common import info, sys_call, check_file
import config
from fname_utils import _get_ext

not_converted_extensions = ['.mgf', '.mzXML', '.mzdata']
extension_after_conversion = '.mgf'


def get_spectra_fpaths(all_fpaths):
    spectra_fpaths = []
    for entry in all_fpaths:
        if isdir(entry):
            spectra_fpaths += [join(path, fpath) for (path, dirs, files) in os.walk(entry)
                               for fpath in files if is_spectrum_file(join(path, fpath))]
        elif is_spectrum_file(entry):
            spectra_fpaths.append(entry)
    spectra_fpaths.sort()
    return spectra_fpaths


def is_spectrum_file(fpath):
    allowed_extensions = map(lambda x: x.lower(), config.allowed_spectrum_extensions)

    if not isfile(fpath):
        return False
    basename, ext = splitext(fpath)
    if ext.lower() in allowed_extensions:
        return True
    return False


def spectra_need_convertion(fpath):
    basename, ext = splitext(fpath)
    if ext.lower() in map(lambda x: x.lower(), not_converted_extensions):
        return False
    return True


def spectra_final_ext(fpath):
    ext = _get_ext(fpath)
    lowercase_not_converted_extensions = list(map(lambda x: x.lower(), not_converted_extensions))
    if ext.lower() in lowercase_not_converted_extensions:
        return not_converted_extensions[lowercase_not_converted_extensions.index(ext.lower())]
    return extension_after_conversion


def convert_to_mgf(external_tools_dir, input_fpath, mgf_fpath):
    def __need_fixing(input_fpath):
        if splitext(input_fpath)[1] != '.mzXML':
            return False
        with open(input_fpath) as f:
            for line in f:
                if 'nameValue' in line:
                    return True
            return False

    def __fix(input_fpath):
        fixed_fpath = mgf_fpath + '.fix.mzXML'
        sys_call('grep -v nameValue {input_fpath} > {fixed_fpath}'.format(**locals()), silent=False)
        return fixed_fpath

    def __is_running(fpath, size):
        if isfile(fpath) and size < getsize(fpath):
            return True
        return False

    def __safe_getsize(fpath):
        return getsize(fpath) if isfile(fpath) else 0

    # main part
    if __need_fixing(input_fpath):
        input_fpath = __fix(input_fpath)

    out_dirname = dirname(mgf_fpath)
    out_basename = basename(mgf_fpath)
    msconvert_fpath = join(external_tools_dir, 'pwiz', config.platform_name, 'msconvert')
    command = ' {msconvert_fpath} "{input_fpath}" --mgf --filter "titleMaker SCAN=<ScanNumber> MSLEVEL=<MsLevel>" ' \
              '-o "{out_dirname}" --outfile "{out_basename}" '.format(**locals())

    new_env = dict(os.environ)
    new_env['LC_ALL'] = 'C'
    init_timeout = 3
    tries = 4
    with open(os.devnull, 'w') as devnull:
        for i in range(tries):
            if isfile(mgf_fpath):
                os.remove(mgf_fpath)
            try_msg = ''
            if i != 0:
                try_msg = ' (try: %d)' % (i + 1)
            info("  Converting input spectra{try_msg}: {command}".format(**locals()))
            proc = subprocess.Popen(shlex.split(command), stdout=devnull, env=new_env)
            mgf_size = 0
            killed = False
            while proc.poll() is None:
                timeout = init_timeout
                for i in range(tries):
                    time.sleep(timeout)
                    if __is_running(mgf_fpath, mgf_size):
                        break
                    timeout *= 2
                    if not isfile(mgf_fpath) and proc.poll() is not None:
                        return False  # error on conversion
                else:
                    pid = proc.pid
                    proc.terminate()
                    try:
                        time.sleep(init_timeout)
                        os.kill(pid, 0)
                        proc.kill()
                    except OSError:
                        pass  # terminated gracefully by "terminate()"
                    killed = True
                    break
                mgf_size = __safe_getsize(mgf_fpath)
            if not killed:
                return __safe_getsize(mgf_fpath) != 0
        return False


def get_scans_from_mgf(fpath):
    scan = ""
    with open(fpath) as f:
        for line in f:
            if line.startswith('BEGIN IONS'):
                if scan:
                    yield scan
                scan = ""
            scan += line
    if scan:
        yield scan


def get_max_scan_len(mgf_fpath):
    max_len = 0
    for scan in get_scans_from_mgf(mgf_fpath):
        len = scan.count('\n')
        if len > max_len:
            max_len = len
    return max_len


def get_max_charge(mgf_fpath):
    MAX_POSSIBLE_CHARGE = 3
    max_charge = 0
    with open(mgf_fpath) as f:
        for line in f:
            if line.startswith('CHARGE='):
                if line.split('CHARGE=')[1].strip().isdigit():
                    max_charge = max(max_charge, int(line.split('=')[1].strip()))
                    if max_charge == MAX_POSSIBLE_CHARGE:
                        return str(max_charge)
    return str(max_charge) if max_charge else None


def get_scan_by_id(fpath, scan_id, external_tools_dir, work_dir):
    if not check_file(fpath):
        return None
    if splitext(fpath)[1] != '.mgf':
        mgf_fpath = join(work_dir, basename(fpath) + '.mgf')
        if not check_file(mgf_fpath):
            if not convert_to_mgf(external_tools_dir, fpath, mgf_fpath):
                return None
        fpath = mgf_fpath
    with open(fpath) as f:
        cur_scan = []
        cur_scan_id = None
        for line in f:
            if line.strip() == 'BEGIN IONS':
                if cur_scan_id == scan_id:
                    return cur_scan
                cur_scan = []
                if cur_scan_id is None:
                    cur_scan_id = '1'
                elif cur_scan_id.isdigit():
                    cur_scan_id = str(int(cur_scan_id) + 1)
            cur_scan.append(line)
            if 'SCAN' in line:
                # SCAN=5
                # SCANS=5
                # TITLE=SCAN=5 MSLEVEL=2
                cur_scan_id = line[line.find('SCAN') + 1:].split('=')[1].split()[0]
        if cur_scan_id == scan_id:
            return cur_scan
    return None