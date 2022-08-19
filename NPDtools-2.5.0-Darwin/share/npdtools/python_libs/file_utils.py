import os
import errno
from os.path import isdir, isfile, getsize
from log_utils import error


def verify_file(fpath, description=''):
    if not check_file(fpath):
        error(fpath + (' (' + description + ' file)' if description else '') + ' is empty or does not exist!')


def verify_dir(dirpath, description=''):
    if not isdir(dirpath):
        error(dirpath + (' (' + description + ' dir)' if description else '') + ' does not exist!')


def check_file(fpath):
    return isfile(fpath) and getsize(fpath) > 0


def remove_if_exists(fpath):
    # http://stackoverflow.com/questions/10840533/most-pythonic-way-to-delete-a-file-which-may-not-exist
    try:
        os.remove(fpath)
    except OSError as e:
        if e.errno != errno.ENOENT:  # errno.ENOENT = no such file or directory
            raise  # re-raise exception if a different error occured


def count_lines(fpath):
    with open(fpath) as f:
        return len(f.readlines())


def is_complete_log(dereplication_result_fpath):  # our tools' logs specifics: should end with DONE
    if not isfile(dereplication_result_fpath) or getsize(dereplication_result_fpath) <= 0:
        return False
    with open(dereplication_result_fpath) as f:
        return f.read().rstrip().endswith("DONE")