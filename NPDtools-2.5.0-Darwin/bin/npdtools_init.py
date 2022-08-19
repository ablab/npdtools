#!/usr/bin/env python

#####################################################################################
# Copyright (c) 2015-2017 Saint Petersburg State University, St. Petersburg, Russia
# Copyright (c) 2015-2017 University of California San Diego, La Jolla, CA, USA
# All Rights Reserved
# See file LICENSE for details.
#####################################################################################

from os.path import abspath, dirname, realpath, isdir, join

# developers configuration
configs_dir = abspath(dirname(realpath(__file__)))
bin_dir = join(configs_dir, 'bin')
python_modules_dir = None
shell_scripts_dir = None
external_tools_dir = None
docs_dir = None


def init():
    global configs_dir
    global bin_dir
    global python_modules_dir
    global shell_scripts_dir
    global external_tools_dir
    global docs_dir

    # check whether we are in "user" configuration (bin + share dirs) or in "developer" one (default setting)
    install_prefix = dirname(configs_dir)
    possible_bin_dir = join(install_prefix, 'bin')
    possible_configs_dir = join(install_prefix, 'share', 'npdtools')
    if isdir(possible_bin_dir) and isdir(possible_configs_dir):  # "user" configuration
        bin_dir = possible_bin_dir
        configs_dir = possible_configs_dir
    if not isdir(bin_dir):  # ProteoSAFe mode
        bin_dir = configs_dir

    python_modules_dir = join(configs_dir, 'python_libs')
    shell_scripts_dir = join(configs_dir, 'shell_scripts')
    external_tools_dir = join(configs_dir, 'external_tools')
    docs_dir = join(configs_dir, 'docs')
