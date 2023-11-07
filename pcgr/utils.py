#!/usr/bin/env python
import sys
import subprocess
import shutil
import logging
import os
import time
import errno
import platform
import string
import random


def getlogger(logger_name):
    logger = logging.getLogger(logger_name)
    logger.setLevel(logging.DEBUG)
    # create console handler and set level to debug
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.DEBUG)
    # add ch to logger
    logger.addHandler(ch)
    # create formatter
    formatter = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s", "20%y-%m-%d %H:%M:%S")
    # add formatter to ch
    ch.setFormatter(formatter)
    return logger

def error_message(message, logger):
    logger.error("")
    logger.error(message)
    logger.error("")
    sys.exit(1)

def warn_message(message, logger):
    logger.warning(message)
    
def random_string(length):
    letters = string.ascii_lowercase
    result_str = ''.join(random.choice(letters) for i in range(length))
    return result_str

def check_subprocess(logger, command, debug):
    if debug:
        logger.info(command)
    try:
        output = subprocess.check_output(
            str(command), stderr=subprocess.STDOUT, shell=True)
        if len(output) > 0:
            print(str(output.decode()).rstrip())
    except subprocess.CalledProcessError as e:
        print(e.output.decode())
        exit(0)

def script_path(env, bin_script):
    """Returns e.g. /path/conda/envs/{env}/{bin_script}
    """
    prefix = conda_env_path(env)
    return os.path.join(prefix, bin_script)

def conda_prefix():
    return os.path.normpath(os.environ.get('CONDA_PREFIX'))

def conda_prefix_basename():
    cp = conda_prefix()          # /conda/envs/FOO
    return(os.path.basename(cp)) # FOO

def conda_prefix_dirname():
    cp = conda_prefix()       # /path/to/conda/envs/FOO
    return(os.path.dirname(cp)) # /path/to/conda/envs

def conda_env_path(env):
    """Construct absolute path to a conda env
    using the current activated env as a prefix.
    e.g. /path/to/conda/envs/{env}
    """
    env_dir = conda_prefix_dirname()  # /path/to/conda/envs
    env_path = os.path.join(env_dir, env) # /path/to/conda/envs/{env}
    return env_path

def get_loftee_dir():
    pcgr_conda_env = conda_prefix_basename()
    return script_path(pcgr_conda_env, "share/loftee")

def get_pcgr_bin():
    """Return abs path to e.g. conda/env/pcgr/bin
    """
    return os.path.dirname(os.path.realpath(sys.executable))

def perl_cmd():
    """Retrieve abs path to locally installed conda Perl or first in PATH,
    e.g. conda/env/pcgr/bin/perl
    """
    perl = shutil.which(os.path.join(get_pcgr_bin(), "perl"))
    if perl:
        return perl
    else:
        return shutil.which("perl")

def get_perl_exports():
    """Environmental exports to use conda installed perl.
    """
    perl_path = os.path.normpath(perl_cmd())      # /conda/env/pcgr/bin/perl
    perl_path_parent = os.path.dirname(perl_path) # /conda/env/pcgr/bin
    out = f"unset PERL5LIB && export PATH={perl_path_parent}:\"$PATH\""
    return out

def is_integer(n):
    try:
        float(n)
    except ValueError:
        return False
    else:
        return float(n).is_integer()

def get_cpsr_version():
    # use pcgrr's Rscript to grab cpsr's R pkg version
    # NOTE: this relies on the pcgrr conda env being named the default pcgrr. If
    # you're using a different name (e.g. foo_pcgrr), it returns UNKNOWN.
    rscript = script_path("pcgrr", "bin/Rscript")
    if not os.path.exists(rscript):
        return "UNKNOWN - check your pcgrr conda env is named 'pcgrr'"
    v_cmd = f"{rscript} -e 'x <- paste0(\"cpsr \", as.character(packageVersion(\"cpsr\"))); cat(x, \"\n\")'"
    return subprocess.check_output(v_cmd, shell=True).decode("utf-8")

# https://stackoverflow.com/a/10840586/2169986
def remove(filename):
    try:
        os.remove(filename)
    except OSError as e:
        if e.errno != errno.ENOENT: # errno.ENOENT = no such file or directory
            raise # re-raise exception if a different error occurred


# from https://github.com/bcbio/bcbio-nextgen/blob/master/bcbio/utils.py
def safe_makedir(dname):
    """Make a directory if it does not exist, handling concurrent race conditions.
    """
    if not dname:
        return dname
    num_tries = 0
    max_tries = 5
    while not os.path.exists(dname):
        # we could get an error here if multiple processes are creating
        # the directory at the same time. Grr, concurrency.
        try:
            os.makedirs(dname)
        except OSError:
            if num_tries > max_tries:
                raise
            num_tries += 1
            time.sleep(2)
    return dname

def sort_bed(unsorted_bed_fname: str, sorted_bed_fname: str, debug = False, logger = None):
    
    ## Sort regions in target BED
    if os.path.exists(unsorted_bed_fname) and os.stat(unsorted_bed_fname).st_size != 0:
        cmd_sort_custom_bed1 = 'egrep \'^[0-9]\' ' + str(unsorted_bed_fname) + \
            ' | sort -k1,1n -k2,2n -k3,3n > ' + str(sorted_bed_fname)
        cmd_sort_custom_bed2 = 'egrep -v \'^[0-9]\' ' + str(unsorted_bed_fname) + \
            ' | egrep \'^[XYM]\' | sort -k1,1 -k2,2n -k3,3n >> ' + str(sorted_bed_fname)

        check_subprocess(logger, cmd_sort_custom_bed1, debug)
        check_subprocess(logger, cmd_sort_custom_bed2, debug)
        if not debug:
            remove(str(unsorted_bed_fname))
    else:
        err_msg = 'File ' + str(unsorted_bed_fname) + ' does not exist or is empty'
        error_message(err_msg, logger)


def check_file_exists(fname: str, logger = None) -> bool:
    err = 0
    if not os.path.isfile(fname):
        err = 1
    else:
        if os.stat(fname).st_size == 0:
            err = 1
    if err == 1:
        err_msg = f"File {fname} does not exist or has zero size"
        error_message(err_msg, logger)
    
    return(True)

def check_tabix_file(fname: str, logger = None) -> bool:
    tabix_file = os.path.abspath(f'{fname}.tbi')
    if not os.path.exists(tabix_file):
        err_msg = f"Tabix file (i.e. '.gz.tbi') is not present for the bgzipped VCF input file ({fname}" + \
            "). Please make sure your input VCF is properly compressed and indexed (bgzip + tabix)"
        error_message(err_msg, logger)
    else:
        ## check file size is more than zero
        check_file_exists(tabix_file)
    return(True)
