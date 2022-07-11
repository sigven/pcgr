#!/usr/bin/env python
import sys
import subprocess
import shutil
import logging
import os
import errno
import platform


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

def conda_env_path(env):
    """Construct absolute path to a conda env
    using the current activated env as a prefix.
    e.g. /path/to/conda/envs/{env}
    """
    cp = os.path.normpath(os.environ.get('CONDA_PREFIX'))  # /path/to/conda/envs/FOO
    env_dir = os.path.dirname(cp)                          # /path/to/conda/envs
    env_path = os.path.join(env_dir, env)                  # /path/to/conda/envs/{env}
    return env_path

def get_loftee_dir():
    return script_path("pcgr", "share/loftee")

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
    rscript = script_path("pcgrr", "bin/Rscript")
    v_cmd = f"{rscript} -e 'x <- paste0(\"cpsr \", as.character(packageVersion(\"cpsr\"))); cat(x, \"\n\")'"
    return subprocess.check_output(v_cmd, shell=True).decode("utf-8")

# https://stackoverflow.com/a/10840586/2169986
def remove(filename):
    try:
        os.remove(filename)
    except OSError as e:
        if e.errno != errno.ENOENT: # errno.ENOENT = no such file or directory
            raise # re-raise exception if a different error occurred
