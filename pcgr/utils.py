#!/usr/bin/env python
import sys
import subprocess
import logging
import os
import platform

def which(program, env=None):
    """ returns the path to an executable or None if it can't be found"""
    if env is None:
        env = os.environ.copy()

    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in env["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file
    for path in get_all_conda_bins():
        exe_file = os.path.join(path, program)
        if is_exe(exe_file):
            return exe_file
    return None



def pcgr_error_message(message, logger):
    logger.error("")
    logger.error(message)
    logger.error("")
    sys.exit(1)


def pcgr_warn_message(message, logger):
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

def rscript_path(docker_run):
    prefix = conda_env_path('pcgrr', docker_run)
    return os.path.join(prefix, 'bin/Rscript')

def pcgrr_script_path(docker_run):
    prefix = conda_env_path('pcgr', docker_run)
    return os.path.join(prefix, 'bin/pcgrr.R')

def conda_env_path(env, docker_run):
    if docker_run:
        env_path = f'/opt/mambaforge/envs/{env}'
    else:
        cp = os.environ.get('CONDA_PREFIX')
        envdir = os.path.dirname(cp)
        env_path = os.path.join(envdir, env)
    return env_path

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

def get_docker_user_id(docker_user_id):
    logger = getlogger('pcgr-get-OS')
    uid = ''
    if docker_user_id:
        uid = docker_user_id
    elif platform.system() == 'Linux' or platform.system() == 'Darwin' or sys.platform == 'darwin' or sys.platform == 'linux2' or sys.platform == 'linux':
        uid = os.getuid()
    else:
        if platform.system() == 'Windows' or sys.platform == 'win32' or sys.platform == 'cygwin':
            uid = getpass.getuser()

    if uid == '':
        warn_msg = (f'Was not able to get user id/username for logged-in user on the underlying platform '
                    f'(platform.system(): {platform.system()},  sys.platform: {sys.platform}, now running PCGR as root')
        logger.warning(warn_msg)
        uid = 'root'
    return uid

def get_pcgr_bin():
    return os.path.dirname(os.path.realpath(sys.executable))

def perl_cmd():
    """Retrieve path to locally installed conda Perl or first in PATH.
    """
    perl = which(os.path.join(get_pcgr_bin(), "perl"))
    if perl:
        return perl
    else:
        return which("perl")

def get_perl_exports():
    """Environmental exports to use conda installed perl.
    """
    perl_path = os.path.dirname(perl_cmd())
    out = f"unset PERL5LIB && export PATH={perl_path}:\"$PATH\""
    return out
