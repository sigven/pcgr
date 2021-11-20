#!/usr/bin/env python
import sys
import subprocess
import logging


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
