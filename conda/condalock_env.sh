#!/usr/bin/env bash

conda-lock --mamba -f env/yml/pcgr.yml -p osx-64 -p linux-64 --filename-template 'pcgr-{platform}.lock'
conda-lock --mamba -f env/yml/pcgrr.yml -p osx-64 -p linux-64 --filename-template 'pcgrr-{platform}.lock'
