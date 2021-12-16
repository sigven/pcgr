#!/usr/bin/env bash

#PLATFORM="linux-64"
PLATFORM="osx-64"

mamba create -n pcgr --file env/lock/pcgr-${PLATFORM}.lock
mamba create -n pcgrr --file env/lock/pcgrr-${PLATFORM}.lock
