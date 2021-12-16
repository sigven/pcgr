#!/usr/bin/env bash

#PLATFORM="linux-64"
PLATFORM="osx-64"

mamba create -n pcgr --file env/lock/pcgr/pcgr-${PLATFORM}.lock
mamba create -n pcgrr --file env/lock/pcgrr/pcgrr-${PLATFORM}.lock
