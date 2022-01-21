#!/usr/bin/env bash

#PLATFORM="linux-64"
PLATFORM="osx-64"

mamba create --prefix ./env/pcgr --file env/lock/pcgr-${PLATFORM}.lock
mamba create --prefix ./env/pcgrr --file env/lock/pcgrr-${PLATFORM}.lock
