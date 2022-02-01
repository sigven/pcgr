#!/usr/bin/env bash
set -euo pipefail

if [[ "$OSTYPE" == "darwin"* ]]; then
    PLATFORM="osx-64"
else
    PLATFORM="linux-64"
fi

PCGR_CONDA_DIR="../conda"

mamba create --prefix ${PCGR_CONDA_DIR}/env/pcgr --file ${PCGR_CONDA_DIR}/env/lock/pcgr-${PLATFORM}.lock
mamba create --prefix ${PCGR_CONDA_DIR}/env/pcgrr --file ${PCGR_CONDA_DIR}/env/lock/pcgrr-${PLATFORM}.lock

# activate pcgr conda environment
conda activate ${PCGR_CONDA_DIR}/env/pcgr
