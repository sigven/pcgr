#!/usr/bin/env bash
set -euo pipefail

if [[ "$OSTYPE" == "darwin"* ]]; then
    PLATFORM="osx-64"
else
    PLATFORM="linux-64"
fi

PCGR_CONDA_ENV_DIR="../conda/env"

mamba create --prefix ${PCGR_CONDA_ENV_DIR}/pcgr --file ${PCGR_CONDA_ENV_DIR}/lock/pcgr-${PLATFORM}.lock
mamba create --prefix ${PCGR_CONDA_ENV_DIR}/pcgrr --file ${PCGR_CONDA_ENV_DIR}/lock/pcgrr-${PLATFORM}.lock

# activate pcgr conda environment
conda activate ${PCGR_CONDA_ENV_DIR}/pcgr
