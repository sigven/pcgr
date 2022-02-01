#!/usr/bin/env bash
set -euo pipefail

if [[ "$OSTYPE" == "darwin"* ]]; then
    PLATFORM="MacOSX"
else
    PLATFORM="Linux"
fi

MINICONDA_URL="https://repo.continuum.io/miniconda/Miniconda3-latest-${PLATFORM}-x86_64.sh"

wget ${MINICONDA_URL} -O miniconda.sh && chmod +x miniconda.sh
# Follow prompts
bash miniconda.sh

# If you want faster conda package installations, use mamba after installing conda:
# conda install -c conda-forge mamba
