#!/usr/bin/env bash
set -euo pipefail

PCGR_VERSION="0.10.12"
OUTPUT_DIRECTORY="PCGR"

git clone \
    -b v${PCGR_VERSION} \
    --depth 1 \
    https://github.com/sigven/pcgr.git \
    ${OUTPUT_DIRECTORY}
