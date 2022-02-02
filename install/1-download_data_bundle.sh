#!/usr/bin/env bash
set -euo pipefail

#GENOME="grch37"
GENOME="grch38"
BUNDLE="pcgr.databundle.${GENOME}.20220119.tgz"

wget http://insilico.hpc.uio.no/pcgr/${BUNDLE}
gzip -dc ${BUNDLE} | tar xvf -
