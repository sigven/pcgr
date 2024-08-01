GENOME="grch38" # or "grch37"
BUNDLE_VERSION="20240621"
BUNDLE="pcgr_ref_data.${BUNDLE_VERSION}.${GENOME}.tgz"

wget https://insilico.hpc.uio.no/pcgr/${BUNDLE}
gzip -dc ${BUNDLE} | tar xvf -
