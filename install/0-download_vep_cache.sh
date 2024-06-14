GENOME="GRCh38" # or "GRCh37"
VEP_VERSION="112"
CACHE="homo_sapiens_vep_${VEP_VERSION}_${GENOME}.tar.gz"

wget https://ftp.ensembl.org/pub/release-${VEP_VERSION}/variation/indexed_vep_cache/${CACHE}
gzip -dc ${CACHE} | tar xvf -
