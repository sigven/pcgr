#!/usr/bin/env bash
set -x

# For the reference: CONDA env vars: https://conda.io/docs/user-guide/tasks/build-packages/environment-variables.html

# Changing permissions to executables
chmod +x ${SRC_DIR}/src/pcgr/*.py
chmod +x ${SRC_DIR}/*.py
chmod +x ${SRC_DIR}/src/*.R
chmod +x ${SRC_DIR}/src/R/*.tar.gz
# Moving libraries and scripts
mkdir -p ${PREFIX}/lib/R/library
mkdir -p ${PREFIX}/bin
mkdir -p ${SP_DIR}
#mv ${SRC_DIR}/src/R/ $PREFIX/lib/R/library/
mv ${SRC_DIR}/src/pcgr/lib/* ${SP_DIR}/  # python modules
mv ${SRC_DIR}/src/pcgr/*.py ${SRC_DIR}/*.py ${PREFIX}/bin/  # python scripts
mv ${SRC_DIR}/src/*.R ${PREFIX}/bin/  # R scripts

# VCF validator
wget https://github.com/EBIvariation/vcf-validator/releases/download/v0.9.3/vcf_validator_linux -O ${PREFIX}/bin/vcf_validator
chmod +x ${PREFIX}/bin/vcf_validator

## Install assertable package (not found through Conda channels), and pcgrr
R -e "install.packages('${SRC_DIR}/src/R/assertable_0.2.7.tar.gz', lib='${PREFIX}/lib/R/library', repos = NULL)"
R -e "install.packages('${SRC_DIR}/src/R/pcgrr_0.9.1.tar.gz', lib='${PREFIX}/lib/R/library', repos = NULL)"

### Loftee. To make sure same LoF version is used in dockerized and non-dockerized installation.
#   ensembl-vep conda package installs most recent version of LoF automatically, however it doesn't work with the most
#   recent perl 5.26 (see https://github.com/sigven/cpsr/issues/2)
#   Also Loftee for hg38 needs Perl-Bio-BigFile (doesn't come with ensemble-vep, but will come
#   with https://github.com/bioconda/bioconda-recipes/pull/18808 once merged)
mkdir ${PREFIX}/share/loftee
tar -xzf ${SRC_DIR}/src/loftee_1.0.3.tgz -C ${PREFIX}/share/loftee
