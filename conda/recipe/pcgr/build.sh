#!/usr/bin/env bash
set -x

# Conda env vars: https://docs.conda.io/projects/conda-build/en/latest/user-guide/environment-variables.html

# VCF validator
mkdir -p ${PREFIX}/bin
# releases available for mac and linux
PLATFORM_INFERRED=$(uname -s)
[[ $PLATFORM_INFERRED == "Darwin" ]] && VCF_VALIDATOR_SUFFIX="macos" || VCF_VALIDATOR_SUFFIX="linux"
VCF_VALIDATOR_RELEASE="https://github.com/EBIvariation/vcf-validator/releases/download/v0.9.3/vcf_validator_${VCF_VALIDATOR_SUFFIX}"
curl -L $VCF_VALIDATOR_RELEASE -o ${PREFIX}/bin/vcf_validator

chmod +x ${PREFIX}/bin/vcf_validator

### Loftee. To make sure same LoF version is used in dockerized and non-dockerized installation.
#   ensembl-vep conda package installs most recent version of LoF automatically, however it doesn't work with the most
#   recent perl 5.26 (see https://github.com/sigven/cpsr/issues/2)
#   Also Loftee for hg38 needs Perl-Bio-BigFile (doesn't come with ensemble-vep, but will come
#   with https://github.com/bioconda/bioconda-recipes/pull/18808 once merged)
mkdir ${PREFIX}/share/loftee
tar -xzf ${SRC_DIR}/src/loftee_1.0.3.tgz -C ${PREFIX}/share/loftee

$PYTHON setup.py install #--single-version-externally-managed --root=/
#chmod -R o+r $PREFIX/lib/python*/site-packages/*
