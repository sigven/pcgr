#!/usr/bin/env bash
# Install all PCGR dependencies to avoid using the Docker image.
# Suitable for HPC systems lacking Docker, and for debugging.
# Works only on Linux due to unavaliability of few packages for MacOS (ensembl-vep -> perl-bio-db-hts)

set -e
set -o pipefail

THIS_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
SKIP_VALIDATOR=$1

### Setting up conda environment
MINICONDA_DIR=${THIS_DIR}/miniconda
ENV_NAME=pcgr

# Install conda if needed:
if [ ! -x "$(command -v conda)" ]; then
    if [ ! -d ${MINICONDA_DIR} ] ; then
        echo "conda executble not in PATH; install conda under ${MINICONDA_DIR}"
        wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ${THIS_DIR}/miniconda.sh
        # on macos: wget https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh -O ${THIS_DIR}/miniconda.sh
        bash ${THIS_DIR}/miniconda.sh -b -p ${MINICONDA_DIR}
    fi
    . ${MINICONDA_DIR}/etc/profile.d/conda.sh
fi

if [ ! -z ${CONDA_DEFAULT_ENV} ] && [ $(basename ${CONDA_DEFAULT_ENV}) = ${ENV_NAME} ] ; then
    echo "Installing into existing environment ${CONDA_DEFAULT_ENV}"
    if [[ ${CONDA_DEFAULT_ENV} = *"/"* ]]; then
        PARAM="-p"
    else
        PARAM="-n"
    fi
    conda env update ${PARAM} ${CONDA_DEFAULT_ENV} --file ${THIS_DIR}/conda_environment.yml
elif [ -d ${MINICONDA_DIR}/envs/${ENV_NAME} ] ; then
    echo "Activating and upding existing environment ${ENV_NAME}"
    conda activate ${ENV_NAME}
    conda env update ${PARAM} ${CONDA_DEFAULT_ENV} --file ${THIS_DIR}/conda_environment.yml
else
    echo "Creating new environment ${ENV_NAME}"
    conda env create -n ${ENV_NAME} --file ${THIS_DIR}/conda_environment.yml
    conda activate ${ENV_NAME}
fi

# Create a loader. Usage: `source load_pcgr.sh`
cat <<EOT > ${THIS_DIR}/load_pcgr.sh
. ${MINICONDA_DIR}/etc/profile.d/conda.sh
conda activate ${ENV_NAME}
EOT

### Few more R packages that are not on yet conda

# The "options(unzip = )" hack to address the install_github issue under conda https://github.com/r-lib/devtools/issues/1722
export TAR=/bin/tar  # to avoid "/bin/gtar: not found"

R -e "library(devtools); options(unzip = '$(which unzip)'); devtools::install_github('mjkallen/rlogging')"
#R -e "library(devtools); options(unzip = '$(which unzip)'); devtools::install_github('kent37/summarywidget', dependencies=FALSE)"  # added to conda env
R -e "options(unzip = '$(which unzip)'); install.packages('configr', dependencies = T, repos = 'http://cran.us.r-project.org')"
#R -e "install.packages('data.tree', dependencies = T, repos = 'http://cran.us.r-project.org')"  # doesn't work
R -e "library(devtools); options(unzip = '$(which unzip)'); devtools::install_github('AdeelK93/collapsibleTree', dependencies=FALSE)"  # to avoid re-installing conda's Rcpp and others
R -e "library(devtools); options(unzip = '$(which unzip)'); devtools::install_github('Francescojm/CELLector', dependencies=FALSE)"

# ggplot2 is available in package repositories, but dev version is recommended
R -e "library(devtools); options(unzip = '$(which unzip)'); devtools::install_github('hadley/ggplot2', dependencies=FALSE)"
# This one is local
SRC_DIR=${THIS_DIR}/../src
R -e "library(devtools); devtools::install('${SRC_DIR}/R/pcgrr', dependencies=FALSE)"
# Access to src scripts
chmod +x ${SRC_DIR}/pcgr/*.py
chmod +x ${SRC_DIR}/*.R

# Install VEP plugins:
#vep_install --AUTO p --PLUGINS miRNA,LoF --NO_HTSLIB --NO_UPDATE --NO_BIOPERL
# Updating conda-installed loftee with most recent code
#git clone https://github.com/konradjk/loftee ${CONDA_PREFIX}/share/loftee_repo
#rsync -trv --exclude ".git" ${CONDA_PREFIX}/share/loftee_repo/ ${CONDA_PREFIX}/share/ensembl-vep-94.4-0/
#rm -rf ${CONDA_PREFIX}/share/loftee_repo

if [ -z ${SKIP_VALIDATOR} ] ; then
    # Install the EBI vcf validator
    wget https://github.com/EBIvariation/vcf-validator/releases/download/v0.6/vcf_validator -O ${CONDA_PREFIX}/bin/vcf_validator
    chmod +x ${CONDA_PREFIX}/bin/vcf_validator
fi

echo ""
echo "--------------------------------------------------------------------"
echo "PCGR installation complete. To load, source ${THIS_DIR}/load_pcgr.sh"
