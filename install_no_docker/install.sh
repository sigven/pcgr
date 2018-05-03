# Install all PCGR dependencies to avoid using the Docker image.
# Suitable for HPC systems lacking Docker, and for debugging.
# Works only for Linux due to unavaliability of few packages for MacOS (ensembl-vep -> perl-bio-db-hts)

SRC_DIR=../src

# Install conda
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
bash miniconda.sh -b -p ./miniconda
export PATH=$(pwd)/miniconda/bin:$PATH

# Make pcgr environment
conda env create --file conda_environment.yml
source activate pcgr

# Few more R packages that are not on yet conda (TODO: add them).
# The "options(unzip = )" hack to address the install_github issue under conda https://github.com/r-lib/devtools/issues/1722
R -e "library(devtools); options(unzip = '$(which unzip)'); devtools::install_github('mjkallen/rlogging')"
R -e "library(devtools); options(unzip = '$(which unzip)'); devtools::install_github('kent37/summarywidget')"
# This one is local
R -e "library(devtools); devtools::install('${SRC_DIR}/R/pcgrr')"

# Install VEP separately (doesn't work when within the envirnoment file, for some reason):
conda install -c bioconda -y ensembl-vep
vep_install -a ap -g miRNA -l -n

# Access to src scripts
chmod +x ${SRC_DIR}/pcgr/*.py
chmod +x ${SRC_DIR}/*.R
export PATH=$PATH:$(pwd)/${SRC_DIR}/pcgr
export PYTHONPATH=$(pwd)/${SRC_DIR}/lib:${PYTHONPATH}
