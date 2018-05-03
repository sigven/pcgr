# Install conda
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
bash miniconda.sh -b -p ./miniconda
export PATH=$(pwd)/miniconda/bin:$PATH

# Make pcgr environment
conda env create --file conda_environment.yml
source activate pcgr

# Few more R packages that are not on yet conda (TODO: add them):
R -e "library(devtools); devtools::install_github('mjkallen/rlogging'); devtools::install_github('kent37/summarywidget')"

# Install VEP separately (doesn't work when within the envirnoment file, for some reason):
conda install -c bioconda -y ensembl-vep
vep_install -a ap -g miRNA -l -n

# Access to src scripts
chmod +x ../src/pcgr/*.py
export PATH=$PATH:$(pwd)/../src/pcgr
export PYTHONPATH=$(pwd)/../src/pcgr/lib:${PYTHONPATH}
