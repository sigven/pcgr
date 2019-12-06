#!/usr/bin/env bash
set -x

R -e "remove.packages('crosstalk')"
export TAR=$(which tar)
R -e "library(devtools); options(unzip = 'internal'); devtools::install_github('rstudio/crosstalk', dependencies=FALSE)"
