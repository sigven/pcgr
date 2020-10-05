#!/usr/bin/env bash
set -x

R -e "remove.packages('crosstalk')"
export TAR=$(which tar)
R -e "library(devtools); options(unzip = 'internal'); devtools::install_github('rstudio/crosstalk', ref = '68b0b617ee82e5d6b738f26106933254ce5ede53', dependencies=FALSE)"
