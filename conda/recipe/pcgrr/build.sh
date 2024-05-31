#!/bin/bash

export DISABLE_AUTOBREW=1
${R} -e "install.packages('remotes', repos = 'https://cloud.r-project.org/', lib = '${PREFIX}/lib/R/library')"
${R} -e "remotes::install_github(repo = 'caravagnalab/CNAqc', ref = '274cde9', lib = '${PREFIX}/lib/R/library')"
${R} CMD INSTALL --build . ${R_ARGS}
