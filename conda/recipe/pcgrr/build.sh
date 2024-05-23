#!/bin/bash

export DISABLE_AUTOBREW=1
${R} -e "install.packages('devtools', repos = 'https://cloud.r-project.org/', lib = '${PREFIX}/lib/R/library')"
${R} -e "install.packages('openxlsx2', repos = 'https://cloud.r-project.org/', lib = '${PREFIX}/lib/R/library')"
${R} -e "devtools::install_github(repo = 'caravagnalab/CNAqc', ref = '274cde9', lib = '${PREFIX}/lib/R/library')"
${R} CMD INSTALL --build . ${R_ARGS}
