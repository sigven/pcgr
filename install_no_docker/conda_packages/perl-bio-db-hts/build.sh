#!/bin/bash

# export HTSLIB_DIR=/Users/vsaveliev/miniconda3/lib/python3.6/site-packages/pysam/include/htslib

perl ./Build.PL --extra_compiler_flags "-I$PREFIX/include"
perl ./Build
# Make sure this goes in site
perl ./Build install --installdirs site
