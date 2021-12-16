#!/usr/bin/env bash

# Build conda package for pcgr or pcgrr

RECIPE="pcgr"
#RECIPE="pcgrr"

docker run -it --rm \
  -v $(PWD):/home/pcgr/conda \
  quay.io/condaforge/linux-anvil-comp7:latest \
  conda mambabuild /home/pcgr/conda/recipe/${RECIPE} \
    --output-folder /home/pcgr/conda/recipe/out \
    -c conda-forge \
    -c bioconda \
    -c defaults
