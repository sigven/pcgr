Conda Packaging
===============

Test conda building locally with Docker (and plain MacOS).

R package
---------

Build `noarch` R package:

```
cd path/to/pcgr

docker run -it --rm \
  -v $(PWD):/home/pcgr_root \
  -v $(PWD)/conda/out:/home/condabuild_out \
  quay.io/condaforge/linux-anvil-comp7:latest \
  conda mambabuild /home/pcgr_root/conda/pcgrr \
    --output-folder /home/condabuild_out \
    -c conda-forge \
    -c bioconda
```

- result will be output to `conda/out/noarch/r-pcgrr-0.9.2-r41_0.tar.bz2`):
- should take around 7-8 min (MacOS Big Sur 2017)

Python package
--------------
The `vcf_validator` repo has a binary for Linux and Mac (and no conda release),
so need to build separately for each:

### Linux

```
docker run -it --rm \
  -v $(PWD):/home/pcgr_root \
  -v $(PWD)/conda/out:/home/condabuild_out \
  quay.io/condaforge/linux-anvil-comp7:latest \
  conda mambabuild /home/pcgr_root/conda/pcgr \
    --output-folder /home/condabuild_out \
    -c conda-forge \
    -c bioconda \
    -c /home/condabuild_out # grabs the built R package from previous step!
```

- result will be output to `conda/out/linux-64/pcgr-0.9.2-py36r41_0.tar.bz2`):
- should take around 10 min (MacOS Big Sur 2017)

### MacOS

```
# create conda build environment
mamba create -n cbuild -c conda-forge conda-build conda-verify anaconda-client boa

# activate that env
conda activate cbuild

# now use conda mambabuild which uses boa that is much faster than conda-build
conda mambabuild conda/pcgr \
  --output-folder $(PWD)/conda/out2 \
  -c conda-forge \
  -c bioconda \
  -c defaults \
  -c $(PWD)/conda/out # grabs the already built (noarch) R package
```


- result will be output to `conda/out2/osx-64/pcgr-0.9.2-py36r41_0.tar.bz2`
- should take around 13 min (MacOS Big Sur 2017)
- expect lots of warnings and noise in the screen - you've been warned!
