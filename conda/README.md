Conda Packaging
===============

R package
---------

- Test conda building locally with Docker
  (result will be output to `$(PWD)/conda/out/noarch/r-pcgrr-0.9.2-r40_0.tar.bz2`):

```
cd path/to/pcgr

docker run -it --rm \
  -v $(PWD):/home/pcgr_root \
  -v $(PWD)/conda/out:/home/condabuild_out \
  quay.io/condaforge/linux-anvil-comp7:latest \
  conda mambabuild /home/pcgr_root/conda/pcgrr --output-folder /home/condabuild_out -c conda-forge -c bioconda
```

- Should take around 7-8 min (MacOS Big Sur 2017)
