## No-docker install

This is an alternative way to install PCGR that does not require Docker on your machine. Run the following commands:

```
# Install dependencies
bash -x install_no_docker/install.sh

# Install reference data for your genome build
pip install gdown
gdown https://drive.google.com/uc?id=1MREECbygW47ttJySgfibBpX7romBrb_Q -O - | tar xvzf - # grch37
gdown https://drive.google.com/uc?id=1Xsw0WcKPnWgJDolQfrZATU5suEFJ5BKG -O - | tar xvzf - # grch38
```

The first script will install all dependencies from package repositories like Conda and CRAN. If you don't have conda 
in your path, it'll download miniconda and create a fresh environment. If you are already in an active environment, 
it will attempt to install everything into it.

This does not guarantee a successful install: due to ongoing updates of the packages in public repositories, you might 
end up with issues due to version conflics, or unsupported versions for some packages, or missing packages for your 
system. So try to stick to the dockerized version if possible.

After installing all dependencies, you will need to run download data bundles, either manually following the links in 
the [README](https://github.com/sigven/pcgr), or using googledrive CLI as shown above.

### Running

The installation script will create the `load_pcgr.sh` script, that just loads the pcgr conda environment. 
Use it before running PCGR:

```
source install_no_docker/load_pcgr.sh
```

Then run with `--no-docker` flag:

```
./pcgr.py --input_vcf examples/tumor_sample.BRCA.vcf.gz . test_out grch37 examples/pcgr_conf.BRCA.toml \
    tumor_sample.BRCA --no-docker
```
