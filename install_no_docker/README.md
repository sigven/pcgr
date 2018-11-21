## Installing PCGR using conda

This is an alternative installation approach that does not require Docker on your machine. At the moment it works only in linux systems.

First you need to download conda package manager. Get it with:

```
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
bash miniconda.sh -b -p ./miniconda
```

Run the following to add conda into your PATH. You can even put that into your `~/.bashrc` or `~/.zshrc` to avoid re-running this in the future:

```
. ./miniconda/etc/profile.d/conda.sh
```

Now create a new environment and install "pcgr" conda package into it:

```
conda create -n pcgr -c conda-forge -c bioconda -c pcgr pcgr
conda activate pcgr
```

Or alternatively build the package from source, which is useful if you need to use the development code from the repository:

```
conda install conda-build
conda build -c conda-forge -c bioconda -c pcgr -c defaults install_no_docker/conda_package/pcgr
conda install --use-local $CHANNELS pcgr
```

Finally you can download reference data bundle for your genome build and you are all set:

```
gdown https://drive.google.com/uc?id=1OL5C994HDaeadASz7KzMhPoXfdSiyhNy -O - | tar xvzf - # grch37
gdown https://drive.google.com/uc?id=1CZNc87E0K5AK2RDSNU57FqLp0H1skpUh -O - | tar xvzf - # grch38
```

There is a chance you'll encounter errors during the installation. Due to ongoing updates of the packages in public repositories, some packages might end up conflicting with each other or missing for your system. So try to stick to the dockerized version of PCGR whenever possible.

## Running condarized PCGR

Activate your environment with:

```
conda activate pcgr
```

Run PCGR with `--no-docker` flag. The first agument ("path to pcgr") now doesn't have to contain anything but a `data` directory that you downloaded.

```
pcgr.py --input_vcf examples/tumor_sample.BRCA.vcf.gz . test_out grch37 examples/pcgr_conf.BRCA.toml \
    tumor_sample.BRCA --no-docker
```
