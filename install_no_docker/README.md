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
conda build -c conda-forge -c bioconda -c pcgr -c defaults install_no_docker/conda_pkg/pcgr
conda install --use-local $CHANNELS pcgr
```

Finally you can download reference data bundle for your genome build and you are all set:

```
conda activate pcgr
wget https://raw.githubusercontent.com/circulosmeos/gdown.pl/master/gdown.pl
perl gdown.pl "https://drive.google.com/uc?id=1TdYagetk-l__aYBsaZJHJvYFStDnIEcq" grch37.tar.gz
perl gdown.pl "https://drive.google.com/uc?id=1wpVqlgY5jBKkQaTAxzgf0rgxKxOEJgj-" grch38.tar.gz
tar -xzf grch37.tar.gz  # will extract into ./data/grch37/
tar -xzf grch38.tar.gz  # will extract into ./data/grch38/
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
