### Installation of PCGR using Conda

This is an alternative installation approach that does not require Docker on your machine. At the moment it works only for Linux systems.

A prerequisite is that you have Conda installed. First download the Conda package manager. Get it with:

```
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
bash miniconda.sh -b -p ./miniconda
```

Run the following to add Conda into your PATH. You can even put that into your `~/.bashrc` or `~/.zshrc` to avoid re-running this in the future:

```
. ./miniconda/etc/profile.d/conda.sh
```

#### Alternative 1
Create a new environment (`-n pcgr`), install the _pcgr_ Conda package directly from Anaconda Cloud:

```
conda create -n pcgr -c conda-forge -c bioconda -c pcgr pcgr

```

#### Alternative 2
Build the _pcgr_ package from source, which is useful if you need to use the development code from the repository:

```
conda install conda-build
export CHANNELS="-c conda-forge -c bioconda -c pcgr -c defaults"
conda build $CHANNELS install_no_docker/conda_pkg/pcgr
conda install --use-local $CHANNELS pcgr
```

For both alternatives you also need to download the reference data bundle for your genome build:

```
wget http://insilico.hpc.uio.no/pcgr/pcgr.databundle.grch37.20210627.tgz -O grch37.tgz
wget http://insilico.hpc.uio.no/pcgr/pcgr.databundle.grch38.20210627.tgz -O grch38.tgz
tar -xzf grch37.tar.gz  # will extract into ./data/grch37/
tar -xzf grch38.tar.gz  # will extract into ./data/grch38/
```

There is a chance you'll encounter errors during the installation. Due to ongoing updates of the packages in public repositories, some packages might end up conflicting with each other or missing for your system. So try to stick to the dockerized version of PCGR whenever possible.

### Running condarized PCGR

Activate your environment with:

```
conda activate pcgr
```

Run PCGR with `--no-docker` flag. The `--pcgr_dir` argument now doesn't have to contain anything but a `data` directory that you downloaded.

```
pcgr.py
--pcgr_dir .
--output_dir test_out
--sample_id tumor_sample.BRCA
--genome_assembly grch37
--conf examples/example_BRCA.toml
--input_vcf examples/tumor_sample.BRCA.vcf.gz
--tumor_site 9
--input_cna examples/tumor_sample.BRCA.cna.tsv
--tumor_purity 0.9
--tumor_ploidy 2.0
--include_trials
--assay WES
--estimate_signatures
--estimate_msi_status
--estimate_tmb
--no_vcf_validate
--no-docker
```
If you encounter errors with VEP, you may need to unset/reset PERL5LIB (e.g. `export PERL5LIB=""`), see the [following issue](https://github.com/bioconda/bioconda-recipes/issues/4390)
