## No-docker install

This is an alternative way to install PCGR that does not require Docker on your machine. Instead of pulling the image, 
run the following script:

```
cd install_no_docker
bash install.sh
```

This script installs all dependencies from package repositories like Conda and CRAN. If you don't have conda in your 
path, it'll download miniconda and create a fresh environment. If you are already in an active environment, it will 
attempt to install everything into it.

This does not guarantee a successful install: due to ongoing updates of the packages in public repositories, you might 
end up with issues due to version conflics, or unsupported versions for some packages, or missing packages for your 
system. So try to stick to the dockerized version if possible.

After installing all dependencies, you will still need to download and uncompress data bundles, according to the 
[README](https://github.com/sigven/pcgr). You can do it from the command line with the gdrive CLI:

```
cd ..
git clone https://github.com/circulosmeos/gdown.pl
gdown.pl/gdown.pl https://drive.google.com/file/d/1cGBAmAh5t4miIeRrrd0zHsPCFToOr0Lf/view pcgr.databundle.grch37.tgz
gzip -dc pcgr.databundle.grch37.tgz | tar xvf -
```

### Loading

The installation script will create the `load_pcgr.sh` script, that just loads the pcgr conda environment. 
Use it before running PCGR:

```
source install_no_docker/load_pcgr.sh
./pcgr.py --input_vcf examples/tumor_sample.BRCA.vcf.gz --no-docker . test_out grch37 examples/pcgr_conf.BRCA.toml tumor_sample.BRCA
```