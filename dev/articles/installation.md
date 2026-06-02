# Installation

## Data

PCGR requires the following data to run successfully:

1.  **Reference bundle** - containing data from multiple knowledge
    resources, including information on molecular biomarkers, targeted
    cancer therapies, variant frequencies etc. Key datasets include
    [CIViC](https://civicdb.org),
    [CGI](https://www.cancergenomeinterpreter.org/biomarkers), [Open
    Targets Platform](https://platform.opentargets.org),
    [TCGA](https://portal.gdc.cancer.gov/), and
    [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/).
2.  **Ensembl VEP data cache** - needed for variant annotation with VEP
    (Variant Effect Predictor)
3.  **User-supplied sample-specific inputs** - e.g. somatic variant
    calls in VCF format

PCGR supports both the GRCh37 and GRCh38 human genome assemblies. All
the data above need to match the chosen assembly.

### 1. Reference Bundle

Reference bundles are generated semi-automatically (by the PCGR author)
and are versioned based on their release date. Keep in mind that the
bundles support only certain Ensembl VEP versions, and importantly also
specific software versions of PCGR. Upgrading the PCGR software without
upgrading the bundle specified here is thus *not* a recommended
installation strategy. The latest (**v20260508**) genome-specific
bundles can be downloaded directly from below (size: ~5G):

| Assembly | Download Link |
|:---|:---|
| GRCh38 | <https://insilico.hpc.uio.no/pcgr/pcgr_ref_data.20260508.grch38.tgz> |
| GRCh37 | <https://insilico.hpc.uio.no/pcgr/pcgr_ref_data.20260508.grch37.tgz> |

**Tip 1**: The `data/grch3x/.PCGR_BUNDLE_VERSION` file within the
downloaded bundle indicates the bundle version for reporting purposes.

**Tip 2**: The `data/grch3x/data_overview.grch3x.html` file provides a
report with an overview and statistics of the key resources included in
the reference bundle.

#### Bash Example

    BUNDLE_VERSION="20260508"

``` bash
GENOME="grch38" # or "grch37"
BUNDLE="pcgr_ref_data.${BUNDLE_VERSION}.${GENOME}.tgz"
wget https://insilico.hpc.uio.no/pcgr/${BUNDLE}
gzip -dc ${BUNDLE} | tar xvf -

mkdir ${BUNDLE_VERSION}
mv data/ ${BUNDLE_VERSION}
```

### 2. VEP Cache

[VEP](https://www.ensembl.org/info/docs/tools/vep/index.html) requires a
data cache which is available from the Ensembl [FTP
site](https://ftp.ensembl.org/pub/release-115/variation/indexed_vep_cache/)
(search there for files starting with `homo_sapiens_vep_`). The latest
Ensembl VEP version we support is **v115**.

#### Bash Example

    VEP_VERSION="115"

``` bash
GENOME="GRCh38" # or "GRCh37"
CACHE="homo_sapiens_vep_${VEP_VERSION}_${GENOME}.tar.gz"

wget https://ftp.ensembl.org/pub/release-${VEP_VERSION}/variation/indexed_vep_cache/${CACHE}
gzip -dc ${CACHE} | tar xvf -
```

**Important**: PCGR needs to be pointed to the *parent* directory
containing the downloaded `homo_sapiens/115_GRCh3x/` cache.
Historically, this parent directory has been named `.vep`, but the name
of this directory can be arbitrarily set.

### 3. Sample Inputs

See the [Inputs](https://sigven.github.io/pcgr/dev/articles/input.md)
article.

------------------------------------------------------------------------

## Software

The PCGR workflow can be installed with any of the following:

- A.
  [Conda](https://conda.io/projects/conda/en/latest/user-guide/getting-started.html)
  (Linux only),
- B. [Docker](https://docs.docker.com/), or
- C.
  [Singularity/Apptainer](https://apptainer.org/docs/user/latest/index.html).

### A. Conda

There is Conda support for Linux machines. The following process can
take anywhere from 10 up to 40 minutes when installing from scratch,
mostly depending on the user’s and server’s internet connection. Most of
the time is spent on downloading the `{BSgenome.Hsapiens.UCSC.hg19}` and
`{BSgenome.Hsapiens.UCSC.hg38}` R packages (which happens at the very
end of the conda environment creation).

    PCGR_VERSION="2.2.5.9020"

``` bash
# set up variables
PCGR_REPO="https://raw.githubusercontent.com/sigven/pcgr/v${PCGR_VERSION}/conda/env/lock/"
PLATFORM="linux"
# create conda envs in local directory
mkdir pcgr_conda
conda create --prefix ./pcgr_conda/pcgr --file ${PCGR_REPO}/pcgr-${PLATFORM}-64.lock
conda create --prefix ./pcgr_conda/pcgrr --file ${PCGR_REPO}/pcgrr-${PLATFORM}-64.lock
# you need to specify the directory of the conda env when using --prefix
conda activate ./pcgr_conda/pcgr
# test that it works
pcgr --version
pcgr --help
```

### B. Docker

The PCGR Docker image is available from the GitHub Container Registry at
<https://github.com/sigven/pcgr/pkgs/container/pcgr>. Pull the latest
**v2.2.5.9020** image with:

    docker pull ghcr.io/sigven/pcgr:2.2.5.9020

#### Example Run

In the example below, we assume that you have the following directory
structure:

- `/Users/you/dir1` - the directory you are working in
- `/Users/you/dir1/VEP_cache` - the directory where you have downloaded
  the [VEP cache](#vep-cache) (**Note**: this directory should contain
  the `homo_sapiens` directory as its child directory
- `/Users/you/dir1/20260508` - the directory where you have downloaded
  the PCGR reference bundle (**Note**: this directory should contain a
  directory named `data` as its child directory)
- `/Users/you/dir1/pcgr_inputs` - the directory where you have your
  input files (e.g. VCF files)
- `/Users/you/dir1/pcgr_outputs` - the directory where you want to store
  the PCGR outputs

&nbsp;

    PCGR_VERSION="2.2.5.9020"
    BUNDLE_VERSION="20260508"

``` bash
docker container run -it --rm \
    -v /Users/you/dir1/VEP_cache:/mnt/.vep \
    -v /Users/you/dir1/${BUNDLE_VERSION}:/mnt/bundle \
    -v /Users/you/dir1/pcgr_inputs:/mnt/pcgr_inputs \
    -v /Users/you/dir1/pcgr_outputs:/mnt/pcgr_outputs \
    ghcr.io/sigven/pcgr:${PCGR_VERSION} \
    pcgr \
      --input_vcf "/mnt/pcgr_inputs/T001-BRCA.grch38.vcf.gz" \
      --vep_dir "/mnt/.vep" \
      --refdata_dir "/mnt/bundle" \
      --output_dir "/mnt/pcgr_outputs" \
      --genome_assembly "grch38" \
      --sample_id "SAMPLE_B" \
      --tumor_dp_tag "TDP" \
      --tumor_af_tag "TAF" \
      --assay "WGS" \
      --vcf2maf
```

**NOTE**: If you need to run the Docker-based version of PCGR as a
non-root user, you may need to explicitly add options for quarto to work
properly, i.e. `--env "XDG_CACHE_HOME=/tmp/quarto_cache_home"` (same as
for Singularity/Apptainer below, see also [issue
\#246](https://github.com/sigven/pcgr/issues/246)).  
  

### C. Singularity/Apptainer

The PCGR Singularity/Apptainer image is available on [GitHub Container
Registry](https://ghcr.io/sigven/pcgr). Pull the latest **v2.2.5.9020**
image with:

    apptainer pull oras://ghcr.io/sigven/pcgr:2.2.5.9020.singularity

This will download a Singularity Image File (SIF) called
**pcgr_2.2.5.9020.singularity.sif** that can be run with Singularity or
Apptainer.

#### Example Run

    PCGR_VERSION="2.2.5.9020"
    BUNDLE_VERSION="20260508"

``` bash
apptainer exec \
  --writable-tmpfs \
  --no-home \
  --env "XDG_CACHE_HOME=/tmp/quarto_cache_home" \
  -B /Users/you/dir1/VEP_cache:/mnt/.vep \
  -B /Users/you/dir1/${BUNDLE_VERSION}:/mnt/bundle \
  -B /Users/you/dir1/pcgr_inputs:/mnt/pcgr_inputs \
  -B /Users/you/dir1/pcgr_outputs:/mnt/pcgr_outputs \
  pcgr_${PCGR_VERSION}.singularity.sif \
  pcgr \
    --input_vcf "/mnt/pcgr_inputs/T001-BRCA.grch38.vcf.gz" \
    --vep_dir "/mnt/.vep" \
    --refdata_dir "/mnt/bundle" \
    --output_dir "/mnt/pcgr_outputs" \
    --genome_assembly "grch38" \
    --sample_id "SAMPLE_B" \
    --assay "WGS" \
    --tumor_dp_tag "TDP" \
    --tumor_af_tag "TAF" \
    --vcf2maf
```

**Note**: For any Apptainer/Singularity issues not directly related to
PCGR, we would recommend reaching out to the Apptainer community
(e.g. <https://github.com/apptainer/apptainer>) since we have limited
experience with Apptainer/Singularity.
