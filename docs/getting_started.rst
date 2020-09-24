Getting started
---------------

STEP 0: Python
~~~~~~~~~~~~~~

An installation of Python (version *3.6*) is required to run PCGR. Check
that Python is installed by typing ``python --version`` in your terminal
window. In addition, a `Python library <https://github.com/uiri/toml>`__
for parsing configuration files encoded with
`TOML <https://github.com/toml-lang/toml>`__ is needed. To install,
simply run the following command:

::

   pip install toml

STEP 1: Installation of Docker
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. `Install the Docker
   engine <https://docs.docker.com/engine/installation/>`__ on your
   preferred platform

   -  installing `Docker on
      Linux <https://docs.docker.com/engine/installation/linux/>`__
   -  installing `Docker on Mac
      OS <https://docs.docker.com/engine/installation/mac/>`__
   -  NOTE: We have not yet been able to perform enough testing on the
      Windows platform, and we have received feedback that particular
      versions of Docker/Windows do not work with PCGR (an example being
      `mounting of data
      volumes <https://github.com/docker/toolbox/issues/607>`__)

2. Test that Docker is running, e.g. by typing ``docker ps`` or
   ``docker images`` in the terminal window
3. Adjust the computing resources dedicated to the Docker, i.e.:

   -  Memory: minimum 5GB
   -  CPUs: minimum 4
   -  `How to - Mac OS
      X <https://docs.docker.com/docker-for-mac/#advanced>`__

STEP 2: Download PCGR and data bundle
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Development version
^^^^^^^^^^^^^^^^^^^

a. Clone the PCGR GitHub repository (includes run script and default
   configuration file): ``git clone https://github.com/sigven/pcgr.git``

b. Download and unpack the latest data bundles in the PCGR directory

   -  `grch37 data bundle -
      20200920 <http://insilico.hpc.uio.no/pcgr/pcgr.databundle.grch37.20200920.tgz>`__
      (approx 17Gb)
   -  `grch38 data bundle -
      20200920 <http://insilico.hpc.uio.no/pcgr/pcgr.databundle.grch38.20200920.tgz>`__
      (approx 18Gb)
   -  *Unpacking*:
      ``gzip -dc pcgr.databundle.grch37.YYYYMMDD.tgz | tar xvf -``

c. Pull the `PCGR Docker image
   (dev) <https://hub.docker.com/r/sigven/pcgr/>`__ from DockerHub
   (approx 6.8Gb):

-  ``docker pull sigven/pcgr:dev`` (PCGR annotation engine)

Latest release
^^^^^^^^^^^^^^

a. Download and unpack the `latest software release
   (0.9.0rc) <https://github.com/sigven/pcgr/releases/tag/v0.9.0rc>`__

b. Download and unpack the assembly-specific data bundle in the PCGR
   directory

-  `grch37 data bundle -
   20200920 <http://insilico.hpc.uio.no/pcgr/pcgr.databundle.grch37.20200920.tgz>`__
   (approx 17Gb)
-  `grch38 data bundle -
   20200920 <http://insilico.hpc.uio.no/pcgr/pcgr.databundle.grch38.20200920.tgz>`__
   (approx 18Gb)

   -  *Unpacking*:
      ``gzip -dc pcgr.databundle.grch37.YYYYMMDD.tgz | tar xvf -``

   A *data/* folder within the *pcgr-X.X* software folder should now
   have been produced

c. Pull the `PCGR Docker image
   (0.9.0rc) <https://hub.docker.com/r/sigven/pcgr/>`__ from DockerHub
   (approx 6.8Gb):

   -  ``docker pull sigven/pcgr:0.9.0rc`` (PCGR annotation engine)

STEP 3: Input preprocessing
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The PCGR workflow accepts two types of input files:

-  An unannotated, single-sample VCF file (>= v4.2) with called somatic
   variants (SNVs/InDels)
-  A copy number segment file

PCGR can be run with either or both of the two input files present.

-  We **strongly** recommend that the input VCF is compressed and
   indexed using `bgzip <http://www.htslib.org/doc/tabix.html>`__ and
   `tabix <http://www.htslib.org/doc/tabix.html>`__
-  If the input VCF contains multi-allelic sites, these will be subject
   to `decomposition <http://genome.sph.umich.edu/wiki/Vt#Decompose>`__
-  Variants used for reporting should be designated as ‘PASS’ in the VCF
   FILTER column

The tab-separated values file with copy number aberrations **MUST**
contain the following four columns:

-  Chromosome
-  Start
-  End
-  Segment_Mean

Here, *Chromosome*, *Start*, and *End* denote the chromosomal segment,
and **Segment_Mean** denotes the log(2) ratio for a particular segment,
which is a common output of somatic copy number alteration callers. Note
that coordinates must be **one-based** (i.e. chromosomes start at 1, not
0). Below shows the initial part of a copy number segment file that is
formatted correctly according to PCGR’s requirements:

::

   Chromosome  Start   End Segment_Mean
   1 3218329 3550598 0.0024
   1 3552451 4593614 0.1995
   1 4593663 6433129 -1.0277

STEP 4: Configure your PCGR workflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The PCGR software bundle comes with a default configuration file in the
*conf/* folder, to be used as a starting point for runnning the PCGR
workflow. The configuration file, formatted using
`TOML <https://github.com/toml-lang/toml>`__, enables the user to
configure a number of options related to the following:

-  **IMPORTANT**: Designation of VCF INFO tags that denote variant
   depth/allelic fraction

   -  If this is not available/properly set, the report contents will be
      less informative **AND** it will not be possible to preset
      thresholds for depth/allelic fraction

-  Tumor-only analysis options

   -  tick on/off various filtering schemes for exclusion of germline
      variants (for input VCFs coming from tumor-only sequencing assays)

-  VEP/*vcfanno* options

More details about the exact nature of the `configuration
options <http://pcgr.readthedocs.io/en/latest/input.html#pcgr-configuration-file>`__.

STEP 5: Run example
~~~~~~~~~~~~~~~~~~~

::

   usage: pcgr.py -h [options] --input_vcf INPUT_VCF --pcgr_dir PCGR_DIR --output_dir OUTPUT_DIR --genome_assembly  GENOME_ASSEMBLY --conf CONFIG_FILE --sample_id SAMPLE_ID

   Personal Cancer Genome Reporter (PCGR) workflow for clinical interpretation of somatic nucleotide variants and copy number aberration segments

   Required arguments:
   --input_vcf INPUT_VCF
                   VCF input file with somatic variants in tumor sample, SNVs/InDels
   --pcgr_dir PCGR_DIR   PCGR base directory with accompanying data directory, e.g. ~/pcgr-0.9.0
   --output_dir OUTPUT_DIR
                   Output directory
   --genome_assembly {grch37,grch38}
                   Human genome assembly build: grch37 or grch38
   --conf CONFIGURATION_FILE
                   PCGR configuration file in TOML format
   --sample_id SAMPLE_ID
                   Tumor sample/cancer genome identifier - prefix for output files

   Optional arguments:
   --input_cna INPUT_CNA
                   Somatic copy number alteration segments (tab-separated values)
   --logr_gain LOGR_GAIN
                   Log ratio-threshold for regions containing copy number gains/amplifications (default: 0.8)
   --logr_homdel LOGR_HOMDEL
                   Log ratio-threshold for regions containing homozygous deletions (default: -0.8)
   --cna_overlap_pct CNA_OVERLAP_PCT
                   Mean percent overlap between copy number segment and gene transcripts for reporting of gains/losses in tumor suppressor genes/oncogenes, (default: 50)
   --pon_vcf PON_VCF     VCF file with germline calls from Panel of Normals (PON) - i.e. blacklisted variants, (default: None)
   --tumor_site TSITE    Optional integer code to specify primary tumor type/site of query sample,
                   choose any of the following identifiers:
                   1 = Adrenal Gland
                   2 = Ampulla of Vater
                   3 = Biliary Tract
                   4 = Bladder/Urinary Tract
                   5 = Bone
                   6 = Breast
                   7 = Cervix
                   8 = CNS/Brain
                   9 = Colon/Rectum
                   10 = Esophagus/Stomach
                   11 = Eye
                   12 = Head and Neck
                   13 = Kidney
                   14 = Liver
                   15 = Lung
                   16 = Lymphoid
                   17 = Myeloid
                   18 = Ovary/Fallopian Tube
                   19 = Pancreas
                   20 = Peripheral Nervous System
                   21 = Peritoneum
                   22 = Pleura
                   23 = Prostate
                   24 = Skin
                   25 = Soft Tissue
                   26 = Testis
                   27 = Thymus
                   28 = Thyroid
                   29 = Uterus
                   30 = Vulva/Vagina
                   (default: 0 - any tumor type)
   --tumor_purity TUMOR_PURITY
                   Estimated tumor purity (between 0 and 1, (default: None)
   --tumor_ploidy TUMOR_PLOIDY
                   Estimated tumor ploidy (default: None)
   --tumor_dp_min TUMOR_DP_MIN
                   If VCF INFO tag for sequencing depth (tumor) is provided and found, set minimum required depth for inclusion in report (default: 0)
   --tumor_af_min TUMOR_AF_MIN
                   If VCF INFO tag for variant allelic fraction (tumor) is provided and found, set minimum required AF for inclusion in report (default: 0)
   --control_dp_min CONTROL_DP_MIN
                   If VCF INFO tag for sequencing depth (control) is provided and found, set minimum required depth for inclusion in report (default: 0)
   --control_af_max CONTROL_AF_MAX
                   If VCF INFO tag for variant allelic fraction (control) is provided and found, set maximum tolerated AF for inclusion in report (default: 1)
   --target_size_mb TARGET_SIZE_MB
                   For mutational burden analysis - approximate protein-coding target size of sequencing assay (default: 34 Mb (WES/WGS))
   --tumor_only          Input VCF comes from tumor-only sequencing, calls will be filtered for variants of germline origin (set configurations for filtering in .toml file), (default: False)
   --cell_line           Input VCF comes from tumor cell line sequencing (requires --tumor_only), calls will be filtered for variants of germline origin (set configurations for filtering in .toml file), (default: False)
   --assay {WES,WGS,TARGETED}
                   Type of DNA sequencing assay performed for input data (VCF) default: WES
   --include_trials      (Beta) Include relevant ongoing or future clinical trials, focusing on studies with molecularly targeted interventions
   --estimate_tmb        Estimate tumor mutational burden from the total number of somatic mutations and target region size, default: False
   --estimate_msi_status
                   Predict microsatellite instability status from patterns of somatic mutations/indels, default: False
   --estimate_signatures
                   Estimate relative contributions of reference mutational signatures in query sample and detect potential kataegis events), default: False
   --min_mutations_signatures MIN_MUTATIONS_SIGNATURES
                   Minimum number of SNVs required for reconstruction of mutational signatures (SBS) by MutationalPatterns (default: 200, minimum n = 100)
   --all_reference_signatures
                   Use all reference mutational signatures (SBS, n = 67) in signature reconstruction rather than only those already attributed to the tumor type (default: False)
   --force_overwrite     By default, the script will fail with an error if any output file already exists. You can force the overwrite of existing result files by using this flag
   --version             show program's version number and exit
   --basic               Run functional variant annotation on VCF through VEP/vcfanno, omit other analyses (i.e. Tier assignment/MSI/TMB/Signatures etc. and report generation (STEP 4), default: False
   --no_vcf_validate     Skip validation of input VCF with Ensembl's vcf-validator
   --docker-uid DOCKER_USER_ID
                   Docker user ID. default is the host system user ID. If you are experiencing permission errors, try setting this up to root (`--docker-uid root`)
   --no-docker           Run the PCGR workflow in a non-Docker mode (see install_no_docker/ folder for instructions)
   --debug               Print full Docker commands to log, default: False

The *examples* folder contain input VCF files from two tumor samples
sequenced within TCGA (**GRCh37** only). It also contains a PCGR
configuration file customized for these VCFs. A report for a colorectal
tumor case can be generated by running the following command in your
terminal window:

::

   python ~/pcgr-0.9.0rc/pcgr.py
   --pcgr_dir ~/pcgr-0.9.0rc
   --output_dir ~/pcgr-0.9.0rc
   --sample_id tumor_sample.COAD
   --genome_assembly grch37
   --conf ~/pcgr-0.9.0rc/examples/example_COAD.toml
   --input_vcf ~/pcgr-0.9.0rc/examples/tumor_sample.COAD.vcf.gz
   --tumor_site 9
   --input_cna ~/pcgr-0.9.0rc/examples/tumor_sample.COAD.cna.tsv
   --tumor_purity 0.9
   --tumor_ploidy 2.0
   --include_trials
   --assay WES
   --estimate_signatures
   --estimate_msi_status
   --estimate_tmb
   --no_vcf_validate

This command will run the Docker-based PCGR workflow and produce the
following output files in the *examples* folder:

1. **tumor_sample.COAD.pcgr_acmg.grch37.html** - An interactive HTML
   report for clinical interpretation
2. **tumor_sample.COAD.pcgr_acmg.grch37.pass.vcf.gz** - Bgzipped VCF
   file with rich set of annotations for precision oncology
3. **tumor_sample.COAD.pcgr_acmg.grch37.pass.tsv.gz** - Compressed
   vcf2tsv-converted file with rich set of annotations for precision
   oncology
4. **tumor_sample.COAD.pcgr_acmg.grch37.snvs_indels.tiers.tsv** -
   Tab-separated values file with variants organized according to tiers
   of functional relevance
5. **tumor_sample.COAD.pcgr_acmg.grch37.mutational_signatures.tsv** -
   Tab-separated values file with information on contribution of
   mutational signatures
6. **tumor_sample.COAD.pcgr_acmg.grch37.json.gz** - Compressed JSON dump
   of HTML report content
7. **tumor_sample.COAD.pcgr_acmg.grch37.cna_segments.tsv.gz** -
   Compressed tab-separated values file with annotations of gene
   transcripts that overlap with somatic copy number aberrations
