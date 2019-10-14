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
   -  **NOTE**: For query VCF’s that are very large (size to be
      determined) one may need to adjust these settings

STEP 2: Download PCGR and data bundle
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Development version
^^^^^^^^^^^^^^^^^^^

a. Clone the PCGR GitHub repository:
   ``git clone https://github.com/sigven/pcgr.git``
b. Download and unpack the latest data bundles in the PCGR directory

   -  `grch37 data bundle -
      20191013 <https://drive.google.com/open?id=17DUuy4Z_zid6YIvT8SSYqVzrpzZcoS1i>`__
      (approx 16Gb)
   -  `grch38 data bundle -
      20191013 <https://drive.google.com/open?id=1tuKOmExGRCwte4rrP7yKogCxnq6I4fiH>`__
      (approx 17Gb)
   -  *Unpacking*:
      ``gzip -dc pcgr.databundle.grch37.YYYYMMDD.tgz | tar xvf -``

c. Pull the `PCGR Docker image
   (dev) <https://hub.docker.com/r/sigven/pcgr/>`__ from DockerHub
   (approx 5.2Gb):

   -  ``docker pull sigven/pcgr:dev`` (PCGR annotation engine)

Latest release
^^^^^^^^^^^^^^

a. Download and unpack the `latest software release
   (0.8.3) <https://github.com/sigven/pcgr/releases/tag/v0.8.3>`__
b. Download and unpack the assembly-specific data bundle in the PCGR
   directory

-  `grch37 data bundle -
   20191013 <https://drive.google.com/open?id=17DUuy4Z_zid6YIvT8SSYqVzrpzZcoS1i>`__
   (approx 16Gb)
-  `grch38 data bundle -
   20191013 <https://drive.google.com/open?id=1tuKOmExGRCwte4rrP7yKogCxnq6I4fiH>`__
   (approx 17Gb)

   -  *Unpacking*:
      ``gzip -dc pcgr.databundle.grch37.YYYYMMDD.tgz | tar xvf -``

   A *data/* folder within the *pcgr-X.X* software folder should now
   have been produced

c. Pull the `PCGR Docker image
   (0.8.3) <https://hub.docker.com/r/sigven/pcgr/>`__ from DockerHub
   (approx 5.2Gb):

   -  ``docker pull sigven/pcgr:0.8.3`` (PCGR annotation engine)

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

    Chromosome Start   End Segment_Mean
    1 3218329 3550598 0.0024
    1 3552451 4593614 0.1995
    1 4593663 6433129 -1.0277

STEP 4: Configure PCGR
~~~~~~~~~~~~~~~~~~~~~~

The PCGR software bundle comes with a default configuration file in the
*conf/* folder, to be used as a starting point for runnning the PCGR
workflow. The configuration file, formatted using
`TOML <https://github.com/toml-lang/toml>`__, enables the user to
configure a number of options related to the following:

-  Sequencing depth/allelic support (definition of tags + thresholds)
-  MSI prediction
-  Mutational signatures analysis
-  Mutational burden analysis (e.g. target size)
-  VCF to MAF conversion
-  Tumor-only analysis options (i.e. exclusion of germline
   variants/enrichment for somatic calls)
-  VEP/\ *vcfanno* options
-  Log-ratio thresholds for gains/losses in CNA analysis

More details about the exact `usage of the configuration
options <http://pcgr.readthedocs.io/en/latest/input.html#pcgr-configuration-file>`__.

STEP 5: Run example
~~~~~~~~~~~~~~~~~~~

A tumor sample report is generated by calling the Python script
**pcgr.py**, which takes the following arguments and options:

::

   usage: pcgr.py -h [options] <PCGR_DIR> <OUTPUT_DIR> <GENOME_ASSEMBLY> <CONFIG_FILE> <SAMPLE_ID>

   Personal Cancer Genome Reporter (PCGR) workflow for clinical interpretation of somatic nucleotide variants and copy number aberration segments

   positional arguments:
     pcgr_dir              PCGR base directory with accompanying data directory, e.g. ~/pcgr-0.8.3
     output_dir            Output directory
     {grch37,grch38}       Genome assembly build: grch37 or grch38
     configuration_file    PCGR configuration file (TOML format)
     sample_id             Tumor sample/cancer genome identifier - prefix for output files

   optional arguments:
     -h, --help            show this help message and exit
     --input_vcf INPUT_VCF
                           VCF input file with somatic query variants (SNVs/InDels).
     --input_cna INPUT_CNA
                           Somatic copy number alteration segments (tab-separated values)
     --input_cna_plot INPUT_CNA_PLOT
                           Somatic copy number alteration plot
     --pon_vcf PON_VCF     VCF file with germline calls from Panel of Normals (PON) - i.e. blacklisted variants, (default: None)
     --tumor_type TTYPE    Optional integer code to specify tumor type of query,
                            choose any of the following identifiers:
                           1 = Adrenal_Gland_Cancer_NOS
                           2 = Ampullary_Carcinoma_NOS
                           3 = Biliary_Tract_Cancer_NOS
                           4 = Bladder_Urinary_Tract_Cancer_NOS
                           5 = Bone_Cancer_NOS
                           6 = Breast_Cancer_NOS
                           7 = CNS_Brain_Cancer_NOS
                           8 = Cancer_Unknown_Primary_NOS
                           9 = Cervical_Cancer_NOS
                           10 = Colorectal_Cancer_NOS
                           11 = Esophageal_Cancer_NOS
                           12 = Head_And_Neck_Cancer_NOS
                           13 = Kidney_Cancer_NOS
                           14 = Leukemia_NOS
                           15 = Liver_Cancer_NOS
                           16 = Lung_Cancer_NOS
                           17 = Lymphoma_Hodgkin_NOS
                           18 = Lymphoma_Non_Hodgkin_NOS
                           19 = Multiple_Myeloma
                           20 = Ocular_Cancer_NOS
                           21 = Ovarian_Fallopian_Tube_Cancer_NOS
                           22 = Pancreatic_Cancer_NOS
                           23 = Penile_Cancer_NOS
                           24 = Peripheral_Nervous_System_Cancer_NOS
                           25 = Peritoneal_Cancer_NOS
                           26 = Pleural_Cancer_NOS
                           27 = Prostate_Cancer_NOS
                           28 = Skin_Cancer_NOS
                           29 = Soft_Tissue_Cancer_Sarcoma_NOS
                           30 = Stomach_Cancer_NOS
                           31 = Testicular_Cancer_NOS
                           32 = Thymic_Cancer_NOS
                           33 = Thyroid_Cancer_NOS
                           34 = Uterine_Cancer_NOS
                           35 = Vulvar_Vaginal_Cancer_NOS
                           (default: 0 - any tumor type)
     --tumor_purity TUMOR_PURITY
                           Estimated tumor purity (between 0 and 1, (default: None)
     --tumor_ploidy TUMOR_PLOIDY
                           Estimated tumor ploidy (default: None)
     --target_size_mb TARGET_SIZE_MB
                           For mutational burden analysis - approximate protein-coding target size of sequencing assay (default: 34 Mb (WES))
     --tumor_only          Input VCF comes from tumor-only sequencing, calls will be filtered for variants of germline origin (set configurations for filtering in .toml file), (default: False)
     --force_overwrite     By default, the script will fail with an error if any output file already exists. You can force the overwrite of existing result files by using this flag
     --version             show program's version number and exit
     --basic               Run functional variant annotation on VCF through VEP/vcfanno, omit other analyses (i.e. CNA, MSI, report generation etc. (STEP 4)
     --no_vcf_validate     Skip validation of input VCF with Ensembl's vcf-validator
     --docker-uid DOCKER_USER_ID
                           Docker user ID. Default is the host system user ID. If you are experiencing permission errors, try setting this up to root (`--docker-uid root`)
     --no-docker           Run the PCGR workflow in a non-Docker mode (see install_no_docker/ folder for instructions
     --debug               Print full docker commands to log

The *examples* folder contain input files from two tumor samples
sequenced within TCGA (**GRCh37** only). It also contains PCGR
configuration files customized for these samples. A report for a
colorectal tumor case can be generated by running the following command
in your terminal window:

``python pcgr.py --input_vcf ~/pcgr-0.8.3/examples/tumor_sample.COAD.vcf.gz --tumor_type 10``
``--input_cna ~/pcgr-0.8.3/examples/tumor_sample.COAD.cna.tsv --tumor_purity 0.9 --tumor_ploidy 2.0``
``~/pcgr-0.8.3 ~/pcgr-0.8.3/examples grch37 ~/pcgr-0.8.3/examples/examples_COAD.toml tumor_sample.COAD``

This command will run the Docker-based PCGR workflow and produce the
following output files in the *examples* folder:

1. **tumor_sample.COAD.pcgr_acmg.grch37.html** - An interactive HTML
   report for clinical interpretation
2. **tumor_sample.COAD.pcgr_acmg.grch37.pass.vcf.gz (.tbi)** - Bgzipped
   VCF file with rich set of annotations for precision oncology
3. **tumor_sample.COAD.pcgr_acmg.grch37.pass.tsv.gz** - Compressed
   vcf2tsv-converted file with rich set of annotations for precision
   oncology
4. **tumor_sample.COAD.pcgr_acmg.grch37.snvs_indels.tiers.tsv** -
   Tab-separated values file with variants organized according to tiers
   of functional relevance
5. **tumor_sample.COAD.pcgr_acmg.grch37.json.gz** - Compressed JSON dump
   of HTML report content
6. **tumor_sample.COAD.pcgr_acmg.grch37.cna_segments.tsv.gz** -
   Compressed tab-separated values file with annotations of gene
   transcripts that overlap with somatic copy number aberrations
