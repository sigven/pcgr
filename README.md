## Personal Cancer Genome Reporter (PCGR) - variant interpretation for precision cancer medicine

### Contents

- [Overview](#overview)
- [News](#news)
- [Example reports](#example-reports)
- [PCGR Documentation](#documentation)
- [Annotation resources](#annotation-resources-included-in-pcgr---0.9.0)
- [Getting started](#getting-started)
- [FAQ](#faq)
- [Contact](#contact)

### Overview

The Personal Cancer Genome Reporter (PCGR) is a stand-alone software package for functional annotation and translation of individual cancer genomes for precision cancer medicine. Currently, it interprets both somatic SNVs/InDels and copy number aberrations. The software extends basic gene and variant annotations from the [Ensembl’s Variant Effect Predictor (VEP)](http://www.ensembl.org/info/docs/tools/vep/index.html) with oncology-relevant, up-to-date annotations retrieved flexibly through [vcfanno](https://github.com/brentp/vcfanno), and produces interactive HTML reports intended for clinical interpretation.

A few screenshots of the dashboard-type HTML output (new in 0.9.0) is shown below.


![PCGR overview](pcgr_dashboard_views.png)

### News
* _Sep 24th 2020_: **0.9.0rc release**
   * Major data bundle updates (CIViC, ClinVar, CancerMine, UniProt KB, Open Targets Platform, Pfam, DisGeNET, GENCODE)
   * VEP v101
   * _New_: report type output - [flexdashboard](https://rmarkdown.rstudio.com/flexdashboard/) HTML
   * _New_: detection of kataegis events
   * _New_: improved distinction between exact and regional biomarkers in tables
   * _New_: addition of other drug targets in copy number view
   * _New_: Inclusion of annotated molecularly targeted clinical trials (beta)
   * _Changed_: estimation of mutational signatures - using [MutationalPatterns](https://github.com/UMCUGenetics/MutationalPatterns) with the [n = 67 reference collection (COSMIC v3)](https://cancer.sanger.ac.uk/cosmic/signatures)
   * _Changed_: a number of options in the configuration file is now arguments for the main Python script
   * __NOTE__: Non-Dockerized version (Conda-based) is in the making
   * see [CHANGELOG](http://pcgr.readthedocs.io/en/latest/CHANGELOG.html)
* _Nov 18th 2019_: **0.8.4 release**
   * Data bundle updates (CIViC, ClinVar, CancerMine, UniProt)
   * Software updates: VEP 98.3
* _Oct 14th 2019_: **0.8.3 release**
   * Software updates (VEP 98.2)
   * Data bundle updates (CIViC, ClinVar, CancerMine)
   * Bug fixing, see [CHANGELOG](http://pcgr.readthedocs.io/en/latest/CHANGELOG.html#oct-14th-2019)
* _Sep 29th 2019_: **0.8.2 release**
   * Software updates (VEP 97.3, vcfanno 0.3.2)
   * Data bundle updates (CIViC, CancerMine, Open Targets Platform, UniProt KB, GENCODE, ClinVar)
   * [CHANGELOG](http://pcgr.readthedocs.io/en/latest/CHANGELOG.html#sep-29th-2019)
   * Accompanying release of the [Cancer Predisposition Sequencing Reporter](https://github.com/sigven/cpsr)


### Example reports

* [Cervical cancer sample (tumor-only)](http://insilico.hpc.uio.no/pcgr/example_reports/0.9.0rc/TCGA-FU-A3HZ-01A_TO.pcgr_acmg.grch37.flexdb.html)
* [Lung cancer sample (tumor-control)](http://insilico.hpc.uio.no/pcgr/example_reports/0.9.0rc/TCGA-95-7039-01A.pcgr_acmg.grch37.flexdb.html)
* [Breast cancer sample (tumor-control)](http://insilico.hpc.uio.no/pcgr/example_reports/0.9.0rc/TCGA-EW-A1J5-01A.pcgr_acmg.grch37.flexdb.html)
* [Brain cancer sample (tumor-control)](http://insilico.hpc.uio.no/pcgr/example_reports/0.9.0rc/TCGA-14-0866-01B.pcgr_acmg.grch37.flexdb.html)

(to view the rmarkdown-based reports, simply remove _.flexdb._ in the file names for the flexdashboard reports)



[![Build Status](https://travis-ci.org/sigven/pcgr.svg?branch=master)](https://travis-ci.org/sigven/pcgr)

### PCGR documentation

[![Documentation Status](https://readthedocs.org/projects/pcgr/badge/?version=latest)](http://pcgr.readthedocs.io/en/latest/?badge=latest)

**IMPORTANT**: If you use PCGR, please cite the publication:

Sigve Nakken, Ghislain Fournous, Daniel Vodák, Lars Birger Aaasheim, Ola Myklebost, and Eivind Hovig. __Personal Cancer Genome Reporter: variant interpretation report for precision oncology__ (2017). _Bioinformatics_. 34(10):1778–1780. doi:[10.1093/bioinformatics/btx817](https://doi.org/10.1093/bioinformatics/btx817)

### Annotation resources included in PCGR - 0.9.0

* [VEP](http://www.ensembl.org/info/docs/tools/vep/index.html) - Variant Effect Predictor v101 (GENCODE v35/v19 as the gene reference dataset)
* [CIViC](http://civic.genome.wustl.edu) - Clinical interpretations of variants in cancer (September 20th 2020)
* [ClinVar](http://www.ncbi.nlm.nih.gov/clinvar/) - Database of variants with clinical significance (August 2020)
* [DoCM](http://docm.genome.wustl.edu) - Database of curated mutations (v3.2, Apr 2016)
* [CGI](http://www.cancergenomeinterpreter.org/biomarkers) - Cancer Biomarkers database (Jan 17th 2018)
* [DisGeNET](http://www.disgenet.org) - Database of gene-tumor type associations (v7.0, May 2020)
* [Cancer Hotspots](http://cancerhotspots.org) - Resource for statistically significant mutations in cancer (v2 - 2017)
* [dBNSFP](https://sites.google.com/site/jpopgen/dbNSFP) - Database of non-synonymous functional predictions (v4.1, June 2020)
* [TCGA](https://portal.gdc.cancer.gov/) - somatic mutations discovered across 33 tumor type cohorts (The Cancer Genome Atlas (TCGA), release 25, July 2020)
* [CHASMplus](https://karchinlab.github.io/CHASMplus/) - predicted driver mutations across 33 tumor type cohorts in TCGA
* [UniProt/SwissProt KnowledgeBase](http://www.uniprot.org) - Resource on protein sequence and functional information (2020_04, August 2020)
* [Pfam](http://pfam.xfam.org) - Database of protein families and domains (v33, May 2020)
* [Open Targets Platform](https://targetvalidation.org) - Target-disease and target-drug associations  (2020_04, June 2020)
* [ChEMBL](https://www.ebi.ac.uk/chembl/) - Manually curated database of bioactive molecules (v27, May 2020)
* [CancerMine](https://zenodo.org/record/3472758#.XZjCqeczaL4) - Literature-mined database of tumor suppressor genes/proto-oncogenes (v28, September 2020)


### Getting started

#### STEP 0: Python

An installation of Python (version _3.6_) is required to run PCGR. Check that Python is installed by typing `python --version` in your terminal window. In addition, a [Python library](https://github.com/uiri/toml) for parsing configuration files encoded with [TOML](https://github.com/toml-lang/toml) is needed. To install, simply run the following command:

   	pip install toml

#### STEP 1: Installation of Docker

1. [Install the Docker engine](https://docs.docker.com/engine/installation/) on your preferred platform
   - installing [Docker on Linux](https://docs.docker.com/engine/installation/linux/)
   - installing [Docker on Mac OS](https://docs.docker.com/engine/installation/mac/)
   - NOTE: We have not yet been able to perform enough testing on the Windows platform, and we have received feedback that particular versions of Docker/Windows do not work with PCGR (an example being [mounting of data volumes](https://github.com/docker/toolbox/issues/607))
2. Test that Docker is running, e.g. by typing `docker ps` or `docker images` in the terminal window
3. Adjust the computing resources dedicated to the Docker, i.e.:
   - Memory: minimum 5GB
   - CPUs: minimum 4
   - [How to - Mac OS X](https://docs.docker.com/docker-for-mac/#advanced)

#### STEP 2: Download PCGR and data bundle

##### Development version

a. Clone the PCGR GitHub repository (includes run script and default configuration file): `git clone https://github.com/sigven/pcgr.git`

b. Download and unpack the latest data bundles in the PCGR directory
   * [grch37 data bundle - 20200920](http://insilico.hpc.uio.no/pcgr/pcgr.databundle.grch37.20200920.tgz) (approx 17Gb)
   * [grch38 data bundle - 20200920](http://insilico.hpc.uio.no/pcgr/pcgr.databundle.grch38.20200920.tgz) (approx 18Gb)
   * *Unpacking*: `gzip -dc pcgr.databundle.grch37.YYYYMMDD.tgz | tar xvf -`

c. Pull the [PCGR Docker image (*dev*)](https://hub.docker.com/r/sigven/pcgr/) from DockerHub (approx 6.8Gb):
* `docker pull sigven/pcgr:dev` (PCGR annotation engine)

##### Latest release

a. Download and unpack the [latest software release (0.9.0rc)](https://github.com/sigven/pcgr/releases/tag/v0.9.0rc)

b. Download and unpack the assembly-specific data bundle in the PCGR directory
  * [grch37 data bundle - 20200920](http://insilico.hpc.uio.no/pcgr/pcgr.databundle.grch37.20200920.tgz) (approx 17Gb)
  * [grch38 data bundle - 20200920](http://insilico.hpc.uio.no/pcgr/pcgr.databundle.grch38.20200920.tgz) (approx 18Gb)
     * *Unpacking*: `gzip -dc pcgr.databundle.grch37.YYYYMMDD.tgz | tar xvf -`

    A _data/_ folder within the _pcgr-X.X_ software folder should now have been produced

c. Pull the [PCGR Docker image (0.9.0rc)](https://hub.docker.com/r/sigven/pcgr/) from DockerHub (approx 6.8Gb):
   * `docker pull sigven/pcgr:0.9.0rc` (PCGR annotation engine)

#### STEP 3: Input preprocessing

The PCGR workflow accepts two types of input files:

  * An unannotated, single-sample VCF file (>= v4.2) with called somatic variants (SNVs/InDels)
  * A copy number segment file

PCGR can be run with either or both of the two input files present.

* We __strongly__ recommend that the input VCF is compressed and indexed using [bgzip](http://www.htslib.org/doc/tabix.html) and [tabix](http://www.htslib.org/doc/tabix.html)
* If the input VCF contains multi-allelic sites, these will be subject to [decomposition](http://genome.sph.umich.edu/wiki/Vt#Decompose)
* Variants used for reporting should be designated as 'PASS' in the VCF FILTER column

The tab-separated values file with copy number aberrations __MUST__ contain the following four columns:
* Chromosome
* Start
* End
* Segment_Mean

Here, _Chromosome_, _Start_, and _End_ denote the chromosomal segment, and __Segment_Mean__ denotes the log(2) ratio for a particular segment, which is a common output of somatic copy number alteration callers. Note that coordinates must be **one-based** (i.e. chromosomes start at 1, not 0). Below shows the initial part of a copy number segment file that is formatted correctly according to PCGR's requirements:

    Chromosome	Start	End	Segment_Mean
    1 3218329 3550598 0.0024
    1 3552451 4593614 0.1995
    1 4593663 6433129 -1.0277


#### STEP 4: Configure your PCGR workflow

The PCGR software bundle comes with a default configuration file in the *conf/* folder, to be used as a starting point for runnning the PCGR workflow. The configuration file, formatted using [TOML](https://github.com/toml-lang/toml), enables the user to configure a number of options related to the following:

* __IMPORTANT__: Designation of VCF INFO tags that denote variant depth/allelic fraction
    * If this is not available/properly set, the report contents will be less informative __AND__ it will not be possible to preset thresholds for variant depth/allelic fraction
* Tumor-only analysis options
	* tick on/off various filtering schemes for exclusion of germline variants (for input VCFs coming from tumor-only sequencing assays)
* VEP/_vcfanno_ options

See here for more details about the exact [usage of the configuration options](http://pcgr.readthedocs.io/en/latest/input.html#pcgr-configuration-file).

#### STEP 5: Run example

A tumor sample report is generated by calling the Python script __pcgr.py__, which takes the following arguments and options:

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

The _examples_ folder contain input VCF files from two tumor samples sequenced within TCGA (**GRCh37** only). It also contains a PCGR configuration file customized for these VCFs. A report for a colorectal tumor case can be generated by running the following command in your terminal window:

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

This command will run the Docker-based PCGR workflow and produce the following output files in the _examples_ folder:

  1. __tumor_sample.COAD.pcgr_acmg.grch37.html__ - An interactive HTML report for clinical interpretation
  2. __tumor_sample.COAD.pcgr_acmg.grch37.pass.vcf.gz__ - Bgzipped VCF file with rich set of annotations for precision oncology
  3. __tumor_sample.COAD.pcgr_acmg.grch37.pass.tsv.gz__ - Compressed vcf2tsv-converted file with rich set of annotations for precision oncology
  4. __tumor_sample.COAD.pcgr_acmg.grch37.snvs_indels.tiers.tsv__ - Tab-separated values file with variants organized according to tiers of functional relevance
  5. __tumor_sample.COAD.pcgr_acmg.grch37.mutational_signatures.tsv__ - Tab-separated values file with information on contribution of mutational signatures
  5. __tumor_sample.COAD.pcgr_acmg.grch37.json.gz__ - Compressed JSON dump of HTML report content
  6. __tumor_sample.COAD.pcgr_acmg.grch37.cna_segments.tsv.gz__ - Compressed tab-separated values file with annotations of gene transcripts that overlap with somatic copy number aberrations

  ## FAQ

  Frequently asked questions regarding PCGR usage and functionality:

  __1. Why do I get a long list of lines with “ERROR: Line ..” during “STEP 0: validate input data”?__

  _Answer: Your query VCF does not pass the VCF validation check by_ [EBI's vcf-validator](https://github.com/EBIvariation/vcf-validator).
  _Solution: 1) Fix the VCF so that it adheres to the VCF standard, or 2) run PCGR with option `--vcf_no_validate` if you think the formatting problems is not critical to the contents of the VCF._

  __2. I am not sure how to specify depth/allelic fraction in my query VCF. Why cannot PCGR pull out this information automatically from my VCF file?__

  _Answer: This is something that you as a user need to handle yourself. To our knowledge, there is currently no standard way that variant callers format these types of data (allelic fraction/depth, tumor/normal) in the VCF, and this makes it very challenging for PCGR to automatically grab this information from the variety of VCFs produced by different variant callers. Please take a careful look at the example VCF files (`examples` folder) that comes with PCGR for how PCGR expects this information to be formatted, and make sure your VCF is formatted accordingly._

  __3. Is it possible to utilize PCGR for analysis of multiple samples?__

  _Answer: As the name of the tool implies, PCGR was developed for the detailed analysis of individual tumor samples. However, if you take advantage of the different outputs from PCGR, it can also be utilized for analysis of multiple samples. First, make sure your input files are organized per sample (i.e. one VCF file per sample, one CNA file per sample), so that they can be fed directly to PCGR. Now, once all samples have been processed with PCGR, note that all the tab-separated output files (i.e. tiers, mutational signatures, cna segments) contain the sample identifier, which enable them to be aggregated and suitable for a downstream multi-sample analysis._

  _Also note that the compressed JSON output pr. sample run contains __ALL__ information presented in the report. Explore the JSON contents e.g. with the_ [jsonlite package](https://github.com/jeroen/jsonlite) _in R:_

   `report_data <- jsonlite::fromJSON('<sample_id>.pcgr_acmg.grch37.json.gz')`

  _E.g. tiered SNV/InDel output_:

   `head(report_data$content$snv_indel$variant_set$tsv)`

  _Or TMB estimate_:

   `report_data$content$tmb$variant_statistic$tmb_estimate`

  __4. I do not see the expected transcript-specific consequence for a particular variant. In what way is the primary variant consequence established?__

  _Answer: PCGR relies upon_ [VEP](https://www.ensembl.org/info/docs/tools/vep/index.html)  _for consequence prioritization, in which a specific transcript-specific consequence is chosen as the primary variant consequence. In the PCGR configuration file, you may customise how this is chosen by changing the order of criteria applied when choosing a primary consequence block  - parameter_ [vep_pick_order](https://www.ensembl.org/info/docs/tools/vep/script/vep_other.html#pick_options)

  __5. Is it possible to use RefSeq as the underlying gene transcript model in PCGR?__

  _Answer: PCGR uses GENCODE as the primary gene transcript model, but we provide cross-references to corresponding RefSeq transcripts when this is available._

  __6. I have a VCF with structural variants detected in my tumor sample, can PCGR process those as well?__

  _Answer: This is currently not supported as input for PCGR, but is something we want to incorporate in the future._

  __7. I am surprised to see a particular gene being located in TIER 3 for my sample, since I know that this gene is of potential clinical significance in the tumor type I am investigating?__

  _Answer: PCGR classifies variants into tiers of significance through an implementation of_  [published guidelines by ACMG/AMP](https://pcgr.readthedocs.io/en/latest/tier_systems.html). _No manual efforts for individual tumor types are conducted beyond this rule-based scheme. The users need to keep this in mind when interpreting the tier contents of the report._

  __8. Is it possible to see all the invididual cancer subtypes that belong to each of the 30 different tumor sites?__

  _Answer: Yes, see_ [oncotree_ontology_xref.tsv](https://raw.githubusercontent.com/sigven/pcgr/master/oncotree_ontology_xref.tsv)

  __9. Is there any plans to incorporate data from__ [OncoKB](https://www.oncokb.org) __in PCGR?__

  _Answer: No. PCGR relies upon publicly available open-source resources, and further that the PCGR data bundle can be distributed freely to the user community. It is our understanding that_ [OncoKB's terms of use](https://www.oncokb.org/terms) _do not fit well with this strategy._

  __10. Is it possible for the users to update the data bundle to get the most recent versions of all underlying data sources?__

  _Answer: As of now, the data bundle is updated only with each release of PCGR. This is not yet supported, but we want to allow for such updates in the future._



## Contact

sigven AT ifi.uio.no
