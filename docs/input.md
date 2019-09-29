## Input

The PCGR workflow accepts two types of input files:

  * An unannotated, single-sample VCF file (>= v4.2) with called somatic variants (SNVs/InDels)
  * A copy number segment file


PCGR can be run with either or both of the two input files present.

### VCF

* We __strongly__ recommend that the input VCF is compressed and indexed using [bgzip](http://www.htslib.org/doc/tabix.html) and [tabix](http://www.htslib.org/doc/tabix.html)
* If the input VCF contains multi-allelic sites, these will be subject to [decomposition](http://genome.sph.umich.edu/wiki/Vt#Decompose)
* Variants used for reporting should be designated as 'PASS' in the VCF FILTER column

__IMPORTANT NOTE 1__: Considering the VCF output for the [numerous somatic SNV/InDel callers](https://www.biostars.org/p/19104/) that have been developed, we have a experienced a general lack of uniformity and robustness for the representation of somatic variant genotype data (e.g. variant allelic depths (tumor/normal), genotype quality etc.). Variant genotype data found as INFO tags in the input VCF can be specified as optional arguments to the PCGR workflow, which in turn can be used for interactive exploration in the tumor report.

__IMPORTANT NOTE 2__: PCGR generates a number of VCF INFO annotation tags that is appended to the query VCF. We will therefore encourage the users to submit query VCF files that have not been subject to annotations by other means, but rather a VCF file that comes directly from variant calling. If not, there are likely to be INFO tags in the query VCF file that coincide with those produced by PCGR.

### Copy number segments

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

### PCGR configuration file

The cancer genome sequencing report can be flexibly configured in a TOML-formatted configuration file. The default TOML configuration file, with descriptive comments wrt. usage are shown below:

	# PCGR configuration options (TOML).

	[tumor_only]
	## Variant filtering applied for option --tumor_only = true in pcgr.py
	## Several filters can be configured, all as a means to minimize the proportion of germline calls in the raw set derived from tumor-only calling

	## if vcf_tumor_only = true, exclude variants (SNVs/InDels) with minor allele frequency above the following population-specific thresholds
	## 1000 Genomes Project - WGS data
	maf_onekg_eur = 0.002
	maf_onekg_amr = 0.002
	maf_onekg_afr = 0.002
	maf_onekg_sas = 0.002
	maf_onekg_eas = 0.002
	maf_onekg_global = 0.002

	## exclude variants with minor allele frequency above the following population-specific thresholds
	## gnomAD - WES data
	maf_gnomad_nfe = 0.002
	maf_gnomad_amr = 0.002
	maf_gnomad_afr = 0.002
	maf_gnomad_sas = 0.002
	maf_gnomad_eas = 0.002
	maf_gnomad_fin = 0.002
	maf_gnomad_oth = 0.002
	maf_gnomad_global = 0.002

	## Exclude variants occurring in PoN (panel of normals, if provided as VCF)
	exclude_pon = true

	## Exclude likely homozygous germline variants (100% allelic fraction for alternate allele in tumor, very unlikely somatic event)
	exclude_likely_hom_germline = false

	## Exclude likely heterozygous germline variants
	## Must satisfy i) 40-60 % allelic fraction for alternate allele in tumor sample, ii) present in dbSNP + gnomAD, ii) not existing as somatic event in COSMIC/TCGA
	## Note that the application of this filter may be suboptimal for very impure tumors or variants affected by CNAs etc (under these circumstances, the allelic fraction
	## will be skewed (see e.g. discussion in PMID:29249243)
	exclude_likely_het_germline = false

	## Exclude variants found in dbSNP (only those that are NOT found in ClinVar(somatic origin)/DoCM/TCGA/COSMIC)
	exclude_dbsnp_nonsomatic = false

	## exclude all non-exonic variants
	exclude_nonexonic = true


	[allelic_support]
	## Specify INFO tags in input VCF that denotes depth/allelic fraction in tumor and normal sample
	## An additional tag that denotes call confidence (call_conf_tag) can also be specified, which will
	## be used for exploration in the global variant browser. Note that 'tumor_dp_tag' must be of
	## Type=Integer, and 'tumor_af_tag' must be of Type=Float (similarly for normal sample)
	tumor_dp_tag = ""
	tumor_af_tag = ""
	normal_dp_tag = ""
	normal_af_tag = ""
	call_conf_tag = ""

	## set thresholds for tumor/normal depth/allelic fraction, will be applied before report generation
	## requires that 'tumor/normal_dp_tag' and 'tumor/normal_af_tag' are specified above
	tumor_dp_min = 0
	tumor_af_min = 0.0
	normal_dp_min = 0
	normal_af_max = 1.0

	[mutational_burden]
	## Calculate mutational burden (similar to Chalmers et al., Genome Med, 2017)
	mutational_burden = true
	## Size of coding target region in megabases (defaults to exome ~ 36 Mb)
	## Note: this should ideally denote the callable target size (i.e. reflecting variable
	## sequencing depth)
	target_size_mb = 34.0
	## set upper limits to tumor mutational burden tertiles (mutations/Mb)
	tmb_low_limit = 5
	tmb_intermediate_limit = 20
	## tmb_high = tmb > tmb_intermediate_limit

	[cna]
	## log ratio thresholds for determination of copy number gains and homozygous deletions
	logR_gain = 0.8
	logR_homdel = -0.8

	## percent overlap between copy number segment and gene transcripts (average) for reporting of gains/losses in tumor suppressor genes/oncogenes
	cna_overlap_pct = 100

	[msi]
	## Predict microsatellite instability
	msi = true

	[mutational_signatures]
	## Identify relative contribution of 30 known mutational signatures (COSMIC) through the deconstructSigs framework
	mutsignatures = true
	## deconstructSigs option: number of mutational signatures to limit the search to ('signatures.limit' in whichSignatures)
	mutsignatures_signature_limit = 6
	## deconstructSigs option: type of trimer count normalization for inference of known mutational signatures, see explanation at https://github.com/raerose01/deconstructSigs"
	## options = 'default', 'exome', 'genome', 'exome2genome'
	## NOTE: If your data (VCF) is from exome sequencing, 'default' or 'exome2genome' should be used. See https://github.com/raerose01/deconstructSigs/issues/2
	mutsignatures_normalization = "exome2genome"
	## Require a minimum number of mutations for signature estimation
	mutsignatures_mutation_limit = 100
	## deconstructSigs option: discard any signature contributions with a weight less than this amount
	mutsignatures_cutoff = 0.06

	[visual]
	## Choose visual theme of report, any of: "default", "cerulean", "journal", "flatly", "readable", "spacelab", "united", "cosmo", "lumen", "paper", "sandstone", "simplex", or "yeti" (https://bootswatch.com/)
	report_theme = "default"

	[custom_tags]
	## list VCF info tags that should be present in JSON and TSV output
	## tags should be comma separated, i.e. custom_tags = "MUTECT2_FILTER,STRELKA_FILTER"
	custom_tags = ""

	[other]
	## list/do not list noncoding variants
	list_noncoding = true
	## VEP/vcfanno processing options
	n_vcfanno_proc = 4
	n_vep_forks = 4
	## Customise the order of criteria used to pick the primary transcript in VEP (see https://www.ensembl.org/info/docs/tools/vep/script/vep_options.html#opt_pick_order)
	vep_pick_order = "canonical,appris,biotype,ccds,rank,tsl,length"
	## omit intergenic variants during VEP processing
	vep_skip_intergenic = false
	## generate a MAF for input VCF using https://github.com/mskcc/vcf2maf
	vcf2maf = true
