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

Here, _Chromosome_, _Start_, and _End_ denote the chromosomal segment, and __Segment_Mean__ denotes the log(2) ratio for a particular segment, which is a common output of somatic copy number alteration callers. Below shows the initial part of a copy number segment file that is formatted correctly according to PCGR's requirements:

      Chromosome	Start	End	Segment_Mean
      1 3218329 3550598 0.0024
      1 3552451 4593614 0.1995
      1 4593663 6433129 -1.0277

### PCGR configuration file

The cancer genome sequencing report can be flexibly configured in a TOML-formatted configuration file. The default TOML configuration file, with descriptive comments wrt. usage are shown below:

	# PCGR configuration options (TOML).

	[tier_model]
	## tier model for prioritization of SNVs/InDels ("pcgr_acmg" or "pcgr")
	tier_model = "pcgr_acmg"

	[tumor_only]
	## If input VCF contains mix of germline/somatic (variants called with no matching control, i.e. tumor-only) set vcf_tumor_only to true
	vcf_tumor_only = false

	## if vcf_tumor_only = true, exclude variants (SNVs/InDels) with minor allele frequency above the following population-specific thresholds
	## 1000 Genomes Project - WGS data
	maf_onekg_eur = 0.01
	maf_onekg_amr = 0.01
	maf_onekg_afr = 0.01
	maf_onekg_sas = 0.01
	maf_onekg_eas = 0.01
	maf_onekg_global = 0.01

	## exclude variants with minor allele frequency above the following population-specific thresholds
	## gnomAD - WES data
	maf_gnomad_nfe = 0.01
	maf_gnomad_amr = 0.01
	maf_gnomad_afr = 0.01
	maf_gnomad_sas = 0.01
	maf_gnomad_eas = 0.01
	maf_gnomad_fin = 0.01
	maf_gnomad_oth = 0.01
	maf_gnomad_global = 0.01

	## exclude variants found in dbSNP (only those not found as somatic in ClinVar/Docm)
	exclude_dbsnp_nonclinical = true

	## in variant exclusion from dbSNP, set whether found in TCGA should be kept (at desired recurrence level)
	## E.g. keep_known_tcga = true + tcga_recurrence = 2 keeps all TCGA variants (that intersect dbSNP) found in at least two samples
	keep_known_tcga = true
	tcga_recurrence = 2

	## exclude all non protein-coding variants
	exclude_noncoding = true

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
	target_size_mb = 36.0
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
	mutsignatures_normalization = "default"
	## Require a minimum number of mutations for signature estimation
	mutsignatures_mutation_limit = 100
	## deconstructSigs option: discard any signature contributions with a weight less than this amount
	mutsignatures_cutoff = 0.06

	[tumor_type]
	## Choose tumor type/class of input sample
	## Due to partial overlap between some classes, user can set maximum two types
	Adrenal_Gland_Cancer_NOS = false
	Ampullary_Carcinoma_NOS = false
	Biliary_Tract_Cancer_NOS = false
	Bladder_Urinary_Tract_Cancer_NOS = false
	Blood_Cancer_NOS = false
	Bone_Cancer_NOS = false
	Breast_Cancer_NOS = false
	CNS_Brain_Cancer_NOS = false
	Colorectal_Cancer_NOS = true
	Cervical_Cancer_NOS = false
	DNA_Repair_Deficiency_Disorders = false
	Esophageal_Stomach_Cancer_NOS = false
	Head_And_Neck_Cancer_NOS = false
	Hereditary_Cancer_NOS = false
	Kidney_Cancer_NOS = false
	Leukemia_NOS = false
	Liver_Cancer_NOS = false
	Lung_Cancer_NOS = false
	Lymphoma_Hodgkin_NOS = false
	Lymphoma_Non_Hodgkin_NOS = false
	Mesothelioma = false
	Multiple_Myeloma = false
	Ovarian_Fallopian_Tube_Cancer_NOS = false
	Pancreatic_Cancer_NOS = false
	Penile_Cancer_NOS = false
	Peripheral_Nervous_System_Cancer_NOS = false
	Peritoneal_Cancer_NOS = false
	Pleural_Cancer_NOS = false
	Prostate_Cancer_NOS = false
	Skin_Cancer_NOS = false
	Soft_Tissue_Cancer_NOS = false
	Stomach_Cancer_NOS = false
	Testicular_Cancer_NOS = false
	Thymic_Cancer_NOS = false
	Thyroid_Cancer_NOS = false
	Uterine_Cancer_NOS = false
	Vulvar_Vaginal_Cancer_NOS = false

	[visual]
	## Choose visual theme of report, any of: "default", "cerulean", "journal", "flatly", "readable", "spacelab", "united", "cosmo", "lumen", "paper", "sandstone", "simplex", or "yeti" (https://bootswatch.com/)
	report_theme = "default"

	[other]
	## Keep/skip VCF validation by https://github.com/EBIvariation/vcf-validator
	## The vcf-validator checks that the input VCF is properly encoded. Since the
	## vcf-validator is strict, and with error messages that is not always self-explanatory,
	## the users can skip validation if they are confident that the most critical parts of the VCF
	## are properly encoded
	vcf_validation = true
	## list/do not list noncoding variants
	list_noncoding = true
	## VEP/vcfanno processing options
	n_vcfanno_proc = 4
	n_vep_forks = 4
	## omit intergenic variants during VEP processing
	vep_skip_intergenic = false
