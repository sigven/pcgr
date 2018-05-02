
## CHANGELOG

#### 0.6.1 - May 2nd 2018

##### Fixed
 * Bug in tier assignment 'pcgr_acmg' (case for no variants in tier1,2,3)
 * Bug in tier assignment 'pcgr_acmg' (no tumor type specified, evidence items with weak support detected)
 * Bug: duplicated variants in 'Tier 3' resulting from genes encoded with dual roles as tumor suppressor genes/oncogenes
 * Bug: duplicated variants in 'Tier 1/Noncoding variants' resulting from rare cases of noncoding variants occurring in Tier 1 (synonymous variants with biomarker role)

#### 0.6.0 - April 25th 2018

##### Added
 * New argument in pcgr.py
	 * *assembly* (grch37/grch38)
 * New option in pcgr.py
	 * *--basic* - run comprehensive VCF annotation only, skip report generation and additional analyses
 * New sections in HTML report
	  * *Settings and annotation sources* - now also listing key PCGR configuration settings
	  * *Main findings* - Six value boxes indicating the main findings of clinical relevance
 * New configuration options
	 * \[tier_model\](string) - choice between *pcgr_acmg* and *pcgr*
	 * \[mutational_burden\] - set TMB tertile limits
		 * *tmb_low_limit (float)*
		 * *tmb_intermediate_limit (float)*
	 * \[tumor_type\] - choose between 34 tumor types/classes:
		 * *Adrenal_Gland_Cancer_NOS (logical)*
		 * *Ampullary_Carcinoma_NOS (logical)*
		 * *Biliary_Tract_Cancer_NOS (logical)*
		 * *Bladder_Urinary_Tract_Cancer_NOS (logical)*
		 * *Blood_Cancer_NOS (logical)*
		 * *Bone_Cancer_NOS (logical)*
		 * *Breast_Cancer_NOS (logical)*
		 * *CNS_Brain_Cancer_NOS (logical)*
		 * *Colorectal_Cancer_NOS (logical)*
		 * *Cervical_Cancer_NOS (logical)*
		 * *Esophageal_Stomach_Cancer_NOS (logical)*
		 * *Head_And_Neck_Cancer_NOS (logical)*
		 * *Hereditary_Cancer_NOS (logical)*
		 * *Kidney_Cancer_NOS (logical)*
		 * *Leukemia_NOS (logical)*
		 * *Liver_Cancer_NOS (logical)*
		 * *Lung_Cancer_NOS (logical)*
		 * *Lymphoma_Hodgkin_NOS (logical)*
		 * *Lymphoma_Non_Hodgkin_NOS (logical)*
		 * *Ovarian_Fallopian_Tube_Cancer_NOS (logical)*
		 * *Pancreatic_Cancer_NOS (logical)*
		 * *Penile_Cancer_NOS (logical)*
		 * *Peripheral_Nervous_System_Cancer_NOS (logical)*
		 * *Peritoneal_Cancer_NOS (logical)*
		 * *Pleural_Cancer_NOS (logical)*
		 * *Prostate_Cancer_NOS (logical)*
		 * *Skin_Cancer_NOS (logical)*
		 * *Soft_Tissue_Cancer_NOS (logical)*
		 * *Stomach_Cancer_NOS (logical)*
		 * *Testicular_Cancer_NOS (logical)*
		 * *Thymic_Cancer_NOS (logical)*
		 * *Thyroid_Cancer_NOS (logical)*
		 * *Uterine_Cancer_NOS (logical)*
		 * *Vulvar_Vaginal_Cancer_NOS (logical)*
	 * \[mutational_signatures\]
	    * *mutsignatures_cutoff (float)* - discard any signature contributions with a weight less than the cutoff
	 * \[cna\]
	    * *transcript_cna_overlap (float)* - minimum percent overlap between copy number segment and transcripts (average) for tumor suppressor gene/proto-oncogene to be reported
	 * \[allelic_support\]
		 * If input VCF has correctly formatted depth/allelic fraction as INFO tags, users can add thresholds on depth/support that are applied prior to report generation
			 * *tumor_dp_min (integer)* - minimum sequencing depth for variant in tumor sample
			 * *tumor_af_min (float)* - minimum allelic fraction for variant in tumor sample
			 * *normal_dp_min (integer)* - minimum sequencing depth for variant in normal sample
			 * *normal_af_max (float)* - maximum allelic fraction for variant in normal sample
	 * \[visual\]
	     * *report_theme (string)* - visual theme of report (Bootstrap)
	 * \[other\]
		 * *vcf_validation (logical)* - keep/skip VCF validation by [vcf-validator](https://github.com/EBIvariation/vcf-validator)
 * New output file - JSON output of HTML report content
 * New INFO tags of PCGR-annotated VCF
	 * *CANCER_PREDISPOSITION*
	 * *PFAM_DOMAIN*
	 * *TCGA_FREQUENCY*
	 * *TCGA_PANCANCER_COUNT*
	 * *ICGC_PCAWG_OCCURRENCE*
	 * *ICGC_PCAWG_AFFECTED_DONORS*
	 * *CLINVAR_MEDGEN_CUI*
 * New column entries in annotated SNV/InDel TSV file:
	 * *CANCER_PREDISPOSITION*
	 * *ICGC_PCAWG_OCCURRENCE*
	 * *TCGA_FREQUENCY*
 * New column in CNA output
	* *TRANSCRIPTS* - aberration-overlapping transcripts (Ensembl transcript IDs)
	* *MEAN_TRANSCRIPT_CNA_OVERLAP* - Mean overlap (%) betweeen gene transcripts and aberration segment


##### Removed

 * Elements of databundle (now annotated directly through VEP):
	 * dbsnp
	 * gnomad/exac
	 * 1000G project
 * INFO tags of PCGR-annotated VCF
	 * *DBSNPBUILDID*
	 * *DBSNP_VALIDATION*
	 * *DBSNP_SUBMISSIONS*
	 * *DBSNP_MAPPINGSTATUS*
	 * *GWAS_CATALOG_PMID*
	 * *GWAS_CATALOG_TRAIT_URI*
	 * *DOCM_DISEASE*
 * Output files
	 * TSV files with mutational signature results and biomarkers (i.e. *sample_id.pcgr.snvs_indels.biomarkers.tsv* and *sample_id.pcgr.mutational_signatures.tsv*)
	  	* Data can still be retrieved - now from the JSON dump
	 * MAF file
		 * The previous MAF output was generated in a custom fashion, a more accurate MAF output based on https://github.com/mskcc/vcf2maf will be incorporated in the next release


##### Changed
 * HTML report sections
	 * *Tier statistics* and *Variant statistics* are now grouped into the section *Tier and variant statistics*
	 * *Tier 5* is now *Noncoding mutations* (i.e. not considered a tier per se)
	 * Sliders for allelic fraction in the *Global variant browser* are now fixed from 0 to 1 (0.05 intervals)
