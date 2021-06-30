#!/usr/bin/env python

import csv
import re
import argparse
import os
import subprocess
import logging
import sys
import getpass
import platform
from argparse import RawTextHelpFormatter

PCGR_VERSION = "0.9.2"
DB_VERSION = "PCGR_DB_VERSION = 20210627"
VEP_VERSION = "104"
GENCODE_VERSION = "38"
NCBI_BUILD_MAF = "GRCh38"
VEP_ASSEMBLY = "GRCh38"
DOCKER_IMAGE_VERSION = "sigven/pcgr:" + str(PCGR_VERSION)

global debug

tsites = {
		0: "Any",
      1: "Adrenal Gland",
      2: "Ampulla of Vater",
      3: "Biliary Tract",
      4: "Bladder/Urinary Tract",
      5: "Bone",
      6: "Breast",
      7: "Cervix",
      8: "CNS/Brain",
      9: "Colon/Rectum",
      10: "Esophagus/Stomach",
      11: "Eye",
      12: "Head and Neck",
      13: "Kidney",
      14: "Liver",
      15: "Lung",
      16: "Lymphoid",
      17: "Myeloid",
      18: "Ovary/Fallopian Tube",
      19: "Pancreas",
      20: "Peripheral Nervous System",
      21: "Peritoneum",
      22: "Pleura",
      23: "Prostate",
      24: "Skin",
      25: "Soft Tissue",
      26: "Testis",
      27: "Thymus",
      28: "Thyroid",
      29: "Uterus",
      30: "Vulva/Vagina"
	}

def __main__():

   tumor_sites = "1 = Adrenal Gland\n"
   tumor_sites = tumor_sites + "2 = Ampulla of Vater\n"
   tumor_sites = tumor_sites + "3 = Biliary Tract\n"
   tumor_sites = tumor_sites + "4 = Bladder/Urinary Tract\n"
   tumor_sites = tumor_sites + "5 = Bone\n"
   tumor_sites = tumor_sites + "6 = Breast\n"
   tumor_sites = tumor_sites + "7 = Cervix\n"
   tumor_sites = tumor_sites + "8 = CNS/Brain\n"
   tumor_sites = tumor_sites + "9 = Colon/Rectum\n"
   tumor_sites = tumor_sites + "10 = Esophagus/Stomach\n"
   tumor_sites = tumor_sites + "11 = Eye\n"
   tumor_sites = tumor_sites + "12 = Head and Neck\n"
   tumor_sites = tumor_sites + "13 = Kidney\n"
   tumor_sites = tumor_sites + "14 = Liver\n"
   tumor_sites = tumor_sites + "15 = Lung\n"
   tumor_sites = tumor_sites + "16 = Lymphoid\n"
   tumor_sites = tumor_sites + "17 = Myeloid\n"
   tumor_sites = tumor_sites + "18 = Ovary/Fallopian Tube\n"
   tumor_sites = tumor_sites + "19 = Pancreas\n"
   tumor_sites = tumor_sites + "20 = Peripheral Nervous System\n"
   tumor_sites = tumor_sites + "21 = Peritoneum\n"
   tumor_sites = tumor_sites + "22 = Pleura\n"
   tumor_sites = tumor_sites + "23 = Prostate\n"
   tumor_sites = tumor_sites + "24 = Skin\n"
   tumor_sites = tumor_sites + "25 = Soft Tissue\n"
   tumor_sites = tumor_sites + "26 = Testis\n"
   tumor_sites = tumor_sites + "27 = Thymus\n"
   tumor_sites = tumor_sites + "28 = Thyroid\n"
   tumor_sites = tumor_sites + "29 = Uterus\n"
   tumor_sites = tumor_sites + "30 = Vulva/Vagina\n"

   program_description = "Personal Cancer Genome Reporter (PCGR) workflow for clinical interpretation of " + \
      "somatic nucleotide variants and copy number aberration segments"
   program_options = "\n\t--input_vcf <INPUT_VCF>\n\t--pcgr_dir <PCGR_DIR>\n\t--output_dir <OUTPUT_DIR>\n\t--genome_assembly" + \
      " <GENOME_ASSEMBLY>\n\t--sample_id <SAMPLE_ID>"

   parser = argparse.ArgumentParser(description = program_description,
                                    formatter_class=RawTextHelpFormatter, usage="\n\t%(prog)s -h [options] " + str(program_options) + "\n \n")
   parser._action_groups.pop()
   required = parser.add_argument_group("Required arguments")
   #optional_rna = parser.add_argument_group("Bulk RNA-seq and RNA fusion options")
   optional_vcfanno = parser.add_argument_group("vcfanno options")
   optional_vep = parser.add_argument_group("VEP options")
   optional_tmb_msi = parser.add_argument_group("Tumor mutational burden (TMB) and MSI options")
   optional_signatures = parser.add_argument_group("Mutatonal signature options")
   optional_tumor_only = parser.add_argument_group("Tumor-only options")
   optional_allelic_support = parser.add_argument_group("Allelic support options")
   optional_other = parser.add_argument_group("Other options")

   #optional_rna.add_argument("--rna_fusion", dest = "rna_fusion_tumor", help = "File with RNA fusion transcripts detected in tumor (tab-separated values)")
   #optional_rna.add_argument("--rna_expression", dest = "rna_exp_tumor", help = "File with RNA expression levels (bulk) and differential expression status of genes in tumor (tab-separated values)")

   optional_other.add_argument("--input_cna", dest = "input_cna", help = "Somatic copy number alteration segments (tab-separated values)")
   optional_other.add_argument("--logr_gain", type = float, default = 0.8, dest = "logr_gain",help = "Log ratio-threshold for regions containing copy number gains/amplifications (default: %(default)s)")
   optional_other.add_argument("--logr_homdel", type = float, default = -0.8, dest = "logr_homdel",help = "Log ratio-threshold for regions containing homozygous deletions (default: %(default)s)")
   optional_other.add_argument("--cna_overlap_pct", type = float, default = 50, dest = "cna_overlap_pct", help = "Mean percent overlap between copy number segment and gene transcripts for reporting of gains/losses in tumor suppressor genes/oncogenes, (default: %(default)s)")
   optional_other.add_argument("--tumor_site", dest = "tsite", type = int, default = 0, help = "Optional integer code to specify primary tumor type/site of query sample,\n choose any of the following identifiers:\n" + str(tumor_sites) + "(default: %(default)s - any tumor type)")
   optional_other.add_argument("--tumor_purity", type = float, dest = "tumor_purity", help = "Estimated tumor purity (between 0 and 1, (default: %(default)s)")
   optional_other.add_argument("--tumor_ploidy", type = float, dest = "tumor_ploidy", help = "Estimated tumor ploidy (default: %(default)s)")
   
   optional_allelic_support.add_argument("--tumor_dp_tag", dest = "tumor_dp_tag", default="_NA_", help = "Specify VCF INFO tag for sequencing depth (tumor, must be Type=Integer, default: %(default)s")
   optional_allelic_support.add_argument("--tumor_af_tag", dest = "tumor_af_tag", default="_NA_", help = "Specify VCF INFO tag for variant allelic fraction (tumor,  must be Type=Float, default: %(default)s")
   optional_allelic_support.add_argument("--control_dp_tag", dest = "control_dp_tag", default="_NA_", help = "Specify VCF INFO tag for sequencing depth (control, must be Type=Integer, default: %(default)s")
   optional_allelic_support.add_argument("--control_af_tag", dest = "control_af_tag", default="_NA_", help = "Specify VCF INFO tag for variant allelic fraction (control, must be Type=Float, default: %(default)s")
   optional_allelic_support.add_argument("--call_conf_tag", dest = "call_conf_tag", default="_NA_", help = "Specify VCF INFO tag for somatic variant call confidence (must be categorical, e.g. Type=String, default: %(default)s")
   optional_allelic_support.add_argument("--tumor_dp_min", type = int, default = 0, dest = "tumor_dp_min",help = "If VCF INFO tag for sequencing depth (tumor) is specified and found, set minimum required depth for inclusion in report (default: %(default)s)")
   optional_allelic_support.add_argument("--tumor_af_min", type = float, default = 0, dest = "tumor_af_min",help = "If VCF INFO tag for variant allelic fraction (tumor) is specified and found, set minimum required AF for inclusion in report (default: %(default)s)")
   optional_allelic_support.add_argument("--control_dp_min", type = int, default = 0, dest = "control_dp_min",help = "If VCF INFO tag for sequencing depth (control) is specified and found, set minimum required depth for inclusion in report (default: %(default)s)")
   optional_allelic_support.add_argument("--control_af_max", type = float, default = 1, dest = "control_af_max",help = "If VCF INFO tag for variant allelic fraction (control) is specified and found, set maximum tolerated AF for inclusion in report (default: %(default)s)")
   
   optional_tmb_msi.add_argument("--target_size_mb", type = float, default = 34, dest = "target_size_mb",help = "For mutational burden analysis - approximate protein-coding target size in Mb of sequencing assay (default: %(default)s (WES/WGS))")
   optional_tmb_msi.add_argument("--estimate_tmb", action = "store_true", help = "Estimate tumor mutational burden from the total number of somatic mutations and target region size, default: %(default)s")
   optional_tmb_msi.add_argument("--estimate_msi_status", action = "store_true", help = "Predict microsatellite instability status from patterns of somatic mutations/indels, default: %(default)s")
   optional_tmb_msi.add_argument("--tmb_algorithm", dest = "tmb_algorithm", default = "all_coding", choices = ["all_coding","nonsyn"], help = "Method for calculation of TMB, all coding variants (Chalmers et al., Genome Medicine, 2017), or non-synonymous variants only, default: %(default)s")
   
   optional_signatures.add_argument("--estimate_signatures", action = "store_true", help = "Estimate relative contributions of reference mutational signatures in query sample and detect potential kataegis events, default: %(default)s")
   optional_signatures.add_argument("--min_mutations_signatures", type = int, default = 200, dest = "min_mutations_signatures", help = "Minimum number of SNVs required for reconstruction of mutational signatures (SBS) by MutationalPatterns (default: %(default)s, minimum n = 100)")
   optional_signatures.add_argument("--all_reference_signatures", action = "store_true", help = "Use all reference mutational signatures (SBS, n = 67) in signature reconstruction rather than only those already attributed to the tumor type (default: %(default)s)")
   optional_signatures.add_argument('--include_artefact_signatures', action = "store_true", help = "Include sequencing artefacts in the collection of reference signatures (default: %(default)s")
   
   optional_other.add_argument("--cpsr_report", dest = "cpsr_report", help = "CPSR report file (Gzipped JSON - file ending with 'cpsr.<genome_assembly>.json.gz' -  germline report of patient's blood/control sample")
   optional_other.add_argument("--vcf2maf", action = "store_true", help = "Generate a MAF file for input VCF using https://github.com/mskcc/vcf2maf (default: %(default)s)")
   optional_other.add_argument("--show_noncoding", action = "store_true", help = "List non-coding (i.e. non protein-altering) variants in report, default: %(default)s")
   optional_other.add_argument("--assay", dest = "assay", choices = ["WES","WGS","TARGETED"], default = "WES", help = "Type of DNA sequencing assay performed for input data (VCF) default: %(default)s")
   optional_other.add_argument("--include_trials",action = "store_true",help = "(Beta) Include relevant ongoing or future clinical trials, focusing on studies with molecularly targeted interventions")
   optional_other.add_argument("--preserved_info_tags", dest ="preserved_info_tags", default="None", help="Comma-separated string of VCF INFO tags from query VCF that should be kept in PCGR output TSV file")
   optional_other.add_argument("--report_theme",choices = ["default","cerulean","journal","flatly","readable","spacelab","united","cosmo","lumen","paper","sandstone","simplex","yeti"], help="Visual report theme (rmarkdown)", default = 'default')
   optional_other.add_argument('--report_nonfloating_toc', action='store_true', help='Do not float the table of contents (TOC) in output report (rmarkdown), default: %(default)s')
   optional_other.add_argument("--force_overwrite", action = "store_true", help = "By default, the script will fail with an error if any output file already exists. You can force the overwrite of existing result files by using this flag, default: %(default)s")
   optional_other.add_argument("--version", action = "version", version="%(prog)s " + str(PCGR_VERSION))
   optional_other.add_argument("--basic", action = "store_true", help = "Run functional variant annotation on VCF through VEP/vcfanno, omit other analyses (i.e. Tier assignment/MSI/TMB/Signatures etc. and report generation (STEP 4), default: %(default)s")
   optional_other.add_argument("--no_vcf_validate", action = "store_true",help = "Skip validation of input VCF with Ensembl's vcf-validator, default: %(default)s")
   optional_other.add_argument("--docker_uid", dest = "docker_user_id", help = "Docker user ID. default is the host system user ID. If you are experiencing permission errors, try setting this up to root (`--docker-uid root`)")
   optional_other.add_argument("--no_docker", action ="store_true", dest = "no_docker", default=False, help = "Run the PCGR workflow in a non-Docker mode (see install_no_docker/ folder for instructions)")   
   optional_other.add_argument("--debug", action ="store_true", default=False, help = "Print full Docker commands to log, default: %(default)s")
   

   optional_vcfanno.add_argument("--vcfanno_n_proc", default = 4, type = int, help="Number of vcfanno processes (option '-p' in vcfanno), default: %(default)s")
   optional_vep.add_argument("--vep_n_forks", default = 4, type = int, help="Number of forks (option '--fork' in VEP), default: %(default)s")   
   optional_vep.add_argument("--vep_buffer_size", default = 100, type = int, help="Variant buffer size (variants read into memory simultaneously, option '--buffer_size' in VEP) " + \
      "\n- set lower to reduce memory usage, default: %(default)s")
   optional_vep.add_argument("--vep_pick_order", default = "canonical,appris,biotype,ccds,rank,tsl,length,mane", help="Comma-separated string " + \
      "of ordered transcript/variant properties for selection of primary variant consequence\n ( option '--pick_order' in VEP), default: %(default)s")
   optional_vep.add_argument("--vep_no_intergenic", action = "store_true", help="Skip intergenic variants during processing (option '--no_intergenic' in VEP), default: %(default)s")
   optional_vep.add_argument("--vep_regulatory", action = "store_true", help="Add VEP regulatory annotations (option '--regulatory' )or non-coding interpretation, default: %(default)s")

   optional_tumor_only.add_argument("--tumor_only",action = "store_true",help = "Input VCF comes from tumor-only sequencing, calls will be filtered for variants of germline origin, (default: %(default)s)")
   optional_tumor_only.add_argument("--cell_line",action = "store_true",help = "Input VCF comes from tumor cell line sequencing (requires --tumor_only), calls will be filtered for variants of germline origin, (default: %(default)s)")
   optional_tumor_only.add_argument("--pon_vcf", dest = "pon_vcf", help = "VCF file with germline calls from Panel of Normals (PON) - i.e. blacklisted variants, (default: %(default)s)")
   optional_tumor_only.add_argument("--maf_onekg_eur", dest = "maf_onekg_eur", type = float, default = 0.002, help = "Exclude variants in tumor (SNVs/InDels, tumor-only mode) with MAF above the given percent threshold (1000 Genomes Project - European pop, default: %(default)s)")
   optional_tumor_only.add_argument("--maf_onekg_amr", dest = "maf_onekg_amr", type = float, default = 0.002, help = "Exclude variants in tumor (SNVs/InDels, tumor-only mode) with MAF above the given percent threshold (1000 Genomes Project - Ad Mixed American pop, default: %(default)s)")
   optional_tumor_only.add_argument("--maf_onekg_afr", dest = "maf_onekg_afr", type = float, default = 0.002, help = "Exclude variants in tumor (SNVs/InDels, tumor-only mode) with MAF above the given percent threshold (1000 Genomes Project - African pop, default: %(default)s)")
   optional_tumor_only.add_argument("--maf_onekg_eas", dest = "maf_onekg_eas", type = float, default = 0.002, help = "Exclude variants in tumor (SNVs/InDels, tumor-only mode) with MAF above the given percent threshold (1000 Genomes Project - East Asian pop, default: %(default)s)")
   optional_tumor_only.add_argument("--maf_onekg_sas", dest = "maf_onekg_sas", type = float, default = 0.002, help = "Exclude variants in tumor (SNVs/InDels, tumor-only mode) with MAF above the given percent threshold (1000 Genomes Project - South Asian pop, default: %(default)s)")
   optional_tumor_only.add_argument("--maf_onekg_global", dest = "maf_onekg_global", type = float, default = 0.002, help = "Exclude variants in tumor (SNVs/InDels, tumor-only mode) with MAF above the given percent threshold (1000 Genomes Project - global pop, default: %(default)s)")
   optional_tumor_only.add_argument("--maf_gnomad_nfe", dest = "maf_gnomad_nfe", type = float, default = 0.002, help = "Exclude variants in tumor (SNVs/InDels, tumor-only mode) with MAF above the given percent threshold, (gnomAD - European (non-Finnish), default: %(default)s)") 
   optional_tumor_only.add_argument("--maf_gnomad_asj", dest = "maf_gnomad_asj", type = float, default = 0.002, help = "Exclude variants in tumor (SNVs/InDels, tumor-only mode) with MAF above the given percent threshold, (gnomAD - Ashkenazi Jewish, default: %(default)s)") 
   optional_tumor_only.add_argument("--maf_gnomad_fin", dest = "maf_gnomad_fin", type = float, default = 0.002, help = "Exclude variants in tumor (SNVs/InDels, tumor-only mode) with MAF above the given percent threshold, (gnomAD - European (Finnish), default: %(default)s)") 
   optional_tumor_only.add_argument("--maf_gnomad_oth", dest = "maf_gnomad_oth", type = float, default = 0.002, help = "Exclude variants in tumor (SNVs/InDels, tumor-only mode) with MAF above the given percent threshold, (gnomAD - Other, default: %(default)s)") 
   optional_tumor_only.add_argument("--maf_gnomad_amr", dest = "maf_gnomad_amr", type = float, default = 0.002, help = "Exclude variants in tumor (SNVs/InDels, tumor-only mode) with MAF above the given percent threshold, (gnomAD - Latino/Admixed American, default: %(default)s)")
   optional_tumor_only.add_argument("--maf_gnomad_afr", dest = "maf_gnomad_afr", type = float, default = 0.002, help = "Exclude variants in tumor (SNVs/InDels, tumor-only mode) with MAF above the given percent threshold, (gnomAD - African/African-American, default: %(default)s)")
   optional_tumor_only.add_argument("--maf_gnomad_eas", dest = "maf_gnomad_eas", type = float, default = 0.002, help = "Exclude variants in tumor (SNVs/InDels, tumor-only mode) with MAF above the given percent threshold, (gnomAD - East Asian, default: %(default)s)")
   optional_tumor_only.add_argument("--maf_gnomad_sas", dest = "maf_gnomad_sas", type = float, default = 0.002, help = "Exclude variants in tumor (SNVs/InDels, tumor-only mode) with MAF above the given percent threshold, (gnomAD - South Asian, default: %(default)s)")
   optional_tumor_only.add_argument("--maf_gnomad_global", dest = "maf_gnomad_global", type = float, default = 0.002, help = "Exclude variants in tumor (SNVs/InDels, tumor-only mode) with MAF above the given percent threshold, (gnomAD - global population, default: %(default)s)")
   optional_tumor_only.add_argument("--exclude_pon", action = "store_true", help = "Exclude variants occurring in PoN (Panel of Normals, if provided as VCF (--pon_vcf), default: %(default)s)")
   optional_tumor_only.add_argument("--exclude_likely_hom_germline", action = "store_true", help = "Exclude likely homozygous germline variants (100 pct allelic fraction for alternate allele in tumor, very unlikely somatic event, default: %(default)s)")
   optional_tumor_only.add_argument("--exclude_likely_het_germline", action = "store_true", help = "Exclude likely heterozygous germline variants (40-60 pct allelic fraction, AND presence in dbSNP + gnomAD, AND not existing as somatic event in COSMIC/TCGA, default: %(default)s)")
   optional_tumor_only.add_argument("--exclude_dbsnp_nonsomatic", action= "store_true", help= "Exclude variants found in dbSNP (only those that are NOT found in ClinVar(somatic origin)/DoCM/TCGA/COSMIC, defult: %(default)s)")
   optional_tumor_only.add_argument("--exclude_nonexonic", action= "store_true", help= "Exclude non-exonic variants, defult: %(default)s)")

   required.add_argument("--input_vcf", dest = "input_vcf", help = "VCF input file with somatic variants in tumor sample, SNVs/InDels", required = True)
   required.add_argument("--pcgr_dir", dest = "pcgr_dir", help = "PCGR base directory with accompanying data directory, e.g. ~/pcgr-" + str(PCGR_VERSION), required = True)
   required.add_argument("--output_dir", dest = "output_dir", help = "Output directory", required = True)
   required.add_argument("--genome_assembly", dest = "genome_assembly", choices = ["grch37","grch38"], help = "Human genome assembly build: grch37 or grch38", required = True)
   required.add_argument("--sample_id", dest = "sample_id", help = "Tumor sample/cancer genome identifier - prefix for output files", required = True)

   ## Parse required and optional arguments
   args = parser.parse_args()
   arg_dict = vars(args)
   
   logger = getlogger("pcgr-validate-arguments-input")
   print()
   logger.info("PCGR - STEP 0: Validate input data and options")

   ## Required arguments
   ## Check the existence of required arguments
   if arg_dict["pcgr_dir"] is None or not os.path.exists(arg_dict["pcgr_dir"]):
      err_msg = "Required argument '--pcgr_dir' does not exist (" + str(arg_dict["pcgr_dir"]) + "). Type pcgr.py --help to view all options and required arguments"
      pcgr_error_message(err_msg,logger)
   
   if arg_dict["output_dir"] is None or not os.path.exists(arg_dict["output_dir"]):
      err_msg = "Required argument '--output_dir' does not exist (" + str(arg_dict["output_dir"]) + "). Type pcgr.py --help to view all options and required arguments"
      pcgr_error_message(err_msg,logger)
   
   if arg_dict["genome_assembly"] is None:
      err_msg = "Required argument '--genome_assembly' has no/undefined value (" + str(arg_dict["genome_assembly"]) + "). Type pcgr.py --help to view all options and required arguments"
      pcgr_error_message(err_msg,logger)
   
   if arg_dict["input_vcf"] is None:
      err_msg = "Required argument '--input_vcf' does not exist (" + str(arg_dict["input_vcf"]) + "). Type pcgr.py --help to view all options and required arguments"
      pcgr_error_message(err_msg,logger)
   
   if arg_dict["sample_id"] is None:
      err_msg = "Required argument '--sample_id' has no/undefined value (" + str(arg_dict["sample_id"]) + "). Type pcgr.py --help to view all options and required arguments"
      pcgr_error_message(err_msg,logger)
   
   if len(arg_dict["sample_id"]) <= 2 or len(arg_dict["sample_id"]) > 35:
      err_msg = "Sample name identifier (--sample_id) requires a name with more than 2 characters (and less than 35). Current sample identifier: " + str(arg_dict["sample_id"])
      pcgr_error_message(err_msg,logger)

   ## Optional arguments
   ## Check the existence of Docker (if not --no_docker is et)
   global debug
   debug = args.debug
   global DOCKER_IMAGE_VERSION

   if arg_dict["no_docker"]:
      DOCKER_IMAGE_VERSION = None
   else:
      # check that script and Docker image version correspond
      check_docker_command = "docker images -q " + str(DOCKER_IMAGE_VERSION)
      output = subprocess.check_output(str(check_docker_command), stderr=subprocess.STDOUT, shell=True)
      if(len(output) == 0):
         err_msg = "Docker image " + str(DOCKER_IMAGE_VERSION) + " does not exist, pull image from Dockerhub (docker pull " + str(DOCKER_IMAGE_VERSION) + ")"
         pcgr_error_message(err_msg,logger)

   ## check if input is cancer cell line, requires --tumor_only
   if arg_dict["cell_line"] and not arg_dict["tumor_only"]:
      err_msg = "Analysis of cell line (--cell_line) needs option --tumor_only"
      pcgr_error_message(err_msg,logger)

   ## check that tumor primary site/type is set correctly (integer between 0 and 30)
   if arg_dict["tsite"] > max(tsites.keys()) or arg_dict["tsite"] < 0:
      err_msg = "Tumor type code (" + str(arg_dict["tsite"]) + ") needs to be between " + str(min(tsites.keys())) + " and " + str(max(tsites.keys()))
      pcgr_error_message(err_msg,logger)

   ## check that tumor purity and tumor ploidy is set correctly
   if not arg_dict["tumor_purity"] is None:
      if not (arg_dict["tumor_purity"] > 0 and arg_dict["tumor_purity"] <= 1):
         err_msg = "Tumor purity value " + str(arg_dict["tumor_purity"]) + " is not within a valid range [0,1]"
         pcgr_error_message(err_msg,logger)
   
   if not arg_dict["tumor_ploidy"] is None:
      if not arg_dict["tumor_ploidy"] > 0:
         err_msg = "Tumor ploidy value " + str(arg_dict["tumor_ploidy"]) + " is negative"
         pcgr_error_message(err_msg,logger)

   ## check that minimum/maximum depth/allelic fractions are set correctly
   if arg_dict["tumor_dp_min"] < 0:
      err_msg = "Minimum depth tumor ('tumor_dp_min') must be positive, current value is " + str(arg_dict["tumor_dp_min"]) + ")"
      pcgr_error_message(err_msg, logger)

   if arg_dict["tumor_af_min"] < 0 or arg_dict["tumor_af_min"] > 1:
      err_msg = "Minimum AF tumor ('tumor_af_min') must be within the [0,1] range, current value is " + str(arg_dict["tumor_af_min"]) + ")"
      pcgr_error_message(err_msg, logger)
   
   if arg_dict["control_dp_min"] < 0:
      err_msg = "Minimum depth control ('control_dp_min') must be positive, current value is " + str(arg_dict["control_dp_min"]) + ")"
      pcgr_error_message(err_msg, logger)

   if arg_dict["control_af_max"] < 0 or arg_dict["control_af_max"] > 1:
      err_msg = "Maximum AF control ('control_af_max') must be within the [0,1] range, current value is " + str(arg_dict["control_af_max"]) + ")"
      pcgr_error_message(err_msg, logger)
   
   ## Check that coding target size region of sequencing assay is set correctly 
   if arg_dict["target_size_mb"] < 0 or arg_dict["target_size_mb"] > 34:
      err_msg = "Coding target size region in Mb (" + str(arg_dict["target_size_mb"]) + ") is not positive or larger than the likely maximum size of the coding human genome (34 Mb))"
      pcgr_error_message(err_msg,logger)
   if arg_dict["target_size_mb"] < 1:
      warn_msg = "Coding target size region in Mb (" + str(arg_dict["target_size_mb"]) + ") must be greater than 1 Mb for mutational burden estimate to be robust"
      pcgr_warn_message(warn_msg,logger)
   if arg_dict["target_size_mb"] < 34 and arg_dict["assay"] != "TARGETED":
      warn_msg = "Coding target size (" + str(arg_dict["target_size_mb"]) + ") is less than default for WES/WGS (34Mb), assay must be set to 'TARGETED'"
      pcgr_warn_message(warn_msg,logger)

   ## if assay is targeted or mode is Tumor-Only, MSI prediction will not be performed/switched off
   assay_type = "Tumor-Control"
   if arg_dict["estimate_msi_status"] is True and (arg_dict["assay"] == "TARGETED" or arg_dict["tumor_only"] is True):
      if arg_dict["tumor_only"] is True:
         assay_type = "Tumor-Only"
      warn_msg = "MSI status prediction can be applied for WGS/WES tumor-control assays only (query type: " + str(arg_dict["assay"]) + "|" + str(assay_type) + ") - analysis will be omitted"
      pcgr_warn_message(warn_msg, logger)
      arg_dict["estimate_msi_status"] = 0
   
   ## minimum number of mutations required for mutational signature reconstruction cannot be less than 100 (somewhat arbitrary lower threshold, recommended value is 200)
   if arg_dict["min_mutations_signatures"] < 200:
      warn_msg = "Minimum number of mutations required for mutational signature analysis ( n = " + str(arg_dict["min_mutations_signatures"]) + ") is less than the recommended number (n = 200)"
      pcgr_warn_message(warn_msg,logger)
      if arg_dict["min_mutations_signatures"] < 100:
         err_msg = "Minimum number of mutations required for mutational signature analysis ( n = " + str(arg_dict["min_mutations_signatures"]) + ") can not be less than 100"
         pcgr_error_message(err_msg,logger)

   ## if MSI status is to be estimated, mutational burden must be turned on
   if arg_dict["estimate_msi_status"] is True and arg_dict["estimate_tmb"] is False:
      err_msg = "Prediction of MSI status (option '--estimate_msi_status') requires mutational burden analysis ('--estimate_tmb')"
      pcgr_error_message(err_msg,logger)

   if arg_dict["tumor_only"] is True:

      for t in ["exclude_likely_het_germline","exclude_likely_hom_germline"]:
         if arg_dict[t]:
            if arg_dict['tumor_af_tag'] == "_NA_":
               err_msg = "Option '--" + str(t) + "' requires '--tumor_af_tag' option to be set"
               pcgr_error_message(err_msg,logger)

      ## Emit warning if panel-of-normals VCF is not present and exclude_pon is set
      if arg_dict["pon_vcf"] is None and arg_dict["exclude_pon"] is True:
         warn_msg = "Panel-of-normals VCF is NOT provided ('--pon_vcf') - exclusion of calls found in panel-of-normals ('--exclude_pon') will be ignored"
         pcgr_warn_message(warn_msg,logger)
         arg_dict["exclude_pon"] = False

      ## Emit warnings that mutational burden and mutational signatures are less accurate for assays with tumor-only data
      if arg_dict["estimate_tmb"] is True:
         warn_msg = "Estimation of mutational burden in tumor-only mode is suboptimal - results must be interpreted with caution"
         pcgr_warn_message(warn_msg,logger)
      if arg_dict["estimate_signatures"] is True:
         warn_msg = "Estimation of mutational signatures in tumor-only mode is suboptimal - results must be interpreted with caution"
         pcgr_warn_message(warn_msg,logger)

      ## Emit errors when tumor-only filtering thresholds are not properly set
      for pop in ["eur", "afr", "amr", "eas", "sas", "global"]:
         tag = "maf_onekg_" + str(pop)
         if arg_dict[tag]:
            if float(arg_dict[tag]) < 0 or float(arg_dict[tag]) > 1:
               err_msg = "MAF threshold (tumor-only germline filter) for 1000 Genomes Project (pop '" + str(pop).upper() + "') must be within the [0,1} range, current value is " + str(arg_dict[tag])
               pcgr_error_message(err_msg, logger)
      
      for pop in ["nfe", "fin", "amr", "eas", "sas", "asj", "oth", "afr", "global"]:
         tag = "maf_gnomad_" + str(pop)
         if arg_dict[tag]:
            if float(arg_dict[tag]) < 0 or float(arg_dict[tag]) > 1:
               err_msg = "MAF threshold (tumor-only germline filter) for gnomAD (pop '" + str(pop).upper() + "') must be within the [0,1} range, current value is " + str(arg_dict[tag])
               pcgr_error_message(err_msg, logger)

   ## tumor-only is False
   # else:
   #    for t in ["exclude_pon","exclude_likely_het_germline","exclude_likely_hom_germline","exclude_dbsnp_nonsomatic","exclude_nonexonic"]:
   #       if arg_dict[t] is True:
   #          warn_msg = "Option "--" + str(t) + "" requires "--tumor_only" option (not currently set)"
   #          pcgr_warn_message(warn_msg, logger)

   
   ## Emit warning that mutational signature estimation is (likely) not optimal for small/targeted sequencing assays
   if arg_dict["estimate_signatures"] is True and arg_dict["assay"] == "TARGETED":
      warn_msg = "Estimation of mutational signatures (--estimate_signatures) is not optimal for TARGETED sequencing assays - results must be interpreted with caution"
      pcgr_warn_message(warn_msg,logger)

   ## Check that log ratio thresholds for homozygous deletions and amplifications are properly set, and that segment overlap with transcripts are set appropriately
   if arg_dict["logr_homdel"] >= 0:
      err_msg = "Log ratio for homozygous deletions (" + str(arg_dict["logr_homdel"]) + ") should be less than zero"
      pcgr_error_message(err_msg,logger)
   if arg_dict["logr_gain"] <= 0:
      err_msg = "Log ratio for copy number gains/amplifications (" + str(arg_dict["logr_gain"]) + ") should be greater than zero"
      pcgr_error_message(err_msg,logger)
   if arg_dict["cna_overlap_pct"] > 100 or arg_dict["cna_overlap_pct"] <= 0:
      err_msg = "Minimum percent overlap between copy number segment and gene transcript (" + str(arg_dict["cna_overlap_pct"]) + ") should be greater than zero and less than 100"
      pcgr_error_message(err_msg,logger)

   ## VEP options
   if arg_dict["vep_n_forks"] <= 0 or arg_dict["vep_n_forks"] > 4:
      err_msg = "Number of forks that VEP can use during annotation must be above 0 and not more than 4, current value is " + str(arg_dict["vep_n_forks"])
      pcgr_error_message(err_msg,logger)
   
   if arg_dict["vep_buffer_size"] <= 0 or arg_dict["vep_buffer_size"] > 30000:
      err_msg = "Internal VEP buffer size, corresponding to the number of variants that are read in to memory simultaneously, must be above 0 and not more than 30,000, current value is " + str(arg_dict["vep_buffer_size"])
      pcgr_error_message(err_msg,logger)

   ## Check that VEP pick criteria is formatted correctly
   if not arg_dict["vep_pick_order"] is None:
      values = str(arg_dict["vep_pick_order"]).split(",")
      permitted_sources = ["canonical","appris","tsl","biotype","ccds","rank","length","mane"]
      num_permitted_sources = 0
      for v in values:
         if v in permitted_sources:
            num_permitted_sources += 1
               
      if num_permitted_sources != 8:
         err_msg = "Option 'vep_pick_order' = " + str(arg_dict["vep_pick_order"]) + " is formatted incorrectly, should be " + \
            "a comma-separated string of the following values: canonical,appris,tsl,biotype,ccds,rank,length,mane"
         pcgr_error_message(err_msg, logger)


   config_options = {}
   config_options["tumor_purity"] = "NA"
   config_options["tumor_ploidy"] = "NA"
   if not arg_dict["tumor_purity"] is None:
      config_options["tumor_purity"] = float(arg_dict["tumor_purity"])
   if not arg_dict["tumor_ploidy"] is None:
      config_options["tumor_ploidy"] = float(arg_dict["tumor_ploidy"])
   
   config_options["assay"] = arg_dict["assay"]

   config_options["other"] = {}
   config_options["other"]["vcfanno_n_proc"] = int(arg_dict["vcfanno_n_proc"])
   config_options["other"]["vep_buffer_size"] = int(arg_dict["vep_buffer_size"])
   config_options["other"]["vep_pick_order"] = str(arg_dict["vep_pick_order"])
   config_options["other"]["vep_n_forks"] = int(arg_dict["vep_n_forks"])
   config_options["other"]["visual_theme"] = str(arg_dict["report_theme"])
   config_options["other"]["nonfloating_toc"] = int(arg_dict["report_nonfloating_toc"])
   config_options["other"]["vep_no_intergenic"] = int(arg_dict["vep_no_intergenic"])
   config_options["other"]["vep_regulatory"] = int(arg_dict["vep_regulatory"])
   config_options["other"]["vcf2maf"] = int(arg_dict["vcf2maf"])
   config_options["other"]["basic"] = int(arg_dict["basic"])
   config_options["other"]["preserved_info_tags"] = str(arg_dict["preserved_info_tags"])
   config_options["other"]["list_noncoding"] = int(arg_dict["show_noncoding"])
   config_options["other"]["force_overwrite"] = int(arg_dict["force_overwrite"])
   config_options["other"]["no_vcf_validate"] = int(arg_dict["no_vcf_validate"])
   
   config_options["rna"] = {}
   config_options["clinicaltrials"] = {}
   config_options["clinicaltrials"]["run"] = int(arg_dict["include_trials"])

   config_options["msi"] = {}
   config_options["msi"]["run"] = int(arg_dict["estimate_msi_status"])
   config_options["msigs"] = {}
   config_options["msigs"]["run"] = int(arg_dict["estimate_signatures"])
   config_options["msigs"]["mutation_limit"] = int(arg_dict["min_mutations_signatures"])
   config_options["msigs"]["all_reference_signatures"] = int(arg_dict["all_reference_signatures"])
   config_options["msigs"]["include_artefact_signatures"] = int(arg_dict["include_artefact_signatures"])
   
   config_options["tmb"] = {}
   config_options["tmb"]["run"] = int(arg_dict["estimate_tmb"])
   config_options["tmb"]["target_size_mb"] = arg_dict["target_size_mb"]
   config_options["tmb"]["algorithm"] = arg_dict["tmb_algorithm"]
   
   config_options["cna"] = {}
   config_options["cna"]["logR_homdel"] = float(arg_dict["logr_homdel"])
   config_options["cna"]["logR_gain"] = float(arg_dict["logr_gain"])
   config_options["cna"]["cna_overlap_pct"] = float(arg_dict["cna_overlap_pct"])
   
   config_options["allelic_support"] = {}
   config_options["allelic_support"]["tumor_dp_min"] = int(arg_dict["tumor_dp_min"])
   config_options["allelic_support"]["control_dp_min"] = int(arg_dict["control_dp_min"])
   config_options["allelic_support"]["tumor_af_min"] = float(arg_dict["tumor_af_min"])
   config_options["allelic_support"]["control_af_max"] = float(arg_dict["control_af_max"])
   config_options["allelic_support"]["control_dp_tag"] = str(arg_dict["control_dp_tag"])
   config_options["allelic_support"]["control_af_tag"] = str(arg_dict["control_af_tag"])
   config_options["allelic_support"]["tumor_dp_tag"] = str(arg_dict["tumor_dp_tag"])
   config_options["allelic_support"]["tumor_af_tag"] = str(arg_dict["tumor_af_tag"])
   config_options["allelic_support"]["call_conf_tag"] = str(arg_dict["call_conf_tag"])

   config_options["tumor_only"] = {}
   config_options["tumor_only"]["tumor_only"] = int(arg_dict["tumor_only"])
   config_options["tumor_only"]["cell_line"] = int(arg_dict["cell_line"])
   config_options["tumor_only"]["exclude_pon"] = int(arg_dict["exclude_pon"])
   config_options["tumor_only"]["exclude_likely_hom_germline"] = int(arg_dict["exclude_likely_hom_germline"])
   config_options["tumor_only"]["exclude_likely_het_germline"] = int(arg_dict["exclude_likely_het_germline"])
   config_options["tumor_only"]["exclude_dbsnp_nonsomatic"] = int(arg_dict["exclude_dbsnp_nonsomatic"])
   config_options["tumor_only"]["exclude_nonexonic"] = int(arg_dict["exclude_nonexonic"])

   for pop in ["eur", "afr", "amr", "eas", "sas", "global"]:
      tag = "maf_onekg_" + str(pop)
      if arg_dict[tag]:
         config_options["tumor_only"][tag] = float(arg_dict[tag])
      
   for pop in ["nfe", "fin", "amr", "eas", "sas", "asj","oth", "afr", "global"]:
      tag = "maf_gnomad_" + str(pop)
      if arg_dict[tag]:
         config_options["tumor_only"][tag] = float(arg_dict[tag])

   config_options["tumor_type"] = {}
   config_options["tumor_type"]["type"] = str(tsites[arg_dict["tsite"]])

   ## Verify existence of input files
   host_directories = verify_input_files(arg_dict, logger)

   ## Run PCGR workflow ( VEP + vcfanno + summarise + vcf2tsv + HTML report generation)
   run_pcgr(arg_dict, host_directories, config_options)

def pcgr_error_message(message, logger):
   logger.error("")
   logger.error(message)
   logger.error("")
   sys.exit(1)

def pcgr_warn_message(message, logger):
   logger.warning(message)

def verify_input_files(arg_dict, logger):

   """
   Function that 
   1. Checks the input files and directories provided by the user (dictionary arg_dict) and checks for their existence
   2. Checks that the data bundle is of correct date
   """
 
   input_vcf_dir = "NA"
   input_cna_dir = "NA"
   input_rna_fusion_dir = "NA"
   input_cpsr_report_dir = "NA"
   input_rna_expression_dir = "NA"
   input_cna_plot_dir = "NA"
   panel_normal_vcf_dir = "NA"
   db_dir = "NA"
   base_dir = "NA"
   output_dir_full = "NA"
   panel_normal_vcf_basename = "NA"
   input_vcf_basename = "NA"
   input_cna_basename = "NA"
   input_rna_fusion_basename = "NA"
   input_rna_expression_basename = "NA"
   input_cpsr_report_basename = "NA"
   input_cna_plot_basename = "NA"

   arg_dict["rna_fusion_tumor"] = None
   arg_dict["rna_exp_tumor"] = None
   
   ## check that either input vcf or cna segments exist
   if arg_dict["input_vcf"] is None and arg_dict["input_cna"] is None:
      err_msg = "Please specifiy either a VCF input file (--input_vcf) or a copy number segment file (--input_cna)"
      pcgr_error_message(err_msg,logger)
   
   ## check the existence of given output folder
   output_dir_full = os.path.abspath(arg_dict["output_dir"])
   if not os.path.isdir(output_dir_full):
      err_msg = "Output directory (" + str(output_dir_full) + ") does not exist"
      pcgr_error_message(err_msg,logger)
   
   ## check if panel of normal VCF exist
   if not arg_dict["pon_vcf"] is None:
      if not os.path.exists(os.path.abspath(arg_dict["pon_vcf"])):
         err_msg = "Input file (" + str(arg_dict["pon_vcf"]) + ") does not exist"
         pcgr_error_message(err_msg,logger)

      if not (os.path.abspath(arg_dict["pon_vcf"]).endswith(".vcf.gz")):
         err_msg = "Panel of normals VCF file (" + os.path.abspath(arg_dict["pon_vcf"]) + ") does not have the correct file extension (.vcf.gz)"
         pcgr_error_message(err_msg,logger)

      ## check that tabix file exist if bgzipped files is given
      if os.path.abspath(arg_dict["pon_vcf"]).endswith(".vcf.gz"):
         tabix_file = arg_dict["pon_vcf"] + ".tbi"
         if not os.path.exists(os.path.abspath(tabix_file)):
            err_msg = "Tabix file (i.e. '.gz.tbi') is not present for the bgzipped panel of normal VCF file (" + os.path.abspath(arg_dict["pon_vcf"]) + \
               "). Please make sure your the VCF is properly compressed and indexed (bgzip + tabix)"
            pcgr_error_message(err_msg,logger)

      if arg_dict["input_vcf"] is None:
         warn_msg = "Ignoring panel of normal VCF file, --input_vcf missing"
         pcgr_warn_message(warn_msg, logger)
      else:
         panel_normal_vcf_basename = os.path.basename(str(arg_dict["pon_vcf"]))
         panel_normal_vcf_dir = os.path.dirname(os.path.abspath(arg_dict["pon_vcf"]))

   ## check if input vcf exist
   if not arg_dict["input_vcf"] is None:
      if not os.path.exists(os.path.abspath(arg_dict["input_vcf"])):
         err_msg = "Input file (" + str(arg_dict["input_vcf"]) + ") does not exist"
         pcgr_error_message(err_msg,logger)

      if not (os.path.abspath(arg_dict["input_vcf"]).endswith(".vcf") or os.path.abspath(arg_dict["input_vcf"]).endswith(".vcf.gz")):
         err_msg = "VCF input file (" + os.path.abspath(arg_dict["input_vcf"]) + ") does not have the correct file extension (.vcf or .vcf.gz)"
         pcgr_error_message(err_msg,logger)

      ## check that tabix file exist if bgzipped files is given
      if os.path.abspath(arg_dict["input_vcf"]).endswith(".vcf.gz"):
         tabix_file = arg_dict["input_vcf"] + ".tbi"
         if not os.path.exists(os.path.abspath(tabix_file)):
            err_msg = "Tabix file (i.e. '.gz.tbi') is not present for the bgzipped VCF input file (" + os.path.abspath(arg_dict["input_vcf"]) + \
               "). Please make sure your input VCF is properly compressed and indexed (bgzip + tabix)"
            pcgr_error_message(err_msg,logger)

      input_vcf_basename = os.path.basename(str(arg_dict["input_vcf"]))
      input_vcf_dir = os.path.dirname(os.path.abspath(arg_dict["input_vcf"]))

      ## if output vcf exist and overwrite not set
      output_vcf = os.path.join(str(output_dir_full),str(arg_dict["sample_id"])) + ".pcgr_acmg." + str(arg_dict["genome_assembly"])  + ".vcf.gz"
      if os.path.exists(output_vcf) and arg_dict["force_overwrite"] is False:
         err_msg = "Output files (e.g. " + str(output_vcf) + ") already exist - please specify different sample_id or add option --force_overwrite"
         pcgr_error_message(err_msg,logger)
   
   ## check if input cna plot file exist
   # if not arg_dict["input_cna_plot"] is None:
   #    if not os.path.exists(os.path.abspath(arg_dict["input_cna_plot"])):
   #       err_msg = "Input file (" + str(arg_dict["input_cna_plot"]) + ") does not exist"
   #       pcgr_error_message(err_msg,logger)
   #    if not (os.path.abspath(arg_dict["input_cna_plot"]).endswith(".png")):
   #       err_msg = "CNA segment input file (" + os.path.abspath(arg_dict["input_cna_plot"]) + ") does not have the correct file extension (.png)"
   #       pcgr_error_message(err_msg,logger)
   #    if arg_dict["input_cna"] is None:
   #       err_msg = "Input a CNA plot needs to come with a CNA segment file (--input_cna is missing)"
   #       pcgr_error_message(err_msg,logger)
   #    input_cna_plot_basename = os.path.basename(str(arg_dict["input_cna_plot"]))
   #    input_cna_plot_dir = os.path.dirname(os.path.abspath(arg_dict["input_cna_plot"]))

   ## check if input cna segments exist
   if not arg_dict["input_cna"] is None:
      if not os.path.exists(os.path.abspath(arg_dict["input_cna"])):
         err_msg = "Input file (" + str(arg_dict["input_cna"]) + ") does not exist"
         pcgr_error_message(err_msg,logger)
      if not (os.path.abspath(arg_dict["input_cna"]).endswith(".tsv") or os.path.abspath(arg_dict["input_cna"]).endswith(".txt")):
         err_msg = "CNA segment input file (" + os.path.abspath(arg_dict["input_cna"]) + ") does not have the correct file extension (.tsv or .txt)"
         pcgr_error_message(err_msg,logger)
      input_cna_basename = os.path.basename(str(arg_dict["input_cna"]))
      input_cna_dir = os.path.dirname(os.path.abspath(arg_dict["input_cna"]))

      ## if output cna segments exist and overwrite not set
      output_cna_segments = os.path.join(str(output_dir_full), str(arg_dict["sample_id"])) + ".pcgr_acmg." + str(arg_dict["genome_assembly"]) + ".cna_segments.tsv.gz"
      if os.path.exists(output_cna_segments) and arg_dict["force_overwrite"] is False:
         err_msg = "Output files (e.g. " + str(output_cna_segments) + ") already exist - please specify different sample_id or add option --force_overwrite"
         pcgr_error_message(err_msg,logger)
   
   ## check if input rna fusion variants exist
   if not arg_dict["rna_fusion_tumor"] is None:
      if not os.path.exists(os.path.abspath(arg_dict["rna_fusion_tumor"])):
         err_msg = "Input file (" + str(arg_dict["rna_fusion_tumor"]) + ") does not exist"
         pcgr_error_message(err_msg,logger)
      if not (os.path.abspath(arg_dict["rna_fusion_tumor"]).endswith(".tsv") or os.path.abspath(arg_dict["rna_fusion_tumor"]).endswith(".txt")):
         err_msg = "RNA fusion variants file (" + os.path.abspath(arg_dict["rna_fusion_tumor"]) + ") does not have the correct file extension (.tsv or .txt)"
         pcgr_error_message(err_msg,logger)
      input_rna_fusion_basename = os.path.basename(str(arg_dict["rna_fusion_tumor"]))
      input_rna_fusion_dir = os.path.dirname(os.path.abspath(arg_dict["rna_fusion_tumor"]))
   
   ## check if input rna expression exist
   if not arg_dict["rna_exp_tumor"] is None:
      if not os.path.exists(os.path.abspath(arg_dict["rna_exp_tumor"])):
         err_msg = "Input file (" + str(arg_dict["rna_exp_tumor"]) + ") does not exist"
         pcgr_error_message(err_msg,logger)
      if not (os.path.abspath(arg_dict["rna_exp_tumor"]).endswith(".tsv") or os.path.abspath(arg_dict["rna_exp_tumor"]).endswith(".txt")):
         err_msg = "RNA gene expression file (" + os.path.abspath(arg_dict["rna_exp_tumor"]) + ") does not have the correct file extension (.tsv or .txt)"
         pcgr_error_message(err_msg,logger)
      input_rna_expression_basename = os.path.basename(str(arg_dict["rna_exp_tumor"]))
      input_rna_expression_dir = os.path.dirname(os.path.abspath(arg_dict["rna_exp_tumor"]))
   
    ## check if input rna fusion variants exist
   if not arg_dict["cpsr_report"] is None:
      if not os.path.exists(os.path.abspath(arg_dict["cpsr_report"])):
         err_msg = "Input file (" + str(arg_dict["cpsr_report"]) + ") does not exist"
         pcgr_error_message(err_msg,logger)
      if not (os.path.abspath(arg_dict["cpsr_report"]).endswith(".json.gz")):
         err_msg = "CPSR report file (" + os.path.abspath(arg_dict["cpsr_report"]) + ") does not have the correct file extension (.json.gz)"
         pcgr_error_message(err_msg,logger)
      input_cpsr_report_basename = os.path.basename(str(arg_dict["cpsr_report"]))
      input_cpsr_report_dir = os.path.dirname(os.path.abspath(arg_dict["cpsr_report"]))




   ## check the existence of base folder
   base_dir = os.path.abspath(arg_dict["pcgr_dir"])
   if not os.path.isdir(base_dir):
      err_msg = "Base directory (" + str(base_dir) + ") does not exist"
      pcgr_error_message(err_msg,logger)
   
   ## check the existence of data folder within the base folder
   db_dir = os.path.join(os.path.abspath(arg_dict["pcgr_dir"]), "data")
   if not os.path.isdir(db_dir):
      err_msg = "Data directory (" + str(db_dir) + ") does not exist"
      pcgr_error_message(err_msg,logger)
   
   ## check the existence of specified assembly data folder within the base folder
   db_assembly_dir = os.path.join(os.path.abspath(arg_dict["pcgr_dir"]), "data", arg_dict["genome_assembly"])
   if not os.path.isdir(db_assembly_dir):
      err_msg = "Data directory for the specified genome assembly (" + str(db_assembly_dir) + ") does not exist"
      pcgr_error_message(err_msg,logger)
   
   ## check the existence of RELEASE_NOTES (starting from 0.4.0)
   rel_notes_file = os.path.join(os.path.abspath(arg_dict["pcgr_dir"]), "data", arg_dict["genome_assembly"], "RELEASE_NOTES")
   if not os.path.exists(rel_notes_file):
      err_msg = "The PCGR data bundle is outdated - please download the latest data bundle (see github.com/sigven/pcgr for instructions)"
      pcgr_error_message(err_msg,logger)
      
   f_rel_not = open(rel_notes_file,"r")
   compliant_data_bundle = 0
   for line in f_rel_not:
      if DB_VERSION in line:
         compliant_data_bundle = 1
   
   f_rel_not.close()
    
   if compliant_data_bundle == 0:
      err_msg = "The PCGR data bundle is not compliant with the software version - please download the latest software and data bundle (see https://github.com/sigven/pcgr for instructions)"
      pcgr_error_message(err_msg,logger)
   
   host_directories = {}
   host_directories["input_vcf_dir_host"] = input_vcf_dir
   host_directories["input_cna_dir_host"] = input_cna_dir
   host_directories["input_rna_fusion_dir_host"] = input_rna_fusion_dir
   host_directories["input_rna_expression_dir_host"] = input_rna_expression_dir
   host_directories["input_cpsr_report_dir_host"] = input_cpsr_report_dir
   host_directories["input_cna_plot_dir_host"] = input_cna_plot_dir

   host_directories["panel_normal_vcf_dir_host"] = panel_normal_vcf_dir
   host_directories["db_dir_host"] = db_assembly_dir
   host_directories["base_dir_host"] = base_dir
   host_directories["output_dir_host"] = output_dir_full
   host_directories["panel_normal_vcf_basename_host"] = panel_normal_vcf_basename
   host_directories["input_vcf_basename_host"] = input_vcf_basename
   host_directories["input_cna_basename_host"] = input_cna_basename
   host_directories["input_rna_fusion_basename_host"] = input_rna_fusion_basename
   host_directories["input_rna_expression_basename_host"] = input_rna_expression_basename
   host_directories["input_cpsr_report_basename_host"] = input_cpsr_report_basename
   host_directories["input_cna_plot_basename_host"] = input_cna_plot_basename

   
   return host_directories
   

def check_subprocess(logger, command):
   if debug:
      logger.info(command)
   try:
      output = subprocess.check_output(str(command), stderr=subprocess.STDOUT, shell=True)
      if len(output) > 0:
         print (str(output.decode()).rstrip())
   except subprocess.CalledProcessError as e:
      print (e.output.decode())
      exit(0)

def getlogger(logger_name):
   logger = logging.getLogger(logger_name)
   logger.setLevel(logging.DEBUG)

   # create console handler and set level to debug
   ch = logging.StreamHandler(sys.stdout)
   ch.setLevel(logging.DEBUG)

   # add ch to logger
   logger.addHandler(ch)

   # create formatter
   formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s", "20%y-%m-%d %H:%M:%S")

   #add formatter to ch
   ch.setFormatter(formatter)

   return logger

def run_pcgr(arg_dict, host_directories, config_options):
   """
   Main function to run the PCGR workflow
   """

   debug = arg_dict["debug"]
   docker_user_id = arg_dict["docker_user_id"]
   tumor_only = 0
   cell_line = 0
   vcf_validation = 1
   report_nonfloating_toc = 0
   vep_regulatory_annotation = 0
   clinical_trials_set = "OFF"
   msi_prediction_set = "OFF"
   tmb_estimation_set = "OFF"
   msig_estimation_set = "OFF"
   vep_regulatory_annotation = "OFF"
   assay_mode = "Tumor vs. Control"
   if config_options["other"]["nonfloating_toc"]:
      report_nonfloating_toc = 1
   if config_options["other"]["vep_regulatory"] == 1:
      vep_regulatory_annotation = "ON"
   if config_options["clinicaltrials"]["run"]:
      clinical_trials_set = "ON"
   if config_options["msi"]["run"]:
      msi_prediction_set = "ON"
   if config_options["msigs"]["run"]:
      msig_estimation_set = "ON"
   if config_options["tmb"]["run"]:
      tmb_estimation_set = "ON"      
   if config_options["tumor_only"]["tumor_only"]:
      assay_mode = "Tumor-Only"
      tumor_only = 1
      if config_options["tumor_only"]["cell_line"]:
         cell_line = 1
         assay_mode = "Tumor-Only (cell line)"
   if config_options["other"]["no_vcf_validate"]:
      vcf_validation = 0
   ## set basic Docker run commands
   output_vcf = "None"
   output_pass_vcf = "None"
   output_pass_tsv = "None"
   output_maf = "None"
   uid = ""
   
   global GENCODE_VERSION, VEP_ASSEMBLY, NCBI_BUILD_MAF
   if arg_dict["genome_assembly"] == "grch37":
      NCBI_BUILD_MAF = "GRCh37"
      GENCODE_VERSION = "release 19"
      VEP_ASSEMBLY = "GRCh37"
   logger = getlogger("pcgr-get-OS")

   if docker_user_id:
      uid = docker_user_id
   elif platform.system() == "Linux" or platform.system() == "Darwin" or sys.platform == "darwin" or sys.platform == "linux2" or sys.platform == "linux":
      uid = os.getuid()
   else:
      if platform.system() == "Windows" or sys.platform == "win32" or sys.platform == "cygwin":
         uid = getpass.getuser()

   if uid == "":
      logger.warning("Was not able to get user id/username for logged-in user on the underlying platform (platform.system(): " + \
         str(platform.system()) + ", sys.platform: " + str(sys.platform) + "), now running PCGR as root")
      uid = "root"

   vepdb_dir_host = os.path.join(str(host_directories["db_dir_host"]),".vep")
   input_vcf_docker = "None"
   input_cna_docker = "None"
   input_rna_fusion_docker = "None"
   input_rna_expression_docker = "None"
   input_cpsr_report_docker = "None"
   panel_normal_docker = "None"
   docker_cmd_run1 = ""
   docker_cmd_run2 = ""
   docker_cmd_run_end = ""
   ## panel-of-normals annotation
   pon_annotation = 0

   ## Map input files and directories to files/directories within the Docker container
   if DOCKER_IMAGE_VERSION:
      if host_directories["input_vcf_basename_host"] != "NA":
         input_vcf_docker = "/workdir/input_vcf/" + str(host_directories["input_vcf_basename_host"])
      if host_directories["input_cna_basename_host"] != "NA":
         input_cna_docker = "/workdir/input_cna/" + str(host_directories["input_cna_basename_host"])
      if host_directories["input_rna_expression_basename_host"] != "NA":
         input_rna_expression_docker = "/workdir/input_rna_expression/" + str(host_directories["input_rna_expression_basename_host"])
      if host_directories["input_rna_fusion_basename_host"] != "NA":
         input_rna_fusion_docker = "/workdir/input_rna_fusion/" + str(host_directories["input_rna_fusion_basename_host"])
      if host_directories["input_cpsr_report_basename_host"] != "NA":
         input_cpsr_report_docker = "/workdir/input_cpsr/" + str(host_directories["input_cpsr_report_basename_host"])
      if host_directories["panel_normal_vcf_basename_host"] != "NA":
         panel_normal_docker = "/workdir/panel_normal_vcf/" + str(host_directories["panel_normal_vcf_basename_host"])

      vep_volume_mapping = str(vepdb_dir_host) + ":/usr/local/share/vep/data"
      databundle_volume_mapping = str(host_directories["base_dir_host"]) + ":/data"
      input_cna_volume_mapping = str(host_directories["input_cna_dir_host"]) + ":/workdir/input_cna"
      input_vcf_volume_mapping = str(host_directories["input_vcf_dir_host"]) + ":/workdir/input_vcf"
      input_rna_fusion_volume_mapping = str(host_directories["input_rna_fusion_dir_host"]) + ":/workdir/input_rna_fusion"
      input_rna_expression_volume_mapping = str(host_directories["input_rna_expression_dir_host"]) + ":/workdir/input_rna_expression"
      input_cpsr_report_volume_mapping = str(host_directories["input_cpsr_report_dir_host"]) + ":/workdir/input_cpsr"
      output_volume_mapping = str(host_directories["output_dir_host"]) + ":/workdir/output"
      panel_normal_vcf_volume_mapping = str(host_directories["panel_normal_vcf_dir_host"]) + ":/workdir/panel_normal_vcf"

      docker_cmd_run1 = "NA"

      ## VCF file only
      docker_run_basic = "docker run --rm -t -u " + str(uid)
      docker_cmd_run1 = str(docker_run_basic) + " -v=" +  str(databundle_volume_mapping) + " -v=" + str(vep_volume_mapping) + " -v=" + str(output_volume_mapping)
      if host_directories["input_vcf_dir_host"] != "NA" and host_directories["input_cna_dir_host"] == "NA":
         docker_cmd_run1 = docker_cmd_run1  + " -v=" + str(input_vcf_volume_mapping)

      ## CNA file and VCF file provided
      if host_directories["input_vcf_dir_host"] != "NA" and host_directories["input_cna_dir_host"] != "NA":
         docker_cmd_run1 = docker_cmd_run1  + " -v=" + str(input_vcf_volume_mapping) + " -v=" + str(input_cna_volume_mapping)

      ## Panel of normal VCFs provided
      if host_directories["panel_normal_vcf_dir_host"] != "NA":
         docker_cmd_run1 = docker_cmd_run1 + " -v=" + str(panel_normal_vcf_volume_mapping)

      ## RNA fusion variants provided
      if host_directories["input_rna_fusion_dir_host"] != "NA":
         docker_cmd_run1 = docker_cmd_run1 + " -v=" + str(input_rna_fusion_volume_mapping)

      ## RNA expression estimates provided
      if host_directories["input_rna_expression_dir_host"] != "NA": 
         docker_cmd_run1 = docker_cmd_run1 + " -v=" + str(input_rna_expression_volume_mapping)
      
      ## CPSR report provided
      if host_directories["input_cpsr_report_dir_host"] != "NA": 
         docker_cmd_run1 = docker_cmd_run1 + " -v=" + str(input_cpsr_report_volume_mapping)


      docker_cmd_run1 = docker_cmd_run1 + " -w=/workdir/output " + str(DOCKER_IMAGE_VERSION) + " sh -c \""
      
      docker_cmd_run2 = str(docker_run_basic) + " -v=" + str(databundle_volume_mapping) + " -v=" + str(output_volume_mapping)
      if host_directories["panel_normal_vcf_dir_host"] != "NA":
         docker_cmd_run2 = docker_cmd_run2 + " -v=" + str(panel_normal_vcf_volume_mapping)
      docker_cmd_run2 = docker_cmd_run2 + " -w=/workdir/output " + str(DOCKER_IMAGE_VERSION) + " sh -c \""
      docker_cmd_run_end = "\""

      data_dir = "/data"
      output_dir = "/workdir/output"
      vep_dir = "/usr/local/share/vep/data"
      r_scripts_dir = "/"

   ## If running non-Dockerized - specifiy paths for input files and directories
   else:
      if host_directories["input_vcf_basename_host"] != "NA":
         input_vcf_docker = os.path.join(host_directories["input_vcf_dir_host"], host_directories["input_vcf_basename_host"])
      if host_directories["input_cna_basename_host"] != "NA":
         input_cna_docker = os.path.join(host_directories["input_cna_dir_host"], host_directories["input_cna_basename_host"])
      if host_directories["input_rna_fusion_basename_host"] != "NA":
         input_rna_fusion_docker = os.path.join(host_directories["input_rna_fusion_dir_host"], host_directories["input_rna_fusion_basename_host"])
      if host_directories["input_rna_expression_basename_host"] != "NA":
         input_rna_expression_docker = os.path.join(host_directories["input_rna_expression_dir_host"], host_directories["input_rna_expression_basename_host"])
      if host_directories["input_cpsr_report_basename_host"] != "NA":
         input_cpsr_report_docker = os.path.join(host_directories["input_cpsr_report_dir_host"], host_directories["input_cpsr_report_basename_host"])
      if host_directories["panel_normal_vcf_basename_host"] != "NA":
         panel_normal_docker = os.path.join(host_directories["panel_normal_vcf_dir_host"], host_directories["panel_normal_vcf_basename_host"])

      data_dir = host_directories["base_dir_host"]
      output_dir = host_directories["output_dir_host"]
      vep_dir = vepdb_dir_host
      r_scripts_dir = ""

   check_subprocess(logger, docker_cmd_run1.replace("-u " + str(uid), "") + "mkdir -p " + output_dir + docker_cmd_run_end)

   ## PCGR|validate_input - verify that VCF and CNA segment file is of appropriate format
   logger = getlogger("pcgr-validate-arguments-input")
   vcf_validate_command = docker_cmd_run1 + "pcgr_validate_input.py " + data_dir + " " + \
            str(input_vcf_docker) + " " + \
            str(input_cna_docker) + " " + \
            str(input_rna_fusion_docker) + " " + \
            str(input_rna_expression_docker) + " " + \
            str(panel_normal_docker) + " " + \
            str(vcf_validation) + " " + \
            str(tumor_only) + " " + \
            str(arg_dict["genome_assembly"]) + " " + \
            str(config_options["other"]["preserved_info_tags"]) + " " + \
            str(config_options["allelic_support"]["tumor_dp_tag"]) + " " + str(config_options["allelic_support"]["tumor_af_tag"]) + " " + \
            str(config_options["allelic_support"]["control_dp_tag"]) + " " + str(config_options["allelic_support"]["control_af_tag"]) + " " + \
            str(config_options["allelic_support"]["call_conf_tag"]) + " " + \
            str(config_options["tumor_only"]["exclude_likely_hom_germline"]) + " " + \
            str(config_options["tumor_only"]["exclude_likely_het_germline"])
   if not DOCKER_IMAGE_VERSION:
      vcf_validate_command += " --output_dir " + output_dir + docker_cmd_run_end
   else:
      vcf_validate_command += docker_cmd_run_end
   check_subprocess(logger, vcf_validate_command)
   logger.info("Finished")

   ## PCGR|start - Log key information about sample, options and sequencing assay/design
   logger = getlogger("pcgr-start")
   print()
   logger.info("--- Personal Cancer Genome Reporter workflow ----")
   logger.info("Sample name: " + str(arg_dict["sample_id"]))
   if config_options["tumor_type"]["type"] == "Cancer_NOS":
      logger.info("Tumor type: Cancer_NOS (Any tumortype)")
   else:
      logger.info("Tumor type: " + str(config_options["tumor_type"]["type"]))
   logger.info("Sequencing assay - type: " + str(config_options["assay"]))
   logger.info("Sequencing assay - mode: " + str(assay_mode))
   logger.info("Sequencing assay - coding target size: " + str(config_options["tmb"]["target_size_mb"]) + "Mb")
   logger.info("Genome assembly: " + str(arg_dict["genome_assembly"]))
   logger.info("Mutational signature estimation: " + str(msig_estimation_set))
   logger.info("MSI classification: " + str(msi_prediction_set))
   logger.info("Mutational burden estimation: " + str(tmb_estimation_set))
   logger.info("Include molecularly targeted clinical trials (beta): " + str(clinical_trials_set))

   if not input_vcf_docker == "None":

      ## Define temporary output file names
      output_vcf = os.path.join(output_dir, str(arg_dict["sample_id"]) + ".pcgr_acmg."  + str(arg_dict["genome_assembly"]) + ".vcf.gz")
      output_pass_vcf = os.path.join(output_dir, str(arg_dict["sample_id"]) + ".pcgr_acmg."  + str(arg_dict["genome_assembly"]) + ".pass.vcf.gz")
      output_pass_tsv = os.path.join(output_dir, str(arg_dict["sample_id"]) + ".pcgr_acmg."  + str(arg_dict["genome_assembly"]) + ".pass.tsv")
      output_maf = os.path.join(output_dir, str(arg_dict["sample_id"]) + ".pcgr_acmg." + str(arg_dict["genome_assembly"]) + ".tmp.maf")
      output_vcf2maf_log = os.path.join(output_dir, str(arg_dict["sample_id"]) + ".pcgr_acmg." + str(arg_dict["genome_assembly"]) + ".maf.log")
      input_vcf_pcgr_ready = os.path.join(output_dir, re.sub(r"(\.vcf$|\.vcf\.gz$)", ".pcgr_ready.vcf.gz", host_directories["input_vcf_basename_host"]))
      input_vcf_pcgr_ready_uncompressed = os.path.join(output_dir, re.sub(r"(\.vcf$|\.vcf\.gz$)", ".pcgr_ready.vcf", host_directories["input_vcf_basename_host"]))
      vep_vcf = re.sub(r"(\.vcf$|\.vcf\.gz$)", ".vep.vcf", input_vcf_pcgr_ready)
      vep_vcfanno_vcf = re.sub(r"(\.vcf$|\.vcf\.gz$)", ".vep.vcfanno.vcf", input_vcf_pcgr_ready)
      vep_vcfanno_annotated_vcf = re.sub(r"\.vcfanno", ".vcfanno.annotated", vep_vcfanno_vcf) + ".gz"
      vep_vcfanno_annotated_pass_vcf = re.sub(r"\.vcfanno", ".vcfanno.annotated.pass", vep_vcfanno_vcf) + ".gz"

      ## Path for human genome assembly file (FASTA)
      fasta_assembly = os.path.join(vep_dir, "homo_sapiens", str(VEP_VERSION) + "_" + str(VEP_ASSEMBLY), "Homo_sapiens." + str(VEP_ASSEMBLY) + ".dna.primary_assembly.fa.gz")
      
      ## List all VEP flags used when calling VEP
      vep_flags = "--hgvs --af --af_1kg --af_gnomad --variant_class --domains --symbol --protein --ccds --mane " + \
         "--uniprot --appris --biotype --tsl --canonical --gencode_basic --cache --numbers --total_length --allele_number " + \
         "--no_stats --no_escape --xref_refseq --vcf --check_ref --dont_skip --flag_pick_allele --plugin NearestExonJB,max_range=50000"
      vep_options = "--pick_order " + str(config_options["other"]["vep_pick_order"]) + " --force_overwrite --buffer_size " + \
         str(config_options["other"]["vep_buffer_size"]) + " --species homo_sapiens --assembly " + \
         str(VEP_ASSEMBLY) + " --offline --fork " + str(config_options["other"]["vep_n_forks"]) + " " + str(vep_flags) + " --dir " + str(vep_dir)
      vep_options += " --cache_version " + str(VEP_VERSION)
      if config_options["other"]["vep_no_intergenic"] == 1:
         vep_options = vep_options + " --no_intergenic"
      if config_options["other"]["vep_regulatory"] == 1:
         vep_options += " --regulatory"
      if not debug:
         vep_options += " --quiet"
      if debug:
         vep_options += " --verbose"

      ## Compose full VEP command
      vep_main_command = docker_cmd_run1 + "vep --input_file " + str(input_vcf_pcgr_ready) + " --output_file " + str(vep_vcf) + " " + str(vep_options) + " --fasta " + str(fasta_assembly) + docker_cmd_run_end
      vep_bgzip_command = docker_cmd_run1 + "bgzip -f -c " + str(vep_vcf) + " > " + str(vep_vcf) + ".gz" + docker_cmd_run_end
      vep_tabix_command = docker_cmd_run1 + "tabix -f -p vcf " + str(vep_vcf) + ".gz" + docker_cmd_run_end

      ## PCGR|VEP - run consequence annotation with Variant Effect Predictor
      print()
      logger = getlogger("pcgr-vep")
      logger.info("PCGR - STEP 1: Basic variant annotation with Variant Effect Predictor (" + str(VEP_VERSION) + ", GENCODE " + str(GENCODE_VERSION) + ", " + str(arg_dict["genome_assembly"]) + ")")
      logger.info("VEP configuration - one primary consequence block pr. alternative allele (--flag_pick_allele)")
      logger.info("VEP configuration - transcript pick order: " + str(config_options["other"]["vep_pick_order"]))
      logger.info("VEP configuration - transcript pick order: See more at https://www.ensembl.org/info/docs/tools/vep/script/vep_other.html#pick_options")
      logger.info("VEP configuration - skip intergenic: " + str(config_options["other"]["vep_no_intergenic"]))
      logger.info("VEP configuration - regulatory annotation: " + str(vep_regulatory_annotation))
      logger.info("VEP configuration - buffer_size/number of forks: " + str(arg_dict["vep_buffer_size"]) + "/" + str(arg_dict["vep_n_forks"]))

      check_subprocess(logger, vep_main_command)
      check_subprocess(logger, vep_bgzip_command)
      check_subprocess(logger, vep_tabix_command)

      ## PCGR|vcf2maf - if option set, convert VCF to MAF with https://github.com/mskcc/vcf2maf
      if config_options["other"]["vcf2maf"] == 1:
         logger.info("Converting VEP-annotated VCF to MAF with https://github.com/mskcc/vcf2maf")
         vcf2maf_command = str(docker_cmd_run1) + "vcf2maf.pl --inhibit-vep " + \
            "--input-vcf " + str(input_vcf_pcgr_ready_uncompressed) + " --tumor-id " + \
            str(arg_dict["sample_id"]) + " --output-maf " + str(output_maf) + " --ref-fasta " + str(fasta_assembly) + " --filter-vcf 0 --ncbi-build " + \
               str(NCBI_BUILD_MAF) + " > " + str(output_vcf2maf_log) + " 2>&1" + docker_cmd_run_end
         clean_vcf2maf_command = str(docker_cmd_run1) + "rm -f " + str(output_vcf2maf_log) + " " + re.sub(r"(\.vcf$)", ".vep.vcf", input_vcf_pcgr_ready_uncompressed) + " " + docker_cmd_run_end
         check_subprocess(logger, vcf2maf_command)
         check_subprocess(logger, clean_vcf2maf_command)
      logger.info("Finished")

      ## PCGR|vcfanno - annotate VCF against a number of variant annotation resources
      print()
      logger = getlogger("pcgr-vcfanno")
      pcgr_vcfanno_command = str(docker_cmd_run2) + "pcgr_vcfanno.py --num_processes " + str(config_options["other"]["vcfanno_n_proc"]) + \
         " --chasmplus --dbnsfp --docm --clinvar --icgc --civic --cgi --tcga_pcdm --winmsk --simplerepeats --tcga --uniprot --cancer_hotspots --pcgr_onco_xref " + \
            str(vep_vcf) + ".gz " + str(vep_vcfanno_vcf) + " " + os.path.join(data_dir, "data", str(arg_dict["genome_assembly"]))
      if debug:
         pcgr_vcfanno_command += " --debug"
      if panel_normal_docker != "None":
         pon_annotation = 1
         pcgr_vcfanno_command = pcgr_vcfanno_command + " --panel_normal_vcf " + str(panel_normal_docker)
         logger.info("PCGR - STEP 2: Annotation for precision oncology with pcgr-vcfanno")
         logger.info("Annotation sources: Panel-of-Normals, ClinVar, dbNSFP, UniProtKB, cancerhotspots.org, CiVIC, CGI, DoCM, CHASMplus driver mutations, TCGA, ICGC-PCAWG")
      else:
         logger.info("PCGR - STEP 2: Annotation for precision oncology with pcgr-vcfanno")
         logger.info("Annotation sources: ClinVar, dbNSFP, UniProtKB, cancerhotspots.org, CiVIC, CGI, DoCM, CHASMplus driver mutations, TCGA, ICGC-PCAWG")
      pcgr_vcfanno_command = pcgr_vcfanno_command + docker_cmd_run_end
      check_subprocess(logger, pcgr_vcfanno_command)
      logger.info("Finished")

      ## PCGR|pcgr_summarise - expand annotations in VCF file
      print()
      logger = getlogger("pcgr-summarise")
      pcgr_summarise_command = str(docker_cmd_run2) + "pcgr_summarise.py " + str(vep_vcfanno_vcf) + ".gz " + str(pon_annotation) + " " + \
         str(config_options["other"]["vep_regulatory"]) + " " + str(os.path.join(data_dir, "data", str(arg_dict["genome_assembly"]))) + docker_cmd_run_end
      if debug:
         pcgr_summarise_command  += " --debug"
      logger.info("PCGR - STEP 3: Cancer gene annotations with pcgr-summarise")
      check_subprocess(logger, pcgr_summarise_command)

      ## PCGR|clean - move output files and clean up temporary files
      create_output_vcf_command1 = str(docker_cmd_run2) + "mv " + str(vep_vcfanno_annotated_vcf) + " " + str(output_vcf) + docker_cmd_run_end
      create_output_vcf_command2 = str(docker_cmd_run2) + "mv " + str(vep_vcfanno_annotated_vcf) + ".tbi " + str(output_vcf) + ".tbi" + docker_cmd_run_end
      create_output_vcf_command3 = str(docker_cmd_run2) + "mv " + str(vep_vcfanno_annotated_pass_vcf) + " " + str(output_pass_vcf) + docker_cmd_run_end
      create_output_vcf_command4 = str(docker_cmd_run2) + "mv " + str(vep_vcfanno_annotated_pass_vcf) + ".tbi " + str(output_pass_vcf) + ".tbi" + docker_cmd_run_end
      clean_command = str(docker_cmd_run2) + "rm -f " + str(vep_vcf) + "* " + str(vep_vcfanno_annotated_vcf) + " " + \
         str(vep_vcfanno_annotated_pass_vcf) + "* " + str(vep_vcfanno_vcf) + "* " +  str(input_vcf_pcgr_ready_uncompressed) + "* "  + docker_cmd_run_end
      check_subprocess(logger, create_output_vcf_command1)
      check_subprocess(logger, create_output_vcf_command2)
      check_subprocess(logger, create_output_vcf_command3)
      check_subprocess(logger, create_output_vcf_command4)

      ## PCGR|vcf2tsv - convert VCF to TSV with https://github.com/sigven/vcf2tsv
      pcgr_vcf2tsv_command = str(docker_cmd_run2) + "vcf2tsv.py " + str(output_pass_vcf) + " --compress " + str(output_pass_tsv) + docker_cmd_run_end
      logger.info("Converting VCF to TSV with https://github.com/sigven/vcf2tsv")
      check_subprocess(logger, pcgr_vcf2tsv_command)
      if not debug:
         check_subprocess(logger, clean_command)
      logger.info("Finished")

   print()

   ## Generation of HTML reports for VEP/vcfanno-annotated VCF and copy number segment file
   if not arg_dict["basic"]:
      ttype = config_options["tumor_type"]["type"].replace(" ","_").replace("/","@")
      logger = getlogger("pcgr-writer")
      logger.info("PCGR - STEP 4: Generation of output files - variant interpretation report for precision oncology")
      pcgr_report_command = docker_cmd_run1 + os.path.join(r_scripts_dir, "pcgr.R") + " " + \
         output_dir + " " + \
            str(output_pass_tsv) + ".gz" + " " + \
                           input_cna_docker + " " + \
                           input_rna_fusion_docker + " " + \
                           input_rna_expression_docker + " " + \
                           input_cpsr_report_docker + " " + \
                           str(arg_dict["sample_id"]) + " " + \
                           str(PCGR_VERSION) + " " +  \
                           str(arg_dict["genome_assembly"]) + " " + \
                           str(data_dir) + " " + \
                           str(config_options["tumor_purity"]) + " " + \
                           str(config_options["tumor_ploidy"]) + " " + \
                           str(ttype) + " " + \
                           str(config_options["tmb"]["target_size_mb"]) + " " + \
                           str(config_options["assay"]) + " " + \
                           str(tumor_only) +  " " + \
                           str(cell_line) + " " + \
                           str(config_options["tumor_only"]["maf_onekg_afr"]) + " " + \
                           str(config_options["tumor_only"]["maf_onekg_amr"]) + " " + \
                           str(config_options["tumor_only"]["maf_onekg_eas"]) + " " + \
                           str(config_options["tumor_only"]["maf_onekg_eur"]) + " " + \
                           str(config_options["tumor_only"]["maf_onekg_sas"]) + " " + \
                           str(config_options["tumor_only"]["maf_onekg_global"]) + " " + \
                           str(config_options["tumor_only"]["maf_gnomad_afr"]) + " " + \
                           str(config_options["tumor_only"]["maf_gnomad_amr"]) + " " + \
                           str(config_options["tumor_only"]["maf_gnomad_asj"]) + " " + \
                           str(config_options["tumor_only"]["maf_gnomad_eas"]) + " " + \
                           str(config_options["tumor_only"]["maf_gnomad_fin"]) + " " + \
                           str(config_options["tumor_only"]["maf_gnomad_nfe"]) + " " + \
                           str(config_options["tumor_only"]["maf_gnomad_oth"]) + " " + \
                           str(config_options["tumor_only"]["maf_gnomad_sas"]) + " " + \
                           str(config_options["tumor_only"]["maf_gnomad_global"]) + " " + \
                           str(config_options["tumor_only"]["exclude_pon"]) + " " + \
                           str(config_options["tumor_only"]["exclude_likely_hom_germline"]) + " " + \
                           str(config_options["tumor_only"]["exclude_likely_het_germline"]) + " " + \
                           str(config_options["tumor_only"]["exclude_dbsnp_nonsomatic"]) + " " + \
                           str(config_options["tumor_only"]["exclude_nonexonic"]) + " " + \
                           str(config_options["tmb"]["run"]) + " " + \
                           str(config_options["tmb"]["algorithm"]) + " " + \
                           str(config_options["msi"]["run"]) + " " + \
                           str(config_options["msigs"]["run"]) + " " + \
                           str(config_options["msigs"]["mutation_limit"]) + " " + \
                           str(config_options["msigs"]["all_reference_signatures"]) + " " + \
                           str(config_options["msigs"]["include_artefact_signatures"]) + " " + \
                           str(config_options["cna"]["logR_homdel"]) + " " + \
                           str(config_options["cna"]["logR_gain"]) + " " + \
                           str(config_options["cna"]["cna_overlap_pct"]) + " " + \
                           str(config_options["allelic_support"]["tumor_af_min"]) + " " + \
                           str(config_options["allelic_support"]["tumor_dp_min"]) + " " + \
                           str(config_options["allelic_support"]["control_af_max"]) + " " + \
                           str(config_options["allelic_support"]["control_dp_min"]) + " " + \
                           str(config_options["allelic_support"]["tumor_af_tag"]) + " " + \
                           str(config_options["allelic_support"]["tumor_dp_tag"]) + " " + \
                           str(config_options["allelic_support"]["control_af_tag"]) + " " + \
                           str(config_options["allelic_support"]["control_dp_tag"]) + " " + \
                           str(config_options["allelic_support"]["call_conf_tag"]) + " " + \
                           str(config_options["clinicaltrials"]["run"]) + " " + \
                           str(config_options["other"]["vep_n_forks"]) + " " + \
                           str(config_options["other"]["vep_buffer_size"]) + " " + \
                           str(config_options["other"]["vep_no_intergenic"]) + " " + \
                           str(config_options["other"]["vep_pick_order"]) + " " + \
                           str(config_options["other"]["vep_regulatory"]) + " " + \
                           str(config_options["other"]["vcf2maf"]) + " " + \
                           str(config_options["other"]["list_noncoding"]) + " " + \
                           str(config_options["other"]["preserved_info_tags"]) + " " + \
                           str(config_options["other"]["visual_theme"]) + " " + \
                           str(report_nonfloating_toc) + " " + \
                           str(config_options["other"]["no_vcf_validate"]) + docker_cmd_run_end
       
      check_subprocess(logger, pcgr_report_command)
      logger.info("Finished")

   print()



if __name__=="__main__": __main__()

