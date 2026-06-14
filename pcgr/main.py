#!/usr/bin/env python

from pcgr import pcgr_vars, arg_checker, utils, cna
from pcgr.utils import getlogger, check_subprocess, remove_file, random_id_generator, pd_to_csv, error_message
from pcgr.config import populate_config_data, create_config, verify_oncotree_code
from pcgr.maf import update_maf, generate_oncokb_maf_input, add_var_id_to_vcf, append_oncokb_snv_annotations
from pcgr.cna import append_oncokb_fusion_annotations, append_oncokb_cna_annotations
from pcgr.oncogenicity import refine_oncogenicity_with_oncokb
from pcgr.vep import get_vep_command
from pcgr.expression import parse_expression, integrate_variant_expression, correlate_sample_expression
from pcgr.expression import find_expression_outliers, aggregate_tpm_per_cons
from pcgr.variant import clean_annotations, set_allelic_support, append_annotations, calculate_tmb, reduce_variants_for_report

import re
import argparse
import pandas as pd
import yaml
import os
from glob import glob
from argparse import RawTextHelpFormatter


def cli():

    program_description = (f"Personal Cancer Genome Reporter (PCGR) workflow for clinical translation of "
                           f"tumor omics data (SNVs/InDels, CNA, RNA expression, RNA fusions) - version: {pcgr_vars.PCGR_VERSION}")
    program_options = "\n\t[--input_vcf <INPUT_VCF>]\n\t[--input_cna <INPUT_CNA>]\n\t[--input_rna_fusion <INPUT_RNA_FUSION>]\n\t" + \
        "[--input_rna_expression <INPUT_RNA_EXPRESSION>]\n\t[--vep_dir <VEP_DIR>]\n\t--refdata_dir <REFDATA_DIR>\n\t" + \
        "--output_dir <OUTPUT_DIR>\n\t--genome_assembly <GENOME_ASSEMBLY>\n\t--sample_id <SAMPLE_ID>"

    parser = argparse.ArgumentParser(description=program_description,
                                     formatter_class=RawTextHelpFormatter,
                                     usage=f"\n\t%(prog)s -h [options] {program_options} \n\n")
    parser._action_groups.pop()
    required = parser.add_argument_group("Required arguments")
    optional_input = parser.add_argument_group(
        "Molecular input data",
        "At least one input type must be provided")
    optional_sample = parser.add_argument_group("Tumor sample options")
    optional_snv_indel = parser.add_argument_group(
        "Somatic SNV/InDel options",
        "Only applicable when --input_vcf is provided")
    optional_allelic_support = parser.add_argument_group(
        "Allelic support options (SNV/InDel)",
        "Only applicable when --input_vcf is provided")
    optional_tumor_only = parser.add_argument_group(
        "Tumor-only filtering options (SNV/InDel)",
        "Only applicable when --input_vcf is provided with --tumor_only")
    optional_tmb_msi = parser.add_argument_group(
        "TMB and MSI options (SNV/InDel)",
        "Only applicable when --input_vcf is provided")
    optional_signatures = parser.add_argument_group(
        "Mutational signature options (SNV/InDel)",
        "Only applicable when --input_vcf is provided")
    optional_cna = parser.add_argument_group(
        "Somatic CNA analysis options",
        "Only applicable when --input_cna is provided")
    optional_rna = parser.add_argument_group("RNA expression and fusion options")
    optional_germline = parser.add_argument_group("Germline variant options")
    optional_biomarker = parser.add_argument_group("Biomarker and tiering options")
    optional_other = parser.add_argument_group("Other options")

    required.add_argument("--refdata_dir", dest="refdata_dir", 
                          help = "Directory where PCGR reference data bundle was downloaded and unpacked", required=True)
    required.add_argument("--output_dir", dest="output_dir", 
                          help = "Output directory", required=True)
    required.add_argument("--genome_assembly", dest="genome_assembly", choices=["grch37", "grch38"], 
                          help = "Human genome assembly build: grch37 or grch38", required=True)
    required.add_argument("--sample_id", dest="sample_id", 
                          help = "Tumor sample/cancer genome identifier - prefix for output files", required=True)

    optional_input.add_argument("--input_vcf", dest="input_vcf", default=None, 
                                help = "VCF input file with somatic variants in tumor sample, SNVs/InDels")
    optional_input.add_argument("--input_cna", dest = "input_cna", 
                                help = "Somatic copy number alteration segments (tab-separated values)")
    optional_input.add_argument("--input_rna_fusion", dest = "input_rna_fusion", 
                                help = "File with RNA fusion transcripts detected in tumor (tab-separated values)")
    optional_input.add_argument("--input_rna_expression", dest = "input_rna_exp", 
                                help = "File with bulk RNA expression counts (TPM) of transcripts in tumor (tab-separated values)")

    optional_sample.add_argument("--sex", dest="sex", choices=["FEMALE", "MALE", "UNKNOWN"], default="UNKNOWN", 
                                 help = "Sex of cancer case/sample (default: %(default)s)")
    optional_sample.add_argument("--tumor_site", dest="tsite", type=int, default=0, 
                                 help = "Optional integer code to specify primary tumor type/site of query sample,\nchoose any of the following identifiers:\n" + str(pcgr_vars.tumor_sites) + "\n(default: %(default)s - any tumor type)")
    optional_sample.add_argument("--tumor_purity", type=float, dest="tumor_purity", 
                                 help = "Estimated tumor purity (between 0 and 1) (default: %(default)s)")
    optional_sample.add_argument("--tumor_ploidy", type=float, dest="tumor_ploidy", 
                                 help = "Estimated tumor ploidy (default: %(default)s)")

    optional_snv_indel.add_argument("--vep_dir", dest="vep_dir", 
                                    help = "Directory of VEP cache, e.g.  $HOME/.vep (required when --input_vcf is provided)", required=False, default=None)
    optional_snv_indel.add_argument("--assay", dest="assay", default="WES", choices=[ "WGS", "WES","TARGETED"], 
                                    help = "Type of DNA sequencing assay performed for input data (VCF), default: %(default)s")
    optional_snv_indel.add_argument("--effective_target_size_mb", type=float, default=34, dest="effective_target_size_mb", 
                                    help = "Effective target size in Mb (potentially limited by read depth) of sequencing assay (for TMB analysis) (default: %(default)s (WES/WGS))")
    optional_snv_indel.add_argument("--tumor_only", action="store_true", 
                                    help = "Input VCF comes from tumor-only sequencing, calls will be filtered for variants of germline origin, (default: %(default)s)")
    optional_snv_indel.add_argument("--vcf2maf", action="store_true", 
                                    help = "Generate a MAF file for input VCF using https://github.com/mskcc/vcf2maf (default: %(default)s)")
    optional_snv_indel.add_argument("--vcfanno_n_proc", default=4, type=int, 
                                    help = "Number of vcfanno processes (option '-p' in vcfanno), default: %(default)s")
    optional_snv_indel.add_argument("--retained_info_tags", dest="retained_info_tags", default="None", 
                                    help = "Comma-separated string of VCF INFO tags from query VCF that should be kept in PCGR output TSV file")
    optional_snv_indel.add_argument("--ignore_noncoding", action="store_true", 
                                    help = "Ignore non-coding (i.e. non protein-altering) variants in report, default: %(default)s")
    optional_snv_indel.add_argument("--vep_n_forks", default=4, type=int, 
                                    help = "Number of forks (VEP option '--fork'), default: %(default)s")
    optional_snv_indel.add_argument("--vep_buffer_size", default=500, type=int, 
                                    help = "Variant buffer size (variants read into memory simultaneously, VEP option '--buffer_size')\n- set lower to reduce memory usage, default: %(default)s")
    optional_snv_indel.add_argument("--vep_pick_order", default="mane_select,mane_plus_clinical,canonical,biotype,ccds,rank,tsl,appris,length", 
                                    help = "Comma-separated string " + \
        "of ordered transcript/variant properties for selection of primary variant consequence\n(option '--pick_order' in VEP), default: %(default)s")
    optional_snv_indel.add_argument("--vep_no_intergenic", action="store_true", 
                                    help = "Skip intergenic variants during variant annotation (VEP option '--no_intergenic' in VEP), default: %(default)s")
    optional_snv_indel.add_argument("--vep_regulatory", action="store_true", 
                                    help = "Add VEP regulatory annotations (VEP option '--regulatory') or non-coding interpretation, default: %(default)s")
    optional_snv_indel.add_argument("--vep_gencode_basic", action="store_true", 
                                    help = "Consider basic GENCODE transcript set only with Variant Effect Predictor (VEP) (VEP option '--gencode_basic').")

    optional_allelic_support.add_argument("--tumor_dp_tag", dest="tumor_dp_tag", default="_NA_", 
                                          help = "Specify VCF INFO tag for sequencing depth (tumor, must be Type=Integer, default: %(default)s)")
    optional_allelic_support.add_argument("--tumor_af_tag", dest="tumor_af_tag", default="_NA_", 
                                          help = "Specify VCF INFO tag for variant allelic fraction (tumor,  must be Type=Float, default: %(default)s)")
    optional_allelic_support.add_argument("--control_dp_tag", dest="control_dp_tag", default="_NA_", 
                                          help = "Specify VCF INFO tag for sequencing depth (control, must be Type=Integer, default: %(default)s)")
    optional_allelic_support.add_argument("--control_af_tag", dest="control_af_tag", default="_NA_", 
                                          help = "Specify VCF INFO tag for variant allelic fraction (control, must be Type=Float, default: %(default)s)")
    optional_allelic_support.add_argument("--call_conf_tag", dest="call_conf_tag", default="_NA_", 
                                          help = "Specify VCF INFO tag for somatic variant call confidence (must be categorical, e.g. Type=String, default: %(default)s)")
    optional_allelic_support.add_argument("--tumor_dp_min", type=int, default=None, dest="tumor_dp_min", 
                                          help = "If VCF INFO tag for sequencing depth (tumor) is specified and found, set minimum required depth for inclusion in report (default: %(default)s)")
    optional_allelic_support.add_argument("--tumor_ad_min", type=int, default=None, dest="tumor_ad_min",
                                        help=("If VCF INFO tag for sequencing depth (tumor) and variant allelic fraction (tumor) are specified and found, set minimum required\n "
                                              "allelic depth (tumor, i.e. number of reads supporting alternate allele) for inclusion in report (default: %(default)s)"))
    optional_allelic_support.add_argument("--tumor_af_min", type=float, default=None, dest="tumor_af_min", 
                                          help = "If VCF INFO tag for variant allelic fraction (tumor) is specified and found, set minimum required AF for inclusion in report (default: %(default)s)")
    optional_allelic_support.add_argument("--control_dp_min", type=int, default=None, dest="control_dp_min", 
                                          help = "If VCF INFO tag for sequencing depth (control) is specified and found, set minimum required depth for inclusion in report (default: %(default)s)")
    optional_allelic_support.add_argument("--control_ad_max", type=int, default=None, dest="control_ad_max",
                                        help=("If VCF INFO tag for sequencing depth (control) and variant allelic fraction (control) are specified and found, set maximum\n "
                                              "allowed allelic depth (control, i.e. number of reads supporting alternate allele) for inclusion in report (default: %(default)s)"))
    optional_allelic_support.add_argument("--control_af_max", type=float, default=None, dest="control_af_max", 
                                          help = "If VCF INFO tag for variant allelic fraction (control) is specified and found, set maximum tolerated AF for inclusion in report (default: %(default)s)")

    optional_tumor_only.add_argument("--pon_vcf", dest="pon_vcf", help = "VCF file with germline calls from Panel of Normals (PON) - i.e. blacklisted variants, (default: %(default)s)")
    maf_help_msg = "Exclude variants in tumor (SNVs/InDels, tumor-only mode only) with popmax MAF greater than this value"
    optional_tumor_only.add_argument("--gnomad_popmax_af_tolerated", dest="gnomad_popmax_af_tolerated", type=float, default=0.001, help=f"{maf_help_msg}, default: %(default)s)")
    optional_tumor_only.add_argument("--exclude_pon", action="store_true", 
                                     help = "Exclude variants occurring in PoN (Panel of Normals, if provided as VCF (--pon_vcf), default: %(default)s)")
    optional_tumor_only.add_argument("--exclude_likely_hom_germline", action="store_true", 
                                     help = "Exclude likely homozygous germline variants (allelic fraction of 1.0 for alternate allele in tumor - very unlikely somatic event), default: %(default)s)")
    optional_tumor_only.add_argument("--exclude_likely_het_germline", action="store_true", 
                                     help = "Exclude likely heterozygous germline variants (0.4-0.6 allelic fraction, AND presence in dbSNP + gnomAD\n, AND not existing as somatic record in COSMIC OR TCGA, default: %(default)s)")
    optional_tumor_only.add_argument("--exclude_clinvar_germline", action="store_true", 
                                     help = "Exclude variants found in ClinVar (of germline origin only), defult: %(default)s)")
    optional_tumor_only.add_argument("--exclude_dbsnp_nonsomatic", action="store_true", 
                                     help = "Exclude variants found in dbSNP (except for those present in ClinVar (somatic origin), or found in TCGA, or overlapping with COSMIC entries), defult: %(default)s)")
    optional_tumor_only.add_argument("--exclude_nonexonic", action="store_true", 
                                     help = "Exclude non-exonic variants, default: %(default)s)")

    optional_tmb_msi.add_argument("--estimate_tmb", action="store_true", 
                                  help = "Estimate tumor mutational burden from the total number of somatic mutations and target region size, default: %(default)s")
    optional_tmb_msi.add_argument("--tmb_display", dest="tmb_display", default="coding_and_silent", choices=["coding_and_silent", "coding_non_silent", "missense_only"], 
                                  help = "Type of TMB measure to show in report, default: %(default)s")
    optional_tmb_msi.add_argument("--tmb_dp_min", dest="tmb_dp_min", default=None, 
                                  help = "If VCF INFO tag for sequencing depth (tumor) is specified and found, set minimum required sequencing depth for TMB calculation: default: %(default)s")
    optional_tmb_msi.add_argument("--tmb_af_min", dest="tmb_af_min", default=None, 
                                  help = "If VCF INFO tag for allelic fraction (tumor) is specified and found, set minimum required allelic fraction for TMB calculation: default: %(default)s")
    optional_tmb_msi.add_argument("--tmb_ad_min", dest="tmb_ad_min", default=None, 
                                  help=("If VCF INFO tag for sequencing depth (tumor) and allelic fraction (tumor) are specified and found, \n"
                                        "set minimum required allelic depth (tumor, i.e. number of reads supporting alternate allele) for TMB calculation: default: %(default)s"))
    optional_tmb_msi.add_argument("--estimate_msi", action="store_true", 
                                  help = "Predict microsatellite instability status from patterns of somatic mutations/indels, default: %(default)s")

    optional_signatures.add_argument("--estimate_signatures", action="store_true", 
                                     help = "Estimate relative contributions of reference mutational signatures in query sample (re-fitting), default: %(default)s")
    optional_signatures.add_argument("--min_mutations_signatures", type=int, default=200, dest="min_mutations_signatures", 
                                     help = "Minimum number of SNVs required for re-fitting of mutational signatures (SBS) (default: %(default)s, minimum n = 100)")
    optional_signatures.add_argument("--all_reference_signatures", action="store_true", 
                                     help = "Use _all_ reference mutational signatures (SBS) during signature re-fitting rather than only those already attributed to the tumor type (default: %(default)s)")
    optional_signatures.add_argument("--include_artefact_signatures", action="store_true", 
                                     help = "Include sequencing artefacts in the collection of reference signatures (default: %(default)s")
    optional_signatures.add_argument("--prevalence_reference_signatures", type=float, default=0.1, 
                                     help = "Minimum tumor-type prevalence (in percent) of reference signatures to be included in refitting procedure (default: %(default)s)")

    optional_cna.add_argument("--cna_threshold_mode", dest="cna_threshold_mode", default="absolute", choices=["absolute", "relative", "combined"],
                              help=("Thresholding mode applied uniformly to all CNA tiers "
                                    "(amplification, gain, heterozygous deletion), based on either\n "
                                    "absolute total copy number ('absolute'), fold-change over tumor "
                                    "ploidy ('relative'), or a combination of both ('combined') "
                                    "(default: %(default)s)"))
    optional_cna.add_argument("--cna_amp_threshold_absolute", type=int, default=5, dest="cna_amp_threshold_absolute",
                              help = "Absolute total copy number threshold for amplification calls (default: %(default)s)")
    optional_cna.add_argument("--cna_amp_threshold_relative", type=float, default=2.5, dest="cna_amp_threshold_relative",
                              help = "Fold-change over tumor ploidy threshold for amplification calls (default: %(default)s)")
    optional_cna.add_argument("--cna_gain_threshold_absolute", type=int, default=3, dest="cna_gain_threshold_absolute",
                              help = "Absolute total copy number threshold for gain calls (default: %(default)s)")
    optional_cna.add_argument("--cna_gain_threshold_relative", type=float, default=1.5, dest="cna_gain_threshold_relative",
                              help = "Fold-change over tumor ploidy threshold for gain calls (default: %(default)s)")
    optional_cna.add_argument("--cna_del_threshold_absolute", type=int, default=1, dest="cna_del_threshold_absolute",
                              help = "Absolute total copy number threshold for heterozygous deletion calls (default: %(default)s)")
    optional_cna.add_argument("--cna_del_threshold_relative", type=float, default=0.5, dest="cna_del_threshold_relative",
                              help = "Fold-change below tumor ploidy threshold for heterozygous deletion calls (default: %(default)s)")
    optional_cna.add_argument("--cna_transcript_overlap_pct", type = float, default = 50, dest="cna_transcript_overlap_pct",
                              help = "Mean percent overlap between copy number segment and gene transcripts for reporting of gains/losses in tumor suppressor genes/oncogenes, (default: %(default)s)")

    optional_rna.add_argument("--fusion_min_split_reads", type=int, default=3, dest="fusion_min_split_reads", 
                              help = "Minimum number of split reads supporting a fusion event (default: %(default)s, minimum: 2)")
    optional_rna.add_argument('--expression_sim', action='store_true', help = "Compare expression profile of tumor sample to expression profiles of other tumor samples (default: %(default)s)")
    optional_rna.add_argument("--expression_sim_db", dest = "expression_sim_db", default="tcga,depmap,treehouse", 
                              help = "Comma-separated string " + \
        "of databases for used in RNA expression similarity analysis, default: %(default)s")

    optional_germline.add_argument("--input_cpsr", dest="input_cpsr", 
                                   help = "CPSR-classified germline calls (file '<cpsr_sample_id>.cpsr.<genome_assembly>.classification.tsv.gz')")
    optional_germline.add_argument("--input_cpsr_yaml", dest="input_cpsr_yaml", 
                                   help = "CPSR YAML configuration file (file '<cpsr_sample_id>.cpsr.<genome_assembly>.conf.yaml')")
    optional_germline.add_argument("--cpsr_ignore_vus", action="store_true", 
                                   help = "Do not show variants of uncertain significance (VUS) in the germline section of the HTML report (default: %(default)s)")

    optional_biomarker.add_argument("--oncokb_api_token", dest="oncokb_api_token", default = None,
                                    help = "Oncokb API token (default: ")
    optional_biomarker.add_argument("--oncokb_oncotree_code", dest="oncokb_oncotree_code",
                                    help = "OncoKB oncotree code (default: %(default)s)")
    optional_biomarker.add_argument("--oncokb_exclusive", action="store_true",
                                    help = "Limit biomarker reporting to OncoKB only - skip CIViC and CGI sources (default: %(default)s)")
    optional_biomarker.add_argument("--oncokb_maf_query_all", action="store_true",
                                    help = "Query OncoKB for all variant classes, including non-coding variants (IGR, Intron, UTR, flanking regions). "
                                           "By default, non-coding variants are filtered out before OncoKB annotation to reduce processing time. "
                                           "Intended for TARGETED/WES assays only - enabling for WGS may result in very long MafAnnotator.py runtimes (default: %(default)s)")
    optional_other.add_argument("--force_overwrite", action="store_true", 
                                help = "By default, the script will fail with an error if any output file already exists.\nYou can force the overwrite of existing result files by using this flag, default: %(default)s")
    optional_other.add_argument("--no_reporting", action="store_true", 
                                help = "Run functional variant annotation on VCF through VEP/vcfanno, omit other analyses\n(i.e. Tier assignment/MSI/TMB/Signatures etc. and report generation (STEP 4), default: %(default)s")
    optional_other.add_argument("--no_html", action="store_true", 
                                help = "Do not generate HTML report (default: %(default)s)")
    optional_other.add_argument("--debug", action="store_true", 
                                help = "Print full commands to log")
    optional_other.add_argument("--pcgrr_conda", default="pcgrr", 
                                help = "pcgrr conda env name (default: %(default)s)")
    optional_other.add_argument("--version", action="version", version="%(prog)s " + str(pcgr_vars.PCGR_VERSION))

    # Parse required and optional arguments
    args = parser.parse_args()
    arg_dict = vars(args)

    # verify parsed arguments
    print('')
    logger = getlogger('pcgr-verify-args')
    logger.info("PCGR - VERIFY ARGUMENTS SECTION: Verifying input arguments and configuring workflow")
    arg_dict_verified = arg_checker.verify_args(arg_dict)
    logger.info('Finished pcgr-verify_args')

    # create config options
    logger = getlogger('pcgr-create-config')
    conf_options = create_config(arg_dict_verified, workflow = "PCGR", logger = logger)

    # Verify existence of input files, define and check existence of output files
    input_data = arg_checker.verify_input_files(arg_dict_verified)
    output_data = arg_checker.define_output_files(arg_dict_verified)
    # Run PCGR workflow
    run_pcgr(input_data, output_data, conf_options)

def run_pcgr(input_data, output_data, conf_options):
    """
    Main function to run the PCGR workflow
    """

    debug = conf_options['debug']

    # set basic run commands
    #clinical_trials_set = 'ON' if conf_options['clinicaltrials']['run'] else 'OFF'
    msi_prediction_set = 'ON' if conf_options['somatic_snv']['msi']['run'] else 'OFF'
    msig_estimation_set = 'ON' if conf_options['somatic_snv']['mutational_signatures']['run'] else 'OFF'
    tmb_estimation_set = 'ON' if conf_options['somatic_snv']['tmb']['run'] else 'OFF'
    rnaseq_sim_analysis_set = 'ON' if conf_options['expression']['similarity_analysis'] else 'OFF'
    run_vcf2maf = conf_options['other']['vcf2maf']
    assay_mode = 'Tumor vs. Control'
    oncogenicity_annotation = 1
    if conf_options['assay_properties']['vcf_tumor_only']:
        assay_mode = 'Tumor-only'
    NCBI_BUILD_MAF = pcgr_vars.NCBI_BUILD_MAF
    if conf_options['genome_assembly'] == 'grch37':
        NCBI_BUILD_MAF = 'GRCh37'
    logger = getlogger('pcgr-get-OS')


    # Define input and output files
    input_vcf = 'None'
    input_cna = 'None'
    input_rna_fusion = 'None'
    input_rna_expression = 'None'
    input_cpsr_calls = 'None'
    input_cpsr_yaml = 'None'
    pon_vcf = 'None'
   
    # Update input data variables with data provided by user
    if input_data['vcf_basename'] != 'NA':
        input_vcf = os.path.join(input_data['vcf_dir'], input_data['vcf_basename'])
    if input_data['cna_basename'] != 'NA':
        input_cna = os.path.join(input_data['cna_dir'], input_data['cna_basename'])
    if input_data['rna_fusion_basename'] != 'NA':
        input_rna_fusion = os.path.join(input_data['rna_fusion_dir'], input_data['rna_fusion_basename'])
    if input_data['rna_expression_basename'] != 'NA':
        input_rna_expression = os.path.join(input_data['rna_expression_dir'], input_data['rna_expression_basename'])
    if input_data['germline_basename'] != 'NA':
        input_cpsr_calls = os.path.join(input_data['germline_dir'], input_data['germline_basename'])
    if input_data['germline_yaml_basename'] != 'NA':
        input_cpsr_yaml = os.path.join(input_data['germline_yaml_dir'], input_data['germline_yaml_basename'])
    if input_data['pon_vcf_basename'] != 'NA':
        pon_vcf = os.path.join(input_data['pon_vcf_dir'], input_data['pon_vcf_basename'])

    n_vcf_pass_variants = 0
    pon_annotation = 0
    variant_set = pd.DataFrame
    expression_data = None
    output_dir = output_data['dir']
    output_prefix = output_data['prefix']
    utils.safe_makedir(output_dir)

    random_id = random_id_generator(15)
    # Define temporary output files
    input_vcf_validated =             f'{output_prefix}.{random_id}.ready.vcf.gz'
    input_vcf_validated_uncompr =     f'{output_prefix}.{random_id}.ready.vcf'
    vep_vcf =                         f'{output_prefix}.{random_id}.vep.vcf'
    vep_vcfanno_vcf =                 f'{output_prefix}.{random_id}.vep.vcfanno.vcf'
    vep_vcfanno_summarised_vcf =      f'{output_prefix}.{random_id}.vep.vcfanno.summarised.vcf'
    vep_vcfanno_summarised_pass_vcf = f'{output_prefix}.{random_id}.vep.vcfanno.summarised.pass.vcf'
    output_vcf =                      f'{output_prefix}.vcf.gz'
    output_pass_vcf =                 f'{output_prefix}.pass.vcf.gz'
    output_pass_vcf2tsv =             f'{output_prefix}.pass.vcf2tsv.tsv'
    output_pass_vcf2tsv_gz =          f'{output_pass_vcf2tsv}.gz'
    output_pass_tsv =                 f'{output_prefix}.pass.tsv'
    output_pass_tsv_gz =              f'{output_pass_tsv}.gz'
    output_pass_raw_tsv_gz =          f'{output_prefix}.pass.raw.tsv.gz'
    output_tmp_maf =                  f'{output_prefix}.tmp.maf'
    output_maf =                      f'{output_prefix}.maf'
    output_vcf2maf_log =              f'{output_prefix}.maf.log'
    vep_vcf_with_varid =              f'{output_prefix}.varid.vcf'
    oncokb_input_fusion_tsv =         f'{output_prefix}.oncokb_input.fusions.txt'
    oncokb_input_cna_tsv =            f'{output_prefix}.oncokb_input.cna.txt'
    oncokb_input_maf =                f'{output_prefix}.oncokb_input.maf'
    oncokb_output_fusion_tsv =        f'{output_prefix}.oncokb_output.fusion.txt'
    oncokb_output_maf_hgvsg =         f'{output_prefix}.oncokb_output.maf.hgvsg.txt'
    oncokb_output_maf_hgvsp =         f'{output_prefix}.oncokb_output.maf.hgvsp.txt'
    oncokb_output_cna_tsv =           f'{output_prefix}.oncokb_output.cna.txt'
    yaml_fname =                      f'{output_prefix}.conf.yaml'
    tmb_fname =                       f'{output_prefix}.tmb.tsv'

    # PCGR|Generate YAML file - containing configuration options and paths to first-pass annotation of molecular profile datasets
    # - VCF/TSV files (somatic SNVs/InDels)
    # - TSV file (germline SNVs/InDels - CPSR)
    # - TSV files (somatic copy number aberrations)
    # - TSV files (TMB)
    # - TSV files (RNA expression)
    # - TSV files (RNA fusion)
    logger = getlogger('pcgr-write-yaml')

    # update conf_options with paths to output files
    conf_options['output_dir'] = output_dir
    conf_options['output_prefix'] = output_prefix
    if not input_vcf == 'None':
        conf_options['molecular_data']['fname_mut_vcf'] = output_vcf
        conf_options['molecular_data']['fname_mut_tsv'] = output_pass_tsv_gz
        if conf_options['other']['vcf2maf'] == 1:
            conf_options['molecular_data']['fname_maf_tsv'] = output_maf
        if conf_options['somatic_snv']['tmb']['run'] == 1:
            conf_options['molecular_data']['fname_tmb_tsv'] = tmb_fname
    if not input_cna == 'None':
        conf_options['molecular_data']['fname_cna_gene_tsv'] = output_data['cna_gene']
        conf_options['molecular_data']['fname_cna_segment_tsv'] = output_data['cna_segment']
    if not input_cpsr_calls == 'None':
        conf_options['molecular_data']['fname_germline_tsv'] = input_cpsr_calls
    if not input_rna_fusion == 'None':
        conf_options['molecular_data']['fname_rna_fusion_tsv'] = output_data['rna_fusion']
    if not input_cpsr_yaml == 'None':
        conf_options['molecular_data']['fname_germline_yaml'] = input_cpsr_yaml
    if not input_rna_expression == 'None':
        conf_options['molecular_data']['fname_expression_tsv'] = output_data['expression']
        conf_options['molecular_data']['fname_expression_outliers_tsv'] = output_data['expression_outliers']
        if conf_options['expression']['similarity_analysis'] == 1:
            conf_options['molecular_data']['fname_expression_similarity_tsv'] = output_data['expression_similarity']
    ## if oncokb is to be run, add output files for oncokb annotations of SNVs/Indels (MAF, hgvsp + genomic), fusions and CNAs
    if conf_options['oncokb']['run'] == 1:
        ## i want the output, not the input, to be in the YAML file, since the output is what will be used for the oncokb annotation and downstream analyses
        conf_options['molecular_data']['fname_oncokb_output_maf_hgvsg'] = oncokb_output_maf_hgvsg
        conf_options['molecular_data']['fname_oncokb_output_maf_hgvsp'] = oncokb_output_maf_hgvsp
        conf_options['molecular_data']['fname_oncokb_output_fusions'] = oncokb_output_fusion_tsv
        conf_options['molecular_data']['fname_oncokb_output_cna'] = oncokb_output_cna_tsv
            
            
    # make YAML file
    yaml_data = populate_config_data(conf_options, input_data["refdata_assembly_dir"],
                                     workflow = "PCGR", logger = logger)
    genome_assembly = yaml_data['genome_assembly']

    # PCGR|validate_input - verify that VCF and CNA segment file is of appropriate format
    logger = getlogger("pcgr-validate-input-arguments")
    print('----')
    logger.info("PCGR - INPUT VALIDATION SECTION: Validate input data and options")

    vcf_validate_command = (
            f'validate_input_pcgr.py '
            f'{input_data["refdata_assembly_dir"]} '
            f'{input_vcf} '
            f'{input_vcf_validated_uncompr} '
            f'{input_cna} '
            f'{input_rna_fusion} '
            f'{input_rna_expression} '
            f'{input_cpsr_calls} '
            f'{pon_vcf} '
            f'{conf_options["assay_properties"]["vcf_tumor_only"]} '
            f'{yaml_data["sample_id"]} '
            f'{genome_assembly} '
            f'{conf_options["other"]["retained_vcf_info_tags"]} '
            f'{conf_options["somatic_snv"]["allelic_support"]["tumor_dp_tag"]} '
            f'{conf_options["somatic_snv"]["allelic_support"]["tumor_af_tag"]} '
            f'{conf_options["somatic_snv"]["allelic_support"]["control_dp_tag"]} '
            f'{conf_options["somatic_snv"]["allelic_support"]["control_af_tag"]} '
            f'{conf_options["somatic_snv"]["allelic_support"]["call_conf_tag"]} '
            f'{conf_options["somatic_snv"]["tumor_only"]["exclude_likely_hom_germline"]} '
            f'{conf_options["somatic_snv"]["tumor_only"]["exclude_likely_het_germline"]} '
            f'--output_dir {output_dir} '
            f'{"--debug " if debug else ""}'
            f'{"--keep_uncompressed" if run_vcf2maf else ""} '
            )
    check_subprocess(logger, vcf_validate_command, debug)
    logger.info('Finished pcgr-validate-input-arguments')
    print('----')

    # PCGR|start - Log key information about sample, options and sequencing assay/design
    logger = getlogger('pcgr-settings')
    logger.info('--- Personal Cancer Genome Reporter workflow ----')
    logger.info(f'Sample name: {yaml_data["sample_id"]}')
    if conf_options['sample_properties']['site'] == 'Any':
        logger.info('Tumor type: Cancer_NOS (Any tumortype)')
    else:
        logger.info(f'Tumor type: {conf_options["sample_properties"]["site"]}')
    
    ## Log information about which molecular input types are being providd as input
    logger.info(f'Molecular input data - somatic SNVs/InDels: {"OFF" if input_vcf == "None" else "ON"}')
    logger.info(f'Molecular input data - somatic copy number alterations: {"OFF" if input_cna == "None" else "ON"}')
    logger.info(f'Molecular input data - tumor RNA fusion transcripts: {"OFF" if input_rna_fusion == "None" else "ON"}')
    logger.info(f'Molecular input data - tumor RNA expression: {"OFF" if input_rna_expression == "None" else "ON"}')
    logger.info(f'Molecular input data - germline SNVs/InDels (CPSR): {"OFF" if input_cpsr_calls == "None" else "ON"}')
    
    ## Log information about sequencing assay and variant filtering settings (if VCF input is provided)
    if not input_vcf == 'None':
        logger.info(f'Sequencing assay - type: {conf_options["assay_properties"]["type"]}')
        logger.info(f'Sequencing assay - mode: {assay_mode}')
        logger.info((
            f'Sequencing assay - effective (coding) target size: '
            f'{conf_options["assay_properties"]["effective_target_size_mb"]}Mb'))
        if conf_options['somatic_snv']['allelic_support']['tumor_dp_tag'] != '_NA_':
            logger.info((
                f'Variant filtering settings - minimum sequencing depth tumor: '
                f'{conf_options["somatic_snv"]["allelic_support"]["tumor_dp_min"]}'))
        if conf_options['somatic_snv']['allelic_support']['tumor_af_tag'] != '_NA_':
            logger.info((
                f'Variant filtering settings - minimum allelic fraction tumor: '
                f'{conf_options["somatic_snv"]["allelic_support"]["tumor_af_min"]}'))
            logger.info((
                f'Variant filtering settings - minimum allelic depth tumor: '
                f'{conf_options["somatic_snv"]["allelic_support"]["tumor_ad_min"]}'))
        if conf_options['somatic_snv']['allelic_support']['control_dp_tag'] != '_NA_':
            logger.info((
                f'Variant filtering settings - minimum sequencing depth control: '
                f'{conf_options["somatic_snv"]["allelic_support"]["control_dp_min"]}'))
        if conf_options['somatic_snv']['allelic_support']['control_af_tag'] != '_NA_':
            logger.info((
                f'Variant filtering settings - maximum allelic fraction control: '
                f'{conf_options["somatic_snv"]["allelic_support"]["control_af_max"]}'))
            logger.info((
                f'Variant filtering settings - maximum allelic depth control: '
                f'{conf_options["somatic_snv"]["allelic_support"]["control_ad_max"]}'))
    
    logger.info(f'Genome assembly: {genome_assembly}')
    logger.info(f'Mutational signature estimation: {msig_estimation_set}')
    logger.info(f'MSI classification: {msi_prediction_set}')
    logger.info(f'Mutational burden estimation: {tmb_estimation_set}')
    logger.info(f'RNA expression similarity analysis: {rnaseq_sim_analysis_set}')
    #print('----')
    
    #logger.info(f'Include molecularly targeted clinical trials (beta): {clinical_trials_set}')

    if not input_vcf == 'None':
        
        vep_command = get_vep_command(
            file_paths = input_data,
            conf_options = yaml_data,
            input_vcf = input_vcf_validated,
            output_vcf = vep_vcf)
    
        # PCGR|VEP - run consequence annotation with Variant Effect Predictor
        print('----')       
        logger = getlogger('pcgr-vep')
        logger.info(f'PCGR - SNV/INDEL ANALYSIS SECTION - VEP: Variant annotation with Variant Effect Predictor {pcgr_vars.VEP_VERSION}' + \
                    f', GENCODE release {pcgr_vars.GENCODE_VERSION[genome_assembly]}, genome assembly {yaml_data["genome_assembly"]}')
        logger.info('VEP configuration - one primary consequence block pr. alternative gene allele (--flag_pick_allele_gene)')
        logger.info(f'VEP configuration - transcript pick order: {yaml_data["conf"]["vep"]["vep_pick_order"]}')
        logger.info('VEP configuration - transcript pick order: See more at https://www.ensembl.org/info/docs/tools/vep/script/vep_other.html#pick_options')
        logger.info(f'VEP configuration - GENCODE set: {vep_command["gencode_set_in_use"]}')
        logger.info(f'VEP configuration - skip intergenic variants: {"ON" if yaml_data["conf"]["vep"]["vep_no_intergenic"] == 1 else "OFF"}')
        logger.info(f'VEP configuration - regulatory variant annotation: {"ON" if yaml_data["conf"]["vep"]["vep_regulatory"] == 1 else "OFF"}')
        logger.info((
            f'VEP configuration - buffer size/number of forks: '
            f'{yaml_data["conf"]["vep"]["vep_buffer_size"]}/{yaml_data["conf"]["vep"]["vep_n_forks"]}'))
        logger.info(f'VEP - plugins in use: {vep_command["plugins_in_use"]}')

        check_subprocess(logger, vep_command['main'], debug)
        check_subprocess(logger, vep_command['bgzip'], debug)
        check_subprocess(logger, vep_command['tabix'], debug)
        logger.info('Finished pcgr-vep')
        print('----')

        # PCGR|vcf2maf - if option set, convert VCF to MAF with https://github.com/mskcc/vcf2maf
        if run_vcf2maf:

            retained_info_tags = \
                yaml_data["conf"]['somatic_snv']['allelic_support']['control_dp_tag'] + ',' \
                + yaml_data["conf"]['somatic_snv']['allelic_support']['tumor_dp_tag'] + ',' \
                + yaml_data["conf"]['somatic_snv']['allelic_support']['control_af_tag'] + ',' \
                + yaml_data["conf"]['somatic_snv']['allelic_support']['tumor_af_tag']

            retained_info_tags = re.sub(r'_NA_(,)?','', retained_info_tags)
            # Always retain VAR_ID so the MAF can be joined back to VCF-derived data in pcgrr
            retained_info_tags = ('VAR_ID,' + retained_info_tags).rstrip(',')

            logger = getlogger("pcgr-vcf2maf")
            logger.info("PCGR - SNV/INDEL ANALYSIS SECTION - VCF2MAF: vcf2maf conversion of VEP-annotated VCF")

            # VEP writes ##contig header lines with chr-prefixed IDs (from the GRCh38
            # reference cache) even though variant records use bare chromosome names
            # (the chr prefix was stripped during input normalisation). htslib/cyvcf2
            # emits "[W::vcf_parse] Contig '1' is not defined in the header" warnings
            # the moment the file is opened. Fix the header BEFORE any tool touches
            # the file by stripping chr from ##contig IDs.
            vep_vcf_reheadered = vep_vcf.replace('.vcf', '.reheadered.vcf')
            check_subprocess(
                logger,
                f"sed 's/^##contig=<ID=chr/##contig=<ID=/' {vep_vcf} > {vep_vcf_reheadered}",
                debug)

            # Stamp VAR_ID (CHROM_POS_REF_ALT, VCF-format) onto the reheadered VEP VCF
            # before vcf2maf strips anchor bases from indel alleles
            add_var_id_to_vcf(vep_vcf_reheadered, vep_vcf_with_varid, logger)
            vcf2maf_input = vep_vcf_with_varid if os.path.exists(vep_vcf_with_varid) else vep_vcf_reheadered

            update_allelic_support = True
            logger.info('Converting VEP-annotated VCF to MAF with https://github.com/mskcc/vcf2maf')
            vcf2maf_command = (
                    f'vcf2maf.pl --retain-info {retained_info_tags} --inhibit-vep --input-vcf {vcf2maf_input} '
                    f'--tumor-id {yaml_data["sample_id"]} --output-maf {output_tmp_maf} '
                    f'--ref-fasta {vep_command["fasta_assembly"]} '
                    f'--ncbi-build {NCBI_BUILD_MAF} > {output_vcf2maf_log} 2>&1'
                    )
            check_subprocess(logger, vcf2maf_command, debug)
            if not debug:
                remove_file(output_vcf2maf_log)

            ## add information on allelic support in MAF file
            ## (n_depth, n_ref_count, n_alt_count, t_depth, t_ref_count, t_alt_count)
            update_maf(
                maf_tmp_fname = output_tmp_maf,
                maf_fname = output_maf,
                allelic_support_tags = yaml_data["conf"]['somatic_snv']['allelic_support'],
                logger = logger,
                update_allelic_support = update_allelic_support,
                debug = debug
            )

            # Generate OncoKB-compliant MAF input file
            logger.info('Generating OncoKB-API compliant MAF input file')
            generate_oncokb_maf_input(
                maf_fname = output_maf,
                output_fname = oncokb_input_maf,
                oncokb_maf_query_all = bool(conf_options['oncokb']['maf_query_all']),
                logger = logger
            )

            if not debug:
                remove_file(vep_vcf_reheadered)
                remove_file(vep_vcf_with_varid)

            logger.info('Finished pcgr-vep-vcf2maf')
            print('----')

        # PCGR|vcfanno - annotate query VCF against a number of variant annotation tracks (BED/VCF)
        logger = getlogger("pcgr-vcfanno")
        pcgr_vcfanno_command = (
                f'pcgr_vcfanno.py {vep_vcf}.gz {vep_vcfanno_vcf} {input_data["refdata_assembly_dir"]} '
                f'--num_processes {conf_options["other"]["vcfanno_n_proc"]} '
                f'--dbnsfp --clinvar --rmsk --winmsk --simplerepeat '
                f'--tcga --gene_transcript_xref --dbmts --gwas '
                f'{"--debug" if debug else ""}'
                )
        vcfanno_db_src_msg1 = (
                f"Annotation sources I (vcfanno): {'Panel-of-Normals, ' if pon_vcf != 'None' else ''}ClinVar, dbNSFP, "
                f"dbMTS, TCGA, GWAS catalog"
                )
        vcfanno_db_src_msg2 = \
                "Annotation sources II (vcfanno): RepeatMasker, SimpleRepeats, WindowMaskerSDust"

        logger.info("PCGR - SNV/INDEL ANALYSIS SECTION - VCFANNO: Variant annotation for cancer precision medicine with pcgr-vcfanno")
        logger.info(vcfanno_db_src_msg1)
        logger.info(vcfanno_db_src_msg2)
        if pon_vcf != "None":
            pon_annotation = 1
            pcgr_vcfanno_command += f' --pon_vcf {pon_vcf}'
        check_subprocess(logger, pcgr_vcfanno_command, debug)
        logger.info("Finished pcgr-vcfanno")
        print('----')

        # PCGR|pcgr_summarise - expand annotations in VCF file
        logger = getlogger("pcgr-summarise")
        pcgr_summarise_command = (
                f'pcgr_summarise.py {vep_vcfanno_vcf}.gz {vep_vcfanno_summarised_vcf} {pon_annotation} '
                f'{yaml_data["conf"]["vep"]["vep_regulatory"]} {oncogenicity_annotation} '
                f'{yaml_data["conf"]["sample_properties"]["site2"]} {yaml_data["genome_assembly"]} {yaml_data["conf"]["vep"]["vep_pick_order"]} '
                f'{input_data["refdata_assembly_dir"]} --compress_output_vcf '
                f'{"--debug" if debug else ""}'
                )
        summarise_db_src_msg1 = \
                "Annotation sources: cancerhotspots.org, CIViC, Cancer Biomarkers database (CGI), Network of Cancer Genes (NCG)"
        summarise_db_src_msg2 = \
                "Annotation sources: CancerMine, IntOGen, TCGA driver genes"

        logger.info("PCGR - SNV/INDEL ANALYSIS SECTION - SUMMARY: Variant and cancer gene annotations with pcgr-summarise")
        logger.info(summarise_db_src_msg1)
        logger.info(summarise_db_src_msg2)

        logger.info('Variant oncogenicity classification according to ClinGen/CGC/VICC standard operating procedures (Horak et al., Genet Med, 2022)')
        logger.info('Variant biomarker matching (CIViC, CGI) at multiple resolutions (genes, exons, amino acid positions, hgvsp/hgvsc, genomic)')
        logger.info('Tumor suppressor/oncogene annotations based on multiple sources (NCG, CGC, CancerMine)')
        check_subprocess(logger, pcgr_summarise_command, debug)
        #exit(0)

        # PCGR|clean - move output files and clean up temporary files
        os.rename(f'{vep_vcfanno_summarised_vcf}.gz', output_vcf)
        os.rename(f'{vep_vcfanno_summarised_vcf}.gz.tbi', f'{output_vcf}.tbi')
        os.rename(f'{vep_vcfanno_summarised_pass_vcf}.gz', output_pass_vcf)
        os.rename(f'{vep_vcfanno_summarised_pass_vcf}.gz.tbi', f'{output_pass_vcf}.tbi')
        delete_files = (
                glob(f'{vep_vcf}*') +
                glob(f'{vep_vcfanno_summarised_vcf}') +
                glob(f'{vep_vcfanno_summarised_pass_vcf}*') +
                glob(f'{vep_vcfanno_vcf}*') +
                glob(f'{input_vcf_validated_uncompr}*')
                )
        # do not delete if debugging
        if not debug:
            for fn in delete_files:
                remove_file(fn)

        logger.info('Finished pcgr-summarise main command')

        # PCGR|vcf2tsvpy - convert VCF to TSV with https://github.com/sigven/vcf2tsvpy
        pcgr_vcf2tsv_command = f'vcf2tsvpy --input_vcf {output_pass_vcf} --out_tsv {output_pass_vcf2tsv} --compress'
        logger.info("Converting VCF to TSV with https://github.com/sigven/vcf2tsvpy")
        check_subprocess(logger, pcgr_vcf2tsv_command, debug)

        ## Append additional (space-containing) annotations not suitable for VCF INFO
        logger.info("Appending ClinVar traits, official gene names, and protein domain annotations")
        variant_set = \
           append_annotations(
              output_pass_vcf2tsv_gz, refdata_assembly_dir = input_data["refdata_assembly_dir"], logger = logger)
        ## Set allelic support properties (DP_TUMOR, DP_CONTROL, VAF_TUMOR, VAF_CONTROL)
        # Guard: append_annotations can return None; fail early instead of passing None to functions
        if variant_set is None:
            return error_message(
                "Failed to append annotations: variant set is empty (append_annotations returned None).",
                logger
            )
        variant_set = set_allelic_support(variant_set, allelic_support_tags = yaml_data["conf"]['somatic_snv']['allelic_support'], logger = logger)
        ## Clean annotations (formatting etc.)
        variant_set = clean_annotations(variant_set, yaml_data, logger = logger)

        ## Check if AD/DP properties could be detected/pulled from VCFs
        for c in ['DP_TUMOR', 'DP_CONTROL', 'VAF_TUMOR', 'VAF_CONTROL']:
            if c in variant_set.columns:
                if len(variant_set[variant_set[c].isnull()]) == 0:
                    var = str(c).lower() + '_detected'
                    if var in yaml_data['conf']['sample_properties'].keys():
                        yaml_data['conf']['sample_properties'][var] = 1

        ## Write final somatic SNV/InDel variant set to gzipped TSV
        n_vcf_pass_variants = len(variant_set)
        pd_to_csv(
            df = variant_set,
            fname = output_pass_tsv_gz,
            col_sep = "\t"            
        )
        if not debug:
            remove_file(output_pass_vcf2tsv_gz)

        ## Reduce output file for WGS/WES when the number of variants exceeds MAX_VARIANTS_FOR_REPORT
        if yaml_data["conf"]['assay_properties']['type'] in ('WGS', 'WES'):
            reduce_variants_for_report(
                output_pass_tsv_gz=output_pass_tsv_gz,
                output_pass_raw_tsv_gz=output_pass_raw_tsv_gz,
                yaml_data=yaml_data,
                debug=debug,
                logger=logger,
            )

        # PCGR|TMB - calculate TMB in multiple ways - potentially also considering subset of variants (depth and allele frequency filtered)
        if yaml_data['conf']['somatic_snv']['tmb']['run'] == 1:
            logger_tmb = getlogger('pcgr-calculate-tmb')
            calculate_tmb(
                variant_set = variant_set,
                tumor_dp_min = int(yaml_data['conf']['somatic_snv']['tmb']['tmb_dp_min']),
                tumor_af_min = float(yaml_data['conf']['somatic_snv']['tmb']['tmb_af_min']),
                tumor_ad_min = int(yaml_data['conf']['somatic_snv']['tmb']['tmb_ad_min']),
                target_size_mb = float(yaml_data['conf']['assay_properties']['effective_target_size_mb']),
                sample_id = yaml_data['sample_id'],
                tmb_fname = tmb_fname,
                logger = logger_tmb)

        logger = getlogger('pcgr-summarise')
        logger.info('Finished pcgr-summarise')
        print('----')

    # PCGR|Expression - Gene expression (bulk RNA-seq) analyses
    if input_rna_expression != 'None' and yaml_data['molecular_data']['fname_expression_tsv'] != 'None':
        
        logger = getlogger('pcgr-expression-analysis')
        logger.info('PCGR - EXPRESSION ANALYSIS SECTION: Gene expression analysis, similarity analysis, and outlier detection')

        
        expression_data = parse_expression(input_rna_expression, yaml_data, logger=logger)
        if 'transcript' in expression_data.keys():
            if expression_data['transcript'] is not None:
                pd_to_csv(df=expression_data['transcript'],
                          fname=yaml_data['molecular_data']['fname_expression_tsv'],
                          col_sep="\t")
            else:
                if 'gene' in expression_data.keys():
                    if expression_data['gene'] is not None:
                        pd_to_csv(df=expression_data['gene'],
                                  fname=yaml_data['molecular_data']['fname_expression_tsv'],
                                  col_sep="\t")
        
        
        if 'transcript' in expression_data.keys():
            ## check if variant set is not None and non-empty
            if expression_data['transcript'] is not None and variant_set is not None and not variant_set.empty:
                variant_set = aggregate_tpm_per_cons(variant_set, expression_data, logger = logger)

            ## Merge expression data with somatic SNV/InDel variant set
            if expression_data['transcript'] is not None and variant_set is not None and not variant_set.empty:
                variant_set = integrate_variant_expression(
                    variant_set, expression_data, logger = logger)
        
        ## Find expression outliers by comparing with reference datasets
        ## - Preliminary version only compares with TCGA datasets, and only selects specific cohorts (pcgr_vars.SITE_TO_DISEASE)
        logger.info('Identification of expression outliers through comparison with reference datasets')
        expression_outliers = find_expression_outliers(
            expression_data,
            yaml_data,
            protein_coding_only = True,
            logger = logger)
        
        if not expression_outliers.empty:
            ## Write to gzipped TSV file
            pd_to_csv(
                df = expression_outliers,
                fname = yaml_data['molecular_data']['fname_expression_outliers_tsv'],
                col_sep = "\t"                
            )            
        else:
            yaml_data['molecular_data']['fname_expression_outliers_tsv'] = 'None'

        ## Correlate sample gene expression profile with reference samples (expression-based similarity)
        if yaml_data['conf']['expression']['similarity_analysis'] == 1:
            exp_similarity_results = correlate_sample_expression(
                expression_data,
                yaml_data,
                protein_coding_only = True,
                logger = logger)

            ## Aggregate similarity results across sources (TCGA, TreeHouse, Depmap) into single table
            exp_similarity_results_all = pd.DataFrame()
            for source in exp_similarity_results.keys():
                if not exp_similarity_results[source].empty:
                    exp_similarity_results_all = \
                        pd.concat([exp_similarity_results_all,
                                   exp_similarity_results[source]], axis = 0)

            ## Write to gzipped TSV file
            if not exp_similarity_results_all.empty:
                pd_to_csv(
                    df = exp_similarity_results_all,
                    fname = yaml_data['molecular_data']['fname_expression_similarity_tsv'],
                    col_sep = "\t",                    
                )
        
        print('----')
    else:
        logger = getlogger("pcgr-expression-analysis")
        logger.info('PCGR - RNA EXPRESSION ANALYSIS SECTION: Gene expression analysis - OMITTED (no data available)')
        print('----')

    # PCGR|CNA - Annotate allele-specific copy number segments with cytobands, overlapping transcripts, and biomarkers
    if not input_cna == 'None':
        logger = getlogger("pcgr-annotate-cna-segments")
        logger.info('PCGR - CNA ANALYSIS SECTION: Annotation of copy number segments - cytobands, overlapping transcripts, and biomarkers')

        cna_annotation = cna.annotate_cna_segments(
            input_cna_segment_fname = input_cna,
            output_segment_gene_fname = output_data['cna_gene'],
            output_segment_fname = output_data['cna_segment'],
            oncokb_input_fname = oncokb_input_cna_tsv,
            output_dir = output_data['dir'],
            build = yaml_data['genome_assembly'],
            sample_id = yaml_data['sample_id'],
            sex = yaml_data['conf']['sample_properties']['sex'],
            refdata_assembly_dir = input_data['refdata_assembly_dir'],
            amp_threshold_absolute = yaml_data["conf"]['somatic_cna']['amp_threshold_absolute'],
            amp_threshold_relative = yaml_data["conf"]['somatic_cna']['amp_threshold_relative'],
            gain_threshold_absolute = yaml_data["conf"]['somatic_cna']['gain_threshold_absolute'],
            gain_threshold_relative = yaml_data["conf"]['somatic_cna']['gain_threshold_relative'],
            del_threshold_absolute = yaml_data["conf"]['somatic_cna']['del_threshold_absolute'],
            del_threshold_relative = yaml_data["conf"]['somatic_cna']['del_threshold_relative'],
            threshold_mode = yaml_data["conf"]['somatic_cna']['threshold_mode'],
            transcript_overlap_fraction = yaml_data["conf"]['somatic_cna']['cna_transcript_overlap_pct'] / 100,
            tumor_ploidy = yaml_data['conf']['sample_properties']['tumor_ploidy'] if yaml_data['conf']['sample_properties']['tumor_ploidy'] != 'NA' else None,
            tumor_purity = yaml_data['conf']['sample_properties']['tumor_purity'] if yaml_data['conf']['sample_properties']['tumor_purity'] != 'NA' else None,
            expression_data = expression_data,
            logger = logger)
        if cna_annotation['status'] == 0:
            logger.info('Finished pcgr-annotate-cna-segments')
            # Update YAML with the resolved effective amplification threshold and ploidy provenance
            if cna_annotation['amp_threshold_effective'] is not None:
                yaml_data['conf']['somatic_cna']['amp_threshold_effective'] = cna_annotation['amp_threshold_effective']
            yaml_data['conf']['sample_properties']['tumor_ploidy_source'] = cna_annotation['tumor_ploidy_source']
            if cna_annotation['tumor_ploidy'] is not None:
                yaml_data['conf']['sample_properties']['tumor_ploidy'] = cna_annotation['tumor_ploidy']
        else:
            yaml_data['molecular_data']['fname_cna_gene_tsv'] = "None"
            yaml_data['molecular_data']['fname_cna_segment_tsv'] = "None"
        print('----')
    else:
        logger = getlogger("pcgr-annotate-cna-segments")
        logger.info('PCGR - CNA ANALYSIS SECTION: Annotation of copy number segments - cytobands, overlapping transcripts, and biomarkers - OMITTED (no data available)')
        print('----')

    # PCGR|Fusions - Annotate fusions with biomarkers
    if not input_rna_fusion == 'None':
        logger = getlogger("pcgr-annotate-rna-fusion")
        logger.info('PCGR - RNA FUSION ANALYSIS SECTION: Annotation of RNA fusions with biomarkers')
        fusion_annotation = cna.annotate_fusions(
            input_fusion_fname = input_rna_fusion,
            output_fusion_fname = output_data['rna_fusion'],
            oncokb_input_fname = oncokb_input_fusion_tsv, 
            build = yaml_data['genome_assembly'],
            sample_id = yaml_data['sample_id'],           
            refdata_assembly_dir = input_data['refdata_assembly_dir'],            
            logger = logger)
        if fusion_annotation == 0:
            logger.info('Finished pcgr-annotate-rna-fusion')
        else:
            yaml_data['molecular_data']['fname_rna_fusion_tsv'] = "None"
        print('----')
    else:
        logger = getlogger("pcgr-annotate-rna-fusion")
        logger.info('PCGR - RNA FUSION ANALYSIS SECTION: Annotation of RNA fusions with biomarkers - OMITTED (no data available)')
        print('----')

    # PCGR|OncoKB - Annotate variants, CNAs, and fusions with OncoKB
    if conf_options.get('oncokb', {}).get('api_token'):
        from pcgr import biomarker
        logger = getlogger("pcgr-oncokb-annotation")
        logger.info('PCGR - ONCOKB ANNOTATION SECTION: OncoKB annotation of SNVs/InDels, CNAs, and RNA fusions')

        # Resolve oncotree code based on primary site and user input
        oncotree_code = conf_options['oncokb'].get('oncotree_code')
        primary_site = yaml_data['conf']['sample_properties']['site']

        if primary_site == 'Any':
            if oncotree_code is not None:
                logger.warning(
                    f"OncoTree code '{oncotree_code}' is specified but primary tumor site is 'Any' "
                    "- OncoKB annotation will proceed without tumor type specification")
            else:
                logger.info('Primary tumor site is not set (Any) - OncoKB annotation will proceed without tumor type specification')
            oncotree_code = None
        else:
            if oncotree_code is not None:
                oncotree_code = verify_oncotree_code(
                    oncotree_code,
                    primary_site,
                    input_data["refdata_assembly_dir"],
                    logger)
            ## Fallback strategy: if primary site is specified but no OncoTree code provided, attempt to map 
            ## primary site to OncoTree code with SITE_TO_ONCOTREE mapping
            if oncotree_code is None:
                if primary_site in pcgr_vars.SITE_TO_ONCOTREE:
                    oncotree_code = pcgr_vars.SITE_TO_ONCOTREE[primary_site]
                    logger.info(f'Mapped primary site "{primary_site}" to OncoTree code: {oncotree_code}')
                else:
                    logger.warning(f'Primary site "{primary_site}" not found in SITE_TO_ONCOTREE mapping - skipping OncoKB annotation')

            ## update conf with resolved oncotree code (even if None) for record-keeping and downstream use
            yaml_data['conf']['oncokb']['oncotree_code'] = oncotree_code

        # Map genome assembly to OncoKB reference genome format
        reference_genome = "GRCh38" if yaml_data['genome_assembly'] == "grch38" else "GRCh37"

        # Prepare file paths (only pass files that exist)
        maf_file = oncokb_input_maf if os.path.exists(oncokb_input_maf) else None
        fusion_file = oncokb_input_fusion_tsv if os.path.exists(oncokb_input_fusion_tsv) else None
        cna_file = oncokb_input_cna_tsv if os.path.exists(oncokb_input_cna_tsv) else None

        # Check if any input files exist
        if any([maf_file, fusion_file, cna_file]):
            try:
                result = biomarker.run_oncokb_annotator(
                    sample_name = yaml_data['sample_id'],
                    oncokb_token = conf_options['oncokb']['api_token'],
                    oncotree_code = oncotree_code,
                    build = reference_genome,
                    maf_query_file = maf_file,
                    fusion_query_file = fusion_file,
                    cna_query_file = cna_file,
                    output_dir = output_dir,
                    logger = logger
                )
            except Exception as e:
                logger.warning(f'OncoKB annotation failed: {str(e)}')
        else:
            logger.info('No OncoKB input files found - skipping OncoKB annotation')
    else:
        logger = getlogger("pcgr-oncokb-annotation")
        logger.info('PCGR - ONCOKB ANNOTATION SECTION: OncoKB annotation - OMITTED (no API token provided)')
        print('----')

    # Merge OncoKB SNV/InDel annotations into the final variant TSV (always — adds NA columns when OncoKB was not run)
    if os.path.exists(output_pass_tsv_gz):
        if conf_options.get('oncokb', {}).get('api_token'):
            logger = getlogger("pcgr-oncokb-annotation")
            logger.info('Appending OncoKB SNV/InDel annotations to variant TSV')
        append_oncokb_snv_annotations(
            tsv_gz_fname    = output_pass_tsv_gz,
            oncokb_maf_hgvsp = oncokb_output_maf_hgvsp if os.path.exists(oncokb_output_maf_hgvsp) else None,
            oncokb_maf_hgvsg = oncokb_output_maf_hgvsg if os.path.exists(oncokb_output_maf_hgvsg) else None,
            logger          = logger
        )

    # Refine internal oncogenicity classification using OncoKB annotations (only when *_OKB columns are present)
    if os.path.exists(output_pass_tsv_gz):
        if conf_options.get('oncokb', {}).get('api_token'):
            logger = getlogger("pcgr-oncokb-annotation")
            logger.info('Refining ClinGen/VICC/CGC oncogenicity classification with OncoKB annotations')
        refine_oncogenicity_with_oncokb(
            tsv_gz_fname = output_pass_tsv_gz,
            logger       = logger
        )

    # Annotate TSG/LOH CNA records with somatic and germline two-hit candidates.
    # Runs after append_oncokb_snv_annotations so MUTATION_EFFECT_OKB is present in
    # the SNV TSV, enabling LoF detection beyond LOSS_OF_FUNCTION (e.g. OncoKB-curated missense).
    cna_gene_tsv = yaml_data['molecular_data'].get('fname_cna_gene_tsv', 'None')
    if cna_gene_tsv != 'None' and os.path.exists(cna_gene_tsv):
        snv_df_twohit = None
        if os.path.exists(output_pass_tsv_gz):
            snv_df_twohit = pd.read_csv(output_pass_tsv_gz, sep='\t', na_values='.',
                                        compression='gzip', low_memory=False)
        germline_df = None
        if not input_cpsr_calls == 'None':
            try:
                germline_df = pd.read_csv(input_cpsr_calls, sep='\t', na_values='.')
            except Exception as e:
                logger.warning(f"Could not load germline variant calls from {input_cpsr_calls}: {e} - skipping germline two-hit annotation")
        logger = getlogger("pcgr-twohit-candidates")
        logger.info('PCGR - TWO-HIT SECTION: Annotating TSG/LOH CNA records with somatic and germline two-hit candidates')
        cna_df_twohit = pd.read_csv(cna_gene_tsv, sep='\t', na_values='.', low_memory=False)
        cna_df_twohit = cna.append_twohit_candidates(cna_df_twohit, snv_df_twohit, germline_df, logger)
        cna_df_twohit.to_csv(cna_gene_tsv, sep='\t', header=True, index=False)
        logger.info('Finished annotating two-hit candidates')
        print('----')

    # Merge OncoKB fusion annotations into the fusion TSV (always — adds NA columns when OncoKB was not run)
    fusion_tsv = yaml_data['molecular_data'].get('fname_rna_fusion_tsv', 'None')
    if fusion_tsv != 'None' and os.path.exists(fusion_tsv):
        if conf_options.get('oncokb', {}).get('api_token'):
            logger = getlogger("pcgr-oncokb-annotation")
            logger.info('Appending OncoKB fusion annotations to fusion TSV')
        append_oncokb_fusion_annotations(
            tsv_fname            = fusion_tsv,
            oncokb_fusion_output = oncokb_output_fusion_tsv if os.path.exists(oncokb_output_fusion_tsv) else None,
            logger               = logger
        )

    # Merge OncoKB CNA annotations into the gene-level CNA TSV (always — adds NA columns when OncoKB was not run)
    cna_gene_tsv = yaml_data['molecular_data'].get('fname_cna_gene_tsv', 'None')
    if cna_gene_tsv != 'None' and os.path.exists(cna_gene_tsv):
        if conf_options.get('oncokb', {}).get('api_token'):
            logger = getlogger("pcgr-oncokb-annotation")
            logger.info('Appending OncoKB CNA annotations to gene-level CNA TSV')
        append_oncokb_cna_annotations(
            tsv_fname         = cna_gene_tsv,
            oncokb_cna_output = oncokb_output_cna_tsv if os.path.exists(oncokb_output_cna_tsv) else None,
            logger            = logger
        )

    if conf_options.get('oncokb', {}).get('api_token'):
        logger.info('Finished OncoKB annotation and integration')
    print('----')
    

    # Write YAML file with configuration options and paths to annotated molecular profile datasets
    with open(yaml_fname, "w") as outfile:
        outfile.write(yaml.dump(yaml_data))
    outfile.close()

    # PCGR|Report - Generation of Excel workbooks and integrative HTML reports for molecular data interpretation
    ## SNVs/InDels, CNAs, expression, TMB, MSI, mutational signatures
    if not conf_options['other']['no_reporting']:
        logger = getlogger('pcgr-writer')
        logger.info('PCGR - REPORT GENERATION SECTION: Generation of output files - molecular interpretation report for precision cancer medicine')
        # export PATH to R conda env Rscript
        pcgrr_conda = conf_options['pcgrr_conda']
        quarto_env_vars = utils.quarto_evars_path(pcgrr_conda)
        pcgr_conda = utils.conda_prefix_basename()
        rscript = utils.script_path(pcgrr_conda, 'bin/Rscript')
        pcgrr_script = utils.script_path(pcgr_conda, 'bin/pcgrr.R')
        export_pcgrr = utils.pcgrr_conda_env_export(pcgrr_conda)
        pcgr_report_command = (
            f"{export_pcgrr} && {rscript} {pcgrr_script} {yaml_fname} {quarto_env_vars}")

        if debug:
            print(pcgr_report_command)
        if not input_vcf == 'None' and n_vcf_pass_variants == 0:
            logger.info('There are n = 0 PASS variants in the input VCF - skipping report generation')
        else:
            check_subprocess(logger, pcgr_report_command, debug)
        logger.info("Finished PCGR!")
        print('----')

    # Redact OncoKB API token from the YAML config file to avoid credential exposure 
    # in shared output files
    if os.path.exists(yaml_fname):
        with open(yaml_fname) as fh:
            yaml_data_redact = yaml.safe_load(fh)
        if isinstance(yaml_data_redact.get('conf', {}).get('oncokb'), dict):
            yaml_data_redact['conf']['oncokb']['api_token'] = None
        with open(yaml_fname, "w") as fh:
            fh.write(yaml.dump(yaml_data_redact))

    print()


if __name__ == "__main__":
    cli()
