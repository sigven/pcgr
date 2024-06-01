#!/usr/bin/env python

from pcgr import pcgr_vars, arg_checker, utils, cna
from pcgr.utils import getlogger, check_subprocess, remove_file, random_id_generator
from pcgr.config import populate_config_data, create_config
from pcgr.maf import update_maf_allelic_support
from pcgr.vep import get_vep_command
from pcgr.expression import parse_expression, integrate_variant_expression, correlate_sample_expression
from pcgr.expression import find_expression_outliers, aggregate_tpm_per_cons
from pcgr.variant import clean_annotations, set_allelic_support, append_annotations, calculate_tmb

import re
import argparse
import pandas as pd
import yaml
import os
from glob import glob
from argparse import RawTextHelpFormatter


def cli():

    program_description = (f"Personal Cancer Genome Reporter (PCGR) workflow for clinical translation of "
                           f"tumor omics data (SNVs/InDels, CNA, RNA expression) - version: {pcgr_vars.PCGR_VERSION}")
    program_options = "\n\t--input_vcf <INPUT_VCF>\n\t--vep_dir <VEP_DIR>\n\t--refdata_dir <REFDATA_DIR>\n\t" + \
        "--output_dir <OUTPUT_DIR>\n\t--genome_assembly <GENOME_ASSEMBLY>\n\t--sample_id <SAMPLE_ID>"

    parser = argparse.ArgumentParser(description=program_description,
                                     formatter_class=RawTextHelpFormatter,
                                     usage=f"\n\t%(prog)s -h [options] {program_options} \n\n")
    parser._action_groups.pop()
    required = parser.add_argument_group("Required arguments")
    optional_assay = parser.add_argument_group("Sequencing assay options")
    optional_sample = parser.add_argument_group("Tumor sample options")
    optional_allelic_support = parser.add_argument_group("Allelic support options")
    optional_tumor_only = parser.add_argument_group("Tumor-only filtering options")
    optional_vep = parser.add_argument_group("VEP options")
    optional_tmb_msi = parser.add_argument_group("Tumor mutational burden (TMB) and MSI options")
    optional_signatures = parser.add_argument_group("Mutational signature options")
    optional_cna = parser.add_argument_group("Somatic copy number alteration (CNA) data options")
    optional_rna = parser.add_argument_group("Bulk RNA-seq and RNA fusion data options")
    #optional_germline = parser.add_argument_group("Germline variant options")    
    optional_other = parser.add_argument_group("Other options")


    required.add_argument("--input_vcf", dest="input_vcf", help="VCF input file with somatic variants in tumor sample, SNVs/InDels", required=True)
    required.add_argument("--vep_dir", dest="vep_dir", help="Directory of VEP cache, e.g.  $HOME/.vep", required=True)
    required.add_argument("--refdata_dir", dest="refdata_dir", help="Directory where PCGR reference data bundle was downloaded and unpacked", required=True)
    required.add_argument("--output_dir", dest="output_dir", help="Output directory", required=True)
    required.add_argument("--genome_assembly", dest="genome_assembly", choices=["grch37", "grch38"], help="Human genome assembly build: grch37 or grch38", required=True)
    required.add_argument("--sample_id", dest="sample_id", help="Tumor sample/cancer genome identifier - prefix for output files", required=True)

    optional_assay.add_argument("--assay", dest="assay", default="WES", choices=[ "WGS", "WES","TARGETED"], help="Type of DNA sequencing assay performed for input data (VCF), default: %(default)s")
    optional_assay.add_argument("--effective_target_size_mb", type=float, default=34, dest="effective_target_size_mb", help="Effective target size in Mb (potentially limited by read depth) of sequencing assay (for TMB analysis) (default: %(default)s (WES/WGS))")
    optional_assay.add_argument("--tumor_only", action="store_true", help="Input VCF comes from tumor-only sequencing, calls will be filtered for variants of germline origin, (default: %(default)s)")    
    optional_sample.add_argument("--tumor_site", dest="tsite", type=int, default=0, help="Optional integer code to specify primary tumor type/site of query sample,\nchoose any of the following identifiers:\n" + str(pcgr_vars.tumor_sites) + "\n(default: %(default)s - any tumor type)")
    optional_sample.add_argument("--tumor_purity", type=float, dest="tumor_purity", help="Estimated tumor purity (between 0 and 1) (default: %(default)s)")
    optional_sample.add_argument("--tumor_ploidy", type=float, dest="tumor_ploidy", help="Estimated tumor ploidy (default: %(default)s)")

    
    optional_allelic_support.add_argument("--tumor_dp_tag", dest="tumor_dp_tag", default="_NA_", help="Specify VCF INFO tag for sequencing depth (tumor, must be Type=Integer, default: %(default)s")
    optional_allelic_support.add_argument("--tumor_af_tag", dest="tumor_af_tag", default="_NA_", help="Specify VCF INFO tag for variant allelic fraction (tumor,  must be Type=Float, default: %(default)s")
    optional_allelic_support.add_argument("--control_dp_tag", dest="control_dp_tag", default="_NA_", help="Specify VCF INFO tag for sequencing depth (control, must be Type=Integer, default: %(default)s")
    optional_allelic_support.add_argument("--control_af_tag", dest="control_af_tag", default="_NA_", help="Specify VCF INFO tag for variant allelic fraction (control, must be Type=Float, default: %(default)s")
    optional_allelic_support.add_argument("--call_conf_tag", dest="call_conf_tag", default="_NA_", help="Specify VCF INFO tag for somatic variant call confidence (must be categorical, e.g. Type=String, default: %(default)s")
    optional_allelic_support.add_argument("--tumor_dp_min", type=int, default=0, dest="tumor_dp_min", help="If VCF INFO tag for sequencing depth (tumor) is specified and found, set minimum required depth for inclusion in report (default: %(default)s)")
    optional_allelic_support.add_argument("--tumor_af_min", type=float, default=0, dest="tumor_af_min", help="If VCF INFO tag for variant allelic fraction (tumor) is specified and found, set minimum required AF for inclusion in report (default: %(default)s)")
    optional_allelic_support.add_argument("--control_dp_min", type=int, default=0, dest="control_dp_min", help="If VCF INFO tag for sequencing depth (control) is specified and found, set minimum required depth for inclusion in report (default: %(default)s)")
    optional_allelic_support.add_argument("--control_af_max", type=float, default=1, dest="control_af_max", help="If VCF INFO tag for variant allelic fraction (control) is specified and found, set maximum tolerated AF for inclusion in report (default: %(default)s)")

    optional_tumor_only.add_argument("--pon_vcf", dest="pon_vcf", help="VCF file with germline calls from Panel of Normals (PON) - i.e. blacklisted variants, (default: %(default)s)")
    maf_help_msg = "Exclude variants in tumor (SNVs/InDels, tumor-only mode) with MAF > pct"
    optional_tumor_only.add_argument("--maf_gnomad_nfe", dest="maf_gnomad_nfe", type=float, default=0.002, help=f"{maf_help_msg}, (gnomAD - European (non-Finnish), default: %(default)s)")
    optional_tumor_only.add_argument("--maf_gnomad_asj", dest="maf_gnomad_asj", type=float, default=0.002, help=f"{maf_help_msg}, (gnomAD - Ashkenazi Jewish, default: %(default)s)")
    optional_tumor_only.add_argument("--maf_gnomad_fin", dest="maf_gnomad_fin", type=float, default=0.002, help=f"{maf_help_msg}, (gnomAD - European (Finnish), default: %(default)s)")
    optional_tumor_only.add_argument("--maf_gnomad_oth", dest="maf_gnomad_oth", type=float, default=0.002, help=f"{maf_help_msg}, (gnomAD - Other, default: %(default)s)")
    optional_tumor_only.add_argument("--maf_gnomad_amr", dest="maf_gnomad_amr", type=float, default=0.002, help=f"{maf_help_msg}, (gnomAD - Latino/Admixed American, default: %(default)s)")
    optional_tumor_only.add_argument("--maf_gnomad_afr", dest="maf_gnomad_afr", type=float, default=0.002, help=f"{maf_help_msg}, (gnomAD - African/African-American, default: %(default)s)")
    optional_tumor_only.add_argument("--maf_gnomad_eas", dest="maf_gnomad_eas", type=float, default=0.002, help=f"{maf_help_msg}, (gnomAD - East Asian, default: %(default)s)")
    optional_tumor_only.add_argument("--maf_gnomad_sas", dest="maf_gnomad_sas", type=float, default=0.002, help=f"{maf_help_msg}, (gnomAD - South Asian, default: %(default)s)")
    optional_tumor_only.add_argument("--maf_gnomad_global", dest="maf_gnomad_global", type=float, default=0.002, help=f"{maf_help_msg}, (gnomAD - global population, default: %(default)s)")
    optional_tumor_only.add_argument("--exclude_pon", action="store_true", help="Exclude variants occurring in PoN (Panel of Normals, if provided as VCF (--pon_vcf), default: %(default)s)")
    optional_tumor_only.add_argument("--exclude_likely_hom_germline", action="store_true", help="Exclude likely homozygous germline variants (allelic fraction of 1.0 for alternate allele in tumor - very unlikely somatic event), default: %(default)s)")
    optional_tumor_only.add_argument("--exclude_likely_het_germline", action="store_true", help="Exclude likely heterozygous germline variants (0.4-0.6 allelic fraction, AND presence in dbSNP + gnomAD, AND not existing as somatic record in COSMIC OR TCGA, default: %(default)s)")
    optional_tumor_only.add_argument("--exclude_clinvar_germline", action="store_true", help="Exclude variants found in ClinVar (germline variant origin), defult: %(default)s)")
    optional_tumor_only.add_argument("--exclude_dbsnp_nonsomatic", action="store_true", help="Exclude variants found in dbSNP (except for those present in ClinVar (somatic origin) OR TCGA OR COSMIC), defult: %(default)s)")
    optional_tumor_only.add_argument("--exclude_nonexonic", action="store_true", help="Exclude non-exonic variants, default: %(default)s)")

    optional_vep.add_argument("--vep_n_forks", default=4, type=int, help="Number of forks (VEP option '--fork'), default: %(default)s")
    optional_vep.add_argument("--vep_buffer_size", default=500, type=int, help=f"Variant buffer size (variants read into memory simultaneously, VEP option '--buffer_size')\n- set lower to reduce memory usage, default: %(default)s")
    optional_vep.add_argument("--vep_pick_order", default="mane_select,mane_plus_clinical,canonical,appris,tsl,biotype,ccds,rank,length", help=f"Comma-separated string " + \
        "of ordered transcript/variant properties for selection of primary variant consequence\n(option '--pick_order' in VEP), default: %(default)s")
    optional_vep.add_argument("--vep_no_intergenic", action="store_true", help="Skip intergenic variants during variant annotation (VEP option '--no_intergenic' in VEP), default: %(default)s")
    optional_vep.add_argument("--vep_regulatory", action="store_true", help="Add VEP regulatory annotations (VEP option '--regulatory') or non-coding interpretation, default: %(default)s")
    optional_vep.add_argument("--vep_gencode_basic", action="store_true", help = "Consider basic GENCODE transcript set only with Variant Effect Predictor (VEP) (VEP option '--gencode_basic').")

    optional_tmb_msi.add_argument("--estimate_tmb", action="store_true", help="Estimate tumor mutational burden from the total number of somatic mutations and target region size, default: %(default)s")
    optional_tmb_msi.add_argument("--tmb_display", dest="tmb_display", default="coding_and_silent", choices=["coding_and_silent", "coding_non_silent", "missense_only"], help="Type of TMB measure to show in report, default: %(default)s")
    optional_tmb_msi.add_argument("--tmb_dp_min", dest="tmb_dp_min", default=0, help="If VCF INFO tag for sequencing depth (tumor) is specified and found, set minimum required sequencing depth for TMB calculation: default: %(default)s")
    optional_tmb_msi.add_argument("--tmb_af_min", dest="tmb_af_min", default=0, help="If VCF INFO tag for allelic fraction (tumor) is specified and found, set minimum required allelic fraction for TMB calculation: default: %(default)s")
    optional_tmb_msi.add_argument("--estimate_msi", action="store_true", help="Predict microsatellite instability status from patterns of somatic mutations/indels, default: %(default)s")

    optional_signatures.add_argument("--estimate_signatures", action="store_true", help="Estimate relative contributions of reference mutational signatures in query sample (re-fitting), default: %(default)s")
    optional_signatures.add_argument("--min_mutations_signatures", type=int, default=200, dest="min_mutations_signatures", help="Minimum number of SNVs required for re-fitting of mutational signatures (SBS) (default: %(default)s, minimum n = 100)")
    optional_signatures.add_argument("--all_reference_signatures", action="store_true", help="Use _all_ reference mutational signatures (SBS) during signature re-fitting rather than only those already attributed to the tumor type (default: %(default)s)")
    optional_signatures.add_argument("--include_artefact_signatures", action="store_true", help="Include sequencing artefacts in the collection of reference signatures (default: %(default)s")
    optional_signatures.add_argument("--prevalence_reference_signatures", type=float, default=0.1, help="Minimum tumor-type prevalence (in percent) of reference signatures to be included in refitting procedure (default: %(default)s)")

    optional_cna.add_argument("--input_cna", dest="input_cna", help="Somatic copy number alteration segments (tab-separated values)")   
    optional_cna.add_argument("--n_copy_gain", type=int, default=6, dest="n_copy_gain", help="Minimum number of total copy number for segments considered as gains/amplifications (default: %(default)s)")
    optional_cna.add_argument("--cna_overlap_pct", type=float, default=50, dest="cna_overlap_pct", help="Mean percent overlap between copy number segment and gene transcripts for reporting of gains/losses in tumor suppressor genes/oncogenes, (default: %(default)s)")
    
    
    #optional_rna.add_argument("--input_rna_fusion", dest = "input_rna_fusion", help = "File with RNA fusion transcripts detected in tumor (tab-separated values)")
    optional_rna.add_argument("--input_rna_expression", dest = "input_rna_exp", help = "File with bulk RNA expression counts (TPM) of transcripts in tumor (tab-separated values)")
    optional_rna.add_argument('--expression_sim', action='store_true', help="Compare expression profile of tumor sample to known expression profiles (default: %(default)s)")
    optional_rna.add_argument("--expression_sim_db", dest = "expression_sim_db", default="tcga,depmap,treehouse", help=f"Comma-separated string " + \
        "of databases for used in RNA expression similarity analysis, default: %(default)s") 

   
    #optional_germline.add_argument("--input_germline", dest="input_germline", help="CPSR-classified germline calls (file '<sample_id>.cpsr.<genome_assembly>.classification.tsv.gz')")
    #optional_germline.add_argument("--sample_id_germline", dest="sample_id_germline", help="Sample identifier for germline calls - used for verification of input_germline file")
    
    optional_other.add_argument("--vcf2maf", action="store_true", help="Generate a MAF file for input VCF using https://github.com/mskcc/vcf2maf (default: %(default)s)")
    optional_other.add_argument("--vcfanno_n_proc", default=4, type=int, help="Number of vcfanno processes (option '-p' in vcfanno), default: %(default)s")
    optional_other.add_argument("--ignore_noncoding", action="store_true", help="Ignore non-coding (i.e. non protein-altering) variants in report, default: %(default)s")
    #optional_other.add_argument("--include_trials", action="store_true", help="Include relevant ongoing or future clinical trials, focusing on studies with molecularly targeted interventions")
    optional_other.add_argument("--retained_info_tags", dest="retained_info_tags", default="None", help="Comma-separated string of VCF INFO tags from query VCF that should be kept in PCGR output TSV file")
    optional_other.add_argument("--force_overwrite", action="store_true", help="By default, the script will fail with an error if any output file already exists. You can force the overwrite of existing result files by using this flag, default: %(default)s")
    optional_other.add_argument("--version", action="version", version="%(prog)s " + str(pcgr_vars.PCGR_VERSION))
    optional_other.add_argument("--no_reporting", action="store_true", help="Run functional variant annotation on VCF through VEP/vcfanno, omit other analyses (i.e. Tier assignment/MSI/TMB/Signatures etc. and report generation (STEP 4), default: %(default)s")
    optional_other.add_argument("--debug", action="store_true", help="Print full commands to log")
    optional_other.add_argument("--pcgrr_conda", default="pcgrr", help="pcgrr conda env name (default: %(default)s)")

    
    
    # Parse required and optional arguments
    args = parser.parse_args()
    arg_dict = vars(args)

    # verify parsed arguments
    arg_dict_verified = arg_checker.verify_args(arg_dict)
    
    # create config options
    conf_options = create_config(arg_dict_verified, workflow = "PCGR")
    
    # Verify existence of input files, define and check existence of output files
    input_data = arg_checker.verify_input_files(arg_dict_verified)
    output_data = arg_checker.define_output_files(arg_dict_verified)
    
    # Run PCGR workflow
    run_pcgr(input_data, output_data, conf_options)

def run_pcgr(input_data, output_data,conf_options):
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
    input_germline_cpsr = 'None'
    pon_vcf = 'None'
    pon_annotation = 0
    variant_set = pd.DataFrame
    expression_data = None


    # Update input data variables with data provided by user
    if input_data['vcf_basename'] != 'NA':
        input_vcf = os.path.join(input_data['vcf_dir'], input_data['vcf_basename'])
    if input_data['cna_basename'] != 'NA':
        input_cna = os.path.join(input_data['cna_dir'], input_data['cna_basename'])
    if input_data['rna_fusion_basename'] != 'NA':
        input_rna_fusion = os.path.join(input_data['rna_fusion_dir'], input_data['rna_fusion_basename'])
    if input_data['rna_expression_basename'] != 'NA':
        input_rna_expression = os.path.join(input_data['rna_expression_dir'], input_data['rna_expression_basename'])
    #if input_data['germline_basename'] != 'NA':
    #    input_germline_cpsr = os.path.join(input_data['germline_dir'], input_data['germline_basename'])
    if input_data['pon_vcf_basename'] != 'NA':
        pon_vcf = os.path.join(input_data['pon_vcf_dir'], input_data['pon_vcf_basename'])

    if not input_vcf == 'None':
        output_dir = output_data['dir']
        output_prefix = output_data['prefix']        
        check_subprocess(logger, f'mkdir -p {output_dir}', debug)

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
        yaml_fname =                      f'{output_prefix}.conf.yaml'
        tmb_fname =                       f'{output_prefix}.tmb.tsv'

        # PCGR|validate_input - verify that VCF and CNA segment file is of appropriate format
        logger = getlogger("pcgr-validate-input-arguments")
        print('')
        logger.info("PCGR - STEP 0: Validate input data and options")

        vcf_validate_command = (
                f'pcgr_validate_input.py '
                f'{input_data["refdata_assembly_dir"]} '
                f'{input_vcf} '
                f'{input_vcf_validated_uncompr} '
                f'{input_cna} '
                f'{input_rna_fusion} '
                f'{input_rna_expression} '
                f'{pon_vcf} '
                f'{conf_options["assay_properties"]["vcf_tumor_only"]} '
                f'{conf_options["sample_id"]} '
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
        logger.info(f'Sample name: {conf_options["sample_id"]}')
        if conf_options['sample_properties']['site'] == 'Any':
            logger.info('Tumor type: Cancer_NOS (Any tumortype)')
        else:
            logger.info(f'Tumor type: {conf_options["sample_properties"]["site"]}')
        logger.info(f'Sequencing assay - type: {conf_options["assay_properties"]["type"]}')
        logger.info(f'Sequencing assay - mode: {assay_mode}')
        logger.info((
            f'Sequencing assay - effective (coding) target size: '
            f'{conf_options["assay_properties"]["effective_target_size_mb"]}Mb'))
        logger.info((
            f'Variant filtering settings - minimum sequencing depth tumor: '
            f'{conf_options["somatic_snv"]["allelic_support"]["tumor_dp_min"]}'))
        logger.info((
            f'Variant filtering settings - minimum allelic fraction tumor: '
            f'{conf_options["somatic_snv"]["allelic_support"]["tumor_af_min"]}'))
        logger.info((
            f'Variant filtering settings - minimum sequencing depth control: '
            f'{conf_options["somatic_snv"]["allelic_support"]["control_dp_min"]}'))
        logger.info((
            f'Variant filtering settings - maximum allelic fraction control: '
            f'{conf_options["somatic_snv"]["allelic_support"]["control_af_max"]}'))
        logger.info(f'Genome assembly: {conf_options["genome_assembly"]}')
        logger.info(f'Mutational signature estimation: {msig_estimation_set}')
        logger.info(f'MSI classification: {msi_prediction_set}')
        logger.info(f'Mutational burden estimation: {tmb_estimation_set}')
        logger.info(f'RNA expression similarity analysis: {rnaseq_sim_analysis_set}')
        #logger.info(f'Include molecularly targeted clinical trials (beta): {clinical_trials_set}')
        
        # PCGR|Generate YAML file - containing configuration options and paths to annotated molecular profile datasets
        # - VCF/TSV files (SNVs/InDels)
        # - TSV files (copy number aberrations)
        # - TSV files (TMB)
        # - TSV files (RNA expression)
        # - TSV files (RNA fusion) - COMING
        logger = getlogger('pcgr-write-yaml')

        # update conf_options with paths to output files
        conf_options['output_dir'] = output_dir
        conf_options['output_prefix'] = output_prefix
        conf_options['molecular_data']['fname_mut_vcf'] = output_vcf
        conf_options['molecular_data']['fname_mut_tsv'] = output_pass_tsv_gz
        if conf_options['somatic_snv']['tmb']['run'] == 1:
            conf_options['molecular_data']['fname_tmb_tsv'] = tmb_fname
        if not input_cna == 'None': 
            conf_options['molecular_data']['fname_cna_tsv'] = output_data['cna']
        if not input_rna_expression == 'None':
            conf_options['molecular_data']['fname_expression_tsv'] = output_data['expression']
            conf_options['molecular_data']['fname_expression_outliers_tsv'] = output_data['expression_outliers']
            if conf_options['expression']['similarity_analysis'] == 1:
                conf_options['molecular_data']['fname_expression_similarity_tsv'] = output_data['expression_similarity']
                
        # make YAML file
        yaml_data = populate_config_data(conf_options, input_data["refdata_assembly_dir"], 
                                         workflow = "PCGR", logger = logger)
        genome_assembly = yaml_data['genome_assembly']
        
        vep_command = get_vep_command(file_paths = input_data, 
                                conf_options = yaml_data, 
                                input_vcf = input_vcf_validated, 
                                output_vcf = vep_vcf)

        # PCGR|VEP - run consequence annotation with Variant Effect Predictor
        print('----')
        logger = getlogger('pcgr-vep')
        logger.info(f'PCGR - STEP 1: Basic variant annotation with Variant Effect Predictor {pcgr_vars.VEP_VERSION}' + \
                    f', GENCODE release {pcgr_vars.GENCODE_VERSION[genome_assembly]}, genome assembly {yaml_data["genome_assembly"]}')
        logger.info(f'VEP configuration - one primary consequence block pr. alternative gene allele (--flag_pick_allele_gene)')
        logger.info(f'VEP configuration - transcript pick order: {yaml_data["conf"]["vep"]["vep_pick_order"]}')
        logger.info(f'VEP configuration - transcript pick order: See more at https://www.ensembl.org/info/docs/tools/vep/script/vep_other.html#pick_options')
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
            logger = getlogger("pcgr-vcf2maf")
            
            vcf2maf_command_main = "vcf2maf.pl"
            update_allelic_support = False
            if retained_info_tags != "":
                vcf2maf_command_main = f"vcf2maf.pl --retain-info {retained_info_tags} "
                update_allelic_support = True
            
            logger.info('Converting VEP-annotated VCF to MAF with https://github.com/mskcc/vcf2maf')
            vcf2maf_command = (
                    f'{vcf2maf_command_main} --inhibit-vep --input-vcf {vep_vcf} '
                    f'--tumor-id {yaml_data["sample_id"]} --output-maf {output_tmp_maf} '
                    f'--ref-fasta {vep_command["fasta_assembly"]} '
                    f'--ncbi-build {NCBI_BUILD_MAF} > {output_vcf2maf_log} 2>&1'
                    )
            check_subprocess(logger, vcf2maf_command, debug)
            if not debug:                
                remove_file(output_vcf2maf_log)
            
            ## add information on allelic support in MAF file 
            ## (n_depth, n_ref_count, n_alt_count, t_depth, t_ref_count, t_alt_count)   
            update_maf_allelic_support(
                maf_tmp_fname = output_tmp_maf,
                maf_fname = output_maf,
                allelic_support_tags = yaml_data["conf"]['somatic_snv']['allelic_support'],
                logger = logger,
                update_allelic_support = update_allelic_support
            )
            
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
                f"Annotation sources II (vcfanno): RepeatMasker, SimpleRepeats, WindowMaskerSDust, gnomAD non-cancer subset"
        logger.info("PCGR - STEP 2: Variant annotation for cancer precision medicine with pcgr-vcfanno")
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
                f'{yaml_data["conf"]["sample_properties"]["site2"]} {yaml_data["conf"]["vep"]["vep_pick_order"]} '
                f'{input_data["refdata_assembly_dir"]} --compress_output_vcf '
                f'{"--debug" if debug else ""}'
                )
        summarise_db_src_msg1 = \
                f"Annotation sources: cancerhotspots.org, CIViC, Cancer Biomarkers database (CGI), Cancer Gene Census (CGC)"
        summarise_db_src_msg2 = \
                f"Annotation sources: Network of Cancer Genes (NCG), CancerMine, IntOGen, TCGA driver genes"

        logger.info("PCGR - STEP 3: Variant and cancer gene annotations with pcgr-summarise")
        logger.info(summarise_db_src_msg1)
        logger.info(summarise_db_src_msg2)
        
        logger.info('Variant oncogenicity classification according to ClinGen/VICC recommendations (Horak et al., Genet Med, 2022)')
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
        variant_set = set_allelic_support(variant_set, allelic_support_tags = yaml_data["conf"]['somatic_snv']['allelic_support'])
        ## Clean annotations (formatting etc.)
        variant_set = clean_annotations(variant_set, yaml_data, logger = logger)        
        
        ## Check if AD/DP properties could be detected/pulled from VCFs
        for c in ['DP_TUMOR', 'DP_CONTROL', 'VAF_TUMOR', 'VAF_CONTROL']:
            if c in variant_set.columns:
                if len(variant_set[variant_set[c].isnull()]) == 0:
                    var = str(c).lower() + '_detected'
                    if var in yaml_data['conf']['sample_properties'].keys():
                        yaml_data['conf']['sample_properties'][var] = 1
               
        # PCGR|expression1 - 
        ## If gene expression data is provided as input; verify identifiers, annotate, and merge with somatic variant set 
        if input_rna_expression != 'None':
            logger = getlogger('pcgr-gene-expression')
            ## Parse expression data, verify identifiers
            expression_data = parse_expression(
                input_rna_expression, yaml_data["sample_id"], 
                input_data["refdata_assembly_dir"], logger = logger)
            ## Write transcript-level expression data to TSV
            if 'transcript' in expression_data.keys():
                if not expression_data['transcript'] is None:
                    expression_data['transcript'].fillna('.').to_csv(
                        yaml_data['molecular_data']['fname_expression_tsv'], sep = "\t", 
                        compression = "gzip", index = False)
                    
                    variant_set = aggregate_tpm_per_cons(variant_set, expression_data, logger = logger)        
                    #exp_to_cons.fillna('.').to_csv(
                    #    yaml_data['molecular_data']['fname_csq_expression_tsv'], sep = "\t", 
                    #    compression = "gzip", index = False)
            else:
                if 'gene' in expression_data.keys():
                    if not expression_data['gene'] is None:
                        expression_data['gene'].fillna('.').to_csv(
                            yaml_data['molecular_data']['fname_expression_tsv'], sep = "\t", 
                            compression = "gzip", index = False)
            
            ## Merge expression data with somatic SNV/InDel variant set
            variant_set = integrate_variant_expression(
                variant_set, expression_data, logger = logger)
            
            
        ## Write somatic SNV/InDel variant set to TSV
        variant_set.fillna('.').to_csv(output_pass_tsv_gz, sep = "\t", compression = "gzip", index = False)
        if not debug:
            remove_file(output_pass_vcf2tsv_gz)
        
        ## Reduce output file for WGS when the number of variants exceeds MAX_VARIANTS_FOR_REPORT
        if yaml_data["conf"]['assay_properties']['type'] == 'WGS' or yaml_data["conf"]['assay_properties']['type'] == 'WES':
            # check that output file exist
            if os.path.exists(output_pass_tsv_gz):
                # get number of rows/variants annotated, using pandas
                var_data = pd.read_csv(output_pass_tsv_gz, sep = '\t', low_memory = False, na_values='.')
                num_variants_raw = len(var_data)
                if num_variants_raw > pcgr_vars.MAX_VARIANTS_FOR_REPORT:
                    logger.info(f'Number of raw variants in input VCF ({num_variants_raw}) exceeds ' + \
                                f'{pcgr_vars.MAX_VARIANTS_FOR_REPORT} - intergenic/intronic variants will be excluded prior to reporting')

                    # Exclude intronic and intergenic variants prior to analysis with pcgrr (reporting and further analysis)
                    var_data_filtered = var_data[~var_data.CONSEQUENCE.str.contains('^intron') & 
                                                 ~var_data.CONSEQUENCE.str.contains('^intergenic')]
                    num_variants_excluded1 = num_variants_raw - len(var_data_filtered)
                    logger.info(f'Number of intergenic/intronic variants excluded: {num_variants_excluded1}')

                    # Exclude upstream_gene/downstream_gene variants if size of filtered variant set is still above MAX_VARIANTS_FOR_REPORT
                    # TODO: in this case, the TMB calculation will be an underestimate (but still likely huge)
                    var_data_filtered_final = var_data_filtered
                    if len(var_data_filtered) > pcgr_vars.MAX_VARIANTS_FOR_REPORT:
                        var_data_filtered_final = \
                            var_data_filtered[~var_data_filtered.CONSEQUENCE.str.contains('^upstream_gene') & 
                                              ~var_data_filtered.CONSEQUENCE.str.contains('^downstream_gene')]
                        num_variants_excluded2 = len(var_data_filtered) - len(var_data_filtered_final)
                        logger.info(f'Number of upstream_gene/downstream_gene variants excluded: {num_variants_excluded2}')
                        
                    # rename original vcf2tsv (gzipped) to 'raw' filename
                    rename_output_tsv = f'mv {output_pass_tsv_gz} {output_pass_raw_tsv_gz}'
                    check_subprocess(logger, rename_output_tsv, debug)
                    
                    var_data_filtered_final.fillna('.').to_csv(output_pass_tsv_gz, sep = "\t", compression = "gzip", index = False)
                    logger.info(f'Number of variants in final output TSV: {len(var_data_filtered_final)}')
        
        # PCGR|TMB - calculate TMB in multiple ways - potentially also considering subset of variants (depth and allele frequency filtered)
        if yaml_data['conf']['somatic_snv']['tmb']['run'] == 1:
            logger_tmb = getlogger('pcgr-calculate-tmb')
            calculate_tmb(
                variant_set = variant_set,
                tumor_dp_min = int(yaml_data['conf']['somatic_snv']['tmb']['tmb_dp_min']),
                tumor_af_min = float(yaml_data['conf']['somatic_snv']['tmb']['tmb_af_min']),
                target_size_mb = float(yaml_data['conf']['assay_properties']['effective_target_size_mb']),
                sample_id = yaml_data['sample_id'],
                tmb_fname = tmb_fname,
                logger = logger_tmb)
        
        logger = getlogger('pcgr-summarise')
        logger.info('Finished pcgr-summarise')
        print('----')
    
    # PCGR|Expression2 - Gene expression (bulk RNA-seq) analyses
    if not expression_data is None:
        logger = getlogger("pcgr-expression-analysis")
        logger.info('PCGR - STEP 4: Gene expression analysis')
         
        ## Find expression outliers by comparing with reference datasets
        ## - Preliminary version only compares with TCGA datasets, and only selects specific cohorts (pcgr_vars.SITE_TO_DISEASE)
        logger.info('Identification of expression outliers through comparison with reference datasets')
        expression_outliers = find_expression_outliers(
            expression_data,
            yaml_data,
            input_data["refdata_assembly_dir"],
            logger = logger
        )
        if not expression_outliers.empty:
            expression_outliers.to_csv(
                yaml_data['molecular_data']['fname_expression_outliers_tsv'],
                sep = "\t", index = False)
        else:
            yaml_data['molecular_data']['fname_expression_outliers_tsv'] = 'None'
        
        ## Correlate sample gene expression profile with reference samples (expression-based similarity)
        if yaml_data['conf']['expression']['similarity_analysis'] == 1:
            exp_similarity_results = correlate_sample_expression(
                expression_data,
                yaml_data,
                input_data["refdata_assembly_dir"],
                protein_coding_only = True,
                logger = logger)
            
            ## Aggregate similarity results across sources (TCGA, TreeHouse, Depmap) into single table
            exp_similarity_results_all = pd.DataFrame()
            for source in exp_similarity_results.keys():
                if not exp_similarity_results[source].empty:
                    exp_similarity_results_all = \
                        pd.concat([exp_similarity_results_all, 
                                   exp_similarity_results[source]], axis = 0)
            exp_similarity_results_all.fillna('.').to_csv(
                yaml_data['molecular_data']['fname_expression_similarity_tsv'],
                sep = "\t", index = False)
        print('----')
    else:
        logger = getlogger("pcgr-expression-analysis")
        logger.info('PCGR - STEP 4: Gene expression analysis - OMITTED (no data available)')
        print('----')
        
    # PCGR|CNA - Annotate allele-specific copy number segments with cytobands, overlapping transcripts, and biomarkers
    if not input_cna == 'None':
        logger = getlogger("pcgr-annotate-cna-segments")
        logger.info('PCGR - STEP 5: Annotation of copy number segments - cytobands, overlapping transcripts, and biomarkers')
        cna_annotation = cna.annotate_cna_segments(
            output_fname = output_data['cna'], 
            output_dir = output_data['dir'], 
            cna_segment_file = input_cna,
            build = yaml_data['genome_assembly'],
            sample_id = yaml_data['sample_id'],
            refdata_assembly_dir = input_data['refdata_assembly_dir'], 
            n_copy_amplifications = yaml_data["conf"]['somatic_cna']['n_copy_gain'],
            overlap_fraction = 0.5,
            expression_data = expression_data,
            logger = logger)
        if cna_annotation == 0:
            logger.info('Finished pcgr-annotate-cna-segments')
        print('----')
    else:
        logger = getlogger("pcgr-annotate-cna-segments")
        logger.info('PCGR - STEP 5: Annotation of copy number segments - cytobands, overlapping transcripts, and biomarkers - OMITTED (no data available)')
        print('----')
    
    # Write YAML file with configuration options and paths to annotated molecular profile datasets
    with open(yaml_fname, "w") as outfile:
        outfile.write(yaml.dump(yaml_data))
    outfile.close()
            
    # PCGR|Report - Generation of Excel workbooks and integrative HTML reports for molecular data interpretation 
    ## SNVs/InDels, CNAs, expression, TMB, MSI, mutational signatures
    if not conf_options['other']['no_reporting'] and not input_vcf == 'None':
        logger = getlogger('pcgr-writer')
        logger.info('PCGR - STEP 6: Generation of output files - molecular interpretation report for precision cancer medicine')
        # export PATH to R conda env Rscript
        pcgrr_conda = conf_options['pcgrr_conda']
        quarto_env_vars = utils.quarto_evars_path(pcgrr_conda)
        pcgr_conda = utils.conda_prefix_basename()
        rscript = utils.script_path(pcgrr_conda, 'bin/Rscript')
        pcgrr_script = utils.script_path(pcgr_conda, 'bin/pcgrr.R')
        pcgr_report_command = (
                 f"{rscript} {pcgrr_script} {yaml_fname} {quarto_env_vars}")

        if debug:
            print(pcgr_report_command)
        check_subprocess(logger, pcgr_report_command, debug)
        logger.info("Finished PCGR!")
        print('----')

    print()


if __name__ == "__main__":
    cli()
