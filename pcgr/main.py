#!/usr/bin/env python

from pcgr import pcgr_vars, arg_checker, config, utils
from pcgr.utils import getlogger, check_subprocess
import re
import argparse
import pandas
import os
import sys
import getpass
import platform
from glob import glob
from argparse import RawTextHelpFormatter


def cli():

    program_description = (f'Personal Cancer Genome Reporter (PCGR) workflow for clinical interpretation of '
                           f'somatic nucleotide variants and copy number aberration segments')
    program_options = "\n\t--input_vcf <INPUT_VCF>\n\t--pcgr_dir <PCGR_DIR>\n\t--output_dir <OUTPUT_DIR>\n\t--genome_assembly" + \
    " <GENOME_ASSEMBLY>\n\t--sample_id <SAMPLE_ID>"

    parser = argparse.ArgumentParser(description=program_description,
                                     formatter_class=RawTextHelpFormatter,
                                     usage=f'\n\t%(prog)s -h [options] {program_options} \n\n')
    parser._action_groups.pop()
    required = parser.add_argument_group("Required arguments")
    # optional_rna = parser.add_argument_group("Bulk RNA-seq and RNA fusion options")
    optional_vcfanno = parser.add_argument_group("vcfanno options")
    optional_vep = parser.add_argument_group("VEP options")
    optional_tmb_msi = parser.add_argument_group("Tumor mutational burden (TMB) and MSI options")
    optional_signatures = parser.add_argument_group("Mutational signature options")
    optional_tumor_only = parser.add_argument_group("Tumor-only options")
    optional_allelic_support = parser.add_argument_group("Allelic support options")
    optional_other = parser.add_argument_group("Other options")

    # optional_rna.add_argument("--input_rna_fusion", dest = "rna_fusion_tumor", help = "File with RNA fusion transcripts detected in tumor (tab-separated values)")
    # optional_rna.add_argument("--input_rna_expression", dest = "rna_exp_tumor", help = "File with RNA expression levels (bulk) of genes in tumor (tab-separated values)")

    optional_other.add_argument("--input_cna", dest="input_cna", help="Somatic copy number alteration segments (tab-separated values)")
    optional_other.add_argument("--logr_gain", type=float, default=0.8, dest="logr_gain", help="Log ratio-threshold (minimum) for segments containing copy number gains/amplifications (default: %(default)s)")
    optional_other.add_argument("--logr_homdel", type=float, default=-0.8, dest="logr_homdel", help="Log ratio-threshold (maximum) for segments containing homozygous deletions (default: %(default)s)")
    optional_other.add_argument("--cna_overlap_pct", type=float, default=50, dest="cna_overlap_pct", help="Mean percent overlap between copy number segment and gene transcripts for reporting of gains/losses in tumor suppressor genes/oncogenes, (default: %(default)s)")
    optional_other.add_argument("--tumor_site", dest="tsite", type=int, default=0, help="Optional integer code to specify primary tumor type/site of query sample,\n choose any of the following identifiers:\n" + str(pcgr_vars.tumor_sites) + "\n(default: %(default)s - any tumor type)")
    optional_other.add_argument("--tumor_purity", type=float, dest="tumor_purity", help="Estimated tumor purity (between 0 and 1, (default: %(default)s)")
    optional_other.add_argument("--tumor_ploidy", type=float, dest="tumor_ploidy", help="Estimated tumor ploidy (default: %(default)s)")

    optional_allelic_support.add_argument("--tumor_dp_tag", dest="tumor_dp_tag", default="_NA_", help="Specify VCF INFO tag for sequencing depth (tumor, must be Type=Integer, default: %(default)s")
    optional_allelic_support.add_argument("--tumor_af_tag", dest="tumor_af_tag", default="_NA_", help="Specify VCF INFO tag for variant allelic fraction (tumor,  must be Type=Float, default: %(default)s")
    optional_allelic_support.add_argument("--control_dp_tag", dest="control_dp_tag", default="_NA_", help="Specify VCF INFO tag for sequencing depth (control, must be Type=Integer, default: %(default)s")
    optional_allelic_support.add_argument("--control_af_tag", dest="control_af_tag", default="_NA_", help="Specify VCF INFO tag for variant allelic fraction (control, must be Type=Float, default: %(default)s")
    optional_allelic_support.add_argument("--call_conf_tag", dest="call_conf_tag", default="_NA_", help="Specify VCF INFO tag for somatic variant call confidence (must be categorical, e.g. Type=String, default: %(default)s")
    optional_allelic_support.add_argument("--tumor_dp_min", type=int, default=0, dest="tumor_dp_min", help="If VCF INFO tag for sequencing depth (tumor) is specified and found, set minimum required depth for inclusion in report (default: %(default)s)")
    optional_allelic_support.add_argument("--tumor_af_min", type=float, default=0, dest="tumor_af_min", help="If VCF INFO tag for variant allelic fraction (tumor) is specified and found, set minimum required AF for inclusion in report (default: %(default)s)")
    optional_allelic_support.add_argument("--control_dp_min", type=int, default=0, dest="control_dp_min", help="If VCF INFO tag for sequencing depth (control) is specified and found, set minimum required depth for inclusion in report (default: %(default)s)")
    optional_allelic_support.add_argument("--control_af_max", type=float, default=1, dest="control_af_max", help="If VCF INFO tag for variant allelic fraction (control) is specified and found, set maximum tolerated AF for inclusion in report (default: %(default)s)")

    optional_tmb_msi.add_argument("--target_size_mb", type=float, default=34, dest="target_size_mb", help="For mutational burden analysis - approximate protein-coding target size in Mb of sequencing assay (default: %(default)s (WES/WGS))")
    optional_tmb_msi.add_argument("--estimate_tmb", action="store_true", help="Estimate tumor mutational burden from the total number of somatic mutations and target region size, default: %(default)s")
    optional_tmb_msi.add_argument("--estimate_msi_status", action="store_true", help="Predict microsatellite instability status from patterns of somatic mutations/indels, default: %(default)s")
    optional_tmb_msi.add_argument("--tmb_algorithm", dest="tmb_algorithm", default="all_coding", choices=[ "all_coding", "nonsyn"], help="Method for calculation of TMB, all coding variants (Chalmers et al., Genome Medicine, 2017), or non-synonymous variants only, default: %(default)s")

    optional_signatures.add_argument("--estimate_signatures", action="store_true", help="Estimate relative contributions of reference mutational signatures in query sample and detect potential kataegis events, default: %(default)s")
    optional_signatures.add_argument("--min_mutations_signatures", type=int, default=200, dest="min_mutations_signatures", help="Minimum number of SNVs required for reconstruction of mutational signatures (SBS) by MutationalPatterns (default: %(default)s, minimum n = 100)")
    optional_signatures.add_argument("--all_reference_signatures", action="store_true", help="Use all reference mutational signatures (SBS, n = 67) in signature reconstruction rather than only those already attributed to the tumor type (default: %(default)s)")
    optional_signatures.add_argument('--include_artefact_signatures', action="store_true", help="Include sequencing artefacts in the collection of reference signatures (default: %(default)s")
    optional_signatures.add_argument('--prevalence_reference_signatures', type=int, default=5, choices=[1,2,5,10,15,20], help="Minimum tumor-type prevalence (in percent) of reference signatures to be included in refitting procedure (default: %(default)s)")

    optional_other.add_argument("--cpsr_report", dest="cpsr_report", help="CPSR report file (Gzipped JSON - file ending with 'cpsr.<genome_assembly>.json.gz' -  germline report of patient's blood/control sample")
    optional_other.add_argument("--vcf2maf", action="store_true", help="Generate a MAF file for input VCF using https://github.com/mskcc/vcf2maf (default: %(default)s)")
    optional_other.add_argument("--show_noncoding", action="store_true", help="List non-coding (i.e. non protein-altering) variants in report, default: %(default)s")
    optional_other.add_argument("--assay", dest="assay", choices=["WES", "WGS", "TARGETED"], default="WES", help="Type of DNA sequencing assay performed for input data (VCF) default: %(default)s")
    optional_other.add_argument("--include_trials", action="store_true", help="(Beta) Include relevant ongoing or future clinical trials, focusing on studies with molecularly targeted interventions")
    optional_other.add_argument("--preserved_info_tags", dest="preserved_info_tags", default="None", help="Comma-separated string of VCF INFO tags from query VCF that should be kept in PCGR output TSV file")
    optional_other.add_argument("--report_theme", choices=["default", "cerulean", "journal", "flatly", "readable", "spacelab", "united", "cosmo", "lumen", "paper", "sandstone", "simplex", "yeti"], help="Visual report theme (rmarkdown)", default='default')
    optional_other.add_argument('--report_nonfloating_toc', action='store_true', help='Do not float the table of contents (TOC) in output report (rmarkdown), default: %(default)s')
    optional_other.add_argument("--force_overwrite", action="store_true", help="By default, the script will fail with an error if any output file already exists. You can force the overwrite of existing result files by using this flag, default: %(default)s")
    optional_other.add_argument("--version", action="version", version="%(prog)s " + str(pcgr_vars.PCGR_VERSION))
    optional_other.add_argument("--basic", action="store_true", help="Run functional variant annotation on VCF through VEP/vcfanno, omit other analyses (i.e. Tier assignment/MSI/TMB/Signatures etc. and report generation (STEP 4), default: %(default)s")
    optional_other.add_argument("--no_vcf_validate", action="store_true", help="Skip validation of input VCF with Ensembl's vcf-validator, default: %(default)s")
    optional_other.add_argument("--debug", action="store_true", help="Print full commands to log")

    optional_vcfanno.add_argument("--vcfanno_n_proc", default=4, type=int, help="Number of vcfanno processes (option '-p' in vcfanno), default: %(default)s")
    optional_vep.add_argument("--vep_n_forks", default=4, type=int, help="Number of forks (option '--fork' in VEP), default: %(default)s")
    optional_vep.add_argument("--vep_buffer_size", default=100, type=int, help=f"Variant buffer size (variants read into memory simultaneously, option '--buffer_size' in VEP)\n- set lower to reduce memory usage, default: %(default)s")
    optional_vep.add_argument("--vep_pick_order", default="canonical,appris,biotype,ccds,rank,tsl,length,mane", help=f"Comma-separated string of ordered transcript/variant properties for selection of primary variant consequence\n(option '--pick_order' in VEP), default: %(default)s")
    optional_vep.add_argument("--vep_no_intergenic", action="store_true", help="Skip intergenic variants during processing (option '--no_intergenic' in VEP), default: %(default)s")
    optional_vep.add_argument("--vep_regulatory", action="store_true", help="Add VEP regulatory annotations (option '--regulatory' )or non-coding interpretation, default: %(default)s")
    optional_vep.add_argument('--vep_gencode_all', action='store_true', help = "Consider all GENCODE transcripts with Variant Effect Predictor (VEP) (option '--gencode_basic' in VEP is used by default in PCGR).")


    optional_tumor_only.add_argument("--tumor_only", action="store_true", help="Input VCF comes from tumor-only sequencing, calls will be filtered for variants of germline origin, (default: %(default)s)")
    optional_tumor_only.add_argument("--cell_line", action="store_true", help="Input VCF comes from tumor cell line sequencing (requires --tumor_only), calls will be filtered for variants of germline origin, (default: %(default)s)")
    optional_tumor_only.add_argument("--pon_vcf", dest="pon_vcf", help="VCF file with germline calls from Panel of Normals (PON) - i.e. blacklisted variants, (default: %(default)s)")
    maf_help_msg = "Exclude variants in tumor (SNVs/InDels, tumor-only mode) with MAF > pct"
    optional_tumor_only.add_argument("--maf_onekg_eur", dest="maf_onekg_eur", type=float, default=0.002, help=f"{maf_help_msg}, (1000 Genomes Project - European pop, default: %(default)s)")
    optional_tumor_only.add_argument("--maf_onekg_amr", dest="maf_onekg_amr", type=float, default=0.002, help=f"{maf_help_msg}, (1000 Genomes Project - Ad Mixed American pop, default: %(default)s)")
    optional_tumor_only.add_argument("--maf_onekg_afr", dest="maf_onekg_afr", type=float, default=0.002, help=f"{maf_help_msg}, (1000 Genomes Project - African pop, default: %(default)s)")
    optional_tumor_only.add_argument("--maf_onekg_eas", dest="maf_onekg_eas", type=float, default=0.002, help=f"{maf_help_msg}, (1000 Genomes Project - East Asian pop, default: %(default)s)")
    optional_tumor_only.add_argument("--maf_onekg_sas", dest="maf_onekg_sas", type=float, default=0.002, help=f"{maf_help_msg}, (1000 Genomes Project - South Asian pop, default: %(default)s)")
    optional_tumor_only.add_argument("--maf_onekg_global", dest="maf_onekg_global", type=float, default=0.002, help=f"{maf_help_msg}, (1000 Genomes Project - global pop, default: %(default)s)")
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
    optional_tumor_only.add_argument("--exclude_likely_hom_germline", action="store_true", help="Exclude likely homozygous germline variants (100 pct allelic fraction for alternate allele in tumor, very unlikely somatic event, default: %(default)s)")
    optional_tumor_only.add_argument("--exclude_likely_het_germline", action="store_true", help="Exclude likely heterozygous germline variants (40-60 pct allelic fraction, AND presence in dbSNP + gnomAD, AND not existing as somatic event in COSMIC/TCGA, default: %(default)s)")
    optional_tumor_only.add_argument("--exclude_dbsnp_nonsomatic", action="store_true", help="Exclude variants found in dbSNP (only those that are NOT found in ClinVar(somatic origin)/DoCM/TCGA/COSMIC, defult: %(default)s)")
    optional_tumor_only.add_argument("--exclude_nonexonic", action="store_true", help="Exclude non-exonic variants, default: %(default)s)")

    required.add_argument("--input_vcf", dest="input_vcf", help="VCF input file with somatic variants in tumor sample, SNVs/InDels", required=True)
    required.add_argument("--pcgr_dir", dest="pcgr_dir", help="PCGR base directory with accompanying data directory, e.g. ~/pcgr-" + str(pcgr_vars.PCGR_VERSION), required=True)
    required.add_argument("--output_dir", dest="output_dir", help="Output directory", required=True)
    required.add_argument("--genome_assembly", dest="genome_assembly", choices=["grch37", "grch38"], help="Human genome assembly build: grch37 or grch38", required=True)
    required.add_argument("--sample_id", dest="sample_id", help="Tumor sample/cancer genome identifier - prefix for output files", required=True)

    # Parse required and optional arguments
    args = parser.parse_args()
    arg_dict = vars(args)

    # check parsed arguments
    arg_checker.check_args(arg_dict)
    # create config options
    config_options = config.create_config(arg_dict)
    # Verify existence of input files
    pcgr_paths = arg_checker.verify_input_files(arg_dict)
    # Run PCGR workflow (vep, vcfanno, summarise, vcf2tsv, html)
    run_pcgr(pcgr_paths, config_options)

def run_pcgr(pcgr_paths, config_options):
    """
    Main function to run the PCGR workflow
    """

    debug = config_options['debug']

    report_nonfloating_toc = 1 if config_options['other']['nonfloating_toc'] else 0
    vep_regulatory_annotation = 'ON' if config_options['other']['vep_regulatory'] == 1 else 'OFF'
    clinical_trials_set = 'ON' if config_options['clinicaltrials']['run'] else 'OFF'
    msi_prediction_set = 'ON' if config_options['msi']['run'] else 'OFF'
    msig_estimation_set = 'ON' if config_options['msigs']['run'] else 'OFF'
    tmb_estimation_set = 'ON' if config_options['tmb']['run'] else 'OFF'
    vcf_validation = 0 if config_options['other']['no_vcf_validate'] else 1
    run_vcf2maf = config_options['other']['vcf2maf']
    assay_mode = 'Tumor vs. Control'
    tumor_only = 0
    cell_line = 0
    if config_options['tumor_only']['tumor_only']:
        assay_mode = 'Tumor-Only'
        tumor_only = 1
        if config_options['tumor_only']['cell_line']:
            cell_line = 1
            assay_mode = 'Tumor-Only (cell line)'
    # set basic run commands
    output_vcf = 'None'
    output_pass_vcf = 'None'
    output_pass_tsv = 'None'
    output_maf = 'None'
    GENCODE_VERSION = pcgr_vars.GENCODE_VERSION
    NCBI_BUILD_MAF = pcgr_vars.NCBI_BUILD_MAF
    VEP_ASSEMBLY = pcgr_vars.VEP_ASSEMBLY
    MAX_VARIANTS_FOR_REPORT = pcgr_vars.MAX_VARIANTS_FOR_REPORT
    if config_options['genome_assembly'] == 'grch37':
        NCBI_BUILD_MAF = 'GRCh37'
        GENCODE_VERSION = 'release 19'
        VEP_ASSEMBLY = 'GRCh37'
    logger = getlogger('pcgr-get-OS')

    vep_dir = os.path.join(str(pcgr_paths['db_dir']), '.vep')
    input_vcf = 'None'
    input_cna = 'None'
    input_rna_fusion = 'None'
    input_rna_expression = 'None'
    input_cpsr_report = 'None'
    panel_normal = 'None'
    # panel-of-normals annotation
    pon_annotation = 0

    # Specify paths for input files and directories
    if pcgr_paths['input_vcf_basename'] != 'NA':
        input_vcf = os.path.join(pcgr_paths['input_vcf_dir'], pcgr_paths['input_vcf_basename'])
    if pcgr_paths['input_cna_basename'] != 'NA':
        input_cna = os.path.join(pcgr_paths['input_cna_dir'], pcgr_paths['input_cna_basename'])
    if pcgr_paths['input_rna_fusion_basename'] != 'NA':
        input_rna_fusion = os.path.join(pcgr_paths['input_rna_fusion_dir'], pcgr_paths['input_rna_fusion_basename'])
    if pcgr_paths['input_rna_expression_basename'] != 'NA':
        input_rna_expression = os.path.join(pcgr_paths['input_rna_expression_dir'], pcgr_paths['input_rna_expression_basename'])
    if pcgr_paths['input_cpsr_report_basename'] != 'NA':
        input_cpsr_report = os.path.join(pcgr_paths['input_cpsr_report_dir'], pcgr_paths['input_cpsr_report_basename'])
    if pcgr_paths['panel_normal_vcf_basename'] != 'NA':
        panel_normal = os.path.join(pcgr_paths['panel_normal_vcf_dir'], pcgr_paths['panel_normal_vcf_basename'])

    data_dir = pcgr_paths['base_dir']
    output_dir = pcgr_paths['output_dir']

    # PCGR|validate_input - verify that VCF and CNA segment file is of appropriate format
    logger = getlogger("pcgr-validate-input-arguments")
    logger.info("PCGR - STEP 0: Validate input data and options")

    vcf_validate_command = (
            f'pcgr_validate_input.py '
            f'{data_dir} '
            f'{input_vcf} '
            f'{input_cna} '
            f'{input_rna_fusion} '
            f'{input_rna_expression} '
            f'{panel_normal} '
            f'{vcf_validation} '
            f'{tumor_only} '
            f'{config_options["genome_assembly"]} '
            f'{config_options["other"]["preserved_info_tags"]} '
            f'{config_options["allelic_support"]["tumor_dp_tag"]} {config_options["allelic_support"]["tumor_af_tag"]} '
            f'{config_options["allelic_support"]["control_dp_tag"]} {config_options["allelic_support"]["control_af_tag"]} '
            f'{config_options["allelic_support"]["call_conf_tag"]} '
            f'{config_options["tumor_only"]["exclude_likely_hom_germline"]} '
            f'{config_options["tumor_only"]["exclude_likely_het_germline"]} '
            f'--output_dir {output_dir} '
            f'{"--debug " if debug else ""}'
            f'{"--keep_uncompressed" if run_vcf2maf else ""} '
            )
    check_subprocess(logger, vcf_validate_command, debug)
    logger.info('Finished pcgr-validate-input-arguments')
    print('----')

    # PCGR|start - Log key information about sample, options and sequencing assay/design
    logger = getlogger('pcgr-start')
    logger.info('--- Personal Cancer Genome Reporter workflow ----')
    logger.info(f'Sample name: {config_options["sample_id"]}')
    if config_options['tumor_type']['type'] == 'Cancer_NOS':
        logger.info('Tumor type: Cancer_NOS (Any tumortype)')
    else:
        logger.info(f'Tumor type: {config_options["tumor_type"]["type"]}')
    logger.info(f'Sequencing assay - type: {config_options["assay"]}')
    logger.info(f'Sequencing assay - mode: {assay_mode}')
    logger.info(f'Sequencing assay - coding target size: {config_options["tmb"]["target_size_mb"]}Mb')
    logger.info(f'Genome assembly: {config_options["genome_assembly"]}')
    logger.info(f'Mutational signature estimation: {msig_estimation_set}')
    logger.info(f'MSI classification: {msi_prediction_set}')
    logger.info(f'Mutational burden estimation: {tmb_estimation_set}')
    logger.info(f'Include molecularly targeted clinical trials (beta): {clinical_trials_set}')

    if not input_vcf == 'None':
        # Define temporary output file names
        prefix = os.path.join(output_dir, f'{config_options["sample_id"]}.pcgr_acmg.{config_options["genome_assembly"]}')
        output_vcf =             f'{prefix}.vcf.gz'
        output_pass_vcf =        f'{prefix}.pass.vcf.gz'
        output_pass_tsv =        f'{prefix}.pass.tsv'
        output_pass_raw_tsv_gz = f'{prefix}.pass.raw.tsv.gz'
        output_maf =             f'{prefix}.tmp.maf'
        output_vcf2maf_log =     f'{prefix}.maf.log'
        input_vcf_pcgr_ready =   os.path.join(output_dir, re.sub(r"(\.vcf$|\.vcf\.gz$)", ".pcgr_ready.vcf.gz", pcgr_paths["input_vcf_basename"]))
        # needs to be uncompressed for vcf2maf
        input_vcf_pcgr_ready_uncompressed = os.path.join(output_dir, re.sub(r"(\.vcf$|\.vcf\.gz$)", ".pcgr_ready.vcf", pcgr_paths["input_vcf_basename"]))
        vep_vcf = re.sub(r"(\.vcf$|\.vcf\.gz$)", ".vep.vcf.gz", input_vcf_pcgr_ready)
        vep_vcfanno_vcf = re.sub(r"(\.vcf$|\.vcf\.gz$)", ".vep.vcfanno.vcf", input_vcf_pcgr_ready)
        vep_vcfanno_annotated_vcf = re.sub(r"\.vcfanno", ".vcfanno.annotated", vep_vcfanno_vcf) + ".gz"
        vep_vcfanno_annotated_pass_vcf = re.sub(r"\.vcfanno", ".vcfanno.annotated.pass", vep_vcfanno_vcf) + ".gz"
        fasta_assembly = os.path.join(vep_dir, 'homo_sapiens', f'{pcgr_vars.VEP_VERSION}_{VEP_ASSEMBLY}', f'Homo_sapiens.{VEP_ASSEMBLY}.dna.primary_assembly.fa.gz')
        # List all VEP flags used when calling VEP
        vep_flags = (
                f'--hgvs --af --af_1kg --af_gnomad --variant_class --domains --symbol --protein --ccds --mane '
                f'--uniprot --appris --biotype --tsl --canonical --format vcf --cache --numbers --total_length --allele_number '
                f'--no_stats --no_escape --xref_refseq --vcf --check_ref --dont_skip --flag_pick_allele --plugin NearestExonJB,max_range=50000 '
                f'--force_overwrite --species homo_sapiens --offline --compress_output bgzip'
                )
        vep_options = (
                f'--dir {vep_dir} --assembly {VEP_ASSEMBLY} --cache_version {pcgr_vars.VEP_VERSION} '
                f'--fasta {fasta_assembly} --pick_order {config_options["other"]["vep_pick_order"]} '
                f'--buffer_size {config_options["other"]["vep_buffer_size"]} '
                f'--fork {config_options["other"]["vep_n_forks"]} '
                f'{vep_flags} '
                f'{"--verbose" if debug else "--quiet"} '
                )
        gencode_set_in_use = "GENCODE - all transcripts"
        if config_options['other']['vep_no_intergenic'] == 1:
            vep_options += '--no_intergenic '
        if config_options['other']['vep_regulatory'] == 1:
            vep_options += '--regulatory '
        if config_options['other']['vep_gencode_all'] == 0:
            vep_options += '--gencode_basic '
            gencode_set_in_use = "GENCODE - basic transcript set (--gencode_basic)"

        # Compose full VEP command
        vep_main_command = f'{utils.get_perl_exports()} && vep --input_file {input_vcf_pcgr_ready} --output_file {vep_vcf} {vep_options}'
        vep_tabix_command = f'tabix -f -p vcf {vep_vcf}'

        # PCGR|VEP - run consequence annotation with Variant Effect Predictor
        print('----')
        logger = getlogger('pcgr-vep')
        logger.info(f'PCGR - STEP 1: Basic variant annotation with Variant Effect Predictor ({pcgr_vars.VEP_VERSION}, GENCODE {GENCODE_VERSION}, {config_options["genome_assembly"]})')
        logger.info(f'VEP configuration - one primary consequence block pr. alternative allele (--flag_pick_allele)')
        logger.info(f'VEP configuration - transcript pick order: {config_options["other"]["vep_pick_order"]}')
        logger.info(f'VEP configuration - transcript pick order: See more at https://www.ensembl.org/info/docs/tools/vep/script/vep_other.html#pick_options')
        logger.info(f'VEP configuration - GENCODE set: {gencode_set_in_use}')
        logger.info(f'VEP configuration - skip intergenic: {"TRUE" if config_options["other"]["vep_no_intergenic"] else "FALSE"}')
        logger.info(f'VEP configuration - regulatory annotation: {vep_regulatory_annotation}')
        logger.info(f'VEP configuration - buffer_size/number of forks: {config_options["other"]["vep_buffer_size"]}/{config_options["other"]["vep_n_forks"]}')

        check_subprocess(logger, vep_main_command, debug)
        check_subprocess(logger, vep_tabix_command, debug)
        logger.info('Finished pcgr-vep')
        print('----')

        # PCGR|vcf2maf - if option set, convert VCF to MAF with https://github.com/mskcc/vcf2maf
        if run_vcf2maf:
            logger.info('Converting VEP-annotated VCF to MAF with https://github.com/mskcc/vcf2maf')
            vcf2maf_command = (
                    f'vcf2maf.pl --inhibit-vep --input-vcf {input_vcf_pcgr_ready_uncompressed} '
                    f'--tumor-id {config_options["sample_id"]} --output-maf {output_maf} --ref-fasta {fasta_assembly} '
                    f'--ncbi-build {NCBI_BUILD_MAF} > {output_vcf2maf_log} 2>&1'
                    )
            check_subprocess(logger, vcf2maf_command, debug)
            utils.remove(input_vcf_pcgr_ready_uncompressed)
            utils.remove(output_vcf2maf_log)
            logger.info('Finished pcgr-vep-vcf2maf')
            print('----')

        # PCGR|vcfanno - annotate VCF against a number of variant annotation resources
        logger = getlogger("pcgr-vcfanno")
        pcgr_vcfanno_command = (
                f'pcgr_vcfanno.py {vep_vcf} {vep_vcfanno_vcf} {pcgr_paths["db_dir"]} '
                f'--num_processes {config_options["other"]["vcfanno_n_proc"]} '
                f'--chasmplus --dbnsfp --docm --clinvar --icgc --civic --cgi --tcga_pcdm --winmsk --simplerepeats '
                f'--tcga --uniprot --cancer_hotspots --pcgr_onco_xref '
                f'{"--debug " if debug else ""}'
                )
        anno_src_msg = (
                f"Annotation sources: {'Panel-of-Normals, ' if panel_normal != 'None' else ''}ClinVar, dbNSFP, "
                f"UniProtKB, cancerhotspots.org, CiVIC, CGI, DoCM, CHASMplus driver mutations, TCGA, ICGC-PCAWG"
                )
        logger.info("PCGR - STEP 2: Annotation for precision oncology with pcgr-vcfanno")
        logger.info(anno_src_msg)
        if panel_normal != "None":
            pon_annotation = 1
            pcgr_vcfanno_command += f'--panel_normal_vcf {panel_normal}'
        check_subprocess(logger, pcgr_vcfanno_command, debug)
        logger.info("Finished pcgr-vcfanno")
        print('----')

        # PCGR|pcgr_summarise - expand annotations in VCF file
        logger = getlogger("pcgr-summarise")
        pcgr_summarise_command = (
                f'pcgr_summarise.py {vep_vcfanno_vcf}.gz {pon_annotation} '
                f'{config_options["other"]["vep_regulatory"]} '
                f'{pcgr_paths["db_dir"]} '
                f'{"--debug" if debug else ""}'
                )
        logger.info("PCGR - STEP 3: Cancer gene annotations with pcgr-summarise")
        check_subprocess(logger, pcgr_summarise_command, debug)

        # PCGR|clean - move output files and clean up temporary files
        os.rename(vep_vcfanno_annotated_vcf, output_vcf)
        os.rename(f'{vep_vcfanno_annotated_vcf}.tbi', f'{output_vcf}.tbi')
        os.rename(vep_vcfanno_annotated_pass_vcf, output_pass_vcf)
        os.rename(f'{vep_vcfanno_annotated_pass_vcf}.tbi', f'{output_pass_vcf}.tbi')
        delete_files = (
                glob(f'{vep_vcf}*') +
                glob(f'{vep_vcfanno_annotated_vcf}') +
                glob(f'{vep_vcfanno_annotated_pass_vcf}*') +
                glob(f'{vep_vcfanno_vcf}*') +
                glob(f'{input_vcf_pcgr_ready_uncompressed}*')
                )
        # do not delete if debugging
        if not debug:
            for fn in delete_files:
                #print(f"Deleting {fn}")
                utils.remove(fn)

        logger.info('Finished pcgr-summarise main command')

        # PCGR|vcf2tsv - convert VCF to TSV with https://github.com/sigven/vcf2tsv
        pcgr_vcf2tsv_command = f'vcf2tsv.py {output_pass_vcf} --compress {output_pass_tsv}'
        logger.info("Converting VCF to TSV with https://github.com/sigven/vcf2tsv")
        check_subprocess(logger, pcgr_vcf2tsv_command, debug)
        logger.info('Finished pcgr-summarise-vcf2tsv')

        if config_options['assay'] == 'WGS' or config_options['assay'] == 'WES':
            output_pass_tsv_gz = f'{output_pass_tsv}.gz'
            # check that output file exist
            if os.path.exists(output_pass_tsv_gz):
                # get number of rows/variants annotated, using pandas
                var_data = pandas.read_csv(output_pass_tsv_gz, sep = '\t', low_memory = False, header = [1])
                num_variants_raw = len(var_data)
                if num_variants_raw > MAX_VARIANTS_FOR_REPORT:
                    logger.info(f'Number of raw variants in input VCF ({num_variants_raw}) exceeds {MAX_VARIANTS_FOR_REPORT} - intergenic/intronic variants will be excluded prior to reporting')

                    # Exclude intronic and intergenic variants prior to analysis with pcgrr (reporting and further analysis)
                    var_data_filtered = var_data[~var_data.Consequence.str.contains('^intron') & ~var_data.Consequence.str.contains('^intergenic')]
                    num_variants_excluded1 = num_variants_raw - len(var_data_filtered)
                    logger.info(f'Number of intergenic/intronic variants excluded: {num_variants_excluded1}')

                    # Exclude upstream_gene/downstream_gene variants if size of filtered variant set is still above MAX_VARIANTS_FOR_REPORT
                    # TODO: in this case, the TMB calculation will be an underestimate (but still likely huge)
                    var_data_filtered_final = var_data_filtered
                    if len(var_data_filtered) > MAX_VARIANTS_FOR_REPORT:
                        var_data_filtered_final = var_data_filtered[~var_data_filtered.Consequence.str.contains('^upstream_gene') & ~var_data_filtered.Consequence.str.contains('^downstream_gene')]
                        num_variants_excluded2 = len(var_data_filtered) - len(var_data_filtered_final)
                        logger.info(f'Number of upstream_gene/downstream_gene variants excluded: {num_variants_excluded2}')


                    # get vcf2tsv header and pipe to output TSV file
                    get_vcf2tsv_header = f'gzip -dc {output_pass_tsv_gz} | egrep \'^#\' > {output_pass_tsv}'
                    check_subprocess(logger, get_vcf2tsv_header, debug)

                    # rename original vcf2tsv (gzipped) to 'raw' filename
                    rename_output_tsv = f'mv {output_pass_tsv_gz} {output_pass_raw_tsv_gz}'
                    check_subprocess(logger, rename_output_tsv, debug)

                    # append filtered data output to output TSV file
                    var_data_filtered_final.to_csv(output_pass_tsv, sep='\t', encoding='utf-8', mode = 'a', index = False)

                    # gzip filtered output TSV file
                    gzip_filtered_output_tsv = f'gzip -f {output_pass_tsv}'
                    check_subprocess(logger, gzip_filtered_output_tsv, debug)


        logger.info('Finished pcgr-summarise')
        print('----')

    # Generation of HTML reports for VEP/vcfanno-annotated VCF and copy number segment file
    if not config_options['other']['basic']:
        co = config_options
        ttype = co['tumor_type']['type'].replace(' ', '_').replace('/', '@')
        logger = getlogger('pcgr-writer')
        logger.info('PCGR - STEP 4: Generation of output files - variant interpretation report for precision oncology')

        # export PATH to R conda env Rscript
        rscript = utils.script_path('pcgrr', 'bin/Rscript')
        pcgrr_script = utils.script_path('pcgr', 'bin/pcgrr.R')
        pcgr_report_command = (
                f"{rscript} {pcgrr_script} "
                f"{output_dir} "
                f"{output_pass_tsv}.gz "
                f"{input_cna} "
                f"{input_rna_fusion} "
                f"{input_rna_expression} "
                f"{input_cpsr_report} "
                f"{config_options['sample_id']} "
                f"{pcgr_vars.PCGR_VERSION} "
                f"{pcgr_vars.DB_VERSION} "
                f"{config_options['genome_assembly']} "
                f"{data_dir} "
                f"{co['tumor_purity']} "
                f"{co['tumor_ploidy']} "
                f"{ttype} "
                f"{co['tmb']['target_size_mb']} "
                f"{co['assay']} "
                f"{tumor_only} "
                f"{cell_line} "
                f"{co['tumor_only']['maf_onekg_afr']} "
                f"{co['tumor_only']['maf_onekg_amr']} "
                f"{co['tumor_only']['maf_onekg_eas']} "
                f"{co['tumor_only']['maf_onekg_eur']} "
                f"{co['tumor_only']['maf_onekg_sas']} "
                f"{co['tumor_only']['maf_onekg_global']} "
                f"{co['tumor_only']['maf_gnomad_afr']} "
                f"{co['tumor_only']['maf_gnomad_amr']} "
                f"{co['tumor_only']['maf_gnomad_asj']} "
                f"{co['tumor_only']['maf_gnomad_eas']} "
                f"{co['tumor_only']['maf_gnomad_fin']} "
                f"{co['tumor_only']['maf_gnomad_nfe']} "
                f"{co['tumor_only']['maf_gnomad_oth']} "
                f"{co['tumor_only']['maf_gnomad_sas']} "
                f"{co['tumor_only']['maf_gnomad_global']} "
                f"{co['tumor_only']['exclude_pon']} "
                f"{co['tumor_only']['exclude_likely_hom_germline']} "
                f"{co['tumor_only']['exclude_likely_het_germline']} "
                f"{co['tumor_only']['exclude_dbsnp_nonsomatic']} "
                f"{co['tumor_only']['exclude_nonexonic']} "
                f"{co['tmb']['run']} "
                f"{co['tmb']['algorithm']} "
                f"{co['msi']['run']} "
                f"{co['msigs']['run']} "
                f"{co['msigs']['mutation_limit']} "
                f"{co['msigs']['all_reference_signatures']} "
                f"{co['msigs']['include_artefact_signatures']} "
                f"{co['msigs']['prevalence_reference_signatures']} "
                f"{co['cna']['logR_homdel']} "
                f"{co['cna']['logR_gain']} "
                f"{co['cna']['cna_overlap_pct']} "
                f"{co['allelic_support']['tumor_af_min']} "
                f"{co['allelic_support']['tumor_dp_min']} "
                f"{co['allelic_support']['control_dp_min']} "
                f"{co['allelic_support']['control_af_max']} "
                f"{co['allelic_support']['tumor_af_tag']} "
                f"{co['allelic_support']['tumor_dp_tag']} "
                f"{co['allelic_support']['control_af_tag']} "
                f"{co['allelic_support']['control_dp_tag']} "
                f"{co['allelic_support']['call_conf_tag']} "
                f"{co['clinicaltrials']['run']} "
                f"{co['other']['vep_n_forks']} "
                f"{co['other']['vep_buffer_size']} "
                f"{co['other']['vep_no_intergenic']} "
                f"{co['other']['vep_pick_order']} "
                f"{co['other']['vep_regulatory']} "
                f"{co['other']['vep_gencode_all']} "
                f"{co['other']['vcf2maf']} "
                f"{co['other']['list_noncoding']} "
                f"{co['other']['preserved_info_tags']} "
                f"{co['other']['visual_theme']} "
                f"{report_nonfloating_toc} "
                f"{co['other']['no_vcf_validate']}"
                )

        if debug:
            print(pcgr_report_command)
        check_subprocess(logger, pcgr_report_command, debug)
        logger.info("Finished PCGR!")
        print('----')

    print()


if __name__ == "__main__":
    cli()
