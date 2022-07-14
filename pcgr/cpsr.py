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
from glob import glob
from argparse import RawTextHelpFormatter
from pcgr import pcgr_vars, arg_checker, utils
from pcgr.utils import check_subprocess, getlogger, error_message, warn_message

def get_args():

    program_description = "Cancer Predisposition Sequencing Reporter - report of " + \
       "clinically significant cancer-predisposing germline variants"
    program_options = " --input_vcf <INPUT_VCF> --pcgr_dir <PCGR_DIR> --output_dir <OUTPUT_DIR> --genome_assembly " + \
       " <GENOME_ASSEMBLY> --sample_id <SAMPLE_ID>"

    parser = argparse.ArgumentParser(description = program_description,
                                     formatter_class=RawTextHelpFormatter, usage="%(prog)s -h [options] " + str(program_options))
    parser._action_groups.pop()
    required = parser.add_argument_group('Required arguments')
    optional_panel = parser.add_argument_group("Panel options")
    optional_vep = parser.add_argument_group('VEP options')
    optional_vcfanno = parser.add_argument_group('vcfanno options')
    optional_other = parser.add_argument_group('Other options')

    optional_panel.add_argument('--panel_id',dest = "virtual_panel_id",type = str, default = "-1", help="Comma-separated string with identifier(s) of predefined virtual cancer predisposition gene panels,\nchoose any combination of the following identifiers (GEP = Genomics England PanelApp):\n" + str(pcgr_vars.panels))
    optional_panel.add_argument('--custom_list',dest = "custom_list",help="Provide custom list of genes from virtual panel 0 (single-column txt file with Ensembl gene identifiers),\n alternative to predefined panels provided with --panel_id)")
    optional_panel.add_argument('--custom_list_name',dest = "custom_list_name", default="Custom_Panel", help="Set name for custom made panel/list (single word - no whitespace), will be displayed in the report")
    optional_panel.add_argument('--diagnostic_grade_only', action="store_true",help="For panel_id's 1-42 (Genomics England PanelApp) - consider genes with a GREEN status only, default: %(default)s")

    optional_other.add_argument('--force_overwrite', action = "store_true", help='By default, the script will fail with an error if any output file already exists.\n You can force the overwrite of existing result files by using this flag, default: %(default)s')
    optional_other.add_argument('--version', action='version', version=str(utils.get_cpsr_version()))
    optional_other.add_argument('--basic',action="store_true",help="Run functional variant annotation on VCF through VEP/vcfanno, omit Tier assignment/report generation (STEP 4), default: %(default)s")
    optional_other.add_argument('--no_vcf_validate', action = "store_true",help="Skip validation of input VCF with Ensembl's vcf-validator, default: %(default)s")
    optional_other.add_argument('--preserved_info_tags', dest ='preserved_info_tags', default='None', help='Comma-separated string of VCF INFO tags from query VCF that should be kept in CPSR output TSV')
    optional_other.add_argument('--report_theme',choices = ['default','cerulean','journal','flatly','readable','spacelab','united','cosmo','lumen','paper','sandstone','simplex','yeti'], default = 'default', help='Visual report theme (rmarkdown),  default: %(default)s' )
    optional_other.add_argument('--report_nonfloating_toc', action='store_true', help='Do not float the table of contents (TOC) in output HTML report, default: %(default)s')
    optional_other.add_argument('--report_table_display', choices = ['full','light'], default='light', help="Set the level of detail/comprehensiveness in interactive datables of HTML report, very comprehensive (option 'full') or slim/focused ('light'), default: %(default)s")
    optional_other.add_argument('--ignore_noncoding', action='store_true',dest='ignore_noncoding',default=False,help='Do not list non-coding variants in HTML report, default: %(default)s')
    optional_other.add_argument('--secondary_findings', action='store_true',dest='secondary_findings',default=False, help='Include variants found in ACMG-recommended list for secondary findings (v3.0), default: %(default)s')
    optional_other.add_argument('--gwas_findings', action='store_true',dest='gwas_findings',default=False, help='Report overlap with low to moderate cancer risk variants (tag SNPs) identified from genome-wide association studies, default: %(default)s')
    optional_other.add_argument('--gwas_p_value', type = float, default = 0.000005, dest = 'gwas_p_value',help='Required p-value for variants listed as hits from genome-wide association studies, default: %(default)s')
    optional_other.add_argument('--pop_gnomad',choices = ['afr','amr','eas','sas','asj','nfe','fin','global'], default='nfe', help='Population source in gnomAD used for variant frequency assessment (ACMG classification), default: %(default)s')
    optional_other.add_argument('--maf_upper_threshold', type = float, default = 0.9, dest = 'maf_upper_threshold',help='Upper MAF limit (gnomAD global population frequency) for variants to be included in the report, default: %(default)s')
    optional_other.add_argument('--classify_all', action='store_true',dest='classify_all',help='Provide CPSR variant classifications (TIER 1-5) also for variants with existing ClinVar classifications in output TSV, default: %(default)s')
    optional_other.add_argument('--clinvar_ignore_noncancer', action='store_true', help='Ignore (exclude from report) ClinVar-classified variants reported only for phenotypes/conditions NOT related to cancer, default: %(default)s')
    optional_other.add_argument("--debug", action="store_true", help="Print full commands to log")

    optional_vcfanno.add_argument('--vcfanno_n_proc', default = 4, type = int, help="Number of vcfanno processes (option '-p' in vcfanno), default: %(default)s")

    optional_vep.add_argument('--vep_n_forks', default = 4, type = int, help="Number of forks (option '--fork' in VEP), default: %(default)s")
    optional_vep.add_argument('--vep_buffer_size', default = 500, type = int, help="Variant buffer size (variants read into memory simultaneously, option '--buffer_size' in VEP) " + \
       "\n- set lower to reduce memory usage, default: %(default)s")
    #optional_vep.add_argument('--vep_regulatory', action='store_true', help = 'Enable Variant Effect Predictor (VEP) to look for overlap with regulatory regions (option --regulatory in VEP).')
    optional_vep.add_argument('--vep_gencode_all', action='store_true', help = "Consider all GENCODE transcripts with Variant Effect Predictor (VEP) (option '--gencode_basic' in VEP is used by default).")
    optional_vep.add_argument('--vep_pick_order', default = "canonical,appris,biotype,ccds,rank,tsl,length,mane", help="Comma-separated string " + \
       "of ordered transcript properties for primary variant pick\n ( option '--pick_order' in VEP), default: %(default)s")
    optional_vep.add_argument('--vep_no_intergenic', action = "store_true", help="Skip intergenic variants during processing (option '--no_intergenic' in VEP), default: %(default)s")

    required.add_argument('--input_vcf', help='VCF input file with germline query variants (SNVs/InDels).', required = True)
    required.add_argument('--pcgr_dir',help=f"Directory that contains the PCGR data bundle directory, e.g. ~/pcgr-{pcgr_vars.PCGR_VERSION}", required = True)
    required.add_argument('--output_dir',help='Output directory', required = True)
    required.add_argument('--genome_assembly',choices = ['grch37','grch38'], help='Genome assembly build: grch37 or grch38', required = True)
    required.add_argument('--sample_id',help="Sample identifier - prefix for output files", required = True)

    args = parser.parse_args()
    return vars(args)

def main():
    arg_dict = get_args()
    # check parsed arguments
    arg_checker.check_args_cpsr(arg_dict)
    # Verify existence of input files
    cpsr_paths = arg_checker.verify_input_files_cpsr(arg_dict)
    ## Run CPSR workflow
    run_cpsr(arg_dict, cpsr_paths)


def run_cpsr(arg_dict, cpsr_paths):
    """
    Main function to run the CPSR workflow
    """
    debug = arg_dict['debug']
    diagnostic_grade_only = 0
    vcf_validation = 1
    virtual_panel_id = "-1"
    ignore_noncoding = 0
    gwas_findings = 0
    secondary_findings = 0
    classify_all = 0
    clinvar_ignore_noncancer = 0
    report_nonfloating_toc = 0
    vep_no_intergenic = 0
    vep_regulatory = 0
    preserved_info_tags = arg_dict['preserved_info_tags']
    diagnostic_grade_set = "OFF"
    secondary_findings_set = "OFF"
    gwas_findings_set = "OFF"

    if arg_dict['vep_regulatory']:
        vep_regulatory = 1
    if arg_dict["vep_no_intergenic"]:
        vep_no_intergenic = 1
    if arg_dict['clinvar_ignore_noncancer']:
        clinvar_ignore_noncancer = 1
    if arg_dict['classify_all']:
        classify_all = 1
    if arg_dict['gwas_findings']:
        gwas_findings = 1
        gwas_findings_set = "ON"
    if arg_dict['secondary_findings']:
        secondary_findings = 1
        secondary_findings_set = "ON"
    if arg_dict['diagnostic_grade_only']:
        diagnostic_grade_only = 1
        diagnostic_grade_set = "ON"
    if arg_dict['report_nonfloating_toc']:
        report_nonfloating_toc = 1
    if arg_dict['no_vcf_validate']:
        vcf_validation = 0
    if arg_dict['virtual_panel_id'] != "-1":
        virtual_panel_id = arg_dict['virtual_panel_id']
    if arg_dict['custom_list']:
        virtual_panel_id = "-1"
    if arg_dict['ignore_noncoding']:
        ignore_noncoding = 1

    output_vcf = 'None'
    output_pass_vcf = 'None'
    output_pass_tsv = 'None'
    uid = ''
    GENCODE_VERSION = pcgr_vars.GENCODE_VERSION
    VEP_ASSEMBLY = pcgr_vars.VEP_ASSEMBLY
    VEP_VERSION = pcgr_vars.VEP_VERSION
    if arg_dict['genome_assembly'] == 'grch37':
        GENCODE_VERSION = '19'
        VEP_ASSEMBLY = 'GRCh37'

    vepdb_dir = os.path.join(str(cpsr_paths['db_dir']),'.vep')
    input_vcf = 'None'
    input_customlist = 'None'

    if cpsr_paths['input_vcf_basename'] != 'NA':
        input_vcf = os.path.join(cpsr_paths['input_vcf_dir'], cpsr_paths['input_vcf_basename'])
    if cpsr_paths['input_customlist_basename'] != 'NA':
        input_customlist = os.path.join(cpsr_paths['input_customlist_dir'], cpsr_paths['input_customlist_basename'])

    data_dir = cpsr_paths['base_dir']
    output_dir = cpsr_paths['output_dir']
    vep_dir = vepdb_dir

    logger = getlogger('cpsr-validate-input-arguments')
    logger.info("CPSR - STEP 0: Validate input data")
    check_subprocess(logger, f'mkdir -p {output_dir}', debug)

    ## CPSR|Validate input VCF - check formatting, non-overlap with CPSR INFO tags, and whether sample contains any variants in cancer predisposition loci
    vcf_validate_command = (
            f'cpsr_validate_input.py '
            f'{data_dir} '
            f'{input_vcf} '
            f'{input_customlist} '
            f'{preserved_info_tags} '
            f'{vcf_validation} '
            f'{arg_dict["genome_assembly"]} '
            f'{arg_dict["sample_id"]} '
            f'{virtual_panel_id} '
            f'{diagnostic_grade_only} '
            f'--output_dir {output_dir} {"--debug" if debug else ""}'
            )
    check_subprocess(logger, vcf_validate_command, debug)
    logger.info('Finished cpsr-validate-input-arguments')
    print('----')

    ## CPSR|Start - log key information about run
    logger = getlogger("cpsr-start")
    logger.info("--- Cancer Predisposition Sequencing Reporter workflow ----")
    logger.info(f"Sample name: {arg_dict['sample_id']}")
    if not input_customlist == 'None':
        logger.info(f"Virtual gene panel: custom-made list from panel 0: {input_customlist}")
    else:
        #logger.info("Virtual gene panel(s): " + str(pcgr_vars.GE_panels[virtual_panel_id]))
        logger.info(f"Diagnostic-grade genes in virtual panels (GE PanelApp): {diagnostic_grade_set}")
    logger.info(f"Include incidental findings (ACMG recommended list v3.0): {secondary_findings_set}")
    logger.info(f"Include low to moderate cancer risk variants from genome-wide association studies: {gwas_findings_set}")
    logger.info(f"Reference population, germline variant frequencies (gnomAD): {str(arg_dict['pop_gnomad']).upper()}")
    logger.info(f"Genome assembly: {arg_dict['genome_assembly']}")

    if not input_vcf == 'None':

        ## Define input, output and temporary file names
        pcgr_model = 'cpsr'
        output_vcf = os.path.join(output_dir, str(arg_dict['sample_id']) + '.cpsr.' + str(arg_dict['genome_assembly']) + '.vcf.gz')
        output_pass_vcf = os.path.join(output_dir, str(arg_dict['sample_id']) + '.cpsr.' + str(arg_dict['genome_assembly']) + '.pass.vcf.gz')
        output_pass_tsv = os.path.join(output_dir, str(arg_dict['sample_id']) + '.cpsr.' + str(arg_dict['genome_assembly']) + '.pass.tsv')
        input_vcf_cpsr_ready = os.path.join(output_dir, re.sub(r'(\.vcf$|\.vcf\.gz$)','.cpsr_ready_target.vcf.gz', cpsr_paths['input_vcf_basename']))
        input_vcf_cpsr_ready_uncompressed = os.path.join(output_dir, re.sub(r'(\.vcf$|\.vcf\.gz$)','.cpsr_ready_target.vcf', cpsr_paths['input_vcf_basename']))
        vep_vcf = re.sub(r'(\.vcf$|\.vcf\.gz$)','.cpsr_vep.vcf',input_vcf_cpsr_ready)
        vep_vcfanno_vcf = re.sub(r'(\.vcf$|\.vcf\.gz$)','.cpsr_vep.vcfanno.vcf',input_vcf_cpsr_ready)
        vep_vcfanno_annotated_vcf = re.sub(r'\.vcfanno','.vcfanno.annotated',vep_vcfanno_vcf) + '.gz'
        vep_vcfanno_annotated_pass_vcf = re.sub(r'\.vcfanno','.vcfanno.annotated.pass',vep_vcfanno_vcf) + '.gz'
        custom_bed = os.path.join(output_dir, str(arg_dict['sample_id']) + '.' + str(pcgr_model) + '.' + str(arg_dict['genome_assembly']) + '.custom_list.bed')

        ## File names for assembly-specific genome fasta files (VEP)
        fasta_assembly = os.path.join(vep_dir, f"homo_sapiens/{VEP_VERSION}_{VEP_ASSEMBLY}/Homo_sapiens.{VEP_ASSEMBLY}.dna.primary_assembly.fa.gz")
        ancestor_assembly = os.path.join(vep_dir, f"homo_sapiens/{VEP_VERSION}_{VEP_ASSEMBLY}/human_ancestor.fa.gz")

        ## Set all flags used in VEP run
        plugins_in_use = "NearestExonJB, LoF"
        vep_flags = (
            f"--format vcf --vcf --check_ref --flag_pick_allele_gene --hgvs --dont_skip --failed 1 --af --af_1kg --af_gnomad "
            f"--variant_class --domains --symbol --protein --ccds --uniprot --appris --biotype --canonical --cache "
            f"--numbers --total_length --no_stats --allele_number --no_escape --xref_refseq --plugin NearestExonJB,max_range=50000"
            )
        vep_options = (
            f"--pick_order {arg_dict['vep_pick_order']} --force_overwrite --buffer_size {arg_dict['vep_buffer_size']} "
            f"--species homo_sapiens --assembly {VEP_ASSEMBLY} --offline --fork {arg_dict['vep_n_forks']} {vep_flags} --dir {vep_dir} "
            f"--cache_version {VEP_VERSION}"
            )
        gencode_set_in_use = "GENCODE - all transcripts"
        if arg_dict['vep_gencode_all'] == 0:
            vep_options += ' --gencode_basic'
            gencode_set_in_use = "GENCODE - basic transcript set (--gencode_basic)"
        if arg_dict['vep_no_intergenic'] == 1:
            vep_options = vep_options + " --no_intergenic"
        if arg_dict['vep_regulatory'] == 1:
            vep_options = vep_options + " --regulatory"
        if arg_dict['genome_assembly'] == "grch38":
            vep_options = vep_options +  " --mane"
        loftee_dir = utils.get_loftee_dir()
        assert os.path.isdir(loftee_dir), f'LoF VEP plugin is not found in {loftee_dir}. Please make sure you installed pcgr conda package and have corresponding conda environment active.'
        vep_options += f" --plugin LoF,loftee_path:{loftee_dir},human_ancestor_fa:{ancestor_assembly},use_gerp_end_trunc:0 --dir_plugins {loftee_dir}"
        if not debug:
            vep_options += " --quiet"

        ## Compose full VEP command
        vep_main_command = f'{utils.get_perl_exports()} && vep --input_file {input_vcf_cpsr_ready} --output_file {vep_vcf} {vep_options} --fasta {fasta_assembly}'
        vep_bgzip_command = f'bgzip -f {vep_vcf}'
        vep_tabix_command = f'tabix -f -p vcf {vep_vcf}.gz'
        logger = getlogger('cpsr-vep')

        ## CPSR|VEP - run Variant Effect Predictor on query VCF with LoF and NearestExonJB plugins
        logger.info(f"CPSR - STEP 1: Basic variant annotation with Variant Effect Predictor ({VEP_VERSION}, GENCODE {GENCODE_VERSION}, {arg_dict['genome_assembly']})")
        logger.info(f"VEP configuration - one primary consequence block pr. alternative allele (--flag_pick_allele)")
        logger.info(f"VEP configuration - transcript pick order: {arg_dict['vep_pick_order']}")
        logger.info(f"VEP configuration - transcript pick order: See more at https://www.ensembl.org/info/docs/tools/vep/script/vep_other.html#pick_options")
        logger.info(f"VEP configuration - GENCODE set: {gencode_set_in_use}")
        logger.info(f"VEP configuration - skip intergenic: {arg_dict['vep_no_intergenic']}")
        logger.info(f"VEP configuration - look for overlap with regulatory regions: {vep_regulatory}")
        logger.info(f"VEP configuration - plugins in use: {plugins_in_use}")
        logger.info(f"VEP configuration - buffer_size/number of forks: {arg_dict['vep_buffer_size']}/{arg_dict['vep_n_forks']}")
        check_subprocess(logger, vep_main_command, debug)
        check_subprocess(logger, vep_bgzip_command, debug)
        check_subprocess(logger, vep_tabix_command, debug)
        logger.info("Finished cpsr-vep")
        print('----')

        ## CPSR|vcfanno - run vcfanno on query VCF with a number of relevant annotated VCFs
        logger = getlogger('cpsr-vcfanno')
        logger.info("CPSR - STEP 2: Annotation for cancer predisposition with cpsr-vcfanno")
        logger.info("(ClinVar, CIViC, dbNSFP, dbMTS, UniProtKB, cancerhotspots.org, ncER, GERP RS scores, GWAS catalog, gnomAD non-cancer subset)")
        pcgr_vcfanno_command = (
                f"pcgr_vcfanno.py --num_processes {arg_dict['vcfanno_n_proc']} --dbnsfp --clinvar "
                f"--cancer_hotspots --dbmts --ncer --gerp --civic --uniprot --gnomad_cpsr --pcgr_onco_xref "
                f"--gwas --rmsk {vep_vcf}.gz {vep_vcfanno_vcf} {os.path.join(data_dir, 'data', str(arg_dict['genome_assembly']))}"
                )
        check_subprocess(logger, pcgr_vcfanno_command, debug)
        logger.info("Finished cpsr-vcfanno")
        print('----')

        ## CPSR|summarise - expand annotations with separate VCF INFO tags
        logger = getlogger("cpsr-summarise")
        pcgr_summarise_command = (
                f'pcgr_summarise.py {vep_vcfanno_vcf}.gz 0 {vep_regulatory} '
                f'{os.path.join(data_dir, "data", arg_dict["genome_assembly"])} '
                f'--cpsr {"--debug" if debug else ""}'
                )
        logger.info("CPSR - STEP 3: Cancer gene annotations with cpsr-summarise")
        check_subprocess(logger, pcgr_summarise_command, debug)

        ## CPSR|clean - rename output files, remove temporary files
        os.rename(vep_vcfanno_annotated_vcf, output_vcf)
        os.rename(f'{vep_vcfanno_annotated_vcf}.tbi', f'{output_vcf}.tbi')
        os.rename(vep_vcfanno_annotated_pass_vcf, output_pass_vcf)
        os.rename(f'{vep_vcfanno_annotated_pass_vcf}.tbi', f'{output_pass_vcf}.tbi')
        delete_files = (
                glob(f'{vep_vcf}*') +
                glob(f'{vep_vcfanno_annotated_vcf}') +
                glob(f'{vep_vcfanno_annotated_pass_vcf}*') +
                glob(f'{vep_vcfanno_vcf}*') +
                glob(f'{input_vcf_cpsr_ready_uncompressed}*')
                )
        # do not delete if debugging
        if not debug:
            for fn in delete_files:
                #print(f"Deleting {fn}")
                utils.remove(fn)
        logger.info('Finished cpsr-summarise main command')
        ## CPSR|vcf2tsv - perform vcf2tsv conversion on the final annotated VCF file
        cpsr_vcf2tsv_command = f"vcf2tsv.py {output_pass_vcf} --compress {output_pass_tsv}"
        logger.info("Converting VCF to TSV with https://github.com/sigven/vcf2tsv")
        check_subprocess(logger, cpsr_vcf2tsv_command, debug)
        logger.info('Finished cpsr-summarise-vcf2tsv')
    logger.info('Finished cpsr-summarise')
    print('----')

    ## Generation of HTML reports for VEP/vcfanno-annotated VCF file
    if not arg_dict['basic']:
        logger = getlogger('cpsr-writer')
        logger.info("CPSR - STEP 4: Generation of output files - Cancer predisposition sequencing report")

        # export PATH to R conda env Rscript
        rscript = utils.script_path("pcgrr", "bin/Rscript")
        cpsrr_script = utils.script_path('pcgr', 'bin/cpsr.R')
        cpsr_report_command = (
                f"{rscript} {cpsrr_script} "
                f"{output_dir} "
                f"{output_pass_tsv}.gz "
                f"{arg_dict['sample_id']} "
                f"{pcgr_vars.PCGR_VERSION} "
                f"{pcgr_vars.DB_VERSION} "
                f"{arg_dict['genome_assembly']} "
                f"{data_dir} "
                f"{virtual_panel_id} "
                f"{preserved_info_tags} "
                f"{custom_bed} "
                f"{arg_dict['custom_list_name']} "
                f"{arg_dict['report_theme']} "
                f"{arg_dict['report_table_display']} "
                f"{report_nonfloating_toc} "
                f"{gwas_findings} "
                f"{arg_dict['gwas_p_value']} "
                f"{arg_dict['pop_gnomad']} "
                f"{arg_dict['maf_upper_threshold']} "
                f"{arg_dict['vep_pick_order']} "
                f"{arg_dict['vep_n_forks']} "
                f"{arg_dict['vep_buffer_size']} "
                f"{arg_dict['vep_gencode_all']} "
                f"{vep_no_intergenic} "
                f"{vep_regulatory} "
                f"{secondary_findings} "
                f"{classify_all} "
                f"{ignore_noncoding} "
                f"{clinvar_ignore_noncancer} "
                f"{diagnostic_grade_only}"
           )

        if debug:
            print(cpsr_report_command)
        check_subprocess(logger, cpsr_report_command, debug)
        logger.info("Finished CPSR!")
        print('----')
    print()


if __name__=="__main__":
    main()
