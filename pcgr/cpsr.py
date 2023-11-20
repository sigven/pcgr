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
import json
import yaml

from glob import glob
from argparse import RawTextHelpFormatter
from pcgr import pcgr_vars, arg_checker, utils, vep, config, variant
from pcgr.utils import check_subprocess, getlogger, error_message, warn_message
from pcgr.config import populate_config_data

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
    optional_classification = parser.add_argument_group("Variant classification options")
    optional_vep = parser.add_argument_group('VEP options')
    optional_vcfanno = parser.add_argument_group('vcfanno options')
    optional_other = parser.add_argument_group('Other options')

    optional_panel.add_argument('--panel_id',dest = "virtual_panel_id",type = str, default = "-1", help="Comma-separated string with identifier(s) of predefined virtual cancer predisposition gene panels,\nchoose any combination of the following identifiers (GEP = Genomics England PanelApp):\n" + str(pcgr_vars.panels))
    optional_panel.add_argument('--custom_list',dest = "custom_list",help="Provide custom list of genes from virtual panel 0 (single-column .txt/.tsv file with Ensembl gene identifiers),\n alternative to predefined panels provided with --panel_id)")
    optional_panel.add_argument('--custom_list_name',dest = "custom_list_name", default="None", help="Set name for custom made panel/list (single word - no whitespace), will be displayed in the report")
    optional_panel.add_argument('--diagnostic_grade_only', action="store_true",help="For panel_id's 1-42 (Genomics England PanelApp) - consider genes with a GREEN status only, default: %(default)s")

    optional_other.add_argument('--force_overwrite', action = "store_true", help='By default, the script will fail with an error if any output file already exists.\n You can force the overwrite of existing result files by using this flag, default: %(default)s')
    optional_other.add_argument('--version', action='version', version=str(utils.get_cpsr_version()))
    optional_other.add_argument('--no_reporting',action="store_true",help="Run functional variant annotation on VCF through VEP/vcfanno, omit classification/report generation (STEP 4), default: %(default)s")
    optional_other.add_argument('--retained_info_tags', dest ='retained_info_tags', default='None', help='Comma-separated string of VCF INFO tags from query VCF that should be kept in CPSR output TSV')
    optional_other.add_argument('--report_theme',choices = ['default','cerulean','journal','flatly','readable','spacelab','united','cosmo','lumen','paper','sandstone','simplex','yeti'], default = 'default', help='Visual report theme (rmarkdown),  default: %(default)s' )
    optional_other.add_argument('--report_nonfloating_toc', action='store_true', help='Do not float the table of contents (TOC) in output HTML report, default: %(default)s')
    optional_other.add_argument('--report_table_display', choices = ['full','light'], default='light', help="Set the level of detail/comprehensiveness in interactive datables of HTML report, very comprehensive (option 'full') or slim/focused ('light'), default: %(default)s")
    optional_other.add_argument('--ignore_noncoding', action='store_true',dest='ignore_noncoding',default=False,help='Do not list non-coding variants in HTML report, default: %(default)s')
    optional_other.add_argument("--debug", action="store_true", help="Print full commands to log")
    optional_other.add_argument("--pcgrr_conda", default="pcgrr", help="pcgrr conda env name (default: %(default)s)")
    
    optional_classification.add_argument('--secondary_findings', action='store_true',dest='secondary_findings',default=False, help='Include variants found in ACMG-recommended list for secondary findings (v3.2), default: %(default)s')
    optional_classification.add_argument('--gwas_findings', action='store_true',dest='gwas_findings',default=False, help='Report overlap with low to moderate cancer risk variants (tag SNPs) identified from genome-wide association studies, default: %(default)s')    
    optional_classification.add_argument('--pop_gnomad',choices = ['afr','amr','eas','sas','asj','nfe','fin','global'], default='nfe', help='Population source in gnomAD (non-cancer subset) used for variant frequency assessment (ACMG classification), default: %(default)s')
    optional_classification.add_argument('--maf_upper_threshold', type = float, default = 0.9, dest = 'maf_upper_threshold',help='Upper MAF limit (gnomAD global population frequency) for variants to be included in the report, default: %(default)s')
    optional_classification.add_argument('--classify_all', action='store_true',dest='classify_all',help='Provide CPSR variant classifications (TIER 1-5) also for variants with existing ClinVar classifications in output TSV, default: %(default)s')
    optional_classification.add_argument('--clinvar_ignore_noncancer', action='store_true', help='Ignore (exclude from report) ClinVar-classified variants reported only for phenotypes/conditions NOT related to cancer, default: %(default)s')
    
    optional_vcfanno.add_argument('--vcfanno_n_proc', default = 4, type = int, help="Number of vcfanno processes (option '-p' in vcfanno), default: %(default)s")

    optional_vep.add_argument('--vep_n_forks', default = 4, type = int, help="Number of forks (option '--fork' in VEP), default: %(default)s")
    optional_vep.add_argument('--vep_buffer_size', default = 500, type = int, help="Variant buffer size (variants read into memory simultaneously, option '--buffer_size' in VEP) " + \
       "\n- set lower to reduce memory usage, default: %(default)s")
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
    arg_dict['show_noncoding'] = int(not arg_dict['ignore_noncoding'])
    # check parsed arguments
    arg_checker.check_args_cpsr(arg_dict)
    # Verify existence of input files
    cpsr_paths = arg_checker.verify_input_files_cpsr(arg_dict)
    
    # create configuration/options
    conf_options = config.create_config(arg_dict, workflow = "CPSR")
    
    ## Run CPSR workflow
    run_cpsr(conf_options, cpsr_paths)


def run_cpsr(conf_options, cpsr_paths):
    """
    Main function to run the CPSR workflow
    """
    
    prefix = os.path.join(
        cpsr_paths['output_dir'], f'{conf_options["sample_id"]}.cpsr.{conf_options["genome_assembly"]}')

    debug = conf_options['debug']
    vep_skip_intergenic_set = 'ON' if conf_options['vep']['vep_no_intergenic'] == 1 else 'OFF'
    output_vcf = 'None'
    output_pass_vcf = 'None'
    output_pass_tsv = 'None'
    uid = ''
    GENCODE_VERSION = pcgr_vars.GENCODE_VERSION
    VEP_VERSION = pcgr_vars.VEP_VERSION
    if conf_options['genome_assembly'] == 'grch37':
        GENCODE_VERSION = '19'

    input_vcf = 'None'
    input_customlist = 'None'
    output_custom_bed = 'None'
    print('')    
    logger = getlogger('cpsr-validate-input-arguments')
    logger.info("CPSR - STEP 0: Validate input data")

    if cpsr_paths['input_vcf_basename'] != 'NA':
        input_vcf = os.path.join(cpsr_paths['input_vcf_dir'], cpsr_paths['input_vcf_basename'])
    if cpsr_paths['input_customlist_basename'] != 'NA' and cpsr_paths['input_customlist_dir'] != 'NA':
        input_customlist = os.path.join(
            cpsr_paths['input_customlist_dir'], cpsr_paths['input_customlist_basename'])
        output_custom_bed =      f'{prefix}.custom_panel.bed'
        if conf_options['gene_panel']['custom_list_name'] == "None":
            warn_msg = "No custom list name provided, use argument '--custom_list_name' to provide name for custom defined panel"
            warn_message(warn_msg, logger)

    data_dir = cpsr_paths['base_dir']
    output_dir = cpsr_paths['output_dir']

    check_subprocess(logger, f'mkdir -p {output_dir}', debug)

    ## CPSR|Validate input VCF - check formatting, non-overlap with CPSR INFO tags, and whether sample contains any variants in cancer predisposition loci
    vcf_validate_command = (
            f'cpsr_validate_input.py '
            f'{data_dir} '
            f'{input_vcf} '
            f'{input_customlist} '
            f'{output_custom_bed} '
            f'{conf_options["other"]["retained_vcf_info_tags"]} '
            f'{conf_options["genome_assembly"]} '
            f'{conf_options["sample_id"]} '
            f'{conf_options["gene_panel"]["panel_id"]} '
            f'{conf_options["gene_panel"]["diagnostic_grade_only"]} '
            f'{conf_options["variant_classification"]["gwas_findings"]} '
            f'{conf_options["variant_classification"]["secondary_findings"]} '
            f'--output_dir {output_dir} {"--debug" if debug else ""}'
            )
    check_subprocess(logger, vcf_validate_command, debug)
 
    logger.info('Finished cpsr-validate-input-arguments')
    print('----')

    ## CPSR|Start - log key information about run
    logger = getlogger("cpsr-start")
    logger.info("--- Cancer Predisposition Sequencing Reporter workflow ----")
    logger.info(f"Sample name: {conf_options['sample_id']}")
    if not input_customlist == 'None':
        logger.info(f"Virtual gene panel: custom-made list from panel 0: {input_customlist}")
    else:
        logger.info(f"Diagnostic-grade genes in virtual panels (GE PanelApp): " + \
                    f"{'ON' if conf_options['gene_panel']['diagnostic_grade_only'] else 'OFF'}")
    logger.info(f"Include incidental findings (ACMG recommended list v3.1): " + \
                f"{'ON' if conf_options['variant_classification']['secondary_findings'] else 'OFF'}")
    logger.info(f"Include low to moderate cancer risk variants from genome-wide association studies: " + \
                f"{'ON' if conf_options['variant_classification']['gwas_findings'] else 'OFF'}")
    logger.info(f"Reference population, germline variant frequencies (gnomAD): {str(conf_options['variant_classification']['pop_gnomad']).upper()}")
    logger.info(f"Genome assembly: {conf_options['genome_assembly']}")
    print('----')


    if not input_vcf == 'None':

        ## Define input, output and temporary file names        
        output_vcf =             f'{prefix}.vcf.gz'
        output_pass_vcf =        f'{prefix}.pass.vcf.gz'
        output_pass_vcf2tsv =    f'{prefix}.pass.vcf2tsv.tsv'
        output_pass_vcf2tsv_gz = f'{output_pass_vcf2tsv}.gz'
        output_pass_tsv =        f'{prefix}.pass.tsv'
        output_pass_tsv_gz =     f'{output_pass_tsv}.gz'
        yaml_fname =             f'{prefix}.conf.yaml'
        
        conf_options['annotated_tsv'] = output_pass_tsv_gz
        conf_options['annotated_vcf'] = output_vcf
        conf_options['output_dir'] = output_dir
        conf_options['gene_panel']['custom_list_bed'] = "None"
        if not input_customlist == 'None':
            conf_options['gene_panel']['custom_list_bed'] = output_custom_bed
        
        yaml_data = populate_config_data(conf_options, data_dir, workflow = "CPSR", logger = logger)
        
        with open(yaml_fname, "w") as outfile:
            outfile.write(yaml.dump(yaml_data))
        outfile.close()
        
        input_vcf_cpsr_ready =   os.path.join(
            output_dir, 
            re.sub(r"(\.vcf$|\.vcf\.gz$)", 
                   ".cpsr_ready_target.vcf.gz", 
                   cpsr_paths["input_vcf_basename"]))
        input_vcf_cpsr_ready_uncompressed = os.path.join(
            output_dir, 
            re.sub(r"(\.vcf$|\.vcf\.gz$)", 
                    ".cpsr_ready_target.vcf", 
                    cpsr_paths["input_vcf_basename"]))
        vep_vcf = re.sub(r"(\.vcf$|\.vcf\.gz$)", ".cpsr_vep.vcf", input_vcf_cpsr_ready)
        vep_vcfanno_vcf = re.sub(r"(\.vcf$|\.vcf\.gz$)", ".cpsr_vep.vcfanno.vcf", input_vcf_cpsr_ready)
        vep_vcfanno_summarized_vcf = re.sub(r"\.vcfanno", ".vcfanno.summarized", vep_vcfanno_vcf)
        vep_vcfanno_summarized_pass_vcf = re.sub(r"\.vcfanno", ".vcfanno.summarized.pass", vep_vcfanno_vcf)

        ## CPSR|VEP - run Variant Effect Predictor on query VCF with LoF and NearestExonJB plugins
        vep_command = vep.get_command(file_paths = cpsr_paths, 
                                      conf_options = yaml_data, 
                                      input_vcf = input_vcf_cpsr_ready, 
                                      output_vcf = vep_vcf)

        logger = getlogger('cpsr-vep')

        logger.info(f"CPSR - STEP 1: Basic variant annotation with Variant Effect Predictor ({VEP_VERSION}, GENCODE {GENCODE_VERSION}, {yaml_data['genome_assembly']})")
        logger.info(f"VEP configuration - one primary consequence block pr. alternative allele (--flag_pick_allele_gene)")
        logger.info(f"VEP configuration - transcript pick order: {yaml_data['conf']['vep']['vep_pick_order']}")
        logger.info(f"VEP configuration - transcript pick order: See more at https://www.ensembl.org/info/docs/tools/vep/script/vep_other.html#pick_options")
        logger.info(f"VEP configuration - GENCODE set: {vep_command['gencode_set_in_use']}")
        logger.info(f'VEP configuration - skip intergenic variants: {vep_skip_intergenic_set}')
        logger.info(f"VEP configuration - look for overlap with regulatory regions: ON")
        logger.info(f"VEP configuration - plugins in use: {vep_command['plugins_in_use']}")
        logger.info(f"VEP configuration - buffer_size/number of forks: {yaml_data['conf']['vep']['vep_buffer_size']}/{yaml_data['conf']['vep']['vep_n_forks']}")
        check_subprocess(logger, vep_command["main"], debug)
        check_subprocess(logger, vep_command["bgzip"], debug)
        check_subprocess(logger, vep_command["tabix"], debug)
        logger.info("Finished cpsr-vep")
        print('----')

        ## CPSR|vcfanno - run vcfanno on query VCF with a number of relevant annotated VCFs
        logger = getlogger('cpsr-vcfanno')
        logger.info("CPSR - STEP 2: Annotation with BED/VCF tracks with cpsr-vcfanno")
        logger.info("(ClinVar, CIViC, dbNSFP, dbMTS, GERP, GWAS catalog, gnomAD non-cancer subset)")
        pcgr_vcfanno_command = (
                f"pcgr_vcfanno.py --num_processes {yaml_data['conf']['other']['vcfanno_n_proc']} --dbnsfp --clinvar "
                f"--dbmts --gerp --tcga --gnomad_non_cancer --gene_transcript_xref "
                f"--gwas --rmsk {vep_vcf}.gz {vep_vcfanno_vcf} {os.path.join(data_dir, 'data', str(yaml_data['genome_assembly']))}"
                )
        check_subprocess(logger, pcgr_vcfanno_command, debug)
        logger.info("Finished cpsr-vcfanno")
        print('----')

        ## CPSR|summarise - expand annotations with separate VCF INFO tags
        logger = getlogger("cpsr-summarise")
        cpsr_summarise_command = (
                f'pcgr_summarise.py {vep_vcfanno_vcf}.gz {vep_vcfanno_summarized_vcf} 0 '
                f'{yaml_data["conf"]["vep"]["vep_regulatory"]} 0 '
                f'Any {yaml_data["conf"]["vep"]["vep_pick_order"]} '
                f'{cpsr_paths["db_dir"]} --compress_output_vcf '
                f'--cpsr_yaml {yaml_fname} '
                f'--cpsr {"--debug" if debug else ""}'
                )
        logger.info("CPSR - STEP 3: Cancer gene annotations with cpsr-summarise")
        check_subprocess(logger, cpsr_summarise_command, debug)

        ## CPSR|clean - rename output files, remove temporary files
        os.rename(f'{vep_vcfanno_summarized_vcf}.gz', output_vcf)
        os.rename(f'{vep_vcfanno_summarized_vcf}.gz.tbi', f'{output_vcf}.tbi')
        os.rename(f'{vep_vcfanno_summarized_pass_vcf}.gz', output_pass_vcf)
        os.rename(f'{vep_vcfanno_summarized_pass_vcf}.gz.tbi', f'{output_pass_vcf}.tbi')
        delete_files = (
                glob(f'{vep_vcf}*') +
                glob(f'{vep_vcfanno_summarized_vcf}*') +
                glob(f'{vep_vcfanno_summarized_pass_vcf}*') +
                glob(f'{vep_vcfanno_vcf}*') +
                glob(f'{input_vcf_cpsr_ready_uncompressed}*')
                )
        # do not delete if debugging
        if not debug:
            for fn in delete_files:
                utils.remove(fn)
        logger.info('Finished cpsr-summarise main command')
        
        # CPSR|vcf2tsvpy - convert VCF to TSV with https://github.com/sigven/vcf2tsvpy
        cpsr_vcf2tsv_command = f'vcf2tsvpy --input_vcf {output_pass_vcf} --out_tsv {output_pass_vcf2tsv} --compress'
        logger.info("Converting VCF to TSV with https://github.com/sigven/vcf2tsvpy")        
        check_subprocess(logger, cpsr_vcf2tsv_command, debug) 
        logger.info('Finished cpsr-summarise-vcf2tsv')       
        logger.info("Appending ClinVar and protein domain (Pfam) annotations")        
        variant_set = \
           variant.append_annotations(
              output_pass_vcf2tsv_gz, pcgr_db_dir = cpsr_paths["db_dir"], logger = logger)
        variant_set.rename(columns = {'Consequence':'CONSEQUENCE'}, inplace = True)
        variant_set['EFFECT_PREDICTIONS'] = variant_set['EFFECT_PREDICTIONS'].str.replace("\\.&|\\.$", "NA&", regex = True)
        variant_set['EFFECT_PREDICTIONS'] = variant_set['EFFECT_PREDICTIONS'].str.replace("&$", "", regex = True)
        variant_set['EFFECT_PREDICTIONS'] = variant_set['EFFECT_PREDICTIONS'].str.replace("&", ", ", regex = True)
        variant_set.loc[variant_set['CLINVAR_CONFLICTED'] == 1, "CLINVAR_CONFLICTED"] = True
        variant_set.loc[variant_set['CLINVAR_CONFLICTED'] != 1, "CLINVAR_CONFLICTED"] = False
        variant_set['GENOMIC_CHANGE'] = variant_set['CHROM'].astype(str) + ":g." + variant_set['POS'].astype(str) + \
            variant_set['REF'].astype(str) + ">" + variant_set['ALT'].astype(str)
        variant_set['GENOME_VERSION'] = yaml_data['genome_assembly']        
        variant_set.to_csv(output_pass_tsv_gz, sep="\t", compression="gzip", index=False)
        utils.remove(output_pass_vcf2tsv_gz)
                
        logger.info('Finished cpsr-summarise')
        
        # CPSR|Generate YAML file - containing configuration options and paths to annotated molecular profile dataset
        # - VCF/TSV files (SNVs/InDels)       
        logger = getlogger('pcgr-write-yaml')
        
    print('----')

    ## Generation of HTML reports for VEP/vcfanno-annotated VCF file
    if not conf_options['other']['no_reporting']:
        logger = getlogger('cpsr-writer')
        logger.info("CPSR - STEP 4: Generation of output files - Cancer predisposition sequencing report")
        
        # export PATH to R conda env Rscript
        pcgrr_conda = conf_options['pcgrr_conda']
        pcgr_conda = utils.conda_prefix_basename()
        rscript = utils.script_path(pcgrr_conda, 'bin/Rscript')
        cpsrr_script = utils.script_path(pcgr_conda, 'bin/cpsr.R')
        
        cpsr_report_command = (
                 f"{rscript} {cpsrr_script} {yaml_fname}")

        if debug:
            print(cpsr_report_command)
        check_subprocess(logger, cpsr_report_command, debug)
        logger.info("Finished CPSR!")
        print('----')
    print()


if __name__=="__main__":
    main()
