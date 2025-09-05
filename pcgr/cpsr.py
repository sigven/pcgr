#!/usr/bin/env python

import re
import argparse
import os
import yaml

from glob import glob
from argparse import RawTextHelpFormatter
from pcgr import pcgr_vars, arg_checker, utils, config, variant
from pcgr.utils import check_subprocess, getlogger, warn_message, remove_file, random_id_generator
from pcgr.config import populate_config_data
from pcgr.vep import get_vep_command


def get_args():

    program_description = "Cancer Predisposition Sequencing Reporter - report of " + \
       "clinically significant cancer-predisposing germline variants"
    program_options = "\n\t--input_vcf <INPUT_VCF>\n\t--vep_dir <VEP_DIR>\n\t--refdata_dir <REFDATA_DIR>\n\t" + \
        "--output_dir <OUTPUT_DIR>\n\t--genome_assembly <GENOME_ASSEMBLY>\n\t--sample_id <SAMPLE_ID>"

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
    optional_panel.add_argument('--diagnostic_grade_only', action="store_true",help="For panel_id's 1-44 (Genomics England PanelApp) - consider genes with a GREEN status only, default: %(default)s")

    optional_other.add_argument('--force_overwrite', action = "store_true", help='By default, the script will fail with an error if any output file already exists.\n You can force the overwrite of existing result files by using this flag, default: %(default)s')
    optional_other.add_argument('--version', action='version', version=str(utils.get_cpsr_version()))
    optional_other.add_argument('--no_reporting',action="store_true",help="Run functional variant annotation on VCF through VEP/vcfanno, omit classification/report generation (STEP 4), default: %(default)s")
    optional_other.add_argument("--no_html", action="store_true", help="Do not generate HTML report (default: %(default)s)")
    optional_other.add_argument('--retained_info_tags', dest ='retained_info_tags', default='None', help='Comma-separated string of VCF INFO tags from query VCF that should be kept in CPSR output TSV')
    optional_other.add_argument('--ignore_noncoding', action='store_true',dest='ignore_noncoding',default=False,help='Ignore non-coding (i.e. non protein-altering) variants in report, default: %(default)s')
    optional_other.add_argument("--debug", action="store_true", help="Print full commands to log")
    optional_other.add_argument("--pcgrr_conda", default="pcgrr", help="pcgrr conda env name (default: %(default)s)")
    
    optional_classification.add_argument('--secondary_findings', action='store_true',dest='secondary_findings',default=False, help='Include variants found in ACMG-recommended list for secondary findings (v3.2), default: %(default)s')
    optional_classification.add_argument('--pgx_findings', action='store_true',dest='pgx_findings',default=False, help='Report overlap with variants associated with chemotherapy toxicity (PgX findings, CPIC), default: %(default)s')
    optional_classification.add_argument('--gwas_findings', action='store_true',dest='gwas_findings',default=False, help='Report overlap with low to moderate cancer risk variants (tag SNPs) identified from genome-wide association studies, default: %(default)s')    
    optional_classification.add_argument('--pop_gnomad',choices = ['afr','amr','eas','sas','asj','nfe','fin','global'], default='nfe', help='Population source in gnomAD (non-cancer subset) used for variant frequency assessment (ACMG classification), default: %(default)s')
    optional_classification.add_argument('--maf_upper_threshold', type = float, default = 0.9, dest = 'maf_upper_threshold',help='Upper MAF limit (gnomAD global population frequency) for variants to be included in the report, default: %(default)s')
    optional_classification.add_argument('--classify_all', action='store_true',dest='classify_all',help='Provide CPSR variant classifications (TIER 1-5) also for variants with existing ClinVar classifications in output TSV, default: %(default)s')
    optional_classification.add_argument('--clinvar_report_noncancer', action='store_true', help='Report also ClinVar-classified variants attributed to phenotypes/conditions NOT directly related to tumor development, default: %(default)s')
    
    optional_vcfanno.add_argument('--vcfanno_n_proc', default = 4, type = int, help="Number of vcfanno processes (option '-p' in vcfanno), default: %(default)s")

    optional_vep.add_argument('--vep_n_forks', default = 4, type = int, help="Number of forks (option '--fork' in VEP), default: %(default)s")
    optional_vep.add_argument('--vep_buffer_size', default = 500, type = int, help="Variant buffer size (variants read into memory simultaneously, option '--buffer_size' in VEP) " + \
       "\n- set lower to reduce memory usage, default: %(default)s")
    optional_vep.add_argument("--vep_gencode_basic", action="store_true", help = "Consider basic GENCODE transcript set only with Variant Effect Predictor (VEP) (option '--gencode_basic' in VEP).")
    optional_vep.add_argument('--vep_pick_order', default = "mane_select,mane_plus_clinical,canonical,biotype,ccds,rank,tsl,appris,length", help="Comma-separated string " + \
       "of ordered transcript properties for primary variant pick\n ( option '--pick_order' in VEP), default: %(default)s")
    optional_vep.add_argument('--vep_no_intergenic', action = "store_true", help="Skip intergenic variants during processing (option '--no_intergenic' in VEP), default: %(default)s")

    required.add_argument('--input_vcf', help='VCF input file with germline query variants (SNVs/InDels).', required = True)
    required.add_argument("--vep_dir", dest="vep_dir", help="Directory of VEP cache, e.g.  $HOME/.vep", required=True)
    required.add_argument('--refdata_dir',help=f"Directory that contains the PCGR/CPSR reference data, e.g. ~/pcgr-data-{pcgr_vars.PCGR_VERSION}", required = True)
    required.add_argument('--output_dir',help='Output directory', required = True)
    required.add_argument('--genome_assembly',choices = ['grch37','grch38'], help='Genome assembly build: grch37 or grch38', required = True)
    required.add_argument('--sample_id',help="Sample identifier - prefix for output files", required = True)

    args = parser.parse_args()
    return vars(args)

def main():
    arg_dict = get_args()
    # check parsed arguments
    arg_dict_verified = arg_checker.verify_args_cpsr(arg_dict)
    # Verify existence of input files
    
    input_data = arg_checker.verify_input_files_cpsr(arg_dict_verified)
    output_data = arg_checker.define_output_files(arg_dict, cpsr = True)
    
    # create configuration/options
    conf_options = config.create_config(arg_dict, workflow = "CPSR")
    
    ## Run CPSR workflow
    run_cpsr(conf_options, input_data, output_data)


def run_cpsr(conf_options, input_data, output_data):
    """
    Main function to run the CPSR workflow
    """
    
    debug = conf_options['debug']
    vep_skip_intergenic_set = 'ON' if conf_options['vep']['vep_no_intergenic'] == 1 else 'OFF'    
    uid = ''
    genome_assembly = str(conf_options['genome_assembly'])
    input_vcf = 'None'
    input_customlist = 'None'
    output_custom_bed = 'None'
    output_dir = output_data['dir']
    output_prefix = output_data['prefix']
    print('')    
    logger = getlogger('cpsr-validate-input-arguments')
    logger.info("CPSR - STEP 0: Validate input data")

    if input_data['input_vcf_basename'] != 'NA':
        input_vcf = os.path.join(input_data['input_vcf_dir'], input_data['input_vcf_basename'])
    if input_data['input_customlist_basename'] != 'NA' and input_data['input_customlist_dir'] != 'NA':
        input_customlist = os.path.join(
            input_data['input_customlist_dir'], input_data['input_customlist_basename'])
        output_custom_bed =  f'{output_prefix}.custom_panel.bed'
        if conf_options['gene_panel']['custom_list_name'] == "None":
            warn_msg = "No custom list name provided, use argument '--custom_list_name' to provide name for custom defined panel"
            warn_message(warn_msg, logger)

    if not input_vcf == 'None':
        check_subprocess(logger, f'mkdir -p {output_dir}', debug)
        
        ## Define input, output and temporary file names
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
        yaml_fname =                      f'{output_prefix}.conf.yaml'

        ## CPSR|Validate input VCF - check formatting, non-overlap with CPSR INFO tags, and whether sample contains any variants in cancer predisposition loci
        vcf_validate_command = (
                f'cpsr_validate_input.py '
                f'{input_data["refdata_assembly_dir"]} '
                f'{input_vcf} '
                f'{input_vcf_validated_uncompr} '
                f'{input_customlist} '
                f'{output_custom_bed} '
                f'{conf_options["other"]["retained_vcf_info_tags"]} '
                f'{conf_options["sample_id"]} '
                f'{conf_options["gene_panel"]["panel_id"]} '
                f'{conf_options["gene_panel"]["diagnostic_grade_only"]} '
                f'{conf_options["variant_classification"]["gwas_findings"]} '
                f'{conf_options["variant_classification"]["secondary_findings"]} '
                f'{conf_options["variant_classification"]["pgx_findings"]} '
                f'--output_dir {output_dir} {"--debug" if debug else ""}'
                )
        check_subprocess(logger, vcf_validate_command, debug)
    
        logger.info('Finished cpsr-validate-input-arguments')
        print('----')

        ## CPSR|Start - log key information about run
        logger = getlogger("cpsr-settings")
        logger.info("--- Cancer Predisposition Sequencing Reporter workflow ----")
        logger.info(f"Sample name: {conf_options['sample_id']}")
        if not input_customlist == 'None':
            logger.info(f"Virtual gene panel: custom-made list from panel 0: {input_customlist}")
        else:
            logger.info("Diagnostic-grade genes in virtual panels (GE PanelApp): " + \
                        f"{'ON' if conf_options['gene_panel']['diagnostic_grade_only'] else 'OFF'}")
        logger.info("Include incidental findings (ACMG recommended list v3.2): " + \
                    f"{'ON' if conf_options['variant_classification']['secondary_findings'] else 'OFF'}")
        logger.info("Include low to moderate cancer risk variants from genome-wide association studies: " + \
                    f"{'ON' if conf_options['variant_classification']['gwas_findings'] else 'OFF'}")
        logger.info("Include pharmacogenetic findings (PgX - variants related to potential toxicity to chemotherapy): " + \
                    f"{'ON' if conf_options['variant_classification']['pgx_findings'] else 'OFF'}")
        logger.info("Reference population, germline variant frequencies (gnomAD - non-cancer subset): " + \
                    f"{str(conf_options['variant_classification']['pop_gnomad']).upper()}")
        logger.info(f"Genome assembly: {conf_options['genome_assembly']}")
        print('----')
        
        conf_options['molecular_data']['fname_mut_vcf'] = output_vcf
        conf_options['molecular_data']['fname_mut_tsv'] = output_pass_tsv_gz
        conf_options['output_dir'] = output_dir
        conf_options['output_prefix'] = output_prefix
        conf_options['gene_panel']['custom_list_bed'] = "None"
        if not input_customlist == 'None':
            conf_options['gene_panel']['custom_list_bed'] = output_custom_bed
        
        ## Write YAML configuration file  - settings, path to files, reference bundle etc
        yaml_data = populate_config_data(conf_options, input_data["refdata_assembly_dir"], workflow = "CPSR", logger = logger)
        with open(yaml_fname, "w") as outfile:
            outfile.write(yaml.dump(yaml_data))
        outfile.close()
        
        ## CPSR|VEP - run Variant Effect Predictor on query VCF with NearestExonJB plugins
        vep_command = get_vep_command(file_paths = input_data, 
                                      conf_options = yaml_data, 
                                      input_vcf = input_vcf_validated, 
                                      output_vcf = vep_vcf)

        logger = getlogger('cpsr-vep')

        logger.info((
            f"CPSR - STEP 1: Basic variant annotation with Variant Effect Predictor (version {pcgr_vars.VEP_VERSION}, "
            f"GENCODE release {pcgr_vars.GENCODE_VERSION[genome_assembly]}, genome assembly {yaml_data['genome_assembly']})"))
        logger.info("VEP configuration - one primary consequence block pr. alternative allele (--flag_pick_allele_gene)")
        logger.info(f"VEP configuration - transcript pick order: {yaml_data['conf']['vep']['vep_pick_order']}")
        logger.info("VEP configuration - transcript pick order: See more at https://www.ensembl.org/info/docs/tools/vep/script/vep_other.html#pick_options")
        logger.info(f"VEP configuration - GENCODE set: {vep_command['gencode_set_in_use']}")
        logger.info(f'VEP configuration - skip intergenic variants: {vep_skip_intergenic_set}')
        logger.info("VEP configuration - look for overlap with regulatory regions: ON")
        logger.info(f"VEP configuration - plugins in use: {vep_command['plugins_in_use']}")
        logger.info((
            f"VEP configuration - buffer size/number of forks: "
            f"{yaml_data['conf']['vep']['vep_buffer_size']}/{yaml_data['conf']['vep']['vep_n_forks']}"))
        check_subprocess(logger, vep_command["main"], debug)
        check_subprocess(logger, vep_command["bgzip"], debug)
        check_subprocess(logger, vep_command["tabix"], debug)
        logger.info("Finished cpsr-vep")
        print('----')

        ## CPSR|vcfanno - run vcfanno on query VCF with a number of relevant annotated VCFs
        logger = getlogger('cpsr-vcfanno')
        logger.info("CPSR - STEP 2: Annotation against BED/VCF tracks with cpsr-vcfanno")
        logger.info("(ClinVar, CIViC, dbNSFP, dbMTS, GERP, GWAS catalog, gnomAD non-cancer subset)")
        pcgr_vcfanno_command = (
                f"pcgr_vcfanno.py --num_processes {yaml_data['conf']['other']['vcfanno_n_proc']} --dbnsfp --clinvar "
                f'{"--debug" if debug else ""} '
                f"--dbmts --gerp --tcga --gnomad_non_cancer --gene_transcript_xref "
                f"--gwas --rmsk {vep_vcf}.gz {vep_vcfanno_vcf} "
                f"{input_data['refdata_assembly_dir']}"
                )
        check_subprocess(logger, pcgr_vcfanno_command, debug)
        logger.info("Finished cpsr-vcfanno")
        print('----')
        #exit(0)

        ## CPSR|summarise - expand annotations with separate VCF INFO tags
        logger = getlogger("cpsr-summarise")
        cpsr_summarise_command = (
                f'pcgr_summarise.py {vep_vcfanno_vcf}.gz {vep_vcfanno_summarised_vcf} 0 '
                f'{yaml_data["conf"]["vep"]["vep_regulatory"]} 0 '
                f'Any {yaml_data["genome_assembly"]} {yaml_data["conf"]["vep"]["vep_pick_order"]} '
                f'{input_data["refdata_assembly_dir"]} --compress_output_vcf '
                f'--cpsr_yaml {yaml_fname} '
                f'--cpsr {"--debug" if debug else ""}'
                )
        logger.info("CPSR - STEP 3: Cancer gene annotations with cpsr-summarise")
        check_subprocess(logger, cpsr_summarise_command, debug)

        ## CPSR|clean - rename output files, remove temporary files
        os.rename(f'{vep_vcfanno_summarised_vcf}.gz', output_vcf)
        os.rename(f'{vep_vcfanno_summarised_vcf}.gz.tbi', f'{output_vcf}.tbi')
        os.rename(f'{vep_vcfanno_summarised_pass_vcf}.gz', output_pass_vcf)
        os.rename(f'{vep_vcfanno_summarised_pass_vcf}.gz.tbi', f'{output_pass_vcf}.tbi')
        delete_files = (
                glob(f'{vep_vcf}*') +
                glob(f'{vep_vcfanno_summarised_vcf}*') +
                glob(f'{vep_vcfanno_summarised_pass_vcf}*') +
                glob(f'{vep_vcfanno_vcf}*') +
                glob(f'{input_vcf_validated_uncompr}*')
                )
        # do not delete if debugging
        if not debug:
            for fn in delete_files:
                remove_file(fn)
        logger.info('Finished cpsr-summarise main command')
        
        # CPSR|vcf2tsvpy - convert VCF to TSV with https://github.com/sigven/vcf2tsvpy
        cpsr_vcf2tsv_command = f'vcf2tsvpy --input_vcf {output_pass_vcf} --out_tsv {output_pass_vcf2tsv} --compress'
        logger.info("Converting VCF to TSV with https://github.com/sigven/vcf2tsvpy")        
        check_subprocess(logger, cpsr_vcf2tsv_command, debug) 
        logger.info('Finished cpsr-summarise-vcf2tsv')       
        logger.info("Appending ClinVar traits, official gene names, and protein domain annotations")        
        variant_set = \
           variant.append_annotations(
              output_pass_vcf2tsv_gz, refdata_assembly_dir = input_data["refdata_assembly_dir"], logger = logger)
        variant_set = variant.clean_annotations(variant_set, yaml_data, logger = logger)
        
        ## If no genotypes are available, set conf['sample_properties']['genotypes_available'] = 1
        if {'GENOTYPE'}.issubset(variant_set.columns):
            if variant_set.loc[variant_set['GENOTYPE'] == '.'].empty and variant_set.loc[variant_set['GENOTYPE'] == 'undefined'].empty:
                yaml_data['conf']['sample_properties']['gt_detected'] = 1
        if {'DP_CONTROL'}.issubset(variant_set.columns):
            if variant_set.loc[variant_set['DP_CONTROL'] == '-1'].empty:
                yaml_data['conf']['sample_properties']['dp_detected'] = 1
        
        yaml_data['conf']['gene_panel']['panel_id'] = re.sub(r',',';', yaml_data['conf']['gene_panel']['panel_id'])
        with open(yaml_fname, "w") as outfile:
            outfile.write(yaml.dump(yaml_data))
        outfile.close()
        
        variant_set.fillna('.').to_csv(output_pass_tsv_gz, sep="\t", compression="gzip", index=False)
        if not debug:
            remove_file(output_pass_vcf2tsv_gz)
                
        logger.info('Finished cpsr-summarise')
                
    print('----')

    ## Generation of HTML reports for VEP/vcfanno-annotated VCF file
    if not conf_options['other']['no_reporting']:
        logger = getlogger('cpsr-writer')
        logger.info("CPSR - STEP 4: Generation of output files - Cancer predisposition sequencing report")
        
        # export PATH to R conda env Rscript
        pcgrr_conda = conf_options['pcgrr_conda']
        quarto_env_vars = utils.quarto_evars_path(pcgrr_conda)
        pcgr_conda = utils.conda_prefix_basename()
        rscript = utils.script_path(pcgrr_conda, 'bin/Rscript')
        cpsrr_script = utils.script_path(pcgr_conda, 'bin/cpsr.R')
        export_pcgrr = utils.pcgrr_conda_env_export(pcgrr_conda)
        ## CPSR|writer - generate HTML report
        cpsr_report_command = (
                 f"{export_pcgrr} && {rscript} {cpsrr_script} {yaml_fname} {quarto_env_vars}")

        if debug:
            print(cpsr_report_command)
        check_subprocess(logger, cpsr_report_command, debug)
        logger.info("Finished CPSR!")
        print('----')
    print()


if __name__=="__main__":
    main()
