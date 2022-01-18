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
#import toml
from argparse import RawTextHelpFormatter

PCGR_VERSION = '0.9.2'
CPSR_VERSION = '0.6.2'
DB_VERSION = 'PCGR_DB_VERSION = 20210627'
VEP_VERSION = '104'
GENCODE_VERSION = '38'
VEP_ASSEMBLY = 'GRCh38'
DOCKER_IMAGE_VERSION = 'sigven/pcgr:' + str(PCGR_VERSION)

global debug
#global VEP_ASSEMBLY

GE_panels = {
		0: "CPSR exploratory cancer predisposition panel (n = 433, Genomics England PanelApp / TCGA Germline Study / Cancer Gene Census / Other)",
      1: "Adult solid tumours cancer susceptibility (Genomics England PanelApp)",
      2: "Adult solid tumours for rare disease (Genomics England PanelApp)",
      3: "Bladder cancer pertinent cancer susceptibility (Genomics England PanelApp)",
      4: "Brain cancer pertinent cancer susceptibility (Genomics England PanelApp)",
      5: "Breast cancer pertinent cancer susceptibility (Genomics England PanelApp)",
      6: "Childhood solid tumours cancer susceptibility (Genomics England PanelApp)",
      7: "Colorectal cancer pertinent cancer susceptibility (Genomics England PanelApp)",
      8: "Endometrial cancer pertinent cancer susceptibility (Genomics England PanelApp)",
      9: "Familial Tumours Syndromes of the central & peripheral Nervous system (Genomics England PanelApp)",
      10: "Familial breast cancer (Genomics England PanelApp)",
      11: "Familial melanoma (Genomics England PanelApp)",
      12: "Familial prostate cancer (Genomics England PanelApp)",
      13: "Familial rhabdomyosarcoma (Genomics England PanelApp)",
      14: "GI tract tumours (Genomics England PanelApp)",
      15: "Genodermatoses with malignancies (Genomics England PanelApp)",
      16: "Haematological malignancies cancer susceptibility (Genomics England PanelApp)",
      17: "Haematological malignancies for rare disease (Genomics England PanelApp)",
      18: "Head and neck cancer pertinent cancer susceptibility (Genomics England PanelApp)",
      19: "Inherited MMR deficiency (Lynch Syndrome) - Genomics England PanelApp",
      20: "Inherited non-medullary thyroid cancer (Genomics England PanelApp)",
      21: "Inherited ovarian cancer (without breast cancer) (Genomics England PanelApp)",
      22: "Inherited pancreatic cancer (Genomics England PanelApp)",
      23: "Inherited polyposis (Genomics England PanelApp)",
      24: "Inherited predisposition to acute myeloid leukaemia (AML) - Genomics England PanelApp",
      25: "Inherited predisposition to GIST (Genomics England PanelApp)",
      26: "Inherited renal cancer (Genomics England PanelApp)",
      27: "Inherited phaeochromocytoma and paraganglioma (Genomics England PanelApp)",
      28: "Melanoma pertinent cancer susceptibility (Genomics England PanelApp)",
      29: "Multiple endocrine tumours (Genomics England PanelApp)",
      30: "Multiple monogenic benign skin tumours (Genomics England PanelApp)",
      31: "Neuroendocrine cancer pertinent cancer susceptibility (Genomics England PanelApp)",
      32: "Neurofibromatosis Type 1 (Genomics England PanelApp)",
      33: "Ovarian cancer pertinent cancer susceptibility (Genomics England PanelApp)",
      34: "Parathyroid Cancer (Genomics England PanelApp)",
      35: "Prostate cancer pertinent cancer susceptibility (Genomics England PanelApp)",
      36: "Renal cancer pertinent cancer susceptibility (Genomics England PanelApp)",
      37: "Rhabdoid tumour predisposition (Genomics England PanelApp)",
      38: "Sarcoma cancer susceptibility (Genomics England PanelApp)",
      39: "Sarcoma susceptbility (Genomics England PanelApp)",
      40: "Thyroid cancer pertinent cancer susceptibility (Genomics England PanelApp)",
      41: "Tumour predisposition - childhood onset (Genomics England PanelApp)",
      42: "Upper gastrointestinal cancer pertinent cancer susceptibility (Genomics England PanelApp)"
	}


def __main__():

   panels = "0 = CPSR exploratory cancer predisposition panel\n (n = 335, Genomics England PanelApp / TCGA Germline Study / Cancer Gene Census / Other)\n"
   panels = panels + "1 = Adult solid tumours cancer susceptibility (Genomics England PanelApp)\n"
   panels = panels + "2 = Adult solid tumours for rare disease (Genomics England PanelApp)\n"
   panels = panels + "3 = Bladder cancer pertinent cancer susceptibility (Genomics England PanelApp)\n"
   panels = panels + "4 = Brain cancer pertinent cancer susceptibility (Genomics England PanelApp)\n"
   panels = panels + "5 = Breast cancer pertinent cancer susceptibility (Genomics England PanelApp)\n"
   panels = panels + "6 = Childhood solid tumours cancer susceptibility (Genomics England PanelApp)\n"
   panels = panels + "7 = Colorectal cancer pertinent cancer susceptibility (Genomics England PanelApp)\n"
   panels = panels + "8 = Endometrial cancer pertinent cancer susceptibility (Genomics England PanelApp)\n"
   panels = panels + "9 = Familial Tumours Syndromes of the central & peripheral Nervous system (Genomics England PanelApp)\n"
   panels = panels + "10 = Familial breast cancer (Genomics England PanelApp)\n"
   panels = panels + "11 = Familial melanoma (Genomics England PanelApp)\n"
   panels = panels + "12 = Familial prostate cancer (Genomics England PanelApp)\n"
   panels = panels + "13 = Familial rhabdomyosarcoma (Genomics England PanelApp)\n"
   panels = panels + "14 = GI tract tumours (Genomics England PanelApp)\n"
   panels = panels + "15 = Genodermatoses with malignancies (Genomics England PanelApp)\n"
   panels = panels + "16 = Haematological malignancies cancer susceptibility (Genomics England PanelApp)\n"
   panels = panels + "17 = Haematological malignancies for rare disease (Genomics England PanelApp)\n"
   panels = panels + "18 = Head and neck cancer pertinent cancer susceptibility (Genomics England PanelApp)\n"
   panels = panels + "19 = Inherited MMR deficiency (Lynch syndrome) - Genomics England PanelApp\n"
   panels = panels + "20 = Inherited non-medullary thyroid cancer (Genomics England PanelApp)\n"
   panels = panels + "21 = Inherited ovarian cancer (without breast cancer) (Genomics England PanelApp)\n"
   panels = panels + "22 = Inherited pancreatic cancer (Genomics England PanelApp)\n"
   panels = panels + "23 = Inherited polyposis (Genomics England PanelApp)\n"
   panels = panels + "24 = Inherited predisposition to acute myeloid leukaemia (AML) - Genomics England PanelApp\n"
   panels = panels + "25 = Inherited predisposition to GIST (Genomics England PanelApp)\n"
   panels = panels + "26 = Inherited renal cancer (Genomics England PanelApp)\n"
   panels = panels + "27 = Inherited phaeochromocytoma and paraganglioma (Genomics England PanelApp)\n"
   panels = panels + "28 = Melanoma pertinent cancer susceptibility (Genomics England PanelApp)\n"
   panels = panels + "29 = Multiple endocrine tumours (Genomics England PanelApp)\n"
   panels = panels + "30 = Multiple monogenic benign skin tumours (Genomics England PanelApp)\n"
   panels = panels + "31 = Neuroendocrine cancer pertinent cancer susceptibility (Genomics England PanelApp)\n"
   panels = panels + "32 = Neurofibromatosis Type 1 (Genomics England PanelApp)\n"
   panels = panels + "33 = Ovarian cancer pertinent cancer susceptibility (Genomics England PanelApp)\n"
   panels = panels + "34 = Parathyroid Cancer (Genomics England PanelApp)\n"
   panels = panels + "35 = Prostate cancer pertinent cancer susceptibility (Genomics England PanelApp)\n"
   panels = panels + "36 = Renal cancer pertinent cancer susceptibility (Genomics England PanelApp)\n"
   panels = panels + "37 = Rhabdoid tumour predisposition (Genomics England PanelApp)\n"
   panels = panels + "38 = Sarcoma cancer susceptibility (Genomics England PanelApp)\n"
   panels = panels + "39 = Sarcoma susceptibility (Genomics England PanelApp)\n"
   panels = panels + "40 = Thyroid cancer pertinent cancer susceptibility (Genomics England PanelApp)\n"
   panels = panels + "41 = Tumour predisposition - childhood onset (Genomics England PanelApp)\n"
   panels = panels + "42 = Upper gastrointestinal cancer pertinent cancer susceptibility (Genomics England PanelApp)\n\n"
   
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

   optional_panel.add_argument('--panel_id',dest = "virtual_panel_id",type = str, default = "-1", help="Comma-separated string with identifier(s) of predefined virtual cancer predisposition gene panels,\n choose any combination of the following identifiers:\n" + str(panels))
   optional_panel.add_argument('--custom_list',dest = "custom_list",help="Provide custom list of genes from virtual panel 0 (single-column txt file with Ensembl gene identifiers),\n alternative to predefined panels provided with --panel_id)")
   optional_panel.add_argument('--custom_list_name',dest = "custom_list_name", default="Custom_Panel", help="Set name for custom made panel/list (single word - no whitespace), will be displayed in the report")
   optional_panel.add_argument('--diagnostic_grade_only', action="store_true",help="For panel_id's 1-42 (Genomics England PanelApp) - consider genes with a GREEN status only, default: %(default)s")


   optional_other.add_argument('--force_overwrite', action = "store_true", help='By default, the script will fail with an error if any output file already exists.\n You can force the overwrite of existing result files by using this flag, default: %(default)s')
   optional_other.add_argument('--version', action='version', version='%(prog)s ' + str(CPSR_VERSION))
   optional_other.add_argument('--basic',action="store_true",help="Run functional variant annotation on VCF through VEP/vcfanno, omit Tier assignment/report generation (STEP 4), default: %(default)s")
   optional_other.add_argument('--no_vcf_validate', action = "store_true",help="Skip validation of input VCF with Ensembl's vcf-validator, default: %(default)s")
   optional_other.add_argument('--docker_uid', dest='docker_user_id', help='Docker user ID. Default is the host system user ID. If you are experiencing permission errors,\n try setting this up to root (`--docker_uid root`), default: %(default)s')
   optional_other.add_argument('--no_docker', action='store_true', dest='no_docker', default=False, help='Run the CPSR workflow in a non-Docker mode, default: %(default)s')
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
   optional_other.add_argument('--classify_all', action='store_true',dest='classify_all',help='Provide CPSR variant classifications (TIER 1-5) also for variants with exising ClinVar classifications in output TSV, default: %(default)s')
   optional_other.add_argument('--clinvar_ignore_noncancer', action='store_true', help='Ignore (exclude from report) ClinVar-classified variants reported only for phenotypes/conditions NOT related to cancer, default: %(default)s')
   optional_other.add_argument('--debug',action='store_true',default=False, help='Print full docker commands to log, default: %(default)s')
   
   optional_vcfanno.add_argument('--vcfanno_n_proc', default = 4, type = int, help="Number of vcfanno processes (option '-p' in vcfanno), default: %(default)s")
   
   optional_vep.add_argument('--vep_n_forks', default = 4, type = int, help="Number of forks (option '--fork' in VEP), default: %(default)s")   
   optional_vep.add_argument('--vep_buffer_size', default = 500, type = int, help="Variant buffer size (variants read into memory simultaneously, option '--buffer_size' in VEP) " + \
      "\n- set lower to reduce memory usage, default: %(default)s")
   #optional_vep.add_argument('--vep_regulatory', action='store_true', help = 'Enable Variant Effect Predictor (VEP) to look for overlap with regulatory regions (option --regulatory in VEP).')
   optional_vep.add_argument('--vep_pick_order', default = "canonical,appris,biotype,ccds,rank,tsl,length,mane", help="Comma-separated string " + \
      "of ordered transcript properties for primary variant pick\n ( option '--pick_order' in VEP), default: %(default)s")
   optional_vep.add_argument('--vep_no_intergenic', action = "store_true", help="Skip intergenic variants during processing (option '--no_intergenic' in VEP), default: %(default)s")
   
   required.add_argument('--input_vcf', help='VCF input file with germline query variants (SNVs/InDels).', required = True)
   required.add_argument('--pcgr_dir',help="Directory that contains the PCGR data bundle directory, e.g. ~/pcgr-" + str(PCGR_VERSION), required = True)
   required.add_argument('--output_dir',help='Output directory', required = True)
   required.add_argument('--genome_assembly',choices = ['grch37','grch38'], help='Genome assembly build: grch37 or grch38', required = True)
   required.add_argument('--sample_id',help="Sample identifier - prefix for output files", required = True)
   
   args = parser.parse_args()
   arg_dict = vars(args)

   logger = getlogger('cpsr-validate-input-arguments')
   print()
   logger.info("STEP 0: Validate input data")

   arg_dict['vep_regulatory'] = True
   ## Required arguments
   ## Check that query VCF is set and exists
   if arg_dict['input_vcf'] is None or not os.path.exists(arg_dict['input_vcf']):
      err_msg = "Required argument '--input_vcf' does not exist (" + str(arg_dict['input_vcf']) + "). Type cpsr.py --help to view all options and required arguments"
      error_message(err_msg,logger)

   ## Check that PCGR directory (with data bundle) is provided and exists
   if arg_dict['pcgr_dir'] is None or not os.path.exists(arg_dict['pcgr_dir']):
      err_msg = "Required argument '--pcgr_dir' does not exist (" + str(arg_dict['pcgr_dir']) + "). Type cpsr.py --help to view all options and required arguments"
      error_message(err_msg,logger)
   
   ## Check that output directory is provided and exists
   if arg_dict['output_dir'] is None or not os.path.exists(arg_dict['output_dir']):
      err_msg = "Required argument '--output_dir' does not exist (" + str(arg_dict['output_dir']) + "). Type cpsr.py --help to view all options and required arguments"
      error_message(err_msg,logger)
   

   ## Check that genome assembly is set
   if arg_dict['genome_assembly'] is None:
      err_msg = "Required argument '--genome_assembly' has no/undefined value (\'" + str(arg_dict['genome_assembly']) + "'). Type cpsr.py --help to view all options and required arguments"
      error_message(err_msg,logger)
   
   ## Check that sample identifier is set and is of appropriate length (minimum two characters)
   if arg_dict['sample_id'] is None:
      err_msg = "Required argument '--sample_id' has no/undefined value (" + str(arg_dict['sample_id']) + "). Type cpsr.py --help to view all options and required arguments"
      error_message(err_msg,logger)

   if len(arg_dict['sample_id']) <= 2:
      err_msg = "Sample name identifier ('--sample_id') requires a name with more than two characters. Current sample identifier: " + str(arg_dict['sample_id'])
      error_message(err_msg,logger)

   ### Optional arguments
   ## Provide virtual_panel_id or a custom list from panel 0
   if arg_dict['virtual_panel_id'] == "-1" and not arg_dict['custom_list']:
      err_msg = 'Provide valid virtual panel identifier(s) through --panel_id (0 - 42) or provide custom list of panel 0 genes (single column text file) through --custom_list'
      error_message(err_msg,logger)

   if arg_dict['custom_list'] and arg_dict['virtual_panel_id'] != "-1":
      err_msg =  "Option --panel_id cannot be used in conjunction with --custom_list"
      error_message(err_msg, logger)

   if arg_dict['maf_upper_threshold'] <= 0 or arg_dict['maf_upper_threshold'] > 1:
      err_msg = 'MAF upper threshold must be greater than 0 and below 1, current value is ' + str(arg_dict['maf_upper_threshold'])
      error_message(err_msg,logger)

   if arg_dict['vcfanno_n_proc'] <= 0 or arg_dict['vcfanno_n_proc'] > 15:
      err_msg = 'Number of processes that vcfanno can use during annotation must be above 0 and not more than 15, current value is ' + str(arg_dict['vcfanno_n_proc'])
      error_message(err_msg,logger)
   
   
   ## Check that panel identifier(s) are set appropriately
   if arg_dict['virtual_panel_id'] != "-1" and not arg_dict['custom_list']:
      if not ',' in arg_dict['virtual_panel_id']:
         if str(arg_dict['virtual_panel_id']).isdigit():
            panel_id = int(arg_dict['virtual_panel_id'])
            if not (panel_id >= 0 and panel_id <= 42):
               err_msg =  'A single panel chosen with \'--panel_id\' must be in the range 0 - 42'
               error_message(err_msg, logger)
         else:
            err_msg =  'A single panel chosen with \'--panel_id\' must be a proper integer - not \'' + str(arg_dict['virtual_panel_id']) + '\''
            error_message(err_msg, logger)
      else:
         panels = str(arg_dict['virtual_panel_id']).split(',')
         for p in panels:
            #p = int(p)
            if str(p).isdigit():
               panel_id = int(p)
               if panel_id < 1 or panel_id > 42:
                  err_msg =  'Multiple panels submitted as comma-separated string with \'--panel_id\' must take values in the range 1 - 42'
                  error_message(err_msg, logger)
            else:
               err_msg =  'Multiple panels submitted as comma-separated string with \'--panel_id\' must contain proper integer values only -  \'' + str(arg_dict['virtual_panel_id']) + '\' contains non-integer entries'
               error_message(err_msg, logger)


   if (arg_dict['custom_list'] or arg_dict['virtual_panel_id'] == "0" ) and arg_dict['diagnostic_grade_only']:
      warn_msg = 'Option \'--diagnostic_grade_only\' applies ONLY to panel identifiers from Genomics England PanelApp - will be ignored'
      warn_message(warn_msg, logger)

   ## VEP options
   if arg_dict['vep_n_forks'] <= 0 or arg_dict['vep_n_forks'] > 4:
      err_msg = 'Number of forks that VEP can use during annotation must be above 0 and not more than 4, current value is ' + str(arg_dict['vep_n_forks'])
      error_message(err_msg,logger)
   
   if arg_dict['vep_buffer_size'] <= 0 or arg_dict['vep_buffer_size'] > 30000:
      err_msg = 'Internal VEP buffer size, corresponding to the number of variants that are read in to memory simultaneously, must be above 0 and not more than 30,000, current value is ' + str(arg_dict['vep_buffer_size'])
      error_message(err_msg,logger)

   ## Check that VEP pick criteria is formatted correctly
   if not arg_dict['vep_pick_order'] is None:
      values = str(arg_dict['vep_pick_order']).split(',')
      permitted_sources = ['canonical','appris','tsl','biotype','ccds','rank','length','mane']
      num_permitted_sources = 0
      for v in values:
         if v in permitted_sources:
            num_permitted_sources += 1
               
      if num_permitted_sources != 8:
         err_msg = "Option 'vep_pick_order' = " + str(arg_dict['vep_pick_order']) + " is formatted incorrectly, should be " + \
            "a comma-separated string of the following values: canonical,appris,tsl,biotype,ccds,rank,length,mane"
         error_message(err_msg, logger)

   ## Check that Docker image contains (if not --no_docker option set)
   global DOCKER_IMAGE_VERSION
   if arg_dict['no_docker']:
      DOCKER_IMAGE_VERSION = None
   else:
      # check that script and Docker image version correspond
      check_docker_command = 'docker images -q ' + str(DOCKER_IMAGE_VERSION)
      output = subprocess.check_output(str(check_docker_command), stderr=subprocess.STDOUT, shell=True)
      if(len(output) == 0):
          err_msg = 'Docker image ' + str(DOCKER_IMAGE_VERSION) + ' does not exist, pull image from Dockerhub (docker pull ' + str(DOCKER_IMAGE_VERSION) + ')'
          error_message(err_msg,logger)
   
   ## Map local input directories and files to internal paths/volumes in container (Docker)
   host_directories = verify_input_files(arg_dict, logger)
   
   ## Run CPSR workflow
   run_cpsr(arg_dict, host_directories)


def error_message(message, logger):
   logger.error('')
   logger.error(message)
   logger.error('')
   exit(0)
   #sys.exit(1)

def warn_message(message, logger):
   logger.warning(message)

def verify_input_files(arg_dict, logger):

   input_vcf_dir = "NA"
   db_dir = "NA"
   base_dir = "NA"
   output_dir_full = "NA"
   input_vcf_basename = "NA"
   input_customlist_basename = "NA"
   input_customlist_dir = "NA"
   
   ## check the existence of given output folder
   output_dir_full = os.path.abspath(arg_dict['output_dir'])
   if not os.path.isdir(output_dir_full):
      err_msg = "Output directory (" + str(output_dir_full) + ") does not exist"
      error_message(err_msg,logger)
   
   ## check if input BED exist
   if not arg_dict['custom_list'] is None:
      if not os.path.exists(os.path.abspath(arg_dict['custom_list'])):
         err_msg = "Input file (" + str(arg_dict['custom_list']) + ") does not exist"
         error_message(err_msg,logger)
      
      input_customlist_basename = os.path.basename(str(arg_dict['custom_list']))
      input_customlist_dir = os.path.dirname(os.path.abspath(arg_dict['custom_list']))

   ## check if input vcf exist
   if not arg_dict['input_vcf'] is None:
      if not os.path.exists(os.path.abspath(arg_dict['input_vcf'])):
         err_msg = "Input file (" + str(arg_dict['input_vcf']) + ") does not exist"
         error_message(err_msg,logger)

      if not (os.path.abspath(arg_dict['input_vcf']).endswith('.vcf') or os.path.abspath(arg_dict['input_vcf']).endswith('.vcf.gz')):
         err_msg = "VCF input file (" + os.path.abspath(arg_dict['input_vcf']) + ") does not have the correct file extension (.vcf or .vcf.gz)"
         error_message(err_msg,logger)

      ## check that tabix file exist if bgzipped files is given
      if os.path.abspath(arg_dict['input_vcf']).endswith('.vcf.gz'):
         tabix_file = arg_dict['input_vcf'] + '.tbi'
         if not os.path.exists(os.path.abspath(tabix_file)):
            err_msg = "Tabix file (i.e. '.gz.tbi') is not present for the bgzipped VCF input file (" + os.path.abspath(arg_dict['input_vcf']) + "). Please make sure your input VCF is properly compressed and indexed (bgzip + tabix)"
            error_message(err_msg,logger)

      input_vcf_basename = os.path.basename(str(arg_dict['input_vcf']))
      input_vcf_dir = os.path.dirname(os.path.abspath(arg_dict['input_vcf']))

      ## if output vcf exist and overwrite not set
      output_vcf = os.path.join(str(output_dir_full),str(arg_dict['sample_id'])) + '.cpsr.' + str(arg_dict['genome_assembly']) + '.vcf.gz'
      if os.path.exists(output_vcf) and arg_dict['force_overwrite'] is False:
         err_msg = "Output files (e.g. " + str(output_vcf) + ") already exist - please specify different sample_id or add option --force_overwrite"
         error_message(err_msg,logger)
   
   ## check the existence of base folder
   base_dir = os.path.abspath(arg_dict['pcgr_dir'])
   if not os.path.isdir(base_dir):
      err_msg = "Base directory (" + str(base_dir) + ") does not exist"
      error_message(err_msg,logger)
   
   ## check the existence of data folder within the base folder
   db_dir = os.path.join(os.path.abspath(arg_dict['pcgr_dir']),'data')
   if not os.path.isdir(db_dir):
      err_msg = "Data directory (" + str(db_dir) + ") does not exist"
      error_message(err_msg,logger)
   
   ## check the existence of specified assembly data folder within the base folder
   db_assembly_dir = os.path.join(os.path.abspath(arg_dict['pcgr_dir']),'data',arg_dict['genome_assembly'])
   if not os.path.isdir(db_assembly_dir):
      err_msg = "Data directory for the specified genome assembly (" + str(db_assembly_dir) + ") does not exist"
      error_message(err_msg,logger)
   
   ## check the existence of RELEASE_NOTES
   rel_notes_file = os.path.join(os.path.abspath(arg_dict['pcgr_dir']),'data',arg_dict['genome_assembly'],'RELEASE_NOTES')
   if not os.path.exists(rel_notes_file):
      err_msg = 'The PCGR data bundle is outdated - please download the latest data bundle (see github.com/sigven/cpsr for instructions)'
      error_message(err_msg,logger)
      
   f_rel_not = open(rel_notes_file,'r')
   compliant_data_bundle = 0
   for line in f_rel_not:
      if DB_VERSION in line:
         compliant_data_bundle = 1
   
   f_rel_not.close()
    
   if compliant_data_bundle == 0:
      err_msg = 'The PCGR data bundle is not compliant with the software version - please download the latest software and data bundle (see https://github.com/sigven/cpsr for instructions)'
      error_message(err_msg,logger)
   
   host_directories = {}
   host_directories['input_vcf_dir_host'] = input_vcf_dir
   host_directories['input_customlist_dir_host'] = input_customlist_dir
   host_directories['db_dir_host'] = db_assembly_dir
   host_directories['base_dir_host'] = base_dir
   host_directories['output_dir_host'] = output_dir_full
   host_directories['input_vcf_basename_host'] = input_vcf_basename
   host_directories['input_customlist_basename_host'] = input_customlist_basename

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

def run_cpsr(arg_dict, host_directories):
   """
   Main function to run the CPSR workflow (Docker/Conda)
   """

   ## get options
   global debug
   debug = arg_dict['debug']
   docker_user_id = arg_dict['docker_user_id']
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
   #vep_regulatory = "OFF"

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

   logger = getlogger('cpsr-validate-input-arguments')

   ## set basic Docker run commands
   output_vcf = 'None'
   output_pass_vcf = 'None'
   output_pass_tsv = 'None'
   uid = ''
   global GENCODE_VERSION, VEP_ASSEMBLY
   if arg_dict['genome_assembly'] == 'grch37':
      GENCODE_VERSION = '19'
      VEP_ASSEMBLY = 'GRCh37'

   if docker_user_id:
      uid = docker_user_id
   elif platform.system() == 'Linux' or platform.system() == 'Darwin' or sys.platform == 'darwin' or sys.platform == 'linux2' or sys.platform == 'linux':
      uid = os.getuid()
   else:
      if platform.system() == 'Windows' or sys.platform == 'win32' or sys.platform == 'cygwin':
         uid = getpass.getuser()

   if uid == '':
      logger.warning('Was not able to get user id/username for logged-in user on the underlying platform (platform.system(): ' + str(platform.system()) + ', sys.platform: ' + str(sys.platform) + '), now running CPSR as root')
      uid = 'root'
      
   vepdb_dir_host = os.path.join(str(host_directories['db_dir_host']),'.vep')
   

   input_vcf_docker = 'None'
   input_customlist_docker = 'None'
   

   ## Determine basic Docker commands

   if DOCKER_IMAGE_VERSION:

      vep_volume_mapping = str(vepdb_dir_host) + ":/usr/local/share/vep/data"
      databundle_volume_mapping = str(host_directories['base_dir_host']) + ":/data"
      input_vcf_volume_mapping = str(host_directories['input_vcf_dir_host']) + ":/workdir/input_vcf"
      input_customlist_volume_mapping = str(host_directories['input_customlist_dir_host']) + ":/workdir/input_custom"
      output_volume_mapping = str(host_directories['output_dir_host']) + ":/workdir/output"

      if host_directories['input_vcf_basename_host'] != 'NA':
         input_vcf_docker = '/workdir/input_vcf/' + str(host_directories['input_vcf_basename_host'])
      if host_directories['input_customlist_basename_host'] != 'NA':
         input_customlist_docker = '/workdir/input_custom/' + str(host_directories['input_customlist_basename_host'])

      docker_command_run1 = "docker run --rm -u " + str(uid) + " -v=" +  str(databundle_volume_mapping) + " -v=" + str(vep_volume_mapping) + " -v=" + str(output_volume_mapping)
      if host_directories['input_vcf_dir_host'] != 'NA':
         docker_command_run1 = docker_command_run1  + " -v=" + str(input_vcf_volume_mapping)
      if host_directories['input_customlist_dir_host'] != 'NA':
         docker_command_run1 = docker_command_run1  + " -v=" + str(input_customlist_volume_mapping)

      docker_command_run1 = docker_command_run1 + " -w=/workdir/output " + str(DOCKER_IMAGE_VERSION) + " sh -c \""
      docker_command_run2 = "docker run --rm -u " + str(uid) + " -v=" + str(databundle_volume_mapping) + " -v=" + str(output_volume_mapping) + " -w=/workdir/output " + str(DOCKER_IMAGE_VERSION) + " sh -c \""
      docker_command_run_end = '\"'

      data_dir = '/data'
      output_dir = '/workdir/output'
      vep_dir = '/usr/local/share/vep/data'
      r_scripts_dir = '/'

   ## If run in no-docker mode, set commands accordingly
   else:
      if host_directories['input_vcf_basename_host'] != 'NA':
         input_vcf_docker = os.path.join(host_directories['input_vcf_dir_host'], host_directories['input_vcf_basename_host'])
      if host_directories['input_customlist_basename_host'] != 'NA':
         input_customlist_docker = os.path.join(host_directories['input_customlist_dir_host'], host_directories['input_customlist_basename_host'])

      docker_command_run1 = ''
      docker_command_run2 = ''
      docker_command_run_end = ''

      data_dir = host_directories['base_dir_host']
      output_dir = host_directories['output_dir_host']
      vep_dir = vepdb_dir_host
      r_scripts_dir = ''

   check_subprocess(logger, docker_command_run1.replace("-u " + str(uid), "") + 'mkdir -p ' + output_dir + docker_command_run_end)

   ## CPSR|Validate input VCF - check formatting, non-overlap with CPSR INFO tags, and whether sample contains any variants in cancer predisposition loci
   vcf_validate_command = docker_command_run1 + "cpsr_validate_input.py" + " " + data_dir + " " + str(input_vcf_docker) + " " + \
      str(input_customlist_docker) + " " + str(preserved_info_tags) + " " + str(vcf_validation) + " " + str(arg_dict['genome_assembly']) + " " + \
      str(arg_dict['sample_id']) +  " " + str(virtual_panel_id) + " " + str(diagnostic_grade_only)
   if debug:
      vcf_validate_command  += ' --debug'
   if not DOCKER_IMAGE_VERSION:
      vcf_validate_command += ' --output_dir ' + output_dir + docker_command_run_end
   else:
      vcf_validate_command += docker_command_run_end
   check_subprocess(logger, vcf_validate_command)
   logger.info('Finished')


   ## CPSR|Start - log key information about run
   logger = getlogger("cpsr-start")
   print()
   logger.info("--- Cancer Predisposition Sequencing Reporter workflow ----")
   logger.info("Sample name: " + str(arg_dict['sample_id']))
   if not input_customlist_docker == 'None':
      logger.info("Virtual gene panel: custom-made list from panel 0: " + str(input_customlist_docker))
   else:
      #logger.info("Virtual gene panel(s): " + str(GE_panels[virtual_panel_id]))
      logger.info("Diagnostic-grade genes in virtual panels (GE PanelApp): " + str(diagnostic_grade_set))
   logger.info("Include incidential findings (ACMG recommended list v3.0): " + str(secondary_findings_set))
   logger.info("Include low to moderate cancer risk variants from genome-wide association studies: " + str(gwas_findings_set))
   logger.info("Reference population, germline variant frequencies (gnomAD): " + str(arg_dict['pop_gnomad']).upper())
   logger.info("Genome assembly: " + str(arg_dict['genome_assembly']))
   
   if not input_vcf_docker == 'None':
      
      ## Define input, output and temporary file names
      pcgr_model = 'cpsr'
      output_vcf = os.path.join(output_dir, str(arg_dict['sample_id']) + '.cpsr.' + str(arg_dict['genome_assembly']) + '.vcf.gz')
      output_pass_vcf = os.path.join(output_dir, str(arg_dict['sample_id']) + '.cpsr.' + str(arg_dict['genome_assembly']) + '.pass.vcf.gz')
      output_pass_tsv = os.path.join(output_dir, str(arg_dict['sample_id']) + '.cpsr.' + str(arg_dict['genome_assembly']) + '.pass.tsv')
      input_vcf_cpsr_ready = os.path.join(output_dir, re.sub(r'(\.vcf$|\.vcf\.gz$)','.cpsr_ready_target.vcf.gz',host_directories['input_vcf_basename_host']))
      input_vcf_cpsr_ready_uncompressed = os.path.join(output_dir, re.sub(r'(\.vcf$|\.vcf\.gz$)','.cpsr_ready_target.vcf',host_directories['input_vcf_basename_host']))
      vep_vcf = re.sub(r'(\.vcf$|\.vcf\.gz$)','.cpsr_vep.vcf',input_vcf_cpsr_ready)
      vep_vcfanno_vcf = re.sub(r'(\.vcf$|\.vcf\.gz$)','.cpsr_vep.vcfanno.vcf',input_vcf_cpsr_ready)
      vep_vcfanno_annotated_vcf = re.sub(r'\.vcfanno','.vcfanno.annotated',vep_vcfanno_vcf) + '.gz'
      vep_vcfanno_annotated_pass_vcf = re.sub(r'\.vcfanno','.vcfanno.annotated.pass',vep_vcfanno_vcf) + '.gz'
      custom_bed = os.path.join(output_dir, str(arg_dict['sample_id']) + '.' + str(pcgr_model) + '.' + str(arg_dict['genome_assembly']) + '.custom_list.bed')

      ## File names for assembly-specific genome fasta files (VEP)
      fasta_assembly = os.path.join(vep_dir, "homo_sapiens", str(VEP_VERSION) + "_" + str(VEP_ASSEMBLY), "Homo_sapiens." + str(VEP_ASSEMBLY) + ".dna.primary_assembly.fa.gz")
      ancestor_assembly = os.path.join(vep_dir, "homo_sapiens", str(VEP_VERSION) + "_" + str(VEP_ASSEMBLY), "human_ancestor.fa.gz")
      
      ## Set all flags used in VEP run
      plugins_in_use = "NearestExonJB, LoF"
      vep_flags = "--format vcf --vcf --check_ref --flag_pick_allele_gene --hgvs --dont_skip --failed 1 --af --af_1kg --af_gnomad " + \
         "--variant_class --domains --symbol --protein --ccds --uniprot --appris --biotype --canonical --gencode_basic --cache " + \
         "--numbers --total_length --no_stats --allele_number --no_escape --xref_refseq --plugin NearestExonJB,max_range=50000"
      vep_options = "--pick_order " + str(arg_dict['vep_pick_order']) + " --force_overwrite --buffer_size " + \
         str(arg_dict['vep_buffer_size']) + " --species homo_sapiens --assembly " + \
         str(VEP_ASSEMBLY) + " --offline --fork " + str(arg_dict['vep_n_forks']) + " " + str(vep_flags) + " --dir " + str(vep_dir)
      vep_options += " --cache_version " + str(VEP_VERSION)
      loftee_dir = '/opt/vep/src/ensembl-vep/modules'
      if arg_dict['vep_no_intergenic'] == 1:
         vep_options = vep_options + " --no_intergenic"
      if arg_dict['vep_regulatory'] == 1:
         vep_options = vep_options + " --regulatory"
      if arg_dict['genome_assembly'] == "grch38":
         vep_options = vep_options +  " --mane"
      if not DOCKER_IMAGE_VERSION:
         conda_prefix = os.path.dirname(os.path.dirname(sys.executable))
         loftee_dir = os.path.join(conda_prefix, 'share', 'loftee')
         assert os.path.isdir(loftee_dir), 'LoF VEP plugin is not found in ' + loftee_dir + '. Please make sure you installed pcgr conda package and have corresponding conda environment active.'
         vep_options += " --plugin LoF,loftee_path:" + loftee_dir + ",human_ancestor_fa:" + str(ancestor_assembly) + ",use_gerp_end_trunc:0 --dir_plugins " + loftee_dir
      else:
         vep_options += " --plugin LoF,loftee_path:" + loftee_dir + ",human_ancestor_fa:" + str(ancestor_assembly)  + ",use_gerp_end_trunc:0 --dir_plugins " + loftee_dir
      if not debug:
         vep_options += " --quiet"
      
      ## Compose full VEP command
      vep_main_command = str(docker_command_run1) + "vep --input_file " + str(input_vcf_cpsr_ready) + " --output_file " + str(vep_vcf) + \
         " " + str(vep_options) + " --fasta " + str(fasta_assembly) + docker_command_run_end
      vep_bgzip_command = str(docker_command_run1) + "bgzip -f " + str(vep_vcf) + docker_command_run_end
      vep_tabix_command = str(docker_command_run1) + "tabix -f -p vcf " + str(vep_vcf) + ".gz" + docker_command_run_end
      logger = getlogger('cpsr-vep')

      ## CPSR|VEP - run Variant Effect Predictor on query VCF with LoF and NearestExonJB plugins
      print()
      logger.info("STEP 1: Basic variant annotation with Variant Effect Predictor (" + str(VEP_VERSION) + ", GENCODE release " + \
         str(GENCODE_VERSION) + ", " + str(arg_dict['genome_assembly']) + ")")
      logger.info("VEP configuration - one primary consequence block pr. alternative allele (--flag_pick_allele)")
      logger.info("VEP configuration - transcript pick order: " + str(arg_dict['vep_pick_order']))
      logger.info("VEP configuration - transcript pick order: See more at https://www.ensembl.org/info/docs/tools/vep/script/vep_other.html#pick_options")
      logger.info("VEP configuration - skip intergenic: " + str(arg_dict['vep_no_intergenic']))
      logger.info("VEP configuration - look for overlap with regulatory regions: " + str(vep_regulatory))
      logger.info("VEP configuration - plugins in use: " + str(plugins_in_use))
      logger.info("VEP configuration - buffer_size/number of forks: " + str(arg_dict['vep_buffer_size']) + '/' + str(arg_dict['vep_n_forks']))
      check_subprocess(logger, vep_main_command)
      check_subprocess(logger, vep_bgzip_command)
      check_subprocess(logger, vep_tabix_command)
      logger.info("Finished")
   
      ## CPSR|vcfanno - run vcfanno on query VCF with a number of relevant annotated VCFs
      print()
      logger = getlogger('cpsr-vcfanno')
      logger.info("STEP 2: Annotation for cancer predisposition with cpsr-vcfanno (ClinVar, CIViC, dbNSFP, dbMTS, UniProtKB, cancerhotspots.org, ncER, GERP RS scores, GWAS catalog, gnomAD non-cancer subset)")
      pcgr_vcfanno_command = str(docker_command_run2) + "pcgr_vcfanno.py --num_processes "  + str(arg_dict['vcfanno_n_proc']) + \
         " --dbnsfp --clinvar --cancer_hotspots --dbmts --ncer --gerp --civic --uniprot --gnomad_cpsr --pcgr_onco_xref --gwas --rmsk " + str(vep_vcf) + ".gz " + \
         str(vep_vcfanno_vcf) + " " + os.path.join(data_dir, "data", str(arg_dict['genome_assembly'])) + docker_command_run_end      
      check_subprocess(logger, pcgr_vcfanno_command)
      logger.info("Finished")

      ## CPSR|summarise - expand annotations with separate VCF INFO tags
      print()
      logger = getlogger("cpsr-summarise")
      pcgr_summarise_command = str(docker_command_run2) + "pcgr_summarise.py " + str(vep_vcfanno_vcf) + ".gz 0 " + \
         str(vep_regulatory) + " " + os.path.join(data_dir, "data", str(arg_dict['genome_assembly'])) + " --cpsr" + docker_command_run_end
      if debug:
         pcgr_summarise_command  += ' --debug'
      logger.info("STEP 3: Cancer gene annotations with cpsr-summarise")
      check_subprocess(logger, pcgr_summarise_command)

      ## CPSR|clean - rename output files, remove temporary files
      create_output_vcf_command1 = str(docker_command_run2) + 'mv ' + str(vep_vcfanno_annotated_vcf) + ' ' + str(output_vcf) + docker_command_run_end
      create_output_vcf_command2 = str(docker_command_run2) + 'mv ' + str(vep_vcfanno_annotated_vcf) + '.tbi ' + str(output_vcf) + '.tbi' + docker_command_run_end
      create_output_vcf_command3 = str(docker_command_run2) + 'mv ' + str(vep_vcfanno_annotated_pass_vcf) + ' ' + str(output_pass_vcf) + docker_command_run_end
      create_output_vcf_command4 = str(docker_command_run2) + 'mv ' + str(vep_vcfanno_annotated_pass_vcf) + '.tbi ' + str(output_pass_vcf) + '.tbi' + docker_command_run_end
      clean_command = str(docker_command_run2) + 'rm -f ' + str(vep_vcf) + '* ' + str(vep_vcfanno_annotated_vcf) + ' ' + \
         str(vep_vcfanno_annotated_pass_vcf) + '* ' + str(vep_vcfanno_vcf) + '* ' +  str(input_vcf_cpsr_ready_uncompressed) + "* " + docker_command_run_end
      check_subprocess(logger, create_output_vcf_command1)
      check_subprocess(logger, create_output_vcf_command2)
      check_subprocess(logger, create_output_vcf_command3)
      check_subprocess(logger, create_output_vcf_command4)

      ## CPSR|vcf2tsv - perform vcf2tsv conversion on the final annotated VCF file
      cpsr_vcf2tsv_command = str(docker_command_run2) + "vcf2tsv.py " + str(output_pass_vcf) + " --compress " + str(output_pass_tsv) + docker_command_run_end
      logger.info("Converting VCF to TSV with https://github.com/sigven/vcf2tsv")
      check_subprocess(logger, cpsr_vcf2tsv_command)
      if not debug:
         check_subprocess(logger, clean_command)
      logger.info("Finished")

  
   print()
   
   ## Generation of HTML reports for VEP/vcfanno-annotated VCF file
   if not arg_dict['basic']: 
      logger = getlogger('cpsr-writer')
      logger.info("STEP 4: Generation of output files - Cancer predisposition sequencing report")
      cpsr_report_command = docker_command_run1 + os.path.join(r_scripts_dir, "cpsr.R") + " " + output_dir + " " + \
         str(output_pass_tsv) + ".gz " +  \
         str(arg_dict['sample_id'])  + " " + \
         str(PCGR_VERSION) + " " + \
         str(CPSR_VERSION) + " " + \
         str(arg_dict['genome_assembly']) + " " + \
         str(data_dir) + " " + \
         str(virtual_panel_id) + " " + \
         str(preserved_info_tags) + " " + \
         str(custom_bed) + " " + \
         str(arg_dict['custom_list_name']) + " " +  \
         str(arg_dict['report_theme']) + " " + \
         str(arg_dict['report_table_display']) + " " + \
         str(report_nonfloating_toc) + " " + \
         str(gwas_findings) + " " + \
         str(arg_dict['gwas_p_value']) + " " + \
         str(arg_dict['pop_gnomad']) + " " + \
         str(arg_dict['maf_upper_threshold']) + " " + \
         str(arg_dict['vep_pick_order']) + " " + \
         str(arg_dict['vep_n_forks']) + " " + \
         str(arg_dict['vep_buffer_size']) + " " + \
         str(vep_no_intergenic) + " " + \
         str(vep_regulatory) + " " + \
         str(secondary_findings) + " " + \
         str(classify_all) + " " + \
         str(ignore_noncoding) + " " + \
         str(clinvar_ignore_noncancer) + " " + \
         str(diagnostic_grade_only) + docker_command_run_end

      check_subprocess(logger, cpsr_report_command)
      logger.info("Finished")
   

   
if __name__=="__main__": __main__()

