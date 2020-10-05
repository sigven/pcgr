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
import toml
from argparse import RawTextHelpFormatter


PCGR_VERSION = 'dev'
DB_VERSION = 'PCGR_DB_VERSION = 20200920'
VEP_VERSION = '101'
GENCODE_VERSION = '35'
NCBI_BUILD_MAF = "GRCh38"
VEP_ASSEMBLY = "GRCh38"
DOCKER_IMAGE_VERSION = 'sigven/pcgr:' + str(PCGR_VERSION)


#global vep_assembly
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
   program_options = "--input_vcf <INPUT_VCF> --pcgr_dir <PCGR_DIR> --output_dir <OUTPUT_DIR> --genome_assembly " + \
      " <GENOME_ASSEMBLY> --conf <CONFIG_FILE> --sample_id <SAMPLE_ID>"

   parser = argparse.ArgumentParser(description = program_description,
                                    formatter_class=RawTextHelpFormatter, usage="%(prog)s -h [options] " + str(program_options))
   parser._action_groups.pop()
   required = parser.add_argument_group('Required arguments')
   optional = parser.add_argument_group('Optional arguments')

   required.add_argument('--input_vcf', dest = "input_vcf", help = 'VCF input file with somatic variants in tumor sample, SNVs/InDels', required = True)
   optional.add_argument('--input_cna', dest = "input_cna", help = 'Somatic copy number alteration segments (tab-separated values)')
   optional.add_argument('--logr_gain', type = float, default = 0.8, dest = "logr_gain",help = 'Log ratio-threshold for regions containing copy number gains/amplifications (default: %(default)s)')
   optional.add_argument('--logr_homdel', type = float, default = -0.8, dest = "logr_homdel",help = 'Log ratio-threshold for regions containing homozygous deletions (default: %(default)s)')
   optional.add_argument('--cna_overlap_pct', type = float, default = 50, dest = "cna_overlap_pct", help = 'Mean percent overlap between copy number segment and gene transcripts for reporting of gains/losses in tumor suppressor genes/oncogenes, (default: %(default)s)')
   #optional.add_argument('--input_cna_plot', dest = "input_cna_plot", help = 'Somatic copy number alteration plot, default: %(default)s')
   optional.add_argument('--pon_vcf', dest = "pon_vcf", help = "VCF file with germline calls from Panel of Normals (PON) - i.e. blacklisted variants, (default: %(default)s)")
   optional.add_argument('--tumor_site', dest = "tsite", type = int, default = 0, help = "Optional integer code to specify primary tumor type/site of query sample,\n choose any of the following identifiers:\n" + str(tumor_sites) + "(default: %(default)s - any tumor type)")
   optional.add_argument('--tumor_purity', type = float, dest = "tumor_purity", help = "Estimated tumor purity (between 0 and 1, (default: %(default)s)")
   optional.add_argument('--tumor_ploidy', type = float, dest = "tumor_ploidy", help = "Estimated tumor ploidy (default: %(default)s)")
   optional.add_argument('--tumor_dp_min', type = int, default = 0, dest = "tumor_dp_min",help = "If VCF INFO tag for sequencing depth (tumor) is provided and found, set minimum required depth for inclusion in report (default: %(default)s)")
   optional.add_argument('--tumor_af_min', type = float, default = 0, dest = "tumor_af_min",help = "If VCF INFO tag for variant allelic fraction (tumor) is provided and found, set minimum required AF for inclusion in report (default: %(default)s)")
   optional.add_argument('--control_dp_min', type = int, default = 0, dest = "control_dp_min",help = "If VCF INFO tag for sequencing depth (control) is provided and found, set minimum required depth for inclusion in report (default: %(default)s)")
   optional.add_argument('--control_af_max', type = float, default = 1, dest = "control_af_max",help = "If VCF INFO tag for variant allelic fraction (control) is provided and found, set maximum tolerated AF for inclusion in report (default: %(default)s)")
   optional.add_argument('--target_size_mb', type = float, default = 34, dest = "target_size_mb",help = "For mutational burden analysis - approximate protein-coding target size of sequencing assay (default: %(default)s Mb (WES/WGS))")
   optional.add_argument('--tumor_only',action = "store_true",help = "Input VCF comes from tumor-only sequencing, calls will be filtered for variants of germline origin (set configurations for filtering in .toml file), (default: %(default)s)")
   optional.add_argument('--cell_line',action = "store_true",help = "Input VCF comes from tumor cell line sequencing (requires --tumor_only), calls will be filtered for variants of germline origin (set configurations for filtering in .toml file), (default: %(default)s)")
   optional.add_argument('--assay', dest = 'assay', choices = ['WES','WGS','TARGETED'], default = "WES", help = 'Type of DNA sequencing assay performed for input data (VCF) default: %(default)s')
   optional.add_argument('--include_trials',action = "store_true",help = "(Beta) Include relevant ongoing or future clinical trials, focusing on studies with molecularly targeted interventions")
   optional.add_argument('--estimate_tmb', action = "store_true", help = "Estimate tumor mutational burden from the total number of somatic mutations and target region size, default: %(default)s")
   optional.add_argument('--estimate_msi_status', action = "store_true", help = "Predict microsatellite instability status from patterns of somatic mutations/indels, default: %(default)s")
   optional.add_argument('--estimate_signatures', action = "store_true", help = "Estimate relative contributions of reference mutational signatures in query sample and detect potential kataegis events), default: %(default)s")
   optional.add_argument('--tmb_algorithm', dest = "tmb_algorithm", default = "all_coding", choices = ['all_coding','nonsyn'], help = "Method for calculation of TMB, all coding variants (Chalmers et al., Genome Medicine, 2017), or non-synonymous variants only, default: %(default)s")
   optional.add_argument('--min_mutations_signatures', type = int, default = 200, dest = "min_mutations_signatures", help = "Minimum number of SNVs required for reconstruction of mutational signatures (SBS) by MutationalPatterns (default: %(default)s, minimum n = 100)")
   optional.add_argument('--all_reference_signatures', action = "store_true", help = "Use all reference mutational signatures (SBS, n = 67) in signature reconstruction rather than only those already attributed to the tumor type (default: %(default)s)")
   optional.add_argument('--force_overwrite', action = "store_true", help = 'By default, the script will fail with an error if any output file already exists. You can force the overwrite of existing result files by using this flag')
   optional.add_argument('--version', action = 'version', version='%(prog)s ' + str(PCGR_VERSION))
   optional.add_argument('--basic', action = "store_true", help = "Run functional variant annotation on VCF through VEP/vcfanno, omit other analyses (i.e. Tier assignment/MSI/TMB/Signatures etc. and report generation (STEP 4), default: %(default)s")
   optional.add_argument('--no_vcf_validate', action = "store_true",help = "Skip validation of input VCF with Ensembl's vcf-validator")
   optional.add_argument('--docker-uid', dest = 'docker_user_id', help = 'Docker user ID. default is the host system user ID. If you are experiencing permission errors, try setting this up to root (`--docker-uid root`)')
   optional.add_argument('--no-docker', action ='store_true', dest = 'no_docker', default=False, help = 'Run the PCGR workflow in a non-Docker mode (see install_no_docker/ folder for instructions)')   
   optional.add_argument('--debug', action ='store_true', default=False, help = 'Print full Docker commands to log, default: %(default)s')
   required.add_argument('--pcgr_dir', dest = 'pcgr_dir', help = 'PCGR base directory with accompanying data directory, e.g. ~/pcgr-' + str(PCGR_VERSION), required = True)
   required.add_argument('--output_dir', dest = 'output_dir', help = 'Output directory', required = True)
   required.add_argument('--genome_assembly', dest = 'genome_assembly', choices = ['grch37','grch38'], help = 'Human genome assembly build: grch37 or grch38', required = True)
   required.add_argument('--conf', dest = 'configuration_file', help = 'PCGR configuration file in TOML format', required = True)
   required.add_argument('--sample_id', dest = 'sample_id', help = "Tumor sample/cancer genome identifier - prefix for output files", required = True)

   ## Parse required and optional arguments
   args = parser.parse_args()
   arg_dict = vars(args)

   logger = getlogger("pcgr-validate-arguments-input")
   print()
   logger.info("PCGR - STEP 0: Validate input data and options")

   ## Required arguments
   ## Check the existence of required arguments
   if arg_dict['pcgr_dir'] is None or not os.path.exists(arg_dict['pcgr_dir']):
      err_msg = "Required argument --pcgr_dir has no/undefined value (" + str(arg_dict['pcgr_dir']) + "). Type pcgr.py --help to view all options and required arguments"
      pcgr_error_message(err_msg,logger)
   
   if arg_dict['output_dir'] is None or not os.path.exists(arg_dict['output_dir']):
      err_msg = "Required argument --output_dir has no/undefined value (" + str(arg_dict['output_dir']) + "). Type pcgr.py --help to view all options and required arguments"
      pcgr_error_message(err_msg,logger)
   
   if arg_dict['configuration_file'] is None or not os.path.exists(arg_dict['configuration_file']):
      err_msg = "Required argument --conf has no/undefined value (" + str(arg_dict['configuration_file']) + "). Type pcgr.py --help to view all options and required arguments"
      pcgr_error_message(err_msg,logger)
   
   if arg_dict['genome_assembly'] is None:
      err_msg = "Required argument --genome_assembly has no/undefined value (" + str(arg_dict['genome_assembly']) + "). Type pcgr.py --help to view all options and required arguments"
      pcgr_error_message(err_msg,logger)
   
   if arg_dict['input_vcf'] is None:
      err_msg = "Required argument --input_vcf has no/undefined value (" + str(arg_dict['input_vcf']) + "). Type pcgr.py --help to view all options and required arguments"
      pcgr_error_message(err_msg,logger)
   
   if arg_dict['sample_id'] is None:
      err_msg = "Required argument --sample_id has no/undefined value (" + str(arg_dict['sample_id']) + "). Type pcgr.py --help to view all options and required arguments"
      pcgr_error_message(err_msg,logger)
   
   if len(arg_dict['sample_id']) <= 2 or len(arg_dict['sample_id']) > 35:
      err_msg = "Sample name identifier (--sample_id) requires a name with more than 2 characters (and less than 35). Current sample identifier: " + str(arg_dict['sample_id'])
      pcgr_error_message(err_msg,logger)

   ## Optional arguments
   ## Check the existence of Docker (if not --no_docker is et)
   global debug
   debug = args.debug
   global DOCKER_IMAGE_VERSION

   if arg_dict['no_docker']:
      DOCKER_IMAGE_VERSION = None
   else:
      # check that script and Docker image version correspond
      check_docker_command = 'docker images -q ' + str(DOCKER_IMAGE_VERSION)
      output = subprocess.check_output(str(check_docker_command), stderr=subprocess.STDOUT, shell=True)
      if(len(output) == 0):
         err_msg = 'Docker image ' + str(DOCKER_IMAGE_VERSION) + ' does not exist, pull image from Dockerhub (docker pull ' + str(DOCKER_IMAGE_VERSION) + ')'
         pcgr_error_message(err_msg,logger)

   
   ## read PCGR configuration file
   config_options = {}
   if os.path.exists(arg_dict['configuration_file']):
      config_options = read_config_options(arg_dict, logger)
   else:
      err_msg = "PCGR configuration file " + str(arg_dict['configuration_file']) + " does not exist - exiting"
      pcgr_error_message(err_msg,logger)

   ## check if input is cancer cell line, requires --tumor_only
   if arg_dict['cell_line'] and not arg_dict['tumor_only']:
      err_msg = "Analysis of cell line (--cell_line) needs option --tumor_only"
      pcgr_error_message(err_msg,logger)

   ## check that tumor primary site/type is set correctly (integer between 0 and 30)
   if arg_dict['tsite'] > max(tsites.keys()) or arg_dict['tsite'] < 0:
      err_msg = "Tumor type code (" + str(arg_dict['tsite']) + ") needs to be between " + str(min(tsites.keys())) + " and " + str(max(tsites.keys()))
      pcgr_error_message(err_msg,logger)

   ## check that tumor purity and tumor ploidy is set correctly
   if not arg_dict['tumor_purity'] is None:
      if not (arg_dict['tumor_purity'] > 0 and arg_dict['tumor_purity'] <= 1):
         err_msg = "Tumor purity value " + str(arg_dict['tumor_purity']) + " is not within a valid range [0,1]"
         pcgr_error_message(err_msg,logger)
   
   if not arg_dict['tumor_ploidy'] is None:
      if not arg_dict['tumor_ploidy'] > 0:
         err_msg = "Tumor ploidy value " + str(arg_dict['tumor_ploidy']) + " is negative"
         pcgr_error_message(err_msg,logger)

   ## check that minimum/maximum depth/allelic fractions are set correctly
   if arg_dict['tumor_dp_min'] < 0:
      err_msg = "Minimum depth tumor ('tumor_dp_min') must be positive, current value is " + str(arg_dict['tumor_dp_min']) + ")"
      pcgr_error_message(err_msg, logger)

   if arg_dict['tumor_af_min'] < 0 or arg_dict['tumor_af_min'] > 1:
      err_msg = "Minimum AF tumor ('tumor_af_min') must be within the [0,1] range, current value is " + str(arg_dict['tumor_af_min']) + ")"
      pcgr_error_message(err_msg, logger)
   
   if arg_dict['control_dp_min'] < 0:
      err_msg = "Minimum depth control ('control_dp_min') must be positive, current value is " + str(arg_dict['control_dp_min']) + ")"
      pcgr_error_message(err_msg, logger)

   if arg_dict['control_af_max'] < 0 or arg_dict['control_af_max'] > 1:
      err_msg = "Maximum AF control ('control_af_max') must be within the [0,1] range, current value is " + str(arg_dict['control_af_max']) + ")"
      pcgr_error_message(err_msg, logger)
   
   ## Check that coding target size region of sequencing assay is set correctly 
   if arg_dict['target_size_mb'] < 0 or arg_dict['target_size_mb'] > 34:
      err_msg = "Coding target size region in Mb (" + str(arg_dict['target_size_mb']) + ") is not positive or larger than the likely maximum size of the coding human genome (34 Mb))"
      pcgr_error_message(err_msg,logger)
   if arg_dict['target_size_mb'] < 1:
      warn_msg = "Coding target size region in Mb (" + str(arg_dict['target_size_mb']) + ") must be greater than 1 Mb for mutational burden estimate to be robust"
      pcgr_warn_message(warn_msg,logger)

   ## if assay is targeted or mode is Tumor-Only, MSI prediction will not be performed/switched off
   assay_type = "Tumor-Control"
   if arg_dict['estimate_msi_status'] is True and (arg_dict['assay'] == 'TARGETED' or arg_dict['tumor_only'] is True):
      if arg_dict['tumor_only'] is True:
         assay_type = "Tumor-Only"
      warn_msg = "MSI status prediction can be applied for WGS/WES tumor-control assays only (query type: " + str(arg_dict['assay']) + "|" + str(assay_type) + ") - analysis will be omitted"
      pcgr_warn_message(warn_msg, logger)
      arg_dict['estimate_msi_status'] = 0
   
   ## minimum number of mutations required for mutational signature reconstruction cannot be less than 100 (somewhat arbitrary lower threshold, recommended value is 200)
   if arg_dict['min_mutations_signatures'] < 200:
      warn_msg = "Minimum number of mutations required for mutational signature analysis ( n = " + str(arg_dict['min_mutations_signatures']) + ") is less than the recommended number (n = 200)"
      pcgr_warn_message(warn_msg,logger)
      if arg_dict['min_mutations_signatures'] < 100:
         err_msg = "Minimum number of mutations required for mutational signature analysis ( n = " + str(arg_dict['min_mutations_signatures']) + ") can not be less than 100"
         pcgr_err_message(err_msg,logger)

   ## if MSI status is to be estimated, mutational burden must be turned on
   if arg_dict['estimate_msi_status'] is True and arg_dict['estimate_tmb'] is False:
      err_msg = "Prediction of MSI status (option '--estimate_msi_status') requires mutational burden analysis ('--estimate_tmb')"
      pcgr_error_message(err_msg,logger)

   ## Emit warnings that mutational burden and mutational signatures are less accurate for assays with tumor-only data
   if arg_dict['tumor_only'] is True:
      if arg_dict['estimate_tmb'] is True:
         warn_msg = 'Estimation of mutational burden in tumor-only mode is suboptimal - results must be interpreted with caution'
         pcgr_warn_message(warn_msg,logger)
      if arg_dict['estimate_signatures'] is True:
         warn_msg = 'Estimation of mutational signatures in tumor-only mode is suboptimal - results must be interpreted with caution'
         pcgr_warn_message(warn_msg,logger)
   
   ## Emit warning that mutational signature estimation is (likely) not optimal for small/targeted sequencing assays
   if arg_dict['estimate_signatures'] is True and arg_dict['assay'] == 'TARGETED':
      warn_msg = 'Estimation of mutational signatures (--estimate_signatures) is not optimal for TARGETED sequencing assays - results must be interpreted with caution'
      pcgr_warn_message(warn_msg,logger)

   ## Check that log ratio thresholds for homozygous deletions and amplifications are properly set, and that segment overlap with transcripts are set appropriately
   if arg_dict['logr_homdel'] >= 0:
      err_msg = "Log ratio for homozygous deletions (" + str(arg_dict['logr_homdel']) + ") should be less than zero"
      pcgr_error_message(err_msg,logger)
   if arg_dict['logr_gain'] <= 0:
      err_msg = "Log ratio for copy number gains/amplifications (" + str(arg_dict['logr_gain']) + ") should be greater than zero"
      pcgr_error_message(err_msg,logger)
   if arg_dict['cna_overlap_pct'] > 100 or arg_dict['cna_overlap_pct'] <= 0:
      err_msg = "Minimum percent overlap between copy number segment and gene transcript (" + str(arg_dict['cna_overlap_pct']) + ") should be greater than zero and less than 100"
      pcgr_error_message(err_msg,logger)

   config_options['tumor_purity'] = "NA"
   config_options['tumor_ploidy'] = "NA"
   if not arg_dict['tumor_purity'] is None:
      config_options['tumor_purity'] = float(arg_dict['tumor_purity'])
   if not arg_dict['tumor_ploidy'] is None:
      config_options['tumor_ploidy'] = float(arg_dict['tumor_ploidy'])
   config_options['assay'] = arg_dict['assay']
   config_options['msi'] = {}
   config_options['msi']['run'] = int(arg_dict['estimate_msi_status'])
   config_options['msigs'] = {}
   config_options['msigs']['run'] = int(arg_dict['estimate_signatures'])
   config_options['msigs']['mutation_limit'] = int(arg_dict['min_mutations_signatures'])
   config_options['msigs']['all_reference_signatures'] = int(arg_dict['all_reference_signatures'])
   config_options['tmb'] = {}
   config_options['tmb']['run'] = int(arg_dict['estimate_tmb'])
   config_options['tmb']['target_size_mb'] = arg_dict['target_size_mb']
   config_options['tmb']['algorithm'] = arg_dict['tmb_algorithm']
   config_options['cna'] = {}
   config_options['cna']['logR_homdel'] = float(arg_dict['logr_homdel'])
   config_options['cna']['logR_gain'] = float(arg_dict['logr_gain'])
   config_options['cna']['cna_overlap_pct'] = float(arg_dict['cna_overlap_pct'])
   #config_options['allelic_support'] = {}
   if not config_options['allelic_support'] is None:
      config_options['allelic_support']['tumor_dp_min'] = int(arg_dict['tumor_dp_min'])
      config_options['allelic_support']['control_dp_min'] = int(arg_dict['control_dp_min'])
      config_options['allelic_support']['tumor_af_min'] = float(arg_dict['tumor_af_min'])
      config_options['allelic_support']['control_af_max'] = float(arg_dict['control_af_max'])

   ## Verify existence of input files
   host_directories = verify_input_files(arg_dict, logger)

   ## Run PCGR workflow ( VEP + vcfanno + summarise + vcf2tsv + HTML report generation)
   run_pcgr(arg_dict, host_directories, config_options)


def read_config_options(arg_dict, logger):
   
   ## read default options
   pcgr_config_options = {}
   pcgr_config_file_default = os.path.join(arg_dict['pcgr_dir'], 'data', str(arg_dict['genome_assembly']), 'pcgr_configuration_default.toml')
   if not os.path.exists(pcgr_config_file_default):
      err_msg = "Default PCGR configuration file " + str(pcgr_config_file_default) + " does not exist - exiting"
      pcgr_error_message(err_msg,logger)
   try:
      pcgr_config_options = toml.load(pcgr_config_file_default)
   except (IndexError,TypeError):
      err_msg = 'Configuration file ' + str(pcgr_config_file_default) + ' is not formatted correctly'
      pcgr_error_message(err_msg, logger)

   ## override with options set by the users
   try:
      user_options = toml.load(arg_dict['configuration_file'])
   except (IndexError,TypeError):
      err_msg = 'Configuration file ' + str(arg_dict['configuration_file']) + ' is not formatted correctly'
      pcgr_error_message(err_msg, logger)

   for section in pcgr_config_options:
      if section in user_options:
         for var in pcgr_config_options[section]:
            if not var in user_options[section]:
               continue
            if isinstance(pcgr_config_options[section][var],bool) and not isinstance(user_options[section][var],bool):
               err_msg = 'Configuration value ' + str(user_options[section][var]) + ' for ' + str(var) + ' cannot be parsed properly (expecting boolean)'
               pcgr_error_message(err_msg, logger)
            if isinstance(pcgr_config_options[section][var],int) and not isinstance(user_options[section][var],int):
                  err_msg = 'Configuration value \"' + str(user_options[section][var]) + '\" for ' + str(var) + ' cannot be parsed properly (expecting integer)'
                  pcgr_error_message(err_msg, logger)
            if isinstance(pcgr_config_options[section][var],float) and (not isinstance(user_options[section][var],float) and not isinstance(user_options[section][var],int)):
               err_msg = 'Configuration value ' + str(user_options[section][var]) + ' for ' + str(var) + ' cannot be parsed properly (expecting float)'
               pcgr_error_message(err_msg, logger)
            if isinstance(pcgr_config_options[section][var],str) and not isinstance(user_options[section][var],str):
               err_msg = 'Configuration value ' + str(user_options[section][var]) + ' for ' + str(var) + ' cannot be parsed properly (expecting string)'
               pcgr_error_message(err_msg, logger)
            normalization_options = ['default','exome','genome','exome2genome']
            theme_options = ['default', 'cerulean', 'journal', 'flatly', 'readable', 'spacelab', 'united', 'cosmo', 'lumen', 'paper', 'sandstone', 'simplex','yeti']
            if var == 'report_theme' and not str(user_options[section][var]) in theme_options:
               err_msg = 'Configuration value ' + str(user_options[section][var]) + ' for ' + str(var) + \
                  ' cannot be parsed properly (expecting \'default\', \'cerulean\', \'journal\', \'flatly\', \'readable\', \'spacelab\', \'united\', \'cosmo\', \'lumen\', \'paper\', \'sandstone\', \'simplex\',or \'yeti\')'
               pcgr_error_message(err_msg, logger)
            if var.startswith('maf_'):
               if user_options['tumor_only'][var] < 0 or user_options[section][var] > 1:
                  err_msg = "MAF value: " + str(var) + " must be within the [0,1] range, current value is " + str(user_options[section][var]) + ")"
                  pcgr_error_message(err_msg,logger)
            if var == 'vep_pick_order':
               values = str(user_options['other'][var]).split(',')
               permitted_sources = ['canonical','appris','tsl','biotype','ccds','rank','length','mane']
               num_permitted_sources = 0
               for v in values:
                  if v in permitted_sources:
                     num_permitted_sources += 1
               
               if num_permitted_sources != 8:
                  err_msg = "Configuration value vep_pick_order = " + str(user_options['other']['vep_pick_order']) + " is formatted incorrectly should be a comma-separated string of the following values: canonical,appris,tsl,biotype,ccds,rank,length, mane"
                  pcgr_error_message(err_msg, logger)

            pcgr_config_options[section][var] = user_options[section][var]
   
   pcgr_config_options['tumor_type'] = {}
   pcgr_config_options['tumor_type']['type'] = str(tsites[arg_dict['tsite']])

   return pcgr_config_options


def pcgr_error_message(message, logger):
   logger.error('')
   logger.error(message)
   logger.error('')
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
   input_cna_plot_dir = "NA"
   input_conf_dir = "NA"
   panel_normal_vcf_dir = "NA"
   db_dir = "NA"
   base_dir = "NA"
   output_dir_full = "NA"
   panel_normal_vcf_basename = "NA"
   input_vcf_basename = "NA"
   input_cna_basename = "NA"
   input_cna_plot_basename = "NA"
   input_conf_basename = "NA"
   
   ## check that either input vcf or cna segments exist
   if arg_dict['input_vcf'] is None and arg_dict['input_cna'] is None:
      err_msg = "Please specifiy either a VCF input file (--input_vcf) or a copy number segment file (--input_cna)"
      pcgr_error_message(err_msg,logger)
   
   ## check the existence of given output folder
   output_dir_full = os.path.abspath(arg_dict['output_dir'])
   if not os.path.isdir(output_dir_full):
      err_msg = "Output directory (" + str(output_dir_full) + ") does not exist"
      pcgr_error_message(err_msg,logger)
   
   ## check if panel of normal VCF exist
   if not arg_dict['pon_vcf'] is None:
      if not os.path.exists(os.path.abspath(arg_dict['pon_vcf'])):
         err_msg = "Input file (" + str(arg_dict['pon_vcf']) + ") does not exist"
         pcgr_error_message(err_msg,logger)

      if not (os.path.abspath(arg_dict['pon_vcf']).endswith('.vcf.gz')):
         err_msg = "Panel of normals VCF file (" + os.path.abspath(arg_dict['pon_vcf']) + ") does not have the correct file extension (.vcf.gz)"
         pcgr_error_message(err_msg,logger)

      ## check that tabix file exist if bgzipped files is given
      if os.path.abspath(arg_dict['pon_vcf']).endswith('.vcf.gz'):
         tabix_file = arg_dict['pon_vcf'] + '.tbi'
         if not os.path.exists(os.path.abspath(tabix_file)):
            err_msg = "Tabix file (i.e. '.gz.tbi') is not present for the bgzipped panel of normal VCF file (" + os.path.abspath(arg_dict['pon_vcf']) + \
               "). Please make sure your the VCF is properly compressed and indexed (bgzip + tabix)"
            pcgr_error_message(err_msg,logger)

      if arg_dict['input_vcf'] is None:
         warn_msg = "Ignoring panel of normal VCF file, --input_vcf missing"
         pcgr_warn_message(warn_msg, logger)
      else:
         panel_normal_vcf_basename = os.path.basename(str(arg_dict['pon_vcf']))
         panel_normal_vcf_dir = os.path.dirname(os.path.abspath(arg_dict['pon_vcf']))

   ## check if input vcf exist
   if not arg_dict['input_vcf'] is None:
      if not os.path.exists(os.path.abspath(arg_dict['input_vcf'])):
         err_msg = "Input file (" + str(arg_dict['input_vcf']) + ") does not exist"
         pcgr_error_message(err_msg,logger)

      if not (os.path.abspath(arg_dict['input_vcf']).endswith('.vcf') or os.path.abspath(arg_dict['input_vcf']).endswith('.vcf.gz')):
         err_msg = "VCF input file (" + os.path.abspath(arg_dict['input_vcf']) + ") does not have the correct file extension (.vcf or .vcf.gz)"
         pcgr_error_message(err_msg,logger)

      ## check that tabix file exist if bgzipped files is given
      if os.path.abspath(arg_dict['input_vcf']).endswith('.vcf.gz'):
         tabix_file = arg_dict['input_vcf'] + '.tbi'
         if not os.path.exists(os.path.abspath(tabix_file)):
            err_msg = "Tabix file (i.e. '.gz.tbi') is not present for the bgzipped VCF input file (" + os.path.abspath(arg_dict['input_vcf']) + \
               "). Please make sure your input VCF is properly compressed and indexed (bgzip + tabix)"
            pcgr_error_message(err_msg,logger)

      input_vcf_basename = os.path.basename(str(arg_dict['input_vcf']))
      input_vcf_dir = os.path.dirname(os.path.abspath(arg_dict['input_vcf']))

      ## if output vcf exist and overwrite not set
      output_vcf = os.path.join(str(output_dir_full),str(arg_dict['sample_id'])) + '.pcgr_acmg.' + str(arg_dict['genome_assembly'])  + '.vcf.gz'
      if os.path.exists(output_vcf) and arg_dict['force_overwrite'] is False:
         err_msg = "Output files (e.g. " + str(output_vcf) + ") already exist - please specify different sample_id or add option --force_overwrite"
         pcgr_error_message(err_msg,logger)
   
   if not arg_dict['configuration_file'] is None:
      if not os.path.exists(os.path.abspath(arg_dict['configuration_file'])):
         err_msg = "Input file (" + str(arg_dict['configuration_file']) + ") does not exist"
         pcgr_error_message(err_msg,logger)

      if not os.path.abspath(arg_dict['configuration_file']).endswith('.toml'):
         err_msg = "Configuration file (" + os.path.abspath(arg_dict['configuration_file']) + ") does not have the correct file extension (.toml)"
         pcgr_error_message(err_msg,logger)

      input_conf_basename = os.path.basename(str(arg_dict['configuration_file']))
      input_conf_dir = os.path.dirname(os.path.abspath(arg_dict['configuration_file']))
   
   ## check if input cna plot file exist
   # if not arg_dict['input_cna_plot'] is None:
   #    if not os.path.exists(os.path.abspath(arg_dict['input_cna_plot'])):
   #       err_msg = "Input file (" + str(arg_dict['input_cna_plot']) + ") does not exist"
   #       pcgr_error_message(err_msg,logger)
   #    if not (os.path.abspath(arg_dict['input_cna_plot']).endswith('.png')):
   #       err_msg = "CNA segment input file (" + os.path.abspath(arg_dict['input_cna_plot']) + ") does not have the correct file extension (.png)"
   #       pcgr_error_message(err_msg,logger)
   #    if arg_dict['input_cna'] is None:
   #       err_msg = "Input a CNA plot needs to come with a CNA segment file (--input_cna is missing)"
   #       pcgr_error_message(err_msg,logger)
   #    input_cna_plot_basename = os.path.basename(str(arg_dict['input_cna_plot']))
   #    input_cna_plot_dir = os.path.dirname(os.path.abspath(arg_dict['input_cna_plot']))

   ## check if input cna segments exist
   if not arg_dict['input_cna'] is None:
      if not os.path.exists(os.path.abspath(arg_dict['input_cna'])):
         err_msg = "Input file (" + str(arg_dict['input_cna']) + ") does not exist"
         pcgr_error_message(err_msg,logger)
      if not (os.path.abspath(arg_dict['input_cna']).endswith('.tsv') or os.path.abspath(arg_dict['input_cna']).endswith('.txt')):
         err_msg = "CNA segment input file (" + os.path.abspath(arg_dict['input_cna']) + ") does not have the correct file extension (.tsv or .txt)"
         pcgr_error_message(err_msg,logger)
      input_cna_basename = os.path.basename(str(arg_dict['input_cna']))
      input_cna_dir = os.path.dirname(os.path.abspath(arg_dict['input_cna']))

      ## if output cna segments exist and overwrite not set
      output_cna_segments = os.path.join(str(output_dir_full), str(arg_dict['sample_id'])) + '.pcgr_acmg.' + str(arg_dict['genome_assembly']) + '.cna_segments.tsv.gz'
      if os.path.exists(output_cna_segments) and arg_dict['force_overwrite'] is False:
         err_msg = "Output files (e.g. " + str(output_cna_segments) + ") already exist - please specify different sample_id or add option --force_overwrite"
         pcgr_error_message(err_msg,logger)
   
   ## check the existence of base folder
   base_dir = os.path.abspath(arg_dict['pcgr_dir'])
   if not os.path.isdir(base_dir):
      err_msg = "Base directory (" + str(base_dir) + ") does not exist"
      pcgr_error_message(err_msg,logger)
   
   ## check the existence of data folder within the base folder
   db_dir = os.path.join(os.path.abspath(arg_dict['pcgr_dir']), 'data')
   if not os.path.isdir(db_dir):
      err_msg = "Data directory (" + str(db_dir) + ") does not exist"
      pcgr_error_message(err_msg,logger)
   
   ## check the existence of specified assembly data folder within the base folder
   db_assembly_dir = os.path.join(os.path.abspath(arg_dict['pcgr_dir']), 'data', arg_dict['genome_assembly'])
   if not os.path.isdir(db_assembly_dir):
      err_msg = "Data directory for the specified genome assembly (" + str(db_assembly_dir) + ") does not exist"
      pcgr_error_message(err_msg,logger)
   
   ## check the existence of RELEASE_NOTES (starting from 0.4.0)
   rel_notes_file = os.path.join(os.path.abspath(arg_dict['pcgr_dir']), 'data', arg_dict['genome_assembly'], 'RELEASE_NOTES')
   if not os.path.exists(rel_notes_file):
      err_msg = 'The PCGR data bundle is outdated - please download the latest data bundle (see github.com/sigven/pcgr for instructions)'
      pcgr_error_message(err_msg,logger)
      
   f_rel_not = open(rel_notes_file,'r')
   compliant_data_bundle = 0
   for line in f_rel_not:
      if DB_VERSION in line:
         compliant_data_bundle = 1
   
   f_rel_not.close()
    
   if compliant_data_bundle == 0:
      err_msg = 'The PCGR data bundle is not compliant with the software version - please download the latest software and data bundle (see https://github.com/sigven/pcgr for instructions)'
      pcgr_error_message(err_msg,logger)
   
   host_directories = {}
   host_directories['input_vcf_dir_host'] = input_vcf_dir
   host_directories['input_cna_dir_host'] = input_cna_dir
   host_directories['input_cna_plot_dir_host'] = input_cna_plot_dir
   host_directories['input_conf_dir_host'] = input_conf_dir
   host_directories['panel_normal_vcf_dir_host'] = panel_normal_vcf_dir
   host_directories['db_dir_host'] = db_assembly_dir
   host_directories['base_dir_host'] = base_dir
   host_directories['output_dir_host'] = output_dir_full
   host_directories['panel_normal_vcf_basename_host'] = panel_normal_vcf_basename
   host_directories['input_vcf_basename_host'] = input_vcf_basename
   host_directories['input_cna_basename_host'] = input_cna_basename
   host_directories['input_cna_plot_basename_host'] = input_cna_plot_basename
   host_directories['input_conf_basename_host'] = input_conf_basename

   
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

   debug = arg_dict['debug']
   docker_user_id = arg_dict['docker_user_id']
   tumor_only = 0
   cell_line = 0
   vcf_validation = 1
   include_trials = 0
   clinical_trials_set = "OFF"
   msi_prediction_set = "OFF"
   tmb_estimation_set = "OFF"
   msig_estimation_set = "OFF"
   if arg_dict['include_trials']:
      include_trials = 1
      clinical_trials_set = "ON"
   if config_options['msi']['run']:
      msi_prediction_set = "ON"
   if config_options['msigs']['run']:
      msig_estimation_set = "ON"
   if config_options['tmb']['run']:
      tmb_estimation_set = "ON"
   if arg_dict['cell_line']:
      cell_line = 1
   if arg_dict['tumor_only']:
      tumor_only = 1
   if arg_dict['no_vcf_validate']:
      vcf_validation = 0
   ## set basic Docker run commands
   output_vcf = 'None'
   output_pass_vcf = 'None'
   output_pass_tsv = 'None'
   output_maf = 'None'
   uid = ''
   
   global GENCODE_VERSION, VEP_ASSEMBLY, NCBI_BUILD_MAF
   if arg_dict['genome_assembly'] == 'grch37':
      NCBI_BUILD_MAF = 'GRCh37'
      GENCODE_VERSION = 'release 19'
      VEP_ASSEMBLY = 'GRCh37'
   logger = getlogger('pcgr-get-OS')

   if docker_user_id:
      uid = docker_user_id
   elif platform.system() == 'Linux' or platform.system() == 'Darwin' or sys.platform == 'darwin' or sys.platform == 'linux2' or sys.platform == 'linux':
      uid = os.getuid()
   else:
      if platform.system() == 'Windows' or sys.platform == 'win32' or sys.platform == 'cygwin':
         uid = getpass.getuser()

   if uid == '':
      logger.warning('Was not able to get user id/username for logged-in user on the underlying platform (platform.system(): ' + \
         str(platform.system()) + ', sys.platform: ' + str(sys.platform) + '), now running PCGR as root')
      uid = 'root'

   vepdb_dir_host = os.path.join(str(host_directories['db_dir_host']),'.vep')
   input_vcf_docker = 'None'
   input_cna_docker = 'None'
   input_cna_plot_docker = 'None'
   input_conf_docker = 'None'
   panel_normal_docker = 'None'
   docker_cmd_run1 = ''
   docker_cmd_run2 = ''
   docker_cmd_run_end = ''
   ## panel-of-normals annotation
   pon_annotation = 0

   ## Map input files and directories to files/directories within the Docker container
   if DOCKER_IMAGE_VERSION:
      if host_directories['input_vcf_basename_host'] != 'NA':
         input_vcf_docker = '/workdir/input_vcf/' + str(host_directories['input_vcf_basename_host'])
      if host_directories['input_cna_basename_host'] != 'NA':
         input_cna_docker = '/workdir/input_cna/' + str(host_directories['input_cna_basename_host'])
      if host_directories['input_cna_plot_basename_host'] != 'NA':
         input_cna_plot_docker = '/workdir/input_cna_plot/' + str(host_directories['input_cna_plot_basename_host'])
      if host_directories['input_conf_basename_host'] != 'NA':
         input_conf_docker = '/workdir/input_conf/' + str(host_directories['input_conf_basename_host'])
      if host_directories['panel_normal_vcf_basename_host'] != 'NA':
         panel_normal_docker = '/workdir/panel_normal_vcf/' + str(host_directories['panel_normal_vcf_basename_host'])

      vep_volume_mapping = str(vepdb_dir_host) + ":/usr/local/share/vep/data"
      databundle_volume_mapping = str(host_directories['base_dir_host']) + ":/data"
      input_cna_volume_mapping = str(host_directories['input_cna_dir_host']) + ":/workdir/input_cna"
      input_vcf_volume_mapping = str(host_directories['input_vcf_dir_host']) + ":/workdir/input_vcf"
      input_conf_volume_mapping = str(host_directories['input_conf_dir_host']) + ":/workdir/input_conf"
      output_volume_mapping = str(host_directories['output_dir_host']) + ":/workdir/output"
      input_cna_plot_volume_mapping = str(host_directories['input_cna_plot_dir_host']) + ":/workdir/input_cna_plot"
      panel_normal_vcf_volume_mapping = str(host_directories['panel_normal_vcf_dir_host']) + ":/workdir/panel_normal_vcf"

      docker_cmd_run1 = 'NA'

      ## VCF file only
      docker_run_basic = "docker run --rm -t -u " + str(uid)
      docker_cmd_run1 = str(docker_run_basic) + " -v=" +  str(databundle_volume_mapping) + " -v=" + str(vep_volume_mapping) + " -v=" + str(input_conf_volume_mapping) + " -v=" + str(output_volume_mapping)
      if host_directories['input_vcf_dir_host'] != 'NA' and host_directories['input_cna_dir_host'] == 'NA':
         docker_cmd_run1 = docker_cmd_run1  + " -v=" + str(input_vcf_volume_mapping)

      ## CNA file and VCF file provided
      if host_directories['input_vcf_dir_host'] != 'NA' and host_directories['input_cna_dir_host'] != 'NA':
         docker_cmd_run1 = docker_cmd_run1  + " -v=" + str(input_vcf_volume_mapping) + " -v=" + str(input_cna_volume_mapping)

      ## CNA plot provided
      if host_directories['input_cna_plot_dir_host'] != "NA":
         docker_cmd_run1 = docker_cmd_run1 + " -v=" + str(input_cna_plot_volume_mapping)
      
      ## Panel of normal VCFs provided
      if host_directories['panel_normal_vcf_dir_host'] != "NA":
         docker_cmd_run1 = docker_cmd_run1 + " -v=" + str(panel_normal_vcf_volume_mapping)

      docker_cmd_run1 = docker_cmd_run1 + " -w=/workdir/output " + str(DOCKER_IMAGE_VERSION) + " sh -c \""
      
      docker_cmd_run2 = str(docker_run_basic) + " -v=" + str(databundle_volume_mapping) + " -v=" + str(output_volume_mapping)
      if host_directories['panel_normal_vcf_dir_host'] != "NA":
         docker_cmd_run2 = docker_cmd_run2 + " -v=" + str(panel_normal_vcf_volume_mapping)
      docker_cmd_run2 = docker_cmd_run2 + " -w=/workdir/output " + str(DOCKER_IMAGE_VERSION) + " sh -c \""
      docker_cmd_run_end = '\"'

      data_dir = '/data'
      output_dir = '/workdir/output'
      vep_dir = '/usr/local/share/vep/data'
      r_scripts_dir = '/'

   ## If running non-Dockerized - specifiy paths for input files and directories
   else:
      if host_directories['input_vcf_basename_host'] != 'NA':
         input_vcf_docker = os.path.join(host_directories['input_vcf_dir_host'], host_directories['input_vcf_basename_host'])
      if host_directories['input_cna_basename_host'] != 'NA':
         input_cna_docker = os.path.join(host_directories['input_cna_dir_host'], host_directories['input_cna_basename_host'])
      if host_directories['input_cna_plot_basename_host'] != 'NA':
         input_cna_plot_docker = os.path.join(host_directories['input_cna_plot_dir_host'], host_directories['input_cna_plot_basename_host'])
      if host_directories['input_conf_basename_host'] != 'NA':
         input_conf_docker = os.path.join(host_directories['input_conf_dir_host'], host_directories['input_conf_basename_host'])
      if host_directories['panel_normal_vcf_basename_host'] != 'NA':
         panel_normal_docker = os.path.join(host_directories['panel_normal_vcf_dir_host'], host_directories['panel_normal_vcf_basename_host'])

      data_dir = host_directories['base_dir_host']
      output_dir = host_directories['output_dir_host']
      vep_dir = vepdb_dir_host
      r_scripts_dir = ''

   check_subprocess(logger, docker_cmd_run1.replace("-u " + str(uid), "") + 'mkdir -p ' + output_dir + docker_cmd_run_end)


   ## PCGR|validate_input - verify that VCF and CNA segment file is of appropriate format
   logger = getlogger('pcgr-validate-arguments-input')
   vcf_validate_command = docker_cmd_run1 + "pcgr_validate_input.py " + data_dir + " " + str(input_vcf_docker) + " " + \
      str(input_cna_docker) + " " + str(input_conf_docker) + " " + str(panel_normal_docker) + " " + str(vcf_validation) + \
         " " + str(tumor_only) + " " + str(arg_dict['genome_assembly'])
   if not DOCKER_IMAGE_VERSION:
      vcf_validate_command += ' --output_dir ' + output_dir + docker_cmd_run_end
   else:
      vcf_validate_command += docker_cmd_run_end
   check_subprocess(logger, vcf_validate_command)
   logger.info('Finished')

   ## PCGR|start - Log key information about sample, options and sequencing assay/design
   logger = getlogger("pcgr-start")
   print()
   assay_mode = "Tumor vs. Control"
   if tumor_only == 1:
      assay_mode = "Tumor-Only"
      if cell_line == 1:
         assay_mode = "Tumor-Only (cell line)"
   logger.info("--- Personal Cancer Genome Reporter workflow ----")
   logger.info("Sample name: " + str(arg_dict['sample_id']))
   if config_options['tumor_type']['type'] == "Cancer_NOS":
      logger.info("Tumor type: Cancer_NOS (Any tumortype)")
   else:
      logger.info("Tumor type: " + str(config_options['tumor_type']['type']))
   logger.info("Sequencing assay - type: " + str(arg_dict['assay']))
   logger.info("Sequencing assay - mode: " + str(assay_mode))
   logger.info("Sequencing assay - coding target size: " + str(config_options['tmb']['target_size_mb']) + "Mb")
   logger.info("Genome assembly: " + str(arg_dict['genome_assembly']))
   logger.info("Mutational signature estimation: " + str(msig_estimation_set))
   logger.info("MSI classification: " + str(msi_prediction_set))
   logger.info("Mutational burden estimation: " + str(tmb_estimation_set))
   logger.info("Include molecularly targeted clinical trials (beta): " + str(clinical_trials_set))

   if not input_vcf_docker == 'None':

      ## Define temporary output file names
      output_vcf = os.path.join(output_dir, str(arg_dict['sample_id']) + '.pcgr_acmg.'  + str(arg_dict['genome_assembly']) + '.vcf.gz')
      output_pass_vcf = os.path.join(output_dir, str(arg_dict['sample_id']) + '.pcgr_acmg.'  + str(arg_dict['genome_assembly']) + '.pass.vcf.gz')
      output_pass_tsv = os.path.join(output_dir, str(arg_dict['sample_id']) + '.pcgr_acmg.'  + str(arg_dict['genome_assembly']) + '.pass.tsv')
      output_maf = os.path.join(output_dir, str(arg_dict['sample_id']) + '.pcgr_acmg.' + str(arg_dict['genome_assembly']) + '.tmp.maf')
      output_vcf2maf_log = os.path.join(output_dir, str(arg_dict['sample_id']) + '.pcgr_acmg.' + str(arg_dict['genome_assembly']) + '.maf.log')
      input_vcf_pcgr_ready = os.path.join(output_dir, re.sub(r'(\.vcf$|\.vcf\.gz$)', '.pcgr_ready.vcf.gz', host_directories['input_vcf_basename_host']))
      input_vcf_pcgr_ready_uncompressed = os.path.join(output_dir, re.sub(r'(\.vcf$|\.vcf\.gz$)', '.pcgr_ready.vcf', host_directories['input_vcf_basename_host']))
      vep_vcf = re.sub(r'(\.vcf$|\.vcf\.gz$)', '.vep.vcf', input_vcf_pcgr_ready)
      vep_vcfanno_vcf = re.sub(r'(\.vcf$|\.vcf\.gz$)', '.vep.vcfanno.vcf', input_vcf_pcgr_ready)
      vep_vcfanno_annotated_vcf = re.sub(r'\.vcfanno', '.vcfanno.annotated', vep_vcfanno_vcf) + '.gz'
      vep_vcfanno_annotated_pass_vcf = re.sub(r'\.vcfanno', '.vcfanno.annotated.pass', vep_vcfanno_vcf) + '.gz'

      ## Path for human genome assembly file (FASTA)
      fasta_assembly = os.path.join(vep_dir, "homo_sapiens", str(VEP_VERSION) + "_" + str(VEP_ASSEMBLY), "Homo_sapiens." + str(VEP_ASSEMBLY) + ".dna.primary_assembly.fa.gz")
      
      ## List all VEP flags used when calling VEP
      vep_flags = "--hgvs --af --af_1kg --af_gnomad --variant_class --domains --symbol --protein --ccds --mane " + \
         "--uniprot --appris --biotype --tsl --canonical --gencode_basic --cache --numbers --total_length --allele_number " + \
         "--no_stats --no_escape --xref_refseq --vcf --check_ref --dont_skip --flag_pick_allele --plugin NearestExonJB,max_range=50000"
      vep_options = "--pick_order " + str(config_options['other']['vep_pick_order']) + " --force_overwrite --species homo_sapiens --assembly " + \
         str(VEP_ASSEMBLY) + " --offline --fork " + str(config_options['other']['n_vep_forks']) + " " + str(vep_flags)  + " --dir " + vep_dir
      vep_options += " --cache_version " + str(VEP_VERSION)
      if config_options['other']['vep_skip_intergenic'] == 1:
         vep_options = vep_options + " --no_intergenic"
      if not debug:
         vep_options += " --quiet"

      ## Compose full VEP command
      vep_main_command = docker_cmd_run1 + "vep --input_file " + str(input_vcf_pcgr_ready) + " --output_file " + str(vep_vcf) + " " + str(vep_options) + " --fasta " + str(fasta_assembly) + docker_cmd_run_end
      vep_bgzip_command = docker_cmd_run1 + "bgzip -f -c " + str(vep_vcf) + " > " + str(vep_vcf) + '.gz' + docker_cmd_run_end
      vep_tabix_command = docker_cmd_run1 + "tabix -f -p vcf " + str(vep_vcf) + ".gz" + docker_cmd_run_end

      ## PCGR|VEP - run consequence annotation with Variant Effect Predictor
      print()
      logger = getlogger('pcgr-vep')
      logger.info("PCGR - STEP 1: Basic variant annotation with Variant Effect Predictor (" + str(VEP_VERSION) + ", GENCODE " + str(GENCODE_VERSION) + ", " + str(arg_dict['genome_assembly']) + ")")
      logger.info("VEP configuration - one primary consequence block pr. alternative allele (--flag_pick_allele)")
      logger.info("VEP configuration - transcript pick order: " + str(config_options['other']['vep_pick_order']))
      logger.info("VEP configuration - transcript pick order: See more at https://www.ensembl.org/info/docs/tools/vep/script/vep_other.html#pick_options")
      logger.info("VEP configuration - skip intergenic: " + str(config_options['other']['vep_skip_intergenic']))
      check_subprocess(logger, vep_main_command)
      check_subprocess(logger, vep_bgzip_command)
      check_subprocess(logger, vep_tabix_command)

      ## PCGR|vcf2maf - if option set, convert VCF to MAF with https://github.com/mskcc/vcf2maf
      if config_options['other']['vcf2maf'] == 1:
         logger.info('Converting VEP-annotated VCF to MAF with https://github.com/mskcc/vcf2maf')
         vcf2maf_command = str(docker_cmd_run1) + "vcf2maf.pl --input-vcf " + str(input_vcf_pcgr_ready_uncompressed) + " --tumor-id " + \
            str(arg_dict['sample_id']) + " --output-maf " + str(output_maf) + " --ref-fasta " + str(fasta_assembly) + " --filter-vcf 0 --ncbi-build " + \
               str(NCBI_BUILD_MAF) + " > " + str(output_vcf2maf_log) + " 2>&1" + docker_cmd_run_end
         clean_vcf2maf_command = str(docker_cmd_run1) + "rm -f " + str(output_vcf2maf_log) + " " + re.sub(r'(\.vcf$)', '.vep.vcf', input_vcf_pcgr_ready_uncompressed) + " " + docker_cmd_run_end
         check_subprocess(logger, vcf2maf_command)
         check_subprocess(logger, clean_vcf2maf_command)
      logger.info("Finished")

      ## PCGR|vcfanno - annotate VCF against a number of variant annotation resources
      print()
      logger = getlogger('pcgr-vcfanno')
      pcgr_vcfanno_command = str(docker_cmd_run2) + "pcgr_vcfanno.py --num_processes " + str(config_options['other']['n_vcfanno_proc']) + \
         " --chasmplus --dbnsfp --docm --clinvar --icgc --civic --cgi --tcga_pcdm --winmsk --simplerepeats --tcga --uniprot --cancer_hotspots --pcgr_onco_xref " + \
            str(vep_vcf) + ".gz " + str(vep_vcfanno_vcf) + " " + os.path.join(data_dir, "data", str(arg_dict['genome_assembly']))
      if debug:
         pcgr_vcfanno_command += " --debug"
      if panel_normal_docker != 'None':
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
      pcgr_summarise_command = str(docker_cmd_run2) + "pcgr_summarise.py " + str(vep_vcfanno_vcf) + ".gz " + str(pon_annotation) + " " + str(os.path.join(data_dir, "data", str(arg_dict['genome_assembly']))) + docker_cmd_run_end
      if debug:
         pcgr_summarise_command  += ' --debug'
      logger.info("PCGR - STEP 3: Cancer gene annotations with pcgr-summarise")
      check_subprocess(logger, pcgr_summarise_command)

      ## PCGR|clean - move output files and clean up temporary files
      create_output_vcf_command1 = str(docker_cmd_run2) + 'mv ' + str(vep_vcfanno_annotated_vcf) + ' ' + str(output_vcf) + docker_cmd_run_end
      create_output_vcf_command2 = str(docker_cmd_run2) + 'mv ' + str(vep_vcfanno_annotated_vcf) + '.tbi ' + str(output_vcf) + '.tbi' + docker_cmd_run_end
      create_output_vcf_command3 = str(docker_cmd_run2) + 'mv ' + str(vep_vcfanno_annotated_pass_vcf) + ' ' + str(output_pass_vcf) + docker_cmd_run_end
      create_output_vcf_command4 = str(docker_cmd_run2) + 'mv ' + str(vep_vcfanno_annotated_pass_vcf) + '.tbi ' + str(output_pass_vcf) + '.tbi' + docker_cmd_run_end
      clean_command = str(docker_cmd_run2) + 'rm -f ' + str(vep_vcf) + '* ' + str(vep_vcfanno_annotated_vcf) + ' ' + \
         str(vep_vcfanno_annotated_pass_vcf) + '* ' + str(vep_vcfanno_vcf) + '* ' +  str(input_vcf_pcgr_ready_uncompressed) + "* "  + docker_cmd_run_end
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
   if not arg_dict['basic']:
      ttype = config_options['tumor_type']['type'].replace(" ","_").replace("/","@")
      logger = getlogger('pcgr-writer')
      logger.info("PCGR - STEP 4: Generation of output files - variant interpretation report for precision oncology")
      pcgr_report_command = (docker_cmd_run1 + os.path.join(r_scripts_dir, "pcgr.R") + " " + output_dir + " " + str(output_pass_tsv) + ".gz" + " " + \
                           input_cna_docker + " " + str(arg_dict['sample_id']) + " " + input_conf_docker + " " + str(PCGR_VERSION) + " " +  \
                           str(arg_dict['genome_assembly']) + " " + data_dir + " " + \
                           str(input_cna_plot_docker) + " " + str(config_options['tumor_purity']) + " " + \
                           str(config_options['tumor_ploidy']) + " " + str(config_options['assay']) + " " + str(tumor_only) +  " " + \
                           str(config_options['tmb']['run']) + " " + str(config_options['tmb']['algorithm']) + " " + \
                           str(config_options['msi']['run']) + " " + str(config_options['msigs']['run']) + " " + \
                           str(config_options['tmb']['target_size_mb']) + " " + str(config_options['cna']['logR_homdel']) + " " + \
                           str(config_options['cna']['logR_gain']) + " " + str(config_options['cna']['cna_overlap_pct']) + " "  + \
                           str(config_options['msigs']['mutation_limit']) + " " + str(config_options['msigs']['all_reference_signatures']) + " " + \
                           str(config_options['allelic_support']['tumor_af_min']) + " " + str(config_options['allelic_support']['tumor_dp_min']) + " "  + \
                           str(config_options['allelic_support']['control_af_max']) + " " + str(config_options['allelic_support']['control_dp_min']) + " " + \
                           str(cell_line) + " " + str(include_trials) + " " + str(ttype) + docker_cmd_run_end)
      check_subprocess(logger, pcgr_report_command)
      logger.info("Finished")

   print()



if __name__=="__main__": __main__()

