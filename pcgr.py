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

version = '0.6.0'

def __main__():
   
   parser = argparse.ArgumentParser(description='Personal Cancer Genome Reporter (PCGR) workflow for clinical interpretation of somatic nucleotide variants and copy number aberration segments',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
   parser.add_argument('--input_vcf', dest = "input_vcf", help='VCF input file with somatic query variants (SNVs/InDels).')
   parser.add_argument('--input_cna', dest = "input_cna",help='Somatic copy number alteration segments (tab-separated values)')
   parser.add_argument('--force_overwrite', action = "store_true", help='By default, the script will fail with an error if any output file already exists. You can force the overwrite of existing result files by using this flag')
   parser.add_argument('--version', action='version', version='%(prog)s ' + str(version))
   parser.add_argument('--basic',action="store_true",help="Run functional variant annotation on VCF through VEP/vcfanno, omit other analyses (i.e. CNA, MSI, report generation etc. (STEP 4)")
   parser.add_argument('pcgr_dir',help='PCGR base directory with accompanying data directory, e.g. ~/pcgr-0.6.0')
   parser.add_argument('output_dir',help='Output directory')
   parser.add_argument('genome_assembly',choices = ['grch37','grch38'], help='Genome assembly build: grch37 or grch38')
   parser.add_argument('configuration_file',help='PCGR configuration file (TOML format)')
   parser.add_argument('sample_id',help="Tumor sample/cancer genome identifier - prefix for output files")
   
   docker_image_version = 'sigven/pcgr:' + str(version)
   args = parser.parse_args()
   
   overwrite = 0
   if args.force_overwrite is True:
      overwrite = 1
   
   # check that script and Docker image version correspond
   check_docker_command = 'docker images -q ' + str(docker_image_version)
   output = subprocess.check_output(str(check_docker_command), stderr=subprocess.STDOUT, shell=True)
   logger = getlogger('pcgr-validate-config')
   
   if(len(output) == 0):
      err_msg = 'Docker image ' + str(docker_image_version) + ' does not exist, pull image from Dockerhub (docker pull ' + str(docker_image_version) + ')'
      pcgr_error_message(err_msg,logger)
   
   config_options = {}
   if os.path.exists(args.configuration_file):
      config_options = read_config_options(args.configuration_file, args.pcgr_dir, args.genome_assembly, logger)
   else:
      err_msg = "PCGR configuration file " + str(args.configuration_file) + " does not exist - exiting"
      pcgr_error_message(err_msg,logger)

   logger = getlogger('pcgr-check-files')
   host_directories = verify_input_files(args.input_vcf, args.input_cna, args.configuration_file, config_options, args.pcgr_dir, args.output_dir, args.sample_id, args.genome_assembly, overwrite, logger)

   run_pcgr(host_directories, docker_image_version, config_options, args.sample_id, args.genome_assembly, version, args.basic)


def read_config_options(configuration_file, pcgr_dir, genome_assembly, logger):
   
   ## read default options
   pcgr_config_options = {}
   pcgr_configuration_file_default = os.path.join(pcgr_dir,'data',str(genome_assembly),'pcgr_configuration_somatic_default.toml')
   if not os.path.exists(pcgr_configuration_file_default):
      err_msg = "Default PCGR configuration file " + str(pcgr_configuration_file_default) + " does not exist - exiting"
      pcgr_error_message(err_msg,logger)
   try:
      pcgr_config_options = toml.load(pcgr_configuration_file_default)
   except (IndexError,TypeError):
      err_msg = 'Configuration file ' + str(configuration_file) + ' is not formatted correctly'
      pcgr_error_message(err_msg, logger)

   ## override with options set by the users
   try:
      user_options = toml.load(configuration_file)
   except (IndexError,TypeError):
      err_msg = 'Configuration file ' + str(configuration_file) + ' is not formatted correctly'
      pcgr_error_message(err_msg, logger)

   tumor_types = []
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
            if section == 'tumor_type':
               if user_options[section][var]:
                  tumor_types.append(var)
            tier_options = ['pcgr','pcgr_acmg']
            normalization_options = ['default','exome','genome','exome2genome']
            theme_options = ['default', 'cerulean', 'journal', 'flatly', 'readable', 'spacelab', 'united', 'cosmo', 'lumen', 'paper', 'sandstone', 'simplex','yeti']
            if var == 'mutsignatures_normalization' and not str(user_options[section][var]) in normalization_options:
               err_msg = 'Configuration value ' + str(user_options[section][var]) + ' for ' + str(var) + ' cannot be parsed properly (expecting \'default\', \'exome\', \'genome\', or \'exome2genome\')'
               pcgr_error_message(err_msg, logger)
            if var == 'mutsignatures_cutoff' and (float(user_options[section][var]) > 1 or float(user_options[section][var]) < 0) :
               err_msg = 'Configuration value ' + str(user_options[section][var]) + " must be within the [0,1] range"
               pcgr_error_message(err_msg, logger)
            if var == 'mutsignatures_signature_limit':
               if int(user_options[section][var]) < 1 or int(user_options[section][var]) > 30:
                  err_msg = "Number of mutational signatures in search space ('mutsignatures_signature_limit') must be positive and not more than 30 (retrieved value: " + str(user_options[section][var]) + ")"
                  pcgr_error_message(err_msg,logger)
            if var == 'tier_model' and not str(user_options[section][var]) in tier_options:
               err_msg = 'Configuration value ' + str(user_options[section][var]) + ' for ' + str(var) + ' cannot be parsed properly (expecting \'pcgr\', or \'pcgr_acmg\')'
               pcgr_error_message(err_msg, logger)
            if var == 'report_theme' and not str(user_options[section][var]) in theme_options:
               err_msg = 'Configuration value ' + str(user_options[section][var]) + ' for ' + str(var) + ' cannot be parsed properly (expecting \'default\', \'cerulean\', \'journal\', \'flatly\', \'readable\', \'spacelab\', \'united\', \'cosmo\', \'lumen\', \'paper\', \'sandstone\', \'simplex\',or \'yeti\')'
               pcgr_error_message(err_msg, logger)
            if var.startswith('maf_'):
               if user_options['tumor_only'][var] < 0 or user_options[section][var] > 1:
                  err_msg = "MAF value: " + str(var) + " must be within the [0,1] range, current value is " + str(user_options[section][var]) + ")"
                  pcgr_error_message(err_msg,logger)
            if var == 'min_af_tumor':
               if user_options['allelic_support'][var] < 0 or user_options[section][var] > 1:
                  err_msg = "Minimum AF tumor: " + str(var) + " must be within the [0,1] range, current value is " + str(user_options[section][var]) + ")"
                  pcgr_error_message(err_msg,logger)
            if var == 'max_af_normal':
               if user_options['allelic_support'][var] < 0 or user_options[section][var] > 1:
                  err_msg = "Maximum AF normal: " + str(var) + " must be within the [0,1] range, current value is " + str(user_options[section][var]) + ")"
                  pcgr_error_message(err_msg,logger)
            if var == 'target_size_mb':
               if user_options[section][var] < 0 or user_options['mutational_burden'][var] > 50:
                  err_msg = "Target size region in Mb (" + str(user_options[section][var]) + ") is not positive or larger than the likely maximum size of the coding human genome (50 Mb))"
                  pcgr_error_message(err_msg,logger)
               if user_options[section][var] < 1:
                  warn_msg = "Target size region in Mb (" + str(user_options[section][var]) + ") must be greater than 1 for mutational burden estimate to be robust (ignoring TMB calculation)"
                  pcgr_warn_message(warn_msg,logger)
                  pcgr_config_options[section]['mutational_burden'] = False
            if var == 'logR_homdel':
               if user_options['cna'][var] > 0:
                  err_msg = "Log ratio for homozygous deletions (" + str(user_options[section][var]) + ") should be less than zero"
                  pcgr_error_message(err_msg,logger)
            if var == 'cna_overlap_pct':
               if user_options['cna'][var] > 100 or user_options['cna'][var] <= 0:
                  err_msg = "Minimum percent overlap between copy number segment and gene transcript (" + str(user_options[section][var]) + ") should be greater than zero and less than 100"
                  pcgr_error_message(err_msg,logger)
            if var == 'logR_gain':
               if user_options['cna'][var] < 0:
                  err_msg = "Log ratio for copy number amplifications (" + str(user_options[section][var]) + ") should be greater than zero"
                  pcgr_error_message(err_msg,logger)
            pcgr_config_options[section][var] = user_options[section][var]

   if len(tumor_types) > 2:
      err_msg = "Two many tumor types (", str(",".join(tumor_types)) + ")  set to True - limit is set to two"
      pcgr_error_message(err_msg,logger)
   if pcgr_config_options['msi']['msi'] == 1 and pcgr_config_options['mutational_burden']['mutational_burden'] == 0:
      err_msg = "Prediction of MSI status (msi = true) requires mutational burden/target_size input (mutational_burden = true) - this is currently set as false"
      pcgr_error_message(err_msg,logger)
   if pcgr_config_options['tumor_only']['vcf_tumor_only'] == 1:
      if pcgr_config_options['msi']['msi'] == 1:
         warn_msg = 'Prediction of MSI status is not perfomed in tumor-only mode (vcf_tumor_only = true)'
         pcgr_warn_message(warn_msg,logger)
      if pcgr_config_options['mutational_burden']['mutational_burden'] == 1:
         warn_msg = 'Estimation of mutational burden is not performed in tumor-only mode (vcf_tumor_only = true)'
         pcgr_warn_message(warn_msg,logger)
      if pcgr_config_options['mutational_signatures']['mutsignatures'] == 1:
         warn_msg = 'Estimation of mutational signatures is not perfomed in tumor-only mode (vcf_tumor_only = true)'
         pcgr_warn_message(warn_msg,logger)
   #print(str(pcgr_config_options))
   return pcgr_config_options


def pcgr_error_message(message, logger):
   logger.error('')
   logger.error(message)
   logger.error('')
   exit(0)

def pcgr_warn_message(message, logger):
   logger.warning(message)

def verify_input_files(input_vcf, input_cna, configuration_file, pcgr_config_options, base_pcgr_dir, output_dir, sample_id, genome_assembly, overwrite, logger):
   """
   Function that checks the input files and directories provided by the user and checks for their existence
   """
 
   input_vcf_dir = "NA"
   input_cna_dir = "NA"
   input_conf_dir = "NA"
   db_dir = "NA"
   base_dir = "NA"
   output_dir_full = "NA"
   input_vcf_basename = "NA"
   input_cna_basename = "NA"
   input_conf_basename = "NA"
   
   ## check that either input vcf or cna segments exist
   if input_vcf is None and input_cna is None:
      err_msg = "Please specifiy either a VCF input file (--input_vcf) or a copy number segment file (--input_cna)"
      pcgr_error_message(err_msg,logger)
   
   ## if msi or mutsignatures is set to True, input VCF is needed
   if input_vcf is None:
      if pcgr_config_options['tumor_only']['vcf_tumor_only'] == 1:
         err_msg = "Tumor-only mode ('vcf_tumor_only = true') requires a VCF input file (--input_vcf) - this is currently missing"
         pcgr_error_message(err_msg,logger)
      if pcgr_config_options['msi']['msi'] == 1:
         err_msg = "Prediction of MSI status (msi = true) requires a VCF input file (--input_vcf) - this is currently missing"
         pcgr_error_message(err_msg,logger)
      if pcgr_config_options['mutational_burden']['mutational_burden'] == 1:
         err_msg = "Calculation of mutational burden (mutational_burden = true) requires a VCF input file (--input_vcf) - this is currently missing"
         pcgr_error_message(err_msg,logger)
   
      if pcgr_config_options['mutational_signatures']['mutsignatures'] == 1:
         err_msg = "Identification of mutational signatures (mutsignatures = true) requires a VCF input file (--input_vcf) - this is currently missing"
         pcgr_error_message(err_msg,logger)
   
   
   ## check the existence of given output folder
   output_dir_full = os.path.abspath(output_dir)
   if not os.path.isdir(output_dir_full):
      err_msg = "Output directory (" + str(output_dir_full) + ") does not exist"
      pcgr_error_message(err_msg,logger)
   
   ## check if input vcf exist
   if not input_vcf is None:
      if not os.path.exists(os.path.abspath(input_vcf)):
         err_msg = "Input file (" + str(input_vcf) + ") does not exist"
         pcgr_error_message(err_msg,logger)

      if not (os.path.abspath(input_vcf).endswith('.vcf') or os.path.abspath(input_vcf).endswith('.vcf.gz')):
         err_msg = "VCF input file (" + os.path.abspath(input_vcf) + ") does not have the correct file extension (.vcf or .vcf.gz)"
         pcgr_error_message(err_msg,logger)

      ## check that tabix file exist if bgzipped files is given
      if os.path.abspath(input_vcf).endswith('.vcf.gz'):
         tabix_file = input_vcf + '.tbi'
         if not os.path.exists(os.path.abspath(tabix_file)):
            err_msg = "Tabix file (i.e. '.gz.tbi') is not present for the bgzipped VCF input file (" + os.path.abspath(input_vcf) + "). Please make sure your input VCF is properly compressed and indexed (bgzip + tabix)"
            pcgr_error_message(err_msg,logger)

      input_vcf_basename = os.path.basename(str(input_vcf))
      input_vcf_dir = os.path.dirname(os.path.abspath(input_vcf))

      ## if output vcf exist and overwrite not set
      output_vcf = os.path.join(str(output_dir_full),str(sample_id)) + '.' + str(pcgr_config_options['tier_model']['tier_model']) + '.vcf.gz'
      if os.path.exists(output_vcf) and overwrite == 0:
         err_msg = "Output files (e.g. " + str(output_vcf) + ") already exist - please specify different sample_id or add option --force_overwrite"
         pcgr_error_message(err_msg,logger)
   
   if not configuration_file is None:
      if not os.path.exists(os.path.abspath(configuration_file)):
         err_msg = "Input file (" + str(configuration_file) + ") does not exist"
         pcgr_error_message(err_msg,logger)

      if not os.path.abspath(configuration_file).endswith('.toml'):
         err_msg = "Configuration file (" + os.path.abspath(configuration_file) + ") does not have the correct file extension (.toml)"
         pcgr_error_message(err_msg,logger)

      input_conf_basename = os.path.basename(str(configuration_file))
      input_conf_dir = os.path.dirname(os.path.abspath(configuration_file))
   
   
   ## check if input cna segments exist
   if not input_cna is None:
      if not os.path.exists(os.path.abspath(input_cna)):
         err_msg = "Input file (" + str(input_cna) + ") does not exist"
         pcgr_error_message(err_msg,logger)
      if not (os.path.abspath(input_cna).endswith('.tsv') or os.path.abspath(input_cna).endswith('.txt')):
         err_msg = "CNA segment input file (" + os.path.abspath(input_cna) + ") does not have the correct file extension (.tsv or .txt)"
         pcgr_error_message(err_msg,logger)
      input_cna_basename = os.path.basename(str(input_cna))
      input_cna_dir = os.path.dirname(os.path.abspath(input_cna))

      ## if output cna segments exist and overwrite not set
      output_cna_segments = os.path.join(str(output_dir_full),str(sample_id)) + '.' + str(pcgr_config_options['tier_model']['tier_model']) + '.cna_segments.tsv.gz'
      if os.path.exists(output_cna_segments) and overwrite == 0:
         err_msg = "Output files (e.g. " + str(output_cna_segments) + ") already exist - please specify different sample_id or add option --force_overwrite"
         pcgr_error_message(err_msg,logger)
   
   ## check the existence of base folder
   base_dir = os.path.abspath(base_pcgr_dir)
   if not os.path.isdir(base_dir):
      err_msg = "Base directory (" + str(base_dir) + ") does not exist"
      pcgr_error_message(err_msg,logger)
   
   ## check the existence of data folder within the base folder
   db_dir = os.path.join(os.path.abspath(base_pcgr_dir),'data')
   if not os.path.isdir(db_dir):
      err_msg = "Data directory (" + str(db_dir) + ") does not exist"
      pcgr_error_message(err_msg,logger)
   
   ## check the existence of specified assembly data folder within the base folder
   db_assembly_dir = os.path.join(os.path.abspath(base_pcgr_dir),'data',genome_assembly)
   if not os.path.isdir(db_assembly_dir):
      err_msg = "Data directory for the specified genome assembly (" + str(db_assembly_dir) + ") does not exist"
      pcgr_error_message(err_msg,logger)
   
   ## check the existence of RELEASE_NOTES (starting from 0.4.0)
   rel_notes_file = os.path.join(os.path.abspath(base_pcgr_dir),'data',genome_assembly,'RELEASE_NOTES')
   if not os.path.exists(rel_notes_file):
      err_msg = 'The PCGR data bundle is outdated - please download the latest data bundle (see github.com/sigven/pcgr for instructions)'
      pcgr_error_message(err_msg,logger)
      
   f_rel_not = open(rel_notes_file,'r')
   compliant_data_bundle = 0
   for line in f_rel_not:
      version_check = 'PCGR_DB_VERSION = 20180422'
      if version_check in line:
         compliant_data_bundle = 1
   
   f_rel_not.close()
    
   if compliant_data_bundle == 0:
      err_msg = 'The PCGR data bundle is not compliant with the software version - please download the latest software and data bundle (see https://github.com/sigven/pcgr for instructions)'
      pcgr_error_message(err_msg,logger)
   
   host_directories = {}
   host_directories['input_vcf_dir_host'] = input_vcf_dir
   host_directories['input_cna_dir_host'] = input_cna_dir
   host_directories['input_conf_dir_host'] = input_conf_dir
   host_directories['db_dir_host'] = db_assembly_dir
   host_directories['base_dir_host'] = base_dir
   host_directories['output_dir_host'] = output_dir_full
   host_directories['input_vcf_basename_host'] = input_vcf_basename
   host_directories['input_cna_basename_host'] = input_cna_basename
   host_directories['input_conf_basename_host'] = input_conf_basename

   
   return host_directories
   

def check_subprocess(command):
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

def run_pcgr(host_directories, docker_image_version, config_options, sample_id, genome_assembly, version, basic):
   """
   Main function to run the PCGR workflow using Docker
   """
   
   ## set basic Docker run commands
   output_vcf = 'None'
   output_pass_vcf = 'None'
   output_pass_tsv = 'None'
   uid = ''
   vep_version = '92'
   gencode_version = 'release 28'
   ncbi_build_maf = 'GRCh38'
   if genome_assembly == 'grch37':
      ncbi_build_maf = 'GRCh37'
      gencode_version = 'release 19'
   logger = getlogger('pcgr-get-OS')
   if platform.system() == 'Linux' or platform.system() == 'Darwin' or sys.platform == 'darwin' or sys.platform == 'linux2' or sys.platform == 'linux':
      uid = os.getuid()
   else:
      if platform.system() == 'Windows' or sys.platform == 'win32' or sys.platform == 'cygwin':
         uid = getpass.getuser()
   
   if uid == '':
      logger.warning('Was not able to get user id/username for logged-in user on the underlying platform (platform.system(): ' + str(platform.system()) + ', sys.platform: ' + str(sys.platform) + '), now running PCGR as root')
      uid = 'root'
   
   vepdb_dir_host = os.path.join(str(host_directories['db_dir_host']),'.vep')

   input_vcf_docker = 'None'
   input_cna_docker = 'None'
   input_conf_docker = 'None'
   
   if host_directories['input_vcf_basename_host'] != 'NA':
      input_vcf_docker = '/workdir/input_vcf/' + str(host_directories['input_vcf_basename_host'])
   if host_directories['input_cna_basename_host'] != 'NA':
      input_cna_docker = '/workdir/input_cna/' + str(host_directories['input_cna_basename_host'])
   if host_directories['input_conf_basename_host'] != 'NA':
      input_conf_docker = '/workdir/input_conf/' + str(host_directories['input_conf_basename_host'])

   docker_command_run1 = 'NA'
   if host_directories['input_vcf_dir_host'] != 'NA' and host_directories['input_cna_dir_host'] == 'NA':
      docker_command_run1 = "docker run --rm -t -u " + str(uid) + " -v=" + str(host_directories['base_dir_host']) + ":/data -v=" + str(vepdb_dir_host) + ":/usr/local/share/vep/data -v=" + str(host_directories['input_vcf_dir_host']) + ":/workdir/input_vcf -v=" + str(host_directories['output_dir_host']) + ":/workdir/output -w=/workdir/output " + str(docker_image_version) + " sh -c \""
      
      if host_directories['input_conf_dir_host'] != 'NA':
         docker_command_run1 = "docker run --rm -t -u " + str(uid) + " -v=" + str(host_directories['base_dir_host']) + ":/data -v=" + str(vepdb_dir_host) + ":/usr/local/share/vep/data -v=" + str(host_directories['input_vcf_dir_host']) + ":/workdir/input_vcf -v=" + str(host_directories['input_conf_dir_host']) + ":/workdir/input_conf -v=" + str(host_directories['output_dir_host']) + ":/workdir/output -w=/workdir/output " + str(docker_image_version) + " sh -c \""

   if host_directories['input_vcf_dir_host'] == 'NA' and host_directories['input_cna_dir_host'] != 'NA':
      docker_command_run1 = "docker run --rm -t -u " + str(uid) + " -v=" + str(host_directories['base_dir_host']) + ":/data -v=" + str(vepdb_dir_host) + ":/usr/local/share/vep/data -v=" + str(host_directories['input_cna_dir_host']) + ":/workdir/input_cna -v=" + str(host_directories['input_conf_dir_host']) + ":/workdir/input_conf -v=" + str(host_directories['output_dir_host']) + ":/workdir/output -w=/workdir/output " + str(docker_image_version) + " sh -c \""

   if host_directories['input_vcf_dir_host'] != 'NA' and host_directories['input_cna_dir_host'] != 'NA':
      docker_command_run1 = "docker run --rm -t -u " + str(uid) + " -v=" + str(host_directories['base_dir_host']) + ":/data -v=" + str(vepdb_dir_host) + ":/usr/local/share/vep/data -v=" + str(host_directories['input_cna_dir_host']) + ":/workdir/input_cna -v=" + str(host_directories['input_vcf_dir_host']) + ":/workdir/input_vcf -v=" + str(host_directories['output_dir_host']) + ":/workdir/output -w=/workdir/output " + str(docker_image_version) + " sh -c \""
      
      if host_directories['input_conf_dir_host'] != 'NA':
         docker_command_run1 = "docker run --rm -t -u " + str(uid) + " -v=" + str(host_directories['base_dir_host']) + ":/data -v=" + str(vepdb_dir_host) + ":/usr/local/share/vep/data -v=" + str(host_directories['input_cna_dir_host']) + ":/workdir/input_cna -v=" + str(host_directories['input_vcf_dir_host']) + ":/workdir/input_vcf -v=" + str(host_directories['input_conf_dir_host']) + ":/workdir/input_conf -v=" + str(host_directories['output_dir_host']) + ":/workdir/output -w=/workdir/output " + str(docker_image_version) + " sh -c \""
   docker_command_run2 = "docker run --rm -t -u " + str(uid) + " -v=" + str(host_directories['base_dir_host']) + ":/data -v=" + str(host_directories['output_dir_host']) + ":/workdir/output -w=/workdir " + str(docker_image_version) + " sh -c \""
   
   
   ## verify VCF and CNA segment file
   logger = getlogger('pcgr-validate-input')
   logger.info("STEP 0: Validate input data")
   vcf_validate_command = str(docker_command_run1) + "pcgr_validate_input.py /data " + str(input_vcf_docker) + " " + str(input_cna_docker) + " " + str(input_conf_docker) + " " + str(genome_assembly) + "\""
   check_subprocess(vcf_validate_command)
   ## Log tumor type of query genome
   logger.info('Finished')
   
   if not input_vcf_docker == 'None':
      
      ## Define input, output and temporary file names
      output_vcf = '/workdir/output/' + str(sample_id) + '.' + str(config_options['tier_model']['tier_model']) + '.vcf.gz'
      #output_maf = '/workdir/output/' + str(sample_id) + '.pcgr.maf'
      output_pass_vcf = '/workdir/output/' + str(sample_id) + '.' + str(config_options['tier_model']['tier_model']) + '.pass.vcf.gz'
      output_pass_tsv = '/workdir/output/' + str(sample_id) + '.' + str(config_options['tier_model']['tier_model']) + '.pass.tsv'
      input_vcf_pcgr_ready = '/workdir/output/' + re.sub(r'(\.vcf$|\.vcf\.gz$)','.pcgr_ready.vcf.gz',host_directories['input_vcf_basename_host'])
      input_vcf_pcgr_ready_uncompressed = '/workdir/output/' + re.sub(r'(\.vcf$|\.vcf\.gz$)','.pcgr_ready.vcf',host_directories['input_vcf_basename_host'])
      vep_vcf = re.sub(r'(\.vcf$|\.vcf\.gz$)','.pcgr_vep.vcf',input_vcf_pcgr_ready)
      vep_vcfanno_vcf = re.sub(r'(\.vcf$|\.vcf\.gz$)','.pcgr_vep.vcfanno.vcf',input_vcf_pcgr_ready)
      vep_tmp_vcf = vep_vcf + '.tmp'
      vep_vcfanno_annotated_vcf = re.sub(r'\.vcfanno','.vcfanno.annotated',vep_vcfanno_vcf) + '.gz'
      vep_vcfanno_annotated_pass_vcf = re.sub(r'\.vcfanno','.vcfanno.annotated.pass',vep_vcfanno_vcf) + '.gz'

      fasta_assembly = "/usr/local/share/vep/data/homo_sapiens/92_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz"
      vep_assembly = 'GRCh37'
      if genome_assembly == 'grch38':
         vep_assembly = 'GRCh38'
         fasta_assembly = "/usr/local/share/vep/data/homo_sapiens/92_GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
      vep_options = "--vcf --check_ref --flag_pick_allele --force_overwrite --species homo_sapiens --assembly " + str(vep_assembly) + " --offline --fork " + str(config_options['other']['n_vep_forks']) + " --hgvs --dont_skip --failed 1 --af --af_1kg --af_gnomad --variant_class --regulatory --domains --symbol --protein --ccds --uniprot --appris --biotype --canonical --gencode_basic --cache --numbers --total_length --allele_number --no_escape --xref_refseq --dir /usr/local/share/vep/data"
      if config_options['other']['vep_skip_intergenic'] == 1:
         vep_options = vep_options + " --no_intergenic"
      vep_main_command = str(docker_command_run1) + "vep --input_file " + str(input_vcf_pcgr_ready) + " --output_file " + str(vep_tmp_vcf) + " " + str(vep_options) + " --fasta " + str(fasta_assembly) + "\""
      vep_sed_command =  str(docker_command_run1) + "sed -r 's/:p\.[A-Z]{1}[a-z]{2}[0-9]+=//g' " + str(vep_tmp_vcf) + " > " + str(vep_vcf) + "\""
      vep_bgzip_command = str(docker_command_run1) + "bgzip -f " + str(vep_vcf) + "\""
      vep_tabix_command = str(docker_command_run1) + "tabix -f -p vcf " + str(vep_vcf) + ".gz" + "\""
      logger = getlogger('pcgr-vep')

      #logger = getlogger('pcgr-vcf2maf')
      #vcf2maf_command = str(docker_command_run1) + "perl /usr/local/bin/vcf2maf.pl --input-vcf " + str(input_vcf_pcgr_ready_uncompressed) + " --output-maf " + str(output_maf) + " --vep-path /home/vep/src/ensembl-vep --vep-data /usr/local/share/vep/data --ref-fasta " + str(fasta_assembly) + " --ncbi_build " + str(ncbi_build_maf) + "\""
      #print(str(vcf2maf_command))
      #check_subprocess(vcf2maf_command)
   
      print()
      logger.info("STEP 1: Basic variant annotation with Variant Effect Predictor (" + str(vep_version) + ", GENCODE " + str(gencode_version) + ", " + str(genome_assembly) + ")")
      check_subprocess(vep_main_command)
      check_subprocess(vep_sed_command)
      check_subprocess(vep_bgzip_command)
      check_subprocess(vep_tabix_command)
      logger.info("Finished")
      #eturn
   
      ## vcfanno command
      print()
      logger = getlogger('pcgr-vcfanno')
      logger.info("STEP 2: Annotation for precision oncology with pcgr-vcfanno (ClinVar, dbNSFP, UniProtKB, cancerhotspots.org, CiVIC, CBMDB, DoCM, TCGA, ICGC-PCAWG, IntoGen_drivers)")
      pcgr_vcfanno_command = str(docker_command_run2) + "pcgr_vcfanno.py --num_processes "  + str(config_options['other']['n_vcfanno_proc']) + " --dbnsfp --docm --clinvar --icgc --civic --cbmdb --intogen_driver_mut --tcga --uniprot --cancer_hotspots --pcgr_onco_xref " + str(vep_vcf) + ".gz " + str(vep_vcfanno_vcf) + " /data/data/" + str(genome_assembly) + "\""
      check_subprocess(pcgr_vcfanno_command)
      logger.info("Finished")
   
      ## summarise command
      print()
      logger = getlogger("pcgr-summarise")
      pcgr_summarise_command = str(docker_command_run2) + "pcgr_summarise.py " + str(vep_vcfanno_vcf) + ".gz /data/data/" + str(genome_assembly) + "\""
      logger.info("STEP 3: Cancer gene annotations with pcgr-summarise")
      check_subprocess(pcgr_summarise_command)
      
      create_output_vcf_command1 = str(docker_command_run2) + 'mv ' + str(vep_vcfanno_annotated_vcf) + ' ' + str(output_vcf) + "\""
      create_output_vcf_command2 = str(docker_command_run2) + 'mv ' + str(vep_vcfanno_annotated_vcf) + '.tbi ' + str(output_vcf) + '.tbi' + "\""
      create_output_vcf_command3 = str(docker_command_run2) + 'mv ' + str(vep_vcfanno_annotated_pass_vcf) + ' ' + str(output_pass_vcf) + "\""
      create_output_vcf_command4 = str(docker_command_run2) + 'mv ' + str(vep_vcfanno_annotated_pass_vcf) + '.tbi ' + str(output_pass_vcf) + '.tbi' + "\""
      clean_command = str(docker_command_run2) + 'rm -f ' + str(vep_vcf) + '* ' + str(vep_vcfanno_annotated_vcf) + ' ' + str(vep_vcfanno_annotated_pass_vcf) + '* ' + str(vep_vcfanno_vcf) + '* ' +  str(input_vcf_pcgr_ready_uncompressed) + "* "  + "\""
      check_subprocess(create_output_vcf_command1)
      check_subprocess(create_output_vcf_command2)
      check_subprocess(create_output_vcf_command3)
      check_subprocess(create_output_vcf_command4)
      pcgr_vcf2tsv_command = str(docker_command_run2) + "vcf2tsv.py " + str(output_pass_vcf) + " --compress " + str(output_pass_tsv) + "\""
      logger.info("Converting VCF to TSV with https://github.com/sigven/vcf2tsv")
      check_subprocess(pcgr_vcf2tsv_command)
      check_subprocess(clean_command)
      logger.info("Finished")

      #return
  
   print()
   
   ## Generation of HTML reports for VEP/vcfanno-annotated VCF and copy number segment file
   if not basic: 
      logger = getlogger('pcgr-writer')
      logger.info("STEP 4: Generation of output files - variant interpretation report for precision oncology")
      pcgr_report_command = str(docker_command_run1) + "/pcgr.R /workdir/output " + str(output_pass_tsv) + ".gz " + str(input_cna_docker) + " "  + str(sample_id)  + " " + str(input_conf_docker) + " " + str(version) + " " + str(genome_assembly) + "\""
      #print(pcgr_report_command)
      check_subprocess(pcgr_report_command)
      logger.info("Finished")

   print()
   
   
   
if __name__=="__main__": __main__()

