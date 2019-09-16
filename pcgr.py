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

pcgr_version = '0.8.2'
db_version = 'PCGR_DB_VERSION = 20190905'
vep_version = '97'
global vep_assembly

def __main__():
   
   parser = argparse.ArgumentParser(description='Personal Cancer Genome Reporter (PCGR) workflow for clinical interpretation of somatic nucleotide variants and copy number aberration segments',formatter_class=argparse.ArgumentDefaultsHelpFormatter, usage="%(prog)s [options] <PCGR_DIR> <OUTPUT_DIR> <GENOME_ASSEMBLY> <CONFIG_FILE> <SAMPLE_ID>" )
   parser.add_argument('--input_vcf', dest = "input_vcf", help='VCF input file with somatic query variants (SNVs/InDels).')
   parser.add_argument('--input_cna', dest = "input_cna",help='Somatic copy number alteration segments (tab-separated values)')
   parser.add_argument('--input_cna_plot', dest = "input_cna_plot",help='Somatic copy number alteration plot')
   parser.add_argument('--pon_vcf',dest = "pon_vcf", help="VCF file with germline calls from Panel of Normals (PON) - i.e. blacklist variants")
   parser.add_argument('--tumor_purity',type = float, dest="tumor_purity", help="Estimated tumor purity (between 0 and 1)")
   parser.add_argument('--tumor_ploidy',type = float, dest="tumor_ploidy", help="Estimated tumor ploidy")
   #parser.add_argument('--tumor_only',action = "store_true",help="Input VCF comes from tumor-only sequencing, calls will be filtered for variants of germline origin")
   #parser.add_argument('--tumor_dp_tag',dest="tumor_dp_tag",help="VCF INFO tag that denotes total sequencing depth at variant site in tumor sample")
   #parser.add_argument('--tumor_af_tag',dest="tumor_af_tag",help="VCF INFO tag that denotes allelic fraction of alternate allele in tumor sample")
   #parser.add_argument('--normal_dp_tag',dest="normal_dp_tag",help="VCF INFO tag that denotes total sequencing depth at variant site in normal (control) sample")
   #parser.add_argument('--normal_af_tag',dest="normal_af_tag",help="VCF INFO tag that denotes allelic fraction of alternate allele in normal (control) sample")
   #parser.add_argument('--tumor_dp_min',dest="tumor_dp_min",default=0,type=int,help="Minimum sequencing depth at variant site in tumor sample (variant filtering criteria, applied to callset before report generation)")
   #parser.add_argument('--tumor_af_min',dest="tumor_af_min",default=0,type=float,help="Minimum allelic fraction of alternate allele in tumor sample (variant filtering criteria, applied to callset before report generation")
   #parser.add_argument('--normal_dp_min',dest="normal_dp_min",default=0,type=int,help="Minimum sequencing depth at variant site in normal (control) sample (variant filtering criteria, applied to callset before report generation)")
   #parser.add_argument('--normal_af_max',dest="normal_af_max",default=1,type=float,help="Maximum value of allelic fraction of alternate allele i normal (control) sample, defaults to 1 (no filtering)")
   #parser.add_argument('--call_conf_tag',dest="call_conf_tag",help="VCF INFO tag that denotes variant call confidence (e.g. categorical variable)")
   #parser.add_argument('--skip_msi',action = "store_true",help="Skip variant-based prediction of MSI status")
   #parser.add_argument('--skip_tmb',action = "store_true",help="Skip estimation of mutational burden")
   #parser.add_argument('--skip_msig',action = "store_true",help="Skip mutational signature analysis")
   #parser.add_argument('--vcf2maf',action = "store_true",help="Generate MAF file of input VCF using vcf2maf")
   #parser.add_argument('--tumor_type',dest="tumor_type",help="Tumor type of query cancer genome, as integer code (make comma-separated if > 1)\nBALLE\nBALLE\n")
   #parser.add_argument('--skip_vcf_validation', action = "store_true",help="Skip validation of input VCF with Ensembl's vcf-validator")
   #parser.add_argument('--target_size_mb',dest="target_size_mb",default=36,type=int,help="Size of coding target region in in megabases (defaults to exome ~ 36Mb")
   parser.add_argument('--force_overwrite', action = "store_true", help='By default, the script will fail with an error if any output file already exists. You can force the overwrite of existing result files by using this flag')
   parser.add_argument('--version', action='version', version='%(prog)s ' + str(pcgr_version))
   parser.add_argument('--basic',action="store_true",help="Run functional variant annotation on VCF through VEP/vcfanno, omit other analyses (i.e. CNA, MSI, report generation etc. (STEP 4)")
   parser.add_argument('--no_vcf_validate', action = "store_true",help="Skip validation of input VCF with Ensembl's vcf-validator")
   parser.add_argument('--docker-uid', dest='docker_user_id', help='Docker user ID. Default is the host system user ID. If you are experiencing permission errors, try setting this up to root (`--docker-uid root`)')
   parser.add_argument('--no-docker', action='store_true', dest='no_docker', default=False, help='Run the PCGR workflow in a non-Docker mode (see install_no_docker/ folder for instructions')
   parser.add_argument('pcgr_dir',help='PCGR base directory with accompanying data directory, e.g. ~/pcgr-0.8.2')
   parser.add_argument('output_dir',help='Output directory')
   parser.add_argument('genome_assembly',choices = ['grch37','grch38'], help='Genome assembly build: grch37 or grch38')
   parser.add_argument('configuration_file',help='PCGR configuration file (TOML format), in conf/ folder')
   parser.add_argument('sample_id',help="Tumor sample/cancer genome identifier - prefix for output files")

   docker_image_version = 'sigven/pcgr:' + str(pcgr_version)
   args = parser.parse_args()
   
   overwrite = 0
   if args.force_overwrite is True:
      overwrite = 1

   logger = getlogger('pcgr-validate-config')

   if args.no_docker:
      docker_image_version = None
   else:
      # check that script and Docker image version correspond
      check_docker_command = 'docker images -q ' + str(docker_image_version)
      output = subprocess.check_output(str(check_docker_command), stderr=subprocess.STDOUT, shell=True)
      if(len(output) == 0):
         err_msg = 'Docker image ' + str(docker_image_version) + ' does not exist, pull image from Dockerhub (docker pull ' + str(docker_image_version) + ')'
         pcgr_error_message(err_msg,logger)

   config_options = {}
   if os.path.exists(args.configuration_file):
      config_options = read_config_options(args.configuration_file, args.pcgr_dir, args.genome_assembly, logger)
   else:
      err_msg = "PCGR configuration file " + str(args.configuration_file) + " does not exist - exiting"
      pcgr_error_message(err_msg,logger)

   tumor_properties = {}
   tumor_properties['tumor_purity'] = "NA"
   tumor_properties['tumor_ploidy'] = "NA"
   if not args.tumor_purity is None:
      if args.tumor_purity > 0 and args.tumor_purity <= 1:
         tumor_properties['tumor_purity'] = str(args.tumor_purity)
      else:
         err_msg = "Tumor purity value " + str(args.tumor_purity) + " is not within a valid range [0,1]"
         pcgr_error_message(err_msg,logger)
   
   if not args.tumor_ploidy is None:
      if args.tumor_ploidy > 0:
         tumor_properties['tumor_ploidy'] = str(args.tumor_ploidy)
      else:
         err_msg = "Tumor ploidy value " + str(args.tumor_ploidy) + " is negative"
         pcgr_error_message(err_msg,logger)
   

   logger = getlogger('pcgr-check-files')
   host_directories = verify_input_files(args.input_vcf, args.input_cna, args.input_cna_plot, args.pon_vcf, args.configuration_file, config_options, args.pcgr_dir, args.output_dir, args.sample_id, args.genome_assembly, overwrite, logger)

   run_pcgr(host_directories, docker_image_version, config_options, args.sample_id, args.genome_assembly, tumor_properties, pcgr_version, args.basic, args.no_vcf_validate, docker_user_id=args.docker_user_id)


def read_config_options(configuration_file, pcgr_dir, genome_assembly, logger):
   
   ## read default options
   pcgr_config_options = {}
   pcgr_config_file_default = os.path.join(pcgr_dir, 'data', str(genome_assembly), 'pcgr_configuration_default.toml')
   if not os.path.exists(pcgr_config_file_default):
      err_msg = "Default PCGR configuration file " + str(pcgr_config_file_default) + " does not exist - exiting"
      pcgr_error_message(err_msg,logger)
   try:
      pcgr_config_options = toml.load(pcgr_config_file_default)
   except (IndexError,TypeError):
      err_msg = 'Configuration file ' + str(configuration_file) + ' is not formatted correctly'
      pcgr_error_message(err_msg, logger)

   ## override with options set by the users
   try:
      user_options = toml.load(configuration_file)
   except (IndexError,TypeError):
      err_msg = 'Configuration file ' + str(configuration_file) + ' is not formatted correctly'
      pcgr_error_message(err_msg, logger)
   
   valid_tumor_types = ['Adrenal_Gland_Cancer_NOS','Ampullary_Carcinoma_NOS','Biliary_Tract_Cancer_NOS','Bladder_Urinary_Tract_Cancer_NOS',
                        'Bone_Cancer_NOS','Breast_Cancer_NOS','Cancer_Unknown_Primary_NOS','Cervical_Cancer_NOS','CNS_Brain_Cancer_NOS',
                        'Colorectal_Cancer_NOS','Esophageal_Cancer_NOS','Head_And_Neck_Cancer_NOS','Kidney_Cancer','Leukemia_NOS',
                        'Liver_Cancer_NOS','Lung_Cancer_NOS','Lymphoma_Hodgkin_NOS','Lymphoma_Non_Hodgkin_NOS','Multiple_Myeloma_NOS',
                        'Ovarian_Fallopian_Tube_Cancer_NOS','Pancreatic_Cancer_NOS','Penile_Cancer_NOS','Peripheral_Nervous_System_Cancer_NOS',
                        'Peritoneal_Cancer_NOS','Pleural_Cancer_NOS','Prostate_Cancer_NOS','Skin_Cancer_NOS','Soft_Tissue_Cancer_Sarcoma_NOS',
                        'Stomach_Cancer_NOS','Testicular_Cancer_NOS','Thymic_Cancer_NOS','Thyroid_Cancer_NOS','Uterine_Cancer_NOS',
                        'Vulvar_Vaginal_Cancer_NOS','']

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
            if section == 'tumor_type' and var == 'type':
               if not str(user_options[section][var]) in valid_tumor_types:
                  err_msg = 'Configuration value for tumor type (' + str(user_options[section][var]) + ') is not a valid type'
                  pcgr_error_message(err_msg, logger)
            #tier_options = ['pcgr','pcgr_acmg']
            normalization_options = ['default','exome','genome','exome2genome']
            theme_options = ['default', 'cerulean', 'journal', 'flatly', 'readable', 'spacelab', 'united', 'cosmo', 'lumen', 'paper', 'sandstone', 'simplex','yeti']
            if var == 'mutsignatures_normalization' and not str(user_options[section][var]) in normalization_options:
               err_msg = 'Configuration value ' + str(user_options[section][var]) + ' for ' + str(var) + \
                  ' cannot be parsed properly (expecting \'default\', \'exome\', \'genome\', or \'exome2genome\')'
               pcgr_error_message(err_msg, logger)
            if var == 'mutsignatures_cutoff' and (float(user_options[section][var]) > 1 or float(user_options[section][var]) < 0) :
               err_msg = 'Configuration value ' + str(user_options[section][var]) + " must be within the [0,1] range"
               pcgr_error_message(err_msg, logger)
            if var == 'mutsignatures_signature_limit':
               if int(user_options[section][var]) < 1 or int(user_options[section][var]) > 30:
                  err_msg = "Number of mutational signatures in search space ('mutsignatures_signature_limit') must be positive and not more than 30 (retrieved value: " + str(user_options[section][var]) + ")"
                  pcgr_error_message(err_msg,logger)
            # if var == 'tier_model' and not str(user_options[section][var]) in tier_options:
            #    err_msg = 'Configuration value ' + str(user_options[section][var]) + ' for ' + str(var) + \
            #       ' cannot be parsed properly (expecting \'pcgr\', or \'pcgr_acmg\')'
            #   pcgr_error_message(err_msg, logger)
            if var == 'report_theme' and not str(user_options[section][var]) in theme_options:
               err_msg = 'Configuration value ' + str(user_options[section][var]) + ' for ' + str(var) + \
                  ' cannot be parsed properly (expecting \'default\', \'cerulean\', \'journal\', \'flatly\', \'readable\', \'spacelab\', \'united\', \'cosmo\', \'lumen\', \'paper\', \'sandstone\', \'simplex\',or \'yeti\')'
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

   #if pcgr_config_options['tumor_type']['type'] == '':
      #err_msg = "Tumor type not defined - please specify a tumor type in the configuration file ([tumor_type] section)"
      #pcgr_error_message(err_msg,logger)
   if pcgr_config_options['msi']['msi'] == 1 and pcgr_config_options['mutational_burden']['mutational_burden'] == 0:
      err_msg = "Prediction of MSI status (msi = true) requires mutational burden/target_size input (mutational_burden = true) - this is currently set as false"
      pcgr_error_message(err_msg,logger)
   if pcgr_config_options['tumor_only']['vcf_tumor_only'] == 1:
      if pcgr_config_options['msi']['msi'] == 1:
         warn_msg = 'Prediction of MSI status in tumor-only mode is not performed - valid for tumor-control data only'
         pcgr_warn_message(warn_msg,logger)
      if pcgr_config_options['mutational_burden']['mutational_burden'] == 1:
         warn_msg = 'Estimation of mutational burden in tumor-only mode is suboptimal - results must be interpreted with caution'
         pcgr_warn_message(warn_msg,logger)
      if pcgr_config_options['mutational_signatures']['mutsignatures'] == 1:
         warn_msg = 'Estimation of mutational signatures in tumor-only mode is not performed - valid for tumor-control data only'
         pcgr_warn_message(warn_msg,logger)
         #err_msg = 'Estimation of mutational signatures in tumor-only mode is suboptimal - results must be interpreted with caution (vcf_tumor_only = true)'
         #pcgr_error_message(warn_msg,logger)


   #print(str(pcgr_config_options))
   return pcgr_config_options


def pcgr_error_message(message, logger):
   logger.error('')
   logger.error(message)
   logger.error('')
   exit(0)

def pcgr_warn_message(message, logger):
   logger.warning(message)

def verify_input_files(input_vcf, input_cna, input_cna_plot, panel_normal_vcf, configuration_file, pcgr_config_options, base_pcgr_dir, output_dir, sample_id, genome_assembly, overwrite, logger):
   """
   Function that checks the input files and directories provided by the user and checks for their existence
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
   
   ## check if panel of normal VCF exist
   if not panel_normal_vcf is None:
      if not os.path.exists(os.path.abspath(panel_normal_vcf)):
         err_msg = "Input file (" + str(panel_normal_vcf) + ") does not exist"
         pcgr_error_message(err_msg,logger)

      if not (os.path.abspath(panel_normal_vcf).endswith('.vcf.gz')):
         err_msg = "Panel of normals VCF file (" + os.path.abspath(panel_normal_vcf) + ") does not have the correct file extension (.vcf.gz)"
         pcgr_error_message(err_msg,logger)

      ## check that tabix file exist if bgzipped files is given
      if os.path.abspath(panel_normal_vcf).endswith('.vcf.gz'):
         tabix_file = panel_normal_vcf + '.tbi'
         if not os.path.exists(os.path.abspath(tabix_file)):
            err_msg = "Tabix file (i.e. '.gz.tbi') is not present for the bgzipped panel of normal VCF file (" + os.path.abspath(input_vcf) + \
               "). Please make sure your the VCF is properly compressed and indexed (bgzip + tabix)"
            pcgr_error_message(err_msg,logger)

      if input_vcf is None:
         warn_msg = "Ignoring panel of normal VCF file, --input_vcf missing"
         pcgr_warn_message(warn_msg, logger)
      else:
         panel_normal_vcf_basename = os.path.basename(str(panel_normal_vcf))
         panel_normal_vcf_dir = os.path.dirname(os.path.abspath(panel_normal_vcf))

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
            err_msg = "Tabix file (i.e. '.gz.tbi') is not present for the bgzipped VCF input file (" + os.path.abspath(input_vcf) + \
               "). Please make sure your input VCF is properly compressed and indexed (bgzip + tabix)"
            pcgr_error_message(err_msg,logger)

      input_vcf_basename = os.path.basename(str(input_vcf))
      input_vcf_dir = os.path.dirname(os.path.abspath(input_vcf))

      ## if output vcf exist and overwrite not set
      output_vcf = os.path.join(str(output_dir_full),str(sample_id)) + '.pcgr_acmg.' + str(genome_assembly)  + '.vcf.gz'
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
   
   ## check if input cna plot file exist
   if not input_cna_plot is None:
      if not os.path.exists(os.path.abspath(input_cna_plot)):
         err_msg = "Input file (" + str(input_cna_plot) + ") does not exist"
         pcgr_error_message(err_msg,logger)
      if not (os.path.abspath(input_cna_plot).endswith('.png')):
         err_msg = "CNA segment input file (" + os.path.abspath(input_cna_plot) + ") does not have the correct file extension (.png)"
         pcgr_error_message(err_msg,logger)
      if input_cna is None:
         err_msg = "Input a CNA plot needs to come with a CNA segment file (--input_cna is missing)"
         pcgr_error_message(err_msg,logger)
      input_cna_plot_basename = os.path.basename(str(input_cna_plot))
      input_cna_plot_dir = os.path.dirname(os.path.abspath(input_cna_plot))

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
      output_cna_segments = os.path.join(str(output_dir_full), str(sample_id)) + '.pcgr_acmg.' + str(genome_assembly) + '.cna_segments.tsv.gz'
      if os.path.exists(output_cna_segments) and overwrite == 0:
         err_msg = "Output files (e.g. " + str(output_cna_segments) + ") already exist - please specify different sample_id or add option --force_overwrite"
         pcgr_error_message(err_msg,logger)
   
   ## check the existence of base folder
   base_dir = os.path.abspath(base_pcgr_dir)
   if not os.path.isdir(base_dir):
      err_msg = "Base directory (" + str(base_dir) + ") does not exist"
      pcgr_error_message(err_msg,logger)
   
   ## check the existence of data folder within the base folder
   db_dir = os.path.join(os.path.abspath(base_pcgr_dir), 'data')
   if not os.path.isdir(db_dir):
      err_msg = "Data directory (" + str(db_dir) + ") does not exist"
      pcgr_error_message(err_msg,logger)
   
   ## check the existence of specified assembly data folder within the base folder
   db_assembly_dir = os.path.join(os.path.abspath(base_pcgr_dir), 'data', genome_assembly)
   if not os.path.isdir(db_assembly_dir):
      err_msg = "Data directory for the specified genome assembly (" + str(db_assembly_dir) + ") does not exist"
      pcgr_error_message(err_msg,logger)
   
   ## check the existence of RELEASE_NOTES (starting from 0.4.0)
   rel_notes_file = os.path.join(os.path.abspath(base_pcgr_dir), 'data', genome_assembly, 'RELEASE_NOTES')
   if not os.path.exists(rel_notes_file):
      err_msg = 'The PCGR data bundle is outdated - please download the latest data bundle (see github.com/sigven/pcgr for instructions)'
      pcgr_error_message(err_msg,logger)
      
   f_rel_not = open(rel_notes_file,'r')
   compliant_data_bundle = 0
   for line in f_rel_not:
      if db_version in line:
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
   

def check_subprocess(command):
   #print(command)
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

def run_pcgr(host_directories, docker_image_version, config_options, sample_id, genome_assembly, tumor_properties, pcgr_version, basic, no_vcf_validate, docker_user_id=None):
   """
   Main function to run the PCGR workflow using Docker
   """

   ## set basic Docker run commands
   output_vcf = 'None'
   output_pass_vcf = 'None'
   output_pass_tsv = 'None'
   output_maf = 'None'
   uid = ''
   gencode_version = 'release 31'
   ncbi_build_maf = 'GRCh38'
   vep_assembly = 'GRCh38'
   if genome_assembly == 'grch37':
      ncbi_build_maf = 'GRCh37'
      gencode_version = 'release 19'
      vep_assembly = 'GRCh37'
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
   vcf_validation = 1
   if no_vcf_validate:
      vcf_validation = 0

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

   if docker_image_version:
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

      ## CNA file only
      if host_directories['input_vcf_dir_host'] == 'NA' and host_directories['input_cna_dir_host'] != 'NA':
         docker_cmd_run1 = docker_cmd_run1  + " -v=" + str(input_cna_volume_mapping)

      ## CNA file and VCF file provided
      if host_directories['input_vcf_dir_host'] != 'NA' and host_directories['input_cna_dir_host'] != 'NA':
         docker_cmd_run1 = docker_cmd_run1  + " -v=" + str(input_vcf_volume_mapping) + " -v=" + str(input_cna_volume_mapping)

      ## CNA plot provided
      if host_directories['input_cna_plot_dir_host'] != "NA":
         docker_cmd_run1 = docker_cmd_run1 + " -v=" + str(input_cna_plot_volume_mapping)
      
      ## Panel of normal VCFs provided
      if host_directories['panel_normal_vcf_dir_host'] != "NA":
         docker_cmd_run1 = docker_cmd_run1 + " -v=" + str(panel_normal_vcf_volume_mapping)

      docker_cmd_run1 = docker_cmd_run1 + " -w=/workdir/output " + str(docker_image_version) + " sh -c \""
      
      docker_cmd_run2 = str(docker_run_basic) + " -v=" + str(databundle_volume_mapping) + " -v=" + str(output_volume_mapping)
      if host_directories['panel_normal_vcf_dir_host'] != "NA":
         docker_cmd_run2 = docker_cmd_run2 + " -v=" + str(panel_normal_vcf_volume_mapping)
      docker_cmd_run2 = docker_cmd_run2 + " -w=/workdir/output " + str(docker_image_version) + " sh -c \""
      docker_cmd_run_end = '\"'

      data_dir = '/data'
      output_dir = '/workdir/output'
      vep_dir = '/usr/local/share/vep/data'
      r_scripts_dir = '/'

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

   check_subprocess(docker_cmd_run1.replace("-u " + str(uid), "") + 'mkdir -p ' + output_dir + docker_cmd_run_end)

   logger = getlogger("pcgr-start")
   logger.info("--- Personal Cancer Genome Reporter workflow ----")
   logger.info("Sample name: " + str(sample_id))
   if config_options['tumor_type']['type'] == "":
      logger.info("Tumor type: Cancer_NOS (Any tumortype)")
   else:
      logger.info("Tumor type: " + str(config_options['tumor_type']['type']))
   logger.info("Tumor-only: " + str(config_options['tumor_only']['vcf_tumor_only']))
   print()

   #logger.info("Finished")

   ## verify VCF and CNA segment file
   logger = getlogger('pcgr-validate-input')
   logger.info("STEP 0: Validate input data")
   vcf_validate_command = docker_cmd_run1 + "pcgr_validate_input.py " + data_dir + " " + str(input_vcf_docker) + " " + str(input_cna_docker) + " " + str(input_conf_docker) + " " + str(panel_normal_docker) + " " + str(vcf_validation) + " " + str(genome_assembly)
   if not docker_image_version:
      vcf_validate_command += ' --output_dir ' + output_dir + docker_cmd_run_end
   else:
      vcf_validate_command += docker_cmd_run_end
   #print(str(vcf_validate_command))
   check_subprocess(vcf_validate_command)
   logger.info('Finished')

   #Valid criteria are: [ canonical appris tsl biotype ccds rank length ]. e.g.:
   #config_options['other']['vep_pick_order'] = canonical,appris,biotype,ccds,rank,tsl,length
   if not input_vcf_docker == 'None':

      ## Define input, output and temporary file names
      output_vcf = os.path.join(output_dir, str(sample_id) + '.' + str(config_options['tier_model']['tier_model']) + '.' + str(genome_assembly) + '.vcf.gz')
      output_pass_vcf = os.path.join(output_dir, str(sample_id) + '.' + str(config_options['tier_model']['tier_model']) + '.' + str(genome_assembly) + '.pass.vcf.gz')
      output_pass_tsv = os.path.join(output_dir, str(sample_id) + '.' + str(config_options['tier_model']['tier_model']) + '.' + str(genome_assembly) + '.pass.tsv')
      output_maf = os.path.join(output_dir, str(sample_id) + '.' + str(genome_assembly) + '.tmp.maf')
      output_vcf2maf_log = os.path.join(output_dir, str(sample_id) + '.' + str(genome_assembly) + '.maf.log')
      input_vcf_pcgr_ready = os.path.join(output_dir, re.sub(r'(\.vcf$|\.vcf\.gz$)', '.pcgr_ready.vcf.gz', host_directories['input_vcf_basename_host']))
      input_vcf_pcgr_ready_uncompressed = os.path.join(output_dir, re.sub(r'(\.vcf$|\.vcf\.gz$)', '.pcgr_ready.vcf', host_directories['input_vcf_basename_host']))
      vep_vcf = re.sub(r'(\.vcf$|\.vcf\.gz$)', '.vep.vcf', input_vcf_pcgr_ready)
      vep_vcfanno_vcf = re.sub(r'(\.vcf$|\.vcf\.gz$)', '.vep.vcfanno.vcf', input_vcf_pcgr_ready)
      vep_vcfanno_annotated_vcf = re.sub(r'\.vcfanno', '.vcfanno.annotated', vep_vcfanno_vcf) + '.gz'
      vep_vcfanno_annotated_pass_vcf = re.sub(r'\.vcfanno', '.vcfanno.annotated.pass', vep_vcfanno_vcf) + '.gz'

      #pick_order = "canonical,appris,tsl,biotype,ccds,rank,length" 
      #pick_order = "biotype,canonical,appris,tsl,ccds,rank,length"
      fasta_assembly = os.path.join(vep_dir, "homo_sapiens", str(vep_version) + "_" + str(vep_assembly), "Homo_sapiens." + str(vep_assembly) + ".dna.primary_assembly.fa.gz")
      #vep_flags = "--vcf --quiet --check_ref --flag_pick_allele_gene --hgvs --dont_skip --af --af_1kg --af_gnomad " + \
      #    "--variant_class --domains --symbol --protein --ccds --uniprot --appris --biotype --canonical --gencode_basic --cache " + \
      #    "--numbers --total_length --no_stats --allele_number --no_escape --xref_refseq"
      vep_flags = "--hgvs --af --af_1kg --af_gnomad --variant_class --domains --symbol --protein --ccds " + \
         "--uniprot --appris --biotype --canonical --gencode_basic --cache --numbers --total_length --allele_number " + \
         "--no_stats --no_escape --xref_refseq --vcf --quiet --check_ref --dont_skip --flag_pick_allele"
      vep_options = "--pick_order " + str(config_options['other']['vep_pick_order']) + " --force_overwrite --species homo_sapiens --assembly " \
         + str(vep_assembly) + " --offline --fork " + str(config_options['other']['n_vep_forks']) + " " + str(vep_flags)  + " --dir " + vep_dir
      vep_options += " --cache_version " + str(vep_version)
      if config_options['other']['vep_skip_intergenic'] == 1:
         vep_options = vep_options + " --no_intergenic"
      vep_main_command = docker_cmd_run1 + "vep --input_file " + str(input_vcf_pcgr_ready) + " --output_file " + str(vep_vcf) + " " + str(vep_options) + " --fasta " + str(fasta_assembly) + docker_cmd_run_end
      vep_bgzip_command = docker_cmd_run1 + "bgzip -f -c " + str(vep_vcf) + " > " + str(vep_vcf) + '.gz' + docker_cmd_run_end
      vep_tabix_command = docker_cmd_run1 + "tabix -f -p vcf " + str(vep_vcf) + ".gz" + docker_cmd_run_end

      #print(str(vep_main_command))
      ## vep commands
      print()
      logger = getlogger('pcgr-vep')
      logger.info("STEP 1: Basic variant annotation with Variant Effect Predictor (" + str(vep_version) + ", GENCODE " + str(gencode_version) + ", " + str(genome_assembly) + ")")
      check_subprocess(vep_main_command)
      check_subprocess(vep_bgzip_command)
      check_subprocess(vep_tabix_command)
      #exit(1)

      if config_options['other']['vcf2maf'] == 1:
         logger.info('Converting VEP-annotated VCF to MAF with https://github.com/mskcc/vcf2maf')
         vcf2maf_command = str(docker_cmd_run1) + "vcf2maf.pl --input-vcf " + str(input_vcf_pcgr_ready_uncompressed) + " --tumor-id " + str(sample_id) + " --output-maf " + str(output_maf) + " --ref-fasta " + str(fasta_assembly) + " --filter-vcf 0 --ncbi-build " + str(ncbi_build_maf) + " > " + str(output_vcf2maf_log) + " 2>&1" + docker_cmd_run_end
         clean_vcf2maf_command = str(docker_cmd_run1) + "rm -f " + str(output_vcf2maf_log) + " " + re.sub(r'(\.vcf$)', '.vep.vcf', input_vcf_pcgr_ready_uncompressed) + " " + docker_cmd_run_end
         check_subprocess(vcf2maf_command)
         check_subprocess(clean_vcf2maf_command)

      logger.info("Finished")

      ## vcfanno command
      print()
      logger = getlogger('pcgr-vcfanno')
      pcgr_vcfanno_command = str(docker_cmd_run2) + "pcgr_vcfanno.py --num_processes " + str(config_options['other']['n_vcfanno_proc']) + " --chasmplus --dbnsfp --docm --clinvar --icgc --civic --cbmdb --tcga_pcdm --winmsk --simplerepeats --tcga --uniprot --cancer_hotspots --pcgr_onco_xref " + str(vep_vcf) + ".gz " + str(vep_vcfanno_vcf) + " " + os.path.join(data_dir, "data", str(genome_assembly))
      if panel_normal_docker != 'None':
         pon_annotation = 1
         pcgr_vcfanno_command = pcgr_vcfanno_command + " --panel_normal_vcf " + str(panel_normal_docker)
         logger.info("STEP 2: Annotation for precision oncology with pcgr-vcfanno (Panel-of-Normals, ClinVar, dbNSFP, UniProtKB, cancerhotspots.org, CiVIC, CBMDB, DoCM, CHASMplus driver mutations, TCGA - putative driver mutations/recurrence, ICGC-PCAWG)")
      else:
         logger.info("STEP 2: Annotation for precision oncology with pcgr-vcfanno (ClinVar, dbNSFP, UniProtKB, cancerhotspots.org, CiVIC, CBMDB, DoCM, CHASMplus driver mutations, TCGA - putative driver mutations/recurrence, ICGC-PCAWG)")
      pcgr_vcfanno_command = pcgr_vcfanno_command + docker_cmd_run_end
      check_subprocess(pcgr_vcfanno_command)
      logger.info("Finished")
      #return

      ## summarise command
      print()
      logger = getlogger("pcgr-summarise")
      pcgr_summarise_command = str(docker_cmd_run2) + "pcgr_summarise.py " + str(vep_vcfanno_vcf) + ".gz " + str(pon_annotation) + " " + str(os.path.join(data_dir, "data", str(genome_assembly))) + docker_cmd_run_end
      logger.info("STEP 3: Cancer gene annotations with pcgr-summarise")
      check_subprocess(pcgr_summarise_command)

      create_output_vcf_command1 = str(docker_cmd_run2) + 'mv ' + str(vep_vcfanno_annotated_vcf) + ' ' + str(output_vcf) + docker_cmd_run_end
      create_output_vcf_command2 = str(docker_cmd_run2) + 'mv ' + str(vep_vcfanno_annotated_vcf) + '.tbi ' + str(output_vcf) + '.tbi' + docker_cmd_run_end
      create_output_vcf_command3 = str(docker_cmd_run2) + 'mv ' + str(vep_vcfanno_annotated_pass_vcf) + ' ' + str(output_pass_vcf) + docker_cmd_run_end
      create_output_vcf_command4 = str(docker_cmd_run2) + 'mv ' + str(vep_vcfanno_annotated_pass_vcf) + '.tbi ' + str(output_pass_vcf) + '.tbi' + docker_cmd_run_end
      clean_command = str(docker_cmd_run2) + 'rm -f ' + str(vep_vcf) + '* ' + str(vep_vcfanno_annotated_vcf) + ' ' + str(vep_vcfanno_annotated_pass_vcf) + '* ' + str(vep_vcfanno_vcf) + '* ' +  str(input_vcf_pcgr_ready_uncompressed) + "* "  + docker_cmd_run_end
      check_subprocess(create_output_vcf_command1)
      check_subprocess(create_output_vcf_command2)
      check_subprocess(create_output_vcf_command3)
      check_subprocess(create_output_vcf_command4)

      ## vcf2tsv command
      pcgr_vcf2tsv_command = str(docker_cmd_run2) + "vcf2tsv.py " + str(output_pass_vcf) + " --compress " + str(output_pass_tsv) + docker_cmd_run_end
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
      pcgr_report_command = (docker_cmd_run1 + os.path.join(r_scripts_dir, "pcgr.R") + " " + output_dir + " " + str(output_pass_tsv) + ".gz" + " " + input_cna_docker + " " + str(sample_id) + " " + input_conf_docker + " " + str(pcgr_version) + " " + genome_assembly + " " + data_dir + " " + str(input_cna_plot_docker) + " " + str(tumor_properties['tumor_purity']) + " " + str(tumor_properties['tumor_ploidy']) + docker_cmd_run_end)
      #print(pcgr_report_command)
      check_subprocess(pcgr_report_command)
      logger.info("Finished")

   print()



if __name__=="__main__": __main__()

