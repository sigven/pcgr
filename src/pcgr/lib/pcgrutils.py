#!/usr/bin/env python

import os,re,sys
import csv
import logging
import gzip
import toml

csv.field_size_limit(500 * 1024 * 1024)

def read_infotag_file(vcf_info_tags_tsv):
   """
   Function that reads a VCF info tag file that denotes annotation tags produced by PCGR.
   An example of the VCF info tag file is the following:
   
   tag	number	type	description category
   Consequence	.	String	"Impact modifier for the consequence type (picked by VEP's --flag_pick_allele option)."   vep
   
   A dictionary is returned, with the tag as the key, and the full dictionary record as the value
   """
   info_tag_xref = {} ##dictionary returned
   if not os.path.exists(vcf_info_tags_tsv):
      return info_tag_xref
   tsvfile = open(vcf_info_tags_tsv, 'r')
   reader = csv.DictReader(tsvfile, delimiter='\t')
   for row in reader:
      if not row['tag'] in info_tag_xref:
         info_tag_xref[row['tag']] = row
   
   return info_tag_xref


def pcgr_error_message(message, logger):
   logger.error('')
   logger.error(message)
   logger.error('')
   exit(1)

def pcgr_warn_message(message, logger):
   logger.warning(message)

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
                  tumor_types.append(str(var))
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
            if var == 'cna_overlap_pct':
               if user_options['cna'][var] > 100 or user_options['cna'][var] <= 0:
                  err_msg = "Minimum percent overlap between copy number segment and gene transcript (" + str(user_options[section][var]) + ") should be greater than zero and less than 100"
                  pcgr_error_message(err_msg,logger)
            if var == 'logR_homdel':
               if user_options['cna'][var] > 0:
                  err_msg = "Log ratio for homozygous deletions (" + str(user_options[section][var]) + ") should be less than zero"
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
#    if pcgr_config_options['tumor_only']['vcf_tumor_only'] == 1:
#       if pcgr_config_options['msi']['msi'] == 1:
#          warn_msg = 'Prediction of MSI status is not perfomed in tumor-only mode (vcf_tumor_only = true)'
#          pcgr_warn_message(warn_msg,logger)
#       if pcgr_config_options['mutational_burden']['mutational_burden'] == 1:
#          warn_msg = 'Estimation of mutational burden is not performed in tumor-only mode (vcf_tumor_only = true)'
#          pcgr_warn_message(warn_msg,logger)
#       if pcgr_config_options['mutational_signatures']['mutsignatures'] == 1:
#          warn_msg = 'Estimation of mutational signatures is not perfomed in tumor-only mode (vcf_tumor_only = true)'
#          pcgr_warn_message(warn_msg,logger)

      

   return pcgr_config_options

