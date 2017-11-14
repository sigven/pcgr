#!/usr/bin/env python

import os,re,sys
import csv
import logging
import gzip
from bx.intervals.intersection import IntervalTree
import toml

csv.field_size_limit(500 * 1024 * 1024)

def read_infotag_file(vcf_info_tags_tsv):
   """
   Function that reads a VCF info tag file that denotes annotation tags produced by PCGR.
   An example of the VCF info tag file is the following:
   
   tag	number	type	description
   Consequence	.	String	"Impact modifier for the consequence type (picked by VEP's --flag_pick_allele option)."
   
   A dictionary is returned, with the tag as the key, and the full dictionary record as the value
   """
   info_tag_xref = {} ##dictionary returned
   if not os.path.exists(vcf_info_tags_tsv):
      return info_tag_xref
   with open(vcf_info_tags_tsv, 'rb') as tsvfile:
      reader = csv.DictReader(tsvfile, delimiter='\t')
      for rec in reader:
         if not info_tag_xref.has_key(rec['tag']):
            info_tag_xref[rec['tag']] = rec
   
   return info_tag_xref


def pcgr_error_message(message, logger):
   logger.error('')
   logger.error(message)
   logger.error('')
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



def read_config_options(configuration_file, logger):
   pcgr_config_options = {}
   for section in ['tumor_only','allelic_support','mutational_burden','cna','msi','mutational_signatures','other']:
      pcgr_config_options[section] = {}
   for tag in ['maf_onekg_eur','maf_onekg_amr','maf_onekg_afr','maf_onekg_sas','maf_onekg_eas','maf_onekg_global','maf_gnomad_nfe','maf_gnomad_amr','maf_gnomad_afr','maf_gnomad_sas','maf_gnomad_eas','maf_gnomad_global']:
      pcgr_config_options['tumor_only'][tag] = 0.01
   pcgr_config_options['tumor_only']['exclude_dbsnp_nonclinical'] = 1
   pcgr_config_options['tumor_only']['keep_known_tcga'] = 0
   pcgr_config_options['tumor_only']['tcga_recurrence'] = 2
   pcgr_config_options['tumor_only']['vcf_tumor_only'] = 0
   pcgr_config_options['tumor_only']['exclude_noncoding'] = 1
   for tag in ['normal_dp_tag','normal_af_tag','tumor_dp_tag','tumor_af_tag','call_conf_tag']:
      pcgr_config_options['allelic_support'][tag] = "_na"
   pcgr_config_options['mutational_burden']['target_size_mb'] = 40
   pcgr_config_options['cna']['logR_gain'] = 0.8
   pcgr_config_options['cna']['logR_homdel'] = -0.8
   pcgr_config_options['mutational_signatures']['mutsignatures'] = 1
   pcgr_config_options['mutational_signatures']['mutsignatures_signature_limit'] = 6
   pcgr_config_options['mutational_signatures']['mutsignatures_mutation_limit'] = 50
   pcgr_config_options['mutational_signatures']['mutsignatures_normalization'] = "default"
   pcgr_config_options['msi']['msi'] = 1
   pcgr_config_options['other']['n_vcfanno_proc'] = 4
   pcgr_config_options['other']['n_vep_forks'] = 4
   pcgr_config_options['other']['list_noncoding'] = 1

   try:
      toml_options = toml.load(configuration_file)
   except IndexError,TypeError:
      err_msg = 'Configuration file ' + str(configuration_file) + ' is not formatted correctly'
      pcgr_error_message(err_msg, logger)
   
   float_tags = ['logR_homdel','logR_gain','target_size_mb','maf_onekg_eur','maf_onekg_amr','maf_onekg_afr','maf_onekg_sas','maf_onekg_eas','maf_onekg_global','maf_gnomad_nfe','maf_gnomad_amr','maf_gnomad_afr','maf_gnomad_sas','maf_gnomad_eas','maf_gnomad_global']
   boolean_tags = ['mutsignatures','msi','keep_known_tcga','exclude_dbsnp_nonclinical','exclude_noncoding','vcf_tumor_only','list_noncoding']
   integer_tags = ['n_vcfanno_proc','n_vep_forks','mutsignatures_signature_limit','mutsignatures_mutation_limit','tcga_recurrence']
   string_tags = ['normal_dp_tag','normal_af_tag','tumor_dp_tag','tumor_af_tag','call_conf_tag','mutsignatures_normalization']
   for section in ['tumor_only','allelic_support','mutational_burden','cna','msi','mutational_signatures','other']:
      if toml_options.has_key(section):
         for t in float_tags:
            if toml_options[section].has_key(t):
               if not isinstance(toml_options[section][t],float) and not isinstance(toml_options[section][t],int):
                  err_msg = 'Configuration value ' + str(toml_options[section][t]) + ' for ' + str(t) + ' cannot be parsed properly (expecting float)'
                  pcgr_error_message(err_msg, logger)
               pcgr_config_options[section][t] = toml_options[section][t]
         for t in boolean_tags:
            if toml_options[section].has_key(t):
               if not isinstance(toml_options[section][t],bool):
                  err_msg = 'Configuration value ' + str(toml_options[section][t]) + ' for ' + str(t) + ' cannot be parsed properly (expecting true/false)'
                  pcgr_error_message(err_msg, logger)
               pcgr_config_options[section][t] = int(toml_options[section][t])
         for t in integer_tags:
            if toml_options[section].has_key(t):
               if not isinstance(toml_options[section][t],int):
                  err_msg = 'Configuration value ' + str(toml_options[section][t]) + ' for ' + str(t) + ' cannot be parsed properly (expecting integer)'
                  pcgr_error_message(err_msg, logger)
               pcgr_config_options[section][t] = toml_options[section][t]
         
         for t in string_tags:
            if toml_options[section].has_key(t):
               if not isinstance(toml_options[section][t],basestring):
                  err_msg = 'Configuration value "' + str(toml_options[section][t]) + '" for ' + str(t) + ' cannot be parsed properly (expecting string)'
                  pcgr_error_message(err_msg, logger)
               normalization_options = ['default','exome','genome','exome2genome']
               if tag == 'mutsignatures_normalization' and not str(toml_options[section][t]).encode('utf-8') in normalization_options:
                  err_msg = 'Configuration value ' + str(toml_options[section][t]) + ' for ' + str(t) + ' cannot be parsed properly (expecting \'default\', \'exome\', \'genome\', or \'exome2genome\')'
                  pcgr_error_message(err_msg, logger)
               pcgr_config_options[section][t] = str(toml_options[section][t]).encode('utf-8')
   
   ## check that msig_n is greater than zero and less than 30
   if pcgr_config_options['mutational_signatures']['mutsignatures_signature_limit'] < 0 or pcgr_config_options['mutational_signatures']['mutsignatures_signature_limit'] > 30:
      err_msg = "Number of mutational signatures in search space ('mutsignatures_signature_limit') must be positive and not more than 30 (retrieved value: " + pcgr_config_options['mutational_signatures']['mutsignatures_signature_limit'] + ")"
      pcgr_error_message(err_msg,logger)

   for t in float_tags:
      if t.startswith('maf_'):
         if pcgr_config_options['tumor_only'][t] < 0 or pcgr_config_options['tumor_only'][t] > 1:
            err_msg = "MAF value: " + str(t) + " must be within the [0,1] range, current value " + pcgr_config_options['tumor_only'][t] + ")"
            pcgr_error_message(err_msg,logger)
      if t == 'target_size_mb':
         if pcgr_config_options['mutational_burden'][t] < 0 or pcgr_config_options['mutational_burden'][t] > 50:
            err_msg = "Coding target size value (" + pcgr_config_options['mutational_burden'][t] + ") is not positive or larger than the likely maximum protein-coding size of the human genome (~50Mb))"
            pcgr_error_message(err_msg,logger)
      if t == 'logR_homdel':
         if pcgr_config_options['cna'][t] > 0:
            err_msg = "Log ratio for homozygous deletions (" + pcgr_config_options['cna'][t] + ") should be less than zero"
            pcgr_error_message(err_msg,logger)
      if t == 'logR_gain':
         if pcgr_config_options['cna'][t] < 0:
            err_msg = "Log ratio for copy number amplifications (" + pcgr_config_options['cna'][t] + ") should be greater than zero"
            pcgr_error_message(err_msg,logger)
   
   return pcgr_config_options

