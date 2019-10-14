#!/usr/bin/env python

import csv
import re
import argparse
import os
import subprocess
import logging
import sys
import annoutils
import pandas as np
import toml
from cyvcf2 import VCF

def __main__():
   
   parser = argparse.ArgumentParser(description='Verify input data for PCGR')
   parser.add_argument('pcgr_dir',help='Docker location of PCGR base directory with accompanying data directory, e.g. /data')
   parser.add_argument('input_vcf', help='VCF input file with query variants (SNVs/InDels)')
   parser.add_argument('input_cna', help='Somatic copy number query segments (tab-separated values)')
   parser.add_argument('configuration_file', help='Configuration file (TOML-formatted, e.g. pcgr_conf.toml)')
   parser.add_argument('panel_normal_vcf',help="VCF file with germline calls from panel of normals")
   parser.add_argument('vcf_validation',type=int, default=0,choices=[0,1],help="Perform VCF validation with Ensembl's vcf-validator")
   parser.add_argument('tumor_only',type=int, default=0,choices=[0,1],help="Tumor only sequencing")
   parser.add_argument('genome_assembly',help='grch37 or grch38')
   parser.add_argument('--output_dir', dest='output_dir', help='Output directory', default='/workdir/output')
   args = parser.parse_args()

   ret = validate_pcgr_input(args.pcgr_dir, args.input_vcf, args.input_cna, args.configuration_file, args.panel_normal_vcf, args.vcf_validation, args.tumor_only,args.genome_assembly, args.output_dir)
   if ret != 0:
      sys.exit(1)

def pcgr_error_message(message, logger):
   logger.error('')
   logger.error(message)
   logger.error('')
   return -1

def is_valid_cna(cna_segment_file, logger):
   """
   Function that checks whether the CNA segment file adheres to the correct format
   """
   cna_reader = csv.DictReader(open(cna_segment_file,'r'), delimiter='\t')
   ## check that required columns are present
   if not ('Chromosome' in cna_reader.fieldnames and 'Segment_Mean' in cna_reader.fieldnames and 'Start' in cna_reader.fieldnames and 'End' in cna_reader.fieldnames): 
      err_msg = "Copy number segment file (" + str(cna_segment_file) + ") is missing required column(s): 'Chromosome', 'Start', 'End', or  'Segment_Mean'\n. Column names present in file: " + str(cna_reader.fieldnames)
      return pcgr_error_message(err_msg, logger)
   
   cna_dataframe = np.read_csv(cna_segment_file, sep="\t")
   if not cna_dataframe['Start'].dtype.kind in 'i': ## check that 'Start' is of type integer
      err_msg = '\'Start\' column of copy number segment file contains non-integer values'
      return pcgr_error_message(err_msg, logger)
   if not cna_dataframe['End'].dtype.kind in 'i': ## check that 'End' is of type integer
      err_msg = '\'End\' column of copy number segment file contains non-integer values'
      return pcgr_error_message(err_msg, logger)

   if not cna_dataframe['Segment_Mean'].dtype.kind in 'if': ## check that 'Segment_Mean' is of type integer/float
      err_msg = '\'Segment_Mean\' column of copy number segment file contains non-numerical values'
      return pcgr_error_message(err_msg, logger)

   for rec in cna_reader:
      if int(rec['End']) < int(rec['Start']): ## check that 'End' is always greather than 'Start'
         err_msg = 'Detected wrongly formatted chromosomal segment - \'Start\' is greater than \'End\' (' + str(rec['Chromosome']) + ':' + str(rec['Start']) + '-' + str(rec['End']) + ')'
         return pcgr_error_message(err_msg, logger)
      if int(rec['End']) < 1 or int(rec['Start']) < 1: ## check that 'Start' and 'End' is always non-negative
         err_msg = 'Detected wrongly formatted chromosomal segment - \'Start\' or \'End\' is less than or equal to zero (' + str(rec['Chromosome']) + ':' + str(rec['Start']) + '-' + str(rec['End']) + ')'
         return pcgr_error_message(err_msg, logger)
   logger.info('Copy number segment file (' + str(cna_segment_file) + ') adheres to the correct format')
   return 0

def is_valid_vcf(input_vcf, output_dir, logger):
   """
   Function that reads the output file of EBIvariation/vcf-validator and reports potential errors and validation status 
   """
   
   logger.info('Validating VCF file with EBIvariation/vcf-validator (version 0.6)')
   vcf_validation_output_file = os.path.join(output_dir, re.sub(r'(\.vcf$|\.vcf\.gz$)','.vcf_validator_output',os.path.basename(input_vcf)))
   command_v42 = 'vcf_validator --input ' + str(input_vcf) + ' > ' + str(vcf_validation_output_file) + ' 2>&1'
   if input_vcf.endswith('.gz'):
      command_v42 = 'bgzip -dc ' + str(input_vcf) + ' | vcf_validator  > ' + str(vcf_validation_output_file) + ' 2>&1'
   os.system(command_v42)
   
   #is_valid_vcf = -1
   validation_results = {}
   validation_results['validation_status'] = 0
   validation_results['error_messages'] = []
   if os.path.exists(vcf_validation_output_file):
      f = open(vcf_validation_output_file,'r')
      for line in f:
         if not re.search(r' \(warning\)$|Reading from ',line.rstrip()): ## ignore warnings
            if line.startswith('Line '):
               validation_results['error_messages'].append('ERROR: ' + line.rstrip())
            if 'the input file is valid' in line.rstrip(): ## valid VCF
               validation_results['validation_status'] = 1
            if 'the input file is not valid' in line.rstrip():  ## non-valid VCF
               validation_results['validation_status'] = 0
      f.close()
      os.system('rm -f ' + str(vcf_validation_output_file))
   else:
      err_msg = str(vcf_validation_output_file) + ' does not exist'
      return pcgr_error_message(err_msg, logger)
   
   if validation_results['validation_status'] == 0:
      error_string_42 = '\n'.join(validation_results['error_messages'])
      validation_status = 'According to the VCF specification, the VCF file (' + str(input_vcf) + ') is NOT valid'
      err_msg = validation_status + ':\n' + str(error_string_42)
      return pcgr_error_message(err_msg, logger)
   else:
      validation_status = 'According to the VCF specification, the VCF file ' + str(input_vcf) + ' is valid'
      logger.info(validation_status)
   return 0

def check_existing_vcf_info_tags(input_vcf, pcgr_directory, genome_assembly, logger):
   
   """
   Function that compares the INFO tags in the query VCF and the INFO tags generated by PCGR
   If any coinciding tags, an error will be returned
   """
   
   pcgr_infotags_desc = annoutils.read_infotag_file(os.path.join(pcgr_directory,'data',genome_assembly, 'pcgr_infotags.tsv'))
         
   vcf = VCF(input_vcf)
   logger.info('Checking if existing INFO tags of query VCF file coincide with PCGR INFO tags')
   ret = 1
   for e in vcf.header_iter():
      header_element = e.info()
      if 'ID' in header_element.keys() and 'HeaderType' in header_element.keys():
         if header_element['HeaderType'] == 'INFO':
            if header_element['ID'] in pcgr_infotags_desc.keys():
               err_msg = 'INFO tag ' + str(header_element['ID']) + ' in the query VCF coincides with a VCF annotation tag produced by PCGR - please remove or rename this tag in your query VCF'
               return pcgr_error_message(err_msg, logger)
            if header_element['ID'] == 'DP_TUMOR' or header_element['ID'] == 'AF_TUMOR' or header_element['ID'] == 'AF_NORMAL' or header_element['ID'] == 'DP_NORMAL' or header_element['ID'] == 'CALL_CONFIDENCE':
               err_msg = 'INFO tag ' + str(header_element['ID']) + ' in the query VCF coincides with a VCF annotation tag produced by PCGR - please remove or rename this tag in your query VCF'
               return pcgr_error_message(err_msg, logger)
   
   logger.info('No query VCF INFO tags coincide with PCGR INFO tags')
   return ret

def validate_panel_normal_vcf(vcf, logger):
   
   """
   Function that checks the INFO tags in the panel of normal VCF for the presense of 'PANEL_OF_NORMAL' (logical tag)
   If any coinciding tags, an error will be returned
   """
            
   vcf = VCF(vcf)
   ret = -1
   for e in vcf.header_iter():
      header_element = e.info()
      if 'ID' in header_element.keys() and 'HeaderType' in header_element.keys():
         if header_element['HeaderType'] == 'INFO' and header_element['Type'] == 'Flag':
            if header_element['ID'] == 'PANEL_OF_NORMALS':
               logger.info('Found \'PANEL_OF_NORMALS\' INFO flag in the VCF header section of the of panel of normals VCF file')
               ret = 1


   if ret == -1:
      err_msg = 'INFO flag \'PANEL_OF_NORMALS\' is missing from the panel of normal VCF header'
      return pcgr_error_message(err_msg, logger)

   return ret


def check_format_ad_dp_tags(vcf, pcgr_directory, config_options, tumor_only, logger):
   
   found_taf_tag = 0
   found_tdp_tag = 0
   found_naf_tag = 0
   found_ndp_tag = 0
   found_call_conf_tag = 0
   
   tumor_dp_tag = config_options['allelic_support']['tumor_dp_tag']
   tumor_af_tag = config_options['allelic_support']['tumor_af_tag']
   control_dp_tag = config_options['allelic_support']['control_dp_tag']
   control_af_tag = config_options['allelic_support']['control_af_tag']
   call_conf_tag = config_options['allelic_support']['call_conf_tag']
   
   annoutils.detect_reserved_info_tag(tumor_dp_tag,'tumor_dp_tag', logger)
   annoutils.detect_reserved_info_tag(control_dp_tag,'control_dp_tag', logger)
   annoutils.detect_reserved_info_tag(tumor_af_tag,'tumor_af_tag', logger)
   annoutils.detect_reserved_info_tag(control_af_tag,'control_af_tag', logger)
   annoutils.detect_reserved_info_tag(call_conf_tag,'call_conf_tag', logger)
   
   for e in vcf.header_iter():
      header_element = e.info()
      if 'ID' in header_element.keys() and 'HeaderType' in header_element.keys():
         if header_element['HeaderType'] == 'INFO':
            if header_element['ID'] == tumor_dp_tag:
               if header_element['Type'] == 'Integer':
                  logger.info('Found INFO tag for tumor variant sequencing depth (tumor_dp_tag ' + str(tumor_dp_tag) + ') in input VCF')
                  found_tdp_tag = 1
               else:
                  err_msg = 'INFO tag for tumor variant sequencing depth (tumor_dp_tag ' + str(tumor_dp_tag) + ') is not correctly specified in input VCF (Type=' + str(header_element['Type']) + '), should be Type=Integer'
                  return pcgr_error_message(err_msg, logger)
            if header_element['ID'] == tumor_af_tag:
               if header_element['Type'] == 'Float':
                  logger.info('Found INFO tag for tumor variant allelic fraction (tumor_af_tag ' + str(tumor_af_tag) + ') in input VCF')
                  found_taf_tag = 1
               else:
                  err_msg = 'INFO tag for tumor variant allelic fraction (tumor_af_tag ' + str(tumor_af_tag) + ') is not correctly specified in input VCF (Type=' + str(header_element['Type']) + '), should be Type=Float'
                  return pcgr_error_message(err_msg, logger)
            if header_element['ID'] == control_dp_tag:
               if header_element['Type'] == 'Integer':
                  logger.info('Found INFO tag for normal/control variant sequencing depth (control_dp_tag ' + str(control_dp_tag) + ') in input VCF')
                  found_ndp_tag = 1
               else:
                  err_msg = 'INFO tag for normal/control variant sequencing depth (control_dp_tag ' + str(control_dp_tag) + ') is not correctly specified in input VCF (Type=' + str(header_element['Type']) + '), should be Type=Integer'
                  return pcgr_error_message(err_msg, logger)
            if header_element['ID'] == control_af_tag:
               if header_element['Type'] == 'Float':
                  logger.info('Found INFO tag for normal/control allelic fraction (control_af_tag ' + str(control_af_tag) + ') in input VCF')
                  found_naf_tag = 1
               else:
                  err_msg = 'INFO tag for for normal/control allelic fraction (control_af_tag ' + str(control_af_tag) + ') is not correctly specified in input VCF (Type=' + str(header_element['Type']) + ') should be Type=Float'
                  return pcgr_error_message(err_msg, logger)
            if header_element['ID'] == call_conf_tag:
               if header_element['Type'] == 'String':
                  logger.info('Found INFO tag for variant call confidence (call_conf_tag ' + str(call_conf_tag) + ') in input VCF')
                  found_call_conf_tag = 1
               else:
                  err_msg = 'INFO tag for variant call confidence (call_conf_tag) is not correctly specified in input VCF (Type=' + str(header_element['Type']) + '), should be Type=String'
                  return pcgr_error_message(err_msg, logger)
   
   
   if call_conf_tag != '' and found_call_conf_tag == 0:
      logger.warn('Could not find the specified call_conf_tag (' + str(call_conf_tag) + ') in INFO column of input VCF')
   if tumor_dp_tag != '' and found_tdp_tag == 0:
      logger.warn('Could not find the specified tumor_dp_tag (' + str(tumor_dp_tag) + ') in INFO column of input VCF')
   if tumor_af_tag != '' and found_taf_tag == 0:
      logger.warn('Could not find the specified tumor_af_tag (' + str(tumor_af_tag) + ') in INFO column of input VCF')
   if control_dp_tag != '' and found_ndp_tag == 0:
      logger.warn('Could not find the specified control_dp_tag (' + str(control_dp_tag) + ') in INFO column of input VCF')
   if control_af_tag != '' and found_naf_tag == 0:
      logger.warn('Could not find the specified control_af_tag (' + str(control_af_tag) + ') in INFO column of input VCF')
   
   if config_options['tumor_only']['exclude_likely_hom_germline'] is True and tumor_only == 1 and found_taf_tag == 0:
      logger.warn('Could not find the specified tumor_af_tag (' + str(tumor_af_tag) + ') in INFO column of input VCF - filtering of homozygous germline variants in tumor-only mode will be ignored')

   if config_options['tumor_only']['exclude_likely_het_germline'] is True and tumor_only == 1 and found_taf_tag == 0:
      logger.warn('Could not find the specified tumor_af_tag (' + str(tumor_af_tag) + ') in INFO column of input VCF - filtering of heterozygous germline variants in tumor-only mode will be ignored')


   if found_tdp_tag == 1 and found_taf_tag == 0:
      logger.warn('BOTH \' tumor_dp_tag\' AND \' tumor_af_tag\' need to be specified for use in tumor report (\'tumor_af_tag\' is missing)')
   
   if found_tdp_tag == 0 and found_taf_tag == 1:
      logger.warn('BOTH \'tumor_dp_tag\' AND \'tumor_af_tag\' need to be specified for use in tumor report (\'tumor_dp_tag\' is missing)')
   
   if found_ndp_tag == 1 and found_naf_tag == 0:
      logger.warn('BOTH \'control_dp_tag\' AND \'control_af_tag\' need to be specified for use in tumor report (\'control_af_tag\' is missing)')
   
   if found_ndp_tag == 0 and found_naf_tag == 1:
      logger.warn('BOTH \'control_dp_tag\' AND \'control_af_tag\' need to be specified for use in tumor report (\'control_dp_tag\' is missing)')
   
   ## if filtering turned on for AF-based tumor-only filtering, return error if TVAF not defined

   return 0


def simplify_vcf(input_vcf, vcf, output_dir, logger):
   
   """
   Function that performs tre things on the validated input VCF:
   1. Strip of any genotype data
   2. If VCF have variants with multiple alternative alleles ("multiallelic", e.g. 'A,T'), these are decomposed into variants with a single alternative allele
   3. Final VCF file is sorted and indexed (bgzip + tabix)
   """
   
   input_vcf_pcgr_ready = os.path.join(output_dir, re.sub(r'(\.vcf$|\.vcf\.gz$)','.pcgr_ready.tmp.vcf',os.path.basename(input_vcf)))
   input_vcf_pcgr_ready_decomposed = os.path.join(output_dir, re.sub(r'(\.vcf$|\.vcf\.gz$)','.pcgr_ready.vcf',os.path.basename(input_vcf)))
   
   multiallelic_alt = 0
   for rec in vcf:
      POS = rec.start + 1
      alt = ",".join(str(n) for n in rec.ALT)
      if len(rec.ALT) > 1:
         logger.warning("Multiallelic site detected:" + str(rec.CHROM) + '\t' + str(POS) + '\t' + str(rec.REF) + '\t' + str(alt))
         multiallelic_alt = 1
   command_vcf_sample_free1 = 'egrep \'^##\' ' + str(input_vcf) + ' | egrep -v \'^##FORMAT=\' > ' + str(input_vcf_pcgr_ready)
   command_vcf_sample_free2 = 'egrep \'^#CHROM\' ' + str(input_vcf) + ' | cut -f1-8 >> ' + str(input_vcf_pcgr_ready)
   command_vcf_sample_free3 = 'egrep -v \'^#\' ' + str(input_vcf) + ' | sed \'s/^chr//\' | cut -f1-8 | egrep \'^[0-9]\' | sort -k1,1n -k2,2n -k4,4 -k5,5 >> ' + str(input_vcf_pcgr_ready)
   command_vcf_sample_free4 = 'egrep -v \'^#\' ' + str(input_vcf) + ' | sed \'s/^chr//\' | cut -f1-8 | egrep -v \'^[0-9]\' | egrep \'^[XYM]\' | sort -k1,1 -k2,2n -k4,4 -k5,5 >> ' + str(input_vcf_pcgr_ready)
   command_vcf_sample_free5 = 'egrep -v \'^#\' ' + str(input_vcf) + ' | sed \'s/^chr//\' | cut -f1-8 | egrep -v \'^[0-9]\' | egrep -v \'^[XYM]\' | sort -k1,1 -k2,2n -k4,4 -k5,5 >> ' + str(input_vcf_pcgr_ready)
   if input_vcf.endswith('.gz'):
      command_vcf_sample_free1 = 'bgzip -dc ' + str(input_vcf) + ' | egrep \'^##\' | egrep -v \'^##FORMAT=\' > ' + str(input_vcf_pcgr_ready)
      command_vcf_sample_free2 = 'bgzip -dc ' + str(input_vcf) + ' | egrep \'^#CHROM\' | cut -f1-8 >> ' + str(input_vcf_pcgr_ready)
      command_vcf_sample_free3 = 'bgzip -dc ' + str(input_vcf) + ' | egrep -v \'^#\' | sed \'s/^chr//\' | cut -f1-8 | egrep \'^[0-9]\' | sort -k1,1n -k2,2n -k4,4 -k5,5 >> ' + str(input_vcf_pcgr_ready)
      command_vcf_sample_free4 = 'bgzip -dc ' + str(input_vcf) + ' | egrep -v \'^#\' | sed \'s/^chr//\' | cut -f1-8 | egrep -v \'^[0-9]\' | egrep \'^[XYM]\' | sort -k1,1 -k2,2n -k4,4 -k5,5 >> ' + str(input_vcf_pcgr_ready)
      command_vcf_sample_free5 = 'bgzip -dc ' + str(input_vcf) + ' | egrep -v \'^#\' | sed \'s/^chr//\' | cut -f1-8 | egrep -v \'^[0-9]\' | egrep -v \'^[XYM]\' | sort -k1,1 -k2,2n -k4,4 -k5,5 >> ' + str(input_vcf_pcgr_ready)
   
   os.system(command_vcf_sample_free1)
   os.system(command_vcf_sample_free2)
   os.system(command_vcf_sample_free3)
   os.system(command_vcf_sample_free4)
   os.system(command_vcf_sample_free5)

   if multiallelic_alt == 1:
      logger.info('Decomposing multi-allelic sites in input VCF file using \'vt decompose\'')
      command_decompose = 'vt decompose -s ' + str(input_vcf_pcgr_ready) + ' > ' + str(input_vcf_pcgr_ready_decomposed) + ' 2> ' + os.path.join(output_dir, 'decompose.log')
      os.system(command_decompose)
   else:
      command_copy = 'cp ' + str(input_vcf_pcgr_ready) + ' ' + str(input_vcf_pcgr_ready_decomposed)
      os.system(command_copy)
   
   

   os.system('bgzip -cf ' + str(input_vcf_pcgr_ready_decomposed) + ' > ' + str(input_vcf_pcgr_ready_decomposed) + '.gz')
   os.system('tabix -p vcf ' + str(input_vcf_pcgr_ready_decomposed) + '.gz')

   if os.path.exists(input_vcf_pcgr_ready_decomposed + '.gz') and os.path.getsize(input_vcf_pcgr_ready_decomposed + '.gz') > 0:
      vcf = VCF(input_vcf_pcgr_ready_decomposed + '.gz')
      i = 0
      for rec in vcf:
         i = i + 1
      if len(vcf.seqnames) == 0 or i == 0:
         logger.info('')
         logger.info("Input VCF contains NO valid variants after VCF cleaning - quitting workflow")
         logger.info('')
         exit(1)

   os.system('rm -f ' + str(input_vcf_pcgr_ready) + ' ' + os.path.join(output_dir, 'decompose.log'))

def validate_pcgr_input(pcgr_directory, input_vcf, input_cna, configuration_file, panel_normal_vcf, vcf_validation, tumor_only, genome_assembly, output_dir):
   """
   Function that reads the input files to PCGR (VCF file and Tab-separated values file with copy number segments) and performs the following checks:
   1. Check that VCF file is properly formatted (according to EBIvariation/vcf-validator - VCF v4.2) - optional (vcf_validation in config file)
   2. Check that no INFO annotation tags in the query VCF coincides with those generated by PCGR
   3. Check that provided columns for tumor/normal coverage and allelic depths are found in VCF
   4. Check that if VCF have variants with multiple alternative alleles (e.g. 'A,T') run vt decompose
   5. Check that panel-of-normals VCF adheres to the required format (PANEL_OF_NORMALS INFO tag in header)
   6. Check that copy number segment file has required columns and correct data types (and range)
   7. Any genotype data from VCF input file is stripped, and the resulting VCF file is sorted and indexed (bgzip + tabix) 
   """
   logger = annoutils.getlogger('pcgr-validate-input')
   config_options = annoutils.read_config_options(configuration_file, pcgr_directory, genome_assembly, logger, wflow = 'pcgr')

   if panel_normal_vcf == "None" and tumor_only == 1 and config_options['tumor_only']['exclude_pon'] is True:
      logger.warn('Panel-of-normals VCF is not present - exclusion of calls found in panel-of-normals will be ignored')

   if not input_vcf == 'None':
      if vcf_validation == 1:
         valid_vcf = is_valid_vcf(input_vcf, output_dir, logger)
         if valid_vcf == -1:
            return -1
      else:
         logger.info('Skipping validation of VCF file - as provided by option --no_vcf_validate')
      tag_check = check_existing_vcf_info_tags(input_vcf, pcgr_directory, genome_assembly, logger)
      if tag_check == -1:
         return -1
      
      vcf = VCF(input_vcf)
      allelic_support_check = check_format_ad_dp_tags(vcf, pcgr_directory, config_options, tumor_only, logger)
      if allelic_support_check == -1:
         return -1
      
      simplify_vcf(input_vcf, vcf, output_dir, logger)
   
   if not panel_normal_vcf == "None":
      valid_panel_normals = validate_panel_normal_vcf(panel_normal_vcf, logger)
      if valid_panel_normals == -1:
         return -1
      
   if not input_cna == 'None':
      valid_cna = is_valid_cna(input_cna, logger)
      if valid_cna == -1:
         return -1
   
   return 0
   
if __name__=="__main__": __main__()

