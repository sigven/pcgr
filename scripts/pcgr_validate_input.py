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
from cyvcf2 import VCF

def __main__():

   parser = argparse.ArgumentParser(description='Verify input data for PCGR')
   parser.add_argument('pcgr_dir',help='Docker location of PCGR base directory with accompanying data directory, e.g. /data')
   parser.add_argument('input_vcf', help='VCF input file with somatic (tumor) query variants (SNVs/InDels)')
   parser.add_argument('input_cna', help='Somatic (tumor) copy number query segments (tab-separated values)')
   parser.add_argument('input_rna_fusion', help='Somatic (tumor) RNA fusion variants (tab-separated values)')
   parser.add_argument('input_rna_exp', help='Somatic/ (tumor) gene expression estimates (tab-separated values)')
   parser.add_argument('panel_normal_vcf',help="VCF file with germline calls from panel of normals")
   parser.add_argument('vcf_validation',type=int, default=0,choices=[0,1],help="Perform VCF validation with Ensembl's vcf-validator")
   parser.add_argument('tumor_only',type=int, default=0,choices=[0,1],help="Tumor only sequencing")
   parser.add_argument('genome_assembly',help='grch37 or grch38')
   parser.add_argument('preserved_info_tags', help="Comma-separated string of custom VCF INFO tags to be kept in PCGR output")
   parser.add_argument('tumor_dp_tag', help='VCF INFO tag that denotes tumor sequencing depth')
   parser.add_argument('tumor_af_tag', help='VCF INFO tag that denotes tumor variant allelic fraction')
   parser.add_argument('control_dp_tag', help='VCF INFO tag that denotes control sequencing depth')
   parser.add_argument('control_af_tag', help='VCF INFO tag that denotes control variant allelic fraction')
   parser.add_argument('call_conf_tag', help='VCF INFO tag that denotes somatic variant call confidence')
   parser.add_argument('exclude_hom_germline', help='Logical indicating if homozygote germline calls are to be filtered based on allelic fraction')
   parser.add_argument('exclude_het_germline', help='Logical indicating if heterozygote germline calls are to be filtered based on allelic fraction')

   parser.add_argument('--output_dir', dest='output_dir', help='Output directory', default='/workdir/output')
   args = parser.parse_args()

   ret = validate_pcgr_input(args.pcgr_dir,
                              args.input_vcf,
                              args.input_cna,
                              args.input_rna_fusion,
                              args.input_rna_exp,
                              args.tumor_dp_tag,
                              args.tumor_af_tag,
                              args.control_dp_tag,
                              args.control_af_tag,
                              args.call_conf_tag,
                              args.exclude_hom_germline,
                              args.exclude_het_germline,
                              args.panel_normal_vcf,
                              args.preserved_info_tags,
                              args.vcf_validation,
                              args.tumor_only,
                              args.genome_assembly,
                              args.output_dir)
   if ret != 0:
      sys.exit(1)

def pcgr_error_message(message, logger):
   logger.error('')
   logger.error(message)
   logger.error('')
   return -1


def check_subprocess(command):
   try:
      output = subprocess.check_output(str(command), stderr=subprocess.STDOUT, shell=True)
      if len(output) > 0:
         print (str(output.decode()).rstrip())
   except subprocess.CalledProcessError as e:
      print (e.output.decode())
      exit(0)

def check_preserved_vcf_info_tags(input_vcf, preserved_tags, logger):

   """
   Function that compares the INFO tags in the query VCF and preserved INFO tags set by the user as retained in PCGR output TSV
   If any preserved tag is not in query VCF, an error will be returned
   """

   tags = str(preserved_tags).split(',')
   info_elements_query_vcf = []

   vcf = VCF(input_vcf)
   logger.info('Checking if existing INFO tags of query VCF file matches preserved tags set by the user')
   ret = 1
   for e in vcf.header_iter():
      header_element = e.info()
      if 'ID' in header_element.keys() and 'HeaderType' in header_element.keys():
         if header_element['HeaderType'] == 'INFO':
            info_elements_query_vcf.append(header_element['ID'])


   for t in tags:
      if not t in info_elements_query_vcf:
         err_msg = "Preserved INFO tag '" + str(t) + "' not found among INFO tags in query VCF - make sure preserved VCF INFO tags are set correctly"
         return annoutils.error_message(err_msg, logger)
      else:
         logger.info("Preserved INFO tag '" + str(t) + "' detected among INFO tags in query VCF")

   return ret

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
   if cna_dataframe.empty is True:
      err_msg = 'Copy number segment file is empty - contains NO segments'
      return pcgr_error_message(err_msg, logger)
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


def is_valid_rna_fusion(rna_fusion_file, logger):
   """
   Function that checks whether the RNA fusion transcript file adheres to the correct format
   """
   rna_fusion_reader = csv.DictReader(open(rna_fusion_file,'r'), delimiter='\t')
   ## check that required columns are present
   if not ('GeneA' in rna_fusion_reader.fieldnames and 'GeneB' in rna_fusion_reader.fieldnames and 'Confidence' in rna_fusion_reader.fieldnames):
      err_msg = "RNA fusion file (" + str(rna_fusion_file) + ") is missing required column(s): 'Gene1', 'Gene2', or  'Confidence'\n. Column names present in file: " + str(rna_fusion_reader.fieldnames)
      return pcgr_error_message(err_msg, logger)

   rna_fusion_dataframe = np.read_csv(rna_fusion_file, sep="\t")
   if rna_fusion_dataframe.empty is True:
      err_msg = 'RNA fusion file is empty - contains NO fusions'
      return pcgr_error_message(err_msg, logger)
   if not rna_fusion_dataframe['Gene1'].dtype.kind in 'O': ## check that 'Gene1' is of type object
      err_msg = "'Gene1' column of RNA fusion file cannot not be of type '" + str(rna_fusion_dataframe['Gene1'].dtype) + "'"
      return pcgr_error_message(err_msg, logger)
   if not rna_fusion_dataframe['Gene2'].dtype.kind in 'O': ## check that 'Gene2' is of type object
      err_msg = "'Gene2' column of RNA fusion file cannot not be of type '" + str(rna_fusion_dataframe['Gene2'].dtype) + "'"
      return pcgr_error_message(err_msg, logger)
   if not rna_fusion_dataframe['Confidence'].dtype.kind in 'O': ## check that 'Confidence' is of type object
      err_msg = "'Confidence' column of RNA fusion file cannot not be of type '" + str(rna_fusion_dataframe['Confidence'].dtype) + "'"
      return pcgr_error_message(err_msg, logger)

   observed_variants = {}
   for rec in rna_fusion_reader:
      if not (rec['Confidence'] == 'high' or rec['Confidence'] == 'medium' or rec['Confidence'] == 'low'): ## check that 'Confidence' column harbor permitted values
         err_msg = "Confidence column contains non-permitted values - only 'high','medium', or 'low' permitted. Value entered was " + str(rec['Confidence'])
         return pcgr_error_message(err_msg, logger)

      variant_key = str(rec['Gene1']) + "_" + str(rec['Gene2'])
      if variant_key in observed_variants.keys():
         err_msg = "Duplicate entry in RNA fusion variants: " + str(variant_key) + " is found in multiple rows"
         return pcgr_error_message(err_msg, logger)
      observed_variants[variant_key] = 1


   logger.info('RNA fusion file (' + str(rna_fusion_file) + ') adheres to the correct format')
   return 0

def is_valid_rna_expression(rna_exp_file, logger):
   """
   Function that checks whether the RNA expression file adheres to the correct format
   """
   rna_exp_reader = csv.DictReader(open(rna_exp_file,'r'), delimiter='\t')
   ## check that required columns are present
   if not ('Gene' in rna_exp_reader.fieldnames and 'TPM' in rna_exp_reader.fieldnames and 'Log2FC' in rna_exp_reader.fieldnames and 'PAdj' in rna_exp_reader.fieldnames and 'DiffExp' in rna_exp_reader.fieldnames):
      err_msg = "RNA fusion file (" + str(rna_exp_file) + ") is missing required column(s): 'Gene', 'TPM', 'Log2FC','PAdj', or 'DiffExp'\n. Column names present in file: " + str(rna_exp_reader.fieldnames)
      return pcgr_error_message(err_msg, logger)

   rna_exp_dataframe = np.read_csv(rna_exp_file, sep="\t")
   if rna_exp_dataframe.empty is True:
      err_msg = 'RNA gene expression file is empty - contains NO gene expression estimates'
      return pcgr_error_message(err_msg, logger)
   if not rna_exp_dataframe['Gene'].dtype.kind in 'O': ## check that 'Gene' is of type object
      err_msg = "'Gene' column of RNA expression file cannot not be of type '" + str(rna_exp_dataframe['Gene'].dtype) + "'"
      return pcgr_error_message(err_msg, logger)
   if not rna_exp_dataframe['TPM'].dtype.kind in 'if': ## check that 'TPM' is of type object
      err_msg = "'TPM' column of RNA expression file cannot not be of type '" + str(rna_exp_dataframe['TPM'].dtype) + "'"
      return pcgr_error_message(err_msg, logger)
   if not rna_exp_dataframe['Log2FC'].dtype.kind in 'if': ## check that 'LogFC' is of type object
      err_msg = "'Log2FC' column of RNA expression file cannot not be of type '" + str(rna_exp_dataframe['Log2FC'].dtype) + "'"
      return pcgr_error_message(err_msg, logger)
   if not rna_exp_dataframe['PAdj'].dtype.kind in 'if': ## check that 'PAdj' is of type object
      err_msg = "'TPM' column of RNA expression file cannot not be of type '" + str(rna_exp_dataframe['PAdj'].dtype) + "'"
      return pcgr_error_message(err_msg, logger)
   if not rna_exp_dataframe['DiffExp'].dtype.kind in 'O': ## check that 'DiffExp' is of type object
      err_msg = "'DiffExp' column of RNA expression file cannot not be of type '" + str(rna_exp_dataframe['DiffExp'].dtype) + "'"
      return pcgr_error_message(err_msg, logger)

   for rec in rna_exp_reader:
      if not (rec['DiffExp'] == 'over' or rec['DiffExp'] == 'under' or rec['DiffExp'] == 'NS'): ## check that 'DiffExp' column harbors permitted values
         err_msg = "Confidence column contains non-permitted values - only 'over','under', or 'NS' permitted. Value entered was " + str(rec['DiffExp'])
         return pcgr_error_message(err_msg, logger)

      if not (rec['TPM'] >= 0):
         err_msg = "'TPM' column cannot contain negative values - value was " + str(rec['TPM'])
         return pcgr_error_message(err_msg, logger)
      if not (rec['PAdj'] >= 0):
         err_msg = "'PAdj' column (adjusted p-value from differential expression testing) cannot contain negative values - value was " + str(rec['PAdj'])
         return pcgr_error_message(err_msg, logger)


   logger.info('RNA expression file (' + str(rna_exp_file) + ') adheres to the correct format')
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


def check_format_ad_dp_tags(vcf,
                           tumor_dp_tag,
                           tumor_af_tag,
                           control_dp_tag,
                           control_af_tag,
                           call_conf_tag,
                           exclude_hom_germline,
                           exclude_het_germline,
                           tumor_only,
                           logger):

   """
   Function that checks whether the INFO tags specified for depth/allelic fraction are correctly formatted in the VCF header (i.e. Type)
   """

   found_taf_tag = 0
   found_tdp_tag = 0
   found_naf_tag = 0
   found_ndp_tag = 0
   found_call_conf_tag = 0

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


   if call_conf_tag != '_NA_' and found_call_conf_tag == 0:
      logger.warning(f"Could not find the specified call_conf_tag ('{call_conf_tag}') in INFO column of input VCF")
   if tumor_dp_tag != '_NA_' and found_tdp_tag == 0:
      logger.warning(f"Could not find the specified tumor_dp_tag ('{tumor_dp_tag}') in INFO column of input VCF")
   if tumor_af_tag != '_NA_' and found_taf_tag == 0:
      logger.warning(f"Could not find the specified tumor_af_tag ('{tumor_af_tag}') in INFO column of input VCF")
   if control_dp_tag != '_NA_' and found_ndp_tag == 0:
      logger.warning(f"Could not find the specified control_dp_tag ('{control_dp_tag}') in INFO column of input VCF")
   if control_af_tag != '_NA_' and found_naf_tag == 0:
      logger.warning(f"Could not find the specified control_af_tag ('{control_af_tag}') in INFO column of input VCF")

   if exclude_hom_germline is True and tumor_only == 1 and found_taf_tag == 0:
      logger.warning(f"Could not find the specified tumor_af_tag ('{tumor_af_tag}') in INFO column of input VCF - filtering of homozygous germline variants in tumor-only mode will be ignored")

   if exclude_het_germline is True and tumor_only == 1 and found_taf_tag == 0:
      logger.warning(f"Could not find the specified tumor_af_tag ('{tumor_af_tag}') in INFO column of input VCF - filtering of heterozygous germline variants in tumor-only mode will be ignored")


   if found_tdp_tag == 1 and found_taf_tag == 0:
      logger.warning('BOTH \' tumor_dp_tag\' AND \' tumor_af_tag\' need to be specified for use in tumor report (\'tumor_af_tag\' is missing)')

   if found_tdp_tag == 0 and found_taf_tag == 1:
      logger.warning('BOTH \'tumor_dp_tag\' AND \'tumor_af_tag\' need to be specified for use in tumor report (\'tumor_dp_tag\' is missing)')

   if found_ndp_tag == 1 and found_naf_tag == 0:
      logger.warning('BOTH \'control_dp_tag\' AND \'control_af_tag\' need to be specified for use in tumor report (\'control_af_tag\' is missing)')

   if found_ndp_tag == 0 and found_naf_tag == 1:
      logger.warning('BOTH \'control_dp_tag\' AND \'control_af_tag\' need to be specified for use in tumor report (\'control_dp_tag\' is missing)')

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

   check_subprocess(command_vcf_sample_free1)
   check_subprocess(command_vcf_sample_free2)
   check_subprocess(command_vcf_sample_free3)
   check_subprocess(command_vcf_sample_free4)
   check_subprocess(command_vcf_sample_free5)

   if multiallelic_alt == 1:
      logger.info('Decomposing multi-allelic sites in input VCF file using \'vt decompose\'')
      command_decompose = 'vt decompose -s ' + str(input_vcf_pcgr_ready) + ' > ' + str(input_vcf_pcgr_ready_decomposed) + ' 2> ' + os.path.join(output_dir, 'decompose.log')
      check_subprocess(command_decompose)
   else:
      command_copy = 'cp ' + str(input_vcf_pcgr_ready) + ' ' + str(input_vcf_pcgr_ready_decomposed)
      check_subprocess(command_copy)

   check_subprocess('bgzip -cf ' + str(input_vcf_pcgr_ready_decomposed) + ' > ' + str(input_vcf_pcgr_ready_decomposed) + '.gz')
   check_subprocess('tabix -p vcf ' + str(input_vcf_pcgr_ready_decomposed) + '.gz')

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

   check_subprocess('rm -f ' + str(input_vcf_pcgr_ready) + ' ' + os.path.join(output_dir, 'decompose.log'))

def validate_pcgr_input(pcgr_directory,
                        input_vcf,
                        input_cna,
                        input_rna_fusion,
                        input_rna_expression,
                        tumor_dp_tag,
                        tumor_af_tag,
                        control_dp_tag,
                        control_af_tag,
                        call_conf_tag,
                        exclude_hom_germline,
                        exclude_het_germline,
                        panel_normal_vcf,
                        preserved_info_tags,
                        vcf_validation,
                        tumor_only,
                        genome_assembly,
                        output_dir):
   """
   Function that reads the input files to PCGR (VCF file and Tab-separated values file with copy number segments) and performs the following checks:
   1. Check that no INFO annotation tags in the query VCF coincides with those generated by PCGR
   2. Check that provided columns for tumor/normal coverage and allelic depths are found in VCF
   3. Check that provided preserved VCF columns are present in VCF file
   4. Check that if VCF have variants with multiple alternative alleles (e.g. 'A,T') run vt decompose
   5. Check that panel-of-normals VCF adheres to the required format (PANEL_OF_NORMALS INFO tag in header)
   6. Any genotype data from VCF input file is stripped, and the resulting VCF file is sorted and indexed (bgzip + tabix)
   7. Check that copy number segment file has required columns and correct data types (and range)
   8. Check that RNA fusion variant file has required columns and correct data types
   9. Check that RNA expression file has required columns and correct data types
   """
   logger = annoutils.getlogger('pcgr-validate-arguments-input')

   # if panel_normal_vcf == "None" and tumor_only == 1 and config_options['tumor_only']['exclude_pon'] is True:
   #    logger.warning('Panel-of-normals VCF is not present - exclusion of calls found in panel-of-normals will be ignored')

   if not input_vcf == 'None':

      ## Perform VCF validation if this option is set
      if vcf_validation == 1:
         logger.info('Skipping validation of VCF file (deprecated as of Dec 2021)')
      else:
         logger.info('Skipping validation of VCF file as provided by option --no_vcf_validate')

      ## Check that VCF does not contain INFO tags that will be appended with PCGR annotation
      tag_check = check_existing_vcf_info_tags(input_vcf, pcgr_directory, genome_assembly, logger)
      if tag_check == -1:
         return -1

      if preserved_info_tags != "None":
         custom_check = check_preserved_vcf_info_tags(input_vcf, preserved_info_tags, logger)
         if custom_check == -1:
            return -1

      ## Check whether specified tags for depth/allelic fraction is properly defined in VCF
      vcf = VCF(input_vcf)
      allelic_support_check = check_format_ad_dp_tags(vcf, tumor_dp_tag, tumor_af_tag, control_dp_tag,
                                                      control_af_tag, call_conf_tag, exclude_hom_germline,
                                                      exclude_het_germline, tumor_only, logger)
      if allelic_support_check == -1:
         return -1

      ## Simplify VCF - remove multiallelic variants
      simplify_vcf(input_vcf, vcf, output_dir, logger)


   ## Validate panel-of-normals VCF is this is provided
   if not panel_normal_vcf == "None":
      valid_panel_normals = validate_panel_normal_vcf(panel_normal_vcf, logger)
      if valid_panel_normals == -1:
         return -1

   ## Check whether file with copy number aberration segments is properly formatted
   if not input_cna == 'None':
      valid_cna = is_valid_cna(input_cna, logger)
      if valid_cna == -1:
         return -1

   ## Check whether file with RNA fusion variants is properly formatted
   if not input_rna_fusion == 'None':
      valid_rna_fusion = is_valid_rna_fusion(input_rna_fusion, logger)
      if valid_rna_fusion == -1:
         return -1

   ## Check whether file with RNA fusion variants is properly formatted
   if not input_rna_expression == 'None':
      valid_rna_expression = is_valid_rna_expression(input_rna_expression, logger)
      if valid_rna_expression == -1:
         return -1

   return 0

if __name__=="__main__": __main__()

