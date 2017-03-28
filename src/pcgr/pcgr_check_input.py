#!/usr/bin/env python

import csv
import re
import argparse
import os
import subprocess
import logging
import sys
import pandas as np
from cyvcf2 import VCF

def __main__():
   
   parser = argparse.ArgumentParser(description='Verify input data for PCGR')
   parser.add_argument('input_vcf', help='VCF input file with somatic query variants (SNVs/InDels)')
   parser.add_argument('input_cna_segments', help='Somatic copy number query segments (tab-separated values)')   
   args = parser.parse_args()
   
   verify_input(args.input_vcf, args.input_cna_segments)

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


def is_valid_cna_segment_file(cna_segment_file, logger):
   cna_dataframe = np.read_csv(cna_segment_file, sep="\t")
   if not cna_dataframe['Start'].dtype.kind in 'i':
      logger.error('')
      logger.error('\'Start\' column of copy number segment file contains non-integer values')
      logger.error('')
      return -1
   if not cna_dataframe['End'].dtype.kind in 'i':
      logger.error('')
      logger.error('\'End\' column of copy number segment file contains non-integer values')
      logger.error('')
      return -1
   if not cna_dataframe['Segment_Mean'].dtype.kind in 'if':
      logger.error('')
      logger.error('\'Segment_Mean\' column of copy number segment file contains non-numerical values')
      logger.error('')
      return -1
   logger.info('Copy number segment file ' + str(cna_segment_file) + ') adheres to the correct format')
   return 0



def is_valid_vcf(vcf_validator_output_file):
   
   valid_vcf = -1
   f = open(vcf_validator_output_file,'r')
   error_messages = []
   for line in f:
      if not re.search(r' \(warning\)$|^Reading from ',line.rstrip()):
         if line.startswith('Line '):
            error_messages.append(line.rstrip())
         if line.endswith('the input file is valid'):
            valid_vcf = 1
         if line.endswith('the input file is not valid'):
            valid_vcf = 0
   f.close()
   
   ret = {}
   ret['error_messages'] = error_messages
   ret['validation_status'] = valid_vcf
   return ret


def verify_input(input_vcf, input_cna_segments):
   
   logger = getlogger('pcgr-check-input')
   input_vcf_pcgr_ready = re.sub(r'(\.vcf$|\.vcf\.gz$)','.pcgr_ready.vcf',input_vcf)
   if not input_vcf == 'None':
      logger.info('Validating VCF file with EBIvariation/vcf-validator')
      vcf_validation_output_file = re.sub(r'(\.vcf$|\.vcf\.gz$)','.vcf_validator_output',input_vcf)
      command_v42 = 'vcf_validator --input ' + str(input_vcf) + ' --version v4.2 > ' + str(vcf_validation_output_file)
      if input_vcf.endswith('.gz'):
         command_v42 = 'bgzip -dc ' + str(input_vcf) + ' | vcf_validator --version v4.2 > ' + str(vcf_validation_output_file)
      
      os.system(command_v42)
      validation_results = is_valid_vcf(vcf_validation_output_file)
      
      if not validation_results['validation_status']:
         error_string_42 = '\n'.join(validation_results['error_messages'])
         validation_status = 'VCF file is NOT valid according to v4.2 specification'
         logger.error(validation_status + ':\n' + str(error_string_42))
         return -1
      else:
         validation_status = 'VCF file is valid according to v4.2 specification'
         logger.info(validation_status)
   
      if validation_results['validation_status']:
         multiallelic_alt = 0
         vcf = VCF(input_vcf)
         for rec in vcf:
            chrom = rec.CHROM
            if chrom.startswith('chr'):
               error_message_chrom = "'chr' must be stripped from chromosome names: " + str(rec.CHROM)
               logger.error(error_message_chrom)
               return -1
            POS = rec.start + 1
            alt = ",".join(str(n) for n in rec.ALT)
            if len(rec.ALT) > 1:
               error_message_multiallelic = "Multiallelic site detected:" + str(rec.CHROM) + '\t' + str(POS) + '\t' + str(rec.REF) + '\t' + str(alt)
               logger.error(error_message_multiallelic)
               multiallelic_alt = 1
               return -1
         command_vcf_sample_free1 = 'egrep \'^#\' ' + str(input_vcf) + ' > ' + str(input_vcf_pcgr_ready)
         command_vcf_sample_free2 = 'egrep -v \'^#\' ' + str(input_vcf) + ' | cut -f1-8 > ' + str(input_vcf_pcgr_ready)
         if input_vcf.endswith('.gz'):
            command_vcf_sample_free1 = 'bgzip -dc ' + str(input_vcf) + ' | egrep \'^#\' > ' + str(input_vcf_pcgr_ready)
            command_vcf_sample_free2 = 'bgzip -dc ' + str(input_vcf) + ' | egrep -v \'^#\' | cut -f1-8 > ' + str(input_vcf_pcgr_ready)
         os.system(command_vcf_sample_free1)
         os.system(command_vcf_sample_free2)
         os.system('bgzip -f ' + str(input_vcf_pcgr_ready))
         os.system('tabix -p vcf ' + str(input_vcf_pcgr_ready) + '.gz')
      
   if not input_cna_segments == 'None':
      with open(input_cna_segments, 'rb') as tsvfile:
         cna_reader = csv.DictReader(tsvfile, delimiter='\t')
         if not ('Chromosome' in cna_reader.fieldnames and 'Segment_Mean' in cna_reader.fieldnames and 'Start' in cna_reader.fieldnames and 'End' in cna_reader.fieldnames):
            error_message_cnv = "Copy number segment file (" + str(input_cna_segments) + ") is missing required column(s): 'Chromosome', 'Start', 'End', or  'Segment_Mean'"
            logger.error(error_message_cnv)
            return -1
         else:
            ret = is_valid_cna_segment_file(input_cna_segments, logger)
            return ret
   
   return 0
   
if __name__=="__main__": __main__()

