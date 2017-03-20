#!/usr/bin/env python

import argparse
import pcgr
import cyvcf
import os
import utils

logger = pcgr.getlogger('pcgr-flatten-multisample')

def __main__():
   
   parser = argparse.ArgumentParser(description='Flatten multisample VCF file to single sample VCF files')
   parser.add_argument('vcf_file', help='Multi-sample, bgzipped VCF file with annotated query variants (SNVs/InDels)')
   parser.add_argument('output_postfix', help='Postfix of output file (e.g. \'mutect.pass.annotated.vcf\')')
   args = parser.parse_args()

   flatten_vcf(args.vcf_file,args.output_postfix)

def flatten_vcf(query_vcf, output_postfix):

   logger.info('Read query VCF using cyvcf - ' + str(query_vcf))
   vcf_reader = cyvcf.Reader(open(query_vcf, 'r'))

   k = 0
   for sample_id in vcf_reader.samples:

      file_prefix = str(sample_id) + '.' + str(output_postfix)
      if os.path.exists(file_prefix):
         continue

      logger.info('Generating sample-specific VCF for ' + str(sample_id))

      sample_index_column = 10 + k
      command0 = 'bgzip -dc ' + str(query_vcf) + ' | egrep \'^##\'  > ' + str(file_prefix)
      os.system(command0)
      command1 = 'bgzip -dc ' + str(query_vcf) + ' | egrep \'^#CHROM\' | cut -f1-9,' + str(sample_index_column) + ' >> ' + str(file_prefix) 
      os.system(command1)
      command2 = 'bgzip -dc ' + str(query_vcf) + ' | egrep -v \'^#\' | cut -f1-9,' + str(sample_index_column) + ' | egrep -v \'\.\/\.:\' >> ' + str(file_prefix)
      os.system(command2)
      k = k + 1


if __name__=="__main__": __main__()
