#!/usr/bin/env python

import argparse
import pcgrutils
import cyvcf
import os

logger = pcgrutils.getlogger('pcgr-flatten-multisample')

def __main__():
   
   parser = argparse.ArgumentParser(description='Flatten multisample VCF file to single sample VCF files')
   parser.add_argument('vcf_file', help='Multi-sample, bgzipped VCF file with annotated query variants (SNVs/InDels)')
   parser.add_argument('output_postfix', help='Postfix of output file (e.g. \'mutect.pass.annotated.vcf\')')
   args = parser.parse_args()

   flatten_vcf(args.vcf_file,args.output_postfix)

def flatten_vcf(query_vcf, output_postfix):

   logger.info('Read query VCF using cyvcf - ' + str(query_vcf))
   vcf_reader = cyvcf.Reader(open(query_vcf, 'r'))

   #info_tags = info_tags_from_format.split(',')

   k = 0
   for sample_id in vcf_reader.samples:

      file_prefix = str(sample_id) + '.' + str(output_postfix)
      if os.path.exists(file_prefix):
         continue

      file_prefix_tmp = file_prefix + '.tmp'

      logger.info('Generating sample-specific VCF for ' + str(sample_id))

      sample_index_column = 10 + k
      command0 = 'bgzip -dc ' + str(query_vcf) + ' | egrep \'^##\'  > ' + str(file_prefix_tmp)
      os.system(command0)
      command1 = 'bgzip -dc ' + str(query_vcf) + ' | egrep \'^#CHROM\' | cut -f1-9,' + str(sample_index_column) + ' >> ' + str(file_prefix_tmp) 
      os.system(command1)
      command2 = 'bgzip -dc ' + str(query_vcf) + ' | egrep -v \'^#\' | cut -f1-9,' + str(sample_index_column) + ' | egrep -v \'\.\/\.:\' >> ' + str(file_prefix_tmp)
      os.system(command2)
      k = k + 1
      
      info_tags = ['DPT','DPC','ADT','ADC']
      
      printed_tvaf_header = 0
      printed_cvaf_header = 0
      sorted_headers = 0
      
      if os.path.exists(file_prefix_tmp):
         vcf_lines = []
         f = open(file_prefix_tmp,'r')
         for line in f:
            printed = 0
            if line.startswith('##FORMAT=<ID='):
               for tag in info_tags:
                  if line.startswith('##FORMAT=<ID=' + str(tag)):
                     #vcf_lines.append(line.rstrip())
                     vcf_lines.append('##INFO=<ID=' + str(tag) + ',' + ','.join(line.rstrip().split(',')[1:]))
                     if printed_tvaf_header == 0:
                        vcf_lines.append('##INFO=<ID=TVAF,Number=1,Type=Float,Description="Variant allelic fraction (VAF, for alternate allele) in tumor sample">')
                        printed_tvaf_header = 1
                     if printed_cvaf_header == 0:
                        vcf_lines.append('##INFO=<ID=CVAF,Number=1,Type=Float,Description="Variant allelic fraction (VAF, for alternate allele) in control sample">')
                        #vcf_lines.append('##INFO=<ID=CVAF,Number=1,Type=Integer,Description="Variant allelic fraction (VAF, for alternate allele) in control sample">')
                        printed_cvaf_header = 1
                     printed = 1
            if line.startswith('#'):
               if printed == 0 and not line.startswith('##FORMAT'):
                  if line.startswith('#CHROM'):
                     vcf_lines.append('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO')
                  else:
                     vcf_lines.append(line.rstrip())
            else:

               var_info = line.strip().split('\t')
               info_data = var_info[7].split(';')
               format_tags = var_info[8].split(':')
               format_data = var_info[9].split(':')
               
               i = 0
               #print str(format_tags)
               read_depth_tumor = 0
               read_depth_normal = 0
               for format_tag in format_tags:
                  if format_tag in info_tags:
                     if format_tag == 'ADT':
                        read_support_tumor = format_data[i].split(',')
                        if read_support_tumor[0] != '.' and read_support_tumor[1] != '.':
                           tvaf = 0
                           if read_depth_tumor > 0:
                              tvaf = round(float(float(read_support_tumor[1]) / read_depth_tumor),3)
                           info_data.append('TVAF=' + str(tvaf))

                     if format_tag == 'ADC':
                        read_support_normal = format_data[i].split(',')
                        if read_support_normal[0] != '.' and read_support_normal[1] != '.':
                           cvaf = 0
                           if read_depth_normal > 0:
                              cvaf = round(float(float(read_support_normal[1]) / read_depth_normal),3)
                           info_data.append('CVAF=' + str(cvaf))
                           #info_data.append('CVAF=.')
                     
                     if format_tag == 'DPC':
                        if format_data[i] != '.':
                           read_depth_normal = int(format_data[i])
                     if format_tag == 'DPT':
                        if format_data[i] != '.':
                           read_depth_tumor = int(format_data[i])
                           
                     
                     info_data.append(format_tag + '=' + str(format_data[i]))
                  i = i + 1
               
               var_entry = '\t'.join(var_info[0:6]) + '\tPASS\t' + ';'.join(info_data)
               vcf_lines.append(var_entry)
         f.close()
         
         os.system('rm -f ' + str(file_prefix_tmp))
         f = open(file_prefix,'w')
         f.write('\n'.join(vcf_lines))
         f.write('\n')
         f.close()
                     
               
               
         


if __name__=="__main__": __main__()
