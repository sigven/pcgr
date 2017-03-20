#!/usr/bin/env python

import argparse
import vcf
import re
#import cyvcf
from cyvcf2 import VCF, Writer
import numpy as np



def __main__():
   parser = argparse.ArgumentParser(description='Convert a VCF file with genomic variants to a file with tab-separated values (TSV). Sample genotypes (if present) will be printed on separate lines. Note: Only non-rejected calls/genotypes will be printed in output', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
   parser.add_argument('query_vcf', help='Bgzipped input VCF file with query variants (SNVs/InDels)')
   parser.add_argument('out_tsv', help='Output TSV file with one line pr non-rejected sample genotype (Variant, genotype and annotation data as tab-separated values)')
   parser.add_argument("--skip_info_column",action = "store_true", help="Ignore tags and values in INFO column")
   args = parser.parse_args()
   
   vcf2tsv(args.query_vcf, args.out_tsv, args.skip_info_column)

def vcf2tsv(query_vcf, out_tsv, skip_info_column):
   
   vcf = VCF(query_vcf)
   out = open(out_tsv,'w')
   
   fixed_columns = ['CHROM','POS','ID','REF','ALT','QUAL','FILTER']
   samples = vcf.samples
   info_columns = []
   format_columns = []
   sample_columns = []
   column_types = {}
   
   if len(samples) > 0:
      sample_columns.append('VCF_SAMPLE_ID')
   
   for e in vcf.header_iter():
      header_element = e.info()
      if 'ID' in header_element.keys() and 'HeaderType' in header_element.keys():
         if header_element['HeaderType'] == 'INFO' or header_element['HeaderType'] == 'FORMAT':
            column_types[header_element['ID']] = header_element['Type']
         if header_element['HeaderType'] == 'INFO':
            #if 'AF_EXAC' in header_element['ID']:
            info_columns.append(header_element['ID'])
         if header_element['HeaderType'] == 'FORMAT':
            if len(sample_columns) > 0:
               format_columns.append(header_element['ID'])
   header_line = '\t'.join(fixed_columns) + '\t' + '\t'.join(sorted(info_columns))
   if len(sample_columns) > 0:
      header_line = '\t'.join(fixed_columns) + '\t' + '\t'.join(sorted(info_columns)) + '\t' + '\t'.join(sample_columns) + '\t' + '\t'.join(format_columns)
   
   #out.write(str(header_line) + '\n')
   for rec in vcf:
      rec_id = '.'
      rec_qual = '.'
      rec_filter = '.'
      alt = ",".join(str(n) for n in rec.ALT)
      if not rec.ID is None:
         rec_id = str(rec.ID)
      if not rec.QUAL is None:
         rec_qual = str("{0:.2f}".format(rec.QUAL))
      rec_filter = str(rec.FILTER)
      if rec.FILTER is None:
         rec_filter = 'PASS'
      if type(rec.FILTER) is list:
         if len(rec.FILTER) == 0:
            rec_filter = 'PASS'
         elif len(rec.FILTER) == 1:
            rec_filter = str(rec.FILTER[0])
         else:
            rec_filter = str(';'.join(str(n) for n in rec.FILTER))
      pos = int(rec.start) + 1
      fixed_fields_string = str(rec.CHROM) + '\t' + str(pos) + '\t' + str(rec_id) + '\t' + str(rec.REF) + '\t' + str(alt) + '\t' + str(rec_qual) + '\t' + str(rec_filter)
      
      variant_info = rec.INFO
      var_info_column = []
      if skip_info_column is False:
         for info_field in sorted(info_columns):
            if column_types[info_field] == 'Flag':
               if variant_info.get(info_field) is None:
                  var_info_column.append('False')
               else:
                  var_info_column.append('True')
            elif column_types[info_field] == 'Float' or column_types[info_field] == 'Integer' or column_types[info_field] == 'String':
               if type(variant_info.get(info_field)) is list:
                  var_info_column.append(",".join(str(n) for n in variant_info.get(info_field)))
               else:
                  if variant_info.get(info_field) is None:
                     var_info_column.append('NA')
                  else:
                     if column_types[info_field] == 'Float':
                        val = str("{0:.7f}".format(variant_info.get(info_field)))
                        var_info_column.append(val)
                     else:
                        var_info_column.append(str(variant_info.get(info_field)))
      
      var_id = str(rec.CHROM) + '_' + str(pos) + '_' + str(rec.REF) + '_' + str(alt)
      print str(var_id) + '\t' + str(rec.format('AL','Integer')) + " " + str(rec.genotypes)
      #if len(samples) > 0:
         #i = 0
         #genotypes = list(np.array(rec.gt_types))
         #ref_depths = list(np.array(rec.gt_ref_depths))
         #alt_depths = list(np.array(rec.gt_alt_depths))
         #quals = list(np.array(rec.gt_quals))
         
         #for sample in samples:
            #genotype = int(genotypes[i])
            #ref_depth = ref_depths[i]
            #alt_depth = alt_depths[i]
            #dp = ref_depth + alt_depth
            #print type(genotype)
            #gt_string = 'het'
            #if genotype == 0:
               #gt_string = 'hom_ref'
            #if genotype == 2:
               #gt_string = 'unknown'
            #if genotype == 3:
               #gt_string = 'hom_alt'
            
            #
            #if gt_string == 'het' or gt_string == 'hom_alt':
               #print str(sample) + '\t' + str(var_id) + '\t' + str(rec_filter) + '\t' + str(gt_string) + '\t' + str(genotype) + '\t' + str(ref_depth) + '\t' + str(alt_depth)
            #i = i + 1
      
      # tmp = rec.FORMAT
      # genotypes = list(np.array(rec.gt_types))
      # ref_depths = list(np.array(rec.gt_ref_depths))
      # alt_depths = list(np.array(rec.gt_alt_depths))
      # print str(rec.ID)
      # print len(ref_depths)
      # print str(vcf.samples)
      # print str(genotypes)
      # print str(ref_depths)
      # print str(alt_depths)
      #print str(rec.genotypes)
      #print str(rec.gt_types)
      #print str(rec.gt_ref_depths)
      #print str(rec.gt_alt_depths)
      #print str(rec.format(tmp[0]))
      #print str(rec.FORMAT)
      #exit(0)
      # if len(samples) > 0:
      #    print str
      #    for sample in samples:
      #       genotype_values = []
      #       for format_item in format_columns:
      #          genotype_values.append(rec.format(format_item))
            
            
            
      
      #out.write(str(fixed_fields_string) + '\t' + '\t'.join(var_info_column) + '\n')
   out.close()

if __name__=="__main__": __main__()


   
