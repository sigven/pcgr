#!/usr/bin/env python

import argparse
import vcf
import re
from cyvcf2 import VCF, Writer
import numpy as np

def __main__():
   parser = argparse.ArgumentParser(description='Convert a VCF file with genomic variants to a file with tab-separated values (TSV). One entry (TSV line) per sample genotype', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
   parser.add_argument('query_vcf', help='Bgzipped input VCF file with query variants (SNVs/InDels)')
   parser.add_argument('out_tsv', help='Output TSV file with one line pr non-rejected sample genotype (Variant, genotype and annotation data as tab-separated values)')
   parser.add_argument("--skip_info_data",action = "store_true", help="Skip printing of data in INFO column")
   parser.add_argument("--skip_genotype_data", action="store_true", help="Skip printing of genotype_data (FORMAT columns)")
   parser.add_argument("--keep_rejected_calls", action="store_true", help="Print data for rejected calls")
   args = parser.parse_args()
   
   vcf2tsv(args.query_vcf, args.out_tsv, args.skip_info_data, args.skip_genotype_data, args.keep_rejected_calls)

def strip_null_values(dat):
   if dat == '-2147483648' or dat == -2147483648:
      return '.'
   elements = dat.split(',')
   if len(elements) == 1:
      if elements[0] == '-2147483648':
         return '.'
      else:
         return dat
   if len(elements) > 1:
      stripped_for_null_values = [item for item in elements if (item != '-2147483647' and item != '-2147483648')]
      #stripped_for_null_values = [item for item in stripped_for_null_values if item != '-2147483648']
      if len(stripped_for_null_values) == 0:
         return '.'
      return ','.join(stripped_for_null_values)
         

def vcf2tsv(query_vcf, out_tsv, skip_info_data, skip_genotype_data, keep_rejected_calls):
   
   vcf = VCF(query_vcf, gts012 = True)
   out = open(out_tsv,'w')
   
   fixed_columns_header = ['CHROM','POS','ID','REF','ALT','QUAL','FILTER']
   samples = vcf.samples
   info_columns_header = []
   format_columns_header = []
   sample_columns_header = []
   column_types = {}
   gt_present_header = 0
   
   if len(samples) > 0:
      sample_columns_header.append('VCF_SAMPLE_ID')
   
   for e in vcf.header_iter():
      header_element = e.info()
      if 'ID' in header_element.keys() and 'HeaderType' in header_element.keys():
         if header_element['HeaderType'] == 'INFO' or header_element['HeaderType'] == 'FORMAT':
            column_types[header_element['ID']] = header_element['Type']
         if header_element['HeaderType'] == 'INFO':
            if skip_info_data is False:
               info_columns_header.append(header_element['ID'])
         if header_element['HeaderType'] == 'FORMAT':
            if len(sample_columns_header) > 0 and skip_genotype_data is False:
               if header_element['ID'] != 'GT':
                  format_columns_header.append(header_element['ID'])
               else:
                  gt_present_header = 1

   header_line = '\t'.join(fixed_columns_header)
   if skip_info_data is False:
      header_line = '\t'.join(fixed_columns_header) + '\t' + '\t'.join(sorted(info_columns_header))
      if len(sample_columns_header) > 0:
         if skip_genotype_data is False:
            header_line = '\t'.join(fixed_columns_header) + '\t' + '\t'.join(sorted(info_columns_header)) + '\t' + '\t'.join(sample_columns_header) + '\t' + '\t'.join(sorted(format_columns_header)) + '\tGT'
         else:
            header_line = '\t'.join(fixed_columns_header) + '\t' + '\t'.join(sorted(info_columns_header))
   else:
      if len(sample_columns_header) > 0:
         if skip_genotype_data is False:
            header_line = '\t'.join(fixed_columns_header) + '\t' + '\t'.join(sample_columns_header) + '\t' + '\t'.join(sorted(format_columns_header)) + '\tGT'
         else:
            header_line = '\t'.join(fixed_columns_header)

   out.write(str(header_line) + '\n')
   #print str(header_line)
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
      vcf_info_data = []
      if skip_info_data is False:
         for info_field in sorted(info_columns_header):
            if column_types[info_field] == 'Flag':
               if variant_info.get(info_field) is None:
                  vcf_info_data.append('False')
               else:
                  vcf_info_data.append('True')
            elif column_types[info_field] == 'Float' or column_types[info_field] == 'Integer' or column_types[info_field] == 'String':
               if type(variant_info.get(info_field)) is list:
                  vcf_info_data.append(",".join(str(n) for n in variant_info.get(info_field)))
               else:
                  if variant_info.get(info_field) is None:
                     vcf_info_data.append('.')
                  else:
                     if column_types[info_field] == 'Float':
                        val = str("{0:.7f}".format(variant_info.get(info_field)))
                        vcf_info_data.append(val)
                     else:
                        vcf_info_data.append(str(variant_info.get(info_field)))
      
      #dictionary, with sample names as keys, values being genotype data (dictionary with format tags as keys)
      vcf_sample_genotype_data = {}
      if len(samples) > 0 and skip_genotype_data is False:
         if gt_present_header == 1:
            gt_cyvcf = rec.gt_types
            i = 0
            while i < len(samples):
               vcf_sample_genotype_data[samples[i]] = {}
               gt = './.'
               if gt_cyvcf[i] == 0:
                  gt = '0/0'
               if gt_cyvcf[i] == 1:
                  gt = '0/1'
               if gt_cyvcf[i] == 2:
                  gt = '1/1'
               vcf_sample_genotype_data[samples[i]]['GT'] = gt
               i = i + 1
               
      for format_tag in sorted(format_columns_header):
         if len(samples) > 0 and skip_genotype_data is False:
            sample_dat = rec.format(format_tag)
            dim = sample_dat.shape
            j = 0
            ## sample-wise
            while j < dim[0]:
               if sample_dat[j].size > 1:
                  d = ','.join(str(e) for e in np.ndarray.tolist(sample_dat[j]))
                  dat = strip_null_values(d)
                  #print str(format_tag) + '\t' + str(d) + '\t' + str(dat)
                  if vcf_sample_genotype_data.has_key(samples[j]):
                     vcf_sample_genotype_data[samples[j]][format_tag] = dat
               else:
                  d = str(sample_dat[j][0])
                  dat = strip_null_values(d)
                  #print str(format_tag) + '\t' + str(d) + '\t' + str(dat)
                  if vcf_sample_genotype_data.has_key(samples[j]):
                     vcf_sample_genotype_data[samples[j]][format_tag] = dat
               j = j + 1
      
      tsv_elements = []
      tsv_elements.append(fixed_fields_string)
      if skip_info_data is False:
         if len(sample_columns_header) > 0:
            if skip_genotype_data is False:
               tsv_elements.append('\t'.join(vcf_info_data))
               for s in sorted(vcf_sample_genotype_data.keys()):
                  sample = s
                  line_elements = []
                  line_elements.extend(tsv_elements)
                  line_elements.append(sample)
                  gt_tag = '.'
                  for tag in sorted(vcf_sample_genotype_data[sample].keys()):
                     if tag != 'GT':
                        line_elements.append(vcf_sample_genotype_data[sample][tag])
                     else:
                        gt_tag = vcf_sample_genotype_data[sample][tag]
                  line_elements.append(gt_tag)
                  if gt_tag == './.':
                     if not keep_rejected_calls is False:
                        out.write('\t'.join(line_elements) + '\n')
                  else:
                     out.write('\t'.join(line_elements) + '\n')
            else:
               tsv_elements.append('\t'.join(vcf_info_data))
               line_elements = []
               line_elements.extend(tsv_elements)
               out.write('\t'.join(line_elements) + '\n')
      else:
         if len(sample_columns_header) > 0:
            if skip_genotype_data is False:
               for s in sorted(vcf_sample_genotype_data.keys()):
                  sample = s
                  line_elements = []
                  line_elements.extend(tsv_elements)
                  line_elements.append(sample)
                  gt_tag = '.'
                  for tag in sorted(vcf_sample_genotype_data[sample].keys()):
                     if tag != 'GT':
                        line_elements.append(vcf_sample_genotype_data[sample][tag])
                     else:
                        gt_tag = vcf_sample_genotype_data[sample][tag]
                  line_elements.append(gt_tag)
                  #print str(line_elements)
                  out.write('\t'.join(line_elements) + '\n')
            else:
               line_elements = []
               line_elements.extend(tsv_elements)
               line_elements = tsv_elements
               out.write('\t'.join(line_elements) + '\n')
       
   out.close()

if __name__=="__main__": __main__()


   
