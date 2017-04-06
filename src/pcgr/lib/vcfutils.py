#!/usr/bin/env python
#Filename: vcfutils.py

import os
import datetime
import vcf
import cyvcf
import sys
import re
import random
import traceback
from operator import itemgetter, attrgetter
import pcgr

logger = pcgr.getlogger('vcfutils')

def print_vcf_meta(outfile,vcf_reader, print_sample_data=0, sample_id = None):
   out = open(outfile,'a')
   info_strings = []
   filter_strings = []
   format_strings = []
   meta_strings = []
   
   ### Meta information
   if not vcf_reader.metadata.has_key('fileformat'):
      meta_strings.append("##fileformat=VCFv4.2")
   else:
      meta_strings.append("##fileformat=" + str(vcf_reader.metadata['fileformat']))
   if not vcf_reader.metadata.has_key('reference'):
      meta_strings.append('##reference=GRCh37')
   for metatag in sorted(vcf_reader.metadata.keys()):
      if metatag != 'fileformat' and metatag != 'filedate':
         if type(vcf_reader.metadata[metatag]) is list:
            metaline = "##" + str(metatag) + '=' + str(vcf_reader.metadata[metatag][0])
         else:
            metaline = "##" + str(metatag) + '=' + str(vcf_reader.metadata[metatag])
         meta_strings.append(metaline)
    
   now = datetime.datetime.now()
   if not vcf_reader.metadata.has_key('filedate'):
      meta_strings.append('##filedate=' + now.strftime("%Y-%m-%d"))
    
   for info_tag in sorted(vcf_reader.infos.keys()):
      if str(vcf_reader.infos[info_tag][1]) == 'None':
         info_line = '##INFO=<ID=' + info_tag + ',' + 'Number=.,Type=' + str(vcf_reader.infos[info_tag][2]) + ',' + 'Description="' + str(vcf_reader.infos[info_tag][3]) + '">'
      else:
         info_line = '##INFO=<ID=' + info_tag + ',' + 'Number=' + str(vcf_reader.infos[info_tag][1]) + ',' + 'Type=' + str(vcf_reader.infos[info_tag][2]) + ',' + 'Description="' + str(vcf_reader.infos[info_tag][3]) + '">'
      info_strings.append(info_line)
   for filter_tag in sorted(vcf_reader.filters.keys()):
      filter_line = '##FILTER=<ID=' + filter_tag + ',' + 'Description="' + vcf_reader.filters[filter_tag][1] + '">'
      filter_strings.append(filter_line)
   if print_sample_data == 1:
      for format_tag in sorted(vcf_reader.formats.keys()):
         if str(vcf_reader.formats[format_tag][1]) == 'None':
            format_line = '##FORMAT=<ID=' + format_tag + ',' + 'Number=.,Type=' + str(vcf_reader.formats[format_tag][2]) + ',' + 'Description="' + str(vcf_reader.formats[format_tag][3]) + '">'
            format_strings.append(format_line)
         else:
            format_line = '##FORMAT=<ID=' + format_tag + ',' + 'Number=' + str(vcf_reader.formats[format_tag][1]) + ',' + 'Type=' + str(vcf_reader.formats[format_tag][2]) + ',' + 'Description="' + str(vcf_reader.formats[format_tag][3]) + '">'
            format_strings.append(format_line)
            
   if len(meta_strings) > 0:
      out.write('\n'.join(meta_strings))
      out.write('\n')
   if len(format_strings) > 0:
      out.write('\n'.join(format_strings))
      out.write('\n')
   if len(filter_strings) > 0:
      out.write('\n'.join(filter_strings))
      out.write('\n')
   if len(info_strings) > 0:
      out.write('\n'.join(info_strings))
      out.write('\n')
   out.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO')
   if print_sample_data == 1:
      if len(vcf_reader.formats.keys()) > 0:
         out.write('\tFORMAT\t')
      if len(vcf_reader.samples) > 0:
         if sample_id is None:
            out.write('\t'.join(vcf_reader.samples))
         else:
            if sample_id in vcf_reader.samples:
               out.write(sample_id)
   out.write('\n')
   out.close()

def print_vcf_content(outfile,vcf_lines):
   if len(vcf_lines) > 0:
      out = open(outfile,'a')
      out.write('\n'.join(vcf_lines))
      out.write('\n')

def get_vcf_fixed_columns(rec):
   alt = ",".join(str(n) for n in rec.ALT)
   qual = '.'
   rec_id = str(rec.ID)
   if rec.ID is None:
      rec_id = '.'
   if not rec.QUAL is None:
      qual = str(rec.QUAL)
   rec_filter = str(rec.FILTER)
   if type(rec.FILTER) is list:
      if len(rec.FILTER) == 0:
         rec_filter = 'PASS'
      elif len(rec.FILTER) == 1:
         rec_filter = str(rec.FILTER[0])
      else:
         rec_filter = str(';'.join(str(n) for n in rec.FILTER))
   else:
      if rec.FILTER is None:
         rec_filter = 'PASS'
   fixed_fields_string = str(rec.CHROM) + '\t' + str(rec.POS) + '\t' + str(rec_id) + '\t' + str(rec.REF) + '\t' + str(alt) + '\t' + str(qual) + '\t' + str(rec_filter)
   
   return fixed_fields_string


def get_info_value(rec, keyword, vcf_reader, existing_info_tags):

	if vcf_reader.infos[keyword][2] == 'Flag':
		if rec.INFO[keyword] == 1 or rec.INFO[keyword] == True:
			existing_info_tags[keyword] = True
	elif vcf_reader.infos[str(keyword)][2] == 'Float' or vcf_reader.infos[str(keyword)][2] == 'Integer':
		if type(rec.INFO[str(keyword)]) is list:
			existing_info_tags[keyword] = ",".join(str(n) for n in rec.INFO[str(keyword)])
		else:
			if re.search(r'^\[(.|\s)+\]$',str(rec.INFO[keyword])):
				tmp = re.sub(r'\[|\]|\'|\s{1,}','',rec.INFO[keyword])
				existing_info_tags[keyword] = str(tmp)
			else:
				existing_info_tags[keyword] = str(rec.INFO[str(keyword)])
	else:
		if type(rec.INFO[keyword]) is list:
			all_vals = []
			for m in rec.INFO[str(keyword)]:
				if m is None:
					all_vals.append('.')
				else:
					all_vals.append(m)
			
			existing_info_tags[keyword] = str(','.join(all_vals))
		else:
			## For some reason(?) some lists are parsed into string objects by PyVCF
			if re.search(r'^\[(.|\s)+\]$',str(rec.INFO[keyword])):
				tmp = re.sub(r'\[|\]|\'|\s{1,}','',rec.INFO[keyword])
				existing_info_tags[keyword] = str(tmp)
			else:
				existing_info_tags[keyword] = str(rec.INFO[str(keyword)])