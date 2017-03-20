#!/usr/bin/env python

import argparse
import vcf
import re
import cyvcf


def __main__():
   parser = argparse.ArgumentParser(description='Convert a VCF file with genomic variants to a file with tab-separated values (TSV). Sample genotypes (if present) will be printed on separate lines. Note: Only non-rejected calls/genotypes will be printed in output', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
   parser.add_argument('query_vcf', help='Bgzipped input VCF file with query variants (SNVs/InDels)')
   parser.add_argument('out_tsv', help='Output TSV file with one line pr non-rejected sample genotype (Variant, genotype and annotation data as tab-separated values)')
   args = parser.parse_args()
   
   vcf2tsv(args.query_vcf, args.out_tsv)


def get_tsv_header(vcf_reader):
   fixedcolumns = ['CHROM','POS','ID','REF','ALT','QUAL','FILTER']
   infocolumns = []
   samplecolumns = []
   for keyw in sorted(vcf_reader.infos.keys()):
      infocolumns.append(keyw.upper())
   if len(vcf_reader.formats.keys()) > 0:
      samplecolumns.append('VCF_SAMPLE_ID')
      for keyw in sorted(vcf_reader.formats.keys()):
         samplecolumns.append(keyw.upper())

   header_line = '\t'.join(fixedcolumns) + '\t' + '\t'.join(infocolumns)
   if len(samplecolumns) > 0:
      header_line = '\t'.join(fixedcolumns) + '\t' + '\t'.join(infocolumns) + '\t' + '\t'.join(samplecolumns)
   return header_line

def vcf2tsv(query_vcf, out_tsv):
   
   vcf_reader = cyvcf.Reader(open(query_vcf, 'r'))
   f_out = open(out_tsv,'w')
   header_line = get_tsv_header(vcf_reader)
   f_out.write(header_line + '\n')
   for rec in vcf_reader:
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
      infocolumn = []
      for keyw in sorted(vcf_reader.infos.keys()):
         if rec.INFO.has_key(keyw):
            #print str(keyw) + '\t' + str(rec.INFO[keyw])
            if vcf_reader.infos[keyw][2] == 'Flag':
               if rec.INFO[keyw] == 1 or rec.INFO[keyw] == True:
                  infocolumn.append('True')
               else:
                  infocolumn.append('False')
            elif vcf_reader.infos[str(keyw)][2] == 'Float' or vcf_reader.infos[str(keyw)][2] == 'Integer':
               if type(rec.INFO[str(keyw)]) is list:
                  infocolumn.append(",".join(str(n) for n in rec.INFO[str(keyw)]))
               else:
                  if re.search(r'^\[(.|\s)+\]$',str(rec.INFO[keyw])):
                     tmp = re.sub(r'\[|\]|\'|\s{1,}','',rec.INFO[keyw])
                     infocolumn.append(str(tmp))
                  else:
                     infocolumn.append(str(rec.INFO[str(keyw)]))
            else:
               if type(rec.INFO[keyw]) is list:
                  all_vals = []
                  for m in rec.INFO[str(keyw)]:
                     if m is None:
                        all_vals.append('NA')
                     else:
                        all_vals.append(m)
                  
                  infocolumn.append(str(','.join(all_vals)))
               else:
                  ## For some reason(?) some lists are parsed into string objects by PyVCF
                  if re.search(r'^\[(.|\s)+\]$',str(rec.INFO[keyw])):
                     #print str(rec.INFO[keyw])
                     tmp = re.sub(r'\[|\]|\'|\s{1,}','',rec.INFO[keyw])
                     #print str(tmp)
                     infocolumn.append(str(tmp))
                  else:
                     infocolumn.append(str(rec.INFO[str(keyw)]))
         else:
            infocolumn.append('NA')
      infostring = '\t'.join(infocolumn)
      i = 0
      if len(vcf_reader.samples) > 0:
         while i < len(vcf_reader.samples):
            genotype_values = []
            skip_null_sample = False
            depth_tumor_non_zero = False
            depth_control_non_zero = False
            for elem in sorted(rec.FORMAT.split(':')):
               p = rec.samples[int(i)][str(elem)]
               if isinstance(p,list):
                  genotype_values.append(','.join(map(str ,p)))
               else:
                  if not p:
                     if str(elem) == 'GT':
                        skip_null_sample = True
                        genotype_values.append('./.')
                     else:
                        if p is None:
                           p = '.'
                        genotype_values.append(str(p))
                  else:
                     if elem == 'GT' and p == './.':
                        skip_null_sample = True
                     genotype_values.append(str(p))
            sample_genotype_data = '\t'.join(genotype_values)
            if skip_null_sample is False:
               tsvline = fixed_fields_string + '\t' + infostring + '\t' + str(vcf_reader.samples[i]) + '\t' + sample_genotype_data + '\n'
               f_out.write(tsvline)
            i = i + 1
      else:
         tsvline = fixed_fields_string + '\t' + infostring + '\n'
         f_out.write(tsvline)
   f_out.close()

if __name__=="__main__": __main__()


   
