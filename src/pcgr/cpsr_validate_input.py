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

logger = logging.getLogger('cpsr-validate-input-arguments')
global debug

def __main__():

   parser = argparse.ArgumentParser(description='Verify input data for CPSR')
   parser.add_argument('pcgr_dir',help='Docker location of PCGR base directory with accompanying data directory, e.g. /data')
   parser.add_argument('input_vcf', help='VCF input file with query variants (SNVs/InDels)')
   parser.add_argument('custom_list',help='Custom text file indicating target genes form panel 0 for screening and reporting')
   parser.add_argument('configuration_file', help='Configuration file (TOML-formatted, e.g. pcgr_conf.toml)')
   parser.add_argument('vcf_validation',type=int, default=0,choices=[0,1], help="Perform VCF validation with Ensembl's vcf-validator")
   parser.add_argument('genome_assembly',help='grch37 or grch38')
   parser.add_argument('sample_id',help='CPSR sample_name')
   parser.add_argument('virtual_panel_id',type=int,help='virtual panel identifier')
   parser.add_argument('diagnostic_grade_only', type=int, default=0, choices=[0,1], help="Green virtual panels only (Genomics England PanelApp)")
   parser.add_argument('--output_dir', dest='output_dir', help='Output directory', default='/workdir/output')
   parser.add_argument('--debug',action='store_true',default=False, help='Print full docker commands to log')
   args = parser.parse_args()

   global debug
   debug = args.debug

   ret = validate_cpsr_input(args.pcgr_dir, args.input_vcf, args.custom_list, args.configuration_file, args.vcf_validation, args.genome_assembly, args.sample_id, args.virtual_panel_id, args.diagnostic_grade_only, args.output_dir)
   if ret != 0:
      sys.exit(1)

def check_subprocess(command):
   if debug:
      logger.info(command)
   try:
      output = subprocess.check_output(str(command), stderr=subprocess.STDOUT, shell=True)
      if len(output) > 0:
         print (str(output.decode()).rstrip())
   except subprocess.CalledProcessError as e:
      print (e.output.decode())
      exit(0)


def is_valid_custom_bed(bed_file, logger):
   """
   Function that checks whether the custom panel (BED) adheres to the correct format
   """
   
   
   bed_reader = csv.DictReader(open(bed_file,'r'), delimiter='\t')
   for row in bed_reader:
      if len(row) != 4:
         err_msg = 'BED file with custom screening regions must contain four columns: \'Chromosome\', \'Start\',\'End\',\'GeneSymbol\' - found entry containing ' + len(row) + ' columns'
         return annoutils.error_message(err_msg, logger)
   
   bed_reader = csv.DictReader(open(bed_file,'r'), delimiter='\t', fieldnames=['Chromosome','Start','End','Symbol'])
   
   bed_dataframe = np.read_csv(bed_file, usecols = [0,1,2,3], sep="\t",names=["Chromosome", "Start", "End","Symbol"])
   if not bed_dataframe['Start'].dtype.kind in 'i': ## check that 'Start' is of type integer
      err_msg = '\'Start\' column of BED file (custom panel) contains non-integer values'
      return annoutils.error_message(err_msg, logger)
   if not bed_dataframe['End'].dtype.kind in 'i': ## check that 'End' is of type integer
      err_msg = '\'End\' column of BED file (custom panel) contains non-integer values'
      return annoutils.error_message(err_msg, logger)

   for rec in bed_reader:
      if int(rec['End']) < int(rec['Start']): ## check that 'End' is always greather than 'Start'
         err_msg = 'Detected wrongly formatted BED segment - \'Start\' is greater than \'End\' (' + str(rec['Chromosome']) + ':' + str(rec['Start']) + '-' + str(rec['End']) + ')'
         return annoutils.error_message(err_msg, logger)
      if int(rec['End']) < 1 or int(rec['Start']) < 0: ## check that 'Start' and 'End' is always non-negative
         err_msg = 'Detected wrongly formatted BED segment - \'Start\' or \'End\' is less than or equal to zero (' + str(rec['Chromosome']) + ':' + str(rec['Start']) + '-' + str(rec['End']) + ')'
         return annoutils.error_message(err_msg, logger)
   logger.info('Custom panel BED file (' + str(bed_file) + ') adheres to the correct format (gene symbols not checked)')

   return 0

def get_valid_custom_genelist(genelist_fname, genelist_bed_fname, sample_id, pcgr_dir, output_dir, genome_assembly, logger):
   """
   Function that checks whether the custom genelist contains valid entries from the complete exploratory track
   """
   
   genelist_reader = csv.DictReader(open(genelist_fname,'r'), delimiter='\n', fieldnames=['symbol'])
   exploratory_track_bed = os.path.join(pcgr_dir,'data',genome_assembly, 'virtual_panels',  "0." + genome_assembly + ".bed.gz")
   exploratory_track_tsv = os.path.join(pcgr_dir,'data',genome_assembly, 'cpsr_exploratory_track.tsv')
   genelist_bed_fname_unsorted = genelist_bed_fname + '.tmp_unsorted'

   customlist_symbols = {}
   exploratory_track = []
   exploratory_symbols_all = {}
   valid_custom_symbols = []

   for row in genelist_reader:
      customlist_symbols[str(row['symbol']).strip()] = 1
      
   exploratory_reader = csv.DictReader(open(exploratory_track_tsv, 'r'), delimiter = '\t')

   for row in exploratory_reader:
      exploratory_track.append(dict(row))

   i = 0
   while i < len(exploratory_track):
      exploratory_symbols_all[exploratory_track[i]['symbol']] = 1
      i = i + 1

   for g in customlist_symbols.keys():
      if g in exploratory_symbols_all.keys():
         valid_custom_symbols.append(g)
      else:
         logger.warning("Ignoring custom-provided gene symbol (" + str(g) + ") NOT found in CPSR exploratory track 0")
         logger.warning("Choose only symbols from this set: https://github.com/sigven/cpsr/predisposition.md")
   all_valid_str = ', '.join(sorted(valid_custom_symbols))

   logger.info('Detected n = ' + str(len(valid_custom_symbols)) + ' valid targets in custom-provided gene list file (--custom_list)):')
   logger.info(all_valid_str)

   if len(valid_custom_symbols) == 0:
      logger.info('')
      logger.info("NO valid gene symbols from panel 0 in custom-provided genelist - exiting")
      logger.info('')
      exit(1)

   logger.info('Creating BED file with custom target genes: ' + str(genelist_bed_fname))
   for g in valid_custom_symbols:
      cmd_target_regions_bed = 'bgzip -dc ' + str(exploratory_track_bed) + ' | egrep \'\|' + g + '\|\' >> ' + str(genelist_bed_fname_unsorted)
      check_subprocess(cmd_target_regions_bed)

   if os.path.exists(genelist_bed_fname_unsorted) and os.stat(genelist_bed_fname_unsorted).st_size != 0:
      cmd_sort_custom_bed1 = 'egrep \'^[0-9]\' ' + str(genelist_bed_fname_unsorted) + ' | sort -k1,1n -k2,2n -k3,3n > ' + str(genelist_bed_fname)
      cmd_sort_custom_bed2 = 'egrep -v \'^[0-9]\' ' + str(genelist_bed_fname_unsorted) + ' | egrep \'^[XYM]\' | sort -k1,1 -k2,2n -k3,3n >> ' + str(genelist_bed_fname)
      cmd_rm_unsorted = 'rm -f ' + str(genelist_bed_fname_unsorted)

      check_subprocess(cmd_sort_custom_bed1)
      check_subprocess(cmd_sort_custom_bed2)
      check_subprocess(cmd_rm_unsorted)

   return 0


def check_existing_vcf_info_tags(input_vcf, pcgr_directory, genome_assembly, logger):

   """
   Function that compares the INFO tags in the query VCF and the INFO tags generated by PCGR
   If any coinciding tags, an error will be returned
   """

   pcgr_infotags_desc = annoutils.read_infotag_file(os.path.join(pcgr_directory,'data',genome_assembly, 'cpsr_infotags.tsv'))

   vcf = VCF(input_vcf)
   logger.info('Checking if existing INFO tags of query VCF file coincide with CPSR INFO tags')
   ret = 1
   for e in vcf.header_iter():
      header_element = e.info()
      if 'ID' in header_element.keys() and 'HeaderType' in header_element.keys():
         if header_element['HeaderType'] == 'INFO':
            if header_element['ID'] in pcgr_infotags_desc.keys():
               err_msg = 'INFO tag ' + str(header_element['ID']) + ' in the query VCF coincides with a VCF annotation tag produced by CPSR - please remove or rename this tag in your query VCF'
               return annoutils.error_message(err_msg, logger)

   logger.info('No query VCF INFO tags coincide with CPSR INFO tags')
   return ret

def simplify_vcf(input_vcf, vcf, custom_bed, pcgr_directory, genome_assembly, virtual_panel_id, diagnostic_grade_only, output_dir, logger):

   """
   Function that performs tre things on the validated input VCF:
   1. Strip of any genotype data
   2. If VCF have variants with multiple alternative alleles ("multiallelic", e.g. 'A,T'), these are decomposed into variants with a single alternative allele
   3. Filters against predisposition loci (virtual panel id or custom target)
   4. Final VCF file is sorted and indexed (bgzip + tabix)
   """

   input_vcf_cpsr_ready = os.path.join(output_dir, re.sub(r'(\.vcf$|\.vcf\.gz$)','.cpsr_ready.tmp.vcf', os.path.basename(input_vcf)))
   input_vcf_cpsr_ready_decomposed = os.path.join(output_dir, re.sub(r'(\.vcf$|\.vcf\.gz$)','.cpsr_ready.vcf', os.path.basename(input_vcf)))
   input_vcf_cpsr_ready_decomposed_target = os.path.join(output_dir, re.sub(r'(\.vcf$|\.vcf\.gz$)','.cpsr_ready_target.vcf', os.path.basename(input_vcf)))

   multiallelic_alt = 0
   for rec in vcf:
      POS = rec.start + 1
      alt = ",".join(str(n) for n in rec.ALT)
      if len(rec.ALT) > 1:
         logger.warning("Multiallelic site detected:" + str(rec.CHROM) + '\t' + str(POS) + '\t' + str(rec.REF) + '\t' + str(alt))
         multiallelic_alt = 1
   command_vcf_sample_free1 = 'egrep \'^##\' ' + str(input_vcf) + ' > ' + str(input_vcf_cpsr_ready)
   command_vcf_sample_free2 = 'egrep \'^#CHROM\' ' + str(input_vcf) + ' >> ' + str(input_vcf_cpsr_ready)
   command_vcf_sample_free3 = 'egrep -v \'^#\' ' + str(input_vcf) + ' | sed \'s/^chr//\' | egrep \'^[0-9]\' | sort -k1,1n -k2,2n -k4,4 -k5,5 >> ' + str(input_vcf_cpsr_ready)
   command_vcf_sample_free4 = 'egrep -v \'^#\' ' + str(input_vcf) + ' | sed \'s/^chr//\' | egrep -v \'^[0-9]\' | egrep \'^[XYM]\' | sort -k1,1 -k2,2n -k4,4 -k5,5 >> ' + str(input_vcf_cpsr_ready)
   command_vcf_sample_free5 = 'egrep -v \'^#\' ' + str(input_vcf) + ' | sed \'s/^chr//\' | egrep -v \'^[0-9]\' | egrep -v \'^[XYM]\' | sort -k1,1 -k2,2n -k4,4 -k5,5 >> ' + str(input_vcf_cpsr_ready)
   
   if input_vcf.endswith('.gz'):
      command_vcf_sample_free1 = 'bgzip -dc ' + str(input_vcf) + ' | egrep \'^##\' > ' + str(input_vcf_cpsr_ready)
      command_vcf_sample_free2 = 'bgzip -dc ' + str(input_vcf) + ' | egrep \'^#CHROM\'  >> ' + str(input_vcf_cpsr_ready)
      command_vcf_sample_free3 = 'bgzip -dc ' + str(input_vcf) + ' | egrep -v \'^#\' | sed \'s/^chr//\' | egrep \'^[0-9]\' | sort -k1,1n -k2,2n -k4,4 -k5,5 >> ' + str(input_vcf_cpsr_ready)
      command_vcf_sample_free4 = 'bgzip -dc ' + str(input_vcf) + ' | egrep -v \'^#\' | sed \'s/^chr//\' | egrep -v \'^[0-9]\' | egrep \'^[XYM]\' | sort -k1,1 -k2,2n -k4,4 -k5,5 >> ' + str(input_vcf_cpsr_ready)
      command_vcf_sample_free5 = 'bgzip -dc ' + str(input_vcf) + ' | egrep -v \'^#\' | sed \'s/^chr//\' | egrep -v \'^[0-9]\' | egrep -v \'^[XYM]\' | sort -k1,1 -k2,2n -k4,4 -k5,5 >> ' + str(input_vcf_cpsr_ready)

   check_subprocess(command_vcf_sample_free1)
   check_subprocess(command_vcf_sample_free2)
   check_subprocess(command_vcf_sample_free3)
   check_subprocess(command_vcf_sample_free4)
   check_subprocess(command_vcf_sample_free5)
  
   if multiallelic_alt == 1:
      logger.info('Decomposing multi-allelic sites in input VCF file using \'vt decompose\'')
      command_decompose = 'vt decompose -s ' + str(input_vcf_cpsr_ready) + ' > ' + str(input_vcf_cpsr_ready_decomposed) + ' 2> ' + os.path.join(output_dir, 'decompose.log')
      check_subprocess(command_decompose)
   else:
      command_copy = 'cp ' + str(input_vcf_cpsr_ready) + ' ' + str(input_vcf_cpsr_ready_decomposed)
      check_subprocess(command_copy)


   if not custom_bed == 'None':
      logger.info('Limiting variant set to user-defined screening loci (custom list from panel 0)')
      if os.path.exists(custom_bed) and os.stat(custom_bed).st_size != 0:
         target_variants_intersect_cmd = "bedtools intersect -wa -u -header -a " + str(input_vcf_cpsr_ready_decomposed) + " -b " + str(custom_bed) + " > " + str(input_vcf_cpsr_ready_decomposed_target)
         check_subprocess(target_variants_intersect_cmd)
      else:
         logger.info('Custom BED file has a filesize of zero or does not exist')
   else:
      logger.info('Limiting variant set to cancer predisposition loci - virtual panel_id ' + str(virtual_panel_id))
      target_bed = os.path.join(pcgr_directory,'data',genome_assembly, 'virtual_panels', str(virtual_panel_id) + "." + genome_assembly + ".bed.gz")
      if diagnostic_grade_only == 1 and virtual_panel_id != 0:
         target_bed = os.path.join(pcgr_directory, 'data', genome_assembly, 'virtual_panels', str(virtual_panel_id) + "." + genome_assembly + ".GREEN.bed.gz")
      if os.path.exists(target_bed):
         target_variants_intersect_cmd = 'bedtools intersect -wa -u -header -a ' + str(input_vcf_cpsr_ready_decomposed) + ' -b ' + str(target_bed) + ' > ' + str(input_vcf_cpsr_ready_decomposed_target)
         check_subprocess(target_variants_intersect_cmd)
      else:
         logger.info('')
         logger.info("Virtual gene panel - BED file (" + str(target_bed) + ")  not found - exiting")
         logger.info('')
         exit(1)
   
   check_subprocess('bgzip -cf ' + str(input_vcf_cpsr_ready_decomposed_target) + ' > ' + str(input_vcf_cpsr_ready_decomposed_target) + '.gz')
   check_subprocess('tabix -p vcf ' + str(input_vcf_cpsr_ready_decomposed_target) + '.gz')
   if not debug:
      check_subprocess('rm -f ' + str(input_vcf_cpsr_ready) + ' ' + str(input_vcf_cpsr_ready_decomposed) + ' ' + os.path.join(output_dir, 'decompose.log'))

   if os.path.exists(input_vcf_cpsr_ready_decomposed_target + '.gz') and os.path.getsize(input_vcf_cpsr_ready_decomposed_target + '.gz') > 0:
      vcf = VCF(input_vcf_cpsr_ready_decomposed_target + '.gz')
      i = 0
      for rec in vcf:
         i = i + 1
      if len(vcf.seqnames) == 0 or i == 0:
         logger.info('')
         logger.info("Query VCF contains NO variants within the selected cancer predisposition geneset or ACMG-recommended genes for secondary findings - quitting workflow")
         logger.info('')
         exit(1)

def validate_cpsr_input(pcgr_directory, input_vcf, custom_list_fname, configuration_file, vcf_validation, genome_assembly, sample_id, virtual_panel_id, diagnostic_grade_only, output_dir):
   """
   Function that reads the input files to CPSR (VCF file) and performs the following checks:
   1. Check that VCF file is properly formatted (according to EBIvariation/vcf-validator - VCF v4.2) - optional (vcf_validation in config file)
   2. Check that no INFO annotation tags in the query VCF coincides with those generated by CPSR
   3. Check that if VCF have variants with multiple alternative alleles (e.g. 'A,T') run vt decompose
   4. Check that VCF contains a single sample column 
   5. The resulting VCF file is sorted and indexed (bgzip + tabix)
   """
   logger = annoutils.getlogger('cpsr-validate-input')

   custom_list_bed_fname = 'None'
   if not custom_list_fname == 'None':
      custom_list_bed_fname = os.path.join(output_dir, sample_id + '.cpsr.' + genome_assembly + '.custom_list.bed') 
      get_valid_custom_genelist(custom_list_fname, custom_list_bed_fname, sample_id, pcgr_directory, output_dir, genome_assembly, logger)

   #config_options = annoutils.read_config_options(configuration_file, pcgr_directory, genome_assembly, logger, wflow = 'cpsr')
   if not input_vcf == 'None':
      if vcf_validation == 1:
         valid_vcf = annoutils.is_valid_vcf(input_vcf, output_dir, logger, debug)
         if valid_vcf == -1:
            return -1
      else:
         logger.info('Skipping validation of VCF file - as provided by option --no_vcf_validate')

      tag_check = check_existing_vcf_info_tags(input_vcf, pcgr_directory, genome_assembly, logger)
      if tag_check == -1:
         return -1

      vcf = VCF(input_vcf)
      samples = vcf.samples
      if len(samples) > 1:
         err_msg = "Query VCF contains more than one sample column (" + ', '.join(samples) + ") - CPSR expects a germline VCF with a single sample column - exiting"
         return annoutils.error_message(err_msg, logger)
      simplify_vcf(input_vcf, vcf, custom_list_bed_fname, pcgr_directory, genome_assembly, virtual_panel_id, diagnostic_grade_only, output_dir, logger)

   return 0

if __name__=="__main__": __main__()
