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
from pcgr import utils
from pcgr.utils import error_message, check_subprocess

def __main__():

   parser = argparse.ArgumentParser(description='Verify input data for CPSR')
   parser.add_argument('pcgr_dir',help='Docker location of PCGR base directory with accompanying data directory, e.g. /data')
   parser.add_argument('input_vcf', help='VCF input file with query variants (SNVs/InDels)')
   parser.add_argument('custom_list',help='Custom text file indicating target genes form panel 0 for screening and reporting')
   parser.add_argument('preserved_info_tags',help='Comma-separated string of VCF INFO tags in query VCF to be retained for output')
   parser.add_argument('vcf_validation',type=int, default=0,choices=[0,1], help="Perform VCF validation with Ensembl's vcf-validator")
   parser.add_argument('genome_assembly',help='grch37 or grch38')
   parser.add_argument('sample_id',help='CPSR sample_name')
   parser.add_argument('virtual_panel_id',type=str,help='virtual panel identifier(s)')
   parser.add_argument('diagnostic_grade_only', type=int, default=0, choices=[0,1], help="Green virtual panels only (Genomics England PanelApp)")
   parser.add_argument('--output_dir', dest='output_dir', help='Output directory')
   parser.add_argument("--debug", action="store_true", help="Print full commands to log")
   args = parser.parse_args()

   ret = validate_cpsr_input(args.pcgr_dir,
                             args.input_vcf,
                             args.custom_list,
                             args.preserved_info_tags,
                             args.vcf_validation,
                             args.genome_assembly,
                             args.sample_id,
                             args.virtual_panel_id,
                             args.diagnostic_grade_only,
                             args.output_dir,
                             args.debug)
   if ret != 0:
      sys.exit(1)

# def is_valid_custom_bed(bed_file, logger):
#    """
#    Function that checks whether the custom panel (BED) adheres to the correct format
#    """
#    bed_reader = csv.DictReader(open(bed_file,'r'), delimiter='\t')
#    for row in bed_reader:
#       if len(row) != 4:
#          err_msg = 'BED file with custom screening regions must contain four columns: \'Chromosome\', \'Start\',\'End\',\'GeneSymbol\' - found entry containing ' + len(row) + ' columns'
#          return error_message(err_msg, logger)
#    bed_reader = csv.DictReader(open(bed_file,'r'), delimiter='\t', fieldnames=['Chromosome','Start','End','Symbol'])
#    bed_dataframe = np.read_csv(bed_file, usecols = [0,1,2,3], sep="\t",names=["Chromosome", "Start", "End","Symbol"])
#    if not bed_dataframe['Start'].dtype.kind in 'i': ## check that 'Start' is of type integer
#       err_msg = '\'Start\' column of BED file (custom panel) contains non-integer values'
#       return error_message(err_msg, logger)
#    if not bed_dataframe['End'].dtype.kind in 'i': ## check that 'End' is of type integer
#       err_msg = '\'End\' column of BED file (custom panel) contains non-integer values'
#       return error_message(err_msg, logger)
#    for rec in bed_reader:
#       if int(rec['End']) < int(rec['Start']): ## check that 'End' is always greather than 'Start'
#          err_msg = 'Detected wrongly formatted BED segment - \'Start\' is greater than \'End\' (' + str(rec['Chromosome']) + ':' + str(rec['Start']) + '-' + str(rec['End']) + ')'
#          return error_message(err_msg, logger)
#       if int(rec['End']) < 1 or int(rec['Start']) < 0: ## check that 'Start' and 'End' is always non-negative
#          err_msg = 'Detected wrongly formatted BED segment - \'Start\' or \'End\' is less than or equal to zero (' + str(rec['Chromosome']) + ':' + str(rec['Start']) + '-' + str(rec['End']) + ')'
#          return error_message(err_msg, logger)
#    logger.info('Custom panel BED file (' + str(bed_file) + ') adheres to the correct format (gene symbols not checked)')

#    return 0

def get_valid_custom_genelist(genelist_fname, genelist_bed_fname, pcgr_dir, genome_assembly, logger, debug):
   """
   Function that checks whether the custom genelist contains valid entries from the complete exploratory track
   """
   genelist_reader = csv.DictReader(open(genelist_fname,'r'), delimiter='\n', fieldnames=['ensembl_gene_id'])
   superpanel_track_bed = os.path.join(pcgr_dir, "data", genome_assembly, "virtual_panels",  "0." + genome_assembly + ".bed.gz")
   superpanel_track_tsv = os.path.join(pcgr_dir, "data", genome_assembly, "virtual_panels", "cpsr_superpanel." + genome_assembly + ".tsv")
   genelist_bed_fname_unsorted = genelist_bed_fname + '.tmp_unsorted'

   customlist_identifiers = {}
   superpanel_track = []
   superpanel_identifiers_all = {}
   valid_custom_identifiers = []
   valid_custom_symbols = []

   for row in genelist_reader:
      if not re.match(r'^ENSG[0-9]{1,}$',str(row['ensembl_gene_id']).rstrip()):
         err_msg = "Custom list of genes from CPSR superpanel (panel 0) should be provided as Ensembl gene identifiers, '" + str(row['ensembl_gene_id']) + "' is not a valid identifier"
         return error_message(err_msg, logger)
      else:
         customlist_identifiers[str(row['ensembl_gene_id']).strip()] = 1

   superpanel_reader = csv.DictReader(open(superpanel_track_tsv, 'r'), delimiter = '\t')

   for row in superpanel_reader:
      superpanel_track.append(dict(row))
   #superpanel_track = list(set(superpanel_track))

   i = 0
   while i < len(superpanel_track):
      superpanel_identifiers_all[superpanel_track[i]['ensembl_gene_id']] = superpanel_track[i]['symbol']
      i = i + 1

   for g in customlist_identifiers.keys():
      if g in superpanel_identifiers_all.keys():
         valid_custom_identifiers.append(g)
         valid_custom_symbols.append(superpanel_identifiers_all[g])
      else:
         logger.warning("Ignoring custom-provided gene identifier (" + str(g) + ") NOT found in CPSR superpanel (panel 0)")
         logger.warning("Choose only Ensembl gene identifiers from this set in data bundle: data/" + str(genome_assembly) + "/virtual_panels/cpsr_superpanel." + str(genome_assembly) + ".tsv")
   all_valid_custom_geneset = ', '.join(sorted(valid_custom_symbols))

   logger.info('Detected n = ' + str(len(valid_custom_identifiers)) + ' valid targets in custom-provided gene list file (--custom_list)):')
   logger.info(all_valid_custom_geneset)

   if len(valid_custom_identifiers) == 0:
      logger.info('')
      logger.info("NO valid gene identifiers from panel 0 in custom-provided genelist - exiting")
      logger.info('')
      exit(1)

   ## Add secondary findings genes to target BED
   cmd_secondary_regions_bed = 'bgzip -dc ' + str(superpanel_track_bed) + ' | egrep \'\|ACMG_SF30\|\' > ' + str(genelist_bed_fname_unsorted)
   check_subprocess(logger, cmd_secondary_regions_bed, debug)

   ## Add GWAS hits to target BED
   cmd_gwas_regions_bed = 'bgzip -dc ' + str(superpanel_track_bed) + ' | egrep \'rs[0-9]{3,}\|\' >> ' + str(genelist_bed_fname_unsorted)
   check_subprocess(logger, cmd_gwas_regions_bed, debug)

   ## Add custom set genes to target BED
   logger.info('Creating BED file with custom target genes: ' + str(genelist_bed_fname))
   for g in valid_custom_identifiers:
      cmd_target_regions_bed = 'bgzip -dc ' + str(superpanel_track_bed) + ' | egrep \'\|' + g + '\|\' >> ' + str(genelist_bed_fname_unsorted)
      check_subprocess(logger, cmd_target_regions_bed, debug)

   ## Sort regions in target BED
   if os.path.exists(genelist_bed_fname_unsorted) and os.stat(genelist_bed_fname_unsorted).st_size != 0:
      cmd_sort_custom_bed1 = 'egrep \'^[0-9]\' ' + str(genelist_bed_fname_unsorted) + ' | sort -k1,1n -k2,2n -k3,3n > ' + str(genelist_bed_fname)
      cmd_sort_custom_bed2 = 'egrep -v \'^[0-9]\' ' + str(genelist_bed_fname_unsorted) + ' | egrep \'^[XYM]\' | sort -k1,1 -k2,2n -k3,3n >> ' + str(genelist_bed_fname)
      cmd_rm_unsorted = 'rm -f ' + str(genelist_bed_fname_unsorted)

      check_subprocess(logger, cmd_sort_custom_bed1, debug)
      check_subprocess(logger, cmd_sort_custom_bed2, debug)
      if not debug:
         check_subprocess(logger, cmd_rm_unsorted, debug)
   #else:
      #print('balle')

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
               return error_message(err_msg, logger)

   logger.info('No query VCF INFO tags coincide with CPSR INFO tags')
   return ret

def check_preserved_vcf_info_tags(input_vcf, preserved_info_tags, logger):

   """
   Function that compares the INFO tags in the query VCF and preserved INFO tags set by the user as retained in CPSR output TSV
   If any preserved tag is not in query VCF, an error will be returned
   """

   tags = str(preserved_info_tags).split(',')
   info_elements_query_vcf = []

   vcf = VCF(input_vcf)
   logger.info('Checking if existing INFO tags of query VCF file matches preserved INFO tags set by the user')
   ret = 1
   for e in vcf.header_iter():
      header_element = e.info()
      if 'ID' in header_element.keys() and 'HeaderType' in header_element.keys():
         if header_element['HeaderType'] == 'INFO':
            info_elements_query_vcf.append(header_element['ID'])


   for t in tags:
      if not t in info_elements_query_vcf:
         err_msg = "Preserved INFO tag '" + str(t) + "' not found among INFO tags in query VCF - make sure preserved VCF INFO tags are set correctly"
         return error_message(err_msg, logger)
      else:
         logger.info("Preserved INFO tag '" + str(t) + "' detected among INFO tags in query VCF")

   return ret

def simplify_vcf(input_vcf, vcf, custom_bed, pcgr_directory, genome_assembly, virtual_panel_id, sample_id, diagnostic_grade_only, output_dir, logger, debug):

   """
   Function that performs four separate checks/filters on the validated input VCF:
   1. Remove/Strip off any genotype data (not needed for annotation)
   2. If VCF have variants with multiple alternative alleles ("multiallelic", e.g. 'A,T'), these are decomposed into variants with a single alternative allele
   3. Filters against predisposition loci (virtual panel id or custom target)
   4. Final VCF file is sorted and indexed (bgzip + tabix)
   """

   input_vcf_cpsr_ready = os.path.join(output_dir, re.sub(r'(\.vcf$|\.vcf\.gz$)','.cpsr_ready.tmp.vcf', os.path.basename(input_vcf)))
   input_vcf_cpsr_ready_decomposed = os.path.join(output_dir, re.sub(r'(\.vcf$|\.vcf\.gz$)','.cpsr_ready.vcf', os.path.basename(input_vcf)))
   input_vcf_cpsr_ready_decomposed_target = os.path.join(output_dir, re.sub(r'(\.vcf$|\.vcf\.gz$)','.cpsr_ready_target.vcf', os.path.basename(input_vcf)))
   virtual_panels_tmp_bed = os.path.join(output_dir, "virtual_panels_all." + str(sample_id) + ".tmp.bed")
   virtual_panels_bed = os.path.join(output_dir, "virtual_panels_all." + str(sample_id) + ".bed")

   multiallelic_list = list()
   for rec in vcf:
      POS = rec.start + 1
      alt = ",".join(str(n) for n in rec.ALT)
      if len(rec.ALT) > 1:
         variant_id = f"{rec.CHROM}:{POS}_{rec.REF}->{alt}"
         multiallelic_list.append(variant_id)

   is_gzipped = True if input_vcf.endswith('.gz') else False
   cat_vcf = f"bgzip -dc {input_vcf}" if is_gzipped else "cat {input_vcf}"

   command_vcf_sample_free1 = f'{cat_vcf} | egrep \'^##\' > {input_vcf_cpsr_ready}'
   command_vcf_sample_free2 = f'{cat_vcf} | egrep \'^#CHROM\' >> {input_vcf_cpsr_ready}'
   command_vcf_sample_free3 = f'{cat_vcf} | egrep -v \'^#\' | sed \'s/^chr//\' | egrep \'^[0-9]\' | sort -k1,1n -k2,2n -k4,4 -k5,5 >> {input_vcf_cpsr_ready}'
   command_vcf_sample_free4 = f'{cat_vcf} | egrep -v \'^#\' | sed \'s/^chr//\' | egrep -v \'^[0-9]\' | egrep \'^[XYM]\' | sort -k1,1 -k2,2n -k4,4 -k5,5 >> {input_vcf_cpsr_ready}'
   command_vcf_sample_free5 = f'{cat_vcf} | egrep -v \'^#\' | sed \'s/^chr//\' | egrep -v \'^[0-9]\' | egrep -v \'^[XYM]\' | sort -k1,1 -k2,2n -k4,4 -k5,5 >> {input_vcf_cpsr_ready}'

   check_subprocess(logger, command_vcf_sample_free1, debug)
   check_subprocess(logger, command_vcf_sample_free2, debug)
   check_subprocess(logger, command_vcf_sample_free3, debug)
   check_subprocess(logger, command_vcf_sample_free4, debug)
   check_subprocess(logger, command_vcf_sample_free5, debug)

   if multiallelic_list:
      logger.warning(f"There were {len(multiallelic_list)} multiallelic sites detected. Showing (up to) the first 100:")
      print('----')
      print(', '.join(multiallelic_list[:100]))
      print('----')
      logger.info('Decomposing multi-allelic sites in input VCF file using \'vt decompose\'')
      command_decompose = f'vt decompose -s {input_vcf_cpsr_ready} > {input_vcf_cpsr_ready_decomposed}  2> {os.path.join(output_dir, "decompose.log")}'
      check_subprocess(logger, command_decompose, debug)
   else:
      command_copy = f'cp {input_vcf_cpsr_ready} {input_vcf_cpsr_ready_decomposed}'
      check_subprocess(logger, command_copy, debug)


   if not custom_bed == 'None':
      logger.info('Limiting variant set to user-defined screening loci (custom list from panel 0)')
      if os.path.exists(custom_bed) and os.stat(custom_bed).st_size != 0:
         target_variants_intersect_cmd = "bedtools intersect -wa -u -header -a " + str(input_vcf_cpsr_ready_decomposed) + " -b " + str(custom_bed) + " > " + str(input_vcf_cpsr_ready_decomposed_target)
         check_subprocess(logger, target_variants_intersect_cmd, debug)
      else:
         logger.info('Custom BED file has a filesize of zero or does not exist')
   else:
      logger.info('Limiting variant set to cancer predisposition loci - virtual panel id(s): ' + str(virtual_panel_id))

      ## Concatenate all panel BEDs to one big virtual panel BED, sort and make unique
      panel_ids = str(virtual_panel_id).split(',')
      for pid in panel_ids:
         target_bed_gz = os.path.join(pcgr_directory,'data',genome_assembly, 'virtual_panels', str(pid) + "." + genome_assembly + ".bed.gz")
         if diagnostic_grade_only == 1 and virtual_panel_id != 0:
            logger.info('Considering diagnostic-grade only genes in panel ' + str(pid) + ' - (GREEN status in Genomics England PanelApp)')
            target_bed_gz = os.path.join(pcgr_directory, 'data', genome_assembly, 'virtual_panels', str(pid) + "." + genome_assembly + ".GREEN.bed.gz")
         check_subprocess(logger, f'bgzip -dc {target_bed_gz} >> {virtual_panels_tmp_bed}', debug)

      ## sort the collection of virtual panels
      if os.path.exists(virtual_panels_tmp_bed) and os.stat(virtual_panels_tmp_bed).st_size != 0:
         cmd_sort_virtual_panel_bed1 = 'egrep \'^[0-9]\' ' + str(virtual_panels_tmp_bed) + ' | sort -k1,1n -k2,2n -k3,3n | uniq > ' + str(virtual_panels_bed)
         cmd_sort_virtual_panel_bed2 = 'egrep -v \'^[0-9]\' ' + str(virtual_panels_tmp_bed) + ' | egrep \'^[XYM]\' | sort -k1,1 -k2,2n -k3,3n | uniq >> ' + str(virtual_panels_bed)
         cmd_rm_unsorted_vp = 'rm -f ' + str(virtual_panels_tmp_bed)
         check_subprocess(logger, cmd_sort_virtual_panel_bed1, debug)
         check_subprocess(logger, cmd_sort_virtual_panel_bed2, debug)
         if not debug:
            check_subprocess(logger, cmd_rm_unsorted_vp, debug)

         if os.path.exists(virtual_panels_bed):
            target_variants_intersect_cmd = 'bedtools intersect -wa -u -header -a ' + str(input_vcf_cpsr_ready_decomposed) + ' -b ' + str(virtual_panels_bed) + ' > ' + str(input_vcf_cpsr_ready_decomposed_target)
            check_subprocess(logger, target_variants_intersect_cmd, debug)


   check_subprocess(logger, f'bgzip -cf {input_vcf_cpsr_ready_decomposed_target} > {input_vcf_cpsr_ready_decomposed_target}.gz', debug)
   check_subprocess(logger, f'tabix -p vcf {input_vcf_cpsr_ready_decomposed_target}.gz', debug)
   if not debug:
      check_subprocess(logger, f'rm -f {input_vcf_cpsr_ready} {virtual_panels_bed} {input_vcf_cpsr_ready_decomposed} {os.path.join(output_dir, "decompose.log")}', debug)

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

def validate_cpsr_input(pcgr_directory, input_vcf, custom_list_fname, preserved_info_tags, vcf_validation, genome_assembly, sample_id, virtual_panel_id, diagnostic_grade_only, output_dir, debug):
   """
   Function that reads the input files to CPSR (VCF file) and performs the following checks:
   1. Check that VCF file is properly formatted (according to EBIvariation/vcf-validator - VCF v4.2) - optional (vcf_validation in config file)
   2. Check that no INFO annotation tags in the query VCF coincides with those generated by CPSR
   3. Check that custom VCF INFO tags set by user as retained for output is found in query VCF
   4. Check that if VCF have variants with multiple alternative alleles (e.g. 'A,T') run vt decompose
   5. Check that VCF contains a single sample column
   6. The resulting VCF file is sorted and indexed (bgzip + tabix)
   """
   logger = utils.getlogger('cpsr-validate-input-arguments')

   custom_list_bed_fname = 'None'
   if not custom_list_fname == 'None':
      logger.info('Establishing BED track with custom list of genes from panel 0')
      custom_list_bed_fname = os.path.join(output_dir, sample_id + '.cpsr.' + genome_assembly + '.custom_list.bed')
      get_valid_custom_genelist(custom_list_fname, custom_list_bed_fname, pcgr_directory, genome_assembly, logger, debug)

   #config_options = annoutils.read_config_options(configuration_file, pcgr_directory, genome_assembly, logger, wflow = 'cpsr')
   if not input_vcf == 'None':
      if vcf_validation == 1:
         logger.info('Skipping validation of VCF file (deprecated as of Dec 2021)')
      else:
         logger.info('Skipping validation of VCF file as provided by option --no_vcf_validate')

      tag_check = check_existing_vcf_info_tags(input_vcf, pcgr_directory, genome_assembly, logger)
      if tag_check == -1:
         return -1

      if preserved_info_tags != "None":
         custom_check = check_preserved_vcf_info_tags(input_vcf, preserved_info_tags, logger)
         if custom_check == -1:
            return -1

      vcf = VCF(input_vcf)
      samples = vcf.samples
      if len(samples) > 1:
         err_msg = "Query VCF contains more than one sample column (" + ', '.join(samples) + ") - CPSR expects a germline VCF with a single sample column - exiting"
         return error_message(err_msg, logger)
      simplify_vcf(input_vcf, vcf, custom_list_bed_fname, pcgr_directory, genome_assembly, virtual_panel_id, sample_id, diagnostic_grade_only, output_dir, logger, debug)

   return 0

if __name__=="__main__":
    __main__()
