#!/usr/bin/env python

import argparse
from cyvcf2 import VCF
import random
import annoutils
import os
import re
import sys

logger = annoutils.getlogger('pcgr-vcfanno')


def __main__():
   parser = argparse.ArgumentParser(description='Run brentp/vcfanno - annotate a VCF file against multiple VCF files in parallel', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
   parser.add_argument('query_vcf', help='Bgzipped input VCF file with query variants (SNVs/InDels)')
   parser.add_argument('out_vcf', help='Output VCF file with appended annotations from multiple VCF files')
   parser.add_argument('pcgr_db_dir', help='PCGR data directory')
   parser.add_argument('--num_processes', help="Number of processes vcfanno can use during annotation", default=4)
   parser.add_argument("--docm",action = "store_true", help="Annotate VCF with annotations from Database of Curated Mutations")
   parser.add_argument("--intogen_driver_mut",action = "store_true", help="Annotate VCF with predicted cancer driver mutations from IntoGen's Catalog of Driver Mutations")
   parser.add_argument("--clinvar",action = "store_true", help="Annotate VCF with annotations from ClinVar")
   parser.add_argument("--dbnsfp",action = "store_true", help="Annotate VCF with annotations from database of non-synonymous functional predictions")
   parser.add_argument("--tcga",action = "store_true", help="Annotate VCF with variant frequencies from the The Cancer Genome Atlas")
   parser.add_argument("--civic",action = "store_true", help="Annotate VCF with annotations from the Clinical Interpretation of Variants in Cancer database")
   parser.add_argument("--cbmdb",action = "store_true", help="Annotate VCF with annotations from the Cancer bioMarkers database")
   parser.add_argument("--icgc",action = "store_true", help="Annotate VCF with known variants found in the ICGC-PCAWG sequencing project")
   parser.add_argument("--cancer_hotspots",action = "store_true", help="Annotate VCF with mutation hotspots from cancerhotspots.org")
   parser.add_argument("--uniprot",action = "store_true", help="Annotate VCF with protein functional features from the UniProt Knowledgebase")
   parser.add_argument("--pcgr_onco_xref",action = "store_true", help="Annotate VCF with transcript annotations from PCGR (targeted drugs, protein complexes, cancer gene associations, etc)")
   
   args = parser.parse_args()
   query_info_tags = get_vcf_info_tags(args.query_vcf)
   vcfheader_file = args.out_vcf + '.tmp.' + str(random.randrange(0,10000000)) + '.header.txt'
   conf_fname = args.out_vcf + '.tmp.conf.toml'
   print_vcf_header(args.query_vcf, vcfheader_file, chromline_only = False)
   run_vcfanno(args.num_processes, args.query_vcf, query_info_tags, vcfheader_file, args.pcgr_db_dir, conf_fname, args.out_vcf, args.docm, args.intogen_driver_mut, args.clinvar, args.tcga, args.dbnsfp, args.civic, args.cbmdb, args.icgc, args.uniprot, args.cancer_hotspots, args.pcgr_onco_xref)


def prepare_vcfanno_configuration(vcfanno_data_directory, conf_fname, vcfheader_file, logger, datasource_info_tags, query_info_tags, datasource):
   for t in datasource_info_tags:
      if t in query_info_tags:
         logger.warning("Query VCF has INFO tag " + str(t) + ' - this is also present in the ' + str(datasource) + ' VCF/BED annotation file. This tag will be overwritten if not renamed in the query VCF')
   append_to_conf_file(datasource, datasource_info_tags, vcfanno_data_directory, conf_fname)
   append_to_vcf_header(vcfanno_data_directory, datasource, vcfheader_file)

def run_vcfanno(num_processes, query_vcf, query_info_tags, vcfheader_file, pcgr_db_directory, conf_fname, output_vcf, docm, intogen_driver_mut, clinvar, tcga, dbnsfp, civic, cbmdb, icgc, uniprot, cancer_hotspots, pcgr_onco_xref):
   """
   Function that annotates a VCF file with vcfanno against a user-defined set of germline and somatic VCF files
   """
   
   civic_info_tags = ["CIVIC_ID","CIVIC_ID_2"]
   cbmdb_info_tags = ["CBMDB_ID"]
   icgc_info_tags = ["ICGC_PCAWG_OCCURRENCE","ICGC_PCAWG_AFFECTED_DONORS"]
   docm_info_tags = ["DOCM_PMID"]
   tcga_info_tags = ["TCGA_FREQUENCY","TCGA_PANCANCER_COUNT"]
   intogen_driver_mut_info_tags = ["INTOGEN_DRIVER_MUT"]
   clinvar_info_tags = ["CLINVAR_MSID","CLINVAR_PMID","CLINVAR_CLNSIG","CLINVAR_VARIANT_ORIGIN","CLINVAR_CONFLICTED","CLINVAR_MEDGEN_CUI","CLINVAR_MEDGEN_CUI_SOMATIC","CLINVAR_CLNSIG_SOMATIC","CLINVAR_PMID_SOMATIC","CLINVAR_ALLELE_ID","CLINVAR_HGVSP"]
   cancer_hotspots_info_tags = ["CANCER_MUTATION_HOTSPOT"]
   dbnsfp_info_tags = ["DBNSFP"]
   uniprot_info_tags = ["UNIPROT_FEATURE"]
   pcgr_onco_xref_info_tags = ["PCGR_ONCO_XREF"]
      
   if icgc is True:
      prepare_vcfanno_configuration(pcgr_db_directory, conf_fname, vcfheader_file, logger, icgc_info_tags, query_info_tags, "icgc")
   if clinvar is True:
      prepare_vcfanno_configuration(pcgr_db_directory, conf_fname, vcfheader_file, logger, clinvar_info_tags, query_info_tags, "clinvar")
   if dbnsfp is True:
      prepare_vcfanno_configuration(pcgr_db_directory, conf_fname, vcfheader_file, logger, dbnsfp_info_tags, query_info_tags, "dbnsfp")
   if cbmdb is True:
      prepare_vcfanno_configuration(pcgr_db_directory, conf_fname, vcfheader_file, logger, cbmdb_info_tags, query_info_tags, "cbmdb")
   if tcga is True:
      prepare_vcfanno_configuration(pcgr_db_directory, conf_fname, vcfheader_file, logger, tcga_info_tags, query_info_tags, "tcga")
   if civic is True:
      prepare_vcfanno_configuration(pcgr_db_directory, conf_fname, vcfheader_file, logger, civic_info_tags, query_info_tags, "civic")
   if cancer_hotspots is True:
      prepare_vcfanno_configuration(pcgr_db_directory, conf_fname, vcfheader_file, logger, cancer_hotspots_info_tags, query_info_tags, "cancer_hotspots")
   if uniprot is True:
      prepare_vcfanno_configuration(pcgr_db_directory, conf_fname, vcfheader_file, logger, uniprot_info_tags, query_info_tags, "uniprot")
   if docm is True:
      prepare_vcfanno_configuration(pcgr_db_directory, conf_fname, vcfheader_file, logger, docm_info_tags, query_info_tags, "docm")
   if intogen_driver_mut is True:
      prepare_vcfanno_configuration(pcgr_db_directory, conf_fname, vcfheader_file, logger, intogen_driver_mut_info_tags, query_info_tags, "intogen_driver_mut")
   if pcgr_onco_xref is True:
      prepare_vcfanno_configuration(pcgr_db_directory, conf_fname, vcfheader_file, logger, pcgr_onco_xref_info_tags, query_info_tags, "pcgr_onco_xref")
   
   out_vcf_vcfanno_unsorted1 = output_vcf + '.tmp.unsorted.1'
   #out_vcf_vcfanno_unsorted2 = output_vcf + '.tmp.unsorted.2'
   #out_vcf_vcfanno_sorted = output_vcf + '.tmp.sorted.1'
   query_prefix = re.sub(r'\.vcf.gz$','',query_vcf)
   print_vcf_header(query_vcf, vcfheader_file, chromline_only = True)
   command1 = "vcfanno -p=" + str(num_processes) + " " + str(conf_fname) + " " + str(query_vcf) + " > " + str(out_vcf_vcfanno_unsorted1) + " 2> " + str(query_prefix) + '.vcfanno.log'
   os.system(command1)
   
   os.system('cat ' + str(vcfheader_file) + ' > ' + str(output_vcf))
   os.system('cat ' + str(out_vcf_vcfanno_unsorted1) + ' | grep -v \'^#\' >> ' + str(output_vcf))
   os.system('rm -f ' + str(output_vcf) + '.tmp*')
   os.system('bgzip -f ' + str(output_vcf))
   os.system('tabix -f -p vcf ' + str(output_vcf) + '.gz')
   return 0
   
def append_to_vcf_header(pcgr_db_directory, datasource, vcfheader_file):
   """
   Function that appends the VCF header information for a given 'datasource' (containing INFO tag formats/descriptions, and datasource version)
   """
   vcf_info_tags_file = str(pcgr_db_directory) + '/' + str(datasource) + '/' + str(datasource) + '.vcfanno.vcf_info_tags.txt'
   os.system('cat ' + str(vcf_info_tags_file) + ' >> ' + str(vcfheader_file))


def append_to_conf_file(datasource, datasource_info_tags, pcgr_db_directory, conf_fname):
   """
   Function that appends data to a vcfanno conf file ('conf_fname') according to user-defined ('datasource'). The datasource defines the set of tags that will be appended during annotation
   """
   fh = open(conf_fname,'a')
   if datasource != 'civic' and datasource != 'uniprot' and datasource != 'pcgr_onco_xref':
      fh.write('[[annotation]]\n')
      fh.write('file="' + str(pcgr_db_directory) + '/' + str(datasource) + '/' + str(datasource) + '.vcf.gz"\n')
      fields_string = 'fields = ["' + '","'.join(datasource_info_tags) + '"]'
      ops = ['concat'] * len(datasource_info_tags)
      ops_string = 'ops=["' + '","'.join(ops) + '"]'
      fh.write(fields_string + '\n')
      fh.write(ops_string + '\n\n')
   else:
      if datasource == 'uniprot' or datasource == 'pcgr_onco_xref':
         fh.write('[[annotation]]\n')
         fh.write('file="' + str(pcgr_db_directory) + '/' + str(datasource) + '/' + str(datasource) + '.bed.gz"\n')
         fh.write('columns=[4]\n')
         names_string = 'names=["' + '","'.join(datasource_info_tags) + '"]'
         fh.write(names_string +'\n')
         fh.write('ops=["concat"]\n\n')
      if datasource == 'civic':
         fh.write('[[annotation]]\n')
         fh.write('file="' + str(pcgr_db_directory) + '/' + str(datasource) + '/' + str(datasource) + '.bed.gz"\n')
         fh.write('columns=[4]\n')
         fh.write('names=["CIVIC_ID_2"]\n')
         fh.write('ops=["concat"]\n\n')

         fh.write('[[annotation]]\n')
         fh.write('file="' + str(pcgr_db_directory) + '/' + str(datasource) + '/' + str(datasource) + '.vcf.gz"\n')
         fh.write('fields = ["CIVIC_ID"]\n')
         fh.write('ops=["concat"]\n\n')
   fh.close()
   return

def get_vcf_info_tags(vcffile):
   vcf = VCF(vcffile)
   info_tags = {}
   for e in vcf.header_iter():
      header_element = e.info()
      if 'ID' in header_element.keys() and 'HeaderType' in header_element.keys():
         if header_element['HeaderType'] == 'INFO':
            info_tags[str(header_element['ID'])] = 1
   
   return info_tags


def print_vcf_header(query_vcf, vcfheader_file, chromline_only = False):
   if chromline_only == True:
      os.system('bgzip -dc ' + str(query_vcf) + ' | egrep \'^#\' | egrep \'^#CHROM\' >> ' + str(vcfheader_file))
   else:
      os.system('bgzip -dc ' + str(query_vcf) + ' | egrep \'^#\' | egrep -v \'^#CHROM\' > ' + str(vcfheader_file))

if __name__=="__main__": __main__()
