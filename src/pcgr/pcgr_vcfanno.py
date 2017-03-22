#!/usr/bin/env python

import argparse
import cyvcf
import vcfutils
import random
import pcgr
import os
import re
import sys

logger = pcgr.getlogger('pcgr-vcfanno')


def __main__():
   parser = argparse.ArgumentParser(description='Run brentp/vcfanno - annotate a VCF file against multiple VCF files in parallel', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
   parser.add_argument('query_vcf', help='Bgzipped input VCF file with query variants (SNVs/InDels)')
   parser.add_argument('out_vcf', help='Output VCF file with appended annotations from multiple VCF files')
   parser.add_argument('pcgr_dir', help='PCGR base directory')
   parser.add_argument('--num_processes', help="Number of processes vcfanno can use during annotation", default=4)
   parser.add_argument("--cosmic",action = "store_true", help="Annotate VCF with annotations from Catalogue of somatic mutations in cancer")
   parser.add_argument("--icgc",action = "store_true", help="Annotate VCF with annotations from all simple somatic mutations in ICGC projects")
   parser.add_argument("--exac",action = "store_true", help="Annotate VCF with allele frequencies from Exome Aggregation Consortium")
   parser.add_argument("--gnomad",action = "store_true", help="Annotate VCF with allele frequencies (exome) from genome Aggregation Database")
   parser.add_argument("--docm",action = "store_true", help="Annotate VCF with annotations from Database of Curated Mutations")
   parser.add_argument("--intogen_driver_mut",action = "store_true", help="Annotate VCF with predicted cancer driver mutations from IntoGen's Catalog of Driver Mutations")
   parser.add_argument("--clinvar",action = "store_true", help="Annotate VCF with annotations from ClinVar")
   parser.add_argument("--dbsnp",action = "store_true", help="Annotate VCF with annotations from database of short genetic variations")
   parser.add_argument("--dbnsfp",action = "store_true", help="Annotate VCF with annotations from database of non-synonymous functional predictions")
   parser.add_argument("--oneKG",action = "store_true", help="Annotate VCF with allele frequencies from the 1000 Genome Project")
   parser.add_argument("--civic",action = "store_true", help="Annotate VCF with annotations from the Clinical Interpretation of Variants in Cancer database")
   parser.add_argument("--cbmdb",action = "store_true", help="Annotate VCF with annotations from the Cancer bioMarkers database")
      
   args = parser.parse_args()

   query_info_tags = get_vcf_info_tags(args.query_vcf)
   vcfheader_file = args.out_vcf + '.tmp.' + str(random.randrange(0,10000000)) + '.header.txt'
   conf_fname = args.out_vcf + '.tmp.conf.toml'
   print_vcf_header(args.query_vcf, vcfheader_file, chromline_only = False)
   run_vcfanno(args.num_processes, args.query_vcf, query_info_tags, vcfheader_file, args.pcgr_dir, conf_fname, args.out_vcf, args.cosmic, args.icgc, args.exac, args.docm, args.intogen_driver_mut, args.clinvar, args.dbsnp, args.dbnsfp, args.oneKG, args.civic, args.cbmdb, args.gnomad)

def run_vcfanno(num_processes, query_vcf, query_info_tags, vcfheader_file, pcgr_directory, conf_fname, output_vcf, cosmic, icgc, exac, docm, intogen_driver_mut, clinvar, dbsnp, dbnsfp, oneKG, civic, cbmdb, gnomad):
   
   pcgr_db_directory = pcgr_directory + '/data'
   if cosmic is True:
      cosmic_tags = ["COSMIC_MUTATION_ID","COSMIC_COUNT_GW","COSMIC_DRUG_RESISTANCE","COSMIC_FATHMM_PRED","COSMIC_CANCER_TYPE_ALL","COSMIC_CANCER_TYPE_GW","COSMIC_CODON_COUNT_GW","COSMIC_CODON_FRAC_GW","COSMIC_SITE_HISTOLOGY","COSMIC_CONSEQUENCE","COSMIC_VARTYPE","COSMIC_SAMPLE_SOURCE"]
      for t in cosmic_tags:
         if query_info_tags.has_key(t):
            logger.warning("Query VCF has INFO tag " + str(t) + ' - this is also present in the COSMIC VCF annotation file. This tag will be overwritten if not renamed in the query VCF')
      append_to_conf_file("cosmic",pcgr_db_directory, conf_fname)
      append_to_vcf_header(pcgr_db_directory, "cosmic", vcfheader_file)
   if icgc is True:
      icgc_tags = ["ICGC_PROJECTS"]
      for t in icgc_tags:
         if query_info_tags.has_key(t):
            logger.warning("Query VCF has INFO tag " + str(t) + ' - this is also present in the ICGC VCF annotation file. This tag will be overwritten if not renamed in the query VCF')
      append_to_conf_file("icgc",pcgr_db_directory, conf_fname)
      append_to_vcf_header(pcgr_db_directory, "icgc", vcfheader_file)
   if civic is True:
      civic_tags = ["CIVIC_ID"]
      for t in civic_tags:
         if query_info_tags.has_key(t):
            logger.warning("Query VCF has INFO tag " + str(t) + ' - this is also present in the CIVIC VCF annotation file. This tag will be overwritten if not renamed in the query VCF')
      append_to_conf_file("civic",pcgr_db_directory, conf_fname)
      append_to_vcf_header(pcgr_db_directory, "civic", vcfheader_file)
   
   if cbmdb is True:
      cbmdb_tags = ["CBMDB_ID"]
      for t in cbmdb_tags:
         if query_info_tags.has_key(t):
            logger.warning("Query VCF has INFO tag " + str(t) + ' - this is also present in the CBMDB VCF annotation file. This tag will be overwritten if not renamed in the query VCF')
      append_to_conf_file("cbmdb",pcgr_db_directory, conf_fname)
      append_to_vcf_header(pcgr_db_directory, "cbmdb", vcfheader_file)

   if exac is True:
      exac_tags = ["AMR_AF_EXAC","AFR_AF_EXAC","NFE_AF_EXAC","FIN_AF_EXAC","OTH_AF_EXAC","GLOBAL_AF_EXAC","EAS_AF_EXAC","SAS_AF_EXAC"]
      for t in exac_tags:
         if query_info_tags.has_key(t):
            logger.warning("Query VCF has INFO tag " + str(t) + ' - this is also present in the ExAC VCF annotation file. This tag will be overwritten if not renamed in the query VCF')
      append_to_conf_file("exac",pcgr_db_directory, conf_fname)
      append_to_vcf_header(pcgr_db_directory, "exac", vcfheader_file)
   
   if gnomad is True:
      gnomad_tags = ["AMR_AF_GNOMAD","AFR_AF_GNOMAD","NFE_AF_GNOMAD","FIN_AF_GNOMAD","ASJ_AF_GNOMAD","OTH_AF_GNOMAD","GLOBAL_AF_GNOMAD","EAS_AF_GNOMAD","SAS_AF_GNOMAD"]
      for t in gnomad_tags:
         if query_info_tags.has_key(t):
            logger.warning("Query VCF has INFO tag " + str(t) + ' - this is also present in the GNOMAD VCF annotation file. This tag will be overwritten if not renamed in the query VCF')
      append_to_conf_file("gnomad",pcgr_db_directory, conf_fname)
      append_to_vcf_header(pcgr_db_directory, "gnomad", vcfheader_file)
   
   if docm is True:
      docm_tags = ["DOCM_DISEASE","DOCM_PMID"]
      for t in docm_tags:
         if query_info_tags.has_key(t):
            logger.warning("Query VCF has INFO tag " + str(t) + ' - this is also present in the DoCM VCF annotation file. This tag will be overwritten if not renamed in the query VCF')
      append_to_conf_file("docm",pcgr_db_directory, conf_fname)
      append_to_vcf_header(pcgr_db_directory, "docm", vcfheader_file)
   if intogen_driver_mut is True:
      intogen_driver_mut_tags = ["INTOGEN_DRIVER_MUT"]
      for t in intogen_driver_mut_tags:
         if query_info_tags.has_key(t):
            logger.warning("Query VCF has INFO tag " + str(t) + ' - this is also present in the Intogen Driver Mutations VCF annotation file. This tag will be overwritten if not renamed in the query VCF')
      append_to_conf_file("intogen_driver_mut",pcgr_db_directory, conf_fname)
      append_to_vcf_header(pcgr_db_directory, "intogen_driver_mut", vcfheader_file)

   if clinvar is True:
      clinvar_tags = ["CLINVAR_MSID","CLINVAR_PMIDS","CLINVAR_SIG","CLINVAR_VARIANT_ORIGIN"]
      for t in clinvar_tags:
         if query_info_tags.has_key(t):
            logger.warning("Query VCF has INFO tag " + str(t) + ' - this is also present in the ClinVar VCF annotation file. This tag will be overwritten if not renamed in the query VCF')
      append_to_conf_file("clinvar",pcgr_db_directory, conf_fname)
      append_to_vcf_header(pcgr_db_directory, "clinvar", vcfheader_file)
   if dbsnp is True:
      dbsnp_tags = ["GWAS_CATALOG_PMID","GWAS_CATALOG_TRAIT_URI","DBSNPRSID", "DBSNPBUILDID", "DBSNP_VALIDATION","DBSNP_MAPPINGSTATUS","DBSNP_SUBMISSIONS"]
      for t in dbsnp_tags:
         if query_info_tags.has_key(t):
            logger.warning("Query VCF has INFO tag " + str(t) + ' - this is also present in the dbSNP VCF annotation file. This tag will be overwritten if not renamed in the query VCF')
      append_to_conf_file("dbsnp",pcgr_db_directory, conf_fname)
      append_to_vcf_header(pcgr_db_directory, "dbsnp", vcfheader_file)
   if dbnsfp is True:
      dbnsfp_tags = ["DBNSFP"]
      for t in dbnsfp_tags:
         if query_info_tags.has_key(t):
            logger.warning("Query VCF has INFO tag " + str(t) + ' - this is also present in the dbNSFP VCF annotation file. This tag will be overwritten if not renamed in the query VCF')
      append_to_conf_file("dbnsfp",pcgr_db_directory, conf_fname)
      append_to_vcf_header(pcgr_db_directory, "dbnsfp", vcfheader_file)
   if oneKG is True:
      oneKG_tags = ["EAS_AF_1KG","EUR_AF_1KG","AMR_AF_1KG","AFR_AF_1KG","SAS_AF_1KG","GLOBAL_AF_1KG"]
      for t in oneKG_tags:
         if query_info_tags.has_key(t):
            logger.warning("Query VCF has INFO tag " + str(t) + ' - this is also present in the 1000Genomes Project VCF annotation file. This tag will be overwritten if not renamed in the query VCF')
      append_to_conf_file("oneKG",pcgr_db_directory, conf_fname)
      append_to_vcf_header(pcgr_db_directory, "oneKG", vcfheader_file)

   out_vcf_vcfanno_unsorted1 = output_vcf + '.tmp.unsorted.1'
   out_vcf_vcfanno_unsorted2 = output_vcf + '.tmp.unsorted.2'
   out_vcf_vcfanno_sorted = output_vcf + '.tmp.sorted.1'
   query_prefix = re.sub('\.vcf.gz$','',query_vcf)
   print_vcf_header(query_vcf, vcfheader_file, chromline_only = True)
   command1 = "vcfanno -p=" + str(num_processes) + " " + str(conf_fname) + " " + str(query_vcf) + " > " + str(out_vcf_vcfanno_unsorted1) + " 2> " + str(query_prefix) + '.vcfanno.log'
   os.system(command1)
   
   os.system('cat ' + str(vcfheader_file) + ' > ' + str(output_vcf))
   os.system('cat ' + str(out_vcf_vcfanno_unsorted1) + ' | grep -v \'^#\' >> ' + str(output_vcf))
   os.system('rm -f ' + str(output_vcf) + '.tmp*')
   os.system('bgzip -f ' + str(output_vcf))
   os.system('tabix -f -p vcf ' + str(output_vcf) + '.gz')
   return 0
   
def append_to_vcf_header(pcgr_db_directory, dbsource, vcfheader_file):
   
   vcf_info_tags_file = str(pcgr_db_directory) + '/' + str(dbsource) + '/' + str(dbsource) + '.vcfanno.vcf_info_tags.txt'
   os.system('cat ' + str(vcf_info_tags_file) + ' >> ' + str(vcfheader_file))


def append_to_conf_file(dbsource, pcgr_db_directory, conf_fname):
   fh = open(conf_fname,'a')
   if dbsource != 'civic':
      fh.write('[[annotation]]\n')
      fh.write('file="' + str(pcgr_db_directory) + '/' + str(dbsource) + '/' + str(dbsource) + '.vcf.gz"\n')
   if dbsource == 'cosmic':
      fh.write('fields = ["COSMIC_MUTATION_ID","COSMIC_DRUG_RESISTANCE","COSMIC_FATHMM_PRED","COSMIC_CANCER_TYPE_ALL","COSMIC_CANCER_TYPE_GW", "COSMIC_COUNT_GW", "COSMIC_CODON_COUNT_GW","COSMIC_CODON_FRAC_GW","COSMIC_SITE_HISTOLOGY","COSMIC_CONSEQUENCE","COSMIC_VARTYPE","COSMIC_SAMPLE_SOURCE"]\n')
      fh.write('ops=["concat","concat","concat","concat","concat","concat","concat","concat","concat","concat","concat","concat"]\n\n')
   if dbsource == 'dbsnp':
      fh.write('fields = ["GWAS_CATALOG_PMID","GWAS_CATALOG_TRAIT_URI","DBSNPRSID", "DBSNPBUILDID", "DBSNP_VALIDATION","DBSNP_MAPPINGSTATUS","DBSNP_SUBMISSIONS"]\n')
      fh.write('ops=["concat","concat", "concat", "concat", "concat", "concat","concat"]\n\n')
   if dbsource == 'dbnsfp':
      
      fh.write('fields = ["DBNSFP"]\n')
      fh.write('ops=["concat"]\n\n')
      
   if dbsource == 'exac':
      fh.write('fields = ["AMR_AF_EXAC","AFR_AF_EXAC","NFE_AF_EXAC","FIN_AF_EXAC","OTH_AF_EXAC","GLOBAL_AF_EXAC","EAS_AF_EXAC","SAS_AF_EXAC"]\n')
      fh.write('ops=["concat","concat","concat","concat","concat","concat","concat","concat"]\n\n')
   
   if dbsource == 'gnomad':
      fh.write('fields = ["AMR_AF_GNOMAD","AFR_AF_GNOMAD","NFE_AF_GNOMAD","FIN_AF_GNOMAD","ASJ_AF_GNOMAD","OTH_AF_GNOMAD","GLOBAL_AF_GNOMAD","EAS_AF_GNOMAD","SAS_AF_GNOMAD"]\n')
      fh.write('ops=["concat","concat","concat","concat","concat","concat","concat","concat","concat"]\n\n')
   
   if dbsource == 'intogen_driver_mut':
      fh.write('fields = ["INTOGEN_DRIVER_MUT"]\n')
      fh.write('ops=["concat"]\n\n')
      
   if dbsource == 'cbmdb':
      fh.write('fields = ["CBMDB_ID"]\n')
      fh.write('ops=["concat"]\n\n')
      
   if dbsource == 'oneKG':
      fh.write('fields = ["EAS_AF_1KG","EUR_AF_1KG","AMR_AF_1KG","AFR_AF_1KG","SAS_AF_1KG","GLOBAL_AF_1KG"]\n')
      fh.write('ops=["concat","concat","concat","concat","concat","concat"]\n\n')

   if dbsource == 'docm':
      fh.write('fields = ["DOCM_DISEASE","DOCM_PMID"]\n')
      fh.write('ops=["concat","concat"]\n\n')
   if dbsource == 'civic':
      fh.write('[[annotation]]\n')
      fh.write('file="' + str(pcgr_db_directory) + '/' + str(dbsource) + '/' + str(dbsource) + '.bed.gz"\n')
      fh.write('columns=[4]\n')
      fh.write('names=["CIVIC_ID_2"]\n')
      fh.write('ops=["concat"]\n\n')

      fh.write('[[annotation]]\n')
      fh.write('file="' + str(pcgr_db_directory) + '/' + str(dbsource) + '/' + str(dbsource) + '.vcf.gz"\n')
      fh.write('fields = ["CIVIC_ID"]\n')
      fh.write('ops=["concat"]\n\n')
      
   if dbsource == 'clinvar':
      fh.write('fields = ["CLINVAR_MSID","CLINVAR_PMIDS","CLINVAR_SIG","CLINVAR_VARIANT_ORIGIN"]\n')
      fh.write('ops=["concat","concat","concat","concat"]\n\n')
   if dbsource == 'icgc':
      fh.write('fields = ["ICGC_PROJECTS"]\n')
      fh.write('ops=["concat"]\n\n')
   fh.close()
   return

def get_vcf_info_tags(vcffile):
   vcf_reader = cyvcf.Reader(open(vcffile, 'r'))
   info_tags = {}
   for info_tag in sorted(vcf_reader.infos.keys()):
      info_tags[str(info_tag)] = 1
   
   return info_tags


def print_vcf_header(query_vcf, vcfheader_file, chromline_only = False):
   if chromline_only == True:
      os.system('bgzip -dc ' + str(query_vcf) + ' | egrep \'^#\' | egrep \'^#CHROM\' >> ' + str(vcfheader_file))
   else:
      os.system('bgzip -dc ' + str(query_vcf) + ' | egrep \'^#\' | egrep -v \'^#CHROM\' > ' + str(vcfheader_file))

if __name__=="__main__": __main__()
