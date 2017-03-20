#!/usr/bin/env python

import csv
import re
import argparse
import os
import subprocess
import logging
import sys

def __main__():
   
   parser = argparse.ArgumentParser(description='Personal Cancer Genome Reporter (PCGR) workflow for clinical interpretation of somatic nucleotide variants and copy number segments')
   parser.add_argument('--input_vcf', dest="input_vcf",help='VCF input file with somatic query variants (SNVs/InDels)')
   parser.add_argument('--input_cnv_segments', dest="input_cnv_segments",help='Somatic copy number query segments (tab-separated values)')
   parser.add_argument('--no_html_report', action = "store_true",help='Skip generation of HTML reports, output will consist of annotated VCF/TSV files')
   parser.add_argument('--logR_threshold_amplification', dest = "logR_threshold_amplification",help='Log2 ratio treshold for copy number amplification')
   parser.add_argument('--logR_threshold_loss', dest = "logR_treshold_loss",help='Log2 ratio treshold for homozygous deletion')
   parser.add_argument('pcgr_directory',help='PCGR base directory')
   parser.add_argument('working_directory',help="Working directory")
   parser.add_argument('sample_id',help="Tumor sample/cancer genome identifier - prefix for output files")
   
   args = parser.parse_args()
   
   ##TODO: implement method validate_arguments()
   ## 1. check that input files exist and are > 0 in size
   ## 2. verify that VCF is properly formatted (vcf-validate?), IMPORTANT: no multi-allelic calls
   ## 3. Chromosomes should not have 'chr' prefix

   ## 3. make custom check of CNV segment file
   ## 4. ensure that either input_vcf or input_cnv is not NULL

   run_pcgr(args.input_vcf, args.input_cnv_segments, args.pcgr_directory, args.working_directory, args.no_html_report, args.sample_id)

def getlogger(logger_name):
   logger = logging.getLogger(logger_name)
   logger.setLevel(logging.DEBUG)

   # create console handler and set level to debug
   ch = logging.StreamHandler(sys.stdout)
   ch.setLevel(logging.DEBUG)

   # add ch to logger
   logger.addHandler(ch)
   
   # create formatter
   formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s", "20%y-%m-%d %H:%M:%S")
   
   #add formatter to ch
   ch.setFormatter(formatter)
   
   return logger

def run_pcgr(input_vcf, input_cnv_segments, pcgr_directory, working_directory, no_html_report, sample_id):
   
   
   
   ## MOVE TO validate_arguments
   vepdb_directory = pcgr_directory + "/data/.vep"
   if not os.path.exists(vepdb_directory):
   ## logger - error - PCGR data directory not installed or EMPTY (TODO)
      return 
   
   output_vcf = str(sample_id) + '.pcgr.vcf.gz'
   vep_vcf = re.sub(r'(\.vcf$|\.vcf\.gz$)','.pcgr_vep.vcf',input_vcf)
   vep_vcfanno_vcf = re.sub(r'(\.vcf$|\.vcf\.gz$)','.pcgr_vep.vcfanno.vcf',input_vcf)
   vep_vcf_gz = vep_vcf + '.gz'
   vep_tmp_vcf = vep_vcf + '.tmp'
   vep_vcfanno_annotated_vcf = re.sub(r'\.vcfanno','.vcfanno.annotated',vep_vcfanno_vcf) + '.gz'

   docker_command_run1 = "docker run -t -v=" + str(vepdb_directory) + ":/usr/local/share/vep/data -v=" + str(working_directory) + ":/workdir -w=/workdir sigven/pcgr:latest sh -c \""
   docker_command_run2 = "docker run -t -v=" + str(pcgr_directory) + ":/data -v=" + str(working_directory) + ":/workdir -w=/workdir sigven/pcgr:latest sh -c \""
   docker_command_run3 = "docker run -t -v=" + str(working_directory) + ":/workdir -w=/workdir sigven/pcgr:latest sh -c \""
   
   print
   logger = getlogger('pcgr-check-input')
   logger.info("STEP 0: Validate input data")
   vcf_validate_command = str(docker_command_run3) + "pcgr_check_input.py " + str(input_vcf) + " " + str(input_cnv_segments) + "\""
   subprocess.check_call(str(vcf_validate_command), shell=True)
   logger.info('Finished')

   
   vep_main_command = str(docker_command_run1) + "variant_effect_predictor.pl --input_file " + str(input_vcf) + " --output_file " + str(vep_tmp_vcf) + " --vcf --check_ref --flag_pick_allele --force_overwrite --species homo_sapiens --assembly GRCh37 --offline --fork 4 --no_progress --variant_class --regulatory --domains --shift_hgvs 1 --hgvs --symbol --protein --ccds --uniprot --appris --biotype --canonical --gencode_basic --cache --numbers --check_alleles --total_length --allele_number --no_escape --xref_refseq --dir /usr/local/share/vep/data --fasta /usr/local/share/vep/data/homo_sapiens/85_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz\""
   vep_sed_command =  str(docker_command_run1) + "sed -r 's/\(p\.=\)|polyphen=2.2.2 |ClinVar=201507 |dbSNP=144 |ESP=20141103|sift=sift5.2.2 |COSMIC=71 | HGMD-PUBLIC=20152//g' " + str(vep_tmp_vcf) + " > " + str(vep_vcf) + "\""
   vep_bgzip_command = str(docker_command_run1) + "bgzip -f " + str(vep_vcf) + "\""
   vep_tabix_command = str(docker_command_run1) + "tabix -f -p vcf " + str(vep_vcf) + ".gz" + "\""
   logger = getlogger('pcgr-vep')


   print
   logger.info("STEP 1: Basic variant annotation with Variant Effect Predictor (v85, GENCODE, GRCh37)")
   subprocess.check_call(str(vep_main_command), shell=True)
   subprocess.check_call(str(vep_sed_command), shell=True)
   subprocess.check_call(str(vep_bgzip_command), shell=True)
   subprocess.check_call(str(vep_tabix_command), shell=True)

   print
   logger = getlogger('pcgr-vcfanno')
   logger.info("STEP 2: Annotation for precision oncology with pcgr-vcfanno (ClinVar,dbSNP,dbNSFP,1000Genomes Project,ExAC,CiVIC,CBMDB,DoCM,COSMIC,Intogen_driver_mutations,ICGC)")
   pcgr_vcfanno_command = str(docker_command_run2) + "pcgr_vcfanno.py --dbsnp --dbnsfp --oneKG --docm --clinvar --exac --icgc --civic --cbmdb --intogen_driver_mut --cosmic " + str(vep_vcf_gz) + " " + str(vep_vcfanno_vcf) + " /data\""
   subprocess.check_call(str(pcgr_vcfanno_command), shell=True)
   logger.info("Finished")
   
   print
   logger = getlogger("pcgr-summarise")
   pcgr_summarise_command = str(docker_command_run2) + "pcgr_summarise.py " + str(vep_vcfanno_vcf) + ".gz /data\""
   logger.info("STEP 3: Cancer gene annotations with pcgr-summarise")
   subprocess.check_call(str(pcgr_summarise_command), shell=True)
   logger.info("Finished")
   
   #print str(vep_vcfanno_annotated_vcf)
   create_output_vcf_command1 = str(docker_command_run3) + 'mv ' + str(vep_vcfanno_annotated_vcf) + ' ' + str(output_vcf) + "\""
   create_output_vcf_command2 = str(docker_command_run3) + 'mv ' + str(vep_vcfanno_annotated_vcf) + '.tbi ' + str(output_vcf) + '.tbi' + "\""
   subprocess.check_call(create_output_vcf_command1, shell=True)
   subprocess.check_call(create_output_vcf_command2, shell=True)
  
   print
   logger = getlogger('pcgr-writer')
   logger.info("STEP 4: Generation of output files")
      
   #pcgr_report_command = str(docker_command_run3) + "/pcgr.R /workdir " + str(output_vcf) + " " + str(input_cnv_segments) + " "  + str(sample_id)  + " " + str(no_html_report) + "\""
   pcgr_report_command = str(docker_command_run2) + "/pcgr_v2.R /workdir " + str(output_vcf) + " " + str(input_cnv_segments) + " "  + str(sample_id)  + " " + str(no_html_report) + "\""

   subprocess.check_call(str(pcgr_report_command), shell=True)
   logger.info("Finished")
   
   clean_command = str(docker_command_run3) + 'rm -f ' + str(vep_vcf) + '*' + "\""
   subprocess.check_call(str(clean_command), shell=True)
   
if __name__=="__main__": __main__()

