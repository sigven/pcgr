#!/usr/bin/env python

import csv
import re
import argparse
import os
import subprocess
import logging
import sys

def __main__():
   
   parser = argparse.ArgumentParser(description='Personal Cancer Genome Reporter (PCGR) workflow for clinical interpretation of somatic nucleotide variants and copy number segments',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
   parser.add_argument('--input_vcf', dest="input_vcf",help='VCF input file with somatic query variants (SNVs/InDels)')
   parser.add_argument('--input_cna_segments', dest="input_cna_segments",help='Somatic copy number alteration segments (tab-separated values)')
   parser.add_argument('--logR_threshold_amplification', dest = "logR_threshold_amplification",default=0.8, help='Log(2) ratio treshold for copy number amplification')
   parser.add_argument('--logR_threshold_homozygous_deletion', dest = "logR_threshold_homozygous_deletion", default=-0.8, help='Log(2) ratio treshold for homozygous deletion')
   parser.add_argument('--num_vcfanno_processes', dest = "num_vcfanno_processes", default=4, help='Number of processes used during vcfanno annotation')
   parser.add_argument('--num_vep_forks', dest = "num_vep_forks", default=4, help='Number of forks (--forks) used during VEP annotation')
   parser.add_argument('pcgr_directory',help='PCGR base directory')
   parser.add_argument('working_directory',help="Working directory")
   parser.add_argument('sample_id',help="Tumor sample/cancer genome identifier - prefix for output files")
   
   args = parser.parse_args()
   
   if(args.input_vcf is None and args.input_cna_segments is None):
      print 
      print "ERROR: Please specifiy either a VCF input file (--input_vcf) or a copy number segment file (input_cna_segments)"
      print
      return
   
   fname_vcf = str(args.working_directory) + '/' + str(args.input_vcf)
   fname_cna = str(args.working_directory) + '/' + str(args.input_cna_segments)
   if not os.path.exists(fname_vcf) or os.path.getsize(fname_vcf) == 0:
      print 
      print "ERROR: VCF input file (", args.input_vcf,") is empty or does not exist"
      print
      return
   
   if not os.path.exists(fname_cna) or os.path.getsize(fname_cna) == 0:
      print 
      print "ERROR: Input copy number segment file (", args.input_cna_segments,") is empty or does not exist"
      print
      return

   run_pcgr(args.input_vcf, args.input_cna_segments, args.logR_threshold_amplification, args.logR_threshold_homozygous_deletion, args.num_vcfanno_processes, args.num_vep_forks, args.pcgr_directory, args.working_directory, args.sample_id)

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

def run_pcgr(input_vcf, input_cna_segments, logR_threshold_amplification, logR_threshold_homozygous_deletion, num_vcfanno_processes, num_vep_forks, pcgr_directory, working_directory, sample_id):
   
   ## MOVE TO validate_arguments
   vepdb_directory = pcgr_directory + "/data/.vep"
   if not os.path.exists(vepdb_directory):
   ## logger - error - PCGR data directory not installed or EMPTY (TODO)
     print str(vepdb_directory) + ' does not exist' 
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
   vcf_validate_command = str(docker_command_run3) + "pcgr_check_input.py " + str(input_vcf) + " " + str(input_cna_segments) + "\""
   subprocess.check_call(str(vcf_validate_command), shell=True)
   logger.info('Finished')

   
   vep_main_command = str(docker_command_run1) + "variant_effect_predictor.pl --input_file " + str(input_vcf) + " --output_file " + str(vep_tmp_vcf) + " --vcf --check_ref --flag_pick_allele --force_overwrite --species homo_sapiens --assembly GRCh37 --offline --fork " + str(num_vep_forks) + " --no_progress --variant_class --regulatory --domains --shift_hgvs 1 --hgvs --symbol --protein --ccds --uniprot --appris --biotype --canonical --gencode_basic --cache --numbers --check_alleles --total_length --allele_number --no_escape --xref_refseq --dir /usr/local/share/vep/data --fasta /usr/local/share/vep/data/homo_sapiens/85_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz\""
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
   logger.info("STEP 2: Annotation for precision oncology with pcgr-vcfanno (ClinVar, dbSNP, dbNSFP, 1000Genomes Project, ExAC, gnomAD, CiVIC, CBMDB, DoCM, COSMIC, Intogen_driver_mutations, ICGC)")
   pcgr_vcfanno_command = str(docker_command_run2) + "pcgr_vcfanno.py --num_processes " + str(num_vcfanno_processes) + " --dbsnp --dbnsfp --oneKG --docm --clinvar --exac --gnomad --icgc --civic --cbmdb --intogen_driver_mut --cosmic " + str(vep_vcf_gz) + " " + str(vep_vcfanno_vcf) + " /data\""
   subprocess.check_call(str(pcgr_vcfanno_command), shell=True)
   logger.info("Finished")
   
   print
   logger = getlogger("pcgr-summarise")
   pcgr_summarise_command = str(docker_command_run2) + "pcgr_summarise.py " + str(vep_vcfanno_vcf) + ".gz /data\""
   logger.info("STEP 3: Cancer gene annotations with pcgr-summarise")
   subprocess.check_call(str(pcgr_summarise_command), shell=True)
   logger.info("Finished")
   
   create_output_vcf_command1 = str(docker_command_run3) + 'mv ' + str(vep_vcfanno_annotated_vcf) + ' ' + str(output_vcf) + "\""
   create_output_vcf_command2 = str(docker_command_run3) + 'mv ' + str(vep_vcfanno_annotated_vcf) + '.tbi ' + str(output_vcf) + '.tbi' + "\""
   subprocess.check_call(create_output_vcf_command1, shell=True)
   subprocess.check_call(create_output_vcf_command2, shell=True)
  
   print
   logger = getlogger('pcgr-writer')
   logger.info("STEP 4: Generation of output files")
   pcgr_report_command = str(docker_command_run2) + "/pcgr.R /workdir " + str(output_vcf) + " " + str(input_cna_segments) + " "  + str(sample_id)  + " " + str(logR_threshold_amplification) + " " + str(logR_threshold_homozygous_deletion) + "\""
   subprocess.check_call(str(pcgr_report_command), shell=True)
   logger.info("Finished")
   
   clean_command = str(docker_command_run3) + 'rm -f ' + str(vep_vcf) + '*' + "\""
   subprocess.check_call(str(clean_command), shell=True)
   
if __name__=="__main__": __main__()

