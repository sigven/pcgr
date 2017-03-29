#!/usr/bin/env python

import csv
import re
import argparse
import os
import subprocess
import logging
import sys
import ntpath

def __main__():
   
   parser = argparse.ArgumentParser(description='Personal Cancer Genome Reporter (PCGR) workflow for clinical interpretation of somatic nucleotide variants and copy number aberration segments',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
   parser.add_argument('--input_vcf', dest="input_vcf",help='VCF input file with somatic query variants (SNVs/InDels)')
   parser.add_argument('--input_cna_segments', dest="input_cna_segments",help='Somatic copy number alteration segments (tab-separated values)')
   parser.add_argument('--logR_threshold_amplification', dest = "logR_threshold_amplification",default=0.8, type=float, help='Log(2) ratio treshold for calling copy number amplifications in HTML report')
   parser.add_argument('--logR_threshold_homozygous_deletion', dest = "logR_threshold_homozygous_deletion", default=-0.8, type=float, help='Log(2) ratio treshold for calling homozygous deletions in HTML report')
   parser.add_argument('--num_vcfanno_processes', dest = "num_vcfanno_processes", default=4, type=int, help='Number of processes used during vcfanno annotation')
   parser.add_argument('--num_vep_forks', dest = "num_vep_forks", default=4, type=int, help='Number of forks (--forks option in VEP) used during VEP annotation')
   parser.add_argument('--force_overwrite', action = "store_true", help='By default, the script will fail with an error if any output file already exists. You can force the overwrite of existing result files by using this flag')
   parser.add_argument('pcgr_directory',help='Full path of PCGR base directory')
   parser.add_argument('working_directory',help="Full path of working directory - directory with input/output files")
   parser.add_argument('sample_id',help="Tumor sample/cancer genome identifier - prefix for output files")
   
   args = parser.parse_args()
   
   overwrite = 0
   if args.force_overwrite is True:
      overwrite = 1
   
   fname_cna = None
   fname_vcf = None
   
   logger = getlogger('pcgr-prepare')
   
   pcgr_db_directory = str(args.pcgr_directory) + "/data"
   if not os.path.exists(pcgr_db_directory):
      logger.error("") 
      logger.error("PCGR data directory (" + str(pcgr_db_directory) + ") does not exist")
      logger.error("")
      return
   
   if(args.input_vcf is None and args.input_cna_segments is None):
      logger.error("")
      logger.error("Please specifiy either a VCF input file (--input_vcf) or a copy number segment file (--input_cna_segments)")
      logger.error("")
      return
   
   if not args.input_vcf is None:
      fname_vcf_full = str(args.working_directory) + '/' + path_leaf(str(args.input_vcf))
      fname_vcf = path_leaf(str(args.input_vcf))
      if not (fname_vcf.endswith('.vcf') or fname_vcf.endswith('.vcf.gz')):
         logger.error('')
         logger.error('Input VCF does not have the correct file extension (.vcf or .vcf.gz)')
         logger.error('')
         return
      output_vcf = str(args.working_directory) + '/' + str(args.sample_id) + '.pcgr.vcf.gz'
      if not os.path.exists(fname_vcf_full) or os.path.getsize(fname_vcf_full) == 0:
         logger.error("") 
         logger.error("VCF input file (" + fname_vcf_full +") is empty or does not exist")
         logger.error("")
         return
      if os.path.exists(output_vcf) and overwrite == 0:
         logger.error('')
         logger.error("Output files (e.g. " + str(output_vcf) + ") already exist - please specify different sample_id or add option --force_overwrite")
         logger.error('')
         return
         
   if not args.input_cna_segments is None:
      fname_cna_full = str(args.working_directory) + '/' + path_leaf(str(args.input_cna_segments))
      fname_cna = path_leaf(str(args.input_cna_segments))
      output_cna = str(args.working_directory) + '/' + str(args.sample_id) + '.pcgr.cna_segments.tsv.gz'
      if not os.path.exists(fname_cna_full) or os.path.getsize(fname_cna_full) == 0:
         logger.error("") 
         logger.error("Input copy number segment file (" + fname_cna_full + ") is empty or does not exist")
         logger.error("") 
         return
      if os.path.exists(output_cna) and overwrite == 0:
         logger.error('')
         logger.error("Output files (e.g. " + str(output_cna) + ") already exist - please specify different sample_id or add option --force_overwrite")
         logger.error('')
         return

   run_pcgr(fname_vcf, fname_cna, args.logR_threshold_amplification, args.logR_threshold_homozygous_deletion, args.num_vcfanno_processes, args.num_vep_forks, args.pcgr_directory, args.working_directory, args.sample_id, overwrite)

def path_leaf(path):
   head, tail = ntpath.split(path)
   return tail or ntpath.basename(head)


def check_subprocess(command):
   try:
      output = subprocess.check_output(str(command), stderr=subprocess.STDOUT, shell=True)
      if len(output) > 0:
         print str(output).rstrip()
   except subprocess.CalledProcessError as e:
      print e.output
      exit(0)

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

def run_pcgr(input_vcf, input_cna_segments, logR_threshold_amplification, logR_threshold_homozygous_deletion, num_vcfanno_processes, num_vep_forks, pcgr_directory, working_directory, sample_id, overwrite):
   
   ## set basic Docker run commands
   output_vcf = 'None'
   vepdb_directory = str(pcgr_directory) + '/data/.vep'
   docker_command_run1 = "docker run -t -v=" + str(vepdb_directory) + ":/usr/local/share/vep/data -v=" + str(working_directory) + ":/workdir -w=/workdir sigven/pcgr:latest sh -c \""
   docker_command_run2 = "docker run -t -v=" + str(pcgr_directory) + ":/data -v=" + str(working_directory) + ":/workdir -w=/workdir sigven/pcgr:latest sh -c \""
   docker_command_run3 = "docker run -t -v=" + str(working_directory) + ":/workdir -w=/workdir sigven/pcgr:latest sh -c \""
   
   ## verify VCF and CNA segment file
   logger = getlogger('pcgr-check-input')
   logger.info("STEP 0: Validate input data")
   vcf_validate_command = str(docker_command_run3) + "pcgr_check_input.py " + str(input_vcf) + " " + str(input_cna_segments) + "\""
   check_subprocess(vcf_validate_command)
   logger.info('Finished')

   if not input_vcf is None:
      
      ## Run VEP + vcfanno + summarise
      output_vcf = str(sample_id) + '.pcgr.vcf.gz'
      input_vcf_pcgr_ready = re.sub(r'(\.vcf$|\.vcf\.gz$)','.pcgr_ready.vcf.gz',input_vcf)
      vep_vcf = re.sub(r'(\.vcf$|\.vcf\.gz$)','.pcgr_vep.vcf',input_vcf_pcgr_ready)
      vep_vcfanno_vcf = re.sub(r'(\.vcf$|\.vcf\.gz$)','.pcgr_vep.vcfanno.vcf',input_vcf_pcgr_ready)
      vep_vcf_gz = vep_vcf + '.gz'
      vep_tmp_vcf = vep_vcf + '.tmp'
      vep_vcfanno_annotated_vcf = re.sub(r'\.vcfanno','.vcfanno.annotated',vep_vcfanno_vcf) + '.gz'
      vep_main_command = str(docker_command_run1) + "variant_effect_predictor.pl --input_file " + str(input_vcf_pcgr_ready) + " --output_file " + str(vep_tmp_vcf) + " --vcf --check_ref --flag_pick_allele --force_overwrite --species homo_sapiens --assembly GRCh37 --offline --fork " + str(num_vep_forks) + " --no_progress --variant_class --regulatory --domains --shift_hgvs 1 --hgvs --symbol --protein --ccds --uniprot --appris --biotype --canonical --gencode_basic --cache --numbers --check_alleles --total_length --allele_number --no_escape --xref_refseq --dir /usr/local/share/vep/data --fasta /usr/local/share/vep/data/homo_sapiens/85_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz\""
      vep_sed_command =  str(docker_command_run1) + "sed -r 's/\(p\.=\)|polyphen=2.2.2 |ClinVar=201507 |dbSNP=144 |ESP=20141103|sift=sift5.2.2 |COSMIC=71 | HGMD-PUBLIC=20152//g' " + str(vep_tmp_vcf) + " > " + str(vep_vcf) + "\""
      vep_bgzip_command = str(docker_command_run1) + "bgzip -f " + str(vep_vcf) + "\""
      vep_tabix_command = str(docker_command_run1) + "tabix -f -p vcf " + str(vep_vcf) + ".gz" + "\""
      logger = getlogger('pcgr-vep')
   
      print
      logger.info("STEP 1: Basic variant annotation with Variant Effect Predictor (v85, GENCODE, GRCh37)")
      check_subprocess(vep_main_command)
      check_subprocess(vep_sed_command)
      check_subprocess(vep_bgzip_command)
      check_subprocess(vep_tabix_command)
   
      print
      logger = getlogger('pcgr-vcfanno')
      logger.info("STEP 2: Annotation for precision oncology with pcgr-vcfanno (ClinVar, dbSNP, dbNSFP, 1000Genomes Project, ExAC, gnomAD, CiVIC, CBMDB, DoCM, COSMIC, Intogen_driver_mutations, ICGC)")
      pcgr_vcfanno_command = str(docker_command_run2) + "pcgr_vcfanno.py --num_processes " + str(num_vcfanno_processes) + " --dbsnp --dbnsfp --oneKG --docm --clinvar --exac --gnomad --icgc --civic --cbmdb --intogen_driver_mut --cosmic " + str(vep_vcf_gz) + " " + str(vep_vcfanno_vcf) + " /data\""
      check_subprocess(pcgr_vcfanno_command)
      logger.info("Finished")
      
      print
      logger = getlogger("pcgr-summarise")
      pcgr_summarise_command = str(docker_command_run2) + "pcgr_summarise.py " + str(vep_vcfanno_vcf) + ".gz /data\""
      logger.info("STEP 3: Cancer gene annotations with pcgr-summarise")
      check_subprocess(pcgr_summarise_command)
      logger.info("Finished")
      
      create_output_vcf_command1 = str(docker_command_run3) + 'mv ' + str(vep_vcfanno_annotated_vcf) + ' ' + str(output_vcf) + "\""
      create_output_vcf_command2 = str(docker_command_run3) + 'mv ' + str(vep_vcfanno_annotated_vcf) + '.tbi ' + str(output_vcf) + '.tbi' + "\""
      clean_command = str(docker_command_run3) + 'rm -f ' + str(vep_vcf) + '* ' +  str(input_vcf_pcgr_ready) + "* " + str(input_vcf) + "*_validator_output" + "\""
      check_subprocess(create_output_vcf_command1)
      check_subprocess(create_output_vcf_command2)
      check_subprocess(clean_command)
  
   print
   
   ## Generation of HTML reports for VEP/vcfanno-annotated VCF and copy number segment file
   logger = getlogger('pcgr-writer')
   logger.info("STEP 4: Generation of output files")
   pcgr_report_command = str(docker_command_run2) + "/pcgr.R /workdir " + str(output_vcf) + " " + str(input_cna_segments) + " "  + str(sample_id)  + " " + str(logR_threshold_amplification) + " " + str(logR_threshold_homozygous_deletion) + "\""
   check_subprocess(pcgr_report_command)
   logger.info("Finished")
   
   
   
if __name__=="__main__": __main__()

