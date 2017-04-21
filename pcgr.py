#!/usr/bin/env python

import csv
import re
import argparse
import os
import subprocess
import logging
import sys

version = '0.3.4'

def __main__():
   
   parser = argparse.ArgumentParser(description='Personal Cancer Genome Reporter (PCGR) workflow for clinical interpretation of somatic nucleotide variants and copy number aberration segments',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
   parser.add_argument('--input_vcf', dest="input_vcf",help='VCF input file with somatic query variants (SNVs/InDels)')
   parser.add_argument('--input_cna_segments', dest="input_cna_segments",help='Somatic copy number alteration segments (tab-separated values)')
   parser.add_argument('--logR_threshold_amplification', dest = "logR_threshold_amplification",default=0.8, type=float, help='Log(2) ratio treshold for calling copy number amplifications in HTML report')
   parser.add_argument('--logR_threshold_homozygous_deletion', dest = "logR_threshold_homozygous_deletion", default=-0.8, type=float, help='Log(2) ratio treshold for calling homozygous deletions in HTML report')
   parser.add_argument('--num_vcfanno_processes', dest = "num_vcfanno_processes", default=4, type=int, help='Number of processes used during vcfanno annotation')
   parser.add_argument('--num_vep_forks', dest = "num_vep_forks", default=4, type=int, help='Number of forks (--forks option in VEP) used during VEP annotation')
   parser.add_argument('--force_overwrite', action = "store_true", help='By default, the script will fail with an error if any output file already exists. You can force the overwrite of existing result files by using this flag')
   parser.add_argument('--version', action='version', version='%(prog)s ' + str(version))
   parser.add_argument('pcgr_dir',help='PCGR base directory with accompanying data directory, e.g. ~/pcgr-0.3.1')
   parser.add_argument('output_dir',help='Output directory')
   parser.add_argument('sample_id',help="Tumor sample/cancer genome identifier - prefix for output files")
   
   docker_image_version = 'sigven/pcgr:' + str(version)
   args = parser.parse_args()
   
   overwrite = 0
   if args.force_overwrite is True:
      overwrite = 1
   
   ## check that script and Docker image version correspond
   check_docker_command = 'docker images -q ' + str(docker_image_version)
   output = subprocess.check_output(str(check_docker_command), stderr=subprocess.STDOUT, shell=True)
   if(len(output) == 0):
      logger.error('')
      logger.error('Docker image ' + str(docker_image_version) + ' does not exist, pull image from Dockerhub (docker pull ' + str(docker_image_version) + ')')
      logger.error('')
      exit(0)
     
   logger = getlogger('pcgr-check-files')
   host_directories = check_input_files(args.input_vcf, args.input_cna_segments, args.pcgr_dir, args.output_dir, overwrite, args.sample_id, logger)

   run_pcgr(host_directories, docker_image_version, args.logR_threshold_amplification, args.logR_threshold_homozygous_deletion, args.num_vcfanno_processes, args.num_vep_forks, args.sample_id)


def check_input_files(input_vcf, input_cna_segments, base_pcgr_dir, output_dir, overwrite, sample_id, logger):
   """
   Function that checks the input files and directories provided by the user and checks for their existence
   """
   
   input_vcf_dir = "NA"
   input_cna_dir = "NA"
   db_dir = "NA"
   base_dir = "NA"
   output_dir_full = "NA"
   input_vcf_basename = "NA"
   input_cna_basename = "NA"
   
   ## check that either input vcf or cna segments exist
   if(input_vcf is None and input_cna_segments is None):
      logger.error("")
      logger.error("Please specifiy either a VCF input file (--input_vcf) or a copy number segment file (--input_cna_segments)")
      logger.error("")
      exit(0)
   
   ## check the existence of given output folder
   output_dir_full = os.path.abspath(output_dir)
   if not os.path.isdir(output_dir_full):
      logger.error("")
      logger.error("Output directory (" + str(output_dir_full) + ") does not exist")
      logger.error("")
      exit(0)
   
   ## check if input vcf exist
   if not input_vcf is None:
      if not os.path.exists(os.path.abspath(input_vcf)):
         logger.error("")
         logger.error("Input file (" + str(input_vcf) + ") does not exist")
         logger.error("")
         exit(0)
      if not (os.path.abspath(input_vcf).endswith('.vcf') or os.path.abspath(input_vcf).endswith('.vcf.gz')):
         logger.error('')
         logger.error("VCF input file (" + os.path.abspath(input_vcf) + ") does not have the correct file extension (.vcf or .vcf.gz)")
         logger.error('')
         exit(0)

      input_vcf_basename = os.path.basename(str(input_vcf))
      input_vcf_dir = os.path.dirname(os.path.abspath(input_vcf))

      ## if output vcf exist and overwrite not set
      output_vcf = os.path.join(str(output_dir_full),str(sample_id)) + '.pcgr.vcf.gz'
      if os.path.exists(output_vcf) and overwrite == 0:
         logger.error('')
         logger.error("Output files (e.g. " + str(output_vcf) + ") already exist - please specify different sample_id or add option --force_overwrite")
         logger.error('')
         exit(0)
   
   ## check if input cna segments exist
   if not input_cna_segments is None:
      if not os.path.exists(os.path.abspath(input_cna_segments)):
         logger.error("")
         logger.error("Input file (" + str(input_cna_segments) + ") does not exist")
         logger.error("")
         exit(0)
      if not (os.path.abspath(input_cna_segments).endswith('.tsv') or os.path.abspath(input_cna_segments).endswith('.txt')):
         logger.error('')
         logger.error("CNA segment input file (" + os.path.abspath(input_cna_segments) + ") does not have the correct file extension (.tsv or .txt)")
         logger.error('')
         exit(0)

      input_cna_basename = os.path.basename(str(input_cna_segments))
      input_cna_dir = os.path.dirname(os.path.abspath(input_cna_segments))

      ## if output cna segments exist and overwrite not set
      output_cna_segments = os.path.join(str(output_dir_full),str(sample_id)) + '.pcgr.cna_segments.tsv.gz'
      if os.path.exists(output_cna_segments) and overwrite == 0:
         logger.error('')
         logger.error("Output files (e.g. " + str(output_cna_segments) + ") already exist - please specify different sample_id or add option --force_overwrite")
         logger.error('')
         exit(0)
   
   ## check the existence of base folder
   base_dir = os.path.abspath(base_pcgr_dir)
   if not os.path.isdir(base_dir):
      logger.error("")
      logger.error("Base directory (" + str(base_dir) + ") does not exist")
      logger.error("")
      exit(0)
   
   ## check the existence of data folder within the base folder
   db_dir = os.path.join(os.path.abspath(base_pcgr_dir),'data')
   if not os.path.isdir(db_dir):
      logger.error("")
      logger.error("Data directory (" + str(db_dir) + ") does not exist")
      logger.error("")
      exit(0)
   
   host_directories = {}
   host_directories['input_vcf_dir_host'] = input_vcf_dir
   host_directories['input_cna_dir_host'] = input_cna_dir
   host_directories['db_dir_host'] = db_dir
   host_directories['base_dir_host'] = base_dir
   host_directories['output_dir_host'] = output_dir_full
   host_directories['input_vcf_basename_host'] = input_vcf_basename
   host_directories['input_cna_basename_host'] = input_cna_basename
   
   return host_directories
   

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

def run_pcgr(host_directories, docker_image_version, logR_threshold_amplification, logR_threshold_homozygous_deletion, num_vcfanno_processes, num_vep_forks, sample_id):
   """
   Main function to run the PCGR workflow using Docker
   """
   
   ## set basic Docker run commands
   output_vcf = 'None'
   uid = os.getuid()
   vepdb_dir_host = os.path.join(str(host_directories['base_dir_host']),'data','.vep')

   input_vcf_docker = 'None'
   input_cna_segments_docker = 'None'
   if host_directories['input_vcf_basename_host'] != 'NA':
      input_vcf_docker = '/workdir/input_vcf/' + str(host_directories['input_vcf_basename_host'])
   if host_directories['input_cna_basename_host'] != 'NA':
      input_cna_segments_docker = '/workdir/input_cna/' + str(host_directories['input_cna_basename_host'])
   
   if host_directories['input_vcf_dir_host'] != 'NA':
      docker_command_run1 = "docker run -t -u " + str(uid) + " -v=" + str(host_directories['base_dir_host']) + ":/data -v=" + str(vepdb_dir_host) + ":/usr/local/share/vep/data -v=" + str(host_directories['input_vcf_dir_host']) + ":/workdir/input_vcf -v=" + str(host_directories['output_dir_host']) + ":/workdir/output -w=/workdir/output " + str(docker_image_version) + " sh -c \""
   if host_directories['input_cna_dir_host'] != 'NA':
      docker_command_run1 = "docker run -t -u " + str(uid) + " -v=" + str(host_directories['base_dir_host']) + ":/data -v=" + str(vepdb_dir_host) + ":/usr/local/share/vep/data -v=" + str(host_directories['input_cna_dir_host']) + ":/workdir/input_cna -v=" + str(host_directories['output_dir_host']) + ":/workdir/output -w=/workdir/output " + str(docker_image_version) + " sh -c \""
   if host_directories['input_vcf_dir_host'] != 'NA' and host_directories['input_cna_dir_host'] != 'NA':
      docker_command_run1 = "docker run -t -u " + str(uid) + " -v=" + str(host_directories['base_dir_host']) + ":/data -v=" + str(vepdb_dir_host) + ":/usr/local/share/vep/data -v=" + str(host_directories['input_cna_dir_host']) + ":/workdir/input_cna -v=" + str(host_directories['input_vcf_dir_host']) + ":/workdir/input_vcf -v=" + str(host_directories['output_dir_host']) + ":/workdir/output -w=/workdir/output " + str(docker_image_version) + " sh -c \""
   docker_command_run2 = "docker run -t -u " + str(uid) + " -v=" + str(host_directories['base_dir_host']) + ":/data -v=" + str(host_directories['output_dir_host']) + ":/workdir/output -w=/workdir " + str(docker_image_version) + " sh -c \""
   
   
   ## verify VCF and CNA segment file
   logger = getlogger('pcgr-validate-input')
   logger.info("STEP 0: Validate input data")
   vcf_validate_command = str(docker_command_run1) + "pcgr_check_input.py /data " + str(input_vcf_docker) + " " + str(input_cna_segments_docker) + "\""
   check_subprocess(vcf_validate_command)
   logger.info('Finished')
   
   if not input_vcf_docker == 'None':
      
      ## Define input, output and temporary file names
      output_vcf = '/workdir/output/' + str(sample_id) + '.pcgr.vcf.gz'
      input_vcf_pcgr_ready = '/workdir/output/' + re.sub(r'(\.vcf$|\.vcf\.gz$)','.pcgr_ready.vcf.gz',host_directories['input_vcf_basename_host'])
      vep_vcf = re.sub(r'(\.vcf$|\.vcf\.gz$)','.pcgr_vep.vcf',input_vcf_pcgr_ready)
      vep_vcfanno_vcf = re.sub(r'(\.vcf$|\.vcf\.gz$)','.pcgr_vep.vcfanno.vcf',input_vcf_pcgr_ready)
      vep_tmp_vcf = vep_vcf + '.tmp'
      vep_vcfanno_annotated_vcf = re.sub(r'\.vcfanno','.vcfanno.annotated',vep_vcfanno_vcf) + '.gz'
      
      ## VEP command
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
   
      ## vcfanno command
      print
      logger = getlogger('pcgr-vcfanno')
      logger.info("STEP 2: Annotation for precision oncology with pcgr-vcfanno (ClinVar, dbSNP, dbNSFP, 1000Genomes Project, ExAC, gnomAD, CiVIC, CBMDB, DoCM, COSMIC, Intogen_driver_mutations, ICGC)")
      pcgr_vcfanno_command = str(docker_command_run2) + "pcgr_vcfanno.py --num_processes " + str(num_vcfanno_processes) + " --dbsnp --dbnsfp --oneKG --docm --clinvar --exac --gnomad --icgc --civic --cbmdb --intogen_driver_mut --cosmic " + str(vep_vcf) + ".gz " + str(vep_vcfanno_vcf) + " /data\""
      check_subprocess(pcgr_vcfanno_command)
      logger.info("Finished")
      
      ## summarise command
      print
      logger = getlogger("pcgr-summarise")
      pcgr_summarise_command = str(docker_command_run2) + "pcgr_summarise.py " + str(vep_vcfanno_vcf) + ".gz /data\""
      logger.info("STEP 3: Cancer gene annotations with pcgr-summarise")
      check_subprocess(pcgr_summarise_command)
      logger.info("Finished")
      
      create_output_vcf_command1 = str(docker_command_run2) + 'mv ' + str(vep_vcfanno_annotated_vcf) + ' ' + str(output_vcf) + "\""
      create_output_vcf_command2 = str(docker_command_run2) + 'mv ' + str(vep_vcfanno_annotated_vcf) + '.tbi ' + str(output_vcf) + '.tbi' + "\""
      clean_command = str(docker_command_run2) + 'rm -f ' + str(vep_vcf) + '* ' + str(vep_vcfanno_annotated_vcf) + ' ' + str(vep_vcfanno_vcf) + '* ' +  str(input_vcf_pcgr_ready) + "* "  + "\""
      check_subprocess(create_output_vcf_command1)
      check_subprocess(create_output_vcf_command2)
      check_subprocess(clean_command)
  
   print
   
   ## Generation of HTML reports for VEP/vcfanno-annotated VCF and copy number segment file
   logger = getlogger('pcgr-writer')
   logger.info("STEP 4: Generation of output files")
   pcgr_report_command = str(docker_command_run1) + "/pcgr.R /workdir/output " + str(output_vcf) + " " + str(input_cna_segments_docker) + " "  + str(sample_id)  + " " + str(logR_threshold_amplification) + " " + str(logR_threshold_homozygous_deletion) + "\""
   check_subprocess(pcgr_report_command)
   logger.info("Finished")
   
   
   
if __name__=="__main__": __main__()

