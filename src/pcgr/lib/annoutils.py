#!/usr/bin/env python

import os,re,sys
import csv
import logging
import gzip
import toml
from cyvcf2 import VCF, Writer


csv.field_size_limit(500 * 1024 * 1024)
threeLettertoOneLetterAA = {'Ala':'A','Arg':'R','Asn':'N','Asp':'D','Cys':'C','Glu':'E','Gln':'Q','Gly':'G','His':'H','Ile':'I','Leu':'L','Lys':'K', 'Met':'M','Phe':'F','Pro':'P','Ser':'S','Thr':'T','Trp':'W','Tyr':'Y','Val':'V','Ter':'X'}


def read_infotag_file(vcf_info_tags_tsv):
   """
   Function that reads a VCF info tag file that denotes annotation tags produced by PCGR/CPSR/gvanno.
   An example of the VCF info tag file is the following:
   
   tag	number	type	description category
   Consequence	.	String	"Impact modifier for the consequence type (picked by VEP's --flag_pick_allele option)."   vep
   
   A dictionary is returned, with the tag as the key, and the full dictionary record as the value
   """
   info_tag_xref = {} ##dictionary returned
   if not os.path.exists(vcf_info_tags_tsv):
      return info_tag_xref
   tsvfile = open(vcf_info_tags_tsv, 'r')
   reader = csv.DictReader(tsvfile, delimiter='\t')
   for row in reader:
      if not row['tag'] in info_tag_xref:
         info_tag_xref[row['tag']] = row
   
   return info_tag_xref


def error_message(message, logger):
   logger.error('')
   logger.error(message)
   logger.error('')
   exit(1)

# def verify_config_file(configuration_file, logger):
#    if not os.path.exists(os.path.abspath(configuration_file)):
#       err_msg = "Input file (" + str(configuration_file) + ") does not exist"
#       error_message(err_msg,logger)

#    if not os.path.abspath(configuration_file).endswith('.toml'):
#       err_msg = "Configuration file (" + os.path.abspath(configuration_file) + ") does not have the correct file extension (.toml)"
#       error_message(err_msg,logger)

#    config_host_info = {'basename':os.path.basename(str(configuration_file)), 'dir':os.path.dirname(os.path.abspath(configuration_file))}
#    return config_host_info

# def verify_input_vcf(input_vcf, logger):
#    if not os.path.exists(os.path.abspath(input_vcf)):
#       err_msg = "Input file (" + str(input_vcf) + ") does not exist"
#       error_message(err_msg,logger)

#    if not (os.path.abspath(input_vcf).endswith('.vcf') or os.path.abspath(input_vcf).endswith('.vcf.gz')):
#       err_msg = "VCF input file (" + os.path.abspath(input_vcf) + ") does not have the correct file extension (.vcf or .vcf.gz)"
#       error_message(err_msg,logger)

#    ## check that tabix file exist if bgzipped files is given
#    if os.path.abspath(input_vcf).endswith('.vcf.gz'):
#       tabix_file = input_vcf + '.tbi'
#       if not os.path.exists(os.path.abspath(tabix_file)):
#          err_msg = "Tabix file (i.e. '.gz.tbi') is not present for the bgzipped VCF input file (" + os.path.abspath(input_vcf) + "). Please make sure your input VCF is properly compressed and indexed (bgzip + tabix)"
#          error_message(err_msg,logger)

#    vcf_host_info = {'basename':os.path.basename(str(input_vcf)), 'dir':os.path.dirname(os.path.abspath(input_vcf))}
#    return vcf_host_info

# def verify_somatic_cna(input_cna, logger):
#    if not os.path.exists(os.path.abspath(input_cna)):
#       err_msg = "Input file (" + str(input_cna) + ") does not exist"
#       error_message(err_msg,logger)
#    if not (os.path.abspath(input_cna).endswith('.tsv') or os.path.abspath(input_cna).endswith('.txt')):
#       err_msg = "CNA segment input file (" + os.path.abspath(input_cna) + ") does not have the correct file extension (.tsv or .txt)"
#       error_message(err_msg,logger)
   
#    input_cna_basename = os.path.basename(str(input_cna))
#    input_cna_dir = os.path.dirname(os.path.abspath(input_cna))

#    cna_host_info = {'basename':os.path.basename(str(input_cna_basename)), 'dir':os.path.dirname(os.path.abspath(input_cna_dir))}
#    return cna_host_info

# def verify_output_folder(output_dir, logger):
#    ## check the existence of given output folder
#    output_dir_full = os.path.abspath(output_dir)
#    if not os.path.isdir(output_dir_full):
#       err_msg = "Output directory (" + str(output_dir_full) + ") does not exist"
#       error_message(err_msg,logger)

# def verify_databundle(directory, db_version, genome_assembly, logger, workflow = 'pcgr'):
#    ## check the existence of base folder
#    base_dir = os.path.abspath(directory)
#    if not os.path.isdir(base_dir):
#       err_msg = "Base directory (" + str(base_dir) + ") does not exist"
#       error_message(err_msg,logger)
   
#    ## check the existence of data folder within the base folder
#    db_dir = os.path.join(os.path.abspath(directory), 'data')
#    if not os.path.isdir(db_dir):
#       err_msg = "Data directory (" + str(db_dir) + ") does not exist"
#       error_message(err_msg,logger)
   
#    ## check the existence of specified assembly data folder within the base folder
#    db_assembly_dir = os.path.join(os.path.abspath(db_dir), genome_assembly)
#    if not os.path.isdir(db_assembly_dir):
#       err_msg = "Data directory for the specified genome assembly (" + str(db_assembly_dir) + ") does not exist"
#       error_message(err_msg,logger)
   
#    ## check the existence of RELEASE_NOTES (starting from 0.4.0)
#    rel_notes_file = os.path.join(os.path.abspath(db_assembly_dir), 'RELEASE_NOTES')
#    if not os.path.exists(rel_notes_file):
#       err_msg = 'The PCGR data bundle is outdated - please download the latest data bundle (see github.com/sigven/pcgr for instructions)'
#       if workflow == 'gvanno':
#          err_msg = 'The gvanno data bundle is outdated - please download the latest data bundle (see github.com/sigven/gvanno for instructions)'
#       if workflow == 'cpsr':
#          err_msg = 'The cpsr data bundle is outdated - please download the latest data bundle (see github.com/sigven/cpsr for instructions)'
#       error_message(err_msg,logger)
      
#    release_notes = open(rel_notes_file,'r')
#    compliant_data_bundle = 0
#    for line in release_notes:
#       if db_version in line:
#          compliant_data_bundle = 1
   
#    release_notes.close()

#    if compliant_data_bundle == 0:
#       err_msg = 'The PCGR data bundle is not compliant with the software version - please download the latest software and data bundle (see https://github.com/sigven/pcgr for instructions)'
#       if workflow == 'gvanno':
#           err_msg = 'The gvanno data bundle is not compliant with the software version - please download the latest software and data bundle (see https://github.com/sigven/gvanno for instructions)'
#       if workflow == 'cpsr':
#           err_msg = 'The cpsr data bundle is not compliant with the software version - please download the latest software and data bundle (see https://github.com/sigven/cpsr for instructions)'
#       error_message(err_msg,logger)
#    return

def write_pass_vcf(annotated_vcf, logger):
   
   out_vcf = re.sub(r'\.annotated\.vcf\.gz$','.annotated.pass.vcf',annotated_vcf)
   vcf = VCF(annotated_vcf)
   w = Writer(out_vcf, vcf)

   num_rejected = 0
   num_pass = 0
   for rec in vcf:
      if rec.FILTER is None or rec.FILTER == 'None':
         w.write_record(rec)
         num_pass += 1
      else:
         num_rejected +=1

   vcf.close()
   w.close()
   
   logger.info('Number of non-PASS/REJECTED variant calls: ' + str(num_rejected))
   logger.info('Number of PASSed variant calls: ' + str(num_pass))
   if num_pass == 0:
      logger.warning('There are zero variants with a \'PASS\' filter in the VCF file')
      os.system('bgzip -dc ' + str(annotated_vcf) + ' egrep \'^#\' > ' + str(out_vcf))
   #else:
   os.system('bgzip -f ' + str(out_vcf))
   os.system('tabix -f -p vcf ' + str(out_vcf) + '.gz')

   return

def warn_message(message, logger):
   logger.warning(message)

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


def read_config_options(configuration_file, base_dir, genome_assembly, logger, wflow = 'pcgr'):
   
   ## read default options
   config_options = {}
   configuration_file_default = os.path.join(base_dir,'data',str(genome_assembly),'pcgr_configuration_default.toml')  
   if wflow == 'cpsr':
      configuration_file_default = os.path.join(base_dir,'data',str(genome_assembly),'cpsr_configuration_default.toml')  
   if wflow == 'gvanno':
      configuration_file_default = os.path.join(base_dir,'data',str(genome_assembly),'gvanno_configuration_default.toml')  
   if not os.path.exists(configuration_file_default):
      err_msg = "Default configuration file " + str(configuration_file_default) + " does not exist - exiting"
      error_message(err_msg,logger)
   try:
      config_options = toml.load(configuration_file_default)
   except (IndexError,TypeError):
      err_msg = 'Default configuration file ' + str(configuration_file_default) + ' is not formatted correctly'
      error_message(err_msg, logger)

   ## override with options set by the users
   try:
      user_options = toml.load(configuration_file)
   except (IndexError,TypeError):
      err_msg = 'Configuration file ' + str(configuration_file) + ' is not formatted correctly'
      error_message(err_msg, logger)
   
   tumor_types = []
   for section in config_options:
      if section in user_options:
         for var in config_options[section]:
            if not var in user_options[section]:
               continue
            #print(str(section) + '\t' + str(var))
            if isinstance(config_options[section][var],bool) and not isinstance(user_options[section][var],bool):
               err_msg = 'Configuration value ' + str(user_options[section][var]) + ' for ' + str(var) + ' cannot be parsed properly (expecting boolean)'
               error_message(err_msg, logger)
            if isinstance(config_options[section][var],int) and not isinstance(user_options[section][var],int):
               err_msg = 'Configuration value \"' + str(user_options[section][var]) + '\" for ' + str(var) + ' cannot be parsed properly (expecting integer)'
               error_message(err_msg, logger)
            if isinstance(config_options[section][var],float) and (not isinstance(user_options[section][var],float) and not isinstance(user_options[section][var],int)):
               err_msg = 'Configuration value ' + str(user_options[section][var]) + ' for ' + str(var) + ' cannot be parsed properly (expecting float)'
               error_message(err_msg, logger)
            if isinstance(config_options[section][var],str) and not isinstance(user_options[section][var],str):
               err_msg = 'Configuration value ' + str(user_options[section][var]) + ' for ' + str(var) + ' cannot be parsed properly (expecting string)'
               error_message(err_msg, logger)
            if section == 'tumor_type':
               if user_options[section][var]:
                  tumor_types.append(str(var))
            tier_options = ['pcgr','pcgr_acmg']
            normalization_options = ['default','exome','genome','exome2genome']
            populations_tgp = ['afr','amr','eas','sas','eur','global']
            populations_gnomad = ['afr','amr','eas','sas','nfe','oth','fin','asj','global']
            theme_options = ['default', 'cerulean', 'journal', 'flatly', 'readable', 'spacelab', 'united', 'cosmo', 'lumen', 'paper', 'sandstone', 'simplex','yeti']
            if var == 'mutsignatures_normalization' and not str(user_options[section][var]) in normalization_options:
               err_msg = 'Configuration value ' + str(user_options[section][var]) + ' for ' + str(var) + ' cannot be parsed properly (expecting \'default\', \'exome\', \'genome\', or \'exome2genome\')'
               error_message(err_msg, logger)
            if var == 'mutsignatures_cutoff' and (float(user_options[section][var]) > 1 or float(user_options[section][var]) < 0) :
               err_msg = 'Configuration value ' + str(user_options[section][var]) + " must be within the [0,1] range"
               error_message(err_msg, logger)
            if var == 'mutsignatures_signature_limit':
               if int(user_options[section][var]) < 1 or int(user_options[section][var]) > 30:
                  err_msg = "Number of mutational signatures in search space ('mutsignatures_signature_limit') must be positive and not more than 30 (retrieved value: " + str(user_options[section][var]) + ")"
                  error_message(err_msg,logger)
            if var == 'tier_model' and not str(user_options[section][var]) in tier_options:
               err_msg = 'Configuration value \'' + str(user_options[section][var]) + '\' for ' + str(var) + ' cannot be parsed properly (expecting \'pcgr\', or \'pcgr_acmg\')'
               error_message(err_msg, logger)
            if var == 'pop_gnomad' and not str(user_options[section][var]) in populations_gnomad:
               err_msg = 'Configuration value \'' + str(user_options[section][var]) + '\' for ' + str(var) + ' cannot be parsed properly (expecting \'afr\', \'amr\', \'asj\', \'eas\', \'fin\', \'global\', \'nfe\', \'oth\', or \'sas\')'
               error_message(err_msg, logger)
            if var == 'pop_tgp' and not str(user_options[section][var]) in populations_tgp:
               err_msg = 'Configuration value \'' + str(user_options[section][var]) + '\' for ' + str(var) + ' cannot be parsed properly (expecting \'afr\', \'amr\', \'eas\', \'eur\', \'global\', or \'sas\')'
               error_message(err_msg, logger)
            if var == 'report_theme' and not str(user_options[section][var]) in theme_options:
               err_msg = 'Configuration value \'' + str(user_options[section][var]) + '\' for ' + str(var) + ' cannot be parsed properly (expecting \'default\', \'cerulean\', \'journal\', \'flatly\', \'readable\', \'spacelab\', \'united\', \'cosmo\', \'lumen\', \'paper\', \'sandstone\', \'simplex\', or \'yeti\')'
               error_message(err_msg, logger)
            if var.startswith('maf_'):
               if user_options[section][var] < 0 or user_options[section][var] > 1:
                  err_msg = "MAF value: " + str(var) + " must be within the [0,1] range, current value is " + str(user_options[section][var]) + ")"
                  error_message(err_msg,logger)
            if var == 'min_af_tumor':
               if user_options[section][var] < 0 or user_options[section][var] > 1:
                  err_msg = "Minimum AF tumor: " + str(var) + " must be within the [0,1] range, current value is " + str(user_options[section][var]) + ")"
                  error_message(err_msg,logger)
            if var == 'max_af_normal':
               if user_options[section][var] < 0 or user_options[section][var] > 1:
                  err_msg = "Maximum AF normal: " + str(var) + " must be within the [0,1] range, current value is " + str(user_options[section][var]) + ")"
                  error_message(err_msg,logger)
            if var == 'target_size_mb':
               if user_options[section][var] < 0 or user_options[section][var] > 50:
                  err_msg = "Target size region in Mb (" + str(user_options[section][var]) + ") is not positive or larger than the likely maximum size of the coding human genome (50 Mb))"
                  error_message(err_msg,logger)
               if user_options[section][var] < 1:
                  warn_msg = "Target size region in Mb (" + str(user_options[section][var]) + ") must be greater than 1 for mutational burden estimate to be robust (ignoring TMB calculation)"
                  warn_message(warn_msg,logger)
                  config_options[section]['mutational_burden'] = False
            if var == 'cna_overlap_pct':
               if user_options['cna'][var] > 100 or user_options['cna'][var] <= 0:
                  err_msg = "Minimum percent overlap between copy number segment and gene transcript (" + str(user_options[section][var]) + ") should be greater than zero and less than 100"
                  error_message(err_msg,logger)
            if var == 'logR_homdel':
               if user_options['cna'][var] > 0:
                  err_msg = "Log ratio for homozygous deletions (" + str(user_options[section][var]) + ") should be less than zero"
                  error_message(err_msg,logger)
            if var == 'logR_gain':
               if user_options['cna'][var] < 0:
                  err_msg = "Log ratio for copy number amplifications (" + str(user_options[section][var]) + ") should be greater than zero"
                  error_message(err_msg,logger)
            if var == 'min_majority' and section == 'dbnsfp':
               if user_options['dbnsfp'][var] > 8 or user_options['dbnsfp'][var] < 5:
                  err_msg = "Minimum number of majority votes for consensus calls among dbNSFP predictions should not exceed 8 and should not be less than 5"
                  error_message(err_msg,logger)
            if var == 'max_minority' and section == 'dbnsfp':
               if user_options['dbnsfp'][var] >= config_options['dbnsfp']['min_majority'] or user_options['dbnsfp'][var] > 2 or user_options['dbnsfp'][var] < 0 or (user_options['dbnsfp'][var] + config_options['dbnsfp']['min_majority'] > 8):
                  err_msg = "Maximum number of minority votes for consensus calls among dbNSFP predictions should not exceed 2 (8 algorithms in total) and should be less than min_majority (" + str(config_options['dbnsfp']['min_majority']) + ")"
                  error_message(err_msg,logger)
            
            config_options[section][var] = user_options[section][var]

   if len(tumor_types) > 2:
      err_msg = "Two many tumor types (", str(",".join(tumor_types)) + ")  set to True - limit is set to two"
      error_message(err_msg,logger)
   if wflow == 'pcgr':
      if 'msi' in config_options.keys() and 'mutational_burden' in config_options.keys():
         if config_options['msi']['msi'] == 1 and config_options['mutational_burden']['mutational_burden'] == 0:
            err_msg = "Prediction of MSI status (msi = true) requires mutational burden/target_size input (mutational_burden = true) - this is currently set as false"
            error_message(err_msg,logger)

   return config_options


def threeToOneAA(aa_change):
	
   for three_letter_aa in threeLettertoOneLetterAA.keys():
      aa_change = aa_change.replace(three_letter_aa,threeLettertoOneLetterAA[three_letter_aa])

   aa_change = re.sub(r'[A-Z]{1}fsX([0-9]{1,}|\?)','fs',aa_change)
   return aa_change

def map_variant_effect_predictors(rec, algorithms):
    
   dbnsfp_predictions = map_dbnsfp_predictions(str(rec.INFO.get('DBNSFP')), algorithms)
   if rec.INFO.get('Gene') is None or rec.INFO.get('Consequence') is None:
      return
   gene_id = str(rec.INFO.get('Gene'))
   consequence = str(rec.INFO.get('Consequence'))
     
   dbnsfp_key = ''

   found_key = 0
   if not rec.INFO.get('HGVSp_short') is None and not rec.INFO.get('HGVSp_short') == '.':
      aa_change = str(rec.INFO.get('HGVSp_short'))
      dbnsfp_key = gene_id + ':' + str(aa_change)
      if dbnsfp_key in dbnsfp_predictions:
         found_key = 1
   
   if found_key == 0 and re.search('splice_',consequence):
      dbnsfp_key = gene_id

   if dbnsfp_key != '':
      if dbnsfp_key in dbnsfp_predictions:
         rec.INFO['EFFECT_PREDICTIONS'] = dbnsfp_predictions[dbnsfp_key]
         for algo_pred in rec.INFO['EFFECT_PREDICTIONS'].split('&'):
            if algo_pred.startswith('sift:'):
               rec.INFO['SIFT_DBNSFP'] = str(algo_pred.split(':')[1])
            if algo_pred.startswith('provean:'):
               rec.INFO['PROVEAN_DBNSFP'] = str(algo_pred.split(':')[1])
            if algo_pred.startswith('m-cap:'):
               rec.INFO['M_CAP_DBNSFP'] = str(algo_pred.split(':')[1])
            if algo_pred.startswith('mutpred:'):
               rec.INFO['MUTPRED_DBNSFP'] = str(algo_pred.split(':')[1])
            if algo_pred.startswith('metalr:'):
               rec.INFO['META_LR_DBNSFP'] = str(algo_pred.split(':')[1])
            if algo_pred.startswith('fathmm:'):
               rec.INFO['FATHMM_DBNSFP'] = str(algo_pred.split(':')[1])
            if algo_pred.startswith('fathmm_mkl_coding:'):
               rec.INFO['FATHMM_MKL_DBNSFP'] = str(algo_pred.split(':')[1])
            if algo_pred.startswith('mutationtaster:'):
               rec.INFO['MUTATIONTASTER_DBNSFP'] = str(algo_pred.split(':')[1])
            if algo_pred.startswith('mutationassessor:'):
               rec.INFO['MUTATIONASSESSOR_DBNSFP'] = str(algo_pred.split(':')[1])
            if algo_pred.startswith('splice_site_rf:'):
               rec.INFO['SPLICE_SITE_RF_DBNSFP'] = str(algo_pred.split(':')[1])
            if algo_pred.startswith('splice_site_ada:'):
               rec.INFO['SPLICE_SITE_ADA_DBNSFP'] = str(algo_pred.split(':')[1])


def detect_reserved_info_tag(tag, tag_name, logger):
   reserved_tags = ['AA','AC','AF','AN','BQ','CIGAR','DB','DP','END','H2','H3','MQ','MQ0','NS','SB','SOMATIC','VALIDATED','1000G']
   if tag in reserved_tags:
      err_msg = 'Custom INFO tag (' + str(tag_name) + ') needs another name - ' + str(tag) + ' is a reserved field in the VCF specification (INFO)'
      return error_message(err_msg, logger)
   
   reserved_format_tags = ['GT','DP','FT','GL','GLE','GQ','PL','HQ','PS','PQ','EC','MQ']
   if tag in reserved_format_tags:
      err_msg = 'Custom INFO tag (' + str(tag_name) + ') needs another name - ' + str(tag) + ' is a reserved field in the VCF specification (FORMAT)'
      return error_message(err_msg, logger)

def set_coding_change(rec):
   
   for m in ['HGVSp_short','CDS_CHANGE']:
      rec.INFO[m] = '.'
   if not rec.INFO.get('HGVSc') is None:
      if rec.INFO.get('HGVSc') != '.':
         if 'splice_acceptor_variant' in rec.INFO.get('Consequence') or 'splice_donor_variant' in rec.INFO.get('Consequence'):
            key = str(rec.INFO.get('Consequence')) + ':' + str(rec.INFO.get('HGVSc'))
            rec.INFO['CDS_CHANGE'] = key
   if rec.INFO.get('Amino_acids') is None or rec.INFO.get('Protein_position') is None or rec.INFO.get('Consequence') is None:
      return
   if not rec.INFO.get('Protein_position') is None:
      if rec.INFO.get('Protein_position').startswith('-'):
         return

   protein_change = '.'
   if '/' in rec.INFO.get('Protein_position'):
      protein_position = str(rec.INFO.get('Protein_position').split('/')[0])
      if '-' in protein_position:
         if protein_position.split('-')[0].isdigit():
            rec.INFO['AMINO_ACID_START'] = protein_position.split('-')[0]
         if protein_position.split('-')[1].isdigit():
            rec.INFO['AMINO_ACID_END'] = protein_position.split('-')[1]
      else:
         if protein_position.isdigit():
            rec.INFO['AMINO_ACID_START'] = protein_position
            rec.INFO['AMINO_ACID_END'] = protein_position
   
   if not rec.INFO.get('HGVSp') is None:
      if rec.INFO.get('HGVSp') != '.':
         if ':' in rec.INFO.get('HGVSp'):
            protein_identifier = str(rec.INFO.get('HGVSp').split(':')[0])
            if protein_identifier.startswith('ENSP'):
               protein_change_VEP = str(rec.INFO.get('HGVSp').split(':')[1])
               protein_change = threeToOneAA(protein_change_VEP)
  
   if 'synonymous_variant' in rec.INFO.get('Consequence'):
      protein_change = 'p.' + str(rec.INFO.get('Amino_acids')) + str(protein_position) + str(rec.INFO.get('Amino_acids'))
      if 'stop_lost' in str(rec.INFO.get('Consequence')) and '/' in str(rec.INFO.get('Amino_acids')):
         protein_change = 'p.X' + str(protein_position) + str(rec.INFO.get('Amino_acids')).split('/')[1]
    
   rec.INFO['HGVSp_short'] = protein_change
   exon_number = 'NA'
   if not rec.INFO.get('EXON') is None:
      if rec.INFO.get('EXON') != '.':
         if '/' in rec.INFO.get('EXON'):
            exon_number = str(rec.INFO.get('EXON')).split('/')[0]
  
   if not rec.INFO.get('HGVSc') is None:
      if rec.INFO.get('HGVSc') != '.':
         if protein_change != '.':
            key = str(rec.INFO.get('Consequence')) + ':' + str(rec.INFO.get('HGVSc')) + ':exon' + str(exon_number) + ':' + str(protein_change)
            rec.INFO['CDS_CHANGE'] = key

   return

def map_dbnsfp_predictions(dbnsfp_tag, algorithms):
   
   effect_predictions = {}
   
   for v in dbnsfp_tag.split(','):
   
      dbnsfp_info = v.split('|')
      if len(dbnsfp_info) == 1:
         return effect_predictions
      ref_aa = dbnsfp_info[0]
      alt_aa = dbnsfp_info[1]
      all_ids = dbnsfp_info[4].split('&')
      unique_ids = {}
      for s in all_ids:
         unique_ids[s] = 1
         
      isoform_aa_keys = []
      if ref_aa != '.' and alt_aa != '.' and ref_aa != '' and alt_aa != '':
         aa_pos = dbnsfp_info[6].split('&')
         for pos in aa_pos:
            for gene_id in unique_ids:
               k = str(gene_id) + ':p.' + str(ref_aa) + pos + str(alt_aa)
               isoform_aa_keys.append(k)
      else:
         #continue
         for gene_id in unique_ids:
            isoform_aa_keys.append(gene_id)
   
      algorithm_raw_predictions = {}
   
      i = 7
      v = 0
      
      if len(algorithms) != len(dbnsfp_info[7:]):
         return effect_predictions
      
      while i < len(dbnsfp_info):
         algorithm_raw_predictions[str(algorithms[v]).lower()] = dbnsfp_info[i].split('&')
         i = i + 1
         v = v + 1
      dbnsfp_predictions = {}
      
      for k in isoform_aa_keys:
         if not k in dbnsfp_predictions:
            dbnsfp_predictions[k] = {}
         all_preds = []
         for algo in algorithm_raw_predictions.keys():
            unique_algo_predictions = {}
            for pred in algorithm_raw_predictions[algo]:
               if pred != '':
                  if not pred in unique_algo_predictions:
                     unique_algo_predictions[pred] = 1
               else:
                  unique_algo_predictions['.'] = 1
            if len(unique_algo_predictions.keys()) > 1 and '.' in unique_algo_predictions.keys():
               del unique_algo_predictions['.']
            dbnsfp_predictions[k][algo] = str(algo) + ':' + '|'.join(unique_algo_predictions.keys())  
            all_preds.append(dbnsfp_predictions[k][algo])
         effect_predictions[k] = '&'.join(all_preds)

   return effect_predictions


def vep_dbnsfp_meta_vcf(query_vcf, info_tags_wanted):
   vep_to_pcgr_af = {'gnomAD_AMR_AF':'AMR_AF_GNOMAD','gnomAD_AFR_AF':'AFR_AF_GNOMAD','gnomAD_EAS_AF':'EAS_AF_GNOMAD','gnomAD_NFE_AF':'NFE_AF_GNOMAD','gnomAD_AF':'GLOBAL_AF_GNOMAD',
                     'gnomAD_SAS_AF':'SAS_AF_GNOMAD','gnomAD_OTH_AF':'OTH_AF_GNOMAD','gnomAD_ASJ_AF':'ASJ_AF_GNOMAD','gnomAD_FIN_AF':'FIN_AF_GNOMAD','AFR_AF':'AFR_AF_1KG',
                     'AMR_AF':'AMR_AF_1KG','SAS_AF':'SAS_AF_1KG','EUR_AF':'EUR_AF_1KG','EAS_AF':'EAS_AF_1KG', 'AF':'GLOBAL_AF_1KG'}

   vcf = VCF(query_vcf)
   vep_csq_index2fields = {}
   vep_csq_fields2index = {}
   dbnsfp_prediction_algorithms = []
   for e in vcf.header_iter():
      header_element = e.info()
      if 'ID' in header_element.keys():
         identifier = str(header_element['ID'])
         if identifier == 'CSQ' or identifier == 'DBNSFP':
            description = str(header_element['Description'])
            if 'Format: ' in description:
               subtags = description.split('Format: ')[1].split('|')
               if identifier == 'CSQ':
                  i = 0
                  for t in subtags:
                     v = t
                     if t in vep_to_pcgr_af:
                        v = str(vep_to_pcgr_af[t])
                     if v in info_tags_wanted:
                        vep_csq_index2fields[i] = v
                        vep_csq_fields2index[v] = i
                     i = i + 1
               if identifier == 'DBNSFP':
                  #if len(subtags) > 7:
                     #effect_predictions_description = "Format: " + '|'.join(subtags[7:])
                  i = 7
                  while(i < len(subtags)):
                     dbnsfp_prediction_algorithms.append(str(re.sub(r'((_score)|(_pred))"*$','',subtags[i])))
                     i = i + 1

   vep_dbnsfp_meta_info = {}
   vep_dbnsfp_meta_info['vep_csq_index2fields'] = vep_csq_index2fields
   vep_dbnsfp_meta_info['vep_csq_fields2index'] = vep_csq_fields2index
   vep_dbnsfp_meta_info['dbnsfp_prediction_algorithms'] = dbnsfp_prediction_algorithms

   return vep_dbnsfp_meta_info