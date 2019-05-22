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
   
   valid_tumor_types = ['Adrenal_Gland_Cancer_NOS','Ampullary_Carcinoma_NOS','Biliary_Tract_Cancer_NOS','Bladder_Urinary_Tract_Cancer_NOS',
                        'Bone_Cancer_NOS','Breast_Cancer_NOS','Cancer_Unknown_Primary_NOS','Cervical_Cancer_NOS','CNS_Brain_Cancer_NOS',
                        'Colorectal_Cancer_NOS','Esophageal_Cancer_NOS','Head_And_Neck_Cancer_NOS','Kidney_Cancer','Leukemia_NOS',
                        'Liver_Cancer_NOS','Lung_Cancer_NOS','Lymphoma_Hodgkin_NOS','Lymphoma_Non_Hodgkin_NOS','Multiple_Myeloma_NOS',
                        'Ovarian_Fallopian_Tube_Cancer_NOS','Pancreatic_Cancer_NOS','Penile_Cancer_NOS','Peripheral_Nervous_System_Cancer_NOS',
                        'Peritoneal_Cancer_NOS','Pleural_Cancer_NOS','Prostate_Cancer_NOS','Skin_Cancer_NOS','Soft_Tissue_Cancer_NOS',
                        'Stomach_Cancer_NOS','Testicular_Cancer_NOS','Thymic_Cancer_NOS','Thyroid_Cancer_NOS','Uterine_Cancer_NOS',
                        'Vulvar_Vaginal_Cancer_NOS','']

   for section in config_options:
      if section in user_options:
         for var in config_options[section]:
            if not var in user_options[section]:
               continue
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
            if section == 'tumor_type' and var == 'type':
               if not str(user_options[section][var]) in valid_tumor_types:
                  err_msg('Configuration value for tumor type (' + str(user_options[section][var]) + ') is not a valid type')
                  error_message(err_msg, logger)
            #tier_options = ['pcgr','pcgr_acmg']
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
            # if var == 'tier_model' and not str(user_options[section][var]) in tier_options:
            #    err_msg = 'Configuration value \'' + str(user_options[section][var]) + '\' for ' + str(var) + \
            #       ' cannot be parsed properly (expecting \'pcgr\', or \'pcgr_acmg\')'
            #    error_message(err_msg, logger)
            if var == 'pop_gnomad' and not str(user_options[section][var]) in populations_gnomad:
               err_msg = 'Configuration value \'' + str(user_options[section][var]) + '\' for ' + str(var) + \
                  ' cannot be parsed properly (expecting \'afr\', \'amr\', \'asj\', \'eas\', \'fin\', \'global\', \'nfe\', \'oth\', or \'sas\')'
               error_message(err_msg, logger)
            if var == 'pop_tgp' and not str(user_options[section][var]) in populations_tgp:
               err_msg = 'Configuration value \'' + str(user_options[section][var]) + '\' for ' + str(var) + \
                  ' cannot be parsed properly (expecting \'afr\', \'amr\', \'eas\', \'eur\', \'global\', or \'sas\')'
               error_message(err_msg, logger)
            if var == 'report_theme' and not str(user_options[section][var]) in theme_options:
               err_msg = 'Configuration value \'' + str(user_options[section][var]) + '\' for ' + str(var) + \
                  ' cannot be parsed properly (expecting \'default\', \'cerulean\', \'journal\', \'flatly\', \'readable\', \'spacelab\', \'united\', \'cosmo\', \'lumen\', \'paper\', \'sandstone\', \'simplex\', or \'yeti\')'
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
               if user_options[section][var] < 0 or user_options[section][var] > 34:
                  err_msg = "Target size region in Mb (" + str(user_options[section][var]) + ") is not positive or larger than the likely maximum size of the coding human genome (34 Mb))"
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
            if var == 'vep_pick_order':
               values = str(user_options['other'][var]).split(',')
               permitted_sources = ['canonical','appris','tsl','biotype','ccds','rank','length']
               num_permitted_sources = 0
               for v in values:
                  if v in permitted_sources:
                     num_permitted_sources += 1
               
               if num_permitted_sources < 7:
                  err_msg = "Configuration value vep_pick_order = " + str(user_options['other']['vep_pick_order']) + " is formatted incorrectly should be a comma-separated string of the following values: canonical,appris,tsl,biotype,ccds,rank,length"
                  error_message(err_msg, logger)
            config_options[section][var] = user_options[section][var]

   
   if wflow == 'pcgr':
      #if config_options['tumor_type']['type'] == '':
         #err_msg = "Tumor type not defined - please specify a tumor type in the configuration file ([tumor_type] section)"
         #error_message(err_msg,logger)
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
            if algo_pred.startswith('sift4g:'):
               rec.INFO['SIFT4G_DBNSFP'] = str(algo_pred.split(':')[1])
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
            if algo_pred.startswith('deogen2:'):
               rec.INFO['DEOGEN2_DBNSFP'] = str(algo_pred.split(':')[1])
            if algo_pred.startswith('primateai:'):
               rec.INFO['PRIMATEAI_DBNSFP'] = str(algo_pred.split(':')[1])
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


def assign_cds_exon_intron_annotations(csq_record):
   csq_record['CODING_STATUS'] = 'noncoding'
   csq_record['EXONIC_STATUS'] = 'nonexonic'
   csq_record['SPLICE_DONOR_RELEVANT'] = False
   csq_record['NULL_VARIANT'] = False
   
   coding_csq_pattern = r"^(stop_|start_lost|frameshift_|missense_|splice_donor|splice_acceptor|protein_altering|inframe_)"
   wes_csq_pattern = r"^(stop_|start_lost|frameshift_|missense_|splice_donor|splice_acceptor|inframe_|protein_altering|synonymous)"
   null_pattern = r"^(stop_|frameshift_)"
   if re.match(coding_csq_pattern, str(csq_record['Consequence'])):
      csq_record['CODING_STATUS'] = 'coding'
   
   if re.match(wes_csq_pattern, str(csq_record['Consequence'])):
      csq_record['EXONIC_STATUS'] = 'exonic'
   
   if re.match(null_pattern, str(csq_record['Consequence'])):
      csq_record['NULL_VARIANT'] = True

   if re.match(r"^splice_region_variant", str(csq_record['Consequence'])) and re.search(r'(\+3(A|G)>|\+4G>|\+5G>)', str(csq_record['HGVSc'])):
      csq_record['SPLICE_DONOR_RELEVANT'] = True

   for m in ['HGVSp_short','CDS_CHANGE']:
      csq_record[m] = '.'
   if not csq_record['HGVSc'] is None:
      if csq_record['HGVSc'] != '.':
         if 'splice_acceptor_variant' in csq_record['Consequence'] or 'splice_donor_variant' in csq_record['Consequence']:
            key = str(csq_record['Consequence']) + ':' + str(csq_record['HGVSc'])
            csq_record['CDS_CHANGE'] = key
   if csq_record['Amino_acids'] is None or csq_record['Protein_position'] is None or csq_record['Consequence'] is None:
      return
   if not csq_record['Protein_position'] is None:
      if csq_record['Protein_position'].startswith('-'):
         return

   protein_change = '.'
   if '/' in csq_record['Protein_position']:
      protein_position = str(csq_record['Protein_position'].split('/')[0])
      if '-' in protein_position:
         if protein_position.split('-')[0].isdigit():
            csq_record['AMINO_ACID_START'] = protein_position.split('-')[0]
         if protein_position.split('-')[1].isdigit():
            csq_record['AMINO_ACID_END'] = protein_position.split('-')[1]
      else:
         if protein_position.isdigit():
            csq_record['AMINO_ACID_START'] = protein_position
            csq_record['AMINO_ACID_END'] = protein_position
   
   if not csq_record['HGVSp'] is None:
      if csq_record['HGVSp'] != '.':
         if ':' in csq_record['HGVSp']:
            protein_identifier = str(csq_record['HGVSp'].split(':')[0])
            if protein_identifier.startswith('ENSP'):
               protein_change_VEP = str(csq_record['HGVSp'].split(':')[1])
               protein_change = threeToOneAA(protein_change_VEP)
  
   if 'synonymous_variant' in csq_record['Consequence']:
      protein_change = 'p.' + str(csq_record['Amino_acids']) + str(protein_position) + str(csq_record['Amino_acids'])
      if 'stop_lost' in str(csq_record['Consequence']) and '/' in str(csq_record['Amino_acids']):
         protein_change = 'p.X' + str(protein_position) + str(csq_record['Amino_acids']).split('/')[1]
    
   csq_record['HGVSp_short'] = protein_change
   exon_number = 'NA'
   if not csq_record['EXON'] is None:
      if csq_record['EXON'] != '.':
         if '/' in csq_record['EXON']:
            exon_number = str(csq_record['EXON']).split('/')[0]
            num_exons = str(csq_record['EXON']).split('/')[1]
            if exon_number == num_exons:
               csq_record['LAST_EXON'] = True

   if not csq_record['INTRON'] is None:
      if csq_record['INTRON'] != '.':
         if '/' in csq_record['INTRON']:
            intron_number = str(csq_record['INTRON']).split('/')[0]
            num_introns = str(csq_record['INTRON']).split('/')[1]
            if intron_number == num_introns:
               csq_record['LAST_INTRON'] = True

   if not csq_record['HGVSc'] is None:
      if csq_record['HGVSc'] != '.':
         if protein_change != '.':
            key = str(csq_record['Consequence']) + ':' + str(csq_record['HGVSc']) + ':exon' + str(exon_number) + ':' + str(protein_change)
            csq_record['CDS_CHANGE'] = key

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
      
      #print(str(algorithms))
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


def make_transcript_xref_map(rec, fieldmap, xref_tag = 'PCGR_ONCO_XREF'):
   transcript_xref_map = {}
   if not rec.INFO.get(xref_tag) is None:
      for tref in rec.INFO.get(xref_tag).split(','):
         xrefs = tref.split('|')
         ensembl_transcript_id = str(xrefs[0])
         transcript_xref_map[ensembl_transcript_id] = {}
         for annotation in fieldmap.keys():
            annotation_index = fieldmap[annotation]
            if annotation_index > (len(xrefs) - 1):
               continue
            if xrefs[annotation_index] != '':
               transcript_xref_map[ensembl_transcript_id][annotation] = xrefs[annotation_index]
   
   return(transcript_xref_map)

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
                  i = 7
                  while(i < len(subtags)):
                     dbnsfp_prediction_algorithms.append(str(re.sub(r'((_score)|(_pred))"*$','',subtags[i])))
                     i = i + 1

   vep_dbnsfp_meta_info = {}
   vep_dbnsfp_meta_info['vep_csq_fieldmap'] = {}
   vep_dbnsfp_meta_info['vep_csq_fieldmap']['field2index'] = vep_csq_fields2index
   vep_dbnsfp_meta_info['vep_csq_fieldmap']['index2field'] = vep_csq_index2fields
   vep_dbnsfp_meta_info['dbnsfp_prediction_algorithms'] = dbnsfp_prediction_algorithms

   return vep_dbnsfp_meta_info

def parse_vep_csq(rec, transcript_xref_map, vep_csq_fields_map, logger, pick_only = True, csq_identifier = 'CSQ'):

   all_csq_pick = []
   all_transcript_consequences = []

   for csq in rec.INFO.get(csq_identifier).split(','):
      csq_fields =  csq.split('|')
      ## loop over VEP consequence blocks PICK'ed according to VEP's ranking scheme
      if csq_fields[vep_csq_fields_map['field2index']['PICK']] == "1": ## only consider the primary/picked consequence when expanding with annotation tags
         j = 0
         csq_record = {}

         ## loop over block annotation elements (separated with '|'), and assign their values to the csq_record dictionary object
         while(j < len(csq_fields)):
            if j in vep_csq_fields_map['index2field']:
               if csq_fields[j] != '':
                  csq_record[vep_csq_fields_map['index2field'][j]] = str(csq_fields[j])
                  if vep_csq_fields_map['index2field'][j] == 'Feature':
                     ensembl_transcript_id = str(csq_fields[j])
                     if ensembl_transcript_id in transcript_xref_map:
                        for annotation in transcript_xref_map[ensembl_transcript_id].keys():
                           if annotation != 'SYMBOL':
                              ## assign additional gene/transcript annotations from the custom transcript xref map (PCGR/CPSR) as key,value pairs in the csq_record object
                              csq_record[annotation] = transcript_xref_map[ensembl_transcript_id][annotation]
                     else:
                        logger.warning('Could not find transcript xrefs for ' + str(ensembl_transcript_id))

                  ## Specifically assign PFAM protein domain as a csq_record key
                  if vep_csq_fields_map['index2field'][j] == 'DOMAINS':
                     domain_identifiers = str(csq_fields[j]).split('&')
                     for v in domain_identifiers:
                        if v.startswith('Pfam_domain'):
                           csq_record['PFAM_DOMAIN'] = str(re.sub(r'\.[0-9]{1,}$','',re.sub(r'Pfam_domain:','',v)))
                     
                     csq_record['DOMAINS'] = None
                  ## Assign COSMIC/DBSNP mutation ID's as individual key,value pairs in the csq_record object
                  if vep_csq_fields_map['index2field'][j] == 'Existing_variation':
                     var_identifiers = str(csq_fields[j]).split('&')
                     cosmic_identifiers = []
                     dbsnp_identifiers = []
                     for v in var_identifiers:
                        if v.startswith('COSM'):
                           cosmic_identifiers.append(v)
                        if v.startswith('rs'):
                           dbsnp_identifiers.append(v)
                     if len(cosmic_identifiers) > 0:
                        csq_record['COSMIC_MUTATION_ID'] = '&'.join(cosmic_identifiers)
                     if len(dbsnp_identifiers) > 0:
                        csq_record['DBSNPRSID'] = '&'.join(dbsnp_identifiers)
               else:
                  csq_record[vep_csq_fields_map['index2field'][j]] = None
            j = j + 1
         
         ## Assign coding status, protein change, coding sequence change, last exon/intron status etc
         assign_cds_exon_intron_annotations(csq_record)
         all_csq_pick.append(csq_record)
      symbol = '.'
      if csq_fields[vep_csq_fields_map['field2index']['SYMBOL']] != "":
         symbol = str(csq_fields[vep_csq_fields_map['field2index']['SYMBOL']])
      consequence_entry = (str(csq_fields[vep_csq_fields_map['field2index']['Consequence']]) + ':' +  
         str(symbol) + ':' + 
         str(csq_fields[vep_csq_fields_map['field2index']['Feature_type']]) + ':' + 
         str(csq_fields[vep_csq_fields_map['field2index']['Feature']]) + ':' + 
         str(csq_fields[vep_csq_fields_map['field2index']['BIOTYPE']]))
      all_transcript_consequences.append(consequence_entry)
         
   vep_csq_results = {}
   vep_csq_results['vep_block'] = all_csq_pick
   vep_csq_results['vep_all_csq'] = all_transcript_consequences

   return(vep_csq_results)
      