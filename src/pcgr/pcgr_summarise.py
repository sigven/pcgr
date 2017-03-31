#!/usr/bin/env python

import csv
import re
import argparse
import pcgr
from itertools import izip, imap
import transcript
import uniprot
import cyvcf
import gzip
import dbnsfp
import os
import utils
import vcfutils

logger = pcgr.getlogger('pcgr-gene-annotate')
csv.field_size_limit(500 * 1024 * 1024)
threeLettertoOneLetterAA = {'Ala':'A','Arg':'R','Asn':'N','Asp':'D','Cys':'C','Glu':'E','Gln':'Q','Gly':'G','His':'H','Ile':'I','Leu':'L','Lys':'K', 'Met':'M','Phe':'F','Pro':'P','Ser':'S','Thr':'T','Trp':'W','Tyr':'Y','Val':'V','Ter':'X'}


def __main__():
   
   parser = argparse.ArgumentParser(description='Cancer gene annotations from PCGR pipeline (SNVs/InDels)')
   parser.add_argument('vcf_file', help='VCF file with VEP-annotated query variants (SNVs/InDels)')
   parser.add_argument('pcgr_dir',help='PCGR base directory')
   args = parser.parse_args()

   extend_vcf_annotations(args.vcf_file, args.pcgr_dir)


def threeToOneAA(aa_change):
	
	for three_letter_aa in threeLettertoOneLetterAA.keys():
		aa_change = aa_change.replace(three_letter_aa,threeLettertoOneLetterAA[three_letter_aa])

	aa_change = re.sub(r'[A-Z]{1}fsX([0-9]{1,}|\?)','fs',aa_change)
	return aa_change

def hamming_distance(str1, str2):
   """hamming1(str1, str2): Hamming distance. Count the number of differences
   between equal length strings str1 and str2."""
   assert len(str1) == len(str2)
   return sum(c1 != c2 for c1, c2 in izip(str1, str2))

def index_clinvar(clinvar_tsv_file_path):
   clinvar_xref = {}
   with gzip.open(clinvar_tsv_file_path, 'rb') as tsvfile:
      cv_reader = csv.DictReader(tsvfile, delimiter='\t')
      for rec in cv_reader:
         
         unique_traits = {}
         traits = ''
         if not rec.has_key('trait'):
            traits = rec['all_traits']
         else:
            traits = rec['trait']
         #traits = rec['all_traits']
         for t in traits.split(';'):
            t_lc = str(t).lower()
            unique_traits[t_lc] = 1
         origin = ''
         if not rec.has_key('variant_origin'):
            origin = rec['origin']
         else:
            origin = rec['variant_origin']
         
         traits_curated = ';'.join(unique_traits.keys())
         traits_origin = traits_curated + ' - ' + str(origin)
         
         clinvar_xref[rec['measureset_id']] = {}
         clinvar_xref[rec['measureset_id']]['phenotype_origin'] = traits_origin
         if rec['symbol'] == '-' or rec['symbol'] == 'more than 10':
            rec['symbol'] = 'NA'
         clinvar_xref[rec['measureset_id']]['genesymbol'] = rec['symbol']
         
   return clinvar_xref


def map_variant_effect_predictors(rec, vep_info_tags, protein_info_tags, variant_prediction_tags, algorithms):
	
	dbnsfp_predictions = dbnsfp.map_dbnsfp_predictions(rec.INFO['DBNSFP'], algorithms)
	variant_prediction_tags['EFFECT_PREDICTIONS'] = {}
	
	for alt_allele in vep_info_tags['Feature'].keys():
		variant_prediction_tags['EFFECT_PREDICTIONS'][alt_allele] = '.'
		gene_id = vep_info_tags['Gene'][alt_allele]
		consequence = vep_info_tags['Consequence'][alt_allele]

		dbnsfp_key = ''
		if protein_info_tags.has_key('HGVSp_short'):
			aa_change = protein_info_tags['HGVSp_short'][alt_allele]
			dbnsfp_key = gene_id + ':' + str(aa_change)
		else:
			if re.search(r'splice_site',consequence):
				dbnsfp_key = gene_id
				
		if dbnsfp_key != '':
			if dbnsfp_predictions.has_key(dbnsfp_key):
				variant_prediction_tags['EFFECT_PREDICTIONS'][alt_allele] = dbnsfp_predictions[dbnsfp_key]

def get_protein_change_info(up_xref, up_feature_xref, swissprot_features, vep_info_tags, protein_info_tags):
   
   for alt_allele in vep_info_tags['Feature'].keys():
      tmp = {}
      tmp['UNIPROT_ID'] = ''
      tmp['SEQ_MATCH'] = '.'
      tmp['SYMBOL'] = '.'
      tmp['AA_position'] = '.'
      aa_positions = {}
        
      protein_change = '.'
      
      protein_info_tags['CDS_CHANGE'] = {}
      protein_info_tags['CDS_CHANGE'][alt_allele] = '.'
      if vep_info_tags['HGVSc'][alt_allele] != '':
         if 'splice_acceptor_variant' in vep_info_tags['Consequence'][alt_allele] or 'splice_donor_variant' in vep_info_tags['Consequence'][alt_allele]:
            key = str(vep_info_tags['Consequence'][alt_allele]) + ':' + str(vep_info_tags['HGVSc'][alt_allele])
            protein_info_tags['CDS_CHANGE'][alt_allele] = key

      if vep_info_tags['Amino_acids'][alt_allele] == '.' or vep_info_tags['Protein_position'][alt_allele] == '.' or vep_info_tags['Protein_position'][alt_allele].startswith('-'):
         continue
   
      for m in ['HGVSp_short','CANCER_MUTATION_HOTSPOT','UNIPROT_ID']:
         protein_info_tags[m] = {}
         protein_info_tags[m][alt_allele] = '.'
      for m in ['UNIPROT_FEATURE','PROTEIN_POSITIONS']:
         protein_info_tags[m] = {}
         protein_info_tags[m][alt_allele] = {}


      uniprot.get_uniprot_data_by_transcript(up_xref, vep_info_tags['Feature'][alt_allele], tmp)
      if '/' in vep_info_tags['Protein_position'][alt_allele]:
         tmp['AA_position'] = vep_info_tags['Protein_position'][alt_allele].split('/')[0]
      uniprot.get_domains_features_by_aapos(up_feature_xref, tmp, qtype = 'feature')
      if tmp.has_key('UNIPROT_FEATURE'):
         for feature in tmp['UNIPROT_FEATURE'].split('&'):
            spkey = str(tmp['UNIPROT_ID']) + ':' + str(feature)
            if swissprot_features.has_key(spkey):
               if swissprot_features[spkey]['type_description'].startswith('Disulfide bond'):
                  disulfid_start = int(swissprot_features[spkey]['aa_start'])
                  disulfid_stop = int(swissprot_features[spkey]['aa_stop'])
                  for aa_pos in aa_positions.keys():
                     if aa_pos == disulfid_start or aa_pos == disulfid_stop:
                        protein_info_tags['UNIPROT_FEATURE'][alt_allele][str(tmp['UNIPROT_ID'])	 + ':' + swissprot_features[spkey]['feature_type'] + ':' + str(disulfid_start) + '-' + str(disulfid_stop)] = 1
               else:
                  protein_info_tags['UNIPROT_FEATURE'][alt_allele][str(tmp['UNIPROT_ID']) + ':' + swissprot_features[spkey]['feature_type'] + ':' + str(swissprot_features[spkey]['aa_start']) + '-' + str(swissprot_features[spkey]['aa_stop'])] = 1
      if vep_info_tags['HGVSp'][alt_allele] != '':
         if ':' in vep_info_tags['HGVSp'][alt_allele]:
            protein_identifier = vep_info_tags['HGVSp'][alt_allele].split(':')[0]
            if protein_identifier.startswith('ENSP'):
               protein_change_VEP = vep_info_tags['HGVSp'][alt_allele].split(':')[1]
               protein_change = threeToOneAA(protein_change_VEP)

      if 'synonymous_variant' in vep_info_tags['Consequence'][alt_allele]:
         protein_change = 'p.' + str(vep_info_tags['Amino_acids'][alt_allele]) + str(tmp['AA_position']) + str(vep_info_tags['Amino_acids'][alt_allele])
      if 'stop_lost' in vep_info_tags['Consequence'][alt_allele] and '/' in str(vep_info_tags['Amino_acids'][alt_allele]):
         protein_change = 'p.X' + str(tmp['AA_position']) + str(vep_info_tags['Amino_acids'][alt_allele].split('/')[1])
      if '-' in tmp['AA_position'] and len(tmp['AA_position']) > 1:
         aa_start = unicode(tmp['AA_position'].split('-')[0])
         aa_stop = unicode(tmp['AA_position'].split('-')[1])
         if aa_start.isnumeric() and aa_stop.isnumeric():
            j = int(aa_start)
            while j <= int(aa_stop):
               if not aa_positions.has_key(j):
                  aa_positions[j] = 1
               j = j + 1
      else:
         if tmp['AA_position'] != '-':
            position = unicode(tmp['AA_position'])
            if position.isnumeric():
               if not aa_positions.has_key(int(position)):
                  aa_positions[int(position)] = 1
                  
      protein_info_tags['HGVSp_short'][alt_allele] = protein_change
      protein_info_tags['PROTEIN_POSITIONS'][alt_allele] = aa_positions
   
      if vep_info_tags['Protein_position'][alt_allele] != '' and vep_info_tags['Amino_acids'][alt_allele] == '':
         protein_info_tags['PROTEIN_POSITIONS'][alt_allele] = aa_positions
      
      exon_number = 'NA'
      if vep_info_tags['EXON'][alt_allele] != '':
         if '/' in vep_info_tags['EXON'][alt_allele]:
            exon_number = str(vep_info_tags['EXON'][alt_allele]).split('/')[0]
         
      if vep_info_tags['HGVSc'][alt_allele] != '':
         if protein_change != '.':
            key = str(vep_info_tags['Consequence'][alt_allele]) + ':' + str(vep_info_tags['HGVSc'][alt_allele]) + ':exon' + str(exon_number) + ':' + str(protein_change)
            protein_info_tags['CDS_CHANGE'][alt_allele] = key


def get_gene_data(gene_xref, vep_info_tags, extended_gene_oncorelevance_tags, idtype = 'transcript'):
  
   for alt_allele in vep_info_tags['Feature'].keys():
      annotations = ['gene_biotype','ccds','entrezgene','principal_isoform_flag','ensembl_gene_id','symbol','cancer_census_somatic','oncoscore','cancer_census_germline','intogen_drivers','tsgene','ts_oncogene','antineoplastic_drugs_dgidb']
      gene_values = {}
      for ann in annotations:
         gene_values[ann] = {}
      
      tid = re.sub(r'\.[0-9]{1,}$','',vep_info_tags['Feature'][alt_allele])
      if gene_xref.has_key(tid):
         gene_trans_mappings = gene_xref[tid]
         for k in gene_trans_mappings:
            for ann in annotations:
               if k[ann] != 'NA':
                  gene_values[ann][re.sub(r' ','_',k[ann])] = 1
      else:
         if tid.startswith('NM_') or tid.startswith('ENST'):
            logger.info("Could not find gene cross-reference information for transcript " + str(tid))
   
      for ann in ['ENTREZ_ID','APPRIS','CANCER_CENSUS_SOMATIC','CANCER_CENSUS_GERMLINE','INTOGEN_DRIVER','ANTINEOPLASTIC_DRUG_INTERACTION','TUMOR_SUPPRESSOR','ONCOGENE','ONCOSCORE']:
         extended_gene_oncorelevance_tags[ann] = {}
         extended_gene_oncorelevance_tags[ann][alt_allele] = {}
      
      for v in gene_values['entrezgene'].keys():
         extended_gene_oncorelevance_tags['ENTREZ_ID'][alt_allele][v] = 1
      
      for v in gene_values['principal_isoform_flag'].keys():
         extended_gene_oncorelevance_tags['APPRIS'][alt_allele][v] = 1
      
      for v in gene_values['cancer_census_somatic'].keys():
         extended_gene_oncorelevance_tags['CANCER_CENSUS_SOMATIC'][alt_allele][v] = 1
      
      for v in gene_values['cancer_census_germline'].keys():
         extended_gene_oncorelevance_tags['CANCER_CENSUS_GERMLINE'][alt_allele][v] = 1
      
      for v in gene_values['intogen_drivers'].keys():
         extended_gene_oncorelevance_tags['INTOGEN_DRIVER'][alt_allele][v] = 1
         
      for v in gene_values['tsgene'].keys():
         extended_gene_oncorelevance_tags['TUMOR_SUPPRESSOR'][alt_allele][v] = 1
      
      for v in gene_values['ts_oncogene'].keys():
         extended_gene_oncorelevance_tags['ONCOGENE'][alt_allele][v] = 1
         
      for v in gene_values['antineoplastic_drugs_dgidb'].keys():
         extended_gene_oncorelevance_tags['ANTINEOPLASTIC_DRUG_INTERACTION'][alt_allele][v] = 1
      
      for v in gene_values['oncoscore'].keys():
         extended_gene_oncorelevance_tags['ONCOSCORE'][alt_allele][v] = 1
  

def extend_vcf_annotations(query_vcf,pcgr_directory):
   
   pfam_domain_names = uniprot.index_pfam_names(pcgr_directory + '/data/pfam/pfam.domains.tsv.gz', ignore_versions = True)
   swissprot_features = uniprot.index_uniprot_feature_names(pcgr_directory + '/data/uniprot/uniprot.features.tsv.gz')
   pfam_xref = uniprot.index_pfam(pcgr_directory + '/data/pfam/pfam.uniprot.tsv.gz')
   uniprot_feature_xref = uniprot.index_uniprot_features(pcgr_directory + '/data/uniprot/uniprot.features.tsv.gz')

   gene_xref = None
   up_xref = None
   gene_xref = transcript.index_gene(pcgr_directory + '/data/gene.transcript.onco_xref.GRCh37.tsv.gz', index = 'ensGene_transcript')
   up_xref = uniprot.index_uniprot(pcgr_directory + '/data/uniprot/uniprot.xref.tsv.gz',index = 'ensGene_transcript')
   clinvar_data = index_clinvar(pcgr_directory + '/data/clinvar/clinvar.tsv.gz')
   out_prefix = re.sub(r'\.vcf(\.gz){0,}$','.annotated.vcf',query_vcf)
   cancer_hotspot_xref = utils.index_cancer_hotspots(pcgr_directory + '/data/cancerhotspots.org/cancer_hotspots.tsv')
  
   vep_infotags_desc = utils.read_infotag_file(pcgr_directory + '/data/vep_infotags.tsv')
   pcgr_infotags_desc = utils.read_infotag_file(pcgr_directory + '/data/pcgr_infotags.tsv')

   vcf_reader = cyvcf.Reader(open(query_vcf, 'r'))
   logger.info('Read query file')
   vep_csq_index2fields = {}
   vep_csq_fields2index = {}
   if 'CSQ' in vcf_reader.infos.keys():
      if 'Format: ' in vcf_reader.infos['CSQ'].desc:
         vep_tags = vcf_reader.infos['CSQ'].desc.split('Format: ')[1].split('|')
         i = 0
         for v in vep_tags:
            if vep_infotags_desc.has_key(v):
               vep_csq_index2fields[i] = v
               vep_csq_fields2index[v] = i
            i = i + 1
   else:
      logger.warning('VCF does not have CSQ tag in its meta information lines')
      no_csq = 1
   
   dbnsfp_prediction_algorithms = []
   effect_predictions_desc = ""
   if 'DBNSFP' in vcf_reader.infos.keys():
      if 'Format:' in vcf_reader.infos['DBNSFP'].desc:
         tmp = vcf_reader.infos['DBNSFP'].desc.split('Format:')[1].split('@')
         if len(tmp) > 7:
            effect_predictions_desc = "Format: " + '&'.join(tmp[7:])
         if len(tmp) == 1:
            ## v3.2
            tmp = vcf_reader.infos['DBNSFP'].desc.split('Format:')[1].split('#')
         i = 7
         while(i < len(tmp)):
            dbnsfp_prediction_algorithms.append(str(re.sub(r'((_score)|(_pred))$','',tmp[i])))
            i = i + 1

   if not 'EFFECT_PREDICTIONS' in vcf_reader.infos.keys():
      vcf_reader.infos['EFFECT_PREDICTIONS'] = ['EFFECT_PREDICTIONS','.',str(vcf_reader.infos['DBNSFP'].type), effect_predictions_desc]

   for tag in vep_infotags_desc:
      if not tag in vcf_reader.infos.keys():
         vcf_reader.infos[tag] = [tag,str(vep_infotags_desc[tag]['number']),str(vep_infotags_desc[tag]['type']),str(vep_infotags_desc[tag]['description'])]
   for tag in pcgr_infotags_desc:
      if not tag in vcf_reader.infos.keys():
         vcf_reader.infos[tag] = [tag,str(pcgr_infotags_desc[tag]['number']),str(pcgr_infotags_desc[tag]['type']),str(pcgr_infotags_desc[tag]['description'])]
   
   vcf_content = []
   current_chrom = None
   num_chromosome_records_processed = 0
   header_printed = 0
   for rec in vcf_reader:
      
      if not rec.INFO.has_key('CSQ'):
         variant_id = 'g.chr' + str(rec.CHROM) + ':' + str(rec.POS) + ':' + str(rec.REF) + '>' + str(rec.ALT)
         logger.warning('Variant record ' + str(variant_id) + ' does not have CSQ tag from Variant Effect Predictor - variant will be skipped')
         continue
      chromosome = rec.CHROM
      if current_chrom is None:
         current_chrom = chromosome
         num_chromosome_records_processed = 0
      else:
         if chromosome != current_chrom:
           
            if len(vcf_content) > 0:
               if header_printed == 0:
                  vcfutils.print_vcf_meta(out_prefix, vcf_reader)
                  header_printed = 1
               out = open(out_prefix, 'a')
               out.write('\n'.join(vcf_content))
               out.write('\n')
               out.close()
               vcf_content = []
            
            logger.info('Completed summary of functional annotations for ' + str(num_chromosome_records_processed) + ' variants on chromosome ' + str(current_chrom))
            current_chrom = chromosome
            num_chromosome_records_processed = 0
      num_chromosome_records_processed += 1
      
      fixed_fields_string = vcfutils.get_vcf_fixed_columns(rec)
      #sample_string = vcfutils.get_vcf_sample_columns(rec, vcf_reader)
      vep_info_tags = {}
      existing_info_tags = {}
      extended_gene_oncorelevance_tags = {}
      extended_protein_info_tags = {}
      variant_effect_prediction_tags = {}
      consequence_all = {}
      
      for keyw in sorted(vcf_reader.infos.keys()):
         if keyw == 'CSQ':
            num_picks = 0
            for csq in rec.INFO['CSQ'].split(','):
               csq_fields =  csq.split('|')
               alt_allele = str(csq_fields[0])
               if not consequence_all.has_key(alt_allele):
                  consequence_all[alt_allele] = {}
               if csq_fields[vep_csq_fields2index['PICK']] == "1":
                  num_picks += 1
                  j = 0
                  while(j < len(csq_fields)):
                     if vep_csq_index2fields.has_key(j):
                        vep_info_tags[vep_csq_index2fields[j]] = {}
                        if csq_fields[j] == '':
                           vep_info_tags[vep_csq_index2fields[j]][alt_allele] = '.'
                        else:
                           vep_info_tags[vep_csq_index2fields[j]][alt_allele] = str(csq_fields[j])
                     j = j + 1
               
               symbol = '.'
               if csq_fields[vep_csq_fields2index['SYMBOL']] != "":
                  symbol = str(csq_fields[vep_csq_fields2index['SYMBOL']])
               entry = str(csq_fields[vep_csq_fields2index['Consequence']]) + ':' + str(symbol) + ':' + str(csq_fields[vep_csq_fields2index['Feature_type']]) + ':' + str(csq_fields[vep_csq_fields2index['Feature']]) + ':' + str(csq_fields[vep_csq_fields2index['BIOTYPE']])
               consequence_all[alt_allele][entry] = 1
            vcfutils.get_info_value(rec, keyw, vcf_reader, existing_info_tags)
         else:
            if rec.INFO.has_key(keyw):
               vcfutils.get_info_value(rec, keyw, vcf_reader, existing_info_tags)
      
      vep_info_tags['VEP_ALL_CONSEQUENCE'] = {}
      for alt_allele in consequence_all.keys():
         vep_info_tags['VEP_ALL_CONSEQUENCE'][alt_allele] = ','.join(consequence_all[alt_allele].keys())
      if vep_info_tags.has_key('Protein_position') and vep_info_tags.has_key('Amino_acids'):
         get_protein_change_info(up_xref, uniprot_feature_xref, swissprot_features,vep_info_tags,extended_protein_info_tags)
      get_gene_data(gene_xref, vep_info_tags, extended_gene_oncorelevance_tags)
      if rec.INFO.has_key('DBNSFP'):
         map_variant_effect_predictors(rec, vep_info_tags, extended_protein_info_tags, variant_effect_prediction_tags, dbnsfp_prediction_algorithms)
      all_info_vals = []
      
      if extended_protein_info_tags.has_key('PROTEIN_POSITIONS'):
         utils.map_cancer_hotspots(cancer_hotspot_xref, vep_info_tags, extended_protein_info_tags)
       
      for k in vep_info_tags.keys():
         values = []
         for alt_allele in sorted(vep_info_tags[k].keys()):
            if vep_info_tags[k][alt_allele] != '.':
               values.append(vep_info_tags[k][alt_allele])
         if(len(values) > 0):
            all_info_vals.append(str(k) + '=' + ','.join(values))
      
      for k in variant_effect_prediction_tags.keys():
         values = []
         for alt_allele in sorted(variant_effect_prediction_tags[k].keys()):
            if variant_effect_prediction_tags[k][alt_allele] != '.':
               values.append(variant_effect_prediction_tags[k][alt_allele])
         if(len(values) > 0):
            all_info_vals.append(str(k) + '=' + ','.join(values))
      
      for k in existing_info_tags.keys():
         if k != 'DBNSFP':
            if existing_info_tags[k] == True:
               all_info_vals.append(str(k))
            else:
               all_info_vals.append(str(k) + '=' + str(existing_info_tags[k]))
         
      
      for k in extended_gene_oncorelevance_tags.keys():
         values = []
         for alt_allele in sorted(extended_gene_oncorelevance_tags[k].keys()):
            if len(extended_gene_oncorelevance_tags[k][alt_allele].keys()) != 0:
               values.append('&'.join(extended_gene_oncorelevance_tags[k][alt_allele].keys()))
         if(len(values) > 0):
            if k == 'ONCOGENE' or k == 'TUMOR_SUPPRESSOR':
               all_info_vals.append(str(k))
            else:
               all_info_vals.append(str(k) + '=' + ','.join(values))
         
      for k in extended_protein_info_tags.keys():
         values = []
         if k == 'PROTEIN_POSITIONS':
            continue
         for alt_allele in sorted(extended_protein_info_tags[k].keys()):
            if k == 'UNIPROT_FEATURE':
               if len(extended_protein_info_tags[k][alt_allele].keys()) != 0:
                  values.append('&'.join(extended_protein_info_tags[k][alt_allele].keys()))
            else:
               if extended_protein_info_tags[k][alt_allele] != '.':
                  values.append(extended_protein_info_tags[k][alt_allele])
         if(len(values) > 0):
            all_info_vals.append(str(k) + '=' + ','.join(values))

      info_string = ';'.join(all_info_vals)
      vcfline = fixed_fields_string + '\t' + str(info_string)
      #if sample_string != '.':
         #vcfline = vcfline + '\t' + str(sample_string)
      vcf_content.append(vcfline)
   
   if len(vcf_content) > 0:
      if header_printed == 0:
         vcfutils.print_vcf_meta(out_prefix, vcf_reader, print_sample_data = 0)
         header_printed = 1
      out = open(out_prefix, 'a')
      out.write('\n'.join(vcf_content))
      out.write('\n')
      out.close()
      vcf_content = []
   
      logger.info('Completed summary of functional annotations for ' + str(num_chromosome_records_processed) + ' variants on chromosome ' + str(current_chrom))
   
   os.system('bgzip -f ' + str(out_prefix))
   os.system('tabix -f -p vcf ' + str(out_prefix) + '.gz')

if __name__=="__main__": __main__()


      
