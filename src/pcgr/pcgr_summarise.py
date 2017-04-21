#!/usr/bin/env python

import csv
import re
import argparse
import pcgr
from itertools import izip, imap
import cyvcf
import gzip
import dbnsfp
import os
import pcgrutils
import vcfutils

logger = pcgrutils.getlogger('pcgr-gene-annotate')
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

def map_variant_effect_predictors(rec, vep_info_tags, pcgr_protein_info_tags, variant_prediction_tags, algorithms):
   
   dbnsfp_predictions = dbnsfp.map_dbnsfp_predictions(rec.INFO['DBNSFP'], algorithms)
   variant_prediction_tags['EFFECT_PREDICTIONS'] = {}
   
   for alt_allele in vep_info_tags['Feature'].keys():
      variant_prediction_tags['EFFECT_PREDICTIONS'][alt_allele] = '.'
      gene_id = vep_info_tags['Gene'][alt_allele]
      consequence = vep_info_tags['Consequence'][alt_allele]

      dbnsfp_key = ''
      if pcgr_protein_info_tags.has_key('HGVSp_short'):
         aa_change = pcgr_protein_info_tags['HGVSp_short'][alt_allele]
         dbnsfp_key = gene_id + ':' + str(aa_change)
      else:
         if re.search(r'splice_site',consequence):
            dbnsfp_key = gene_id
            
      if dbnsfp_key != '':
         if dbnsfp_predictions.has_key(dbnsfp_key):
            variant_prediction_tags['EFFECT_PREDICTIONS'][alt_allele] = dbnsfp_predictions[dbnsfp_key]

def set_pcgr_protein_info_tags(up_xref, up_feature_xref, cancer_hotspot_xref, uniprot_feature_names, vep_info_tags, pcgr_protein_info_tags):
   
   """
   Function that adds specific annotation tags of relevance for protein-coding alterations:
   1. Coding sequence change - 'CDS_CHANGE'
   2. UniProt functional feature, e.g. active site - 'UNIPROT_FEATURE'
   3. PFAM protein domain name '
   4. Short version of HGVS, e.g. p.V600E (as opposed to p.Val600Glu)
   5. Cancer mutation hotspots
   """
   
   for alt_allele in vep_info_tags['Feature'].keys():
      tmp = {}
      tmp['UNIPROT_ID'] = ''
      tmp['SEQ_MATCH'] = '.'
      tmp['SYMBOL'] = '.'
      tmp['AA_position'] = '.'
      aa_positions = {}
        
      protein_change = '.'
      
      pcgr_protein_info_tags['CDS_CHANGE'] = {}
      pcgr_protein_info_tags['CDS_CHANGE'][alt_allele] = '.'
      if vep_info_tags['HGVSc'][alt_allele] != '':
         if 'splice_acceptor_variant' in vep_info_tags['Consequence'][alt_allele] or 'splice_donor_variant' in vep_info_tags['Consequence'][alt_allele]:
            key = str(vep_info_tags['Consequence'][alt_allele]) + ':' + str(vep_info_tags['HGVSc'][alt_allele])
            pcgr_protein_info_tags['CDS_CHANGE'][alt_allele] = key

      if vep_info_tags['Amino_acids'][alt_allele] == '.' or vep_info_tags['Protein_position'][alt_allele] == '.' or vep_info_tags['Protein_position'][alt_allele].startswith('-'):
         continue
   
      for m in ['HGVSp_short','CANCER_MUTATION_HOTSPOT','UNIPROT_ID']:
         pcgr_protein_info_tags[m] = {}
         pcgr_protein_info_tags[m][alt_allele] = '.'
      for m in ['UNIPROT_FEATURE','PROTEIN_POSITIONS']:
         pcgr_protein_info_tags[m] = {}
         pcgr_protein_info_tags[m][alt_allele] = {}


      pcgrutils.get_uniprot_data_by_transcript(up_xref, vep_info_tags['Feature'][alt_allele], tmp)
      if '/' in vep_info_tags['Protein_position'][alt_allele]:
         tmp['AA_position'] = vep_info_tags['Protein_position'][alt_allele].split('/')[0]
      pcgrutils.get_domains_features_by_aapos(up_feature_xref, tmp, qtype = 'feature')
      if tmp.has_key('UNIPROT_FEATURE'):
         for feature in tmp['UNIPROT_FEATURE'].split('&'):
            spkey = str(tmp['UNIPROT_ID']) + ':' + str(feature)
            if uniprot_feature_names.has_key(spkey):
               if uniprot_feature_names[spkey]['type_description'].startswith('Disulfide bond'):
                  disulfid_start = int(uniprot_feature_names[spkey]['aa_start'])
                  disulfid_stop = int(uniprot_feature_names[spkey]['aa_stop'])
                  for aa_pos in aa_positions.keys():
                     if aa_pos == disulfid_start or aa_pos == disulfid_stop:
                        pcgr_protein_info_tags['UNIPROT_FEATURE'][alt_allele][str(tmp['UNIPROT_ID'])   + ':' + uniprot_feature_names[spkey]['feature_type'] + ':' + str(disulfid_start) + '-' + str(disulfid_stop)] = 1
               else:
                  pcgr_protein_info_tags['UNIPROT_FEATURE'][alt_allele][str(tmp['UNIPROT_ID']) + ':' + uniprot_feature_names[spkey]['feature_type'] + ':' + str(uniprot_feature_names[spkey]['aa_start']) + '-' + str(uniprot_feature_names[spkey]['aa_stop'])] = 1
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
                  
      pcgr_protein_info_tags['HGVSp_short'][alt_allele] = protein_change
      pcgr_protein_info_tags['PROTEIN_POSITIONS'][alt_allele] = aa_positions
   
      if vep_info_tags['Protein_position'][alt_allele] != '' and vep_info_tags['Amino_acids'][alt_allele] == '':
         pcgr_protein_info_tags['PROTEIN_POSITIONS'][alt_allele] = aa_positions
      
      if pcgr_protein_info_tags.has_key('PROTEIN_POSITIONS'):
         pcgrutils.map_cancer_hotspots(cancer_hotspot_xref, vep_info_tags, pcgr_protein_info_tags)
      
      exon_number = 'NA'
      if vep_info_tags['EXON'][alt_allele] != '':
         if '/' in vep_info_tags['EXON'][alt_allele]:
            exon_number = str(vep_info_tags['EXON'][alt_allele]).split('/')[0]
         
      if vep_info_tags['HGVSc'][alt_allele] != '':
         if protein_change != '.':
            key = str(vep_info_tags['Consequence'][alt_allele]) + ':' + str(vep_info_tags['HGVSc'][alt_allele]) + ':exon' + str(exon_number) + ':' + str(protein_change)
            pcgr_protein_info_tags['CDS_CHANGE'][alt_allele] = key


def set_pcgr_gene_info_tags(transcript_xref, vep_info_tags, pcgr_gene_info_tags):
  
   """
   Function that sets different cancer-relevant annotations (membership in Cancer Cene Census, predicted driver genes etc.) for a given transcript ID
   """
  
   for alt_allele in vep_info_tags['Feature'].keys():
      pcgr_gene_annotations = ['gene_biotype','ccds','entrezgene','principal_isoform_flag','ensembl_gene_id','symbol','cancer_census_somatic','oncoscore','cancer_census_germline','intogen_drivers','tsgene','ts_oncogene','antineoplastic_drugs_dgidb']
      gene_values = {}
      for ann in pcgr_gene_annotations:
         gene_values[ann] = {}
      
      tid = re.sub(r'\.[0-9]{1,}$','',vep_info_tags['Feature'][alt_allele]) ## ignore Ensembl transcript version when annotating gene annotations (should be similar across transcript versions)
      if transcript_xref.has_key(tid):
         transcript_gene_annotations = transcript_xref[tid]
         for tanno in transcript_gene_annotations:
            for ann in pcgr_gene_annotations:
               if tanno[ann] != 'NA':
                  gene_values[ann][re.sub(r' ','_',tanno[ann])] = 1
      else:
         if tid.startswith('NM_') or tid.startswith('ENST'):
            logger.info("Could not find gene cross-reference information for transcript " + str(tid))
   
      for ann in ['ENTREZ_ID','APPRIS','CANCER_CENSUS_SOMATIC','CANCER_CENSUS_GERMLINE','INTOGEN_DRIVER','ANTINEOPLASTIC_DRUG_INTERACTION','TUMOR_SUPPRESSOR','ONCOGENE','ONCOSCORE']:
         pcgr_gene_info_tags[ann] = {}
         pcgr_gene_info_tags[ann][alt_allele] = {}
      
      for v in gene_values['entrezgene'].keys():
         pcgr_gene_info_tags['ENTREZ_ID'][alt_allele][v] = 1
      
      for v in gene_values['principal_isoform_flag'].keys():
         pcgr_gene_info_tags['APPRIS'][alt_allele][v] = 1
      
      for v in gene_values['cancer_census_somatic'].keys():
         pcgr_gene_info_tags['CANCER_CENSUS_SOMATIC'][alt_allele][v] = 1
      
      for v in gene_values['cancer_census_germline'].keys():
         pcgr_gene_info_tags['CANCER_CENSUS_GERMLINE'][alt_allele][v] = 1
      
      for v in gene_values['intogen_drivers'].keys():
         pcgr_gene_info_tags['INTOGEN_DRIVER'][alt_allele][v] = 1
         
      for v in gene_values['tsgene'].keys():
         pcgr_gene_info_tags['TUMOR_SUPPRESSOR'][alt_allele][v] = 1
      
      for v in gene_values['ts_oncogene'].keys():
         pcgr_gene_info_tags['ONCOGENE'][alt_allele][v] = 1
         
      for v in gene_values['antineoplastic_drugs_dgidb'].keys():
         pcgr_gene_info_tags['ANTINEOPLASTIC_DRUG_INTERACTION'][alt_allele][v] = 1
      
      for v in gene_values['oncoscore'].keys():
         pcgr_gene_info_tags['ONCOSCORE'][alt_allele][v] = 1
  

def extend_vcf_annotations(query_vcf, pcgr_directory):
   """
   Function that reads VEP/vcfanno-annotated VCF and extends the VCF INFO column with tags from
   1. CSQ elements within the primary transcript consequence picked by VEP, e.g. SYMBOL, Feature, Gene, Consequence etc.
   2. Cancer-relevant gene annotations, e.g. known oncogenes/tumor suppressors, known antineoplastic drugs interacting with a given protein etc.
   3. Protein-relevant annotations, e.g. cancer hotspot mutations, functional protein features etc.
   4. Variant effect predictions
   """
   
   ## read annotation datasets
   pfam_domain_names = pcgrutils.index_pfam_names(os.path.join(pcgr_directory,'data','pfam','pfam.domains.tsv.gz'), ignore_versions = True)
   uniprot_feature_names = pcgrutils.index_uniprot_feature_names(os.path.join(pcgr_directory,'data','uniprot','uniprot.features.tsv.gz'))
   pfam_xref = pcgrutils.index_pfam(os.path.join(pcgr_directory,'data','pfam','pfam.uniprot.tsv.gz'))
   uniprot_feature_xref = pcgrutils.index_uniprot_features(os.path.join(pcgr_directory,'data','uniprot','uniprot.features.tsv.gz'))
   gene_xref = pcgrutils.index_gene_transcripts(os.path.join(pcgr_directory,'data','gene.transcript.onco_xref.GRCh37.tsv.gz'), index = 'ensGene_transcript')
   up_xref = pcgrutils.index_uniprot(os.path.join(pcgr_directory, 'data','uniprot','uniprot.xref.tsv.gz'), index = 'ensGene_transcript')
   cancer_hotspot_xref = pcgrutils.index_cancer_hotspots(os.path.join(pcgr_directory,'data','cancerhotspots.org','cancer_hotspots.tsv'))

   ##output VCF fname
   out_prefix = re.sub(r'\.vcf(\.gz){0,}$','.annotated.vcf',query_vcf)
  
   ## read VEP and PCGR tags to be appended to VCF file
   vep_infotags_desc = pcgrutils.read_infotag_file(os.path.join(pcgr_directory,'data','vep_infotags.tsv'))
   pcgr_infotags_desc = pcgrutils.read_infotag_file(os.path.join(pcgr_directory,'data','pcgr_infotags.tsv'))

   ## get the annotation elements of the CSQ tag provided by VEP and create two dictionaries for lookup
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
   
   ## retrieve dbnsfp algorithm predictions from header
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

   ## set VCF header information
   for tag in vep_infotags_desc:
      vcf_reader.infos[tag] = [tag,str(vep_infotags_desc[tag]['number']),str(vep_infotags_desc[tag]['type']),str(vep_infotags_desc[tag]['description'])]
   for tag in pcgr_infotags_desc:
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
      pcgr_gene_info_tags = {}
      pcgr_protein_info_tags = {}
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
               if csq_fields[vep_csq_fields2index['PICK']] == "1": ## only consider the primary/picked consequence when expanding with annotation tags
                  num_picks += 1
                  j = 0
                  ## loop over all CSQ elements and set them in the vep_info_tags dictionary (for each alt_allele)
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
         set_pcgr_protein_info_tags(up_xref, uniprot_feature_xref, cancer_hotspot_xref, uniprot_feature_names,vep_info_tags,pcgr_protein_info_tags)
      set_pcgr_gene_info_tags(gene_xref, vep_info_tags, pcgr_gene_info_tags)
      if rec.INFO.has_key('DBNSFP'):
         map_variant_effect_predictors(rec, vep_info_tags, pcgr_protein_info_tags, variant_effect_prediction_tags, dbnsfp_prediction_algorithms)
      all_info_vals = []
      
      ## set info tags in the current variant record from
      ## 1. vep_info tags - All CSQ annotations
      ## 2. variant_effect_prediction_tags - dbNSFP annotations
      ## 3. existing_info_tags - annotations already present in query VCF
      ## 4. pcgr_gene_info_tags - cancer-relevant gene annotations
      ## 5. pcgr_protein_info_tags - protein-relevant gene annotations (hotspots, features etc.)
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
         
      
      for k in pcgr_gene_info_tags.keys():
         values = []
         for alt_allele in sorted(pcgr_gene_info_tags[k].keys()):
            if len(pcgr_gene_info_tags[k][alt_allele].keys()) != 0:
               values.append('&'.join(pcgr_gene_info_tags[k][alt_allele].keys()))
         if(len(values) > 0):
            if k == 'ONCOGENE' or k == 'TUMOR_SUPPRESSOR':
               all_info_vals.append(str(k))
            else:
               all_info_vals.append(str(k) + '=' + ','.join(values))
         
      for k in pcgr_protein_info_tags.keys():
         values = []
         if k == 'PROTEIN_POSITIONS':
            continue
         for alt_allele in sorted(pcgr_protein_info_tags[k].keys()):
            if k == 'UNIPROT_FEATURE':
               if len(pcgr_protein_info_tags[k][alt_allele].keys()) != 0:
                  values.append('&'.join(pcgr_protein_info_tags[k][alt_allele].keys()))
            else:
               if pcgr_protein_info_tags[k][alt_allele] != '.':
                  values.append(pcgr_protein_info_tags[k][alt_allele])
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


      
