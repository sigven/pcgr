#!/usr/bin/env python

import os,re,sys
import csv
import logging
import gzip
from bx.intervals.intersection import IntervalTree

csv.field_size_limit(500 * 1024 * 1024)

def read_infotag_file(vcf_info_tags_tsv):
   """
   Function that reads a VCF info tag file that denotes annotation tags produced by PCGR.
   An example of the VCF info tag file is the following:
   
   tag	number	type	description
   Consequence	.	String	"Impact modifier for the consequence type (picked by VEP's --flag_pick_allele option)."
   
   A dictionary is returned, with the tag as the key, and the full dictionary record as the value
   """
   info_tag_xref = {} ##dictionary returned
   if not os.path.exists(vcf_info_tags_tsv):
      return info_tag_xref
   with open(vcf_info_tags_tsv, 'rb') as tsvfile:
      reader = csv.DictReader(tsvfile, delimiter='\t')
      for rec in reader:
         if not info_tag_xref.has_key(rec['tag']):
            info_tag_xref[rec['tag']] = rec
   
   return info_tag_xref

def index_cancer_hotspots(cancer_hotspot_fname):
   """
   returns a dictionary of dictionaries, with gene symbols and codons as the respective keys
   Each entry is associated with the full dictionary record of actual hotspot annotations (Q-value, variants, tumor types etc)
   """
   hotspot_xref = {} ##dictionary returned
   
   if not os.path.exists(cancer_hotspot_fname):
      return hotspot_xref
   with open(cancer_hotspot_fname, 'rb') as tsvfile:
      ch_reader = csv.DictReader(tsvfile, delimiter='\t', quotechar='#')
      for rec in ch_reader:
         if 'splice' in rec['Codon']:
            continue
         gene = str(rec['Hugo Symbol']).upper()
         codon =  str(re.sub(r'[A-Z]','',rec['Codon']))
         if not hotspot_xref.has_key(gene):
            hotspot_xref[gene] = {}
         hotspot_xref[gene][codon] = rec
   return hotspot_xref


def index_gene_transcripts(gene_transcript_annotations_fname, index = 'ensGene_transcript'):
   """
   Function that parses the gene transcript annotation file and returns a dictionary with Ensembl transcript
   ID as key (e.g. ENST000.. ) and an array of transcripts as elements
   
   each transcript element is a dictionary record with annotation types as keys (e.g. symbol, name etc)
   """
   transcript_xref = {} ##dictionary returned
   with gzip.open(gene_transcript_annotations_fname, 'rb') as tsvfile:
      greader = csv.DictReader(tsvfile, delimiter='\t', quotechar='|')
      for rec in greader:
         if index == 'ensGene_transcript' and not rec['ensembl_transcript_id'].startswith('ENST'):
            continue
         #if not (rec['gencode_tag'] == 'NA' or rec['gencode_tag'] == 'basic'):
            #continue
         if not transcript_xref.has_key(rec['ensembl_transcript_id']):
            transcript_xref[rec['ensembl_transcript_id']] = []
         transcript_xref[rec['ensembl_transcript_id']].append(rec)
   return transcript_xref


def map_cancer_hotspots(cancer_hotspot_xref, vcfrec_vep_info_tags, vcfrec_protein_info_tags):
   """
   Function that matches the annotations of a VCF record (VEP_INFO + PROTEIN_INFO) with the cancer hotspot dictionary in order
   to return a dictionary of cancer mutation hotspots
   """
   
   for alt_allele in vcfrec_vep_info_tags['Feature'].keys():
      hotspot_hits = {}
      symbol = vcfrec_vep_info_tags['SYMBOL'][alt_allele]
      consequence = vcfrec_vep_info_tags['Consequence'][alt_allele]
      if vcfrec_protein_info_tags['PROTEIN_POSITIONS'].has_key(alt_allele):
         if 'missense_variant' in consequence or 'stop_gained' in consequence:
            if cancer_hotspot_xref.has_key(symbol):
               for codon in vcfrec_protein_info_tags['PROTEIN_POSITIONS'][alt_allele].keys():
                  if cancer_hotspot_xref[symbol].has_key(str(codon)):
                     cancer_hotspot_description = str(symbol) + '|' + str(cancer_hotspot_xref[symbol][str(codon)]['Codon']) + '|' + str(cancer_hotspot_xref[symbol][str(codon)]['Q-value'])
                     vcfrec_protein_info_tags['CANCER_MUTATION_HOTSPOT'][alt_allele] = cancer_hotspot_description


def index_uniprot(uniprot_fname, index = 'ensGene_transcript'):
   
   """
   Function that creates a dictionary of UniProt annotations using ENSEMBL gene transcripts as keys
   """
   
   uniprot_xref = {}
   with gzip.open(uniprot_fname, 'rb') as tsvfile:
      upreader = csv.DictReader(tsvfile, delimiter='\t', quotechar='|')
      for rec in upreader:
         if index == 'ensGene_transcript' and not rec['transcript_id'].startswith('ENST'):
            continue
         if not uniprot_xref.has_key(rec['transcript_id']):
            uniprot_xref[rec['transcript_id']] = []
         uniprot_xref[rec['transcript_id']].append(rec)
   return uniprot_xref

def index_pfam(pfam_fname):
 
   """
   Function that adds one interval tree per uniprot ID, individual PFAM domains and their corresponding amino acid positions are
   appended to the tree, enabling rapid lookup and annotation for gene variants/amino acid positions
   """
   # dictionary mapping uniprot ids to interval trees
   pfam = dict()
 
   # parse the UniProt-PFAM annotations file (tsv) and build the interval trees
   with gzip.open(pfam_fname, 'r') as annotations_file:
      reader = csv.DictReader(annotations_file, delimiter='\t', quotechar = '|')
      for rec in reader:
         # one interval tree per uniprot ID
         if rec['uniprot_id'] in pfam:
            tree = pfam[rec['uniprot_id']]
         else:
            # first time we've encountered this chromosome, create an interval tree
            tree = IntervalTree()
            pfam[rec['uniprot_id']] = tree
            
         # index the feature
         tree.add(int(rec['aa_start']), int(rec['aa_stop']), rec['pfam_id']) 
 
   return pfam


def index_uniprot_features(uniprot_feature_fname):
   """
   Function that adds one interval tree per uniprot ID, individual UniProt features and their corresponding amino acid positions (e.g. active sites etc) are
   appended to the tree, enabling rapid lookup and annotation for gene variants/amino acid positions
   """
   
   uniprot_features = dict()
   
   observed_feats = {}
   # parse the UniProt-PFAM annotations file (tsv) and build the interval trees
   with gzip.open(uniprot_feature_fname, 'rb') as annotations_file:
      reader = csv.DictReader(annotations_file, delimiter='\t', quotechar = '|')
      for rec in reader:
         if re.match(r'(CA_BIND|ZN_FING|DNA_BIND|NP_BIND|REGION|MOTIF|ACT_SITE|METAL|BINDING|SITE|MOD_RES|NON_STD|CARBOHYD|DISULFID|CROSSLNK|MUTAGEN)',rec['feature_type']) != None:  
            # one interval tree per uniprot ID
            if rec['uniprot_id'] in uniprot_features:
               tree = uniprot_features[rec['uniprot_id']]
            else:
               # first time we've encountered this chromosome, create an interval tree
               tree = IntervalTree()
               uniprot_features[rec['uniprot_id']] = tree
            
            feat_key = str(rec['uniprot_id']) + '_' + str(rec['key'])
            if not observed_feats.has_key(feat_key):
               # index the feature
               tree.add(int(rec['aa_start']), int(rec['aa_stop']), rec['key'].replace('#',':'))
               observed_feats[feat_key] = 1
 
   return uniprot_features 

def get_uniprot_data_by_transcript(up_xref, transcript_id, csq):
   
   uniprot_mappings = None
   uniprot_ids = {}
   symbols = {}
   seqmatches = {}
   if up_xref.has_key(transcript_id):
      uniprot_mappings = up_xref[transcript_id]
      if not uniprot_mappings is None:
         if len(uniprot_mappings) > 1:
            for k in uniprot_mappings:
               if k['symbol'] == csq['SYMBOL']:
                  uniprot_ids[k['uniprot_id']] = 1
                  seqmatches[k['uniprot_seq_match']] = 1
         else:
            uniprot_ids[uniprot_mappings[0]['uniprot_id']] = 1
            seqmatches[uniprot_mappings[0]['uniprot_seq_match']] = 1
         
         csq['UNIPROT_ID'] = '&'.join(uniprot_ids.keys())
         csq['SEQ_MATCH'] = '&'.join(seqmatches.keys())



def get_domains_features_by_aapos(xref, csq, qtype = 'domain'):
   
   """
   Function that finds protein features/protein domains based on the position given as input (i.e. csq['AA_position']) in the given protein (e.g. csq['UNIPROT_ID'])
   """
   
   if csq['UNIPROT_ID'] != '' and csq['AA_position'] != '' and csq['SEQ_MATCH'] != '.':
      if xref.has_key(csq['UNIPROT_ID']):
         start = None
         stop = None
         if '-' in csq['AA_position']:
            start, stop = csq['AA_position'].split('-')
         else:
            start = csq['AA_position']
            stop = start
      
         if not start is None and not stop is None and (csq['SEQ_MATCH'] == 'IDENTICAL_SEQUENCE' or csq['SEQ_MATCH'] == 'IDENTICAL_LENGTH'):
            if qtype == 'domain':
               csq['DOMAIN'] = '&'.join(xref[csq['UNIPROT_ID']].find(int(start), int(stop)))
            elif qtype == 'feature':
               csq['UNIPROT_FEATURE'] = '&'.join(xref[csq['UNIPROT_ID']].find(int(start), int(stop)))
               

def index_pfam_names(pfam_domains_fname, ignore_versions = False):
   """
   Function that indexes PFAM protein domain nmames, found in 'pfam_domains_fname' with the following format:
   
   pfam_id	url	name
   F02671.20	<a href="http://pfam.xfam.org/family/PF02671.20" target='_blank'>Paired amphipathic helix repeat</a>	Paired amphipathic helix repeat
   PF09810.8	<a href="http://pfam.xfam.org/family/PF09810.8" target='_blank'>Exonuclease V - a 5' deoxyribonuclease</a>	Exonuclease V - a 5' deoxyribonuclease
   
   A dictionary is returned, with pfam_id as the key and the domain name as the value
   """
   pfam_domain_names = {}
   with gzip.open(pfam_domains_fname, 'rb') as tsvfile:
      pfam_reader = csv.DictReader(tsvfile, delimiter='\t', quotechar='|')
      for rec in pfam_reader:
         if ignore_versions is True:
            pfam_domain_names[re.sub(r'\.[0-9]{1,}$','',rec['pfam_id'])] = rec['name']
         else:
            pfam_domain_names[rec['pfam_id']] = rec['name']
   return pfam_domain_names


def index_uniprot_feature_names(sp_features_fname):
   """
   Function that indexes UniProt/KB function features, found in 'sp_features_fname' with the following format:
   
   feature_type	uniprot_id	description	aa_start	aa_stop	key	type_description
   CHAIN	SERC_HUMAN	Phosphoserine aminotransferase	1	370	CHAIN#1#370	Polypeptide chain in the mature protein
   REGION	SERC_HUMAN	Pyridoxal phosphate binding	79	80	REGION#79#80	Region of interest
   
   'key' denotes the feature type and amino acid position(s) of the feature
   
   A dictionary is returned, with uniprot_id + key as the dictionary key, e.g. SERC_HUMAN:REGION:79:80
   and the full dictionary record as the value
   """
   swissprot_features = {}
   with gzip.open(sp_features_fname, 'rb') as tsvfile:
      sp_reader = csv.DictReader(tsvfile, delimiter='\t', quotechar='|')
      for rec in sp_reader:
         sp_key = rec['uniprot_id'] + ':' + re.sub(r'#',':',rec['key'])
         swissprot_features[sp_key] = rec
   return swissprot_features


def index_clinvar(clinvar_tsv_fname):
   clinvar_xref = {}
   with gzip.open(clinvar_tsv_fname, 'rb') as tsvfile:
      cv_reader = csv.DictReader(tsvfile, delimiter='\t')
      for rec in cv_reader:
         
         unique_traits = {}
         traits = ''
         traits = rec['all_traits']
         for t in traits.split(';'):
            t_lc = str(t).lower()
            unique_traits[t_lc] = 1
         origin = ''
         origin = rec['origin']
         
         traits_curated = ';'.join(unique_traits.keys())
         traits_origin = traits_curated + ' - ' + str(origin)
         
         clinvar_xref[rec['measureset_id']] = {}
         clinvar_xref[rec['measureset_id']]['phenotype_origin'] = traits_origin
         if rec['symbol'] == '-' or rec['symbol'] == 'more than 10':
            rec['symbol'] = 'NA'
         clinvar_xref[rec['measureset_id']]['genesymbol'] = rec['symbol']
         
   return clinvar_xref


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

