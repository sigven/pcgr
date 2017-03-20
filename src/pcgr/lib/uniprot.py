#!/usr/bin/env python

import os,re,sys
import logging
import csv
import gzip
from bx.intervals.intersection import IntervalTree

csv.field_size_limit(500 * 1024 * 1024)


def index_uniprot(uniprot_file_path, index = 'refGene_transcript'):
   
   uniprot_xref = {}
   with gzip.open(uniprot_file_path, 'rb') as tsvfile:
      upreader = csv.DictReader(tsvfile, delimiter='\t', quotechar='|')
      for rec in upreader:
         if index == 'refGene_transcript' and not (rec['transcript_id'].startswith('NM_') or rec['transcript_id'].startswith('NR_')):
            continue
         if index == 'ensGene_transcript' and not rec['transcript_id'].startswith('ENST'):
            continue
         if not uniprot_xref.has_key(rec['transcript_id']):
            uniprot_xref[rec['transcript_id']] = []
         uniprot_xref[rec['transcript_id']].append(rec)
   return uniprot_xref

def index_pfam(pfam_file_path):
 
   # dictionary mapping uniprot_ids to interval trees
   pfam = dict()
 
   # parse the UniProt-PFAM annotations file (tsv) and build the interval trees
   with gzip.open(pfam_file_path, 'r') as annotations_file:
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


def index_uniprot_features(uniprot_feature_file_path):
   uniprot_features = dict()
 
 
   observed_feats = {}
   # parse the UniProt-PFAM annotations file (tsv) and build the interval trees
   with gzip.open(uniprot_feature_file_path, 'rb') as annotations_file:
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
					

def index_pfam_names(pfam_domains_file_path, ignore_versions = False):
   
   pfam_domain_names = {}
   with gzip.open(pfam_domains_file_path, 'rb') as tsvfile:
      pfam_reader = csv.DictReader(tsvfile, delimiter='\t', quotechar='|')
      for rec in pfam_reader:
         if ignore_versions is True:
            pfam_domain_names[re.sub(r'\.[0-9]{1,}$','',rec['pfam_id'])] = rec['name']
         else:
            pfam_domain_names[rec['pfam_id']] = rec['name']
   return pfam_domain_names


def index_uniprot_feature_names(sp_features_file_path):
   
   swissprot_features = {}
   with gzip.open(sp_features_file_path, 'rb') as tsvfile:
      sp_reader = csv.DictReader(tsvfile, delimiter='\t', quotechar='|')
      for rec in sp_reader:
         sp_key = rec['uniprot_id'] + ':' + re.sub(r'#',':',rec['key'])
         swissprot_features[sp_key] = rec
   return swissprot_features


