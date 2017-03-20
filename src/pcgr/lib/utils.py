#!/usr/bin/env python

import os,re,sys
import csv

csv.field_size_limit(500 * 1024 * 1024)

def read_infotag_file(vcf_info_tags_tsv):
   info_tag_xref = {}
   if not os.path.exists(vcf_info_tags_tsv):
      return info_tag_xref
   with open(vcf_info_tags_tsv, 'rb') as tsvfile:
      vep_reader = csv.DictReader(tsvfile, delimiter='\t')
      for rec in vep_reader:
         if not info_tag_xref.has_key(rec['tag']):
            info_tag_xref[rec['tag']] = rec
   
   return info_tag_xref

def index_cancer_hotspots(cancer_hotspot_path):
   hotspot_xref = {}
   if not os.path.exists(cancer_hotspot_path):
      return hotspot_xref
   with open(cancer_hotspot_path, 'rb') as tsvfile:
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

def map_cancer_hotspots(cancer_hotspot_xref, vep_info_tags, protein_info_tags):
   
   for alt_allele in vep_info_tags['Feature'].keys():
      hotspot_hits = {}
      symbol = vep_info_tags['SYMBOL'][alt_allele]
      consequence = vep_info_tags['Consequence'][alt_allele]
      if protein_info_tags['PROTEIN_POSITIONS'].has_key(alt_allele):
         if 'missense_variant' in consequence or 'stop_gained' in consequence:
            if cancer_hotspot_xref.has_key(symbol):
               for codon in protein_info_tags['PROTEIN_POSITIONS'][alt_allele].keys():
                  if cancer_hotspot_xref[symbol].has_key(str(codon)):
                     cancer_hotspot_description = str(symbol) + '|' + str(cancer_hotspot_xref[symbol][str(codon)]['Codon']) + '|' + str(cancer_hotspot_xref[symbol][str(codon)]['Q-value'])
                     protein_info_tags['CANCER_MUTATION_HOTSPOT'][alt_allele] = cancer_hotspot_description
