#!/usr/bin/env python

import os,re
import csv
import gzip
from pcgr import annoutils
from pcgr.annoutils import threeToOneAA

def load_mutation_hotspots(hotspots_fname, logger):

   ## load cancer hotspot entries from 'misc/tsv/hotspot/hotspot.tsv.gz' (provided by github.com/sigven/cancerHotspots)

   hotspots = {} ##dictionary to return
   mutation_hotspots = {}
   codon_hotspots = {}
   splice_hotspots = {}
   if not os.path.exists(hotspots_fname):
      logger.info("ERROR: File '" + str(hotspots_fname) + "' does not exist - exiting")
      exit(1)
   
   with gzip.open(hotspots_fname, mode='rt') as f:
      reader = csv.DictReader(f, delimiter='\t')
      for row in reader:
         mutation_hotspots[str(row['entrezgene']) + '-' + row['hgvsp2']] = row
         codon_hotspots[str(row['entrezgene']) + '-' + row['codon']] = row
         if row['hgvsc'] != '.':
            splice_hotspots[str(row['entrezgene']) + '-' + str(row['hgvsc'])] = row
   f.close()

   hotspots['mutation'] = mutation_hotspots
   hotspots['codon'] = codon_hotspots
   hotspots['splice'] = splice_hotspots
   return hotspots


def match_csq_mutation_hotspot(transcript_csq_elements, cancer_hotspots, rec, principal_csq_properties):

   unique_hotspot_mutations = {}
   unique_hotspot_codons = {}

   principal_hgvsp = principal_csq_properties['hgvsp']
   principal_hgvsc = principal_csq_properties['hgvsc']
   principal_entrezgene = principal_csq_properties['entrezgene']
   principal_codon = principal_csq_properties['codon']   
   
   ## loop through all transcript-specific consequences ('csq_elements') for a given variant, check for the presence of
   ## 1. Exonic, protein-altering mutations (dictionary key = entrezgene + hgvsp) that overlap known cancer hotspots (https://github.com/sigven/cancerHotspots)
   ##    - assert whether a potentially hit reflects the principal hgvsp ('by_hgvsp_principal') or an alternative hgvsp ('by_hgvsp_nonprincipal')
   ## 2. Splice site mutations (dictionary key = entrezgene + hgvsc) that overlap known cancer hotspots (NOTE: only a subset have been curated)
   ##    - assert whether a potentially hit reflects the principal hgvsc ('by_hgvsc_principal') or an alternative hgvsc ('by_hgvsc_nonprincipal')

   for csq in transcript_csq_elements:
      (consequence, symbol, entrezgene, hgvsc, hgvsp, exon, feature_type, feature, biotype) = csq.split(':')

      if not bool(re.search(r'^(missense|stop|start|inframe|splice_donor|splice_acceptor|frameshift)', consequence)) is True:
         continue

      hgvsp_short = threeToOneAA(hgvsp)
      hotspot_key_mutation = "."
      codon_match = []
      if entrezgene != "." and hgvsp != ".":
         hotspot_key_mutation = str(entrezgene) + '-' + str(hgvsp_short)
         codon_match = re.findall(r'p.[A-Z][0-9]{1,}',hgvsp_short)

      if entrezgene != "." and (consequence == 'splice_donor_variant' or consequence == 'splice_acceptor_variant'):
         hgvsc_key = re.sub(r'>(A|G|C|T)$','',hgvsc)
         hotspot_key_mutation = str(entrezgene) + '-' + str(hgvsc_key)

      if hotspot_key_mutation == ".":
         continue

      if hotspot_key_mutation in cancer_hotspots['mutation']:
         hotspot_info = cancer_hotspots['mutation'][hotspot_key_mutation]['MUTATION_HOTSPOT2']
         hotspot_info_ttype = cancer_hotspots['mutation'][hotspot_key_mutation]['MUTATION_HOTSPOT_CANCERTYPE']
         unique_hotspot_mutations['exonic|' + str(hotspot_info)] = hotspot_info_ttype

      if hotspot_key_mutation in cancer_hotspots['splice']:
         hotspot_info = cancer_hotspots['splice'][hotspot_key_mutation]['MUTATION_HOTSPOT2']
         hotspot_info_ttype = cancer_hotspots['splice'][hotspot_key_mutation]['MUTATION_HOTSPOT_CANCERTYPE']
         unique_hotspot_mutations['splice|' + str(hotspot_info)] = hotspot_info_ttype

               
      if len(codon_match) > 0:
         hotspot_key_codon = str(entrezgene) + '-' + str(codon_match[0])

         if hotspot_key_codon in cancer_hotspots['codon']:            
            unique_hotspot_codons[str('exonic|') + cancer_hotspots['codon'][hotspot_key_codon]['MUTATION_HOTSPOT2']] =  \
               cancer_hotspots['codon'][hotspot_key_codon]['MUTATION_HOTSPOT_CANCERTYPE']
         
   if len(unique_hotspot_mutations.keys()) > 0:
      if len(unique_hotspot_mutations.keys()) == 1:
         for gene_mutation_key in unique_hotspot_mutations.keys():
            rec.INFO['MUTATION_HOTSPOT'] = gene_mutation_key
            rec.INFO['MUTATION_HOTSPOT_CANCERTYPE'] = unique_hotspot_mutations[gene_mutation_key]

            if gene_mutation_key.split('|')[0] == 'exonic':
               rec.INFO['MUTATION_HOTSPOT_MATCH'] = 'by_hgvsp_nonprincipal'
               hgvsp_candidate = 'p.' + str(gene_mutation_key.split('|')[3]) + str(gene_mutation_key.split('|')[4])
               if hgvsp_candidate == principal_hgvsp:
                  rec.INFO['MUTATION_HOTSPOT_MATCH'] = 'by_hgvsp_principal'
            else:
               rec.INFO['MUTATION_HOTSPOT_MATCH'] = 'by_hgvsc_nonprincipal'
               hgvsc_candidate = re.sub(r'>(A|G|C|T){1,}$', '' , str(gene_mutation_key.split('|')[4]))
               if hgvsc_candidate == principal_hgvsc:
                  rec.INFO['MUTATION_HOTSPOT_MATCH'] = 'by_hgvsc_principal'
      else:
         ## multiple hotspot matches for alternative hgvsp keys
         ## - keep only match against principal HGVSp 
         for hotspot_info in unique_hotspot_mutations.keys():
            if hotspot_info.split('|')[0] == 'exonic':
               hgvsp_candidate = "p." + str(hotspot_info.split('|')[3]) + str(hotspot_info.split('|')[4]) 

               if hgvsp_candidate == principal_hgvsp:
                  rec.INFO['MUTATION_HOTSPOT'] = hotspot_info
                  rec.INFO['MUTATION_HOTSPOT_CANCERTYPE'] = unique_hotspot_mutations[hotspot_info]
                  rec.INFO['MUTATION_HOTSPOT_MATCH'] = 'by_hgvsp_principal'
            else:
               hgvsc_candidate = re.sub(r'>(A|G|C|T){1,}$', '' , str(hotspot_info.split('|')[4]))

               if hgvsc_candidate == principal_hgvsc:
                  rec.INFO['MUTATION_HOTSPOT'] = hotspot_info
                  rec.INFO['MUTATION_HOTSPOT_CANCERTYPE'] = unique_hotspot_mutations[hotspot_info]
                  rec.INFO['MUTATION_HOTSPOT_MATCH'] = 'by_hgvsc_principal'

   else:
      if len(unique_hotspot_codons.keys()) > 0:
         if len(unique_hotspot_codons.keys()) == 1:
            for gene_codon_key in unique_hotspot_codons.keys():

               if '|' in gene_codon_key:

                  codon = str(gene_codon_key.split('|')[3])

                  if codon == principal_codon:
                     rec.INFO['MUTATION_HOTSPOT'] = gene_codon_key
                     rec.INFO['MUTATION_HOTSPOT_CANCERTYPE'] = unique_hotspot_codons[gene_codon_key]
                     rec.INFO['MUTATION_HOTSPOT_MATCH'] = 'by_codon_principal'
                  else:
                     rec.INFO['MUTATION_HOTSPOT'] = gene_codon_key
                     rec.INFO['MUTATION_HOTSPOT_CANCERTYPE'] = unique_hotspot_codons[gene_codon_key]
                     rec.INFO['MUTATION_HOTSPOT_MATCH'] = 'by_codon_nonprincipal'


   return
