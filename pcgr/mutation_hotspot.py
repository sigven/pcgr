#!/usr/bin/env python

import os,re
import csv
import gzip
from pcgr import annoutils
from pcgr.annoutils import threeToOneAA

from typing import Dict
from logging import Logger

def load_mutation_hotspots(hotspots_fname: str, logger: Logger) -> Dict[str, Dict[str, Dict[str, str]]]:
    """
    Load mutation hotspots from a file and create a dictionary of hotspots.

    Parameters:
        hotspots_fname (str): The file path of the mutation hotspots file.
        logger (Logger): The logger object to log messages.

    Returns:
        Dict[str, Dict[str, Dict[str, str]]]: A dictionary containing mutation hotspots categorized by mutation type.
            The dictionary has three keys: 'mutation', 'codon', and 'splice'.
            Each key maps to a sub-dictionary that contains the hotspots for that mutation type.
            The sub-dictionaries are indexed by a combination of gene and mutation identifier.

    Raises:
        SystemExit: If the mutation hotspots file does not exist.

    """
    hotspots: Dict[str, Dict[str, Dict[str, str]]] = {
        'mutation': {},
        'codon': {},
        'splice': {}
    }

    if not os.path.exists(hotspots_fname):
        logger.info(f"ERROR: File '{hotspots_fname}' does not exist - exiting")
        exit(1)

    with gzip.open(hotspots_fname, mode='rt') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            gene = str(row['entrezgene'])
            hgvsp2 = row['hgvsp2']
            codon = row['codon']
            hgvsc = row['hgvsc']

            hotspots['mutation'][gene + '-' + hgvsp2] = row
            hotspots['codon'][gene + '-' + codon] = row
            if hgvsc != '.':
                hotspots['splice'][gene + '-' + hgvsc] = row

    return hotspots


def match_csq_mutation_hotspot(transcript_csq_elements, cancer_hotspots, rec, principal_csq_properties):

   """
   Function that matches consequence entries from VEP (transcript_csq_elements) with entries in mutation hotspots
   """

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
