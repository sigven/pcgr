#!/usr/bin/env python

import os,re,sys
import csv
import gzip
import annoutils

def load_biomarkers(logger, biomarker_variant_fname, biomarker_clinical_fname):

   ## load actionable cancer variants from 'variant.tsv.gz' (provided by github.com/sigven/cancerHotspots)

   variant_biomarkers = {} ##dictionary to return
   for variant_alias_type in ['dbsnp','hgvsp','hgvsc','genomic','exon','other']:
      variant_biomarkers[variant_alias_type] = {}

   if not os.path.exists(biomarker_clinical_fname):
      logger.info("ERROR: File '" + str(biomarker_clinical_fname) + "' does not exist - exiting")
      exit(1)
   

   variant_to_clinical_evidence = {}
   with gzip.open(biomarker_clinical_fname, mode='rt') as f:
      reader = csv.DictReader(f, delimiter='\t')
      for row in reader:
         if not row['variant_id'] in variant_to_clinical_evidence:
            variant_to_clinical_evidence[row['variant_id']] = {}
         variant_to_clinical_evidence[row['variant_id']][row['evidence_id']] = 1

   f.close()

   for vid in variant_to_clinical_evidence.keys():
      variant_to_clinical_evidence[vid] = '&'.join(sorted(variant_to_clinical_evidence[vid].keys()))

   if not os.path.exists(biomarker_variant_fname):
      logger.info("ERROR: File '" + str(biomarker_variant_fname) + "' does not exist - exiting")
      exit(1)
   
   with gzip.open(biomarker_variant_fname, mode='rt') as f:
      reader = csv.DictReader(f, delimiter='\t')
      for row in reader:
         if 'alteration_type' in row.keys() and 'alias_type' in row.keys() and 'variant_exon' in row.keys() and 'entrezgene' in row.keys() and 'variant_id' in row.keys():
            if row['alteration_type'].startswith('MUT') or row['alteration_type'].startswith('CODON'):
               row['clinical_evidence_items_id'] = '.'
               if row['variant_id'] in variant_to_clinical_evidence.keys():
                  row['clinical_evidence_items_id'] = variant_to_clinical_evidence[row['variant_id']]
               entry_alias_type = str(row['alias_type']).replace("_grch37", "")
               entry_alias_type = entry_alias_type.replace("_grch38", "")

               if entry_alias_type == "other":
                  if bool(re.search(r'^((ACTIVATING )?MUTATION|LOSS|START LOSS)$')) is True:
                     varkey = str(row['entrezgene'])
                     if not varkey in variant_biomarkers['other']:
                        variant_biomarkers['other'][varkey] = []
                     variant_biomarkers['other'][varkey].append(row)

               if entry_alias_type == 'exon':
                  exons = row['variant_exon']
                  if not bool(re.search(r'&|-', exons)) is True:
                     for exon in str(exons).split('|'):
                        varkey = str(row['entrezgene']) + '_' + str(exon)
                        if not varkey in variant_biomarkers['exon']:
                           variant_biomarkers['exon'][varkey] = []
                        variant_biomarkers['exon'][varkey].append(row)

               for alias_type in ['dbsnp','hgvsp','hgvsc','genomic']:                 
                  if entry_alias_type == alias_type:
                     varkey = str(row['entrezgene']) + '_' + str(row['variant_alias'])
                     if not varkey in variant_biomarkers[alias_type]:
                        variant_biomarkers[alias_type][varkey] = []
                     variant_biomarkers[alias_type][varkey].append(row)

               

   f.close()

   return variant_biomarkers


def match_csq_biomarker(transcript_csq_elements, variant_biomarkers, rec, principal_hgvsp, principal_hgvsc):

   
   principal_codon = '.'
   if re.match(r'^(p.[A-Z]{1}[0-9]{1,}[A-Za-z]{1,})', principal_hgvsp):
      codon_match = re.findall(r'[A-Z][0-9]{1,}', principal_hgvsp)
      if len(codon_match) == 1:
         principal_codon = codon_match[0]
   
   ## loop through all transcript-specific consequences ('csq_elements') for a given variant, check for the presence of
   ## 1. Exonic, protein-altering mutations (dictionary key = entrezgene + hgvsp) that overlap biomarker entries (alias_type = "hgvsp")
   ##    - assert whether a potentially hit reflects the principal hgvsp ('by_hgvsp_principal') or an alternative hgvsp ('by_hgvsp_nonprincipal')
   ## 2. Splice site mutations (dictionary key = entrezgene + hgvsc) that overlap known cancer hotspots (NOTE: only a subset have been curated)
   ##    - assert whether a potentially hit reflects the principal hgvsc ('by_hgvsc_principal') or an alternative hgvsc ('by_hgvsc_nonprincipal')


   #genomic_var_id = f"g.{rec.CHROM}:{pos}{rec.REF}>{alt_allele}"
   genomic_var_key = f"{rec.CHROM}_{rec.POS}_{rec.REF}_{rec.ALT}"
   if genomic_var_key in variant_biomarkers['genomic']:
      biomarker_data = variant_biomarkers['genomic'][genomic_var_key]

   for csq in transcript_csq_elements:
      (consequence, symbol, entrezgene, hgvsc, hgvsp, feature_type, feature, biotype) = csq.split(':')

      if not bool(re.search(r'^(missense|stop|start|inframe|splice_donor|splice_acceptor|frameshift)', consequence)) is True:
         continue

      hgvsp_short = annoutils.threeToOneAA(hgvsp)
      biomarker_key = "."
      codon_match = []
      if entrezgene != "." and hgvsp != ".":
         biomarker_key = str(entrezgene) + '_' + str(hgvsp_short)
         codon_match = re.findall(r'p.[A-Z][0-9]{1,}',hgvsp_short)

      if entrezgene != "." and (consequence == 'splice_donor_variant' or consequence == 'splice_acceptor_variant'):
         hgvsc_key = re.sub(r'>(A|G|C|T)$','',hgvsc)
         biomarker_key = str(entrezgene) + '-' + str(hgvsc_key)

      if biomarker_key == ".":
         continue

      if biomarker_key in variant_biomarkers['hgvsp']:
         hotspot_info = variant_biomarkers['hgvsp'][biomarker_key]
         #hotspot_info_ttype = cancer_hotspots['mutation'][hotspot_key_mutation]['MUTATION_HOTSPOT_CANCERTYPE']
         #unique_hotspot_mutations['exonic|' + str(hotspot_info)] = hotspot_info_ttype

      # if biomarker_key in variant_biomarkers['splice']:
      #    hotspot_info = cancer_hotspots['splice'][hotspot_key_mutation]['MUTATION_HOTSPOT2']
      #    hotspot_info_ttype = cancer_hotspots['splice'][hotspot_key_mutation]['MUTATION_HOTSPOT_CANCERTYPE']
      #    unique_hotspot_mutations['splice|' + str(hotspot_info)] = hotspot_info_ttype

               
      # if len(codon_match) > 0:
      #    hotspot_key_codon = str(entrezgene) + '-' + str(codon_match[0])

      #    if hotspot_key_codon in cancer_hotspots['codon']:            
      #       unique_hotspot_codons[str('exonic|') + cancer_hotspots['codon'][hotspot_key_codon]['MUTATION_HOTSPOT2']] =  \
      #          cancer_hotspots['codon'][hotspot_key_codon]['MUTATION_HOTSPOT_CANCERTYPE']
         
   # if len(unique_hotspot_mutations.keys()) > 0:
   #    if len(unique_hotspot_mutations.keys()) == 1:
   #       for gene_mutation_key in unique_hotspot_mutations.keys():
   #          rec.INFO['MUTATION_HOTSPOT'] = gene_mutation_key
   #          rec.INFO['MUTATION_HOTSPOT_CANCERTYPE'] = unique_hotspot_mutations[gene_mutation_key]

   #          if gene_mutation_key.split('|')[0] == 'exonic':
   #             rec.INFO['MUTATION_HOTSPOT_MATCH'] = 'by_hgvsp_nonprincipal'
   #             hgvsp_candidate = 'p.' + str(gene_mutation_key.split('|')[3]) + str(gene_mutation_key.split('|')[4])
   #             if hgvsp_candidate == principal_hgvsp:
   #                rec.INFO['MUTATION_HOTSPOT_MATCH'] = 'by_hgvsp_principal'
   #          else:
   #             rec.INFO['MUTATION_HOTSPOT_MATCH'] = 'by_hgvsc_nonprincipal'
   #             hgvsc_candidate = re.sub(r'>(A|G|C|T){1,}$', '' , str(gene_mutation_key.split('|')[4]))
   #             if hgvsc_candidate == principal_hgvsc:
   #                rec.INFO['MUTATION_HOTSPOT_MATCH'] = 'by_hgvsc_principal'
   #    else:
   #       ## multiple hotspot matches for alternative hgvsp keys
   #       ## - keep only match against principal HGVSp 
   #       for hotspot_info in unique_hotspot_mutations.keys():
   #          if hotspot_info.split('|')[0] == 'exonic':
   #             hgvsp_candidate = "p." + str(hotspot_info.split('|')[3]) + str(hotspot_info.split('|')[4]) 

   #             if hgvsp_candidate == principal_hgvsp:
   #                rec.INFO['MUTATION_HOTSPOT'] = hotspot_info
   #                rec.INFO['MUTATION_HOTSPOT_CANCERTYPE'] = unique_hotspot_mutations[hotspot_info]
   #                rec.INFO['MUTATION_HOTSPOT_MATCH'] = 'by_hgvsp_principal'
   #          else:
   #             hgvsc_candidate = re.sub(r'>(A|G|C|T){1,}$', '' , str(hotspot_info.split('|')[4]))

   #             if hgvsc_candidate == principal_hgvsc:
   #                rec.INFO['MUTATION_HOTSPOT'] = hotspot_info
   #                rec.INFO['MUTATION_HOTSPOT_CANCERTYPE'] = unique_hotspot_mutations[hotspot_info]
   #                rec.INFO['MUTATION_HOTSPOT_MATCH'] = 'by_hgvsc_principal'

   # else:
   #    if len(unique_hotspot_codons.keys()) > 0:
   #       if len(unique_hotspot_codons.keys()) == 1:
   #          for gene_codon_key in unique_hotspot_codons.keys():

   #             if '|' in gene_codon_key:

   #                codon = str(gene_codon_key.split('|')[3])

   #                if codon == principal_codon:
   #                   rec.INFO['MUTATION_HOTSPOT'] = gene_codon_key
   #                   rec.INFO['MUTATION_HOTSPOT_CANCERTYPE'] = unique_hotspot_codons[gene_codon_key]
   #                   rec.INFO['MUTATION_HOTSPOT_MATCH'] = 'by_codon_principal'
   #                else:
   #                   rec.INFO['MUTATION_HOTSPOT'] = gene_codon_key
   #                   rec.INFO['MUTATION_HOTSPOT_CANCERTYPE'] = unique_hotspot_codons[gene_codon_key]
   #                   rec.INFO['MUTATION_HOTSPOT_MATCH'] = 'by_codon_nonprincipal'


   return
