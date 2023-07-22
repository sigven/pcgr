#!/usr/bin/env python

import os,re
import csv
import gzip
from pcgr import annoutils
from pcgr.annoutils import threeToOneAA

def load_biomarkers(logger, biomarker_variant_fname, biomarker_clinical_fname):

   ## load actionable cancer variants from 'variant.tsv.gz' (provided by github.com/sigven/cancerHotspots)

   variant_biomarkers = {} ##dictionary to return
   for variant_alias_type in ['dbsnp','hgvsp','hgvsc','genomic','exon','other']:
      variant_biomarkers[variant_alias_type] = {}

   if not os.path.exists(biomarker_clinical_fname):
      logger.info("ERROR: File '" + str(biomarker_clinical_fname) + "' does not exist - exiting")
      exit(1)
   

   variant_to_clinical_evidence = {}
   variant_to_origin = {}
   with gzip.open(biomarker_clinical_fname, mode='rt') as f:
      reader = csv.DictReader(f, delimiter='\t')
      for row in reader:
         if not row['variant_id'] in variant_to_clinical_evidence:
            variant_to_clinical_evidence[row['variant_id']] = {}
         variant_to_clinical_evidence[row['variant_id']][row['evidence_id']] = 1
         if not row['variant_id'] in variant_to_origin:
            variant_to_origin[row['variant_id']] = {}
         variant_to_origin[row['variant_id']][row['variant_origin']] = 1

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
               if row['variant_id'] in variant_to_origin.keys():
                  row['variant_origin'] = '&'.join(variant_to_origin[row['variant_id']].keys())
               entry_alias_type = str(row['alias_type']).replace("_grch37", "")
               entry_alias_type = entry_alias_type.replace("_grch38", "")

               if entry_alias_type == "other":
                  if bool(re.search(r'^((ACTIVATING )?MUTATION|LOSS|START LOSS)$', row['variant_alias'])) is True:
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


def match_csq_biomarker(transcript_csq_elements, variant_biomarkers, rec, principal_csq_properties):

   ## Picked CSQ by VEP
   principal_hgvsp = principal_csq_properties['hgvsp']
   principal_hgvsc = principal_csq_properties['hgvsc']
   principal_entrezgene = principal_csq_properties['entrezgene']
   principal_codon = principal_csq_properties['codon'] 
   principal_exon = principal_csq_properties['exon']   

   biomarker_hits_all = {}
   
   ## Match variant by genomic coordinate
   genomic_var_key = f"{principal_entrezgene}_{rec.CHROM}_{rec.POS}_{rec.REF}_{rec.ALT[0]}"

   if genomic_var_key in variant_biomarkers['genomic'].keys():
      hits_genomic = variant_biomarkers['genomic'][genomic_var_key]

      #print("BALLE")

      for ghit in hits_genomic:
         bkey_genomic = f"{ghit['biomarker_source']}|{ghit['variant_origin']}|{ghit['variant_id']}|{ghit['clinical_evidence_items_id']}"
         if not bkey_genomic in biomarker_hits_all.keys():            
            biomarker_hits_all[bkey_genomic] = {}
         biomarker_hits_all[bkey_genomic]['by_genomic_coord'] = 1
         

   mut_lof_fs = False
   mut_lof = False
   mut_protein = False
   for csq_elem in transcript_csq_elements:
      (consequence, symbol, entrezgene, hgvsc, hgvsp, exon, feature_type, feature, biotype) = csq_elem.split(':')

      if bool(re.search(r'^(missense|stop|start|inframe|splice_donor|protein|splice_acceptor|frameshift)', consequence)) is True:
         mut_protein = True
      if bool(re.search(r'^(stop_gain|splice_donor_variant|splice_acceptor_variant|frameshift)', consequence)) is True:
         mut_lof = True
      if bool(re.search(r'^(frameshift)', consequence)) is True:
         mut_lof_fs = True


      principal_csq_hgvsp = False
      principal_csq_entrezgene = False
      principal_csq_hgvsc = False

      hgvsp_short = threeToOneAA(hgvsp)

      if hgvsp_short == principal_hgvsp:
         principal_csq_hgvsp = True
      if entrezgene == principal_entrezgene:
         principal_csq_entrezgene = True
      if hgvsc == principal_hgvsc:
         principal_csq_hgvsc = True

      codon_match = []
      if entrezgene != "." and hgvsp != ".":
         biomarker_key_hgvsp = str(entrezgene) + '_' + str(hgvsp_short)
         codon_match = re.findall(r'p.[A-Z][0-9]{1,}',hgvsp_short)

         print(biomarker_key_hgvsp)
         ## match biomarkers annotated at the amino acid level - with HGVSp coordinates
         if biomarker_key_hgvsp in variant_biomarkers['hgvsp']:
            hits = variant_biomarkers['hgvsp'][biomarker_key_hgvsp]
            for h in hits:
               bkey = f"{h['biomarker_source']}|{h['variant_origin']}|{h['variant_id']}|{h['clinical_evidence_items_id']}"
               if not bkey in biomarker_hits_all.keys():            
                  biomarker_hits_all[bkey] = {}
               if principal_csq_entrezgene is True:
                  ## principal hgvsp (VEP's picked csq)
                  if principal_csq_hgvsp is True:
                     biomarker_hits_all[bkey]['by_hgvsp_principal'] = 1
                  else:
                     ## nonprincipal hgvsp
                     biomarker_hits_all[bkey]['by_hgvsp_nonprincipal'] = 1
         
         if len(codon_match) > 0:
            biomarker_key_codon = str(entrezgene) + '_' + str(codon_match[0])

            ## match biomarkers annotated as "CODON" only for a given gene
            if biomarker_key_codon in variant_biomarkers['hgvsp']:
               hits_codon = variant_biomarkers['hgvsp'][biomarker_key_codon]
               for chit in hits_codon:
                  bkey2 = f"{chit['biomarker_source']}|{chit['variant_origin']}|{chit['variant_id']}|{chit['clinical_evidence_items_id']}"
                  if not bkey2 in biomarker_hits_all.keys():            
                     biomarker_hits_all[bkey2] = {}
                  if principal_csq_entrezgene is True:
                     ## principal codon (VEP's picked csq)
                     if codon_match[0] == principal_codon:
                        biomarker_hits_all[bkey2]['by_codon_principal'] = 1
                     else:
                        ## nonprincipal codon
                        biomarker_hits_all[bkey2]['by_codon_nonprincipal'] = 1
         

      if entrezgene != "." and (consequence == 'splice_donor_variant' or consequence == 'splice_acceptor_variant'):
         hgvsc_key = re.sub(r'>(A|G|C|T)$','',hgvsc)
         hgvsc_biomarker_key = str(entrezgene) + '_' + str(hgvsc_key)

      
      if entrezgene != "." and principal_csq_entrezgene is True and exon != ".":
         exon_biomarker_key = str(entrezgene) + '_' + str(exon)
         if exon_biomarker_key in variant_biomarkers['exon'].keys():
            hits_exon = variant_biomarkers['exon'][exon_biomarker_key]

            for ehit in hits_exon:
               bkey4 = f"{ehit['biomarker_source']}|{ehit['variant_origin']}|{ehit['variant_id']}|{ehit['clinical_evidence_items_id']}"
               if ehit['variant_consequence'] == "exon_variant" and re.search(r'MUTATION', ehit['variant_alias']) and re.search(r'^(missense|coding|protein)', consequence):
                  if not bkey4 in biomarker_hits_all.keys():            
                     biomarker_hits_all[bkey4] = {}
                  if exon == principal_exon:
                     biomarker_hits_all[bkey4]['by_exon_mut_principal'] = 1
                  else:
                     biomarker_hits_all[bkey4]['by_exon_mut_nonprincipal'] = 1
               
               elif ehit['variant_consequence'] == "inframe_deletion" and re.search(r'^inframe_deletion', consequence):                  
                  if not bkey4 in biomarker_hits_all.keys():            
                     biomarker_hits_all[bkey4] = {}
                  if exon == principal_exon:
                     biomarker_hits_all[bkey4]['by_exon_deletion_principal'] = 1
                  else:
                     biomarker_hits_all[bkey4]['by_exon_deletion_nonprincipal'] = 1
               
               else:
                  if ehit['variant_consequence'] == "inframe_insertion" and re.search(r'^inframe_insertion', consequence):
                     if not bkey4 in biomarker_hits_all.keys():            
                        biomarker_hits_all[bkey4] = {}
                     if exon == principal_exon:
                        biomarker_hits_all[bkey4]['by_exon_insertion_principal'] = 1
                     else:
                        biomarker_hits_all[bkey4]['by_exon_insertion_nonprincipal'] = 1
               


      if entrezgene != "." and principal_csq_entrezgene is True:
         if str(entrezgene) in variant_biomarkers['other'].keys():
            hits_gene = variant_biomarkers['other'][str(entrezgene)]
            for ghit in hits_gene:
               bkey3 = f"{ghit['biomarker_source']}|{ghit['variant_origin']}|{ghit['variant_id']}|{ghit['clinical_evidence_items_id']}"
               ## match biomarkers annotated as "Mutation" only for a given gene - consider only the principal consequence (VEP's picked)
               if ghit['alteration_type'] == 'MUT' and mut_protein is True and (principal_csq_hgvsc is True or principal_csq_hgvsp is True):
                  if not bkey3 in biomarker_hits_all.keys():            
                     biomarker_hits_all[bkey3] = {}
                  biomarker_hits_all[bkey3]['by_gene_mut'] = 1
               
               ## match biomarkers annotated as "Mutation - loss of function" only for a given gene - consider only the principal consequence (VEP's picked)
               elif ghit['alteration_type'] == 'MUT_LOF' and mut_lof is True and (principal_csq_hgvsc is True or principal_csq_hgvsp is True):
                  if not bkey3 in biomarker_hits_all.keys():            
                     biomarker_hits_all[bkey3] = {}
                  biomarker_hits_all[bkey3]['by_gene_mut_lof'] = 1
               
               ## match biomarkers annotated as "Mutation - loss of function - frameshift" only for a given gene - consider only the principal consequence (VEP's picked)
               else:
                  if ghit['alteration_type'] == 'MUT_LOF_FS' and mut_lof_fs is True and (principal_csq_hgvsc is True or principal_csq_hgvsp is True):
                     if not bkey3 in biomarker_hits_all.keys():            
                        biomarker_hits_all[bkey3] = {}
                     biomarker_hits_all[bkey3]['by_gene_mut_lof_fs'] = 1

   
   if len(biomarker_hits_all.keys()) > 0:
      biomarker_var_matches = {}
      for bm in biomarker_hits_all.keys():
         match_types = '&'.join(biomarker_hits_all[bm].keys())
         vkey = f"{bm}|{match_types}"
         biomarker_var_matches[vkey] = 1

      rec.INFO['BIOMARKER_MATCH'] = str(','.join(sorted(biomarker_var_matches.keys())))

   return
