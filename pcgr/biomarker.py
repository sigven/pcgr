#!/usr/bin/env python

import re
import csv
import gzip
import json
import os
import urllib.request
import pandas as pd
import subprocess
import logging

#from pcgr import annoutils
from pcgr.annoutils import threeToOneAA
from pcgr.utils import check_file_exists, error_message
from pcgr import pcgr_vars


def load_all_biomarkers(
   logger: logging.Logger, 
   refdata_assembly_dir: str, 
   biomarker_vartype: str = 'MUT', 
   biomarker_variant_origin: str = 'Somatic'):
   
   ## biomarker evidence
   biomarkers = {}
   actionable_dict = {}
    
   for db in ['cgi','civic']:
      variant_fname = os.path.join(refdata_assembly_dir, 'biomarker','tsv', f"{db}.variant.tsv.gz")
      clinical_fname = os.path.join(refdata_assembly_dir, 'biomarker','tsv', f"{db}.clinical.tsv.gz")
      if biomarker_vartype == 'CNA':
         logger.info(f"Loading CNA biomarker evidence from {db} ..")
      if biomarker_vartype == 'MUT':
         logger.info(f"Loading SNV/InDel biomarker evidence from {db} ..")
      if biomarker_vartype == 'FUSION':
         logger.info(f"Loading fusion biomarker evidence from {db} ..")
      

      biomarkers[db] = load_biomarkers(
         logger, variant_fname, clinical_fname, 
         biomarker_vartype = biomarker_vartype, 
         biomarker_variant_origin = biomarker_variant_origin)
      
      if biomarker_vartype in ['CNA','FUSION']:
         for key in biomarkers[db]['other_gene']:
            biomarker_data = biomarkers[db]['other_gene'][key]
            match_type = 'by_fusion'
            if biomarker_vartype == 'CNA':
               match_type = 'by_cna_segment'
            biomarker_item = (
               f"{db}|{biomarker_data[0]['variant_id']}|" 
               f"{biomarker_data[0]['clinical_evidence_items']}|{match_type}"
            )
            if key not in actionable_dict:               
               actionable_dict[key] = biomarker_item
            else:
               actionable_dict[key] = actionable_dict[key] + ',' + biomarker_item
            
         biomarkers['actionable_df'] = \
            pd.DataFrame(actionable_dict.items(), columns=['aberration_key', 'biomarker_match'])

   return biomarkers


def load_biomarkers(
   logger: logging.Logger, 
   biomarker_variant_fname: str, 
   biomarker_clinical_fname: str, 
   biomarker_vartype: str = 'MUT', 
   biomarker_variant_origin: str = 'Somatic'):
   """
   Loads biomarkers from the given files and returns a dictionary of variant biomarkers.

   Parameters:
   - logger: A logger object for logging messages.
   - biomarker_variant_fname: The file name of the biomarker variant data.
   - biomarker_clinical_fname: The file name of the biomarker clinical data.
   - biomarker_vartype: The type of biomarker variant (e.g., 'MUT','CNA', or 'FUSION')
   - biomarker_variant_origin: The origin of the biomarker variant (e.g., 'Somatic', 'Germline', or 'Both')

   Returns:
   - variant_biomarkers: A dictionary containing variant biomarkers. The keys are variant alias types 
     ('dbsnp', 'hgvsp', 'hgvsc', 'genomic', 'exon', 'other_gene', 'aa_region'), and the values are 
     dictionaries containing variant information.

   Note:
   - The function checks if the biomarker clinical file exists and logs an error message if it doesn't.
   - The function loads actionable cancer variants from the biomarker clinical file.
   - The function populates the 'variant_to_clinical_evidence' dictionary.
   - The function checks if the biomarker variant file exists and logs an error message if it doesn't.
   - The function reads the biomarker variant file and populates the 'variant_biomarkers' dictionary.
   """

   variant_biomarkers = {} ##dictionary to return
   for variant_alias_type in ['dbsnp','hgvsp','hgvsc','genomic','exon','other_gene','aa_region']:
      variant_biomarkers[variant_alias_type] = {}
   check_file_exists(biomarker_clinical_fname, logger)
   
   ## load actionable cancer variants from 'variant.tsv.gz'
   variant_to_clinical_evidence = {}
   with gzip.open(biomarker_clinical_fname, mode='rt') as f:
      reader = csv.DictReader(f, delimiter='\t')
      for row in reader:
         
         clinical_significance = re.sub(r' or ', '/', row['clinical_significance'])
         clinical_significance = re.sub(r'\s+', '_', clinical_significance)
         evidence_level = re.sub(r': .+', '', row['evidence_level'])
         evidence_type = re.sub(r'\s+', '_', row['evidence_type'])
         primary_site = re.sub(r'\s+', '_', row['primary_site'])
         cancer_type = row['cancer_type']
         if primary_site == ".":
            if cancer_type == "Cancer" or cancer_type == "Solid Tumor":
               primary_site = "Any"
            else:
               primary_site = "Undefined"
         variant_origin = row['variant_origin']
         molecular_profile_type = row['molecular_profile_type']
         if variant_origin == ".":
            variant_origin = "Unknown"
         
         if molecular_profile_type != 'Single':
            continue
         if biomarker_variant_origin == 'Somatic' and variant_origin != 'Somatic':
            continue
         if biomarker_variant_origin == 'Germline' and variant_origin != 'Germline':
            continue
         if biomarker_variant_origin == 'Both' and \
            variant_origin != 'Somatic' and \
               variant_origin != 'Germline':
            continue
         evidence_key = row['evidence_id'] + ':' + str(primary_site) + ':' + \
            str(clinical_significance) + ":" + str(evidence_level) + ":" + \
            str(evidence_type) + ":" + str(variant_origin)
         
         if row['variant_id'] not in variant_to_clinical_evidence:
            variant_to_clinical_evidence[row['variant_id']] = {}
         variant_to_clinical_evidence[row['variant_id']][evidence_key] = 1         

   f.close()

   for vid in variant_to_clinical_evidence.keys():
      variant_to_clinical_evidence[vid] = \
         '&'.join(sorted(variant_to_clinical_evidence[vid].keys()))
   check_file_exists(biomarker_variant_fname, logger)
   
   with gzip.open(biomarker_variant_fname, mode='rt') as f:
      reader = csv.DictReader(f, delimiter='\t')
      for row in reader:
         if 'alteration_type' in row.keys() and \
            'alias_type' in row.keys() and \
            'variant_exon' in row.keys() and \
            'entrezgene' in row.keys() and \
            'variant_id' in row.keys():
            if biomarker_vartype == 'MUT' and \
               (row['alteration_type'].startswith('MUT') or row['alteration_type'].startswith('CODON')):
               row['clinical_evidence_items'] = '.'
               if row['variant_id'] in variant_to_clinical_evidence.keys():
                  row['clinical_evidence_items'] = variant_to_clinical_evidence[row['variant_id']]    
               if row['clinical_evidence_items'] == '.':
                  continue           
               entry_alias_type = row['alias_type'].replace("_grch37", "").replace("_grch38", "")
              
               if entry_alias_type == "other_gene":
                  if bool(re.search(r'^((ACTIVATING )?Mutation|Loss|START Loss)$', row['variant_alias'])) is True:
                     varkey = str(row['entrezgene'])
                     if varkey not in variant_biomarkers['other_gene']:
                        variant_biomarkers['other_gene'][varkey] = []
                     variant_biomarkers['other_gene'][varkey].append(row)

               if entry_alias_type == 'exon':
                  exons = row['variant_exon']
                  if bool(re.search(r'&|-', exons)) is False:
                     for exon in str(exons).split('|'):
                        varkey = str(row['entrezgene']) + '_' + str(exon)
                        if varkey not in variant_biomarkers['exon']:
                           variant_biomarkers['exon'][varkey] = []
                        variant_biomarkers['exon'][varkey].append(row)

               for alias_type in ['dbsnp','hgvsp','hgvsc','genomic','aa_region']:                 
                  if entry_alias_type == alias_type:
                     if entry_alias_type != 'aa_region':
                        varkey = str(row['entrezgene']) + '_' + str(row['variant_alias'])
                        if varkey not in variant_biomarkers[alias_type]:
                           variant_biomarkers[alias_type][varkey] = []
                        variant_biomarkers[alias_type][varkey].append(row)
                     else:
                        if entry_alias_type == "aa_region":                           
                           aa_region = re.sub(r"aa(_|-)region:","",row["variant_alias"]).split("-")                           
                           aa_index = int(aa_region[0])
                           while aa_index <= int(aa_region[1]):
                              varkey = str(row['entrezgene']) + '_' + str(aa_index)
                              if varkey not in variant_biomarkers[entry_alias_type].keys():
                                 variant_biomarkers[entry_alias_type][varkey] = []
                              variant_biomarkers[entry_alias_type][varkey].append(row)
                              aa_index = aa_index + 1
            else:
               if biomarker_vartype == 'CNA' and \
                  ((row['alteration_type'].startswith('CNA')) or \
                     bool(re.search(r'^Loss$', row['variant_alias'])) is True):
                  row['clinical_evidence_items'] = '.'
                  if row['variant_id'] in variant_to_clinical_evidence.keys():
                     row['clinical_evidence_items'] = variant_to_clinical_evidence[row['variant_id']]                                 
                  
                  if row['alias_type'] == "other_gene":
                     if bool(re.search(r'^(Amplification|Deletion|Loss)$', row['variant_alias'])) is True:
                        varkey = str(row['entrezgene']) + "_" + \
                           re.sub(r"transcript_","",str(row['variant_consequence']))
                        if bool(re.search(r'^Loss$', row['variant_alias'])) is True: ## currently annotated with "loss_of_function_variant" as variant consequence
                           varkey = str(row['entrezgene']) + '_ablation'
                        if varkey not in variant_biomarkers['other_gene']:
                           variant_biomarkers['other_gene'][varkey] = []
                        del row['variant_exon']
                        del row['gene']
                        del row['alias_type']
                        variant_biomarkers['other_gene'][varkey].append(row)
               if biomarker_vartype == 'FUSION' and \
                  row['alteration_type'].startswith('FUSION') and \
                     bool(re.search(r'transcript_fusion', row['variant_consequence'])) is True:
                  row['clinical_evidence_items'] = '.'
                  if row['variant_id'] in variant_to_clinical_evidence.keys():
                     row['clinical_evidence_items'] = variant_to_clinical_evidence[row['variant_id']]
                  
                  if row['alias_type'] == "other_gene":
                     varkey = str(row['gene']) + "_transcript_fusion"
                     varkey2 = re.sub(r'((::v$)|(^v::))','::',str(row['gene'])) + "_transcript_fusion"
                     if varkey not in variant_biomarkers['other_gene']:
                        variant_biomarkers['other_gene'][varkey] = []
                     if varkey2 not in variant_biomarkers['other_gene']:
                        variant_biomarkers['other_gene'][varkey2] = []
                     del row['entrezgene']
                     del row['symbol']
                     del row['alias_type']
                     del row['variant_exon']
                     variant_biomarkers['other_gene'][varkey].append(row)
                     if varkey != varkey2:
                        variant_biomarkers['other_gene'][varkey2].append(row)

               

   f.close()

   return variant_biomarkers


def match_csq_biomarker(transcript_csq_elements, variant_biomarkers, rec, principal_csq_properties):

   ## Picked CSQ by VEP
   """
   Matches consequence entries from VEP against known variant biomarkers.

   This function identifies variant biomarkers by comparing transcript-specific
   consequence elements to entries in a biomarker database. It performs matching
   at various levels of resolution, including genomic coordinates, HGVSp, codon,
   amino-acid region, HGVSc, exon, and gene level.

   Parameters:
   - transcript_csq_elements: List of consequence entries from VEP for a given variant.
   - variant_biomarkers: Dictionary containing variant biomarker information.
   - rec: Record containing variant information (e.g., chromosomal coordinates).
   - principal_csq_properties: Dictionary of properties for the principal consequence
     selected by VEP, including 'hgvsp', 'hgvsc', 'entrezgene', 'codon', 'exon', and 'lof'.

   Returns:
   - Updates the 'rec.INFO' field with matched biomarker information, if applicable.
   """
   principal_hgvsp = principal_csq_properties['hgvsp']
   principal_hgvsc = principal_csq_properties['hgvsc']
   principal_entrezgene = principal_csq_properties['entrezgene']
   principal_codon = principal_csq_properties['codon'] 
   principal_exon = principal_csq_properties['exon']
   principal_lof = principal_csq_properties['lof'] 

   biomarker_hits_all = {}
   
   ## Match variant by genomic coordinate
   genomic_var_key = f"{principal_entrezgene}_{rec.CHROM}_{rec.POS}_{rec.REF}_{rec.ALT[0]}"

   if genomic_var_key in variant_biomarkers['genomic'].keys():
      hits_genomic = variant_biomarkers['genomic'][genomic_var_key]

      for ghit in hits_genomic:
         bkey_genomic = f"{ghit['biomarker_source']}|{ghit['variant_id']}|{ghit['clinical_evidence_items']}"
         if bkey_genomic not in biomarker_hits_all.keys():            
            biomarker_hits_all[bkey_genomic] = {}
         biomarker_hits_all[bkey_genomic]['by_genomic_coord'] = 1
         

   mut_lof_fs = False
   mut_protein = False
   principal_csq_hgvsp = False
   principal_csq_entrezgene = False
   principal_csq_hgvsc = False
   for csq_elem in transcript_csq_elements:
      (consequence, symbol, entrezgene, hgvsc, hgvsp, exon, feature_type, feature, biotype) = csq_elem.split(':')
      if bool(re.search(r'^(missense|stop|start|inframe|protein|frameshift)', consequence)) is True:
         mut_protein = True
      hgvsp_short = threeToOneAA(hgvsp)

      if hgvsp_short == principal_hgvsp:
         principal_csq_hgvsp = True
      if entrezgene == principal_entrezgene:
         principal_csq_entrezgene = True
      if hgvsc == principal_hgvsc:
         principal_csq_hgvsc = True
         
      if principal_csq_hgvsp is True and principal_csq_entrezgene is True:
         if bool(re.search(r'^(frameshift)', consequence)) is True and principal_lof is True:
            mut_lof_fs = True

      ## Match biomarkers by HGVSp or codon - "exact" or "codon" level resolution
      codon_match = []
      if entrezgene != "." and hgvsp != ".":
         biomarker_key_hgvsp = str(entrezgene) + '_' + str(hgvsp_short)
         codon_match = re.findall(r'p.[A-Z][0-9]{1,}',hgvsp_short)

         ## match biomarkers annotated at the amino acid level - with HGVSp coordinates
         if biomarker_key_hgvsp in variant_biomarkers['hgvsp']:
            hits = variant_biomarkers['hgvsp'][biomarker_key_hgvsp]
            for h in hits:
               bkey = f"{h['biomarker_source']}|{h['variant_id']}|{h['clinical_evidence_items']}"
               if bkey not in biomarker_hits_all.keys():            
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
                  if not chit['alteration_type'] == "CODON":
                     continue
                  bkey2 = f"{chit['biomarker_source']}|{chit['variant_id']}|{chit['clinical_evidence_items']}"
                  if bkey2 not in biomarker_hits_all.keys():            
                     biomarker_hits_all[bkey2] = {}
                  if principal_csq_entrezgene is True:
                     ## principal codon (VEP's picked csq)
                     if codon_match[0] == principal_codon:
                        biomarker_hits_all[bkey2]['by_codon_principal'] = 1
                     else:
                        ## nonprincipal codon
                        biomarker_hits_all[bkey2]['by_codon_nonprincipal'] = 1
      

      ## Match biomarkers by amino-acid (protein) region - "region level" resolution
      if entrezgene != "." and rec.INFO.get('AMINO_ACID_START') is not None and rec.INFO.get('AMINO_ACID_END') is not None:
         
         if principal_csq_entrezgene is True:
            aa_region_start_key = entrezgene + "_" + str(rec.INFO.get('AMINO_ACID_START'))
            aa_region_end_key = entrezgene + "_" + str(rec.INFO.get('AMINO_ACID_END'))
            if aa_region_start_key in variant_biomarkers['aa_region'].keys() and \
               aa_region_end_key in variant_biomarkers['aa_region'].keys():
               for rhit in variant_biomarkers['aa_region'][aa_region_start_key]:
                  
                  if re.search(str(consequence), rhit['variant_consequence']):
                     region_hit = f"{rhit['biomarker_source']}|{rhit['variant_id']}|{rhit['clinical_evidence_items']}"
                     if principal_csq_hgvsc is True:
                        if region_hit not in biomarker_hits_all.keys():
                           biomarker_hits_all[region_hit] = {}
                        biomarker_hits_all[region_hit]['by_aa_region_principal'] = 1
                     


      ## Match biomarkers by HGVSc identifier - "exact" resolution
      if entrezgene != "." and rec.INFO.get('HGVSc') is not None:
         hgvsc_elements = str(rec.INFO.get('HGVSc')).split(':')
         if len(hgvsc_elements) == 2:
            hgvsc_biomarker_key = str(entrezgene) + '_' + str(hgvsc_elements[1]) 
            if hgvsc_biomarker_key in variant_biomarkers['hgvsc'].keys():
               hits_hgvsc = variant_biomarkers['hgvsc'][hgvsc_biomarker_key]
               for hit_hgvsc in hits_hgvsc:
                  hgvsc_hit = f"{hit_hgvsc['biomarker_source']}|{hit_hgvsc['variant_id']}|{hit_hgvsc['clinical_evidence_items']}"
                  if hgvsc_hit not in biomarker_hits_all.keys():
                     biomarker_hits_all[hgvsc_hit] = {}
                  if principal_csq_hgvsc is True:
                     biomarker_hits_all[hgvsc_hit]['by_hgvsc_principal'] = 1
                  else:
                     biomarker_hits_all[hgvsc_hit]['by_hgvsc_nonprincipal'] = 1


      ## Match biomarkers indicated by exon number (and consequence) - "exon level" resolution
      if entrezgene != "." and principal_csq_entrezgene is True and exon != ".":
         exon_biomarker_key = str(entrezgene) + '_' + str(exon)
         if exon_biomarker_key in variant_biomarkers['exon'].keys():
            hits_exon = variant_biomarkers['exon'][exon_biomarker_key]

            for ehit in hits_exon:
               bkey4 = f"{ehit['biomarker_source']}|{ehit['variant_id']}|{ehit['clinical_evidence_items']}"
               exon_type_match = 'NA'
               if ehit['variant_consequence'] == "exon_variant" and \
                  re.search(r'Mutation', ehit['variant_alias']) and \
                     re.search(r'^(missense|coding|protein)', consequence):
                  exon_type_match = 'by_exon_mut'
               elif ehit['variant_consequence'] == "inframe_deletion" and re.search(r'^inframe_deletion', consequence):
                  exon_type_match = 'by_exon_deletion'
               else: 
                  if ehit['variant_consequence'] == "inframe_insertion" and re.search(r'^inframe_insertion', consequence):
                     exon_type_match = 'by_exon_insertion'
               if exon_type_match != 'NA':
                  if bkey4 not in biomarker_hits_all.keys():            
                     biomarker_hits_all[bkey4] = {}
                  if exon == principal_exon:
                     biomarker_hits_all[bkey4][exon_type_match + '_principal'] = 1
                  else:
                     biomarker_hits_all[bkey4][exon_type_match + '_nonprincipal'] = 1
               
      ## Match biomarkers indicated by gene only - "gene level" resolution
      if entrezgene != "." and principal_csq_entrezgene is True:
         if str(entrezgene) in variant_biomarkers['other_gene'].keys():
            hits_gene = variant_biomarkers['other_gene'][str(entrezgene)]
            for ghit in hits_gene:
               if ghit['variant_consequence'] == "wild_type":
                  continue
               bkey3 = f"{ghit['biomarker_source']}|{ghit['variant_id']}|{ghit['clinical_evidence_items']}"
               gene_type_match = 'NA'
               ## match biomarkers annotated as "Mutation" only for a given gene - 
               ## consider only the principal consequence (VEP's picked)
               if (ghit['alteration_type'] == 'MUT' or ghit['alteration_type'] == 'MUT_ONC') and \
                  mut_protein is True and (principal_csq_hgvsc is True or principal_csq_hgvsp is True):
                     
                     ## ensure matching of inframe insertion and deletion biomarkers with VEP's inframe insertion 
                     ## and deletion consequences, respectively;
                  if ghit['variant_consequence'] == "inframe_deletion" and re.search(r'^inframe_deletion', consequence) or \
                     ghit['variant_consequence'] == "frameshift_variant" and re.search(r'^frameshift', consequence) or \
                     ghit['variant_consequence'] == "inframe_insertion" and re.search(r'^inframe_insertion', consequence) or \
                     ghit['variant_consequence'] == "start_lost" and re.search(r'^start_lost', consequence) or \
                     (re.search(r'(missense_|gain_of_function_|protein_altering|gene_variant|transcript_variant)',ghit['variant_consequence']) and \
                        re.search(r'^(missense|coding|protein)', consequence)):
                        #print(f"Matching gene-level biomarker {ghit['variant_id']} {ghit['clinical_evidence_items']} {ghit['variant_consequence']} with VEP consequence {consequence} ..")
                        gene_type_match = 'by_gene_mut'
                      
               ## match biomarkers annotated as "Mutation - loss of function" only for a given gene - 
               ## consider only the principal consequence (VEP's picked)
               elif ghit['alteration_type'] == 'MUT_LOF' and principal_lof is True and \
                  (principal_csq_hgvsc is True or principal_csq_hgvsp is True):
                     gene_type_match = 'by_gene_mut_lof'
               else:
                  ## match biomarkers annotated as "Mutation - loss of function - frameshift" only for a given gene - 
                  ## consider only the principal consequence (VEP's picked)
                  if ghit['alteration_type'] == 'MUT_LOF_FS' and principal_lof is True and mut_lof_fs is True and \
                     (principal_csq_hgvsc is True or principal_csq_hgvsp is True):
                     gene_type_match = 'by_gene_mut_lof_fs'

               if gene_type_match != 'NA':
                  if bkey3 not in biomarker_hits_all.keys():            
                     biomarker_hits_all[bkey3] = {}
                  biomarker_hits_all[bkey3][gene_type_match] = 1

      mut_lof_fs = False
      mut_protein = False
      principal_csq_hgvsp = False
      principal_csq_entrezgene = False
      principal_csq_hgvsc = False
   
   if len(biomarker_hits_all.keys()) > 0:
      biomarker_var_matches = {}
      for bm in biomarker_hits_all.keys():
         match_types = ':'.join(sorted(biomarker_hits_all[bm].keys()))

         ## Aggregate all types of matching performed ('by_exon, by_hgvsp, by_hgvsc' etc.) for a given variant
         vkey = f"{bm}|{match_types}"
         biomarker_var_matches[vkey] = 1

      if rec.INFO.get('BIOMARKER_MATCH') is not None:
         rec.INFO['BIOMARKER_MATCH'] = rec.INFO['BIOMARKER_MATCH'] + ',' + str(','.join(sorted(biomarker_var_matches.keys())))
      else:
         rec.INFO['BIOMARKER_MATCH'] = str(','.join(sorted(biomarker_var_matches.keys())))

   return



def validate_oncokb_input_file(path, required_columns):
   """
   Check whether variant data file exists and contains all required columns.
   """
   if path is None or not os.path.isfile(path):
      return False, f"File not found: {path}"

   with open(path) as f:
      header = f.readline().strip().split("\t")

   missing = required_columns - set(header)
   if missing:
      return False, f"Missing columns in {path}: {missing}"

   return True, None


def filter_oncokb_output(output_path, logger=None):
   """
   Filter OncoKB annotator output to retain only rows where
   GENE_IN_ONCOKB is True. Overwrites the file in place.
   """
   logger.info(f"Filtering OncoKB output: {os.path.basename(output_path)}")
   if not os.path.isfile(output_path):
      return
   df = pd.read_csv(output_path, sep="\t", low_memory=False)
   cols = df.columns.tolist()
   if "GENE_IN_ONCOKB" in cols:
      df = df[
         (df["GENE_IN_ONCOKB"].astype(str).str.strip().str.lower() == "true")
      ]
      df.to_csv(output_path, sep="\t", index=False)


def fetch_oncokb_cancer_genes(oncokb_token: str, logger=None) -> set:
   """Fetch the OncoKB cancer gene list and return a set of Hugo symbols.
   Returns an empty set on any failure so the caller can proceed without filtering."""
   url = "https://www.oncokb.org/api/v1/utils/cancerGeneList"
   req = urllib.request.Request(
      url,
      headers={
         "accept": "application/json",
         "Authorization": f"Bearer {oncokb_token}"
      }
   )
   try:
      with urllib.request.urlopen(req, timeout=30) as response:
         genes = json.loads(response.read().decode())
         symbols = {g["hugoSymbol"] for g in genes if "hugoSymbol" in g}
         if logger:
            logger.info(f"OncoKB cancer gene list: fetched {len(symbols)} genes")
         return symbols
   except Exception as e:
      if logger:
         logger.warning(f"Could not fetch OncoKB cancer gene list: {e} - skipping CNA pre-filter")
      return set()


def run_oncokb_annotator(
   sample_name: str,
   oncokb_token: str,
   oncotree_code: str = None,
   build: str = "GRCh38",
   maf_query_file: str = None,
   fusion_query_file: str = None,
   cna_query_file: str = None,
   output_dir: str = ".",
   logger: logging.Logger = None):

    
   """
   Wrapper for OncoKB variant annotation through the API.

   Workflow
    - At least one input data file must be provided.
    - Validate headers of each input data file.
    - Produce a clinical file that is necessary for each input script (MafAnnotator.py, FusionAnnotator.py, CNAAnnotator.py).
    - Run annotators only for valid inputs.
   """

   if logger is None:
      logger = logging.getLogger("pcgr-oncokb-annotation")
   # Ensure output directory exists
   if not os.path.isdir(output_dir):
      err_msg = f"Output directory {output_dir} does not exist."
      error_message(err_msg, logger)
      
   # 1) Validate that at least one variant data file is present
   input_files = [maf_query_file, fusion_query_file, cna_query_file]
   if not any(f and os.path.isfile(f) for f in input_files):
      err_msg = "At least one valid input file (MAF, fusion, CNA) must be provided for OncoKB annotation."
      error_message(err_msg, logger)
  
   # 2) Validate contents of each input file
   validations = [
      ("MAF", maf_query_file, pcgr_vars.ONCOKB_MAF_REQUIRED_COLS),
      ("FUSION", fusion_query_file, pcgr_vars.ONCOKB_FUSION_REQUIRED_COLS),
      ("CNA", cna_query_file, pcgr_vars.ONCOKB_CNA_REQUIRED_COLS),
   ]

   valid_inputs = {}

   for label, path, required_cols in validations:
      if path:  # if file provided provided
         ok, error = validate_oncokb_input_file(path, required_cols)
         if not ok:
            err_msg = f"{label} file failed validation: {error}"
            error_message(err_msg, logger)
         else:
            valid_inputs[label] = path

   if not valid_inputs:
      err_msg = "No valid input files passed header validation."
      error_message(err_msg, logger)

   # 3) Pre-filter CNA input to OncoKB cancer genes
   if "CNA" in valid_inputs:
      cancer_genes = fetch_oncokb_cancer_genes(oncokb_token, logger)
      if cancer_genes:
         cna_df = pd.read_csv(valid_inputs["CNA"], sep="\t")
         n_before = len(cna_df)
         cna_df = cna_df[cna_df["Hugo_Symbol"].isin(cancer_genes)]
         n_after = len(cna_df)
         logger.info(f"CNA pre-filter: {n_before} → {n_after} events after restricting to OncoKB cancer genes")
         if n_after == 0:
            logger.warning("No CNA events remain after OncoKB cancer gene filter - skipping CNA annotation")
            del valid_inputs["CNA"]
         else:
            cna_df.to_csv(valid_inputs["CNA"], sep="\t", index=False)

   # 4) Create clinical file (only when oncotree_code is provided)
   clinical_file = None
   if oncotree_code is not None:
      clinical_content = (
         "SAMPLE_ID\tONCOTREE_CODE\n"
         f"{sample_name}\t{oncotree_code}\n")
      clinical_file = os.path.join(output_dir, f"{sample_name}_pcgr_oncokb_clinical.txt")
      with open(clinical_file, "w") as f:
         f.write(clinical_content)
   else:
      logger.info("No OncoTree code specified - running OncoKB annotation without tumor type")

   # 5) Annotator script paths
   scripts = {
      "MAF": "MafAnnotator.py",
      "FUSION": "FusionAnnotator.py",
      "CNA": "CnaAnnotator.py"
   }

   # 6) Run annotators for valid inputs
   outputs = {}

   for label, input_path in valid_inputs.items():
      script = scripts[label]

      # For MAF files, run twice with different query types (HGVSp and HGVSg)
      if label == "MAF":
         query_types = ["HGVSp", "HGVSg"]
         for query_type in query_types:
            output_path = os.path.join(
               output_dir,
               f"{sample_name}.pcgr.{str(build).lower()}.oncokb_output.{label.lower()}.{query_type.lower()}.txt"
            )

            # Create log file for this annotator run
            log_file = os.path.join(
               output_dir,
               f"{sample_name}.pcgr.{str(build).lower()}.oncokb.{label.lower()}.{query_type.lower()}.log"
            )

            cmd = [
               script,
               "-i", input_path,
               "-o", output_path,
               "-b", oncokb_token,
               "-q", query_type,
               "-r", build,
               "-d"
            ]
            if clinical_file is not None:
               cmd.extend(["-c", clinical_file])

            logger.info(f"Running OncoKB-{label} annotator with query type {query_type}...")
            logger.info(f"OncoKB {label} {query_type} log: {os.path.basename(log_file)}")

            # Redirect stdout and stderr to log file
            with open(log_file, 'w') as log_fh:
               subprocess.run(cmd, check=True, stdout=log_fh, stderr=subprocess.STDOUT)

            logger.info(f"OncoKB-{label} annotation complete (query type: {query_type}).")
            filter_oncokb_output(output_path, logger)

            # Store both outputs for MAF
            if label not in outputs:
               outputs[label] = {}
            outputs[label][query_type] = output_path
      else:
         # For FUSION and CNA, run once
         output_path = os.path.join(
            output_dir,
            f"{sample_name}.pcgr.{str(build).lower()}.oncokb_output.{label.lower()}.txt"
         )

         # Create log file for this annotator run
         log_file = os.path.join(
            output_dir,
            f"{sample_name}.pcgr.{str(build).lower()}.oncokb.{label.lower()}.log"
         )

         cmd = [
            script,
            "-i", input_path,
            "-o", output_path,
            "-b", oncokb_token,
            "-d"
         ]
         if clinical_file is not None:
            cmd.extend(["-c", clinical_file])

         if label == "CNA":
            cmd.extend(["-f", "individual"])

         logger.info(f"Running OncoKB-{label} annotator...")
         logger.info(f"OncoKB {label} log: {os.path.basename(log_file)}")

         # Redirect stdout and stderr to log file
         with open(log_file, 'w') as log_fh:
            subprocess.run(cmd, check=True, stdout=log_fh, stderr=subprocess.STDOUT)

         logger.info(f"OncoKB-{label} annotation complete.")
         filter_oncokb_output(output_path, logger)

         outputs[label] = output_path

   # 7) Return summary of outputs
   return {
      "clinical_file": clinical_file,
      "outputs": outputs
   }