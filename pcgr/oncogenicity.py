#!/usr/bin/env python

import os
import re
import csv
import logging

from logging import Logger
import gzip
import pandas as pd

from pcgr.annoutils import threeToOneAA
from pcgr import pcgr_vars
from pcgr.pcgr_vars import (
    ONCOGENICITY_THRESHOLDS,
    ONCOGENICITY_CODE_ORDER,
    OKB_ONCOGENICITY_CRITERIA,
)


def assign_oncogenicity_evidence(rec = None, oncogenicity_criteria = None, tumortype = "Any"):

   ### Benign oncogenic effects of somatic variants
   
   #  1) "ONCG_SBVS1"
   ## Very high MAF: > 0.05 in gnomAD - any five general continental populations
   ## AFR/AMR/EAS/NFE/SAS

   # 2) "ONCG_SBS1"
   ## High MAF: > 0.01 in gnomAD - any five general continental populations
   ## AFR/AMR/EAS/NFE/SAS

   # 3) "ONCG_SBP1"
   ## Multiple lines of computational evidence support a benign
   ## effect on the gene or gene product
   ## (conservation, evolutionary, splicing impact, etc. - from dbNSFP

   # 4) ONCG_SBP2"
   ## Synonymous (silent) variant for which splicing prediction
   ## algorithms predict no effect on the splice consensus sequence
   ## nor the creation of a new splice site and the nucleotide is
   ## not highly conserved

   # 5) "ONCG_SBS2"
   ## Well established in invitro/in vivo functional studies show
   ## no oncogenic effects

   ### Oncogenic effects of somatic variants
   
   # 6) "ONCG_OVS1"
   ## Null variant - predicted as LoF
   ## - Nonsense, frameshift, canonical splice sites, initiation codon,
   ##   single-exon/multi-exon deletion
   ## - Tumor suppressor gene

   # 7) "ONCG_OS1"
   ## Same amino acid change as previously established oncogenic variant
   ## (ClinVar - oncogenicity classifications)

   # 8) "ONCG_OS2"
   ## Well established in invitro/in vivo functional studies show
   ## oncogenic effect of the variant

   # 9) "ONCG_OS3"
   ## Located in a mutation hotspot
   ## - >= 50 samples with a somatic variant at the AA position (cancerhotspots.org)
   ## - Same amino acid change in >= 10 samples (cancerhotspots.org)

   # 10) "ONCG_OM1"
   ## Located in a critical and well-established part of a functional domain
   ## (active site of an enzyme)
   ## NOTE: Here we assume that OM1 is applicable if there is overlap with 
   ## oncogenic/predictive evidence (CIViC/CGI)

   # 11) "ONCG_OM2"
   ## Protein length changes as a result of in-frame deletions/insertions in a
   ## known oncogene/tumor suppressor genes or stop-loss variants in a
   ## tumor suppressor gene
   # - Not applicable if OVS1 is applicable

   # 12) "ONCG_OM3"
   ## Missense variant at an amino acid residue where a different missense
   ## variant determined to be oncogenic (using this standard) has been
   ## documented. Amino acid difference from reference amino acid (Grantham distance) 
   ## should be greater or at least approximately the same as for missense 
   ## change determined to be oncogenic.

   # 13) "ONCG_OM4"
   ## Located in a mutation hotspot
   ## - < 50 samples with a somatic variant at the AA position (cancerhotspots.org)
   ## - Same amino acid change in >= 10 samples (cancerhotspots.org)
   ## - Not applicable if OM1 or OM3 is applicable

   # 14) "ONCG_OP1"
   ## All used lines of computational support an oncogenic effect
   ## of a variant (conservation, evolutionary, splicing effect)

   # 15) "ONCG_OP2"
   ## Somatic variant in a gene in a malignancy with a single genetic
   ## etiology. Example:retinoblastoma is caused by bi-allelic
   ## RB1 inactivation.

   # 16) "ONCG_OP3"
   ## Located in a mutation hotspot
   ## - Same amino acid change in < 10 samples (cancerhotspots.org)
   ## - Not applicable if OM1 or OM3 is applicable

   # 17) "ONCG_OP4"
   ## Absent from controls (gnomAD)
   ## - Extremely low MAF


   required_oncogenicity_vars = [
      "Consequence",
      "MUTATION_HOTSPOT",
      "MUTATION_HOTSPOT_CANCERTYPE",
      "KNOWN_ONCOGENIC",
      "KNOWN_ONCOGENIC_SITE",
      "SYMBOL",
      "GRANTHAM_DISTANCE",
      "BIOMARKER_MATCH",
      "ONCOGENE",
      "ONCOGENE_EVIDENCE",
      "EFFECT_PREDICTIONS",
      "TSG",
      "TSG_EVIDENCE",
      "LOSS_OF_FUNCTION",
      "INTRON_POSITION",
      "EXON_POSITION",
      "CODING_STATUS",
      "SPLICE_DONOR_RELEVANT",
      "gnomADe_EAS_AF",
      "gnomADe_NFE_AF",
      "gnomADe_AFR_AF",
      "gnomADe_AMR_AF",
      "gnomADe_SAS_AF",
      "DBNSFP_SIFT",
      "DBNSFP_PROVEAN",
      "DBNSFP_META_RNN",
      "DBNSFP_MUTATIONTASTER",
      "DBNSFP_DEOGEN2",
      "DBNSFP_PRIMATEAI",
      "DBNSFP_MUTATIONASSESSOR",
      "DBNSFP_FATHMM_XF",
      "DBNSFP_M_CAP",
      "DBNSFP_POLYPHEN2_HVAR",
      "DBNSFP_ALPHA_MISSENSE",
      "DBNSFP_PHACTBOOST",
      "DBNSFP_CLINPRED",
      "DBNSFP_LIST_S2",
      "DBNSFP_BAYESDEL_ADDAF",
      "DBNSFP_SPLICE_SITE_ADA",
      "DBNSFP_SPLICE_SITE_RF"]

   variant_data = {}
   vcf_record_info = getattr(rec, "INFO", None)
   for col in required_oncogenicity_vars:
      if vcf_record_info is None or vcf_record_info.get(col) is None:        
         if col == "TSG" or col == "ONCOGENE" or \
            col == "SPLICE_DONOR_RELEVANT" or col == "LOSS_OF_FUNCTION":
            variant_data[col] = False
         else:
            variant_data[col] = None
      else:
         if vcf_record_info.get(col) == '':
            variant_data[col] = True
         else: 
            variant_data[col] = vcf_record_info.get(col)
   
   oncogenicity_criteria = oncogenicity_criteria or {}
   for code in oncogenicity_criteria.keys():
      variant_data[code] = False
   
   if "KNOWN_ONCOGENIC" in variant_data.keys():
      if variant_data['KNOWN_ONCOGENIC'] is not None:         
         variant_data['ONCG_OS1'] = True
   
   for elem in ['N_INSILICO_CALLED',
                'N_INSILICO_DAMAGING',
                'N_INSILICO_TOLERATED',
                'N_INSILICO_SPLICING_NEUTRAL',
                'N_INSILICO_SPLICING_AFFECTED']:
      variant_data[elem] = 0

   for algo in ['SIFT',
                'PROVEAN',
                'META_RNN',
                'ALPHA_MISSENSE',
                'POLYPHEN2_HVAR',
                'PHACTBOOST',
                'CLINPRED',
                'MUTATIONTASTER',
                'DEOGEN2',
                'PRIMATEAI',
                'MUTATIONASSESSOR',
                'FATHMM_XF',
                'M_CAP',
                'LIST_S2',
                'BAYESDEL_ADDAF',
                'SPLICE_SITE_RF',
                'SPLICE_SITE_ADA']:
      col = 'DBNSFP_' + str(algo)
      if col in variant_data.keys():
         if variant_data[col] != '.':
            variant_data['N_INSILICO_CALLED'] += 1
         if variant_data[col] == 'D' or variant_data[col] == 'PD':
            variant_data['N_INSILICO_DAMAGING'] += 1
         if variant_data[col] == 'T':
            variant_data['N_INSILICO_TOLERATED'] += 1
         if (col == 'DBNSFP_SPLICE_SITE_RF' or \
             col == 'DBNSFP_SPLICE_SITE_ADA') and \
            variant_data[col] == 'SN':
            variant_data['N_INSILICO_SPLICING_NEUTRAL'] += 1
         if (col == 'DBNSFP_SPLICE_SITE_RF' or \
             col == 'DBNSFP_SPLICE_SITE_ADA') and \
            variant_data[col] == 'AS':
            variant_data['N_INSILICO_SPLICING_AFFECTED'] += 1

   ## majority of insilico predictors predict damaging effect
   if (variant_data['N_INSILICO_CALLED'] > pcgr_vars.ONCOGENICITY['insilico_pred_min_majority'] and \
      variant_data['N_INSILICO_DAMAGING'] >= pcgr_vars.ONCOGENICITY['insilico_pred_min_majority'] and \
      variant_data['N_INSILICO_TOLERATED'] <= pcgr_vars.ONCOGENICITY['insilico_pred_max_minority'] and \
      variant_data['N_INSILICO_SPLICING_NEUTRAL'] <= 1) or \
         (variant_data['N_INSILICO_CALLED'] > pcgr_vars.ONCOGENICITY['insilico_pred_min_majority'] and \
      variant_data['N_INSILICO_DAMAGING'] >= pcgr_vars.ONCOGENICITY['insilico_pred_min_majority'] and \
      variant_data['N_INSILICO_TOLERATED'] <= pcgr_vars.ONCOGENICITY['insilico_pred_max_minority'] + 1 and \
      "alphamissense:D" in variant_data['EFFECT_PREDICTIONS'] and \
      variant_data['N_INSILICO_SPLICING_NEUTRAL'] <= 1) or \
         variant_data['N_INSILICO_SPLICING_AFFECTED'] == 2:
         variant_data['ONCG_OP1'] = True
   
   ## majority of insilico predictors predict tolerated effect
   if (variant_data['N_INSILICO_CALLED'] > pcgr_vars.ONCOGENICITY['insilico_pred_min_majority'] and \
      variant_data['N_INSILICO_TOLERATED'] >= pcgr_vars.ONCOGENICITY['insilico_pred_min_majority'] and \
      variant_data['N_INSILICO_DAMAGING'] <= pcgr_vars.ONCOGENICITY['insilico_pred_max_minority'] and \
      variant_data['N_INSILICO_SPLICING_AFFECTED'] == 0) or \
         (variant_data['N_INSILICO_CALLED'] > pcgr_vars.ONCOGENICITY['insilico_pred_min_majority'] and \
      variant_data['N_INSILICO_TOLERATED'] >= pcgr_vars.ONCOGENICITY['insilico_pred_min_majority'] and \
      variant_data['N_INSILICO_DAMAGING'] <= pcgr_vars.ONCOGENICITY['insilico_pred_max_minority'] + 1 and \
      "alphamissense:T" in variant_data['EFFECT_PREDICTIONS'] and \
      variant_data['N_INSILICO_SPLICING_AFFECTED'] == 0):
         variant_data['ONCG_SBP1'] = True

   hotspot_mutated_samples = {}
   hotspot_mutated_samples['aa_variant'] = {}
   hotspot_mutated_samples['aa_site'] = {}
   hotspot_mutated_samples['aa_variant']['Any'] = 0
   hotspot_mutated_samples['aa_site']['Any'] = 0

   if variant_data['MUTATION_HOTSPOT_CANCERTYPE'] is not None:
      ttype_samples_in_hotspots = variant_data['MUTATION_HOTSPOT_CANCERTYPE'].split(',')
      for ttype in ttype_samples_in_hotspots:
         ttype_stats = ttype.split('|')
         if len(ttype_stats) == 3:
            ttype_hotspot = ttype_stats[0]
            ttype_hotspot_samples_site = int(ttype_stats[1])
            ttype_hotspot_samples_aa = ttype_hotspot_samples_site
            if ttype_stats[2] != '':
               ttype_hotspot_samples_aa = int(ttype_stats[2])
            
            if ttype_hotspot != "Unknown":
               hotspot_mutated_samples['aa_variant'][ttype_hotspot] = ttype_hotspot_samples_aa
               hotspot_mutated_samples['aa_site'][ttype_hotspot] = ttype_hotspot_samples_site
            
            hotspot_mutated_samples['aa_variant']['Any'] = hotspot_mutated_samples['aa_variant']['Any'] + ttype_hotspot_samples_aa
            hotspot_mutated_samples['aa_site']['Any'] = hotspot_mutated_samples['aa_site']['Any'] + ttype_hotspot_samples_site
   
   
   ## The SOP's for assessment of overlap with mutational hotspots do not specify that occurrences should
   ## match the tumor type (example in paper ignores tumor-type matching when looking for number of samples overlapping with hotspot)
   ttype = "Any"
   #if tumortype != "Any":
   #   ttype = re.sub(r'/', '@', re.sub(r' ','_', tumortype))

   if ttype in hotspot_mutated_samples['aa_variant'].keys() and \
      ttype in hotspot_mutated_samples['aa_site'].keys():
      if hotspot_mutated_samples['aa_variant'][ttype] >= 10 and \
         hotspot_mutated_samples['aa_site'][ttype] >= 50:
            variant_data['ONCG_OS3'] = True
      if hotspot_mutated_samples['aa_variant'][ttype] >= 10 and \
         hotspot_mutated_samples['aa_site'][ttype] < 50:
            variant_data['ONCG_OM4'] = True
      if hotspot_mutated_samples['aa_variant'][ttype] < 10 and \
         hotspot_mutated_samples['aa_variant'][ttype] > 0:
            variant_data['ONCG_OP3'] = True
   

   ## For the OM1 criterion ("Located in a critical and well-established part of a functional domain"), we presently
   ## lack the access to a resource with such information. As a simplified means to gather some evidence in this regard,
   ## we rather base our criteria matching based on existing oncogenicity/prognostic/diagnostic evidence, for which 
   ## variants are presumably located in critical sites of functional domains (in that sense, indirect evidence for OM1)
   if variant_data['BIOMARKER_MATCH'] is not None:
      
      ## Split all biomarker evidence into a list
      biomarker_evidence = variant_data['BIOMARKER_MATCH'].split(',')
     
      for eitem in biomarker_evidence:
                  
         ## Example 'eitem' element:
         ## cgi|659|CGI1077:Pancreas:Sensitivity/Response:C:Predictive:Somatic|by_hgvsp_principal
         ## EID510:Myeloid:Poor_Outcome:B:Prognostic:Somatic|by_gene_mut
         if ('Oncogenic' in eitem or \
            'Poor_Outcome:B:Prognostic' in eitem or \
            'Positive:B:Diagnostic' in eitem or \
            'Positive:A:Diagnostic' in eitem or \
            'Poor_Outcome:A:Prognostic' in eitem) and \
            'Somatic' in eitem and \
            ('by_genomic_coord' in eitem or \
             'by_hgvsc_principal' in eitem or \
             'by_hgvsp_principal' in eitem or \
             'by_codon_principal' in eitem or \
             'by_aa_region_principal' in eitem):               
               ## only applicable if OS3 and OS1 is not set
               if variant_data['ONCG_OS3'] is False and variant_data['ONCG_OS1'] is False:
                  variant_data['ONCG_OM1'] = True
         
         ## Catch prognostic/diagnostic/oncogenic non-coding variants (e.g. TERT) - these will rank at the top
         ## with respect to oncogenicity (altough not classified as oncogenic in the strict sense)
         if ('Oncogenic' in eitem or \
            'Positive:B:Diagnostic' in eitem or \
            'Positive:A:Diagnostic' in eitem or \
            'Poor_Outcome:B:Prognostic' in eitem or \
            'Poor_Outcome:A:Prognostic' in eitem) and \
            variant_data['CODING_STATUS'] == 'noncoding' and \
            'Somatic' in eitem and \
            ('by_genomic_coord' in eitem or \
             'by_hgvsc_principal' in eitem):
               ## only applicable if OS3 and OS1 is not set
               if variant_data['ONCG_OS3'] is False and variant_data['ONCG_OS1'] is False:
                  variant_data['ONCG_OM1'] = True

   all_gnomad_tags = pcgr_vars.GNOMAD_MAIN_EXOME_AF_TAGS + pcgr_vars.GNOMAD_MAIN_GENOME_AF_TAGS  

   ## check if variant has MAF > 0.01 (SBVS1) or > 0.05 in any of five major gnomAD subpopulations (exome or genome set)
   for pop in all_gnomad_tags:
      if pop not in variant_data.keys():
         continue
      if variant_data[pop] is not None:
         if float(variant_data[pop]) >= float(pcgr_vars.ONCOGENICITY['gnomAD_very_common_AF']):
            variant_data["ONCG_SBVS1"] = True
         if float(variant_data[pop]) >= float(pcgr_vars.ONCOGENICITY['gnomAD_common_AF']) and variant_data["ONCG_SBVS1"] is False:
            variant_data["ONCG_SBS1"] = True
   
   gnomad_tags = {}
   approx_zero_pop_freq = {}
      
   for assay in ['exome', 'genome']:
      gnomad_tags[assay] = pcgr_vars.GNOMAD_MAIN_EXOME_AF_TAGS if assay == 'exome' else pcgr_vars.GNOMAD_MAIN_GENOME_AF_TAGS
      approx_zero_pop_freq[assay] = 0
         
      for pop in gnomad_tags[assay]:
         if pop not in variant_data.keys():
            continue
         ## no AF recorded in gnomAD for this population
         if variant_data[pop] is None:
            approx_zero_pop_freq[assay] = approx_zero_pop_freq[assay] + 1
         else:
            ## Extremely low AF for this population
            if float(variant_data[pop]) < float(pcgr_vars.ONCOGENICITY['gnomAD_extremely_rare_AF']):
               approx_zero_pop_freq[assay] = approx_zero_pop_freq[assay] + 1
      ## check if variant is missing or with AF approximately zero in all five major gnomAD subpopulations (exome or genome set)
      if approx_zero_pop_freq[assay] == 5:
         variant_data["ONCG_OP4"] = True
   
   ## check if variant is a loss-of-function variant in a tumor suppressor gene (Cancer Gene Census/CancerMine)
   if "TSG" in variant_data.keys() and \
      "ONCOGENE" in variant_data.keys() and \
      "LOSS_OF_FUNCTION" in variant_data.keys() and \
      "Consequence" in variant_data.keys():

      if variant_data['LOSS_OF_FUNCTION'] is True and variant_data['TSG'] is True:
         variant_data['ONCG_OVS1'] = True
  
      ## check if variant is creating a stop-lost or protein-length change in oncogene/tumor suppressor genes
      if variant_data['ONCG_OVS1'] is False and \
         ((re.match(r'^(inframe_deletion|inframe_insertion)', variant_data['Consequence']) and \
            (variant_data['TSG'] is True or variant_data['ONCOGENE'] is True)) or \
         (re.match(r'^(stop_lost)', variant_data['Consequence']) and \
            variant_data['TSG'] is True)):
            variant_data['ONCG_OM2'] = True
   
   ## check if variant is silent (synonymous|splice) and outside critical splice region
   if "INTRON_POSITION" in variant_data.keys() and \
      "EXON_POSITION" in variant_data.keys() and \
      "DBNSFP_SPLICE_SITE_RF" in variant_data.keys() and \
      "Consequence" in variant_data.keys():
      
      if (int(variant_data['INTRON_POSITION']) < 0 and int(variant_data['INTRON_POSITION']) < -3 or \
         int(variant_data['INTRON_POSITION']) > 0 and int(variant_data['INTRON_POSITION']) > 6 or \
         int(variant_data['EXON_POSITION']) < 0 and int(variant_data['EXON_POSITION']) < -2 or \
         int(variant_data['EXON_POSITION']) > 0 and int(variant_data['EXON_POSITION']) > 1) and \
         variant_data['DBNSFP_SPLICE_SITE_RF'] != "AS" and \
         re.match(r'^(synonymous_variant|splice_region_variant)', variant_data['Consequence']):
            variant_data['ONCG_SBP2'] = True
   
   if "KNOWN_ONCOGENIC_SITE" in variant_data.keys():
      if variant_data['KNOWN_ONCOGENIC_SITE'] is not None:         
         if variant_data['ONCG_OS1'] is False and \
            variant_data['ONCG_OS3'] is False and \
               variant_data['ONCG_OM1'] is False:
            variant_data['ONCG_OM3'] = True

   variant_data['ONCOGENICITY'] = "VUS"
   variant_data["ONCOGENICITY_DOC"] = "."
   variant_data["ONCOGENICITY_CODE"] = "."
   variant_data["ONCOGENICITY_SCORE"] = 0
   onc_score_pathogenic = 0
   onc_score_benign = 0

   #for code in clingen_vicc_ev_codes:
   for code in oncogenicity_criteria.keys():
      if code in variant_data.keys():
         if variant_data[code] is True:            
            score = float(oncogenicity_criteria[code]['score'])
            pole = oncogenicity_criteria[code]['pole']
            if pole == "P":
               onc_score_pathogenic = onc_score_pathogenic + score
            else:
               onc_score_benign = onc_score_benign + score

            if variant_data['ONCOGENICITY_CODE'] == ".":
               variant_data['ONCOGENICITY_CODE'] = code
            else:
               variant_data['ONCOGENICITY_CODE'] = f'{variant_data["ONCOGENICITY_CODE"]}|{code}'

   variant_data['ONCOGENICITY_SCORE'] = onc_score_benign + onc_score_pathogenic
   variant_data['ONCOGENICITY_CODE'] = _sort_oncogenicity_codes(
      str(variant_data['ONCOGENICITY_CODE']))
   variant_data['ONCOGENICITY'] = _classify_oncogenicity_score(
      variant_data['ONCOGENICITY_SCORE'],
      str(variant_data['ONCOGENICITY_CODE'])
   )
   
   for e in ['ONCOGENICITY_SCORE',
             'ONCOGENICITY',
             'ONCOGENICITY_CODE']:
      if vcf_record_info is not None:
         vcf_record_info[e] = variant_data[e]

   return(rec)

def load_oncogenicity_criteria(oncogenic_criteria_fname: str, logger: Logger):
   """
   Load oncogenic criteria from a file and create a dictionary of criteria.
   """
   
   oncogenic_criteria = {}
   if not os.path.exists(oncogenic_criteria_fname):
      logger.info(f"ERROR: File '{oncogenic_criteria_fname}' does not exist - exiting")
      exit(1)

   with open(oncogenic_criteria_fname, mode='rt') as f:
      reader = csv.DictReader(f, delimiter='\t')
      for row in reader:            
         oncogenic_criteria[row['code']] = row
         
   return oncogenic_criteria

def load_oncogenic_variants(oncogenic_variants_fname: str, logger: Logger):
   """
   Load oncogenic variants from a file and create a dictionary of variants.
   """
   
   oncogenic_variants = {}
   if not os.path.exists(oncogenic_variants_fname):
      logger.info(f"ERROR: File '{oncogenic_variants_fname}' does not exist - exiting")
      exit(1)

   with gzip.open(oncogenic_variants_fname, mode='rt') as f:
      reader = csv.DictReader(f, delimiter='\t')
      for row in reader:            
         gene = str(row['entrezgene'])
         if 'oncogenic' not in str(row['oncogenicity']).lower():
            continue         
         oncogenic_variants[str(gene) + '-' + str(row['var_id'])] = row         
         if 'grantham_distance' in row.keys():
            if row['grantham_distance'] in ('NA', '', None):
            #if row['grantham_distance'] == '' or row['grantham_distance'] is None:
               row['grantham_distance'] = -1
            else:
               row['grantham_distance'] = float(row['grantham_distance'])
         if not len(row['hgvsp']) == 0:
            oncogenic_variants[str(gene) + '-' + str(row['hgvsp'])] = row
         if not len(row['hgvs_c']) == 0:
            oncogenic_variants[str(gene) + '-' + str(row['hgvs_c'])] = row
         if not len(row['molecular_consequence']) == 0 and row['molecular_consequence'] == 'missense_variant':
            key = str(gene) + '-oncogenic_codon-' + str(row['codon'])
            if key in oncogenic_variants:
               if row['grantham_distance'] < oncogenic_variants[key]['grantham_distance']:
                  oncogenic_variants[key] = row
            else:
               oncogenic_variants[key] = row
         
   return oncogenic_variants


def match_oncogenic_variants(transcript_csq_elements, oncogenic_variants, rec, principal_csq_properties):

   """
   Function that matches consequence entries from VEP (transcript_csq_elements) agains known oncogenic variants from ClinVar,
   using both genomic coordinate information, HGVSp and HGVSc information.
   """
   
   principal_hgvsp = principal_csq_properties['hgvsp']
   principal_hgvsc = principal_csq_properties['hgvsc'] 
   
   known_oncogenic_matches = {}
   known_oncogenic_sites = {}

   for csq in transcript_csq_elements:
      (consequence, symbol, entrezgene, hgvsc, hgvsp, exon, feature_type, feature, biotype) = csq.split(':')

      if bool(re.search(r'^(missense|stop|start|inframe|splice_donor|intron|splice_acceptor|frameshift|upstream)', consequence)) is False:
         continue

      hgvsc_match = 'by_hgvsc_nonprincipal'
      hgvsp_match = 'by_hgvsp_nonprincipal'
      codon_match = 'by_codon_nonprincipal'
      if hgvsc == principal_hgvsc:
         hgvsc_match = 'by_hgvsc_principal'
      if threeToOneAA(hgvsp) == principal_hgvsp:
         hgvsp_match = 'by_hgvsp_principal'
         codon_match = 'by_codon_principal'
      
      var_id = str(rec.CHROM) + '_' + str(rec.POS) + '_' + str(rec.REF) + '_' + str(','.join(rec.ALT))
      
      oncogenic_varkeys = {}
      if entrezgene != ".":
         oncogenic_varkeys[str(entrezgene) + '-' + str(var_id)] = 'by_genomic_coord'
         if hgvsp != ".":
            hgvsp_short = threeToOneAA(hgvsp)
            oncogenic_varkeys[str(entrezgene) + '-' + str(hgvsp_short)] = hgvsp_match
            if re.search(r'^missense', consequence) and not re.search(r'(del|ins|fs|Ter|X)', hgvsp_short):
               codon = re.match(r'p.[A-Z]{1}[0-9]{1,}', hgvsp_short)
               if codon:
                  oncogenic_varkeys[str(entrezgene) + '-oncogenic_codon-' + str(codon.group(0))] = codon_match
         if hgvsc != ".":
            oncogenic_varkeys[str(entrezgene) + '-' + str(hgvsc)] = hgvsc_match
      
      for oncogenic_varkey in oncogenic_varkeys:
         if oncogenic_varkey in oncogenic_variants and 'oncogenicity' in oncogenic_variants[oncogenic_varkey]:    
            oncogenic_info = 'clinvar|' + oncogenic_variants[oncogenic_varkey]['symbol'] + '|' + \
               str(oncogenic_variants[oncogenic_varkey]['hgvsp']) + '|' + \
               str(oncogenic_variants[oncogenic_varkey]['hgvs_c']) + '|' + \
               str(oncogenic_variants[oncogenic_varkey]['grantham_distance']) + '|' + \
               str(oncogenic_variants[oncogenic_varkey]['oncogenicity'])

            if oncogenic_varkeys[oncogenic_varkey].startswith('by_codon'):
               grantham_distance = rec.INFO.get('GRANTHAM_DISTANCE')
               if grantham_distance is not None:
                  if grantham_distance > 0 and oncogenic_variants[oncogenic_varkey]['grantham_distance'] > 0:                     
                     if float(grantham_distance / float(oncogenic_variants[oncogenic_varkey]['grantham_distance'])) >= 0.8:                  
                        if oncogenic_info not in known_oncogenic_sites:
                           known_oncogenic_sites[oncogenic_info] = []
                        known_oncogenic_sites[oncogenic_info].append(oncogenic_varkeys[oncogenic_varkey])
            else:
               if oncogenic_info not in known_oncogenic_matches:
                  known_oncogenic_matches[oncogenic_info] = []
               known_oncogenic_matches[oncogenic_info].append(oncogenic_varkeys[oncogenic_varkey])

   
   if len(list(known_oncogenic_matches.keys())) >= 1:
      all_matches = []
      for oncogenic_info in list(known_oncogenic_matches.keys()):
         all_matches.append(oncogenic_info + '|' + \
            '&'.join(sorted(set(known_oncogenic_matches[oncogenic_info]))))
      rec.INFO['KNOWN_ONCOGENIC'] = ','.join(all_matches)
   
   if len(list(known_oncogenic_sites.keys())) >= 1:
      all_matches = []
      for oncogenic_info in list(known_oncogenic_sites.keys()):
         all_matches.append(oncogenic_info + '|' + \
            '&'.join(sorted(set(known_oncogenic_sites[oncogenic_info]))))
      rec.INFO['KNOWN_ONCOGENIC_SITE'] = ','.join(all_matches)

   return


def _sort_oncogenicity_codes(code_str: str) -> str:
   """Return the pipe-separated ONCOGENICITY_CODE string in canonical order.
   Codes absent from ONCOGENICITY_CODE_ORDER are appended at the end."""
   if code_str == '.':
      return code_str
   _order_index = {c: i for i, c in enumerate(ONCOGENICITY_CODE_ORDER)}
   codes = [c for c in code_str.split('|') if c]
   codes.sort(key=lambda c: _order_index.get(c, len(ONCOGENICITY_CODE_ORDER)))
   return '|'.join(codes)


def _classify_oncogenicity_score(score: float, code_str: str) -> str:
   """
   Map a numeric oncogenicity score + evidence code string to a classification
   label using the VICC/ClinGen thresholds defined in ONCOGENICITY_THRESHOLDS.
   Mirrors the logic in assign_oncogenicity_evidence, including the edge-case
   for variants one point below the Likely_Oncogenic threshold.
   """
   t = ONCOGENICITY_THRESHOLDS
   classification = "VUS"

   if t['likely_benign_lower'] <= score <= t['likely_benign_upper']:
      classification = "Likely_Benign"
   if t['likely_oncogenic_lower'] <= score <= t['likely_oncogenic_upper']:
      classification = "Likely_Oncogenic"
   # Edge case: one point below LP threshold — still LP when ≥3 positive criteria
   # present (without any benign criteria) or when ONCG_OS1 is in the code string
   if score == t['likely_oncogenic_lower'] - 1 and \
      re.search(r'_SB', code_str) is None and \
      (('|' in code_str and len(code_str.split('|')) >= 3) or
       re.search(r'ONCG_OS1', code_str)):
      classification = "Likely_Oncogenic"
   if score <= t['benign_upper']:
      classification = "Benign"
   if score >= t['oncogenic_lower']:
      classification = "Oncogenic"

   return classification


def refine_oncogenicity_with_oncokb(
      tsv_gz_fname: str,
      logger: logging.Logger = None) -> None:
   """
   Second-pass oncogenicity refinement using OncoKB annotations already
   present as *_OKB columns in the final variant TSV.

   Three groups of OKB-derived criteria are injected where not already captured
   by the primary VICC/ClinGen classification:

   Oncogenic evidence (positive score) — from ONCOGENICITY_OKB:
     ONCG_OS2_A  (+4) — "Oncogenic"        (strong positive)
     ONCG_OS2_B  (+2) — "Likely Oncogenic" (moderate positive)

   LoF evidence in TSGs (positive score) — from MUTATION_EFFECT_OKB:
     ONCG_OVS1_A (+8) — "Loss-of-function"        in TSG (very strong; only if ONCG_OVS1 not set)
     ONCG_OVS1_B (+4) — "Likely Loss-of-function" in TSG (strong;      only if ONCG_OVS1 not set)

   Benign / neutral evidence (negative score) — from ONCOGENICITY_OKB:
     ONCG_SBS2_A  (−4) — "Neutral"        (strong benign)
     ONCG_SBS2_B  (−2) — "Likely Neutral" (moderate benign)

   After injection the score and classification label are recomputed.
   The TSV is overwritten in place.

   Args:
      tsv_gz_fname: Path to the gzipped variant TSV (read + overwritten)
      logger:       Logger instance
   """
   if logger is None:
      logger = logging.getLogger("pcgr-oncogenicity-refine-okb")

   if not os.path.exists(tsv_gz_fname):
      logger.warning(f"Variant TSV not found: {tsv_gz_fname} - skipping OncoKB oncogenicity refinement")
      return

   tsv_df = pd.read_csv(tsv_gz_fname, sep="\t", low_memory=False)

   required = {'ONCOGENICITY', 'ONCOGENICITY_CODE', 'ONCOGENICITY_SCORE',
               'ONCOGENICITY_OKB', 'MUTATION_EFFECT_OKB', 'TSG', 'VAR_ID'}
   missing = required - set(tsv_df.columns)
   if missing:
      logger.warning(f"Missing columns for OncoKB oncogenicity refinement: {missing} - skipping")
      return

   n_refined = 0

   for idx, row in tsv_df.iterrows():
      existing_code   = str(row['ONCOGENICITY_CODE'])   if pd.notna(row['ONCOGENICITY_CODE'])   else '.'
      existing_score  = float(row['ONCOGENICITY_SCORE']) if pd.notna(row['ONCOGENICITY_SCORE']) else 0.0
      oncogenicity_okb    = str(row['ONCOGENICITY_OKB'])    if pd.notna(row['ONCOGENICITY_OKB'])    else None
      mutation_effect_okb = str(row['MUTATION_EFFECT_OKB'])  if pd.notna(row['MUTATION_EFFECT_OKB'])  else None
      is_tsg = row['TSG'] is True or str(row['TSG']).lower() == 'true'
      var_id = str(row['VAR_ID']) if pd.notna(row['VAR_ID']) else None

      if oncogenicity_okb is None and mutation_effect_okb is None:
         continue

      active_codes = set(existing_code.split('|')) if existing_code != '.' else set()
      new_codes = []
      score_delta = 0.0

      # --- Oncogenic evidence (ONCOGENICITY_OKB) ---
      # OS2_A: "Oncogenic" — strong. Skip if OS1 or either OS2 already present.
      if oncogenicity_okb == 'Oncogenic' and \
         'ONCG_OS1'   not in active_codes and \
         'ONCG_OS2_A' not in active_codes and \
         'ONCG_OS2_B' not in active_codes:
         new_codes.append('ONCG_OS2_A')
         score_delta += OKB_ONCOGENICITY_CRITERIA['ONCG_OS2_A']['score']

      # OS2_B: "Likely Oncogenic" — moderate. Skipped if stronger evidence already present.
      elif oncogenicity_okb == 'Likely Oncogenic' and \
           'ONCG_OS1'   not in active_codes and \
           'ONCG_OS2_A' not in active_codes and \
           'ONCG_OS2_B' not in active_codes:
         new_codes.append('ONCG_OS2_B')
         score_delta += OKB_ONCOGENICITY_CRITERIA['ONCG_OS2_B']['score']

      # --- Benign / neutral evidence (ONCOGENICITY_OKB) ---
      # SBS2_A: "Neutral" — strong benign.
      if oncogenicity_okb == 'Neutral' and \
         'ONCG_SBS2_A' not in active_codes and \
         'ONCG_SBS2_B' not in active_codes:
         new_codes.append('ONCG_SBS2_A')
         score_delta += OKB_ONCOGENICITY_CRITERIA['ONCG_SBS2_A']['score']

      # SBS2_B: "Likely Neutral" — moderate benign. Skipped if stronger evidence present.
      elif oncogenicity_okb == 'Likely Neutral' and \
           'ONCG_SBS2_A' not in active_codes and \
           'ONCG_SBS2_B' not in active_codes:
         new_codes.append('ONCG_SBS2_B')
         score_delta += OKB_ONCOGENICITY_CRITERIA['ONCG_SBS2_B']['score']

      # --- LoF evidence in TSGs (MUTATION_EFFECT_OKB) ---
      # Only applied when the gene is a TSG and PCGR's internal OVS1 was not set,
      # covering cases where the internal LoF predictor missed a functionally
      # characterised variant that OncoKB has curated.
      if is_tsg and 'ONCG_OVS1' not in active_codes and \
         'ONCG_OVS1_A' not in active_codes and \
         'ONCG_OVS1_B' not in active_codes:

         # OVS1_A: "Loss-of-function" — very strong (OVS1-equivalent).
         if mutation_effect_okb == 'Loss-of-function':
            new_codes.append('ONCG_OVS1_A')
            score_delta += OKB_ONCOGENICITY_CRITERIA['ONCG_OVS1_A']['score']

         # OVS1_B: "Likely Loss-of-function" — strong (OS-equivalent).
         elif mutation_effect_okb == 'Likely Loss-of-function':
            new_codes.append('ONCG_OVS1_B')
            score_delta += OKB_ONCOGENICITY_CRITERIA['ONCG_OVS1_B']['score']

      if not new_codes:
         continue

      ## log the update for traceability
      #logger.info(f"Refining oncogenicity for variant (TSG = {is_tsg}) at index {var_id} with new "
      #             f"evidence codes {new_codes} (existing {existing_code}) and score delta {score_delta:.1f}")
      

      updated_code  = existing_code + '|' + '|'.join(new_codes) if existing_code != '.' \
                      else '|'.join(new_codes)
      updated_code  = _sort_oncogenicity_codes(updated_code)
      updated_score = existing_score + score_delta
      updated_class = _classify_oncogenicity_score(updated_score, updated_code)

      tsv_df.at[idx, 'ONCOGENICITY_CODE']  = updated_code
      tsv_df.at[idx, 'ONCOGENICITY_SCORE'] = updated_score
      tsv_df.at[idx, 'ONCOGENICITY']       = updated_class
      n_refined += 1

   tsv_df.to_csv(tsv_gz_fname, sep="\t", compression="gzip", index=False)
   logger.info(f"OncoKB oncogenicity refinement complete: {n_refined} / {len(tsv_df)} variants updated")