#!/usr/bin/env python

import os,re,sys
from cyvcf2 import VCF, Writer
import csv

from typing import Dict
from logging import Logger
import gzip

from pcgr.annoutils import threeToOneAA

def assign_oncogenicity_evidence(rec = None, oncogenic_variants = None, tumortype = "Any"):

   clingen_vicc_ev_codes = [
      "ONCG_SBVS1", 
      "ONCG_SBS1", 
      "ONCG_SBP1", 
      "ONCG_SBP2",
      "ONCG_OS1", 
      "ONCG_OS2",
      "ONCG_OS3", 
      "ONCG_OM1",
      "ONCG_OM2",
      "ONCG_OM3",  
      "ONCG_OM4", 
      "ONCG_OP1", 
      "ONCG_OP3",
      "ONCG_OP4",
      "ONCG_OVS1"]

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
   ## oncogenic effect of the variant - CIViC/CGI

   # 9) "ONCG_OS3"
   ## Located in a mutation hotspot
   ## - >= 50 samples with a somatic variant at the AA position (cancerhotspots.org)
   ## - Same amino acid change in >= 10 samples (cancerhotspots.org)

   # 10) "ONCG_OM1"
   ## Located in a critical and well-established part of a functional domain
   ## (active site of an enzyme)
   ## NOTE: Here we assume that OM1 is applicable if there is overlap with oncogenic/predictive evidence (CIViC/CGI)

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
      "TSG",
      "TSG_EVIDENCE",
      "LOSS_OF_FUNCTION",
      "INTRON_POSITION",
      "EXON_POSITION",
      "CODING_STATUS",
      "gnomADe_EAS_AF",
      "gnomADe_NFE_AF",
      "gnomADe_AFR_AF",
      "gnomADe_AMR_AF",
      "gnomADe_SAS_AF",
      "DBNSFP_SIFT",
      "DBNSFP_PROVEAN",
      "DBNSFP_META_RNN",
      "DBNSFP_FATHMM",
      "DBNSFP_MUTATIONTASTER",
      "DBNSFP_DEOGEN2",
      "DBNSFP_PRIMATEAI",
      "DBNSFP_MUTATIONASSESSOR",
      "DBNSFP_FATHMM_MKL",
      "DBNSFP_M_CAP",
      "DBNSFP_LIST_S2",
      "DBNSFP_BAYESDEL_ADDAF",
      "DBNSFP_SPLICE_SITE_ADA",
      "DBNSFP_SPLICE_SITE_RF"]

   variant_data = {}
   for col in required_oncogenicity_vars:
      if rec.INFO.get(col) is None:
         if col == "TSG" or col == "ONCOGENE":
            variant_data[col] = False
         elif col == "LOSS_OF_FUNCTION":
            variant_data['LOSS_OF_FUNCTION'] = False
         else:
            variant_data[col] = None
      else:
         if rec.INFO.get(col) == '':
            variant_data[col] = True
         else: 
            variant_data[col] = rec.INFO.get(col)
   
   for code in clingen_vicc_ev_codes:
      variant_data[code] = False
   
   if "KNOWN_ONCOGENIC" in variant_data.keys():
      if not variant_data['KNOWN_ONCOGENIC'] is None:         
         variant_data['ONCG_OS1'] = True
   
   dbnsfp_minimum_majority = 8
   dbnsfp_maximum_minority = 2
   dbnsfp_minimum_algos_called = dbnsfp_minimum_majority

   variant_data['N_INSILICO_CALLED'] = 0
   variant_data['N_INSILICO_DAMAGING'] = 0
   variant_data['N_INSILICO_TOLERATED'] = 0
   variant_data['N_INSILICO_SPLICING_NEUTRAL'] = 0
   variant_data['N_INSILICO_SPLICING_AFFECTED'] = 0

   for algo in ['SIFT',
                'PROVEAN',
                'META_RNN',
                'ALPHA_MISSENSE',
                'POLYPHEN_HVAR',
                'PHACTBOOOST',
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
         if (col == 'SPLICE_SITE_RF' or \
             col == 'SPLICE_SITE_ADA') and \
            variant_data[col] == 'SN':
            variant_data['N_INSILICO_SPLICING_NEUTRAL'] += 1
         if (col == 'SPLICE_SITE_RF' or \
             col == 'SPLICE_SITE_ADA') and \
            variant_data[col] == 'AS':
            variant_data['N_INSILICO_SPLICING_AFFECTED'] += 1

   if (variant_data['N_INSILICO_CALLED'] > dbnsfp_minimum_algos_called and \
      variant_data['N_INSILICO_DAMAGING'] >= dbnsfp_minimum_majority and \
      variant_data['N_INSILICO_TOLERATED'] <= dbnsfp_maximum_minority and \
      variant_data['N_INSILICO_SPLICING_NEUTRAL'] <= 1) or \
         variant_data['N_INSILICO_SPLICING_AFFECTED'] == 2:
      variant_data['ONCG_OP1'] = True
   
   if variant_data['N_INSILICO_CALLED'] > dbnsfp_minimum_algos_called and \
      variant_data['N_INSILICO_TOLERATED'] >= dbnsfp_minimum_majority and \
      variant_data['N_INSILICO_DAMAGING'] <= dbnsfp_maximum_minority and \
      variant_data['N_INSILICO_SPLICING_AFFECTED'] == 0:
      variant_data['ONCG_SBP1'] = True

   hotspot_mutated_samples = {}
   hotspot_mutated_samples['aa_variant'] = {}
   hotspot_mutated_samples['aa_site'] = {}
   hotspot_mutated_samples['aa_variant']['Any'] = 0
   hotspot_mutated_samples['aa_site']['Any'] = 0

   if not variant_data['MUTATION_HOTSPOT_CANCERTYPE'] is None:
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
   ## we rather base our criteria matching based on existing actionability evidence (predictive/oncogenic), for which 
   ## variants are presumably located in critical sites of functional domains (in that sense, indirect evidence for OM1)
   if not variant_data['BIOMARKER_MATCH'] is None:
      
      ## Split all biomarker evidence into a list
      biomarker_evidence = variant_data['BIOMARKER_MATCH'].split(',')
     
      for eitem in biomarker_evidence:
                  
         ## Example 'eitem' element:
         ## cgi|659|CGI1077:Pancreas:Sensitivity/Response:C:Predictive:Somatic|by_hgvsp_principal
         if ('Predictive' in eitem or 'Oncogenic' in eitem) and \
            'Somatic' in eitem and \
            ('by_genomic_coord' in eitem or \
             'by_hgvsc_principal' in eitem or \
             'by_hgvsp_principal' in eitem or \
             'by_codon_principal' in eitem or \
             'by_aa_region_principal' in eitem):               
               ## only applicable if OS3 is not set
               if variant_data['ONCG_OS3'] is False and variant_data['ONCG_OS1'] is False:
                  variant_data['ONCG_OM1'] = True
         
         ## Catch prognostic/diagnostic non-coding variants (e.g. TERT) - these will rank at the top
         ## with respect to oncogenicity (altough not classified as oncogenic in the strict sense)
         if ('Diagnostic' in eitem or 'Prognostic' in eitem or 'Predictive' in eitem or 'Oncogenic' in eitem) and \
            variant_data['CODING_STATUS'] == 'noncoding' and \
            'Somatic' in eitem and \
            ('by_genomic_coord' in eitem or \
             'by_hgvsc_principal' in eitem):
               ## only applicable if OS3 is not set
               if variant_data['ONCG_OS3'] is False and variant_data['ONCG_OS1'] is False:
                  variant_data['ONCG_OM1'] = True

   if "gnomADe_EAS_AF" in variant_data.keys() and \
      "gnomADe_SAS_AF" in variant_data.keys() and \
      "gnomADe_AMR_AF" in variant_data.keys() and \
      "gnomADe_AFR_AF" in variant_data.keys() and \
      "gnomADe_NFE_AF" in variant_data.keys():

      ## check if variant has MAF > 0.01 (SBVS1) or > 0.05 in any of five major gnomAD subpopulations (exome set)
      for pop in ['gnomADe_SAS_AF','gnomADe_EAS_AF','gnomADe_AMR_AF','gnomADe_AFR_AF','gnomADe_NFE_AF']:
         if not variant_data[pop] is None:
            ## MAF for this population >= 0.05
            if float(variant_data[pop]) >= 0.05:
               variant_data["ONCG_SBVS1"] = True
      for pop in ['gnomADe_SAS_AF','gnomADe_EAS_AF','gnomADe_AMR_AF','gnomADe_AFR_AF','gnomADe_NFE_AF']:
         if not variant_data[pop] is None:
            ## MAF for this population >= 0.01 (< 0.05)
            if float(variant_data[pop]) >= 0.01 and variant_data["ONCG_SBVS1"] is False:
               variant_data["ONCG_SBS1"] = True

      #missing_pop_freq = 0
      approx_zero_pop_freq = 0
      for pop in ['gnomADe_SAS_AF','gnomADe_EAS_AF','gnomADe_AMR_AF','gnomADe_AFR_AF','gnomADe_NFE_AF']:
         ## no MAF recorded in gnomAD for this population
         if variant_data[pop] is None:
            approx_zero_pop_freq = approx_zero_pop_freq + 1
         else:
            ## Very low MAF for this population
            if float(variant_data[pop]) < 0.0001:
               approx_zero_pop_freq = approx_zero_pop_freq + 1
    
      ## check if variant is missing or with MAF approximately zero in all five major gnomAD subpopulations (exome set)
      if approx_zero_pop_freq == 5:
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
      if not variant_data['KNOWN_ONCOGENIC_SITE'] is None:         
         if variant_data['ONCG_OS1'] is False and \
            variant_data['ONCG_OS3'] is False and \
               variant_data['ONCG_OM1'] is False:
            variant_data['ONCG_OM3'] = True

   og_score_data = {}
   og_score_data['code'] = \
      ['ONCG_SBVS1',
       'ONCG_SBS1',
       'ONCG_OP4',
       'ONCG_SBP1',
       'ONCG_SBP2',
       'ONCG_OVS1',
       'ONCG_OS1',
       'ONCG_OS3',
       'ONCG_OM1',
       'ONCG_OM2',
       'ONCG_OM3',
       'ONCG_OM4',
       'ONCG_OP1',
       'ONCG_OP3'
      ]
   og_score_data['category'] = \
      ['clinpop',
       'clinpop', 
       'clinpop',
       'funccomp',
       'funcvar',
       'funcvar',
       'funcvar',
       'funcvar',
       'funcvar',
       'funcvar',
       'funcvar',
       'funcvar',
       'funcvar',
       'funccomp']
   og_score_data['pole'] = \
      ['B','B',
       'P','B',
       'B','P',
       'P','P',
       'P','P',
       'P','P',
       'P','P']
   
   og_score_data['description'] = \
      ['Very high MAF (> 0.05 in gnomAD - any five major continental pops)',
      'High MAF (> 0.01 in gnomAD - any five major continental pops)',
      'Absent from controls (gnomAD) / very low MAF ( < 0.0001 in all five major subpopulations)',
      'Multiple lines (>=8) of computational evidence support a benign variant effect on the gene or gene product - from dbNSFP',
      'Silent and intronic changes outside of the consensus splice site',
      'Null variant - predicted as LoF - in bona fide tumor suppressor gene',
      'Same amino acid change as previously established oncogenic variant - regardless of nucleotide change (ClinVar oncogenicity records)',
      'Located in a mutation hotspot (cancerhotspots.org). >= 50 samples with a  variant at AA position, >= 10 samples with same AA change',
      'Presumably critical site of functional domain - based on indirect evidence from overlap with predictive biomarkers',
      'Protein length changes from in-frame dels/ins in known oncogene/tumor suppressor genes or stop-loss variants in a tumor suppressor gene',
      'Missense variant at an amino acid residue where a different missense variant determined to be oncogenic (using this standard) has been documented. ',
      'Located in a mutation hotspot (cancerhotspots.org). < 50 samples with a variant at AA position, >= 10 samples with same AA change.',
      'Multiple lines (>=8) of computational evidence support of a damaging variant effect on the gene or gene product - from dbNSFP',
      'Located in a mutation hotspot (cancerhotspots.org). < 10 samples with the same amino acid change.'
      ]
   
   og_score_data['score'] = \
      [-8, -4, 1, -1, -1, 8, 4, 4, 2, 2, 2, 2, 1, 1]
   
   i = 0
   oncogenicity_scores = {}
   for code in og_score_data['code']:
      oncogenicity_scores[code] = {}
      for e in ['category','pole','description','score']:
         oncogenicity_scores[code][e] = og_score_data[e][i]
         if e == 'score':
            oncogenicity_scores[code][e] = float(og_score_data[e][i])
      i = i + 1

   variant_data['ONCOGENICITY'] = "VUS"
   variant_data["ONCOGENICITY_DOC"] = "."
   variant_data["ONCOGENICITY_CODE"] = "."
   variant_data["ONCOGENICITY_SCORE"] = "."
   onc_score_pathogenic = 0
   onc_score_benign = 0

   for code in clingen_vicc_ev_codes:

      if code in variant_data.keys():
         if variant_data[code] is True:
            score = oncogenicity_scores[code]['score']
            pole = oncogenicity_scores[code]['pole']
            if pole == "P":
               onc_score_pathogenic = onc_score_pathogenic + score
            if pole == "B":
               onc_score_benign = onc_score_benign + score

            for e in ['ONCOGENICITY_DOC',
                      'ONCOGENICITY_CODE']:

               if variant_data[e] == ".":
                  if e == 'ONCOGENICITY_DOC':
                     variant_data[e] = oncogenicity_scores[code]['description']
                  else:
                     variant_data[e] = code
               else:
                  if e == 'ONCOGENICITY_DOC':
                     variant_data[e] = f'{variant_data[e]}|{oncogenicity_scores[code]["description"]}'
                  else:
                     variant_data[e] = f'{variant_data[e]}|{code}'

   likely_benign_upper_limit = -1
   likely_benign_lower_limit = -6
   benign_upper_limit = -7
   likely_oncogenic_lower_limit = 5
   likely_oncogenic_upper_limit = 8
   oncogenic_lower_limit = 9

   variant_data['ONCOGENICITY_SCORE'] = onc_score_benign + onc_score_pathogenic
   if variant_data['ONCOGENICITY_SCORE'] <= likely_benign_upper_limit and \
      variant_data['ONCOGENICITY_SCORE'] >= likely_benign_lower_limit:
      variant_data['ONCOGENICITY'] = "Likely_Benign"
   if variant_data['ONCOGENICITY_SCORE'] >= likely_oncogenic_lower_limit and \
      variant_data['ONCOGENICITY_SCORE'] <= likely_oncogenic_upper_limit:
      variant_data['ONCOGENICITY'] = "Likely_Oncogenic"
   if variant_data['ONCOGENICITY_SCORE'] <= benign_upper_limit:
      variant_data['ONCOGENICITY'] = "Benign"
   if variant_data['ONCOGENICITY_SCORE'] >= oncogenic_lower_limit:
      variant_data['ONCOGENICITY'] = "Oncogenic"
   
   for e in ['ONCOGENICITY_SCORE',
             'ONCOGENICITY',
             'ONCOGENICITY_CODE']:
      rec.INFO[e] = variant_data[e]

   return(rec)

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
         oncogenic_variants[str(gene) + '-' + str(row['var_id'])] = row
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

      if not bool(re.search(r'^(missense|stop|start|inframe|splice_donor|intron|splice_acceptor|frameshift|upstream)', consequence)) is True:
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
               if not grantham_distance is None:
                  if float(grantham_distance / float(oncogenic_variants[oncogenic_varkey]['grantham_distance'])) > 0.8:                  
                     if not oncogenic_info in known_oncogenic_sites:
                        known_oncogenic_sites[oncogenic_info] = []
                     known_oncogenic_sites[oncogenic_info].append(oncogenic_varkeys[oncogenic_varkey])
            else:
               if not oncogenic_info in known_oncogenic_matches:
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