#!/usr/bin/env python

import os,re,sys
from cyvcf2 import VCF, Writer

def assign_oncogenicity_evidence(rec = None, tumortype = "Any"):

   clingen_vicc_ev_codes = [
      "CLINGEN_VICC_SBVS1", 
      "CLINGEN_VICC_SBS1", 
      "CLINGEN_VICC_SBP1", 
      "CLINGEN_VICC_SBP2",
      "CLINGEN_VICC_OS3", 
      "CLINGEN_VICC_OM2",
      "CLINGEN_VICC_OM3",  
      "CLINGEN_VICC_OM4", 
      "CLINGEN_VICC_OP1", 
      "CLINGEN_VICC_OP3",
      "CLINGEN_VICC_OP4",
      "CLINGEN_VICC_OVS1"]

   ### Benign oncogenic effects of somatic variants
   
   #  1) "CLINGEN_VICC_SBVS1"
   ## Very high MAF: > 0.05 in gnomAD - any 5 general continental populations
   ## AFR/AMR/EAS/NFE/SAS

   # 2) "CLINGEN_VICC_SBS1"
   ## High MAF: > 0.01 in gnomAD - any 5 general continental populations
   ## AFR/AMR/EAS/NFE/SAS

   # 3) "CLINGEN_VICC_SBP1"
   ## Multiple lines of computational evidence support a benign
   ## effect on the gene or gene product
   ## (conservation, evolutionary, splicing impact, etc. - from dbNSFP

   # 4) CLINGEN_VICC_SBP2"
   ## Synonymous (silent) variant for which splicing prediction
   ## algorithms predict no effect on the splice consensus sequence
   ## nor the creation of a new splice site and the nucleotide is
   ## not highly conserved

   # 5) "CLINGEN_VICC_SBS2"
   ## Well established in invitro/in vivo functional studies show
   ## no oncogenic effects

   ## Oncogenic effects of somatic variants
   
   # 6) "CLINGEN_VICC_OVS1"
   ## Null variant - predicted as LoF by LOFTEE
   ## - Nonsense, frameshift, canonical splice sites, initiation codon,
   ##   single-exon/multi-exon deletion
   ## - Tumor suppressor gene

   # 7) "CLINGEN_VICC_OS1"
   ## Same amino acid change as previously established oncogenic variant
   ## NOTE: Not relevant since current implementation does not consider
   ## any criteria where nucleotide change is affecting annotation

   # 8) "CLINGEN_VICC_OS2"
   ## Well established in invitro/in vivo functional studies show
   ## oncogenic effect of the variant

   # 9) "CLINGEN_VICC_OS3"
   ## Located in a mutation hotspot
   ## - >= 50 samples with a somatic variant at the AA position (cancerhotspots.org)
   ## - Same amino acid change in >= 10 samples (cancerhotspots.org)

   # 10) "CLINGEN_VICC_OM1"
   ## Located in a critical and well-established part of a functional domain
   ## (active site of an enzyme)
   ## NOTE: Not used since determination/retrieval of critical protein sites are non-trivial to automate

   # 11) "CLINGEN_VICC_OM2"
   ## Protein length changes as a result of in-frame deletions/insertions in a
   ## known oncogene/tumor suppressor genes or stop-loss variants in a
   ## tumor suppressor gene

   # 12) "CLINGEN_VICC_OM3"
   ## Missense variant at an amino acid residue where a different missense
   ## variant determined to be oncogenic (using this standard) has been
   ## documented. Amino acid difference from reference amino acid should
   ## be greater or at least approximately the same as for missense change
   ## determined to be oncogenic.

   # 13) "CLINGEN_VICC_OM4"
   ## Located in a mutation hotspot
   ## - < 50 samples with a somatic variant at the AA position (cancerhotspots.org)
   ## - Same amino acid change in >= 10 samples (cancerhotspots.org)
   ## - Not applicable if OM1 or OM3 is applicable

   # 14) "CLINGEN_VICC_OP1"
   ## All used lines of computational support an oncogenic effect
   ## of a variant (conservation, evolutionary, splicing effect)

   # 15) "CLINGEN_VICC_OP2"
   ## Somatic variant in a gene in a malignancy with a single genetic
   ## etiology. Example:retinoblastoma is caused by bi-allelic
   ## RB1 inactivation.

   # 16) "CLINGEN_VICC_OP3"
   ## Located in a mutation hotspot
   ## - Same amino acid change in < 10 samples (cancerhotspots.org)
   ## - Not applicable if OM1 or OM3 is applicable

   # 17) "CLINGEN_VICC_OP4"
   ## Absent from controls (gnomAD)
   ## - Extremely low MAF


   required_oncogenicity_vars = [
      "Consequence",
      "MUTATION_HOTSPOT",
      "MUTATION_HOTSPOT_CANCERTYPE",
      "SYMBOL",
      "ONCOGENE",
      "ONCOGENE_EVIDENCE",
      "TSG",
      "TSG_EVIDENCE",
      "LOSS_OF_FUNCTION",
      "INTRON_POSITION",
      "EXON_POSITION",
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
   
   dbnsfp_minimum_majority = 6
   dbnsfp_maximum_minority = 2
   dbnsfp_minimum_algos_called = dbnsfp_minimum_majority

   variant_data['N_INSILICO_CALLED'] = 0
   variant_data['N_INSILICO_DAMAGING'] = 0
   variant_data['N_INSILICO_TOLERATED'] = 0
   variant_data['N_INSILICO_SPLICING_NEUTRAL'] = 0
   variant_data['N_INSILICO_SPLICING_AFFECTED'] = 0

   for algo in ['SIFT','PROVEAN','META_RNN',
                'MUTATIONTASTER','DEOGEN2',
                'PRIMATEAI','MUTATIONASSESSOR',
                'FATHMM_MKL','M_CAP',
                'LIST_S2','BAYESDEL_ADDAF',
                'SPLICE_SITE_RF','SPLICE_SITE_ADA']:
      col = 'DBNSFP_' + str(algo)
      if col in variant_data.keys():
         if variant_data[col] != '.':
            variant_data['N_INSILICO_CALLED'] += 1
         if variant_data[col] == 'D':
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
      variant_data['CLINGEN_VICC_OP1'] = True
   
   if variant_data['N_INSILICO_CALLED'] > dbnsfp_minimum_algos_called and \
      variant_data['N_INSILICO_TOLERATED'] >= dbnsfp_minimum_majority and \
      variant_data['N_INSILICO_DAMAGING'] <= dbnsfp_maximum_minority and \
      variant_data['N_INSILICO_SPLICING_AFFECTED'] == 0:
      variant_data['CLINGEN_VICC_SBP1'] = True

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
   

   ttype = "Any"
   if tumortype != "Any":
      ttype = re.sub(r'/', '@', re.sub(r' ','_', tumortype))

   if ttype in hotspot_mutated_samples['aa_variant'].keys() and \
      ttype in hotspot_mutated_samples['aa_site'].keys():
      if hotspot_mutated_samples['aa_variant'][ttype] >= 10 and \
         hotspot_mutated_samples['aa_site'][ttype] >= 50:
            variant_data['CLINGEN_VICC_OS3'] = True
      if hotspot_mutated_samples['aa_variant'][ttype] >= 10 and \
         hotspot_mutated_samples['aa_site'][ttype] < 50:
            variant_data['CLINGEN_VICC_OM4'] = True
      if hotspot_mutated_samples['aa_variant'][ttype] < 10 and \
         hotspot_mutated_samples['aa_variant'][ttype] > 0:
            variant_data['CLINGEN_VICC_OP3'] = True
   

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
               variant_data["CLINGEN_VICC_SBVS1"] = True
      for pop in ['gnomADe_SAS_AF','gnomADe_EAS_AF','gnomADe_AMR_AF','gnomADe_AFR_AF','gnomADe_NFE_AF']:
         if not variant_data[pop] is None:
            ## MAF for this population >= 0.01 (< 0.05)
            if float(variant_data[pop]) >= 0.01 and variant_data["CLINGEN_VICC_SBVS1"] is False:
               variant_data["CLINGEN_VICC_SBS1"] = True

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
         variant_data["CLINGEN_VICC_OP4"] = True
  
    
   ## check if variant is a loss-of-function variant (LOFTEE) in a tumor suppressor gene (Cancer Gene Census/CancerMine)
   if "TSG" in variant_data.keys() and \
      "ONCOGENE" in variant_data.keys() and \
      "LOSS_OF_FUNCTION" in variant_data.keys() and \
      "Consequence" in variant_data.keys():

      if variant_data['LOSS_OF_FUNCTION'] is True and variant_data['TSG'] is True:
         variant_data['CLINGEN_VICC_OVS1'] = True
  
      ## check if variant is creating a stop-lost or protein-length change in oncogene/tumor suppressor genes
      if variant_data['CLINGEN_VICC_OVS1'] is False and \
         ((re.match(r'^(inframe_deletion|inframe_insertion)', variant_data['Consequence']) and \
            (variant_data['TSG'] is True or variant_data['ONCOGENE'] is True)) or \
         (re.match(r'^(stop_lost)', variant_data['Consequence']) and \
            variant_data['TSG'] is True)):
            variant_data['CLINGEN_VICC_OM2'] = True
   
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
            variant_data['CLINGEN_VICC_SBP2'] = True
         

   og_score_data = {}
   og_score_data['code'] = \
      ['CLINGEN_VICC_SBVS1','CLINGEN_VICC_SBS1',
      'CLINGEN_VICC_SBP1','CLINGEN_VICC_SBP2',
      'CLINGEN_VICC_OVS1','CLINGEN_VICC_OS3',
      'CLINGEN_VICC_OM2','CLINGEN_VICC_OM4',
      'CLINGEN_VICC_OP1','CLINGEN_VICC_OP3',
      'CLINGEN_VICC_OP4']
   og_score_data['category'] = \
      ['clinpop','clinpop',
      'funccomp','funcvar',
      'funcvar','funcvar',
      'funcvar','funcvar',
      'funccomp','funcvar',
      'clinpop']
   og_score_data['pole'] = \
      ['B','B','B','B',
      'P','P','P','P',
      'P','P','P']
   
   og_score_data['description'] = \
      ['Very high MAF (> 0.05 in gnomAD - any five major continental pops)',
      'High MAF (> 0.01 in gnomAD - any five major continental pops)',
      'Multiple lines (>=7) of computational evidence support a benign effect on the gene or gene product - from dbNSFP',
      'Silent and intronic changes outside of the consensus splice site',
      'Null variant - predicted as LoF by LOFTEE - in bona fide tumor suppressor gene',
      'Located in a mutation hotspot (cancerhotspots.org). >= 50 samples with a  variant at AA position, >= 10 samples with same AA change',
      'Protein length changes from in-frame dels/ins in known oncogene/tumor suppressor genes or stop-loss variants in a tumor suppressor gene',
      'Located in a mutation hotspot (cancerhotspots.org). < 50 samples with a variant at AA position, >= 10 samples with same AA change.',
      'Multiple lines (>=7) of computational evidence support a damaging effect on the gene or gene product - from dbNSFP',
      'Located in a mutation hotspot (cancerhotspots.org). < 10 samples with the same amino acid change.',
      'Absent from controls (gnomAD) / very low MAF ( < 0.0001 in all five major subpopulations)']
   
   og_score_data['score'] = \
      [-8, -4, -1, -1, 8, 4, 2, 2, 1, 1, 1]
   
   i = 0
   oncogenicity_scores = {}
   for code in og_score_data['code']:
      oncogenicity_scores[code] = {}
      for e in ['category','pole','description','score']:
         oncogenicity_scores[code][e] = og_score_data[e][i]
         if e == 'score':
            oncogenicity_scores[code][e] = float(og_score_data[e][i])
      i = i + 1

   variant_data['ONCOGENICITY_CLASSIFICATION'] = "VUS"
   variant_data["ONCOGENICITY_CLASSIFICATION_DOC"] = "."
   variant_data["ONCOGENICITY_CLASSIFICATION_CODE"] = "."
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

            for e in ['ONCOGENICITY_CLASSIFICATION_DOC',
                      'ONCOGENICITY_CLASSIFICATION_CODE']:

               if variant_data[e] == ".":
                  if e == 'ONCOGENICITY_CLASSIFICATION_DOC':
                     variant_data[e] = oncogenicity_scores[code]['description']
                  else:
                     variant_data[e] = code
               else:
                  if e == 'ONCOGENICITY_CLASSIFICATION_DOC':
                     variant_data[e] = f'{variant_data[e]}|{oncogenicity_scores[code]["description"]}'
                  else:
                     variant_data[e] = f'{variant_data[e]}|{code}'

   likely_benign_upper_limit = -1
   likely_benign_lower_limit = -6
   benign_upper_limit = -7
   #vus_lower_limit = 0
   #vus_upper_limit = 4
   likely_oncogenic_lower_limit = 5
   likely_oncogenic_upper_limit = 9
   oncogenic_lower_limit = 10

   variant_data['ONCOGENICITY_SCORE'] = onc_score_benign + onc_score_pathogenic
   if variant_data['ONCOGENICITY_SCORE'] <= likely_benign_upper_limit and \
      variant_data['ONCOGENICITY_SCORE'] >= likely_benign_lower_limit:
      variant_data['ONCOGENICITY_CLASSIFICATION'] = "Likely_Benign"
   if variant_data['ONCOGENICITY_SCORE'] >= likely_oncogenic_lower_limit and \
      variant_data['ONCOGENICITY_SCORE'] <= likely_oncogenic_upper_limit:
      variant_data['ONCOGENICITY_CLASSIFICATION'] = "Likely_Oncogenic"
   if variant_data['ONCOGENICITY_SCORE'] <= benign_upper_limit:
      variant_data['ONCOGENICITY_CLASSIFICATION'] = "Benign"
   if variant_data['ONCOGENICITY_SCORE'] >= oncogenic_lower_limit:
      variant_data['ONCOGENICITY_CLASSIFICATION'] = "Oncogenic"
   
   for e in ['ONCOGENICITY_SCORE',
             'ONCOGENICITY_CLASSIFICATION',
             'ONCOGENICITY_CLASSIFICATION_CODE']:
      rec.INFO[e] = variant_data[e]

   return(rec)