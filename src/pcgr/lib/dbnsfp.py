#!/usr/bin/env python

import os,re,sys
import logging
import csv
from bx.intervals.intersection import IntervalTree

def map_effect_abbreviation(abbrev, algo):
   

   if algo == 'splice_site_ada' or algo == 'splice_site_rf' or algo == 'splice_site_global' or algo == 'cadd_phred':
      return abbrev
   if abbrev == 'D':
      if algo == 'polyphen':
         return 'probably_damaging'
      else:
         return 'damaging'
   if abbrev == 'P':
      return 'possibly_damaging'
   if abbrev == 'H' or abbrev == 'M' or abbrev == 'A':
      return 'damaging'
   if abbrev == 'T' or abbrev == 'N' or abbrev == 'L' or abbrev == 'B':
      return 'tolerated'
   if abbrev == 'U':
      return '.'
   
   return '.' 

def map_dbnsfp_predictions(dbnsfp_tag):
   
   effect_predictions = {}
   
   for v in dbnsfp_tag.split(','):
   
      dbnsfp_info = v.split('@')
      if len(dbnsfp_info) != 23:
         continue
      ref_aa = dbnsfp_info[0]
      alt_aa = dbnsfp_info[1]
      all_ids = dbnsfp_info[4].split('___')
      unique_ids = {}
      for s in all_ids:
         unique_ids[s] = 1
         
      isoform_aa_keys = []
      if ref_aa != '.' and alt_aa != '.':
         aa_pos = dbnsfp_info[6].split('___')
         num_aapos = len(aa_pos)
         for pos in aa_pos:
            for gene_id in unique_ids:
               k = str(gene_id) + ':p.' + str(ref_aa) + pos + str(alt_aa)
               isoform_aa_keys.append(k)
      else:
         #continue
         for gene_id in unique_ids:
            isoform_aa_keys.append(gene_id)
   
      algorithm_raw_predictions = {}
   
      algorithm_raw_predictions['sift'] = dbnsfp_info[7].split('___')
      algorithm_raw_predictions['polyphen'] = dbnsfp_info[8].split('___')
      algorithm_raw_predictions['lrt'] = dbnsfp_info[10].split('___')
      algorithm_raw_predictions['mutationtaster'] = dbnsfp_info[11].split('___')
      algorithm_raw_predictions['mutationassessor'] = dbnsfp_info[12].split('___')
      algorithm_raw_predictions['fathmm'] = dbnsfp_info[13].split('___')
      algorithm_raw_predictions['provean'] = dbnsfp_info[14].split('___')
      algorithm_raw_predictions['cadd_phred'] = dbnsfp_info[15].split('___')
      algorithm_raw_predictions['fathmm_mkl'] = dbnsfp_info[16].split('___')
      algorithm_raw_predictions['meta_svm'] = dbnsfp_info[17].split('___')
      algorithm_raw_predictions['meta_lr'] = dbnsfp_info[18].split('___')
      algorithm_raw_predictions['splice_site_ada'] = dbnsfp_info[20].split('___')
      algorithm_raw_predictions['splice_site_rf'] = dbnsfp_info[21].split('___')
      algorithm_raw_predictions['splice_site_global'] = dbnsfp_info[22].split('___')
      dbnsfp_predictions = {}
      
      for k in isoform_aa_keys:
         if not dbnsfp_predictions.has_key(k):
            dbnsfp_predictions[k] = {}
         for algo in algorithm_raw_predictions.keys():
            unique_algo_predictions = {}
            for pred in algorithm_raw_predictions[algo]:
               if pred != '.':
                  if not unique_algo_predictions.has_key(map_effect_abbreviation(pred,algo)):
                     unique_algo_predictions[map_effect_abbreviation(pred,algo)] = 1
               else:
                  unique_algo_predictions['.'] = 1
            
            if len(unique_algo_predictions.keys()) > 1 and '.' in unique_algo_predictions.keys():
               del unique_algo_predictions['.']
            dbnsfp_predictions[k][algo] = '|'.join(unique_algo_predictions.keys())
         pred_string = dbnsfp_predictions[k]['sift'] + '&' + dbnsfp_predictions[k]['polyphen'] + '&' + dbnsfp_predictions[k]['lrt'] + '&' + dbnsfp_predictions[k]['mutationtaster'] + '&' + dbnsfp_predictions[k]['mutationassessor'] + '&' + dbnsfp_predictions[k]['fathmm'] + '&' + dbnsfp_predictions[k]['provean'] + '&' + dbnsfp_predictions[k]['fathmm_mkl'] + '&' + dbnsfp_predictions[k]['cadd_phred'] + '&' + dbnsfp_predictions[k]['meta_svm'] + '&' + dbnsfp_predictions[k]['meta_lr'] + '&' + dbnsfp_predictions[k]['splice_site_ada'] + '&' + dbnsfp_predictions[k]['splice_site_rf'] + '&' + str(dbnsfp_predictions[k]['splice_site_global'])
      
         effect_predictions[k] = pred_string
   
   return effect_predictions
