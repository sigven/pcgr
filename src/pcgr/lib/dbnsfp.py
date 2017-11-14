#!/usr/bin/env python

import os,re,sys
import logging
import csv
from bx.intervals.intersection import IntervalTree

def map_effect_abbreviation(abbrev, algo):
   
   if algo == 'splice_site_ada' or algo == 'splice_site_rf' or algo == 'splice_site_global' or algo == 'cadd_phred' or algo == 'gerp_rs' or algo == 'mutpred':
      return abbrev
   if abbrev == 'D':
      #if algo == 'polyphen2_hdiv' or algo == 'polyphen2_hvar':
         #return 'probably_damaging'
      #else:
      return 'damaging'
   if abbrev == 'P':
      return 'possibly_damaging'
   if abbrev == 'H' or abbrev == 'M' or abbrev == 'A':
      return 'damaging'
      if abbrev == 'M' and algo == 'mutationassessor':
         return 'possibly_damaging'
   if abbrev == 'T' or abbrev == 'N' or abbrev == 'L' or abbrev == 'B':
      return 'tolerated'
   if abbrev == 'U':
      return '.'
   
   return '.' 

def map_dbnsfp_predictions(dbnsfp_tag, algorithms):
   
   effect_predictions = {}
   
   for v in dbnsfp_tag.split(','):
   
      dbnsfp_info = v.split('|')
      if len(dbnsfp_info) == 1:
         dbnsfp_info = v.split('#')
      ref_aa = dbnsfp_info[0]
      alt_aa = dbnsfp_info[1]
      all_ids = dbnsfp_info[4].split('&')
      unique_ids = {}
      for s in all_ids:
         unique_ids[s] = 1
         
      isoform_aa_keys = []
      if ref_aa != '.' and alt_aa != '.':
         aa_pos = dbnsfp_info[6].split('&')
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
   
      i = 7
      v = 0
      
      if len(algorithms) != len(dbnsfp_info[7:]):
         return effect_predictions
      
      while i < len(dbnsfp_info):
         algorithm_raw_predictions[str(algorithms[v]).lower()] = dbnsfp_info[i].split('&')
         i = i + 1
         v = v + 1
      dbnsfp_predictions = {}
      
      for k in isoform_aa_keys:
         if not dbnsfp_predictions.has_key(k):
            dbnsfp_predictions[k] = {}
         all_preds = []
         for algo in algorithm_raw_predictions.keys():
            unique_algo_predictions = {}
            for pred in algorithm_raw_predictions[algo]:
               if pred != '':
                  if not unique_algo_predictions.has_key(map_effect_abbreviation(pred,algo)):
                     unique_algo_predictions[map_effect_abbreviation(pred,algo)] = 1
               else:
                  unique_algo_predictions['.'] = 1
            
            if len(unique_algo_predictions.keys()) > 1 and '.' in unique_algo_predictions.keys():
               del unique_algo_predictions['.']
            dbnsfp_predictions[k][algo] = str(algo) + ':' + '|'.join(unique_algo_predictions.keys())  
            all_preds.append(dbnsfp_predictions[k][algo])
         effect_predictions[k] = '&'.join(all_preds)
   
   return effect_predictions
