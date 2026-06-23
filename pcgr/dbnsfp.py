#!/usr/bin/env python

import re

from cyvcf2 import VCF
from pcgr import pcgr_vars

def map_variant_effect_predictors(rec, algorithms):

    dbnsfp_predictions = map_dbnsfp_predictions(
        str(rec.INFO.get('DBNSFP')), algorithms)
        
    if rec.INFO.get('Gene') is None or rec.INFO.get('Consequence') is None:
        return
    gene_id = str(rec.INFO.get('Gene'))
    consequence = str(rec.INFO.get('Consequence'))

    dbnsfp_key = ''

    found_key = 0
    if rec.INFO.get('HGVSp_short') is not None and not rec.INFO.get('HGVSp_short') == '.':
        aa_change = str(rec.INFO.get('HGVSp_short'))
        dbnsfp_key = gene_id + ':' + str(aa_change)
        if dbnsfp_key in dbnsfp_predictions:
            found_key = 1
    
    if found_key == 0 and re.search('splice_', consequence):
        dbnsfp_key = gene_id

    if dbnsfp_key != '':
        if dbnsfp_key in dbnsfp_predictions:
            rec.INFO['EFFECT_PREDICTIONS'] = dbnsfp_predictions[dbnsfp_key]
            for algo_pred in rec.INFO['EFFECT_PREDICTIONS'].split('&'):
                if algo_pred.split(':')[0] in pcgr_vars.DBNSFP_ALGORITHMS:
                    rec.INFO[pcgr_vars.DBNSFP_ALGORITHMS[algo_pred.split(':')[0]]] = str(algo_pred.split(':')[1])
                


def map_dbnsfp_predictions(dbnsfp_tag, algorithms):

    effect_predictions = {}
    
    for v in dbnsfp_tag.split(','):
        dbnsfp_info = v.split('|')
        if len(dbnsfp_info) == 1:
            return effect_predictions
        ref_aa = dbnsfp_info[0]
        alt_aa = dbnsfp_info[1]
        all_ids = dbnsfp_info[3].split('&')
        unique_ids = {}
        for s in all_ids:
            unique_ids[s] = 1

        isoform_aa_keys = []
        if ref_aa != '.' and alt_aa != '.' and ref_aa != '' and alt_aa != '':
            aa_pos = dbnsfp_info[5].split('&')
            for pos in aa_pos:
                for gene_id in unique_ids:
                    k = str(gene_id) + ':p.' + str(ref_aa) + pos + str(alt_aa)
                    isoform_aa_keys.append(k)
        else:
            # continue
            for gene_id in unique_ids:
                isoform_aa_keys.append(gene_id)

        algorithm_raw_predictions = {}

        i = 6
        v = 0
               
        if len(algorithms) != len(dbnsfp_info[6:]):            
            return effect_predictions

        while i < len(dbnsfp_info):
            algorithm_raw_predictions[str(
                algorithms[v]).lower()] = dbnsfp_info[i].split('&')
            i = i + 1
            v = v + 1
        dbnsfp_predictions = {}

        for k in isoform_aa_keys:
            if k not in dbnsfp_predictions:
                dbnsfp_predictions[k] = {}
            all_preds = []
            for algo in algorithm_raw_predictions.keys():
                unique_algo_predictions = {}
                for pred in algorithm_raw_predictions[algo]:
                    if pred != '':
                        if pred not in unique_algo_predictions:
                            unique_algo_predictions[pred] = 1
                    else:
                        unique_algo_predictions['.'] = 1
                if len(unique_algo_predictions.keys()) > 1 and '.' in unique_algo_predictions.keys():
                    del unique_algo_predictions['.']
                dbnsfp_predictions[k][algo] = str(
                    algo) + ':' + '|'.join(unique_algo_predictions.keys())
                all_preds.append(dbnsfp_predictions[k][algo])
            effect_predictions[k] = '&'.join(all_preds)

    return effect_predictions


def vep_dbnsfp_meta_vcf(query_vcf, info_tags_wanted):
    
    vcf = VCF(query_vcf)
    vep_csq_index2fields = {}
    vep_csq_fields2index = {}
    dbnsfp_prediction_algorithms = []
    for e in vcf.header_iter():
        header_element = e.info()
        if 'ID' in header_element.keys():
            identifier = str(header_element['ID'])
            if identifier == 'CSQ' or identifier == 'DBNSFP':
                description = str(header_element['Description'])
                if 'Format: ' in description:
                    subtags = description.split('Format: ')[1].split('|')
                    if identifier == 'CSQ':
                        i = 0
                        for t in subtags:
                            v = t.replace('"', '')                            
                            if v in info_tags_wanted:
                                vep_csq_index2fields[i] = v
                                vep_csq_fields2index[v] = i
                            i = i + 1
                    if identifier == 'DBNSFP':
                        i = 6
                        while (i < len(subtags)):
                            dbnsfp_prediction_algorithms.append(
                                str(re.sub(r'((_score)|(_pred))"*$', '', subtags[i])))
                            i = i + 1

    vep_dbnsfp_meta_info = {}
    vep_dbnsfp_meta_info['vep_csq_fieldmap'] = {}
    vep_dbnsfp_meta_info['vep_csq_fieldmap']['field2index'] = vep_csq_fields2index
    vep_dbnsfp_meta_info['vep_csq_fieldmap']['index2field'] = vep_csq_index2fields
    vep_dbnsfp_meta_info['dbnsfp_prediction_algorithms'] = dbnsfp_prediction_algorithms

    return vep_dbnsfp_meta_info

