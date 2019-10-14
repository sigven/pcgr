#!/usr/bin/env python

import csv
import re
import argparse
from cyvcf2 import VCF, Writer
import gzip
import os
import annoutils

csv.field_size_limit(500 * 1024 * 1024)

def __main__():
   
   parser = argparse.ArgumentParser(description='Cancer gene annotations from PCGR pipeline (SNVs/InDels)')
   parser.add_argument('vcf_file', help='VCF file with VEP-annotated query variants (SNVs/InDels)')
   parser.add_argument('pon_annotation',default=0,type=int,help='Include Panel of Normals annotation')
   parser.add_argument('pcgr_db_dir',help='PCGR data directory')
   parser.add_argument('--cpsr',action="store_true",help="Aggregate cancer gene annotations for Cancer Predisposition Sequencing Reporter (CPSR)")
   args = parser.parse_args()

   logger = annoutils.getlogger('pcgr-gene-annotate')
   if args.cpsr is True:
      logger = annoutils.getlogger('cpsr-gene-annotate')

   extend_vcf_annotations(args.vcf_file, args.pcgr_db_dir, logger, args.pon_annotation, args.cpsr)


def extend_vcf_annotations(query_vcf, pcgr_db_directory, logger, pon_annotation, cpsr):
   """
   Function that reads VEP/vcfanno-annotated VCF and extends the VCF INFO column with tags from
   1. CSQ elements within the primary transcript consequence picked by VEP, e.g. SYMBOL, Feature, Gene, Consequence etc.
   2. Cancer-relevant gene annotations, e.g. known oncogenes/tumor suppressors, known antineoplastic drugs interacting with a given protein etc.
   3. Protein-relevant annotations, e.g. cancer hotspot mutations, functional protein features etc.
   4. Variant effect predictions
   5. Panel-of-normal (blacklisted variants) annotation
   """

   ## read VEP and PCGR tags to be appended to VCF file
   vcf_infotags_meta = annoutils.read_infotag_file(os.path.join(pcgr_db_directory,'pcgr_infotags.tsv'))
   if cpsr is True:
      vcf_infotags_meta = annoutils.read_infotag_file(os.path.join(pcgr_db_directory,'cpsr_infotags.tsv'))

   out_vcf = re.sub(r'\.vcf(\.gz){0,}$','.annotated.vcf',query_vcf)

   meta_vep_dbnsfp_info = annoutils.vep_dbnsfp_meta_vcf(query_vcf, vcf_infotags_meta)
   dbnsfp_prediction_algorithms = meta_vep_dbnsfp_info['dbnsfp_prediction_algorithms']
   vep_csq_fields_map = meta_vep_dbnsfp_info['vep_csq_fieldmap']
   vcf = VCF(query_vcf)
   for tag in sorted(vcf_infotags_meta):
      if pon_annotation == 0:
         if not tag.startswith('PANEL_OF_NORMALS'):
            vcf.add_info_to_header({'ID': tag, 'Description': str(vcf_infotags_meta[tag]['description']),'Type':str(vcf_infotags_meta[tag]['type']), 'Number': str(vcf_infotags_meta[tag]['number'])})
      else:
         vcf.add_info_to_header({'ID': tag, 'Description': str(vcf_infotags_meta[tag]['description']),'Type':str(vcf_infotags_meta[tag]['type']), 'Number': str(vcf_infotags_meta[tag]['number'])})

   w = Writer(out_vcf, vcf)
   current_chrom = None
   num_chromosome_records_processed = 0
   pcgr_onco_xref_map = {'ENSEMBL_TRANSCRIPT_ID': 0, 'ENSEMBL_GENE_ID':1, 'ENSEMBL_PROTEIN_ID':2, 'SYMBOL':3, 'SYMBOL_ENTREZ':4, 'ENTREZ_ID':5, 'UNIPROT_ID':6, 'APPRIS':7,'UNIPROT_ACC':8,'REFSEQ_MRNA':9,'CORUM_ID':10,'TUMOR_SUPPRESSOR':11,
                        'ONCOGENE':12,'DRIVER':13, 'NETWORK_CG':14,'DISGENET_CUI':15,'CHEMBL_COMPOUND_ID':16,'INTOGEN_DRIVER':17,'TCGA_DRIVER':18,'ONCOSCORE':19, 'MIM_PHENOTYPE_ID':20, 'CANCER_PREDISPOSITION_SOURCE':21, 
                        'CANCER_SUSCEPTIBILITY_CUI':22, 'CANCER_SYNDROME_CUI':23, 'CANCER_PREDISPOSITION_MOI':24, 'CANCER_PREDISPOSITION_MOD':25, 'SIGNALING_PATHWAY':26, 'OPENTARGETS_DISEASE_ASSOCS':27,
                        'OPENTARGETS_TRACTABILITY_COMPOUND:':28, 'OPENTARGETS_TRACTABILITY_ANTIBODY':29, 'GE_PANEL_ID':30, 'ACTIONABLE_TARGET':31}
   
   vcf_info_element_types = {}
   for e in vcf.header_iter():
      header_element = e.info()
      if 'ID' in header_element and 'HeaderType' in header_element and 'Type' in header_element:
         identifier = str(header_element['ID'])
         fieldtype = str(header_element['Type'])
         vcf_info_element_types[identifier] = fieldtype

   for rec in vcf:
      if current_chrom is None:
         current_chrom = str(rec.CHROM)
         num_chromosome_records_processed = 0
      else:
         if str(rec.CHROM) != current_chrom:
            if not current_chrom is None:
               logger.info('Completed summary of functional annotations for ' + str(num_chromosome_records_processed) + ' variants on chromosome ' + str(current_chrom))
            current_chrom = str(rec.CHROM)
            num_chromosome_records_processed = 0
      if rec.INFO.get('CSQ') is None:
         alt_allele = ','.join(rec.ALT)
         pos = rec.start + 1
         variant_id = 'g.' + str(rec.CHROM) + ':' + str(pos) + str(rec.REF) + '>' + alt_allele
         logger.warning('Variant record ' + str(variant_id) + ' does not have CSQ tag from Variant Effect Predictor (vep_skip_intergenic in config set to true?)  - variant will be skipped')
         continue
      csq_record_results = {}
      num_chromosome_records_processed += 1
      pcgr_onco_xref = annoutils.make_transcript_xref_map(rec, pcgr_onco_xref_map, xref_tag = "PCGR_ONCO_XREF")
      csq_record_results = annoutils.parse_vep_csq(rec, pcgr_onco_xref, vep_csq_fields_map, logger, pick_only = True, csq_identifier = 'CSQ')

      vep_csq_records = None 
      if 'vep_all_csq' in csq_record_results:
         rec.INFO['VEP_ALL_CSQ'] = ','.join(csq_record_results['vep_all_csq'])
      if 'vep_block' in csq_record_results:
         vep_csq_records = csq_record_results['vep_block']

         block_idx = 0
         if cpsr is True:
            block_idx = annoutils.get_correct_cpg_transcript(vep_csq_records)
         record = vep_csq_records[block_idx]
         for k in record:
            if k in vcf_info_element_types:
               if vcf_info_element_types[k] == "Flag" and record[k] == "1":
                  rec.INFO[k] = True
               else:
                  if not record[k] is None:
                     rec.INFO[k] = record[k]
      
      if not rec.INFO.get('DBNSFP') is None:
         annoutils.map_variant_effect_predictors(rec, dbnsfp_prediction_algorithms)


      w.write_record(rec)
   w.close()
   if not current_chrom is None:
      logger.info('Completed summary of functional annotations for ' + str(num_chromosome_records_processed) + ' variants on chromosome ' + str(current_chrom))
   vcf.close()

   if os.path.exists(out_vcf):
      if os.path.getsize(out_vcf) > 0:
         os.system('bgzip -f ' + str(out_vcf))
         os.system('tabix -f -p vcf ' + str(out_vcf) + '.gz')
         annotated_vcf = out_vcf + '.gz'
         annoutils.write_pass_vcf(annotated_vcf, logger)
      else:
         annoutils.error_message('No remaining PASS variants found in query VCF - exiting and skipping STEP 4 (pcgr-writer)', logger)
   else:
      annoutils.error_message('No remaining PASS variants found in query VCF - exiting and skipping STEP 4 (pcgr-writer)', logger)

if __name__=="__main__": __main__()


      
