#!/usr/bin/env python

import csv
import re
import argparse
import cyvcf2
import os
import annoutils
from pcgr import utils
from pcgr.utils import error_message, check_subprocess

csv.field_size_limit(500 * 1024 * 1024)

def __main__():
    parser = argparse.ArgumentParser(description='Cancer gene annotations from PCGR pipeline (SNVs/InDels)')
    parser.add_argument('vcf_file', help='VCF file with VEP-annotated query variants (SNVs/InDels)')
    parser.add_argument('pon_annotation',default=0,type=int,help='Include Panel of Normals annotation')
    parser.add_argument('regulatory_annotation',default=0,type=int,help='Inclusion of VEP regulatory annotations (0/1)')
    parser.add_argument('pcgr_db_dir',help='PCGR data directory')
    parser.add_argument('--cpsr',action="store_true",help="Aggregate cancer gene annotations for Cancer Predisposition Sequencing Reporter (CPSR)")
    parser.add_argument("--debug", action="store_true", default=False, help="Print full commands to log, default: %(default)s")
    args = parser.parse_args()

    logger = utils.getlogger('pcgr-gene-annotate')
    if args.cpsr is True:
        logger = utils.getlogger('cpsr-gene-annotate')

    extend_vcf_annotations(args.vcf_file, args.pcgr_db_dir, logger, args.pon_annotation, args.regulatory_annotation, args.cpsr, args.debug)

def extend_vcf_annotations(query_vcf, pcgr_db_dir, logger, pon_annotation, regulatory_annotation, cpsr, debug):
    """
    Function that reads VEP/vcfanno-annotated VCF and extends the VCF INFO column with tags from
    1. CSQ elements within the primary transcript consequence picked by VEP, e.g. SYMBOL, Feature, Gene, Consequence etc.
    2. Cancer-relevant gene annotations (PCGR_ONCO_XREF), e.g. known oncogenes/tumor suppressors, known antineoplastic drugs interacting with a given protein etc.
    3. Protein-relevant annotations, e.g. cancer hotspot mutations, functional protein features etc.
    4. Variant effect predictions
    5. Panel-of-normal (blacklisted variants) annotation

    List of INFO tags to be produced is provided by the 'infotags' files in the pcgr_db_dir
    """

    ## read VEP and PCGR tags to be appended to VCF file
    vcf_infotags_meta = annoutils.read_infotag_file(os.path.join(pcgr_db_dir, 'pcgr_infotags.tsv'))
    if cpsr is True:
        vcf_infotags_meta = annoutils.read_infotag_file(os.path.join(pcgr_db_dir, 'cpsr_infotags.tsv'))
    pcgr_onco_xref_map = annoutils.read_genexref_namemap(os.path.join(pcgr_db_dir, 'pcgr_onco_xref', 'pcgr_onco_xref_namemap.tsv'))


    out_vcf = re.sub(r'\.vcf(\.gz){0,}$','.annotated.vcf',query_vcf)

    meta_vep_dbnsfp_info = annoutils.vep_dbnsfp_meta_vcf(query_vcf, vcf_infotags_meta)
    dbnsfp_prediction_algorithms = meta_vep_dbnsfp_info['dbnsfp_prediction_algorithms']
    vep_csq_fields_map = meta_vep_dbnsfp_info['vep_csq_fieldmap']
    vcf = cyvcf2.VCF(query_vcf)
    for tag in sorted(vcf_infotags_meta):
        if pon_annotation == 0 and regulatory_annotation == 0:
            if not tag.startswith('PANEL_OF_NORMALS') and not tag.startswith('REGULATORY_'):
                vcf.add_info_to_header({'ID': tag, 'Description': str(vcf_infotags_meta[tag]['description']),'Type':str(vcf_infotags_meta[tag]['type']), 'Number': str(vcf_infotags_meta[tag]['number'])})
        elif pon_annotation == 1 and regulatory_annotation == 0:
            if not tag.startswith('REGULATORY_'):
                vcf.add_info_to_header({'ID': tag, 'Description': str(vcf_infotags_meta[tag]['description']),'Type':str(vcf_infotags_meta[tag]['type']), 'Number': str(vcf_infotags_meta[tag]['number'])})
        elif pon_annotation == 0 and regulatory_annotation == 1:
            if not tag.startswith('PANEL_OF_NORMALS'):
                vcf.add_info_to_header({'ID': tag, 'Description': str(vcf_infotags_meta[tag]['description']),'Type':str(vcf_infotags_meta[tag]['type']), 'Number': str(vcf_infotags_meta[tag]['number'])})
        else:
            vcf.add_info_to_header({'ID': tag, 'Description': str(vcf_infotags_meta[tag]['description']),'Type':str(vcf_infotags_meta[tag]['type']), 'Number': str(vcf_infotags_meta[tag]['number'])})

    w = cyvcf2.Writer(out_vcf, vcf)
    current_chrom = None
    num_chromosome_records_processed = 0

    vcf_info_element_types = {}
    for e in vcf.header_iter():
        header_element = e.info()
        if 'ID' in header_element and 'HeaderType' in header_element and 'Type' in header_element:
            identifier = str(header_element['ID'])
            fieldtype = str(header_element['Type'])
            vcf_info_element_types[identifier] = fieldtype

    vars_no_csq = list()
    for rec in vcf:
        if current_chrom is None:
            current_chrom = str(rec.CHROM)
            num_chromosome_records_processed = 0
        else:
            if str(rec.CHROM) != current_chrom:
                if not current_chrom is None:
                    logger.info(f"Completed summary of functional annotations for {num_chromosome_records_processed} variants on chr{current_chrom}")
                current_chrom = str(rec.CHROM)
                num_chromosome_records_processed = 0
        if rec.INFO.get('CSQ') is None:
            alt_allele = ','.join(rec.ALT)
            pos = rec.start + 1
            variant_id = f"g.{rec.CHROM}:{pos}{rec.REF}>{alt_allele}"
            vars_no_csq.append(variant_id)
            continue

        num_chromosome_records_processed += 1
        pcgr_onco_xref = annoutils.make_transcript_xref_map(rec, pcgr_onco_xref_map, xref_tag = "PCGR_ONCO_XREF")

        if regulatory_annotation == 1:
            csq_record_results_all = annoutils.parse_vep_csq(rec, pcgr_onco_xref, vep_csq_fields_map, logger, pick_only = False, csq_identifier = 'CSQ')
            if 'vep_block' in csq_record_results_all:
                vep_csq_records_all = csq_record_results_all['vep_block']
                rec.INFO['REGULATORY_ANNOTATION'] = annoutils.map_regulatory_variant_annotations(vep_csq_records_all)

        csq_record_results_pick = annoutils.parse_vep_csq(rec, pcgr_onco_xref, vep_csq_fields_map, logger, pick_only = True, csq_identifier = 'CSQ')
        vep_csq_records = None
        if 'vep_all_csq' in csq_record_results_pick:
            rec.INFO['VEP_ALL_CSQ'] = ','.join(csq_record_results_pick['vep_all_csq'])
        if 'vep_block' in csq_record_results_pick:
            vep_csq_records = csq_record_results_pick['vep_block']
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
    if vars_no_csq:
        logger.warning(f"There were {len(vars_no_csq)} records with no CSQ tag from VEP (was --vep_no_intergenic flag set?). Skipping them and showing (up to) the first 100:")
        print('----')
        print(', '.join(vars_no_csq[:100]))
        print('----')
    w.close()
    if current_chrom is not None:
        logger.info(f"Completed summary of functional annotations for {num_chromosome_records_processed} variants on chr{current_chrom}")
    vcf.close()

    if os.path.exists(out_vcf):
        if os.path.getsize(out_vcf) > 0:
            check_subprocess(logger, f'bgzip -f {out_vcf}', debug=False)
            check_subprocess(logger, f'tabix -f -p vcf {out_vcf}.gz', debug=False)
            annotated_vcf = f'{out_vcf}.gz'
            annoutils.write_pass_vcf(annotated_vcf, logger)
        else:
            error_message('No remaining PASS variants found in query VCF - exiting and skipping STEP 4', logger)
    else:
        error_message('No remaining PASS variants found in query VCF - exiting and skipping STEP 4', logger)

if __name__=="__main__":
    __main__()
