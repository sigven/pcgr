#!/usr/bin/env python

import csv
import re
import argparse
import cyvcf2
import os

from pcgr import annoutils
from pcgr.annoutils import read_infotag_file, make_transcript_xref_map, read_genexref_namemap, map_regulatory_variant_annotations, parse_vep_csq,  write_pass_vcf
from pcgr import dbnsfp 
from pcgr.dbnsfp import vep_dbnsfp_meta_vcf, map_variant_effect_predictors
from pcgr import oncogenicity
from pcgr.oncogenicity import assign_oncogenicity_evidence
from pcgr import mutation_hotspot
from pcgr.mutation_hotspot import load_mutation_hotspots, match_csq_mutation_hotspot
from pcgr import biomarker
from pcgr.biomarker import load_biomarkers, match_csq_biomarker
from pcgr import utils
from pcgr.utils import error_message, check_subprocess

csv.field_size_limit(500 * 1024 * 1024)

def __main__():
    parser = argparse.ArgumentParser(description='Cancer gene annotations from PCGR pipeline (SNVs/InDels)')
    parser.add_argument('vcf_file', help='VCF file with VEP-annotated query variants (SNVs/InDels)')
    parser.add_argument('pon_annotation',default=0,type=int,help='Include Panel of Normals annotation')
    parser.add_argument('regulatory_annotation',default=0,type=int,help='Inclusion of VEP regulatory annotations (0/1)')
    parser.add_argument('oncogenicity_annotation',default=0,type=int,help='Include oncogenicity annotation (0/1)')
    parser.add_argument('pcgr_db_dir',help='PCGR data directory')
    parser.add_argument('--cpsr',action="store_true",help="Aggregate cancer gene annotations for Cancer Predisposition Sequencing Reporter (CPSR)")
    parser.add_argument("--debug", action="store_true", default=False, help="Print full commands to log, default: %(default)s")
    args = parser.parse_args()

    logger = utils.getlogger('pcgr-gene-annotate')
    if args.cpsr is True:
        logger = utils.getlogger('cpsr-gene-annotate')

    extend_vcf_annotations(args.vcf_file, args.pcgr_db_dir, logger, args.pon_annotation, args.regulatory_annotation, args.oncogenicity_annotation, args.cpsr, args.debug)

def extend_vcf_annotations(query_vcf, pcgr_db_dir, logger, pon_annotation, regulatory_annotation, oncogenicity_annotation, cpsr, debug):
    """
    Function that reads VEP/vcfanno-annotated VCF and extends the VCF INFO column with tags from
    1. CSQ elements within the primary transcript consequence picked by VEP, e.g. SYMBOL, Feature, Gene, Consequence etc.
    2. Cancer-relevant gene annotations (GENE_TRANSCRIPT_XREF), e.g. known oncogenes/tumor suppressors, driver genes etc
    3. Variant effect predictions - dbNSFP
    4. Panel-of-normal (blacklisted variants) annotation

    Moreover, it uses 
    1. Information from VEP's CSQ information (HGVSp/HGVSc) to match known mutation hotspots in cancer
    2. Gene annotations (tumor suppressor, oncogene) and variant annotations (loss-of-function, gnomAD variant frequencies, variant effect predictions) to assign
       variant oncogenicity levels ("Oncogenic", "Likely oncogenic", "VUS", "Likely Benign", "Benign")

    List of VCF INFO tags to be appended is defined by the 'infotags' files in the pcgr_db_dir
    """
    vcf_infotags_pcgr = read_infotag_file(os.path.join(pcgr_db_dir, 'pcgr_vcf_infotags.tsv'))
    vcf_infotags_vep = read_infotag_file(os.path.join(pcgr_db_dir, 'vep_vcf_infotags.tsv'))
    vcf_infotags_meta = {}

    for tag in vcf_infotags_pcgr:
        vcf_infotags_meta[tag] = vcf_infotags_pcgr[tag]
    for tag in vcf_infotags_vep:
        vcf_infotags_meta[tag] = vcf_infotags_vep[tag]

    if cpsr is True:
        vcf_infotags_meta = read_infotag_file(os.path.join(pcgr_db_dir, 'cpsr_vcf_infotags.tsv'))
        ## add gnomad non-cancer
    gene_transcript_xref_map = read_genexref_namemap(os.path.join(pcgr_db_dir, 'gene','tsv','gene_transcript_xref', 'gene_transcript_xref_bedmap.tsv.gz'), logger)
    cancer_hotspots = load_mutation_hotspots(os.path.join(pcgr_db_dir, 'misc','tsv','hotspot', 'hotspot.tsv.gz'), logger)

    biomarkers = {}
    for db in ['cgi','civic']:
        variant_fname = os.path.join(pcgr_db_dir, 'biomarker','tsv', f"{db}.variant.tsv.gz")
        clinical_fname = os.path.join(pcgr_db_dir, 'biomarker','tsv', f"{db}.clinical.tsv.gz")
        biomarkers[db] = load_biomarkers(logger, variant_fname, clinical_fname)

    out_vcf = re.sub(r'\.vcf(\.gz){0,}$','.annotated.vcf',query_vcf)

    meta_vep_dbnsfp_info = vep_dbnsfp_meta_vcf(query_vcf, vcf_infotags_meta)
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
        pcgr_onco_xref = make_transcript_xref_map(rec, gene_transcript_xref_map, xref_tag = "GENE_TRANSCRIPT_XREF")

        if regulatory_annotation == 1:
            csq_record_results_all = parse_vep_csq(rec, pcgr_onco_xref, vep_csq_fields_map, logger, pick_only = False, csq_identifier = 'CSQ')
            if 'picked_gene_csq' in csq_record_results_all:
                vep_csq_records_all = csq_record_results_all['picked_gene_csq']
                rec.INFO['REGULATORY_ANNOTATION'] = map_regulatory_variant_annotations(vep_csq_records_all)

        vep_csq_record_results = parse_vep_csq(rec, pcgr_onco_xref, vep_csq_fields_map, logger, pick_only = True, csq_identifier = 'CSQ')

        vep_csq_records = None

        principal_csq_properties = {}

        principal_csq_properties['hgvsp'] = '.'
        principal_csq_properties['hgvsc'] = '.'
        principal_csq_properties['entrezgene'] = '.'
        principal_csq_properties['exon'] = '.'
        principal_csq_properties['codon'] = '.'

        if 'picked_csq' in vep_csq_record_results:
            csq_record = vep_csq_record_results['picked_csq']
            for k in csq_record:
                if k in vcf_info_element_types:
                    if vcf_info_element_types[k] == "Flag" and csq_record[k] == "1":
                        rec.INFO[k] = True
                    else:
                        if not csq_record[k] is None:
                            rec.INFO[k] = csq_record[k]

                            #print(k)
                            if k == 'HGVSp_short':
                                principal_csq_properties['hgvsp'] = csq_record[k]
                                if re.match(r'^(p.[A-Z]{1}[0-9]{1,}[A-Za-z]{1,})', principal_csq_properties['hgvsp']):
                                    codon_match = re.findall(r'[A-Z][0-9]{1,}', principal_csq_properties['hgvsp'])
                                    if len(codon_match) == 1:
                                        principal_csq_properties['codon'] = codon_match[0]
                     
                            if k == 'HGVSc':
                                principal_csq_properties['hgvsc'] = csq_record[k].split(':')[1]
                            
                            if k == 'ENTREZGENE':
                                principal_csq_properties['entrezgene'] = csq_record[k]
                            
                            if k == 'EXON':
                                if "/" in csq_record[k]:
                                    principal_csq_properties['exon'] = csq_record[k].split('/')[0]
        
        if 'all_csq' in vep_csq_record_results:
            rec.INFO['VEP_ALL_CSQ'] = ','.join(vep_csq_record_results['all_csq'])
            match_csq_mutation_hotspot(vep_csq_record_results['all_csq'], cancer_hotspots, rec, principal_csq_properties)

            for db in ['civic','cgi']:
                match_csq_biomarker(vep_csq_record_results['all_csq'], biomarkers[db], rec, principal_csq_properties)

        if not rec.INFO.get('DBNSFP') is None:
            map_variant_effect_predictors(rec, dbnsfp_prediction_algorithms)
        
        if oncogenicity_annotation == 1:
            assign_oncogenicity_evidence(rec, tumortype = "Any")

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
            write_pass_vcf(annotated_vcf, logger)
        else:
            error_message('No remaining PASS variants found in query VCF - exiting and skipping STEP 4', logger)
    else:
        error_message('No remaining PASS variants found in query VCF - exiting and skipping STEP 4', logger)

if __name__=="__main__":
    __main__()
