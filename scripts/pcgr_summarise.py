#!/usr/bin/env python

import csv
import re
import argparse
import cyvcf2
import os
import yaml
import glob

from glob import glob
from pcgr.annoutils import read_infotag_file, make_transcript_xref_map, read_genexref_namemap, map_regulatory_variant_annotations, write_pass_vcf, load_grantham
from pcgr.vep import parse_vep_csq
from pcgr.dbnsfp import vep_dbnsfp_meta_vcf, map_variant_effect_predictors
from pcgr.oncogenicity import load_oncogenicity_criteria, assign_oncogenicity_evidence, load_oncogenic_variants, match_oncogenic_variants
from pcgr.mutation_hotspot import load_mutation_hotspots, match_csq_mutation_hotspot
from pcgr.biomarker import load_biomarkers, match_csq_biomarker
from pcgr.utils import error_message, check_subprocess, getlogger, remove_file
from pcgr.splice import load_splice_effects

csv.field_size_limit(500 * 1024 * 1024)

def __main__():
    parser = argparse.ArgumentParser(description='Summarise VEP annotations (gene/variant) from PCGR/CPSR pipeline (SNVs/InDels)')
    parser.add_argument('vcf_file_in', help='Bgzipped VCF file with VEP-annotated query variants (SNVs/InDels)')
    parser.add_argument('vcf_file_out', help='Bgzipped VCF file with extended VEP-annotated query variants (SNVs/InDels)')
    parser.add_argument('pon_annotation',default=0,type=int,help='Include Panel of Normals annotation (0/1)')
    parser.add_argument('regulatory_annotation',default=0,type=int,help='Inclusion of VEP regulatory annotations (0/1)')
    parser.add_argument('oncogenicity_annotation',default=0,type=int,help='Include oncogenicity annotation (0/1)')
    parser.add_argument('tumortype', default='Any', help='Primary tumor type of query VCF')
    parser.add_argument('build', default='GRCh38', help='Genome build of query VCF')
    parser.add_argument('vep_pick_order', default="mane,canonical,appris,biotype,ccds,rank,tsl,length", 
                        help="Comma-separated string of ordered transcript/variant properties for selection of primary variant consequence")
    parser.add_argument('refdata_assembly_dir',help='Assembly-specific reference data directory, e.g. "pcgrdb/data/grch38')
    parser.add_argument('--compress_output_vcf', action="store_true", default=False, help="Compress output VCF file")
    parser.add_argument('--cpsr',action="store_true",help="Aggregate cancer gene annotations for Cancer Predisposition Sequencing Reporter (CPSR)")
    parser.add_argument('--cpsr_yaml',dest="cpsr_yaml", default=None, help='YAML file with list of targeted genes by CPSR (custom, or panel-defined)')
    parser.add_argument("--debug", action="store_true", default=False, help="Print full commands to log, default: %(default)s")
    args = parser.parse_args()

    logger = getlogger('pcgr-gene-annotate')
    if args.cpsr is True:
        logger = getlogger('cpsr-gene-annotate')
    
    arg_dict = vars(args)
    
    extend_vcf_annotations(arg_dict, logger)

def extend_vcf_annotations(arg_dict, logger):
    """
    Function that reads VEP/vcfanno-annotated VCF and populates the VCF INFO column with tags from
    0. Information on sequencing depth and genotype (hom_ref, het, hom_alt) - for germline VCFs (CPSR) only 
    1. CSQ elements for the primary (i.e. "picked") gene transcript consequence from VEP, e.g. SYMBOL, Feature, Gene, Consequence etc.
    2. Cancer-relevant gene annotations (GENE_TRANSCRIPT_XREF), e.g. known oncogenes/tumor suppressors, driver genes etc
    3. Variant effect predictions - dbNSFP
    4. Panel-of-normal (blacklisted variants) annotation

    Moreover, it performs two important matching procedures, using 
    5. Information from VEP's CSQ information (HGVSp/HGVSc) to match known mutation hotspots in cancer
    6. Information from VEP's CSQ information (Consequence, HGVSp/HGVSc, exon information etc) and genomic coordinates to match known biomarkers in cancer
    
    Finally, it assesses somatic variant oncogenicity, using
    7. Gene annotations (tumor suppressor, oncogene) and variant annotations (loss-of-function, gnomAD variant frequencies, variant effect predictions).
       Variant oncogenicity levels are provided for all variants using a recommended five-level scheme ("Oncogenic", "Likely oncogenic", "VUS", "Likely Benign", "Benign")
       - Recommended scoring scheme for variant oncogenicity classification outlined by VICC/ClinGen consortia (Horak et al., Genet Med, 2022)

    List of VCF INFO tags appended by this procedure is defined by the 'infotags' files in the refdata_assembly_dir
    """
    
    vcf_infotags = {}
    
    vcf_infotags['other'] = read_infotag_file(os.path.join(arg_dict['refdata_assembly_dir'], 'vcf_infotags_other.tsv'), scope = "pcgr")
    if arg_dict['cpsr'] is True:
        vcf_infotags['other'] = read_infotag_file(os.path.join(arg_dict['refdata_assembly_dir'], 'vcf_infotags_other.tsv'), scope = "cpsr")
    vcf_infotags['vep'] = read_infotag_file(os.path.join(arg_dict['refdata_assembly_dir'], 'vcf_infotags_vep.tsv'), scope = "vep")
    vcf_infotags['other'].update(vcf_infotags['vep'])
    vcf_info_metadata = vcf_infotags['other']
    
    
    ## load CPSR target genes from YAML file
    cpsr_target_genes = {}
    if arg_dict['cpsr'] is True:
        if arg_dict['cpsr_yaml'] is not None:
            if os.path.exists(arg_dict['cpsr_yaml']):
                with open(arg_dict['cpsr_yaml'], 'r') as f:
                    yaml_data = yaml.safe_load(f)
                    if yaml_data is not None:
                        panel_genes = yaml_data['conf']['gene_panel']['panel_genes']
                        for g in panel_genes:
                            if 'entrezgene' in g:                                
                                cpsr_target_genes[str(g['entrezgene'])] = 1
                            
            
        vcf_infotags['cpsr'] = read_infotag_file(os.path.join(arg_dict['refdata_assembly_dir'], 'vcf_infotags_cpsr.tsv'), scope = "cpsr")
        vcf_infotags['other'].update(vcf_infotags['cpsr'])
        vcf_info_metadata.update(vcf_infotags['cpsr'])

    gene_transcript_xref_map = read_genexref_namemap(
        os.path.join(arg_dict['refdata_assembly_dir'], 'gene','tsv','gene_transcript_xref', 'gene_transcript_xref_bedmap.tsv.gz'), logger)
    cancer_hotspots = load_mutation_hotspots(
        os.path.join(arg_dict['refdata_assembly_dir'], 'misc','tsv','hotspot', 'hotspot.tsv.gz'), logger)
    oncogenicity_criteria = load_oncogenicity_criteria(
        os.path.join(arg_dict['refdata_assembly_dir'], 'misc','tsv', 'oncogenicity', 'oncogenicity.tsv'), logger)
    oncogenic_variants = load_oncogenic_variants(os.path.join(arg_dict['refdata_assembly_dir'], 'variant','tsv', 'clinvar','clinvar_oncogenic.tsv.gz'), logger)
    grantham_scores = load_grantham(os.path.join(arg_dict['refdata_assembly_dir'], 'misc','tsv', 'grantham', 'grantham.tsv'), logger)
    splice_effects = load_splice_effects(os.path.join(arg_dict['refdata_assembly_dir'], 'misc','tsv', 'mutsplicedb', 'mutsplicedb.tsv'), logger)

    biomarkers = {}
    for db in ['cgi','civic']:
        variant_fname = os.path.join(arg_dict['refdata_assembly_dir'], 'biomarker','tsv', f"{db}.variant.tsv.gz")
        clinical_fname = os.path.join(arg_dict['refdata_assembly_dir'], 'biomarker','tsv', f"{db}.clinical.tsv.gz")
        if arg_dict['cpsr'] is True:
            biomarkers[db] = load_biomarkers(
                logger, variant_fname, clinical_fname, biomarker_vartype = 'MUT', biomarker_variant_origin = 'Both')
        else:
            biomarkers[db] = load_biomarkers(
                logger, variant_fname, clinical_fname, biomarker_vartype = 'MUT', biomarker_variant_origin = 'Somatic')

    out_vcf = re.sub(r'(\.gz)$','',arg_dict['vcf_file_out'])

    meta_vep_dbnsfp_info = vep_dbnsfp_meta_vcf(arg_dict['vcf_file_in'], vcf_info_metadata)
    dbnsfp_prediction_algorithms = meta_vep_dbnsfp_info['dbnsfp_prediction_algorithms']
    vep_csq_fields_map = meta_vep_dbnsfp_info['vep_csq_fieldmap']
    
    vcf = cyvcf2.VCF(arg_dict['vcf_file_in'],  gts012=True)
    for tag in sorted(vcf_info_metadata):
        if arg_dict['pon_annotation'] == 0 and tag.startswith('PANEL_OF_NORMALS'):
            continue
        if arg_dict['regulatory_annotation'] == 0 and tag.startswith('REGULATORY_'):
            continue
        vcf.add_info_to_header({'ID': tag, 'Description': str(vcf_info_metadata[tag]['description']),
                                'Type':str(vcf_info_metadata[tag]['type']), 
                                'Number': str(vcf_info_metadata[tag]['number'])})
    w = cyvcf2.Writer(arg_dict['vcf_file_out'], vcf)
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
        alt_allele = ','.join(rec.ALT)
        pos = rec.start + 1
        variant_id = f"g.{rec.CHROM}:{pos}{rec.REF}>{alt_allele}"
        
        ## For germline-called VCFs (single-sample), pull out sequencing depth and genotype using cyvcf2 API
        if arg_dict['cpsr'] is True:
            rec.INFO['DP_CONTROL'] = "-1"
            if rec.gt_depths is not None:
                if len(rec.gt_depths) == 1:
                    rec.INFO['DP_CONTROL'] = str(rec.gt_depths[0])
            rec.INFO['GENOTYPE'] = 'undefined'
            if rec.gt_types is not None:
                if len(rec.gt_types) == 1:
                    if int(rec.gt_types[0]) == 0:
                        rec.INFO['GENOTYPE'] = 'hom_ref'
                    else:
                        if int(rec.gt_types[0]) == 1:
                            rec.INFO['GENOTYPE'] = 'het'
                        if int(rec.gt_types[0]) == 2:
                            rec.INFO['GENOTYPE'] = 'hom_alt'
                
        if current_chrom is None:
            current_chrom = str(rec.CHROM)
            num_chromosome_records_processed = 0
        else:
            if str(rec.CHROM) != current_chrom:
                if current_chrom is not None:
                    logger.info(f"Completed summary of functional annotations for {num_chromosome_records_processed} variants on chr{current_chrom}")
                current_chrom = str(rec.CHROM)
                num_chromosome_records_processed = 0
        if rec.INFO.get('CSQ') is None:
            vars_no_csq.append(variant_id)
            continue

        num_chromosome_records_processed += 1
        transcript_xref_map = make_transcript_xref_map(rec, gene_transcript_xref_map, xref_tag = "GENE_TRANSCRIPT_XREF")

        vep_csq_record_results = {}
        if arg_dict['cpsr'] is True:
            vep_csq_record_results = \
                parse_vep_csq(rec, transcript_xref_map, 
                              vep_csq_fields_map, 
                              grantham_scores,
                              arg_dict['vep_pick_order'], 
                              logger,  
                              pick_only = False, 
                              csq_identifier = 'CSQ', 
                              debug = arg_dict['debug'],
                              targets_entrez_gene = cpsr_target_genes)
            if 'picked_gene_csq' in vep_csq_record_results:
                rec.INFO['REGULATORY_ANNOTATION'] = map_regulatory_variant_annotations(
                    vep_csq_record_results['picked_gene_csq'])
        else:
            vep_csq_record_results = \
                parse_vep_csq(rec, 
                              transcript_xref_map, 
                              vep_csq_fields_map, 
                              grantham_scores, 
                              arg_dict['vep_pick_order'], 
                              logger, 
                              pick_only = True, 
                              csq_identifier = 'CSQ', 
                              debug = arg_dict['debug'])

        principal_csq_properties = {}
        for prop in ['hgvsp','hgvsc','entrezgene','exon','codon','lof','refseq_transcript_id']:
            principal_csq_properties[prop] = '.'
        
        if 'picked_csq' in vep_csq_record_results:
            csq_record = vep_csq_record_results['picked_csq']
            for k in csq_record:
                if k in vcf_info_element_types:
                    if vcf_info_element_types[k] == "Flag" and csq_record[k] == "1":
                        rec.INFO[k] = True
                    else:
                        if csq_record[k] is not None:                            
                            if k == 'DOMAINS':
                                filtered_domains = []
                                for dom in csq_record[k].split('&'):
                                    if not dom.startswith('PDB-'):
                                        filtered_domains.append(dom)
                                if(len(filtered_domains) > 0):
                                    csq_record[k] = '&'.join(filtered_domains)
                            
                            rec.INFO[k] = csq_record[k]

                            if k == 'HGVSp_short':
                                principal_csq_properties['hgvsp'] = csq_record[k]
                                if re.match(r'^(p.[A-Z]{1}[0-9]{1,}[A-Za-z]{1,})', principal_csq_properties['hgvsp']):
                                    codon_match = re.findall(r'[A-Z][0-9]{1,}', principal_csq_properties['hgvsp'])
                                    if len(codon_match) == 1:
                                        principal_csq_properties['codon'] = 'p.' + codon_match[0]
                     
                            if k == 'HGVSc':
                                principal_csq_properties['hgvsc'] = csq_record[k].split(':')[1]
                            
                            if k == 'REFSEQ_TRANSCRIPT_ID':
                                principal_csq_properties['refseq_transcript_id'] = csq_record[k]
                            
                            if k == 'ENTREZGENE':
                                principal_csq_properties['entrezgene'] = csq_record[k]
                            
                            if k == 'LOSS_OF_FUNCTION':
                                principal_csq_properties['lof'] = csq_record[k]
                            
                            if k == 'EXON':
                                if "/" in csq_record[k]:
                                    principal_csq_properties['exon'] = csq_record[k].split('/')[0]
                            
        
        splice_key = f"{principal_csq_properties['entrezgene']}_{principal_csq_properties['refseq_transcript_id']}_{principal_csq_properties['hgvsc']}"
        if splice_key in splice_effects:
            if 'NO_EFFECT' not in splice_effects[splice_key]:
                rec.INFO['SPLICE_EFFECT_MUTSPLICEDB'] = splice_effects[splice_key]
                rec.INFO['LOSS_OF_FUNCTION'] = True
                if rec.INFO.get('LOF_FILTER') is not None:
                    rec.INFO['LOF_FILTER'] = '.'
        
        if 'all_csq' in vep_csq_record_results:
            rec.INFO['VEP_ALL_CSQ'] = ','.join(vep_csq_record_results['all_csq'])
            match_csq_mutation_hotspot(vep_csq_record_results['all_csq'], cancer_hotspots, rec, principal_csq_properties)

            for db in ['civic','cgi']:                
                match_csq_biomarker(vep_csq_record_results['all_csq'], biomarkers[db], rec, principal_csq_properties)

        if rec.INFO.get('DBNSFP') is not None:
            map_variant_effect_predictors(rec, dbnsfp_prediction_algorithms)
        
        if arg_dict['oncogenicity_annotation'] == 1:
            match_oncogenic_variants(vep_csq_record_results['all_csq'], oncogenic_variants, rec, principal_csq_properties)
            assign_oncogenicity_evidence(rec, oncogenicity_criteria, tumortype = arg_dict['tumortype'])

        if "GENE_TRANSCRIPT_XREF" in vcf_info_element_types:
            gene_xref_tag = rec.INFO.get('GENE_TRANSCRIPT_XREF')
            if gene_xref_tag is not None:
                del rec.INFO['GENE_TRANSCRIPT_XREF']                
        w.write_record(rec)
    if vars_no_csq:
        logger.warning(f"There were {len(vars_no_csq)} records with no CSQ tag from VEP (was --vep_no_intergenic flag set?). Skipping them and showing (up to) the first 100:")
        logger.info('----')
        logger.info(', '.join(vars_no_csq[:100]))
        logger.info('----')
    w.close()
    if current_chrom is not None:
        logger.info(f"Completed summary of functional annotations for {num_chromosome_records_processed} variants on chr{current_chrom}")
    vcf.close()

    ## get list of all intermediate vcf files
    vcf_prefix = re.sub(r'(\.vep\.vcfanno\.vcf\.gz)$','',arg_dict['vcf_file_in'])
    vcf_prefix_all_files = glob(f'{vcf_prefix}*')

    if os.path.exists(arg_dict['vcf_file_out']):
        if os.path.getsize(arg_dict['vcf_file_out']) > 0:
            if arg_dict['compress_output_vcf'] is True:
                check_subprocess(logger, f'bgzip -f {out_vcf}', debug=arg_dict['debug'])
                check_subprocess(logger, f'tabix -f -p vcf {out_vcf}.gz', debug=arg_dict['debug'])
                summarized_vcf = f'{arg_dict["vcf_file_out"]}.gz'
                vcf_no_pass_variants = write_pass_vcf(summarized_vcf, logger)
                vcf_prefix_all_files = glob(f'{vcf_prefix}*')               
                if bool(vcf_no_pass_variants):
                    logger.info('Cleaning up all intermediate files from STEP 1/2/3')
                    if not arg_dict['debug'] and len(vcf_prefix_all_files) == 13:
                        for f in vcf_prefix_all_files:
                            remove_file(f)
                    error_message('No remaining PASS variants found in query VCF - exiting STEP 3', logger)
        else:
            if not arg_dict['debug'] and len(vcf_prefix_all_files) == 10:
                logger.info('Cleaning up all intermediate files from STEP 1/2/3')
                for f in vcf_prefix_all_files:
                    remove_file(f)
            error_message('No remaining PASS variants found in query VCF - exiting STEP 3', logger)
    else:
        if not arg_dict['debug'] and len(vcf_prefix_all_files) == 10:
            logger.info('Cleaning up all intermediate files from STEP 1/2/3')
            for f in vcf_prefix_all_files:
                remove_file(f)
        error_message('No remaining PASS variants found in query VCF - exiting STEP 3', logger)

if __name__=="__main__":
    __main__()
