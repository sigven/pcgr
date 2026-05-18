#!/usr/bin/env python

import argparse
import os
import sys

from cyvcf2 import VCF

from pcgr import vcf
from pcgr.utils import check_subprocess, remove_file, random_id_generator, getlogger
from pcgr.vcf import check_existing_vcf_info_tags, check_retained_vcf_info_tags
from pcgr.validate import (
    bgzip_sort_filter_vcf, detect_multiallelic_sites, collect_workflow_info_tags,
    is_valid_cna, is_valid_rna_fusion, is_valid_rna_expression,
    is_valid_germline, is_valid_pon_vcf
)


def __main__():

    parser = argparse.ArgumentParser(description='Verify input data for PCGR')
    parser.add_argument('refdata_assembly_dir',help='Assembly-specific directory with reference data for PCGR/CPSR')
    parser.add_argument('input_vcf', help='Somatic (tumor) query variants (SNVs/InDels) - VCF format')
    parser.add_argument('validated_vcf', help="Validated VCF file with somatic (tumor) query variants (SNVs/InDels)")
    parser.add_argument('input_cna', help='Somatic (tumor) copy number query segments (tab-separated values)')
    parser.add_argument('input_rna_fusion', help='Tumor RNA fusion variants (tab-separated values)')
    parser.add_argument('input_rna_exp', help='Tumor gene expression estimates (tab-separated values)')
    parser.add_argument('input_cpsr', help='Classified germline calls from CPSR (tab-separated values)')
    parser.add_argument('panel_normal_vcf',help="VCF file with germline calls from panel of normals")
    parser.add_argument('tumor_only',type=int, default=0,choices=[0,1],help="Tumor only sequencing")
    parser.add_argument('sample_id',help='PCGR sample_name')
    parser.add_argument('build',help='Genome build (grch37/grch38)')
    parser.add_argument('retained_info_tags', help="Comma-separated string of custom VCF INFO tags to be kept in PCGR output")
    parser.add_argument('tumor_dp_tag', help='VCF INFO tag that denotes tumor sequencing depth')
    parser.add_argument('tumor_af_tag', help='VCF INFO tag that denotes tumor variant allelic fraction')
    parser.add_argument('control_dp_tag', help='VCF INFO tag that denotes control sequencing depth')
    parser.add_argument('control_af_tag', help='VCF INFO tag that denotes control variant allelic fraction')
    parser.add_argument('call_conf_tag', help='VCF INFO tag that denotes somatic variant call confidence')
    parser.add_argument('exclude_hom_germline', help='Logical indicating if homozygote germline calls are to be filtered based on allelic fraction')
    parser.add_argument('exclude_het_germline', help='Logical indicating if heterozygote germline calls are to be filtered based on allelic fraction')
    parser.add_argument('--keep_uncompressed', action="store_true", help='Keep uncompressed VCF for vcf2maf.pl')
    parser.add_argument('--output_dir', dest='output_dir', help='Output directory')
    parser.add_argument("--debug", action="store_true", help="Print full commands to log")
    args = parser.parse_args()

    ret = validate_pcgr_input(args.refdata_assembly_dir,
                              args.input_vcf,
                              args.validated_vcf,
                              args.input_cna,
                              args.input_rna_fusion,
                              args.input_rna_exp,
                              args.input_cpsr,
                              args.tumor_dp_tag,
                              args.tumor_af_tag,
                              args.control_dp_tag,
                              args.control_af_tag,
                              args.call_conf_tag,
                              args.exclude_hom_germline,
                              args.exclude_het_germline,
                              args.panel_normal_vcf,
                              args.retained_info_tags,
                              args.tumor_only,
                              args.sample_id,
                              args.build,
                              args.keep_uncompressed,
                              args.output_dir,
                              args.debug)
    if ret != 0:
        sys.exit(1)

def simplify_vcf(input_vcf, validated_vcf, vcf_obj, output_dir, sample_id, keep_uncompressed, logger, debug):
    """
    input_vcf: path to input VCF
    validated_vcf: path to validated VCF
    vcf_obj: parsed cyvcf2 object
    Function that performs the following on the validated input VCF:
    1. Strip of any genotype data
    2. If VCF has variants with multiple alternative alleles ("multiallelic", e.g. 'A,T'),
       these are decomposed into variants with a single alternative allele
    3. Final VCF file is sorted and indexed (bgzip + tabix)
    """

    random_id = random_id_generator(15)

    vcf_filtered     = os.path.join(output_dir, f'{sample_id}.pcgr_validate.filtered.{random_id}.vcf')
    vt_decompose_log = os.path.join(output_dir, f'{sample_id}.pcgr_validate.vt_decompose.{random_id}.log')

    ## Check for multiallelic sites
    multiallelic_list = detect_multiallelic_sites(vcf_obj)

    # Sort, chromosome-filter, strip FORMAT/genotype columns — temp files managed inside
    bgzip_sort_filter_vcf(input_vcf, vcf_filtered, output_dir, sample_id, 'pcgr', logger, debug)

    if multiallelic_list:
        logger.warning(f"There were {len(multiallelic_list)} multiallelic site(s) detected. Showing (up to) the first 100:")
        print('----')
        print(', '.join(multiallelic_list[:100]))
        print('----')
        logger.info('Decomposing multi-allelic sites in input VCF file using \'vt decompose\'')
        command_decompose = f'vt decompose -s {vcf_filtered} > {validated_vcf} 2> {vt_decompose_log}'
        check_subprocess(logger, command_decompose, debug)
    else:
        logger.info('All sites seem to be decomposed - skipping decomposition of multiallelic sites')
        check_subprocess(logger, f'cp {vcf_filtered} {validated_vcf}', debug)

    # need to keep uncompressed copy for vcf2maf.pl if selected
    bgzip_cmd = f"bgzip -cf {validated_vcf} > {validated_vcf}.gz" if keep_uncompressed else f"bgzip -f {validated_vcf}"
    check_subprocess(logger, bgzip_cmd, debug)
    check_subprocess(logger, f'tabix -p vcf {validated_vcf}.gz', debug)

    if os.path.exists(f'{validated_vcf}.gz') and os.path.getsize(f'{validated_vcf}.gz') > 0:
        vcf_check = VCF(f'{validated_vcf}.gz')
        i = 0
        for rec in vcf_check:
            i = i + 1
        if len(vcf_check.seqnames) == 0 or i == 0:
            logger.info('')
            logger.info("Input VCF contains NO valid variants on autosomal/sex chromosomes after VCF cleaning - quitting workflow")
            logger.info('')
            if not debug:
                remove_file(vcf_filtered)
                remove_file(vt_decompose_log)
            exit(1)

    if not debug:
        remove_file(vcf_filtered)
        remove_file(vt_decompose_log)

def validate_pcgr_input(refdata_assembly_dir,
                        input_vcf,
                        validated_vcf,
                        input_cna,
                        input_rna_fusion,
                        input_rna_expression,
                        input_cpsr,
                        tumor_dp_tag,
                        tumor_af_tag,
                        control_dp_tag,
                        control_af_tag,
                        call_conf_tag,
                        exclude_hom_germline,
                        exclude_het_germline,
                        panel_normal_vcf,
                        retained_info_tags,
                        tumor_only,
                        sample_id,
                        build,
                        keep_uncompressed,
                        output_dir,
                        debug):
    """
    Function that checks the format of input files to PCGR
        - VCF file with somatic SNVs/InDels - mandatory
        - Tab-separated values file with somatic copy number segments - optional
        - Tab-separated values file with RNA fusion variants - optional
        - Tab-separated values file with RNA expression values - optional
        - Tab-separated values file with CPSR-classified germline mutations - optional
    Function performs the following checks:
    1. No INFO annotation tags in the input VCF coincides with those generated by PCGR
    2. Provided columns for tumor/normal coverage and allelic depths are found in VCF
    3. Provided retained VCF INFO tags are present in VCF file
    4. If VCF have variants with multiple alternative alleles (e.g. 'A,T') run vt decompose
    5. panel-of-normals VCF adheres to the required format (PANEL_OF_NORMALS INFO tag in header)
    6. Any genotype data from VCF input file is stripped, and the resulting VCF file is sorted and indexed (bgzip + tabix)
    7. Check that copy number segment file has required columns and correct data types (and range)
    8. Check that RNA fusion variant file has required columns and correct data types
    9. Check that RNA expression file has required columns and correct data types
    10. Check that germline mutation file has required columns and correct data types
    """
    logger = getlogger('pcgr-validate-input-arguments')

    if not input_vcf == 'None':

        vcf_object = VCF(input_vcf, gts012=True)

        ## Check that VCF does not already contain INFO tags that will be appended through PCGR annotation
        vcf_tags_pcgr = collect_workflow_info_tags(refdata_assembly_dir, 'pcgr', logger)

        ## Check that no INFO annotation tags in the input VCF coincides with those generated by PCGR
        tag_check = check_existing_vcf_info_tags(vcf_object, vcf_tags_pcgr, logger)
        if tag_check == -1:
            return -1

        if retained_info_tags != "None":
            custom_check = check_retained_vcf_info_tags(vcf_object, retained_info_tags, logger)
            if custom_check == -1:
                return -1

        ## Check whether specified tags for depth/allelic fraction are properly defined in VCF
        allelic_support_check = vcf.check_format_ad_dp_tags(
            input_vcf, tumor_dp_tag, tumor_af_tag, control_dp_tag,
            control_af_tag, call_conf_tag, exclude_hom_germline,
            exclude_het_germline, tumor_only, logger)
        if allelic_support_check == -1:
            return -1

        ## Simplify VCF - remove multiallelic variants
        simplify_vcf(input_vcf, validated_vcf, vcf_object, output_dir, sample_id, keep_uncompressed, logger, debug)


    ## Validate panel-of-normals VCF is provided
    if not panel_normal_vcf == "None":
        valid_panel_normals = is_valid_pon_vcf(panel_normal_vcf, logger)
        if valid_panel_normals == -1:
            return -1

    ## Check whether file with copy number aberration segments is properly formatted
    if not input_cna == 'None':
        valid_cna = is_valid_cna(input_cna, logger)
        if valid_cna == -1:
            return -1

    ## Check whether file with classified germline calls is properly formatted
    if not input_cpsr == 'None':
        valid_germline = is_valid_germline(input_cpsr, build, logger)
        if valid_germline == -1:
            return -1

    ## Check whether file with RNA fusion variants is properly formatted
    if not input_rna_fusion == 'None':
        valid_rna_fusion = is_valid_rna_fusion(input_rna_fusion, logger)
        if valid_rna_fusion == -1:
            return -1

    ## Check whether file with RNA gene expression data is properly formatted
    if not input_rna_expression == 'None':
        valid_rna_expression = is_valid_rna_expression(input_rna_expression, logger)
        if valid_rna_expression == -1:
            return -1

    return 0

if __name__=="__main__":
    __main__()
