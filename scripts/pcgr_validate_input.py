#!/usr/bin/env python

import csv
import re
import argparse
import os
import logging
import sys
import pandas as np
from cyvcf2 import VCF

from pcgr import annoutils, utils, vcf, cna
from pcgr.annoutils import read_infotag_file
from pcgr import utils
from pcgr.utils import error_message, check_subprocess, random_id_generator

def __main__():

    parser = argparse.ArgumentParser(description='Verify input data for PCGR')
    parser.add_argument('pcgr_dir',help='PCGR base directory with accompanying data directory')
    parser.add_argument('input_vcf', help='VCF input file with somatic (tumor) query variants (SNVs/InDels)')
    parser.add_argument('validated_vcf', help="Validated VCF file with somatic (tumor) query variants (SNVs/InDels)")
    parser.add_argument('input_cna', help='Somatic (tumor) copy number query segments (tab-separated values)')
    parser.add_argument('input_rna_fusion', help='Somatic (tumor) RNA fusion variants (tab-separated values)')
    parser.add_argument('input_rna_exp', help='Somatic/ (tumor) gene expression estimates (tab-separated values)')
    parser.add_argument('panel_normal_vcf',help="VCF file with germline calls from panel of normals")
    parser.add_argument('tumor_only',type=int, default=0,choices=[0,1],help="Tumor only sequencing")
    parser.add_argument('genome_assembly',help='grch37 or grch38')
    parser.add_argument('sample_id',help='PCGR sample_name')
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

    ret = validate_pcgr_input(args.pcgr_dir,
                              args.input_vcf,
                              args.validated_vcf,
                              args.input_cna,
                              args.input_rna_fusion,
                              args.input_rna_exp,
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
                              args.genome_assembly,
                              args.sample_id,
                              args.keep_uncompressed,
                              args.output_dir,
                              args.debug)
    if ret != 0:
        sys.exit(1)

def is_valid_rna_fusion(rna_fusion_file, logger):
    """
    Function that checks whether the RNA fusion transcript file adheres to the correct format
    """
    rna_fusion_reader = csv.DictReader(open(rna_fusion_file,'r'), delimiter='\t')
    ## check that required columns are present
    if not ('GeneA' in rna_fusion_reader.fieldnames and 'GeneB' in rna_fusion_reader.fieldnames and 'Confidence' in rna_fusion_reader.fieldnames):
        err_msg = "RNA fusion file (" + str(rna_fusion_file) + ") is missing required column(s): 'Gene1', 'Gene2', or  'Confidence'\n. Column names present in file: " + str(rna_fusion_reader.fieldnames)
        return error_message(err_msg, logger)

    rna_fusion_dataframe = np.read_csv(rna_fusion_file, sep="\t")
    if rna_fusion_dataframe.empty is True:
        err_msg = 'RNA fusion file is empty - contains NO fusions'
        return error_message(err_msg, logger)
    if not rna_fusion_dataframe['Gene1'].dtype.kind in 'O': ## check that 'Gene1' is of type object
        err_msg = "'Gene1' column of RNA fusion file cannot not be of type '" + str(rna_fusion_dataframe['Gene1'].dtype) + "'"
        return error_message(err_msg, logger)
    if not rna_fusion_dataframe['Gene2'].dtype.kind in 'O': ## check that 'Gene2' is of type object
        err_msg = "'Gene2' column of RNA fusion file cannot not be of type '" + str(rna_fusion_dataframe['Gene2'].dtype) + "'"
        return error_message(err_msg, logger)
    if not rna_fusion_dataframe['Confidence'].dtype.kind in 'O': ## check that 'Confidence' is of type object
        err_msg = "'Confidence' column of RNA fusion file cannot not be of type '" + str(rna_fusion_dataframe['Confidence'].dtype) + "'"
        return error_message(err_msg, logger)

    observed_variants = {}
    for rec in rna_fusion_reader:
        if not (rec['Confidence'] == 'high' or rec['Confidence'] == 'medium' or rec['Confidence'] == 'low'): ## check that 'Confidence' column harbor permitted values
            err_msg = "Confidence column contains non-permitted values - only 'high','medium', or 'low' permitted. Value entered was " + str(rec['Confidence'])
            return error_message(err_msg, logger)

        variant_key = str(rec['Gene1']) + "_" + str(rec['Gene2'])
        if variant_key in observed_variants.keys():
            err_msg = "Duplicate entry in RNA fusion variants: " + str(variant_key) + " is found in multiple rows"
            return error_message(err_msg, logger)
        observed_variants[variant_key] = 1


    logger.info('RNA fusion file (' + str(rna_fusion_file) + ') adheres to the correct format')
    return 0

def is_valid_rna_expression(rna_exp_file, logger):
    """
    Function that checks whether the RNA expression file adheres to the correct format
    """
    rna_exp_reader = csv.DictReader(open(rna_exp_file,'r'), delimiter='\t')
    ## check that required columns are present
    if not ('Gene' in rna_exp_reader.fieldnames and 'TPM' in rna_exp_reader.fieldnames and 'Log2FC' in rna_exp_reader.fieldnames and 'PAdj' in rna_exp_reader.fieldnames and 'DiffExp' in rna_exp_reader.fieldnames):
        err_msg = "RNA fusion file (" + str(rna_exp_file) + ") is missing required column(s): 'Gene', 'TPM', 'Log2FC','PAdj', or 'DiffExp'\n. Column names present in file: " + str(rna_exp_reader.fieldnames)
        return error_message(err_msg, logger)

    rna_exp_dataframe = np.read_csv(rna_exp_file, sep="\t")
    if rna_exp_dataframe.empty is True:
        err_msg = 'RNA gene expression file is empty - contains NO gene expression estimates'
        return error_message(err_msg, logger)
    if not rna_exp_dataframe['Gene'].dtype.kind in 'O': ## check that 'Gene' is of type object
        err_msg = "'Gene' column of RNA expression file cannot not be of type '" + str(rna_exp_dataframe['Gene'].dtype) + "'"
        return error_message(err_msg, logger)
    if not rna_exp_dataframe['TPM'].dtype.kind in 'if': ## check that 'TPM' is of type object
        err_msg = "'TPM' column of RNA expression file cannot not be of type '" + str(rna_exp_dataframe['TPM'].dtype) + "'"
        return error_message(err_msg, logger)
    if not rna_exp_dataframe['Log2FC'].dtype.kind in 'if': ## check that 'LogFC' is of type object
        err_msg = "'Log2FC' column of RNA expression file cannot not be of type '" + str(rna_exp_dataframe['Log2FC'].dtype) + "'"
        return error_message(err_msg, logger)
    if not rna_exp_dataframe['PAdj'].dtype.kind in 'if': ## check that 'PAdj' is of type object
        err_msg = "'TPM' column of RNA expression file cannot not be of type '" + str(rna_exp_dataframe['PAdj'].dtype) + "'"
        return error_message(err_msg, logger)
    if not rna_exp_dataframe['DiffExp'].dtype.kind in 'O': ## check that 'DiffExp' is of type object
        err_msg = "'DiffExp' column of RNA expression file cannot not be of type '" + str(rna_exp_dataframe['DiffExp'].dtype) + "'"
        return error_message(err_msg, logger)

    for rec in rna_exp_reader:
        if not (rec['DiffExp'] == 'over' or rec['DiffExp'] == 'under' or rec['DiffExp'] == 'NS'): ## check that 'DiffExp' column harbors permitted values
            err_msg = "Confidence column contains non-permitted values - only 'over','under', or 'NS' permitted. Value entered was " + str(rec['DiffExp'])
            return error_message(err_msg, logger)

        if not (rec['TPM'] >= 0):
            err_msg = "'TPM' column cannot contain negative values - value was " + str(rec['TPM'])
            return error_message(err_msg, logger)
        if not (rec['PAdj'] >= 0):
            err_msg = "'PAdj' column (adjusted p-value from differential expression testing) cannot contain negative values - value was " + str(rec['PAdj'])
            return error_message(err_msg, logger)

    logger.info('RNA expression file (' + str(rna_exp_file) + ') adheres to the correct format')
    return 0


def validate_panel_normal_vcf(vcf, logger):
    """
    Function that checks the INFO tags in the panel of normal VCF for the presense of 'PANEL_OF_NORMAL' (logical tag)
    If any coinciding tags, an error will be returned
    """

    vcf = VCF(vcf)
    ret = -1
    for e in vcf.header_iter():
        header_element = e.info()
        if 'ID' in header_element.keys() and 'HeaderType' in header_element.keys():
            if header_element['HeaderType'] == 'INFO' and header_element['Type'] == 'Flag':
                if header_element['ID'] == 'PANEL_OF_NORMALS':
                    logger.info('Found \'PANEL_OF_NORMALS\' INFO flag in the VCF header section of the of panel of normals VCF file')
                    ret = 1

    if ret == -1:
        err_msg = 'INFO flag \'PANEL_OF_NORMALS\' is missing from the panel of normal VCF header'
        return error_message(err_msg, logger)

    return ret




def simplify_vcf(input_vcf, validated_vcf, vcf, output_dir, sample_id, keep_uncompressed, logger, debug):
    """
    input_vcf: path to input VCF
    validated_vcf: path to validated VCF
    vcf: parsed cyvcf2 object
    Function that performs the following on the validated input VCF:
    1. Strip of any genotype data
    2. If VCF has variants with multiple alternative alleles ("multiallelic", e.g. 'A,T'), 
       these are decomposed into variants with a single alternative allele
    3. Final VCF file is sorted and indexed (bgzip + tabix)
    """

    random_id = random_id_generator(15) 

    temp_files = {}
    temp_files['vcf_1'] = \
        os.path.join(output_dir, f'{sample_id}.pcgr_validate.bcftools.{random_id}_1.vcf')
    temp_files['vcf_2'] = \
        os.path.join(output_dir, f'{sample_id}.pcgr_validate.bcftools.{random_id}_2.vcf.gz')
    temp_files['vcf_3'] = \
        os.path.join(output_dir, f'{sample_id}.pcgr_validate.bftools.{random_id}_3.vcf.gz')
    bcftools_simplify_log = \
        os.path.join(output_dir, f'{sample_id}.pcgr_validate.bcftools.{random_id}.log')
    vt_decompose_log = \
        os.path.join(output_dir, f'{sample_id}.pcgr_validate.vt_decompose.{random_id}.log')

    multiallelic_list = list()
    for rec in vcf:
        POS = rec.start + 1
        alt = ",".join(str(n) for n in rec.ALT)
        if len(rec.ALT) > 1:
            variant_id = f"{rec.CHROM}:{POS}_{rec.REF}->{alt}"
            multiallelic_list.append(variant_id)

    logger.info('Extracting variants on autosomal/sex/mito chromosomes only (1-22,X,Y, M/MT) with bcftools')
    # bgzip + tabix required for sorting
    cmd_vcf1 = f'bcftools view {input_vcf} | bgzip -cf > {temp_files["vcf_2"]} && tabix -p vcf {temp_files["vcf_2"]} && ' + \
        f'bcftools sort --temp-dir {output_dir} -Oz {temp_files["vcf_2"]} > {temp_files["vcf_3"]} 2> {bcftools_simplify_log} && ' + \
        f'tabix -p vcf {temp_files["vcf_3"]}'
    # Keep only autosomal/sex/mito chrom (handle hg38 and hg19), remove FORMAT metadata lines, keep cols 1-8, sub chr prefix
    chrom_to_keep = [str(x) for x in [*range(1,23), 'X', 'Y', 'M', 'MT']]
    chrom_to_keep = ','.join([*['chr' + chrom for chrom in chrom_to_keep], *[chrom for chrom in chrom_to_keep]])
    cmd_vcf2 = f'bcftools view --regions {chrom_to_keep} {temp_files["vcf_3"]} | egrep -v \'^##FORMAT=\' ' + \
        f'| cut -f1-8 | sed \'s/^chr//\' > {temp_files["vcf_1"]}'

    check_subprocess(logger, cmd_vcf1, debug)
    check_subprocess(logger, cmd_vcf2, debug)

    if multiallelic_list:
        logger.warning(f"There were {len(multiallelic_list)} multiallelic sites detected. Showing (up to) the first 100:")
        print('----')
        print(', '.join(multiallelic_list[:100]))
        print('----')
        logger.info('Decomposing multi-allelic sites in input VCF file using \'vt decompose\'')
        command_decompose = f'vt decompose -s {temp_files["vcf_1"]} > {validated_vcf} 2> {vt_decompose_log}'
        check_subprocess(logger, command_decompose, debug)
    else:
        logger.info('All sites seem to be decomposed - skipping decomposition!')
        check_subprocess(logger, f'cp {temp_files["vcf_1"]} {validated_vcf}', debug)

    # need to keep uncompressed copy for vcf2maf.pl if selected
    bgzip_cmd = f"bgzip -cf {validated_vcf} > {validated_vcf}.gz" if keep_uncompressed else f"bgzip -f {validated_vcf}"
    check_subprocess(logger, bgzip_cmd, debug)
    check_subprocess(logger, f'tabix -p vcf {validated_vcf}.gz', debug)

    if os.path.exists(f'{validated_vcf}.gz') and os.path.getsize(f'{validated_vcf}.gz') > 0:
        vcf = VCF(f'{validated_vcf}.gz')
        i = 0
        for rec in vcf:
            i = i + 1
        if len(vcf.seqnames) == 0 or i == 0:
            logger.info('')
            logger.info("Input VCF contains NO valid variants after VCF cleaning - quitting workflow")
            logger.info('')
            exit(1)

    if not debug:
        utils.remove(temp_files["vcf_1"])
        utils.remove(temp_files["vcf_2"])
        utils.remove(temp_files["vcf_3"])
        utils.remove(temp_files["vcf_2"] + str('.tbi'))
        utils.remove(temp_files["vcf_3"] + str('.tbi'))
        utils.remove(bcftools_simplify_log)
        utils.remove(vt_decompose_log)

def validate_pcgr_input(pcgr_directory,
                        input_vcf,
                        validated_vcf,
                        input_cna,
                        input_rna_fusion,
                        input_rna_expression,
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
                        genome_assembly,
                        sample_id,
                        keep_uncompressed,
                        output_dir,
                        debug):
    """
    Function that reads the input files to PCGR (VCF file and Tab-separated values file with copy number segments) and performs the following checks:
    1. No INFO annotation tags in the query VCF coincides with those generated by PCGR
    2. Provided columns for tumor/normal coverage and allelic depths are found in VCF
    3. Provided retained VCF INFO tags are present in VCF file
    4. If VCF have variants with multiple alternative alleles (e.g. 'A,T') run vt decompose
    5. panel-of-normals VCF adheres to the required format (PANEL_OF_NORMALS INFO tag in header)
    6. Any genotype data from VCF input file is stripped, and the resulting VCF file is sorted and indexed (bgzip + tabix)
    7. Check that copy number segment file has required columns and correct data types (and range)
    8. Check that RNA fusion variant file has required columns and correct data types
    9. Check that RNA expression file has required columns and correct data types
    """
    logger = utils.getlogger('pcgr-validate-input-arguments')

    # if panel_normal_vcf == "None" and tumor_only == 1 and config_options['tumor_only']['exclude_pon'] is True:
    #    logger.warning('Panel-of-normals VCF is not present - exclusion of calls found in panel-of-normals will be ignored')

    if not input_vcf == 'None':

        vcf_object = VCF(input_vcf)

        ## Check that VCF does not contain INFO tags that will be appended with PCGR annotation
        vcf_infotags = {}
        vcf_infotags['pcgr'] = read_infotag_file(
            os.path.join(pcgr_directory,'data',genome_assembly, 'vcf_infotags_other.tsv'), scope = "pcgr")
        vcf_infotags['vep'] = read_infotag_file(
            os.path.join(pcgr_directory,'data',genome_assembly, 'vcf_infotags_vep.tsv'), scope = "vep")
        vcf_infotags['pcgr'].update(vcf_infotags['vep'])
        vcf_tags_pcgr = vcf_infotags['pcgr']
        
        tag_check = vcf.check_existing_vcf_info_tags(vcf_object, vcf_tags_pcgr, logger)
        if tag_check == -1:
            return -1

        if retained_info_tags != "None":
            custom_check = vcf.check_retained_vcf_info_tags(vcf_object, retained_info_tags, logger)
            if custom_check == -1:
                return -1

        ## Check whether specified tags for depth/allelic fraction are properly defined in VCF
        
        allelic_support_check = vcf.check_format_ad_dp_tags(
            vcf_object, tumor_dp_tag, tumor_af_tag, control_dp_tag,
            control_af_tag, call_conf_tag, exclude_hom_germline,
            exclude_het_germline, tumor_only, logger)
        if allelic_support_check == -1:
            return -1

        ## Simplify VCF - remove multiallelic variants
        simplify_vcf(input_vcf, validated_vcf, vcf_object, output_dir, sample_id, keep_uncompressed, logger, debug)


    ## Validate panel-of-normals VCF is provided
    if not panel_normal_vcf == "None":
        valid_panel_normals = validate_panel_normal_vcf(panel_normal_vcf, logger)
        if valid_panel_normals == -1:
            return -1

    ## Check whether file with copy number aberration segments is properly formatted
    if not input_cna == 'None':
        valid_cna = cna.is_valid_cna2(input_cna, logger)
        if valid_cna == -1:
            return -1

    ## Check whether file with RNA fusion variants is properly formatted
    if not input_rna_fusion == 'None':
        valid_rna_fusion = is_valid_rna_fusion(input_rna_fusion, logger)
        if valid_rna_fusion == -1:
            return -1

    ## Check whether file with RNA fusion variants is properly formatted
    if not input_rna_expression == 'None':
        valid_rna_expression = is_valid_rna_expression(input_rna_expression, logger)
        if valid_rna_expression == -1:
            return -1

    return 0

if __name__=="__main__":
    __main__()
