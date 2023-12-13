#!/usr/bin/env python

import logging

from pcgr.utils import error_message, warn_message, check_subprocess
from cyvcf2 import VCF
from typing import Union




def get_vcf_info_tags(vcf_fname):
    vcf = VCF(vcf_fname)
    info_tags = {}
    for e in vcf.header_iter():
        header_element = e.info()
        if 'ID' in header_element.keys() and 'HeaderType' in header_element.keys():
            if header_element['HeaderType'] == 'INFO':
                info_tags[str(header_element['ID'])] = 1

    return info_tags


def print_vcf_header(vcf_fname, vcfheader_file, logger, chromline_only=False):
    if chromline_only == True:
        check_subprocess(
            logger, f'bgzip -dc {vcf_fname} | egrep \'^#\' | egrep \'^#CHROM\' >> {vcfheader_file}', debug=False)
    else:
        check_subprocess(
            logger, f'bgzip -dc {vcf_fname} | egrep \'^#\' | egrep -v \'^#CHROM\' > {vcfheader_file}', debug=False)

def detect_reserved_info_tag(tag, tag_name, logger):
    reserved_tags = ['AA', 'AC', 'AF', 'AN', 'BQ', 'CIGAR', 'DB', 'DP', 'END',
                     'H2', 'H3', 'MQ', 'MQ0', 'NS', 'SB', 'SOMATIC', 'VALIDATED', '1000G']
    if tag in reserved_tags:
        err_msg = f'Custom INFO tag ({tag_name}) needs another name - \'{tag}\' is a reserved field in the VCF specification (INFO)'
        return error_message(err_msg, logger)

    reserved_format_tags = ['GT', 'DP', 'FT', 'GL',
                            'GLE', 'GQ', 'PL', 'HQ', 'PS', 'PQ', 'EC', 'MQ']
    if tag in reserved_format_tags:
        err_msg = 'Custom INFO tag ({tag_name}) needs another name - \'{tag}\' is a reserved field in the VCF specification (FORMAT)'
        return error_message(err_msg, logger)

def check_retained_vcf_info_tags(vcf: VCF, retained_info_tags: str, logger: logging.Logger) -> int:

    """
    Function that compares the INFO tags in the query VCF and retained INFO tags set by the user as retained for output
    If any retained tag is not in query VCF, an error will be returned
    """

    tags = str(retained_info_tags).split(',')
    info_elements_query_vcf = []

    #vcf = VCF(input_vcf)
    logger.info('Checking if existing INFO tags of query VCF file matches retained INFO tags set by the user')
    ret = 1
    for e in vcf.header_iter():
        header_element = e.info()
        if 'ID' in header_element.keys() and 'HeaderType' in header_element.keys():
            if header_element['HeaderType'] == 'INFO':
                info_elements_query_vcf.append(header_element['ID'])

    for t in tags:
        if not t in info_elements_query_vcf:
            err_msg = "Retained INFO tag '" + str(t) + "' not found among INFO tags in query VCF - make sure retained VCF INFO tags are set correctly"
            error_message(err_msg, logger)
        else:
            logger.info("Retained INFO tag '" + str(t) + "' detected among INFO tags in query VCF")
        
    return ret



def check_existing_vcf_info_tags(vcf: VCF, populated_infotags: dict, logger) -> int:

    """
    Function that compares the INFO tags in the query VCF and the INFO tags generated by PCGR/CPSR
    If any coinciding tags, an error will be returned
    """

    logger.info('Checking if existing INFO tags of query VCF file coincide with PCGR/CPSR INFO tags')
    ret = 1
    for e in vcf.header_iter():
        header_element = e.info()
        if 'ID' in header_element.keys() and 'HeaderType' in header_element.keys():
            if header_element['HeaderType'] == 'INFO':
                if header_element['ID'] in populated_infotags.keys():
                    err_msg = 'INFO tag ' + str(header_element['ID']) + ' in the query VCF coincides with a VCF ' + \
                        'annotation tag produced by PCGR/CPSR - please remove or rename this tag in your query VCF'
                    return error_message(err_msg, logger)

    logger.info('No query VCF INFO tags coincide with populated INFO tags by PCGR/CPSR')
    return ret



def check_format_ad_dp_tags(vcf: VCF,
                           tumor_dp_tag: str,
                           tumor_af_tag: str,
                           control_dp_tag: str,
                           control_af_tag: str,
                           call_conf_tag: str,
                           exclude_hom_germline: bool,
                           exclude_het_germline: bool,
                           tumor_only: int,
                           logger: logging.Logger) -> Union[int, None]:

    """
    Function that checks whether the INFO tags specified for depth/allelic fraction are correctly formatted in the VCF header (i.e. Type)
    """

    found_taf_tag = 0
    found_tdp_tag = 0
    found_naf_tag = 0
    found_ndp_tag = 0
    found_call_conf_tag = 0

    detect_reserved_info_tag(tumor_dp_tag,'tumor_dp_tag', logger)
    detect_reserved_info_tag(control_dp_tag,'control_dp_tag', logger)
    detect_reserved_info_tag(tumor_af_tag,'tumor_af_tag', logger)
    detect_reserved_info_tag(control_af_tag,'control_af_tag', logger)
    detect_reserved_info_tag(call_conf_tag,'call_conf_tag', logger)

    for e in vcf.header_iter():
        header_element = e.info()
        if 'ID' in header_element.keys() and 'HeaderType' in header_element.keys():
            if header_element['HeaderType'] == 'INFO':
                if header_element['ID'] == tumor_dp_tag:
                    if header_element['Type'] == 'Integer':
                        logger.info(f'Found INFO tag for tumor variant sequencing depth (tumor_dp_tag {tumor_dp_tag}) in input VCF')
                        found_tdp_tag = 1
                    else:
                        err_msg = f'INFO tag for tumor variant sequencing depth (tumor_dp_tag {tumor_dp_tag}) is not correctly specified in input VCF (Type={header_element["Type"]}), should be Type=Integer'
                        return error_message(err_msg, logger)
                if header_element['ID'] == tumor_af_tag:
                    if header_element['Type'] == 'Float':
                        logger.info(f'Found INFO tag for tumor variant allelic fraction (tumor_af_tag {tumor_af_tag}) in input VCF')
                        found_taf_tag = 1
                    else:
                        err_msg = f'INFO tag for tumor variant allelic fraction (tumor_af_tag {tumor_af_tag}) is not correctly specified in input VCF (Type={header_element["Type"]}), should be Type=Float'
                        return error_message(err_msg, logger)
                if header_element['ID'] == control_dp_tag:
                    if header_element['Type'] == 'Integer':
                        logger.info(f'Found INFO tag for normal/control variant sequencing depth (control_dp_tag {control_dp_tag}) in input VCF')
                        found_ndp_tag = 1
                    else:
                        err_msg = f'INFO tag for normal/control variant sequencing depth (control_dp_tag {control_dp_tag}) is not correctly specified in input VCF (Type={header_element["Type"]}), should be Type=Integer'
                        return error_message(err_msg, logger)
                if header_element['ID'] == control_af_tag:
                    if header_element['Type'] == 'Float':
                        logger.info(f'Found INFO tag for normal/control allelic fraction (control_af_tag {control_af_tag}) in input VCF')
                        found_naf_tag = 1
                    else:
                        err_msg = f'INFO tag for for normal/control allelic fraction (control_af_tag {control_af_tag}) is not correctly specified in input VCF (Type={header_element["Type"]}) should be Type=Float'
                        return error_message(err_msg, logger)
                if header_element['ID'] == call_conf_tag:
                    if header_element['Type'] == 'String':
                        logger.info(f'Found INFO tag for variant call confidence (call_conf_tag {call_conf_tag}) in input VCF')
                        found_call_conf_tag = 1
                    else:
                        err_msg = f'INFO tag for variant call confidence (call_conf_tag) is not correctly specified in input VCF (Type={header_element["Type"]}), should be Type=String'
                        return error_message(err_msg, logger)


    if call_conf_tag != '_NA_' and found_call_conf_tag == 0:
        logger.warning(f"Could not find the specified call_conf_tag ('{call_conf_tag}') in INFO column of input VCF")
    if tumor_dp_tag != '_NA_' and found_tdp_tag == 0:
        logger.warning(f"Could not find the specified tumor_dp_tag ('{tumor_dp_tag}') in INFO column of input VCF")
    if tumor_af_tag != '_NA_' and found_taf_tag == 0:
        logger.warning(f"Could not find the specified tumor_af_tag ('{tumor_af_tag}') in INFO column of input VCF")
    if control_dp_tag != '_NA_' and found_ndp_tag == 0:
        logger.warning(f"Could not find the specified control_dp_tag ('{control_dp_tag}') in INFO column of input VCF")
    if control_af_tag != '_NA_' and found_naf_tag == 0:
        logger.warning(f"Could not find the specified control_af_tag ('{control_af_tag}') in INFO column of input VCF")

    if exclude_hom_germline is True and tumor_only == 1 and found_taf_tag == 0:
        logger.warning(f"Could not find the specified tumor_af_tag ('{tumor_af_tag}') in INFO column of input VCF - filtering of homozygous germline variants in tumor-only mode will be ignored")

    if exclude_het_germline is True and tumor_only == 1 and found_taf_tag == 0:
        logger.warning(f"Could not find the specified tumor_af_tag ('{tumor_af_tag}') in INFO column of input VCF - filtering of heterozygous germline variants in tumor-only mode will be ignored")


    if found_tdp_tag == 1 and found_taf_tag == 0:
        logger.warning('BOTH \' tumor_dp_tag\' AND \' tumor_af_tag\' need to be specified for use in tumor report (\'tumor_af_tag\' is missing)')

    if found_tdp_tag == 0 and found_taf_tag == 1:
        logger.warning('BOTH \'tumor_dp_tag\' AND \'tumor_af_tag\' need to be specified for use in tumor report (\'tumor_dp_tag\' is missing)')

    if found_ndp_tag == 1 and found_naf_tag == 0:
        logger.warning('BOTH \'control_dp_tag\' AND \'control_af_tag\' need to be specified for use in tumor report (\'control_af_tag\' is missing)')

    if found_ndp_tag == 0 and found_naf_tag == 1:
        logger.warning('BOTH \'control_dp_tag\' AND \'control_af_tag\' need to be specified for use in tumor report (\'control_dp_tag\' is missing)')

    ## if filtering turned on for AF-based tumor-only filtering, return error if TVAF not defined

    return 0
