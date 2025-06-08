#!/usr/bin/env python

from dgg_engine import pcgr_vars, utils
from dgg_engine.utils import getlogger, error_message, warn_message, check_file_exists, check_tabix_file
import os


def verify_args(arg_dict):

    logger = getlogger("pcgr-validate-arguments-input-a")
    
    # Check the existence of required arguments    
    verify_required_args(arg_dict, logger)

    # Optional arguments
    if arg_dict['expression_sim'] == True and arg_dict['input_rna_exp'] is None:
        warn_msg = f"RNA expression similarity analysis can only be performed if --input_rna_exp is set (--input_rna_exp = {arg_dict['input_rna_exp']})."
        warn_message(warn_msg, logger)
        arg_dict['expression_sim'] = False
    
    if not arg_dict['expression_sim_db'] is None:
        expression_sim_db_sources = str(arg_dict['expression_sim_db']).split(',')
        permitted_exp_sources = pcgr_vars.EXPRESSION_DB_SOURCES
        num_permitted_sources = 0
        for v in expression_sim_db_sources:
            if v in permitted_exp_sources:
                num_permitted_sources += 1

        if num_permitted_sources != len(expression_sim_db_sources):
            err_msg = (
                f"Option '--expression_sim_db' = {str(arg_dict['expression_sim_db'])} is formatted incorrectly, should be " 
                f"a comma-separated string of the following values: '{' '.join(pcgr_vars.EXPRESSION_DB_SOURCES)}'"
            )
            error_message(err_msg, logger)

    # check that tumor primary site/type is set correctly (integer between 0 and 30)
    if arg_dict['tsite'] > max(pcgr_vars.tsites.keys()) or arg_dict['tsite'] < 0:
        err_msg = f"Tumor type code ('--tumor_site' = {arg_dict['tsite']}) must be within [0, {max(pcgr_vars.tsites.keys())}]"
        error_message(err_msg, logger)

    # check that tumor purity and tumor ploidy is set correctly
    if not arg_dict['tumor_purity'] is None:
        if not (arg_dict['tumor_purity'] > 0 and arg_dict['tumor_purity'] <= 1):
            err_msg = f"Tumor purity value ('--tumor_purity' = {arg_dict['tumor_purity']}) must be within (0, 1]"
            error_message(err_msg, logger)

    if not arg_dict['tumor_ploidy'] is None:
        if not arg_dict['tumor_ploidy'] > 0:
            err_msg = f"Tumor ploidy value ('--tumor_ploidy' = {arg_dict['tumor_ploidy']}) must be > 0"
            error_message(err_msg, logger)

    # check that minimum/maximum depth/allelic fractions are set correctly
    if int(arg_dict['tumor_dp_min']) < 0:
        err_msg = f"Minimum sequencing depth tumor ('tumor_dp_min' = {arg_dict['tumor_dp_min']}) must be >= 0"
        error_message(err_msg, logger)

    if float(arg_dict['tumor_af_min']) < 0 or float(arg_dict['tumor_af_min']) > 1:
        err_msg = f"Minimum AF tumor ('tumor_af_min' = {arg_dict['tumor_af_min']}) must be within [0, 1]"
        error_message(err_msg, logger)

    if int(arg_dict['control_dp_min']) < 0:
        err_msg = f"Minimum sequencing depth control ('control_dp_min' = {arg_dict['control_dp_min']}) must be >= 0"
        error_message(err_msg, logger)

    if float(arg_dict['control_af_max']) < 0 or float(arg_dict['control_af_max']) > 1:
        err_msg = f"Maximum AF control ('control_af_max' = {arg_dict['control_af_max']}) must be within [0, 1]"
        error_message(err_msg, logger)
    
    
      # TMB: check that minimum/maximum depth/allelic fractions are set correctly
    if int(arg_dict['tmb_dp_min']) < 0:
        err_msg = f"Minimum sequencing depth tumor - TMB calculation ('tmb_dp_min' = {arg_dict['tmb_dp_min']}) must be >= 0"
        error_message(err_msg, logger)
        
    if int(arg_dict['tmb_dp_min']) > 0 and (int(arg_dict['tmb_dp_min']) < int(arg_dict['tumor_dp_min'])):
        err_msg = f"Minimum sequencing depth (tumor) for TMB calculation ('tmb_dp_min' = {str(arg_dict['tmb_dp_min'])}) must be "
        err_msg += f"greater or equal to minimum sequencing depth tumor {str(arg_dict['tumor_dp_min'])} (i.e. filter for variant inclusion in report)"
        error_message(err_msg, logger)

    if float(arg_dict['tmb_af_min']) < 0 or float(arg_dict['tmb_af_min']) > 1:
        err_msg = f"Minimum AF (tumor) for TMB calculation ('tmb_af_min' = {arg_dict['tmb_af_min']}) must be within [0, 1]"
        error_message(err_msg, logger)
        
    if float(arg_dict['tmb_af_min']) > 0 and (float(arg_dict['tmb_af_min']) < float(arg_dict['tumor_af_min'])):
        err_msg = f"Minimum AF (tumor) for TMB calculation ('tmb_af_min' = {str(arg_dict['tmb_af_min'])}) must be "
        err_msg += f"greater or equal to minimum AF tumor {str(arg_dict['tumor_dp_min'])} (i.e. filter for variant inclusion in report)"
        error_message(err_msg, logger)

    # Check that coding target size region of sequencing assay is set correctly
    if float(arg_dict['effective_target_size_mb']) < 0 or float(arg_dict['effective_target_size_mb']) > float(pcgr_vars.CODING_EXOME_SIZE_MB):
        err_msg = (
            f"Coding target size region in Mb ('--effective_target_size_mb' = {arg_dict['effective_target_size_mb']}) is not "
            f"positive or larger than the approximate size of the coding human genome ({float(pcgr_vars.CODING_EXOME_SIZE_MB)} Mb))")
        error_message(err_msg, logger)
    if float(arg_dict['effective_target_size_mb']) < 1:
        warn_msg = (
            f"Coding target size region in Mb ('--effective_target_size_mb' = {arg_dict['effective_target_size_mb']}) must be "
            "greater than 1 Mb for mutational burden estimate to be robust")
        warn_message(warn_msg, logger)
    if float(arg_dict['effective_target_size_mb']) < float(pcgr_vars.CODING_EXOME_SIZE_MB) and arg_dict['assay'] != 'TARGETED':
        warn_msg = (
            f"Coding target size region in Mb ('--effective_target_size_mb' = {arg_dict['effective_target_size_mb']}) is less than ",
            f"the default for WES/WGS ({float(pcgr_vars.CODING_EXOME_SIZE_MB)} Mb), assay must be set to 'TARGETED'")
        warn_message(warn_msg, logger)

    # if assay is targeted or mode is Tumor-Only, MSI prediction will not be performed/switched off
    assay_type = 'Tumor-Control'
    if arg_dict['estimate_msi'] is True and (arg_dict['assay'] == 'TARGETED' or arg_dict['tumor_only'] is True):
        if arg_dict['tumor_only'] is True:
            assay_type = 'Tumor-Only'
        warn_msg = f"MSI status prediction can be applied for WGS/WES tumor-control assays only (query type: {arg_dict['assay']}|{assay_type}) - analysis will be omitted"
        warn_message(warn_msg, logger)
        arg_dict['estimate_msi'] = 0

    # minimum number of mutations required for mutational signature re-fitting cannot be less than 100 
    # somewhat arbitrary min threshold), recommended value is 200)
    if int(arg_dict['min_mutations_signatures']) < int(pcgr_vars.RECOMMENDED_N_MUT_SIGNATURE):
        if int(arg_dict['min_mutations_signatures']) < int(pcgr_vars.MINIMUM_N_MUT_SIGNATURE):
            err_msg = (
                f"Minimum number of mutations required for mutational signature analysis ('--min_mutations_signatures' = "
                f"{arg_dict['min_mutations_signatures']}) must be >= {pcgr_vars.MINIMUM_N_MUT_SIGNATURE}")
            error_message(err_msg, logger)
        warn_msg = (
            f"Minimum number of mutations required for mutational signature analysis ('--min_mutations_signatures' "
            f"= {arg_dict['min_mutations_signatures']}) is less than the recommended number (n = {pcgr_vars.RECOMMENDED_N_MUT_SIGNATURE})")
        warn_message(warn_msg, logger)
        

    if float(arg_dict['prevalence_reference_signatures']) > int(pcgr_vars.MAX_SIGNATURE_PREVALENCE) or \
        float(arg_dict['prevalence_reference_signatures']) < 0:
        err_msg = (
            f"Prevalence of reference signatures (percent) must be above zero and less than {pcgr_vars.MAX_SIGNATURE_PREVALENCE} "
            f"'--prevalence_reference_signatures' = {arg_dict['prevalence_reference_signatures']}")
        error_message(err_msg, logger)

    # if MSI status is to be estimated, mutational burden must be turned on
    if arg_dict['estimate_msi'] is True and arg_dict['estimate_tmb'] is False:
        err_msg = "Prediction of MSI status ('--estimate_msi') requires mutational burden analysis ('--estimate_tmb')"
        error_message(err_msg, logger)

    if arg_dict['tumor_only'] is True:
        
        if arg_dict['control_dp_tag'] is not None and arg_dict['control_dp_tag'] != "_NA_":
            err_msg = f"Option '--tumor_only' does not allow '--control_dp_tag' option to be set ({arg_dict['control_dp_tag']})"
            error_message(err_msg, logger)
        
        if arg_dict['control_af_tag'] is not None and arg_dict['control_af_tag'] != "_NA_":
            err_msg = f"Option '--tumor_only' does not allow '--control_af_tag' option to be set ({arg_dict['control_af_tag']})"
            error_message(err_msg, logger)
        
        for t in ['exclude_likely_het_germline', 'exclude_likely_hom_germline']:
            if arg_dict[t]:
                if arg_dict['tumor_af_tag'] == "_NA_":
                    err_msg = f"Option '--{t}' requires '--tumor_af_tag' option to be set"
                    error_message(err_msg, logger)

        # Emit warning if panel-of-normals VCF is not present and exclude_pon is set
        if arg_dict['pon_vcf'] is None and arg_dict['exclude_pon'] is True:
            warn_msg = "Panel-of-normals VCF is NOT provided ('--pon_vcf') - exclusion of calls found in panel-of-normals ('--exclude_pon') will be ignored"
            warn_message(warn_msg, logger)
            arg_dict['exclude_pon'] = False

        # Emit warnings that mutational burden and mutational signatures are less accurate for assays with tumor-only data
        if arg_dict['estimate_tmb'] is True:
            warn_msg = "Estimation of mutational burden in tumor-only mode is suboptimal - results must be interpreted with caution"
            warn_message(warn_msg, logger)
        if arg_dict['estimate_signatures'] is True:
            warn_msg = "Estimation of mutational signatures in tumor-only mode is suboptimal - results must be interpreted with caution"
            warn_message(warn_msg, logger)

        for pop in ['nfe', 'fin', 'amr', 'eas', 'sas', 'asj', 'oth', 'afr', 'global']:
            tag = f'maf_gnomad_{pop}'
            if arg_dict[tag]:
                if float(arg_dict[tag]) < 0 or float(arg_dict[tag]) > 1:
                    err_msg = f"MAF threshold (tumor-only germline filter) for gnomAD (pop '{pop.upper()}') must be within the [0, 1] range, current value is {arg_dict[tag]}"
                    error_message(err_msg, logger)

    # Emit warning that mutational signature estimation is (likely) not optimal for small/targeted sequencing assays
    if arg_dict['estimate_signatures'] is True and arg_dict['assay'] == 'TARGETED':
        warn_msg = "Estimation of mutational signatures ('--estimate_signatures') is not optimal for TARGETED sequencing assays - results must be interpreted with caution"
        warn_message(warn_msg, logger)

    # Check that threshold for gains/amplifications are properly set, and that segment overlap with transcripts are set appropriately
    if arg_dict['n_copy_gain'] <= 2:
        err_msg = f"Total copy number threshold for gains/amplifications ('--n_copy_gain' = {arg_dict['n_copy_gain']}) should be > 2"
        error_message(err_msg, logger)
    if arg_dict['cna_overlap_pct'] > 100 or arg_dict['cna_overlap_pct'] <= 0:
        err_msg = f"Minimum percent overlap between copy number segment and gene transcript ('--cna_overlap_pct' = {arg_dict['cna_overlap_pct']}) must be within (0, 100]"
        error_message(err_msg, logger)

    verify_vep_options(arg_dict, logger)

    return(arg_dict)


def define_output_files(arg_dict, cpsr = False):
    """
    Define output files
    """
    logger = getlogger("define-output-files")
    # create output folder (if not already exists)
    output_dir = utils.safe_makedir(os.path.abspath(arg_dict['output_dir']))

    # Define general output prefix
    output_prefix = os.path.join(
        str(output_dir), f"{arg_dict['sample_id']}.pcgr.{arg_dict['genome_assembly']}")
    if cpsr:
        output_prefix = os.path.join(
            str(output_dir), f"{arg_dict['sample_id']}.cpsr.{arg_dict['genome_assembly']}")

    # Define output data files (absolute paths) - DGG Engine   
    output_data = {}
    output_data['prefix'] = output_prefix
    output_data['dir'] = utils.safe_makedir(os.path.abspath(arg_dict['output_dir']))
    output_data['vcf'] = f"{output_prefix}.vcf.gz"
    output_data['vcf_pass'] = f"{output_prefix}.pass.vcf.gz"
    output_data['vcf2tsv'] = f"{output_prefix}.pass.tsv.gz"
    output_data['html'] = f"{output_