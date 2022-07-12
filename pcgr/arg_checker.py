#!/usr/bin/env python

from pcgr import pcgr_vars
from pcgr.utils import getlogger, error_message, warn_message
import os


def check_args(arg_dict):

    logger = getlogger("pcgr-validate-arguments-input-a")
    # Check the existence of required arguments
    if arg_dict['pcgr_dir'] is None or not os.path.exists(arg_dict['pcgr_dir']):
        err_msg = f"Required argument '--pcgr_dir' does not exist ({arg_dict['pcgr_dir']})."
        error_message(err_msg, logger)

    if arg_dict['genome_assembly'] is None:
        err_msg = f"Required argument '--genome_assembly' has no/undefined value ({arg_dict['genome_assembly']})."
        error_message(err_msg, logger)

    if arg_dict['input_vcf'] is None:
        err_msg = f"Required argument '--input_vcf' does not exist ({arg_dict['input_vcf']})."
        error_message(err_msg, logger)

    if arg_dict['sample_id'] is None:
        err_msg = f"Required argument '--sample_id' has no/undefined value ({arg_dict['sample_id']})."
        error_message(err_msg, logger)

    if len(arg_dict['sample_id']) <= 2 or len(arg_dict['sample_id']) > 35:
        err_msg = f"Sample name identifier ('--sample_id' = {arg_dict['sample_id']}) must be between 2 and 35 characters long"
        error_message(err_msg, logger)

    # Optional arguments

    # check if input is cancer cell line, requires --tumor_only
    if arg_dict['cell_line'] and not arg_dict['tumor_only']:
        err_msg = 'Analysis of cell line (--cell_line) needs option --tumor_only'
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
    if arg_dict['tumor_dp_min'] < 0:
        err_msg = f"Minimum depth tumor ('tumor_dp_min' = {arg_dict['tumor_dp_min']}) must be >= 0"
        error_message(err_msg, logger)

    if arg_dict['tumor_af_min'] < 0 or arg_dict['tumor_af_min'] > 1:
        err_msg = f"Minimum AF tumor ('tumor_af_min' = {arg_dict['tumor_af_min']}) must be within [0, 1]"
        error_message(err_msg, logger)

    if arg_dict['control_dp_min'] < 0:
        err_msg = f"Minimum depth control ('control_dp_min' = {arg_dict['control_dp_min']}) must be >= 0"
        error_message(err_msg, logger)

    if arg_dict['control_af_max'] < 0 or arg_dict['control_af_max'] > 1:
        err_msg = f"Maximum AF control ('control_af_max' = {arg_dict['control_af_max']}) must be within [0, 1]"
        error_message(err_msg, logger)

    # Check that coding target size region of sequencing assay is set correctly
    if arg_dict['target_size_mb'] < 0 or arg_dict['target_size_mb'] > 34:
        err_msg = f"Coding target size region in Mb ('--target_size_mb' = {arg_dict['target_size_mb']}) is not positive or larger than the likely maximum size of the coding human genome (34 Mb))"
        error_message(err_msg, logger)
    if arg_dict['target_size_mb'] < 1:
        warn_msg = f"Coding target size region in Mb ('--target_size_mb' = {arg_dict['target_size_mb']}) must be greater than 1 Mb for mutational burden estimate to be robust"
        warn_message(warn_msg, logger)
    if arg_dict['target_size_mb'] < 34 and arg_dict['assay'] != 'TARGETED':
        warn_msg = f"Coding target size region in Mb ('--target_size_mb' = {arg_dict['target_size_mb']}) is less than default for WES/WGS (34Mb), assay must be set to 'TARGETED'"
        warn_message(warn_msg, logger)

    # if assay is targeted or mode is Tumor-Only, MSI prediction will not be performed/switched off
    assay_type = 'Tumor-Control'
    if arg_dict['estimate_msi_status'] is True and (arg_dict['assay'] == 'TARGETED' or arg_dict['tumor_only'] is True):
        if arg_dict['tumor_only'] is True:
            assay_type = 'Tumor-Only'
        warn_msg = f"MSI status prediction can be applied for WGS/WES tumor-control assays only (query type: {arg_dict['assay']}|{assay_type}) - analysis will be omitted"
        warn_message(warn_msg, logger)
        arg_dict['estimate_msi_status'] = 0

    # minimum number of mutations required for mutational signature reconstruction cannot be less than 100 (somewhat arbitrary lower threshold, recommended value is 200)
    if arg_dict['min_mutations_signatures'] < 200:
        warn_msg = f"Minimum number of mutations required for mutational signature analysis ('--min_mutations_signatures' = {arg_dict['min_mutations_signatures']}) is less than the recommended number (n = 200)"
        warn_message(warn_msg, logger)
        if arg_dict['min_mutations_signatures'] < 100:
            err_msg = f"Minimum number of mutations required for mutational signature analysis ('--min_mutations_signatures' = {arg_dict['min_mutations_signatures']}) must be >= 100"
            error_message(err_msg, logger)

    # if MSI status is to be estimated, mutational burden must be turned on
    if arg_dict['estimate_msi_status'] is True and arg_dict['estimate_tmb'] is False:
        err_msg = "Prediction of MSI status ('--estimate_msi_status') requires mutational burden analysis ('--estimate_tmb')"
        error_message(err_msg, logger)

    if arg_dict['tumor_only'] is True:
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

        # Emit errors when tumor-only filtering thresholds are not properly set
        for pop in ['eur', 'afr', 'amr', 'eas', 'sas', 'global']:
            tag = f'maf_onekg_{pop}'
            if arg_dict[tag]:
                if float(arg_dict[tag]) < 0 or float(arg_dict[tag]) > 1:
                    err_msg = f"MAF threshold (tumor-only germline filter) for 1000 Genomes Project (pop '{pop.upper()}') must be within the [0, 1] range, current value is {arg_dict[tag]}"
                    error_message(err_msg, logger)

        for pop in ['nfe', 'fin', 'amr', 'eas', 'sas', 'asj', 'oth', 'afr', 'global']:
            tag = f'maf_gnomad_{pop}'
            if arg_dict[tag]:
                if float(arg_dict[tag]) < 0 or float(arg_dict[tag]) > 1:
                    err_msg = f"MAF threshold (tumor-only germline filter) for gnomAD (pop '{pop.upper()}') must be within the [0, 1] range, current value is {arg_dict[tag]}"
                    error_message(err_msg, logger)

    ## tumor-only is False
    # else:
    #    for t in ["exclude_pon","exclude_likely_het_germline","exclude_likely_hom_germline","exclude_dbsnp_nonsomatic","exclude_nonexonic"]:
    #       if arg_dict[t] is True:
    #          warn_msg = "Option "--" + str(t) + "" requires "--tumor_only" option (not currently set)"
    #          warn_message(warn_msg, logger)

    # Emit warning that mutational signature estimation is (likely) not optimal for small/targeted sequencing assays
    if arg_dict['estimate_signatures'] is True and arg_dict['assay'] == 'TARGETED':
        warn_msg = "Estimation of mutational signatures ('--estimate_signatures') is not optimal for TARGETED sequencing assays - results must be interpreted with caution"
        warn_message(warn_msg, logger)

    # Check that log ratio thresholds for homozygous deletions and amplifications are properly set, and that segment overlap with transcripts are set appropriately
    if arg_dict['logr_homdel'] >= 0:
        err_msg = f"Log ratio for homozygous deletions ('--logr_homdel' = {arg_dict['logr_homdel']}) should be < 0"
        error_message(err_msg, logger)
    if arg_dict['logr_gain'] <= 0:
        err_msg = f"Log ratio for copy number gains/amplifications ('--logr_gain' = {arg_dict['logr_gain']}) should be > 0"
        error_message(err_msg, logger)
    if arg_dict['cna_overlap_pct'] > 100 or arg_dict['cna_overlap_pct'] <= 0:
        err_msg = f"Minimum percent overlap between copy number segment and gene transcript ('--cna_overlap_pct' = {arg_dict['cna_overlap_pct']}) must be within (0, 100]"
        error_message(err_msg, logger)

    # VEP options
    if arg_dict['vep_n_forks'] <= 0 or arg_dict['vep_n_forks'] > 4:
        err_msg = f"Number of forks that VEP can use during annotation ('--vep_n_forks' = {arg_dict['vep_n_forks']}must be within (0, 4]"
        error_message(err_msg, logger)

    if arg_dict['vep_buffer_size'] <= 0 or arg_dict['vep_buffer_size'] > 30000:
        err_msg = f"Internal VEP buffer size, corresponding to the number of variants that are read in to memory simultaneously ('--vep_buffer_size' = {arg_dict['vep_buffer_size']}),  must be within (0, 30000]"
        error_message(err_msg, logger)

    # Check that VEP pick criteria is formatted correctly
    if not arg_dict['vep_pick_order'] is None:
        values = str(arg_dict['vep_pick_order']).split(',')
        permitted_sources = ['canonical', 'appris', 'tsl', 'biotype', 'ccds', 'rank', 'length', 'mane']
        num_permitted_sources = 0
        for v in values:
            if v in permitted_sources:
                num_permitted_sources += 1

        if num_permitted_sources != 8:
            err_msg = (f"'--vep_pick_order' = {arg_dict['vep_pick_order']} is formatted incorrectly, should be "
                "a comma-separated string of the following values: canonical,appris,tsl,biotype,ccds,rank,length,mane")
            error_message(err_msg, logger)
    return


def verify_input_files(arg_dict):
    """
    1. Checks existence of input files/dirs (arg_dict)
    2. Checks that the data bundle is of correct date
    """
    logger = getlogger("pcgr-validate-arguments-input-b")

    input_vcf_dir = 'NA'
    input_cna_dir = 'NA'
    input_rna_fusion_dir = 'NA'
    input_cpsr_report_dir = 'NA'
    input_rna_expression_dir = 'NA'
    input_cna_plot_dir = 'NA'
    panel_normal_vcf_dir = 'NA'
    db_dir = 'NA'
    base_dir = 'NA'
    output_dir_full = 'NA'
    panel_normal_vcf_basename = 'NA'
    input_vcf_basename = 'NA'
    input_cna_basename = 'NA'
    input_rna_fusion_basename = 'NA'
    input_rna_expression_basename = 'NA'
    input_cpsr_report_basename = 'NA'
    input_cna_plot_basename = 'NA'

    arg_dict['rna_fusion_tumor'] = None
    arg_dict['rna_exp_tumor'] = None

    # check that either input vcf or cna segments exist
    if arg_dict['input_vcf'] is None and arg_dict['input_cna'] is None:
        err_msg = 'Please specifiy either a VCF input file (--input_vcf) or a copy number segment file (--input_cna)'
        error_message(err_msg, logger)

    # create output folder (if not already exists)
    output_dir_full = safe_makedir(os.path.abspath(arg_dict['output_dir']))

    # check if panel of normal VCF exist
    if not arg_dict["pon_vcf"] is None:
        if not os.path.exists(os.path.abspath(arg_dict["pon_vcf"])):
            err_msg = "Input file (" + \
                str(arg_dict["pon_vcf"]) + ") does not exist"
            error_message(err_msg, logger)

        if not (os.path.abspath(arg_dict["pon_vcf"]).endswith(".vcf.gz")):
            err_msg = "Panel of normals VCF file (" + os.path.abspath(
                arg_dict["pon_vcf"]) + ") does not have the correct file extension (.vcf.gz)"
            error_message(err_msg, logger)

        # check that tabix file exist if bgzipped files is given
        if os.path.abspath(arg_dict["pon_vcf"]).endswith(".vcf.gz"):
            tabix_file = arg_dict["pon_vcf"] + ".tbi"
            if not os.path.exists(os.path.abspath(tabix_file)):
                err_msg = "Tabix file (i.e. '.gz.tbi') is not present for the bgzipped panel of normal VCF file (" + os.path.abspath(arg_dict["pon_vcf"]) + \
                    "). Please make sure your the VCF is properly compressed and indexed (bgzip + tabix)"
                error_message(err_msg, logger)

        if arg_dict["input_vcf"] is None:
            warn_msg = "Ignoring panel of normal VCF file, --input_vcf missing"
            warn_message(warn_msg, logger)
        else:
            panel_normal_vcf_basename = os.path.basename(
                str(arg_dict["pon_vcf"]))
            panel_normal_vcf_dir = os.path.dirname(
                os.path.abspath(arg_dict["pon_vcf"]))

    # check if input vcf exists
    if not arg_dict["input_vcf"] is None:
        if not os.path.exists(os.path.abspath(arg_dict["input_vcf"])):
            err_msg = f'Input file ({arg_dict["input_vcf"]}) does not exist'
            error_message(err_msg, logger)

        if not (os.path.abspath(arg_dict["input_vcf"]).endswith(".vcf") or os.path.abspath(arg_dict["input_vcf"]).endswith(".vcf.gz")):
            err_msg = f'VCF input file ({os.path.abspath(arg_dict["input_vcf"])}) does not have the correct file extension (.vcf or .vcf.gz)'
            error_message(err_msg, logger)

        # check that tabix file exists if bgzipped file is given
        if os.path.abspath(arg_dict["input_vcf"]).endswith(".vcf.gz"):
            tabix_file = arg_dict["input_vcf"] + ".tbi"
            if not os.path.exists(os.path.abspath(tabix_file)):
                err_msg = "Tabix file (i.e. '.gz.tbi') is not present for the bgzipped VCF input file (" + os.path.abspath(arg_dict["input_vcf"]) + \
                    "). Please make sure your input VCF is properly compressed and indexed (bgzip + tabix)"
                error_message(err_msg, logger)

        input_vcf_basename = os.path.basename(str(arg_dict["input_vcf"]))
        input_vcf_dir = os.path.dirname(os.path.abspath(arg_dict["input_vcf"]))

        # if output vcf exist and overwrite not set
        output_vcf = os.path.join(str(output_dir_full), f"{arg_dict['sample_id']}.pcgr_acmg.{arg_dict['genome_assembly']}.vcf.gz")
        if os.path.exists(output_vcf) and arg_dict["force_overwrite"] is False:
            err_msg = f"Output files (e.g. {output_vcf}) already exist - please specify different sample_id or add option --force_overwrite"
            error_message(err_msg, logger)

    # check if input cna plot file exist
    # if not arg_dict["input_cna_plot"] is None:
    #    if not os.path.exists(os.path.abspath(arg_dict["input_cna_plot"])):
    #       err_msg = "Input file (" + str(arg_dict["input_cna_plot"]) + ") does not exist"
    #       error_message(err_msg,logger)
    #    if not (os.path.abspath(arg_dict["input_cna_plot"]).endswith(".png")):
    #       err_msg = "CNA segment input file (" + os.path.abspath(arg_dict["input_cna_plot"]) + ") does not have the correct file extension (.png)"
    #       error_message(err_msg,logger)
    #    if arg_dict["input_cna"] is None:
    #       err_msg = "Input a CNA plot needs to come with a CNA segment file (--input_cna is missing)"
    #       error_message(err_msg,logger)
    #    input_cna_plot_basename = os.path.basename(str(arg_dict["input_cna_plot"]))
    #    input_cna_plot_dir = os.path.dirname(os.path.abspath(arg_dict["input_cna_plot"]))

    # check if input cna segments exist
    if not arg_dict["input_cna"] is None:
        if not os.path.exists(os.path.abspath(arg_dict["input_cna"])):
            err_msg = "Input file (" + \
                str(arg_dict["input_cna"]) + ") does not exist"
            error_message(err_msg, logger)
        if not (os.path.abspath(arg_dict["input_cna"]).endswith(".tsv") or os.path.abspath(arg_dict["input_cna"]).endswith(".txt")):
            err_msg = "CNA segment input file (" + os.path.abspath(
                arg_dict["input_cna"]) + ") does not have the correct file extension (.tsv or .txt)"
            error_message(err_msg, logger)
        input_cna_basename = os.path.basename(str(arg_dict["input_cna"]))
        input_cna_dir = os.path.dirname(os.path.abspath(arg_dict["input_cna"]))

        # if output cna segments exist and overwrite not set
        output_cna_segments = os.path.join(str(output_dir_full), str(
            arg_dict["sample_id"])) + ".pcgr_acmg." + str(arg_dict["genome_assembly"]) + ".cna_segments.tsv.gz"
        if os.path.exists(output_cna_segments) and arg_dict["force_overwrite"] is False:
            err_msg = "Output files (e.g. " + str(output_cna_segments) + \
                ") already exist - please specify different sample_id or add option --force_overwrite"
            error_message(err_msg, logger)

    # check if input rna fusion variants exist
    if not arg_dict["rna_fusion_tumor"] is None:
        if not os.path.exists(os.path.abspath(arg_dict["rna_fusion_tumor"])):
            err_msg = "Input file (" + \
                str(arg_dict["rna_fusion_tumor"]) + ") does not exist"
            error_message(err_msg, logger)
        if not (os.path.abspath(arg_dict["rna_fusion_tumor"]).endswith(".tsv") or os.path.abspath(arg_dict["rna_fusion_tumor"]).endswith(".txt")):
            err_msg = "RNA fusion variants file (" + os.path.abspath(
                arg_dict["rna_fusion_tumor"]) + ") does not have the correct file extension (.tsv or .txt)"
            error_message(err_msg, logger)
        input_rna_fusion_basename = os.path.basename(
            str(arg_dict["rna_fusion_tumor"]))
        input_rna_fusion_dir = os.path.dirname(
            os.path.abspath(arg_dict["rna_fusion_tumor"]))

    # check if input rna expression exist
    if not arg_dict["rna_exp_tumor"] is None:
        if not os.path.exists(os.path.abspath(arg_dict["rna_exp_tumor"])):
            err_msg = "Input file (" + \
                str(arg_dict["rna_exp_tumor"]) + ") does not exist"
            error_message(err_msg, logger)
        if not (os.path.abspath(arg_dict["rna_exp_tumor"]).endswith(".tsv") or os.path.abspath(arg_dict["rna_exp_tumor"]).endswith(".txt")):
            err_msg = "RNA gene expression file (" + os.path.abspath(
                arg_dict["rna_exp_tumor"]) + ") does not have the correct file extension (.tsv or .txt)"
            error_message(err_msg, logger)
        input_rna_expression_basename = os.path.basename(
            str(arg_dict["rna_exp_tumor"]))
        input_rna_expression_dir = os.path.dirname(
            os.path.abspath(arg_dict["rna_exp_tumor"]))

      # check if input rna fusion variants exist
    if not arg_dict["cpsr_report"] is None:
        if not os.path.exists(os.path.abspath(arg_dict["cpsr_report"])):
            err_msg = "Input file (" + \
                str(arg_dict["cpsr_report"]) + ") does not exist"
            error_message(err_msg, logger)
        if not (os.path.abspath(arg_dict["cpsr_report"]).endswith(".json.gz")):
            err_msg = "CPSR report file (" + os.path.abspath(
                arg_dict["cpsr_report"]) + ") does not have the correct file extension (.json.gz)"
            error_message(err_msg, logger)
        input_cpsr_report_basename = os.path.basename(
            str(arg_dict["cpsr_report"]))
        input_cpsr_report_dir = os.path.dirname(
            os.path.abspath(arg_dict["cpsr_report"]))

    # check the existence of base folder
    base_dir = os.path.abspath(arg_dict["pcgr_dir"])
    if not os.path.isdir(base_dir):
        err_msg = "Base directory (" + str(base_dir) + ") does not exist"
        error_message(err_msg, logger)

    # check the existence of data folder within the base folder
    db_dir = os.path.join(os.path.abspath(arg_dict["pcgr_dir"]), "data")
    if not os.path.isdir(db_dir):
        err_msg = "Data directory (" + str(db_dir) + ") does not exist"
        error_message(err_msg, logger)

    # check the existence of specified assembly data folder within the base folder
    db_assembly_dir = os.path.join(os.path.abspath(
        arg_dict["pcgr_dir"]), "data", arg_dict["genome_assembly"])
    if not os.path.isdir(db_assembly_dir):
        err_msg = "Data directory for the specified genome assembly (" + str(
            db_assembly_dir) + ") does not exist"
        error_message(err_msg, logger)

    # check the existence of RELEASE_NOTES (starting from 0.4.0)
    rel_notes_file = os.path.join(os.path.abspath(
        arg_dict["pcgr_dir"]), "data", arg_dict["genome_assembly"], "RELEASE_NOTES")
    if not os.path.exists(rel_notes_file):
        err_msg = "The PCGR data bundle is outdated - please download the latest data bundle (see github.com/sigven/pcgr for instructions)"
        error_message(err_msg, logger)

    f_rel_not = open(rel_notes_file, "r")
    compliant_data_bundle = 0
    for line in f_rel_not:
        if pcgr_vars.DB_VERSION in line:
            compliant_data_bundle = 1

    f_rel_not.close()

    if compliant_data_bundle == 0:
        err_msg = "The PCGR data bundle is not compliant with the software version - please download the latest software and data bundle (see https://github.com/sigven/pcgr for instructions)"
        error_message(err_msg, logger)

    pcgr_paths = {
      "input_vcf_dir": input_vcf_dir,
      "input_cna_dir": input_cna_dir,
      "input_rna_fusion_dir": input_rna_fusion_dir,
      "input_rna_expression_dir": input_rna_expression_dir,
      "input_cpsr_report_dir": input_cpsr_report_dir,
      "input_cna_plot_dir": input_cna_plot_dir,
      "panel_normal_vcf_dir": panel_normal_vcf_dir,
      "db_dir": db_assembly_dir,
      "base_dir": base_dir,
      "output_dir": output_dir_full,
      "panel_normal_vcf_basename": panel_normal_vcf_basename,
      "input_vcf_basename": input_vcf_basename,
      "input_cna_basename": input_cna_basename,
      "input_rna_fusion_basename": input_rna_fusion_basename,
      "input_rna_expression_basename": input_rna_expression_basename,
      "input_cpsr_report_basename": input_cpsr_report_basename,
      "input_cna_plot_basename": input_cna_plot_basename,
    }

    return pcgr_paths

##---- CPSR ----##

def check_args_cpsr(arg_dict):

    logger = getlogger('cpsr-validate-input-arguments-a')
    arg_dict['vep_regulatory'] = True
    ## Required arguments
    ## Check that query VCF is set and exists
    if arg_dict['input_vcf'] is None or not os.path.exists(arg_dict['input_vcf']):
        err_msg = f"Required argument '--input_vcf' does not exist ({arg_dict['input_vcf']})."
        error_message(err_msg,logger)
    ## Check that PCGR directory (with data bundle) is provided and exists
    if arg_dict['pcgr_dir'] is None or not os.path.exists(arg_dict['pcgr_dir']):
        err_msg = f"Required argument '--pcgr_dir' does not exist ({arg_dict['pcgr_dir']})."
        error_message(err_msg,logger)
    ## Check that genome assembly is set
    if arg_dict['genome_assembly'] is None:
        err_msg = f"Required argument '--genome_assembly' has no/undefined value ({arg_dict['genome_assembly']})."
        error_message(err_msg,logger)
    ## Check that sample identifier is set and is of appropriate length (minimum two characters)
    if arg_dict['sample_id'] is None:
        err_msg = f"Required argument '--sample_id' has no/undefined value ({arg_dict['sample_id']})."
        error_message(err_msg,logger)

    if len(arg_dict['sample_id']) <= 2:
        err_msg = f"Sample name identifier ('--sample_id') requires a name with more than two characters ({arg_dict['sample_id']})."
        error_message(err_msg,logger)

    ### Optional arguments
    ## Provide virtual_panel_id or a custom list from panel 0
    if arg_dict['virtual_panel_id'] == "-1" and not arg_dict['custom_list']:
        err_msg = 'Provide valid virtual panel identifier(s) through --panel_id (0 - 42) or provide custom list of panel 0 genes (single column text file) through --custom_list'
        error_message(err_msg,logger)
    if arg_dict['custom_list'] and arg_dict['virtual_panel_id'] != "-1":
        err_msg =  "Option --panel_id cannot be used in conjunction with --custom_list"
        error_message(err_msg, logger)
    if arg_dict['maf_upper_threshold'] <= 0 or arg_dict['maf_upper_threshold'] > 1:
        err_msg = f"MAF upper threshold must be greater than 0 and below 1, current value is {arg_dict['maf_upper_threshold']}"
        error_message(err_msg,logger)
    if arg_dict['vcfanno_n_proc'] <= 0 or arg_dict['vcfanno_n_proc'] > 15:
        err_msg = f"Number of processes that vcfanno can use during annotation must be above 0 and not more than 15, current value is {arg_dict['vcfanno_n_proc']}."
        error_message(err_msg,logger)

    ## Check that panel identifier(s) are set appropriately
    if arg_dict['virtual_panel_id'] != "-1" and not arg_dict['custom_list']:
        if not ',' in arg_dict['virtual_panel_id']:
            if str(arg_dict['virtual_panel_id']).isdigit():
                panel_id = int(arg_dict['virtual_panel_id'])
                if not (panel_id >= 0 and panel_id <= 42):
                    err_msg =  'A single panel chosen with \'--panel_id\' must be in the range 0 - 42'
                    error_message(err_msg, logger)
            else:
                err_msg =  'A single panel chosen with \'--panel_id\' must be a proper integer - not \'' + str(arg_dict['virtual_panel_id']) + '\''
                error_message(err_msg, logger)
        else:
            panels = str(arg_dict['virtual_panel_id']).split(',')
            for p in panels:
                #p = int(p)
                if str(p).isdigit():
                    panel_id = int(p)
                    if panel_id < 1 or panel_id > 42:
                        err_msg =  'Multiple panels submitted as comma-separated string with \'--panel_id\' must take values in the range 1 - 42'
                        error_message(err_msg, logger)
                else:
                    err_msg =  f"Multiple panels submitted as comma-separated string with '--panel_id' must contain proper integer values only - \'{arg_dict['virtual_panel_id']}\' contains non-integer entries."
                    error_message(err_msg, logger)


    if (arg_dict['custom_list'] or arg_dict['virtual_panel_id'] == "0" ) and arg_dict['diagnostic_grade_only']:
        warn_msg = 'Option \'--diagnostic_grade_only\' applies ONLY to panel identifiers from Genomics England PanelApp - will be ignored'
        warn_message(warn_msg, logger)

    ## VEP options
    if arg_dict['vep_n_forks'] <= 0 or arg_dict['vep_n_forks'] > 4:
        err_msg = f"Number of forks that VEP can use during annotation must be above 0 and not more than 4, current value is {arg_dict['vep_n_forks']}"
        error_message(err_msg,logger)

    if arg_dict['vep_buffer_size'] <= 0 or arg_dict['vep_buffer_size'] > 30000:
        err_msg = "Internal VEP buffer size, corresponding to the number of variants that are read in to memory simultaneously, must be above 0 and not more than 30,000, current value is {arg_dict['vep_buffer_size']}"
        error_message(err_msg,logger)

    ## Check that VEP pick criteria is formatted correctly
    if not arg_dict['vep_pick_order'] is None:
        values = str(arg_dict['vep_pick_order']).split(',')
        permitted_sources = ['canonical','appris','tsl','biotype','ccds','rank','length','mane']
        num_permitted_sources = 0
        for v in values:
            if v in permitted_sources:
                num_permitted_sources += 1

        if num_permitted_sources != 8:
            err_msg = "Option 'vep_pick_order' = " + str(arg_dict['vep_pick_order']) + " is formatted incorrectly, should be " + \
               "a comma-separated string of the following values: canonical,appris,tsl,biotype,ccds,rank,length,mane"
            error_message(err_msg, logger)
    return


def verify_input_files_cpsr(arg_dict):

    logger = getlogger('cpsr-validate-input-arguments-b')
    input_vcf_dir = "NA"
    db_dir = "NA"
    base_dir = "NA"
    output_dir_full = "NA"
    input_vcf_basename = "NA"
    input_customlist_basename = "NA"
    input_customlist_dir = "NA"

    # create output folder (if not already exists)
    output_dir_full = safe_makedir(os.path.abspath(arg_dict['output_dir']))

    ## check if input BED exist
    if not arg_dict['custom_list'] is None:
        if not os.path.exists(os.path.abspath(arg_dict['custom_list'])):
            err_msg = f"Input file ({arg_dict['custom_list']}) does not exist"
            error_message(err_msg,logger)

        input_customlist_basename = os.path.basename(str(arg_dict['custom_list']))
        input_customlist_dir = os.path.dirname(os.path.abspath(arg_dict['custom_list']))

    ## check if input vcf exist
    if not arg_dict['input_vcf'] is None:
        if not os.path.exists(os.path.abspath(arg_dict['input_vcf'])):
            err_msg = f"Input file ({arg_dict['input_vcf']}) does not exist"
            error_message(err_msg,logger)

        if not (os.path.abspath(arg_dict['input_vcf']).endswith('.vcf') or os.path.abspath(arg_dict['input_vcf']).endswith('.vcf.gz')):
            err_msg = f"VCF input file ({os.path.abspath(arg_dict['input_vcf'])}) does not have the correct file extension (.vcf or .vcf.gz)"
            error_message(err_msg,logger)

        ## check that tabix file exist if bgzipped files is given
        if os.path.abspath(arg_dict['input_vcf']).endswith('.vcf.gz'):
            tabix_file = arg_dict['input_vcf'] + '.tbi'
            if not os.path.exists(os.path.abspath(tabix_file)):
                err_msg = "Tabix file (i.e. '.gz.tbi') is not present for the bgzipped VCF input file (" + os.path.abspath(arg_dict['input_vcf']) + "). Please make sure your input VCF is properly compressed and indexed (bgzip + tabix)"
                error_message(err_msg,logger)

        input_vcf_basename = os.path.basename(str(arg_dict['input_vcf']))
        input_vcf_dir = os.path.dirname(os.path.abspath(arg_dict['input_vcf']))

        ## if output vcf exist and overwrite not set
        output_vcf = os.path.join(str(output_dir_full), str(arg_dict['sample_id'])) + '.cpsr.' + str(arg_dict['genome_assembly']) + '.vcf.gz'
        if os.path.exists(output_vcf) and arg_dict['force_overwrite'] is False:
            err_msg = f"Output files (e.g. {output_vcf}) already exist - please specify different sample_id or add option --force_overwrite"
            error_message(err_msg,logger)

    ## check the existence of base folder
    base_dir = os.path.abspath(arg_dict['pcgr_dir'])
    if not os.path.isdir(base_dir):
        err_msg = f"Base directory ({base_dir}) does not exist"
        error_message(err_msg,logger)

    ## check the existence of data folder within the base folder
    db_dir = os.path.join(os.path.abspath(arg_dict['pcgr_dir']), 'data')
    if not os.path.isdir(db_dir):
        err_msg = f"Data directory ({db_dir}) does not exist"
        error_message(err_msg,logger)

    ## check the existence of specified assembly data folder within the base folder
    db_assembly_dir = os.path.join(os.path.abspath(arg_dict['pcgr_dir']), 'data', arg_dict['genome_assembly'])
    if not os.path.isdir(db_assembly_dir):
        err_msg = f"Data directory for the specified genome assembly ({db_assembly_dir}) does not exist"
        error_message(err_msg,logger)

    ## check the existence of RELEASE_NOTES
    rel_notes_file = os.path.join(os.path.abspath(arg_dict['pcgr_dir']), 'data', arg_dict['genome_assembly'],'RELEASE_NOTES')
    if not os.path.exists(rel_notes_file):
        err_msg = 'The PCGR data bundle is outdated - please download the latest data bundle (see github.com/sigven/cpsr for instructions)'
        error_message(err_msg,logger)

    f_rel_not = open(rel_notes_file,'r')
    compliant_data_bundle = 0
    for line in f_rel_not:
        if pcgr_vars.DB_VERSION in line:
            compliant_data_bundle = 1

    f_rel_not.close()

    if compliant_data_bundle == 0:
        err_msg = 'The PCGR data bundle is not compliant with the software version - please download the latest software and data bundle (see https://github.com/sigven/cpsr for instructions)'
        error_message(err_msg,logger)

    cpsr_paths = {
            "input_vcf_dir": input_vcf_dir,
            "input_customlist_dir": input_customlist_dir,
            "db_dir": db_assembly_dir,
            "base_dir": base_dir,
            "output_dir": output_dir_full,
            "input_vcf_basename": input_vcf_basename,
            "input_customlist_basename": input_customlist_basename,
            }

    return cpsr_paths
