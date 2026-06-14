#!/usr/bin/env python

import os
import csv
import pybedtools
import pandas as pd
import numpy as np
import logging
import re

from typing import Optional, cast
from pybedtools import BedTool
from pcgr import pcgr_vars
from pcgr.annoutils import nuclear_chromosomes
from pcgr.utils import error_message, warn_message, check_file_exists, remove_file, pd_to_csv
from pcgr.biomarker import load_all_biomarkers
from pcgr.expression import integrate_variant_expression


def estimate_tumor_ploidy(cna_df: pd.DataFrame,
                           focal_threshold_mb: float = 3.0,
                           high_level_cn_threshold: int = 10,
                           logger: Optional[logging.Logger] = None) -> float:
    """
    Auto-estimate tumor ploidy from CNA segments using weighted median approach.

    Excludes focal high-level amplifications from the calculation, as these represent
    localized oncogenic events rather than genome-wide ploidy state.

    Args:
        cna_df: DataFrame with nMajor, nMinor, Start, End columns
        focal_threshold_mb: Segments smaller than this (in Mb) are considered focal. Default: 3.0 Mb
        high_level_cn_threshold: Copy number above this is considered high-level amplification. Default: 10
        logger: Logger instance

    Returns:
        Estimated ploidy value (rounded to 2 decimal places)
    """
    if logger is None:
        logger = logging.getLogger("pcgr-estimate-ploidy")

    # Calculate total copy number and segment length
    cna_df_copy = cna_df.copy()
    cna_df_copy['total_cn'] = cna_df_copy['nMajor'].astype(int) + cna_df_copy['nMinor'].astype(int)
    cna_df_copy['seg_length'] = cna_df_copy['End'].astype(int) - cna_df_copy['Start'].astype(int)
    cna_df_copy['seg_length_mb'] = cna_df_copy['seg_length'] / 1e6

    # Identify focal high-level amplifications
    focal_highlevel_mask = (
        (cna_df_copy['seg_length_mb'] < focal_threshold_mb) &
        (cna_df_copy['total_cn'] > high_level_cn_threshold)
    )

    num_excluded = focal_highlevel_mask.sum()
    total_segments = len(cna_df_copy)

    # Filter out focal high-level amplifications for ploidy estimation
    cna_df_filtered = cna_df_copy[~focal_highlevel_mask].copy()

    if len(cna_df_filtered) == 0:
        logger.warning(
            "All segments are focal high-level amplifications - using all segments for ploidy estimation"
        )
        cna_df_filtered = cna_df_copy
        num_excluded = 0

    # Weighted average by segment length
    # Note: Could also use weighted median for more robustness, but weighted mean
    # provides a good balance between accuracy and computational efficiency
    estimated_ploidy = np.average(
        cna_df_filtered['total_cn'],
        weights=cna_df_filtered['seg_length']
    )

    # Round to 2 decimal places
    estimated_ploidy = round(estimated_ploidy, 2)

    # Log filtering statistics
    if num_excluded > 0:
        logger.info(
            f"Ploidy estimation: excluded {num_excluded}/{total_segments} segments "
            f"(focal < {focal_threshold_mb} Mb AND CN > {high_level_cn_threshold})"
        )

    # Apply bounds
    if estimated_ploidy < 1.0:
        logger.warning(f"Auto-estimated ploidy ({estimated_ploidy}) is < 1.0 - setting to 1.0")
        estimated_ploidy = 1.0

    if estimated_ploidy > 8.0:
        logger.warning(f"Auto-estimated ploidy ({estimated_ploidy}) is > 8.0 - setting to 8.0")
        estimated_ploidy = 8.0

    return estimated_ploidy


def generate_oncokb_fusion_input(fusion_df: pd.DataFrame,
                                  output_fname: str,
                                  logger: Optional[logging.Logger] = None) -> None:
    """
    Generate OncoKB-compatible fusion input file from annotated fusion DataFrame.

    Args:
        fusion_df: Annotated fusion DataFrame with SAMPLE_ID and FUSION_GENE columns
        output_fname: Path to output OncoKB fusion input file
        logger: Logger instance
    """
    if logger is None:
        logger = logging.getLogger("pcgr-oncokb-fusion")

    if fusion_df.empty:
        logger.warning("Empty fusion dataframe - skipping OncoKB input generation")
        return

    # Create OncoKB input format: Tumor_Sample_Barcode, Fusion, VAR_ID
    # VAR_ID is included for cross-reference (FusionAnnotator does not pass it through;
    # the join back to the fusion TSV uses the normalized Fusion gene name)
    has_varid = 'VAR_ID' in fusion_df.columns
    oncokb_df = pd.DataFrame({
        'Tumor_Sample_Barcode': fusion_df['SAMPLE_ID'],
        'Fusion':               fusion_df['FUSION_GENE'],
        **({"VAR_ID": fusion_df['VAR_ID']} if has_varid else {}),
    })

    # Remove any rows with missing fusion gene
    oncokb_df = oncokb_df[oncokb_df['Fusion'].notna()]

    # Replace -- with - (OncoKB API preferred format) and :: with - as well
    oncokb_df['Fusion'] = oncokb_df['Fusion'].str.replace('--', '-', regex=False)
    oncokb_df['Fusion'] = oncokb_df['Fusion'].str.replace('::', '-', regex=False)

    # Deduplicate on Fusion (keep first VAR_ID for each unique gene pair)
    oncokb_df = oncokb_df.drop_duplicates(subset=['Fusion'])

    if len(oncokb_df) > 0:
        # Write as plain text TSV file
        oncokb_df.to_csv(output_fname, sep="\t", index=False)
        logger.info(f"Generated OncoKB fusion input file: {os.path.basename(output_fname)} ({len(oncokb_df)} fusions)")
    else:
        logger.warning("No valid fusions for OncoKB input - skipping file generation")


def generate_oncokb_cna_input(cna_df: pd.DataFrame,
                               output_fname: str,
                               logger: Optional[logging.Logger] = None) -> None:
    """
    Generate OncoKB-compatible CNA input file from annotated CNA DataFrame.

    Args:
        cna_df: Annotated CNA DataFrame with SAMPLE_ID, SYMBOL, and VARIANT_CLASS columns
        output_fname: Path to output OncoKB CNA input file
        logger: Logger instance
    """
    if logger is None:
        logger = logging.getLogger("pcgr-oncokb-cna")

    if cna_df.empty:
        logger.warning("Empty CNA dataframe - skipping OncoKB input generation")
        return

    # Filter for only annotated CNAs (amplification/deletion events)
    cna_oncokb = cna_df[cna_df['VARIANT_CLASS'].isin(pcgr_vars.VARIANT_CLASS_TO_OKB_CNA.keys())].copy()

    if len(cna_oncokb) == 0:
        logger.warning("No amplification/deletion events found for OncoKB input - skipping file generation")
        return

    # Create OncoKB input format: Tumor_Sample_Barcode, Hugo_Symbol, Copy_Number_Alteration, VAR_ID
    # VAR_ID is included for cross-reference (CnaAnnotator does not pass it through;
    # the join back to the CNA TSV uses Hugo_Symbol + Copy_Number_Alteration)
    has_varid = 'VAR_ID' in cna_oncokb.columns
    oncokb_df = pd.DataFrame({
        'Tumor_Sample_Barcode': cna_oncokb['SAMPLE_ID'],
        'Hugo_Symbol': cna_oncokb['SYMBOL'],
        'Copy_Number_Alteration': cna_oncokb['VARIANT_CLASS'].map(pcgr_vars.VARIANT_CLASS_TO_OKB_CNA),
        **({"VAR_ID": cna_oncokb['VAR_ID']} if has_varid else {}),
    })

    # Remove any rows with missing gene symbol
    oncokb_df = oncokb_df[oncokb_df['Hugo_Symbol'].notna()]

    # Deduplicate on the join key (keep first VAR_ID for each unique gene + alteration type)
    oncokb_df = oncokb_df.drop_duplicates(subset=['Hugo_Symbol', 'Copy_Number_Alteration'])

    if len(oncokb_df) > 0:
        # Write as plain text TSV file
        oncokb_df.to_csv(output_fname, sep="\t", index=False)
        logger.info(f"Generated OncoKB CNA input file: {os.path.basename(output_fname)} ({len(oncokb_df)} events)")
    else:
        logger.warning("No valid CNA events for OncoKB input - skipping file generation")


def append_oncokb_fusion_annotations(
        tsv_fname: str,
        oncokb_fusion_output: Optional[str] = None,
        logger: Optional[logging.Logger] = None) -> None:
    """
    Merge OncoKB fusion annotations into the final PCGR fusion TSV.

    Join key: normalized FUSION_GENE (-- and :: replaced with -) matched
    against the Fusion column in the FusionAnnotator.py output.
    Adds NA columns when OncoKB was not run.

    Args:
        tsv_fname:            Path to the fusion TSV (read + overwritten, uncompressed)
        oncokb_fusion_output: Path to FusionAnnotator.py output file
        logger:               Logger instance
    """
    if logger is None:
        logger = logging.getLogger("pcgr-oncokb-fusion-append")

    okb_source_cols  = list(pcgr_vars.ONCOKB_COLS.keys())
    okb_renamed_cols = list(pcgr_vars.ONCOKB_COLS.values())

    if not os.path.exists(tsv_fname):
        logger.warning(f"Fusion TSV not found: {tsv_fname} - skipping OncoKB fusion merge")
        return

    tsv_df = pd.read_csv(tsv_fname, sep="\t", low_memory=False)

    if 'FUSION_GENE' not in tsv_df.columns:
        logger.warning("FUSION_GENE absent from fusion TSV - cannot merge OncoKB fusion annotations")
        for col in okb_renamed_cols:
            tsv_df[col] = pd.NA
        tsv_df.to_csv(tsv_fname, sep="\t", index=False)
        return

    # Drop any pre-existing OKB columns to make this function idempotent
    tsv_df = tsv_df.drop(columns=[c for c in okb_renamed_cols if c in tsv_df.columns])

    # Normalize FUSION_GENE the same way the OncoKB input was prepared
    tsv_df['_OKB_FUSION'] = (tsv_df['FUSION_GENE']
                              .str.replace('--', '-', regex=False)
                              .str.replace('::', '-', regex=False))

    oncokb_df = None
    if oncokb_fusion_output and os.path.exists(oncokb_fusion_output):
        try:
            df = pd.read_csv(oncokb_fusion_output, sep="\t", low_memory=False)
            if 'Fusion' in df.columns:
                present = [c for c in okb_source_cols if c in df.columns]
                if present:
                    keep = ['Fusion'] + present
                    oncokb_df = df[keep].drop_duplicates(subset='Fusion')
            else:
                logger.warning("Fusion column absent from OncoKB fusion output - cannot merge")
        except Exception as e:
            logger.warning(f"Could not read OncoKB fusion output: {e}")

    if oncokb_df is not None and len(oncokb_df) > 0:
        for col in okb_source_cols:
            if col not in oncokb_df.columns:
                oncokb_df[col] = pd.NA
        oncokb_df = oncokb_df.rename(columns={**pcgr_vars.ONCOKB_COLS, 'Fusion': '_OKB_FUSION'})
        tsv_df = tsv_df.merge(oncokb_df[['_OKB_FUSION'] + okb_renamed_cols],
                              on='_OKB_FUSION', how='left')
        n_hits = tsv_df[okb_renamed_cols[0]].notna().sum()
        logger.info(f"OncoKB fusion annotations merged: {n_hits} / {len(tsv_df)} fusions annotated")
    else:
        logger.info("No OncoKB fusion annotations available - adding NA columns")
        for col in okb_renamed_cols:
            tsv_df[col] = pd.NA

    tsv_df = tsv_df.drop(columns=['_OKB_FUSION'])
    tsv_df.to_csv(tsv_fname, sep="\t", index=False)
    logger.info(f"OncoKB fusion columns written to {os.path.basename(tsv_fname)}")


def append_oncokb_cna_annotations(
        tsv_fname: str,
        oncokb_cna_output: Optional[str] = None,
        logger: Optional[logging.Logger] = None) -> None:
    """
    Merge OncoKB CNA annotations into the PCGR gene-level CNA TSV.

    Join key: SYMBOL + Copy_Number_Alteration (mapped from VARIANT_CLASS).
    Adds NA columns when OncoKB was not run.

    Args:
        tsv_fname:         Path to the gene-level CNA TSV (read + overwritten, uncompressed)
        oncokb_cna_output: Path to CnaAnnotator.py output file
        logger:            Logger instance
    """
    if logger is None:
        logger = logging.getLogger("pcgr-oncokb-cna-append")

    okb_source_cols  = list(pcgr_vars.ONCOKB_COLS.keys())
    okb_renamed_cols = list(pcgr_vars.ONCOKB_COLS.values())

    if not os.path.exists(tsv_fname):
        logger.warning(f"CNA gene TSV not found: {tsv_fname} - skipping OncoKB CNA merge")
        return

    tsv_df = pd.read_csv(tsv_fname, sep="\t", low_memory=False)

    if 'SYMBOL' not in tsv_df.columns or 'VARIANT_CLASS' not in tsv_df.columns:
        logger.warning("SYMBOL or VARIANT_CLASS absent from CNA gene TSV - cannot merge OncoKB CNA annotations")
        for col in okb_renamed_cols:
            tsv_df[col] = pd.NA
        tsv_df.to_csv(tsv_fname, sep="\t", index=False)
        return

    # Drop any pre-existing OKB columns to make this function idempotent
    tsv_df = tsv_df.drop(columns=[c for c in okb_renamed_cols if c in tsv_df.columns])

    # Map VARIANT_CLASS to OncoKB Copy_Number_Alteration for the join
    tsv_df['_OKB_CNA_TYPE'] = tsv_df['VARIANT_CLASS'].map(pcgr_vars.VARIANT_CLASS_TO_OKB_CNA)

    oncokb_df = None
    if oncokb_cna_output and os.path.exists(oncokb_cna_output):
        try:
            df = pd.read_csv(oncokb_cna_output, sep="\t", low_memory=False)
            required = {'Hugo_Symbol', 'Copy_Number_Alteration'}
            if required.issubset(df.columns):
                present = [c for c in okb_source_cols if c in df.columns]
                if present:
                    keep = ['Hugo_Symbol', 'Copy_Number_Alteration'] + present
                    oncokb_df = df[keep].drop_duplicates(
                        subset=['Hugo_Symbol', 'Copy_Number_Alteration'])
            else:
                logger.warning("Hugo_Symbol or Copy_Number_Alteration absent from OncoKB CNA output")
        except Exception as e:
            logger.warning(f"Could not read OncoKB CNA output: {e}")

    if oncokb_df is not None and len(oncokb_df) > 0:
        for col in okb_source_cols:
            if col not in oncokb_df.columns:
                oncokb_df[col] = pd.NA
        oncokb_df = oncokb_df.rename(columns=pcgr_vars.ONCOKB_COLS)
        tsv_df = tsv_df.merge(
            oncokb_df[['Hugo_Symbol', 'Copy_Number_Alteration'] + okb_renamed_cols],
            left_on  = ['SYMBOL', '_OKB_CNA_TYPE'],
            right_on = ['Hugo_Symbol', 'Copy_Number_Alteration'],
            how      = 'left'
        ).drop(columns=['Hugo_Symbol', 'Copy_Number_Alteration', '_OKB_CNA_TYPE'])
        n_hits = tsv_df[okb_renamed_cols[0]].notna().sum()
        logger.info(f"OncoKB CNA annotations merged: {n_hits} / {len(tsv_df)} gene-level events annotated")
    else:
        logger.info("No OncoKB CNA annotations available - adding NA columns")
        tsv_df = tsv_df.drop(columns=['_OKB_CNA_TYPE'])
        for col in okb_renamed_cols:
            tsv_df[col] = pd.NA

    tsv_df.to_csv(tsv_fname, sep="\t", index=False)
    logger.info(f"OncoKB CNA columns written to {os.path.basename(tsv_fname)}")


def annotate_fusions(input_fusion_fname: str,
                     output_fusion_fname: str,
                     oncokb_input_fname: str,
                     sample_id: str,
                     build: str,
                     refdata_assembly_dir: str,
                     logger: Optional[logging.Logger] = None) -> int:
    """
    Annotate gene fusions in a given fusion file.
    Args:
        input_fusion_file (str): Path to the user-provided fusion file.
        output_fusion_fname (str): File name of the annotated output file.
        sample_id: Sample identifier
        build: Genome build version
        refdata_assembly_dir (str): Path to the build-specific PCGR database directory.
        logger (logging.Logger, optional): Logger. Defaults to None.
    Returns:
        int: 0 if successful.
    """
    if logger is None:
        logger = logging.getLogger("pcgr-annotate-rna-fusion")

    check_file_exists(input_fusion_fname, logger = logger)
    fusion_df = pd.read_csv(input_fusion_fname, sep="\t", na_values=".")
    if 'FusionGene' in fusion_df.columns:
        fusion_df['aberration_key'] = (
            fusion_df['FusionGene'].fillna('').str.replace('--', '::', regex=False).apply(
                lambda x: x + "_transcript_fusion" if x else None)
            )

    ## load fusion biomarker evidence
    biomarkers = load_all_biomarkers(
        logger, refdata_assembly_dir, biomarker_vartype = 'FUSION', biomarker_variant_origin = 'Somatic')
        
    if fusion_df.empty is True:
        warn_msg = 'Fusion file is empty - contains NO fusion events'
        warn_message(warn_msg, logger)
        return -1
    
    if 'actionable_df' in biomarkers:
        if 'aberration_key' not in fusion_df.columns:
            warn_msg = "Could not find required column 'aberration_key' in input fusion file - exiting"
            warn_message(warn_msg, logger)
        if fusion_df['aberration_key'].isnull().all() is True:
            warn_msg = "Could not find any valid fusion gene annotations in input fusion file - omitting gene-level annotations"
            warn_message(warn_msg, logger)
            return -1
        fusion_df = fusion_df.merge(
            biomarkers['actionable_df'], left_on=["aberration_key"], right_on=["aberration_key"], how="left")
        fusion_df.loc[fusion_df['biomarker_match'].isnull(),"biomarker_match"] = '.'
    else:
        warn_msg = "Could not find any fusion biomarker evidence in PCGR data bundle - omitting gene-level annotations"
        warn_message(warn_msg, logger)
        fusion_df['biomarker_match'] = '.'
    
    
    fusion_df = fusion_df.drop(columns=['aberration_key'])
    fusion_df = fusion_df.rename(columns={
        "FusionGene": "FUSION_GENE",
        "LeftBreakpoint": "BREAKPOINT_5P",
        "RightBreakpoint": "BREAKPOINT_3P",
        "SplitReads": "SPLIT_READS",
        "biomarker_match": "BIOMARKER_MATCH"
    })
    # Normalise breakpoint chromosome prefix: strip any leading 'chr' so that
    # downstream R code and display URLs always work with bare chromosome names
    # (e.g. "2:42264731" not "chr2:42264731"). The R display layer adds 'chr'
    # back when building UCSC browser links.
    for bp_col in ("BREAKPOINT_5P", "BREAKPOINT_3P"):
        if bp_col in fusion_df.columns:
            fusion_df[bp_col] = fusion_df[bp_col].str.replace(
                r'^chr', '', regex=True)
    if("Score" in fusion_df.columns):
        fusion_df = fusion_df.rename(columns={"Score": "FUSION_SCORE"})
    fusion_df['SAMPLE_ID'] = sample_id
    fusion_df['VARIANT_CLASS'] = 'fusion'
    fusion_df["VAR_ID"] = "FUSION_" + (fusion_df.index + 1).astype(str)

    desired_order = [
        "SAMPLE_ID",
        "VAR_ID",
        "FUSION_GENE",
        "VARIANT_CLASS",
        "BREAKPOINT_5P",
        "BREAKPOINT_3P",
        "SPLIT_READS",
        "FUSION_SCORE",       
        "BIOMARKER_MATCH"
    ]
    fusion_df = fusion_df[[c for c in desired_order if c in fusion_df.columns]]

    pd_to_csv(
        df = fusion_df,
        fname = output_fusion_fname,
        col_sep = "\t")

    # Generate OncoKB input file
    generate_oncokb_fusion_input(fusion_df, oncokb_input_fname, logger)

    return 0

def _twohit_vaf_flag(
    vaf_tumor,
    n_major: int,
    n_minor: int,
    loh_type: str,
    tumor_purity,
    tolerance: float = 0.15) -> str:
    """
    Assess whether a somatic variant's observed VAF is consistent with it
    sitting on the retained (major) allele of a LOH segment.

    Expected VAF formula (clonal variant on retained allele):
        VAF_expected = (p * n_mut) / (p * n_total + 2 * (1 - p))

    For deletion LOH (nMajor=1, nMinor=0): n_mut=1, n_total=1
        → VAF_expected = p / (2 - p)
    For copy-neutral LOH (nMajor=2, nMinor=0): n_mut ambiguous (1 or 2).
        Lower bound uses n_mut=1 → VAF_expected = p / 2.
    For homdel (nMajor=0, nMinor=0): no retained allele → VAF_UNKNOWN.

    Returns one of: 'VAF_CONSISTENT', 'VAF_LOW', 'VAF_UNKNOWN'
    """
    try:
        purity = float(tumor_purity)
    except (TypeError, ValueError):
        return 'VAF_UNKNOWN'
    if pd.isna(purity) or purity <= 0 or purity > 1:
        return 'VAF_UNKNOWN'

    try:
        vaf = float(vaf_tumor)
    except (TypeError, ValueError):
        return 'VAF_UNKNOWN'

    n_total = int(n_major) + int(n_minor)
    if n_total == 0:
        return 'VAF_UNKNOWN'

    denominator = purity * n_total + 2.0 * (1.0 - purity)
    if denominator <= 0:
        return 'VAF_UNKNOWN'

    if loh_type == 'copy-neutral':
        # Use lower bound: n_mut=1 (single copy of mutant on major allele)
        vaf_expected_low = purity / denominator
        consistent = vaf >= (vaf_expected_low - tolerance)
    else:
        # deletion LOH or other: n_mut=1
        vaf_expected = purity / denominator
        consistent = vaf >= (vaf_expected - tolerance)

    return 'VAF_CONSISTENT' if consistent else 'VAF_LOW'


def append_twohit_candidates(cna_df: pd.DataFrame,
                              snv_df_somatic: Optional[pd.DataFrame] = None,
                              snv_df_germline: Optional[pd.DataFrame] = None,
                              logger: Optional[logging.Logger] = None) -> pd.DataFrame:
    """
    For each gene-level CNA record annotated as TSG (TSG == TRUE) with allele-specific LOH
    (deletion or copy-neutral: nMinor=0, nMajor>0), find:
    1. Overlapping somatic LOF SNV/InDels (LOSS_OF_FUNCTION == True) recorded in
       TWOHIT_CANDIDATE_SOMATIC.
    2. Germline pathogenic/likely-pathogenic LOF variants matching by SYMBOL, recorded
       in TWOHIT_CANDIDATE_GERMLINE.

    Somatic matching is by CHROM + SYMBOL + genomic position (SNV POS within the CNA
    segment interval). Segment coordinates are 0-based BED-style, VCF POS is 1-based,
    so the overlap condition is: POS > seg_start AND POS <= seg_end.

    Germline matching is by SYMBOL only, filtered to CLASSIFICATION in
    ['Pathogenic', 'Likely_Pathogenic'] and LOSS_OF_FUNCTION == True.

    Args:
        cna_df: Annotated CNA dataframe (uppercase columns: CHROM, SEGMENT_START, SEGMENT_END,
                SYMBOL, TSG, LOH)
        snv_df_somatic: Somatic SNV/InDel dataframe with CHROM, POS, SYMBOL, LOSS_OF_FUNCTION,
                VAR_ID, CONSEQUENCE, MUTATION_EFFECT_OKB columns
        snv_df_germline: Germline SNV/InDel dataframe (from CPSR) with SYMBOL, CLASSIFICATION,
                LOSS_OF_FUNCTION, VAR_ID, CONSEQUENCE
        logger: Logger instance

    Returns:
        cna_df with TWOHIT_CANDIDATE_SOMATIC and TWOHIT_CANDIDATE_GERMLINE columns added.
        Somatic: comma-separated VAR_ID;CONSEQUENCE;VAF_FLAG
          where VAF_FLAG is one of VAF_CONSISTENT, VAF_LOW, VAF_UNKNOWN
        Germline: comma-separated VAR_ID;CONSEQUENCE
    """
    if logger is None:
        logger = logging.getLogger("pcgr-twohit-candidates")

    cna_df['TWOHIT_CANDIDATE_SOMATIC'] = '.'
    cna_df['TWOHIT_CANDIDATE_GERMLINE'] = '.'

    required_cna_cols = {'CHROM', 'SEGMENT_START', 'SEGMENT_END', 'SYMBOL', 'TSG', 'LOH'}
    missing_cna = required_cna_cols - set(cna_df.columns)
    if missing_cna:
        logger.warning(f"CNA dataframe missing columns for two-hit candidate annotation: {missing_cna} - skipping")
        return cna_df

    #logger.debug(f"Two-hit candidates - CNA TSG unique values: {sorted(cna_df['TSG'].dropna().astype(str).unique().tolist())}")
    #logger.debug(f"Two-hit candidates - CNA LOH unique values: {sorted(cna_df['LOH'].dropna().astype(str).unique().tolist())}")

    _truthy = {True, 'TRUE', 'true', 1, '1'}
    ## TSG with allele-specific LOH (deletion or copy-neutral subtypes: nMinor=0, nMajor>0)
    tsg_loh_mask = cna_df['TSG'].isin(_truthy) & (cna_df['LOH'].isin({'deletion', 'copy-neutral'}))
    tsg_loh_cna = cna_df[tsg_loh_mask]

    if tsg_loh_cna.empty:
        return cna_df

    # --- Somatic two-hit candidates ---
    if snv_df_somatic is not None and isinstance(snv_df_somatic, pd.DataFrame) and not snv_df_somatic.empty:
        required_som_cols = {'CHROM', 'POS', 'SYMBOL', 
                             'LOSS_OF_FUNCTION', 
                             'MUTATION_EFFECT_OKB',
                             'VAR_ID', 'CONSEQUENCE'}
        missing_som = required_som_cols - set(snv_df_somatic.columns)
        if missing_som:
            logger.warning(f"Somatic SNV dataframe missing columns for two-hit candidate annotation: {missing_som} - skipping somatic")
        else:
            lof_mask = snv_df_somatic['LOSS_OF_FUNCTION'].isin(_truthy)
            if 'MUTATION_EFFECT_OKB' in snv_df_somatic.columns:
                okb_lof_mask = snv_df_somatic['MUTATION_EFFECT_OKB'].isin(
                    {'Loss-of-function', 'Likely Loss-of-function'})
                lof_mask = lof_mask | okb_lof_mask
            lof_somatic = snv_df_somatic[lof_mask].copy()
            if lof_somatic.empty:
                logger.info("Two-hit candidates (somatic): no LOF somatic SNV/InDels found")
            else:
                lof_somatic['_var_base'] = (
                    lof_somatic['VAR_ID'].fillna('.').astype(str) + ';' +
                    lof_somatic['CONSEQUENCE'].fillna('.').astype(str)
                )
                lof_somatic['_chrom'] = lof_somatic['CHROM'].astype(str)
                lof_somatic['_symbol'] = lof_somatic['SYMBOL'].astype(str)
                lof_somatic['_pos'] = pd.to_numeric(lof_somatic['POS'], errors='coerce')
                _has_vaf = 'VAF_TUMOR' in lof_somatic.columns

                n_matched_somatic = 0
                for idx, cna_row in tsg_loh_cna.iterrows():
                    chrom = str(cna_row['CHROM'])
                    symbol = str(cna_row['SYMBOL'])
                    seg_start = int(cna_row['SEGMENT_START'])  # 0-based BED
                    seg_end = int(cna_row['SEGMENT_END'])       # 0-based BED, exclusive end
                    loh_type = str(cna_row.get('LOH', '.'))
                    n_major = cna_row.get('CN_MAJOR', cna_row.get('N_MAJOR', 0))
                    n_minor = cna_row.get('CN_MINOR', cna_row.get('N_MINOR', 0))
                    tumor_purity = cna_row.get('TUMOR_PURITY', None)

                    # VCF POS (1-based) overlaps BED interval [seg_start, seg_end) when:
                    # POS > seg_start AND POS <= seg_end
                    overlapping = lof_somatic[
                        (lof_somatic['_chrom'] == chrom) &
                        (lof_somatic['_symbol'] == symbol) &
                        (lof_somatic['_pos'] > seg_start) &
                        (lof_somatic['_pos'] <= seg_end)
                    ].copy()

                    if not overlapping.empty:
                        if _has_vaf:
                            overlapping['_vaf_flag'] = overlapping['VAF_TUMOR'].apply(
                                lambda v: _twohit_vaf_flag(v, n_major, n_minor, loh_type, tumor_purity)
                            )
                        else:
                            overlapping['_vaf_flag'] = 'VAF_UNKNOWN'
                        overlapping['var_string'] = overlapping['_var_base'] + ';' + overlapping['_vaf_flag']
                        cna_df.at[idx, 'TWOHIT_CANDIDATE_SOMATIC'] = ','.join(overlapping['var_string'].tolist())
                        n_matched_somatic += 1

                if n_matched_somatic > 0:
                    logger.info(f"Two-hit candidates (somatic): {n_matched_somatic} TSG LOH gene records annotated with overlapping somatic LoF SNV/InDels")
                else:
                    logger.info("Two-hit candidates (somatic): no overlapping somatic LoF SNV/InDels found for TSG LOH records")

    # --- Germline two-hit candidates ---
    if snv_df_germline is not None and isinstance(snv_df_germline, pd.DataFrame) and not snv_df_germline.empty:
        required_germ_cols = {'SYMBOL', 'CLASSIFICATION', 
                              'LOSS_OF_FUNCTION', 'VAR_ID', 
                              'CONSEQUENCE'}
        missing_germ = required_germ_cols - set(snv_df_germline.columns)
        if missing_germ:
            logger.warning(f"Germline SNV dataframe missing columns for two-hit candidate annotation: {missing_germ} - skipping germline")
        else:
            #logger.debug(f"Two-hit candidates - germline LOSS_OF_FUNCTION unique values: {sorted(snv_df_germline['LOSS_OF_FUNCTION'].dropna().astype(str).unique().tolist())}")
            het_mask = (
                snv_df_germline['GENOTYPE'].str.lower() == 'het'
                if 'GENOTYPE' in snv_df_germline.columns
                else pd.Series(True, index=snv_df_germline.index)
            )
            pathogenic_lof_germline = snv_df_germline[
                snv_df_germline['CLASSIFICATION'].isin(['Pathogenic', 'Likely_Pathogenic']) &
                snv_df_germline['LOSS_OF_FUNCTION'].isin(_truthy) &
                het_mask
            ].copy()

            if pathogenic_lof_germline.empty:
                logger.info("Two-hit candidates (germline): no pathogenic/likely-pathogenic LoF germline variants found")
            else:
                pathogenic_lof_germline['var_string'] = (
                    pathogenic_lof_germline['VAR_ID'].fillna('.').astype(str) + ';' +
                    pathogenic_lof_germline['CONSEQUENCE'].fillna('.').astype(str)
                )
                pathogenic_lof_germline['_symbol'] = pathogenic_lof_germline['SYMBOL'].astype(str)

                n_matched_germline = 0
                for idx, cna_row in tsg_loh_cna.iterrows():
                    symbol = str(cna_row['SYMBOL'])
                    matching = pathogenic_lof_germline[pathogenic_lof_germline['_symbol'] == symbol]

                    if not matching.empty:
                        cna_df.at[idx, 'TWOHIT_CANDIDATE_GERMLINE'] = ','.join(matching['var_string'].tolist())
                        n_matched_germline += 1

                if n_matched_germline > 0:
                    logger.info(f"Two-hit candidates (germline): {n_matched_germline} TSG LOH gene records annotated with pathogenic/LP germline LOF variants")
                else:
                    logger.info("Two-hit candidates (germline): no pathogenic/likely-pathogenic germline LOF variants found for TSG LOH records")

    return cna_df


def _annotate_amplifications(
    df: pd.DataFrame,
    tumor_ploidy: float,
    amp_threshold_absolute: int,
    amp_threshold_relative: float,
    threshold_mode: str,
    logger: logging.Logger,
) -> tuple[pd.DataFrame, float]:
    """
    Mark amplification segments and compute the effective amplification threshold.
    Adds columns: amp_cond, total_cn, fold_change, aberration_key (for amplifications).
    Returns (df, effective_amp_threshold).
    NOTE: total_cn is retained for use by _annotate_gains; caller must drop it afterwards.
    """
    relative_threshold_value = round(tumor_ploidy * amp_threshold_relative, 1)
    df['aberration_key'] = 'nan'
    df['amp_cond'] = False
    df['total_cn'] = df['n_major'] + df['n_minor']
    df['fold_change'] = (df['total_cn'] / tumor_ploidy).round(4)

    if threshold_mode == "absolute":
        logger.info(f"Using absolute threshold mode: amplification = CN >= {amp_threshold_absolute}")
        df.loc[df['total_cn'] >= amp_threshold_absolute, "amp_cond"] = True
    elif threshold_mode == "relative":
        logger.info(
            f"Using relative threshold mode: amplification = CN >= {tumor_ploidy} × {amp_threshold_relative} = {relative_threshold_value}"
        )
        df.loc[df['total_cn'] >= relative_threshold_value, "amp_cond"] = True
    elif threshold_mode == "combined":
        logger.info("Using combined threshold mode (requires BOTH criteria):")
        logger.info(f"Absolute criterion: CN >= {amp_threshold_absolute}")
        logger.info(f"Relative criterion: CN >= {tumor_ploidy} × {amp_threshold_relative} = {relative_threshold_value}")
        df.loc[
            (df['total_cn'] >= amp_threshold_absolute) &
            (df['total_cn'] >= relative_threshold_value), "amp_cond"] = True

    if threshold_mode == "absolute":
        effective_amp_threshold = float(amp_threshold_absolute)
    elif threshold_mode == "relative":
        effective_amp_threshold = float(relative_threshold_value)
    else:
        effective_amp_threshold = float(max(amp_threshold_absolute, relative_threshold_value))
    logger.info(f"Effective amplification threshold: CN >= {effective_amp_threshold}")

    df.loc[df['amp_cond'], 'aberration_key'] = \
        df.loc[df['amp_cond'], 'entrezgene'].astype(str) + '_amplification'

    return df, effective_amp_threshold


def _annotate_gains(
    df: pd.DataFrame,
    tumor_ploidy: float,
    gain_threshold_absolute: int,
    gain_threshold_relative: float,
    threshold_mode: str,
    logger: logging.Logger,
) -> pd.DataFrame:
    """
    Mark gain segments (above gain threshold but below amplification threshold).
    Requires amp_cond and total_cn columns from _annotate_amplifications.
    Adds: gain_cond, aberration_key (for gains).
    Drops: total_cn.
    """
    relative_gain_threshold_value = round(tumor_ploidy * gain_threshold_relative, 1)
    df['gain_cond'] = False

    if threshold_mode == "absolute":
        logger.info(f"Using absolute threshold mode: gain = CN >= {gain_threshold_absolute} (and < amplification threshold)")
        df.loc[df['total_cn'] >= gain_threshold_absolute, "gain_cond"] = True
    elif threshold_mode == "relative":
        logger.info(
            f"Using relative threshold mode: gain = CN >= {tumor_ploidy} × {gain_threshold_relative} = {relative_gain_threshold_value} (and < amplification threshold)"
        )
        df.loc[df['total_cn'] >= relative_gain_threshold_value, "gain_cond"] = True
    elif threshold_mode == "combined":
        logger.info("Using combined threshold mode for gains (requires BOTH criteria, and < amplification threshold):")
        logger.info(f"Absolute criterion: CN >= {gain_threshold_absolute}")
        logger.info(f"Relative criterion: CN >= {tumor_ploidy} × {gain_threshold_relative} = {relative_gain_threshold_value}")
        df.loc[
            (df['total_cn'] >= gain_threshold_absolute) &
            (df['total_cn'] >= relative_gain_threshold_value), "gain_cond"] = True

    # Gains and amplifications are mutually exclusive
    df['gain_cond'] = df['gain_cond'] & ~df['amp_cond']
    df = df.drop(columns=['total_cn'])

    df.loc[df['gain_cond'], 'aberration_key'] = \
        df.loc[df['gain_cond'], 'entrezgene'].astype(str) + '_gain'

    return df


def _annotate_deletions(
    df: pd.DataFrame,
    tumor_ploidy: float,
    del_threshold_absolute: int,
    del_threshold_relative: float,
    threshold_mode: str,
    sex: str,
    logger: logging.Logger,
) -> tuple[pd.DataFrame, int]:
    """
    Mark homozygous, heterozygous, and hemizygous deletion segments.
    Adds: homloss_cond, hetloss_cond, hemloss_cond, aberration_key (for ablations).
    Returns (df, baseline_cn).
    """
    baseline_cn = round(tumor_ploidy)
    relative_del_threshold_value = round(tumor_ploidy * del_threshold_relative, 1)
    _total_cn_s = df['n_major'] + df['n_minor']
    _autosome_mask = df['chromosome'].isin(pcgr_vars.AUTOSOMES)

    # Homozygous deletions — autosomes
    df['homloss_cond'] = False
    df.loc[(_total_cn_s == 0) & _autosome_mask, "homloss_cond"] = True

    # Heterozygous deletions — autosomes (configurable threshold)
    df['hetloss_cond'] = False
    if threshold_mode == "absolute":
        logger.info(f"Using absolute threshold mode: hetdel = 0 < CN <= {del_threshold_absolute}")
        _hetloss_base = (_total_cn_s > 0) & (_total_cn_s <= del_threshold_absolute) & _autosome_mask
    elif threshold_mode == "relative":
        logger.info(
            f"Using relative threshold mode: hetdel = 0 < CN <= {tumor_ploidy} × {del_threshold_relative} = {relative_del_threshold_value}"
        )
        _hetloss_base = (_total_cn_s > 0) & (_total_cn_s <= relative_del_threshold_value) & _autosome_mask
    else:  # combined
        logger.info("Using combined threshold mode for hetdel (requires BOTH criteria):")
        logger.info(f"Absolute criterion: CN <= {del_threshold_absolute}")
        logger.info(f"Relative criterion: CN <= {tumor_ploidy} × {del_threshold_relative} = {relative_del_threshold_value}")
        _hetloss_base = (
            (_total_cn_s > 0) &
            (_total_cn_s <= del_threshold_absolute) &
            (_total_cn_s <= relative_del_threshold_value) &
            _autosome_mask
        )
    df.loc[_hetloss_base, "hetloss_cond"] = True

    # Sex chromosome losses
    df['hemloss_cond'] = False
    if sex == "UNKNOWN":
        warn_message("Argument --sex is 'UNKNOWN' - skipping copy number loss annotations on sex chromosomes", logger)
    elif sex == "MALE":
        df.loc[
            (_total_cn_s == 0) & df['chromosome'].isin(pcgr_vars.SEX_CHROMOSOMES), "hemloss_cond"] = True
    else:  # FEMALE
        _sex_mask = df['chromosome'].isin(pcgr_vars.SEX_CHROMOSOMES)
        df.loc[(_total_cn_s == 0) & _sex_mask, "homloss_cond"] = True
        if threshold_mode == "absolute":
            _sex_hetloss = (_total_cn_s > 0) & (_total_cn_s <= del_threshold_absolute) & _sex_mask
        elif threshold_mode == "relative":
            _sex_hetloss = (_total_cn_s > 0) & (_total_cn_s <= relative_del_threshold_value) & _sex_mask
        else:
            _sex_hetloss = (
                (_total_cn_s > 0) &
                (_total_cn_s <= del_threshold_absolute) &
                (_total_cn_s <= relative_del_threshold_value) &
                _sex_mask
            )
        df.loc[_sex_hetloss, "hetloss_cond"] = True

    # Aberration keys for ablations (used for biomarker matching)
    df.loc[df['homloss_cond'], 'aberration_key'] = \
        df.loc[df['homloss_cond'], 'entrezgene'].astype(str) + '_ablation'
    df.loc[df['hemloss_cond'], 'aberration_key'] = \
        df.loc[df['hemloss_cond'], 'entrezgene'].astype(str) + '_ablation'

    return df, baseline_cn


def _assign_variant_class(df: pd.DataFrame) -> pd.DataFrame:
    """
    Assign VARIANT_CLASS based on the condition flags set by the annotation functions.
    Requires: amp_cond, gain_cond, homloss_cond, hetloss_cond, hemloss_cond.
    """
    df['variant_class'] = 'neutral'
    df.loc[df['amp_cond'], 'variant_class'] = 'amplification'
    df.loc[df['gain_cond'], 'variant_class'] = 'gain'
    df.loc[df['homloss_cond'], 'variant_class'] = 'homdel'
    df.loc[df['hetloss_cond'], 'variant_class'] = 'hetdel'
    df.loc[df['hemloss_cond'], 'variant_class'] = 'hemdel'
    return df


def _annotate_loh(
    df: pd.DataFrame,
    baseline_cn: int,
    sex: str,
    logger: logging.Logger,
) -> pd.DataFrame:
    """
    Annotate Loss of Heterozygosity (LOH) subtypes.
    LOH: n_minor == 0 and n_major > 0 (complete loss of one allele).
    Subtypes: copy_neutral, deletion, amplification (based on n_major vs baseline).
    Adds: loh column.
    """
    df['loh'] = '.'

    loh_base = (df['n_minor'] == 0) & (df['n_major'] > 0)

    # Autosomes
    _auto = loh_base & df['chromosome'].isin(pcgr_vars.AUTOSOMES)
    df.loc[_auto & (df['n_major'] == baseline_cn), 'loh'] = 'copy-neutral'
    df.loc[_auto & (df['n_major'] < baseline_cn),  'loh'] = 'deletion'
    df.loc[_auto & (df['n_major'] > baseline_cn),  'loh'] = 'amplification'

    # Female sex chromosomes (XX baseline same as autosomes)
    if sex == "FEMALE":
        _sex = loh_base & df['chromosome'].isin(pcgr_vars.SEX_CHROMOSOMES)
        df.loc[_sex & (df['n_major'] == baseline_cn), 'loh'] = 'copy-neutral'
        df.loc[_sex & (df['n_major'] < baseline_cn),  'loh'] = 'deletion'
        df.loc[_sex & (df['n_major'] > baseline_cn),  'loh'] = 'amplification'

    loh_segments = df[df['loh'] != '.']
    if len(loh_segments) > 0:
        loh_unique = loh_segments.drop_duplicates(subset=['chromosome', 'segment_start', 'segment_end'])
        logger.info(f"LOH annotation: {len(loh_unique)} unique segments with LOH detected")
        for loh_type, count in loh_unique['loh'].value_counts().items():
            logger.info(f"LOH - {loh_type}: {count} segments")

    return df


def annotate_cna_segments(input_cna_segment_fname: str,
                          output_segment_gene_fname: str, 
                          output_segment_fname: str,
                          oncokb_input_fname: str,
                          output_dir: str, 
                          refdata_assembly_dir: str, 
                          build: str, 
                          sample_id: str,
                          sex: str,
                          amp_threshold_absolute: int = 5,
                          amp_threshold_relative: float = 2.5,
                          gain_threshold_absolute: int = 3,
                          gain_threshold_relative: float = 1.5,
                          del_threshold_absolute: int = 1,
                          del_threshold_relative: float = 0.5,
                          threshold_mode: str = "absolute",
                          transcript_overlap_fraction: float = 0.5,
                          tumor_ploidy: Optional[float] = None,
                          tumor_purity: Optional[float] = None,
                          expression_data: Optional[dict] = None,
                          logger: Optional[logging.Logger] = None) -> dict:
    """
    Annotate copy number aberrations in a given segment file with ploidy-aware thresholds.

    Args:
        input_cna_segment_fname (str): Path to the user-provided allele-specific copy number aberrations segment file.
        output_segment_gene_fname (str): File name of the annotated transcript-level output file.
        output_segment_fname (str): File name of the annotated segment-level output file.
        oncokb_input_fname: str, File name of the OncoKB CNA input file to be generated.
        output_dir (str): Path to the output directory.
        refdata_assembly_dir (str): Path to the build-specific PCGR database directory.
        build (str): Genome assembly build of input data.
        sample_id (str): Sample identifier
        sex (str): Sex of the patient (MALE, FEMALE, or UNKNOWN)
        amp_threshold_absolute (int, optional): Absolute copy number threshold for amplifications. Defaults to 5.
            Used in 'absolute' and 'combined' modes. Valid range: 3-20.
        amp_threshold_relative (float, optional): Relative threshold as fold-change over tumor ploidy. Defaults to 2.5.
            Used in 'relative' and 'combined' modes. Valid range: 1.2-5.0.
            Example: If ploidy=2.0 and relative=2.5, amplification threshold = 2.0 × 2.5 = 5.0
        gain_threshold_absolute (int, optional): Absolute CN threshold for gains. Defaults to 3.
        gain_threshold_relative (float, optional): Fold-change over ploidy threshold for gains. Defaults to 1.5.
        del_threshold_absolute (int, optional): Absolute CN at or below which (but > 0) a segment is a hetdel. Defaults to 1.
        del_threshold_relative (float, optional): Fold-change below ploidy at or below which a segment is a hetdel. Defaults to 0.5.
        threshold_mode (str, optional): Thresholding mode applied uniformly to all CNA tiers. Defaults to "absolute".
            Options:
            - "absolute": use absolute copy number thresholds
            - "relative": use fold-change over tumor ploidy thresholds
            - "combined": both absolute AND relative criteria must be satisfied (most conservative)
        transcript_overlap_fraction (float, optional): Fraction of overlap required for gene annotation. Defaults to 0.5.
        tumor_ploidy (float, optional): Tumor ploidy value. Defaults to None (auto-estimated from segments).
            Valid range: 1.0-8.0. If not provided in 'relative' or 'combined' modes, will be auto-estimated
            using genome-wide weighted median of segment copy numbers.
        tumor_purity (float, optional): Tumor purity value. Defaults to None.
        snv_df_somatic (pd.DataFrame, optional): DataFrame of somatic SNV/InDel variants for two-hit candidate annotation. Defaults to None.
        snv_df_germline (pd.DataFrame, optional): DataFrame of germline SNV/InDel variants (from CPSR) for two-hit candidate annotation. Defaults to None.
        expression_data (dict, optional): Expression data for integration. Defaults to None.
        logger (logging.Logger, optional): Logger instance. Defaults to None.

    Returns:
        dict with keys:
            'status': 0 if successful, -1 if error
            'amp_threshold_effective': the resolved effective amplification threshold (float or None on error)
            'tumor_ploidy': the ploidy value actually used (float or None on error)
            'tumor_ploidy_source': 'provided', 'estimated', or 'NA' on error
    """
    _error_result = {
        'status': -1,
        'amp_threshold_effective': None,
        'tumor_ploidy': tumor_ploidy,
        'tumor_ploidy_source': 'NA',
    }
    
    #logger = getlogger("pcgr-annotate-cna-segments")
    pybedtools.set_tempdir(output_dir)
    if logger is None:
        logger = logging.getLogger("pcgr-annotate-cna-segments")
    
    temp_files = []
    
    # Read user-defined copy number aberrations segment file
    check_file_exists(input_cna_segment_fname, logger = logger)
    cna_query_segment_df = pd.read_csv(input_cna_segment_fname, sep="\t", na_values=".")
    
    ## Check that required columns are present
    if not {'Chromosome','Start','End','nMajor','nMinor'}.issubset(cna_query_segment_df.columns):
        err_msg = f"Could not find required columns in CNA segment file: {input_cna_segment_fname} - exiting."
        error_message(err_msg, logger)
    
    cna_query_segment_df = cna_query_segment_df[['Chromosome', 'Start','End','nMajor','nMinor']]
    
    ## round nMajor and nMinor to integers
    cna_query_segment_df['nMajor'] = cna_query_segment_df['nMajor'].round(0).astype(int)
    cna_query_segment_df['nMinor'] = cna_query_segment_df['nMinor'].round(0).astype(int)

    ## Ensure cna_query_segment_df is a pd.DataFrame
    cna_query_segment_df = cast(pd.DataFrame, cna_query_segment_df)

    for col in ['Chromosome', 'Start', 'End', 'nMajor', 'nMinor']:
        cna_query_segment_df[col] = pd.Series(cna_query_segment_df[col]).astype("string")
    
    ## Remove 'chr' prefix from chromosome names
    cna_query_segment_df['Chromosome'] = \
        cna_query_segment_df['Chromosome'].str.replace("^(C|c)hr", "", regex = True)
    
    ## Only consider chromosome names that are present in the nuclear chromosomes list
    cna_query_segment_df = cna_query_segment_df[cna_query_segment_df['Chromosome'].isin(nuclear_chromosomes)]
    if len(cna_query_segment_df) == 0:
        warn_msg = f"Could not find any CNA query segments listed on nuclear chromosomes: {nuclear_chromosomes} - returning."
        warn_message(warn_msg, logger)
        return _error_result
    
    cna_query_segment_df['segment_id'] = \
        'chr' + cna_query_segment_df['Chromosome'].str.cat(
            cna_query_segment_df['Start'].astype(str), sep = ":").str.cat(
                cna_query_segment_df['End'].astype(str), sep='-')
    
    ## Create Name column of BED file
    cna_query_segment_df["Name"] = \
        cna_query_segment_df["segment_id"].str.cat(
            cna_query_segment_df["nMajor"], sep="|").str.cat(
                cna_query_segment_df["nMinor"], sep = "|")
    cna_query_segment_df = cna_query_segment_df[['Chromosome','Start','End','Name']]
    
    
    ## Check that 'End' of segments do not exceed chromosome lengths
    chromsizes_fname = \
        os.path.join(refdata_assembly_dir, 'chromsize.' + build + '.tsv')
    
    check_file_exists(chromsizes_fname, logger = logger)
    chromsizes = pd.read_csv(chromsizes_fname, sep="\t", skiprows=1, names=['Chromosome', 'ChromLength','CentromereStart','CentromereEnd'])
    chromsizes = chromsizes[['Chromosome','ChromLength']]
    if not isinstance(cna_query_segment_df, pd.DataFrame):
        cna_query_segment_df = pd.DataFrame(cna_query_segment_df)
    
    cna_query_segment_df = cna_query_segment_df.merge(
        chromsizes, left_on=["Chromosome"], right_on=["Chromosome"], how="left")

    # Convert Start and End back to integers for chromosome length comparisons
    cna_query_segment_df['Start'] = cna_query_segment_df['Start'].astype(int)
    cna_query_segment_df['End'] = cna_query_segment_df['End'].astype(int)

    segments_beyond_chromlength = \
        cna_query_segment_df[cna_query_segment_df['End'] > cna_query_segment_df['ChromLength']]
    
    ## Issue warning if segments exceed chromosome lengths
    if segments_beyond_chromlength.empty is False:
        warn_msg = f"Ignoring parts of n = {len(segments_beyond_chromlength)} copy number segments that " + \
            f"exceed the chromosomal lengths of {build}"
        warn_message(warn_msg, logger)
           
    cna_query_segment_df.loc[cna_query_segment_df['End'] > cna_query_segment_df['ChromLength'],"End"] = \
        cna_query_segment_df.loc[cna_query_segment_df['End'] > cna_query_segment_df['ChromLength'],'ChromLength']
        
    cna_query_segment_df = \
        cna_query_segment_df[cna_query_segment_df['Start'] <= cna_query_segment_df['ChromLength']]
    
    cna_query_segment_df = cna_query_segment_df[['Chromosome','Start','End','Name']]    
    
    ## transform cna segments to pybedtools object
    cna_query_segment_bed = pybedtools.BedTool.from_dataframe(cna_query_segment_df)
    temp_files.append(cna_query_segment_bed.fn)
    
    ## annotate segments with cytobands
    cna_query_segment_df = annotate_cytoband(cna_query_segment_bed, output_dir, refdata_assembly_dir, logger)
    
    cna_query_segment_df['chromosome2'] = cna_query_segment_df['chromosome'].astype(object)
    cna_query_segment_df.loc[cna_query_segment_df['chromosome2'] == "X","chromosome2"] = 23
    cna_query_segment_df.loc[cna_query_segment_df['chromosome2'] == "Y","chromosome2"] = 24
    cna_query_segment_df['chromosome2'] = cna_query_segment_df['chromosome2'].astype(int)
        
    cna_query_segment_df = cna_query_segment_df.sort_values(['chromosome2','segment_start'], ascending=True)

    cna_query_segment_df = cna_query_segment_df.drop(columns=['chromosome2'])

    ## annotate with protein-coding transcripts
    cna_query_segment_bed = pybedtools.BedTool.from_dataframe(cna_query_segment_df)
    temp_files.append(cna_query_segment_bed.fn)

    cna_query_segment_df = annotate_transcripts(
       cna_query_segment_bed, output_dir, refdata_assembly_dir, transcript_overlap_fraction=transcript_overlap_fraction, logger=logger)
    
    if cna_query_segment_df.empty is True:
        warn_msg = "Could not find any protein-coding gene annotations for input CNA segments - omitting gene-level annotations"
        warn_message(warn_msg, logger)
        return _error_result
    cna_query_segment_df['segment_length_mb'] = \
        ((cna_query_segment_df['segment_end'] - cna_query_segment_df['segment_start']) / 1e6).astype(float).round(4)
    
    ## load copy-number biomarker evidence
    biomarkers = load_all_biomarkers(
        logger, refdata_assembly_dir, biomarker_vartype = 'CNA', biomarker_variant_origin = 'Somatic')

    ## Determine amplification threshold based on mode
    #if threshold_mode == "relative" or threshold_mode == "combined":
    # Auto-estimate ploidy if not provided
    tumor_ploidy_source = 'provided'
    effective_amplification_threshold = None
    if tumor_ploidy is None or tumor_ploidy == 'NA':
        # Need to read the original segment file to get nMajor/nMinor before annotation
        temp_segment_df = pd.read_csv(input_cna_segment_fname, sep="\t", na_values=".")
        tumor_ploidy = estimate_tumor_ploidy(temp_segment_df, logger=logger)
        tumor_ploidy_source = 'estimated'
        logger.info(f"Auto-estimated tumor ploidy: {tumor_ploidy} (genome-wide weighted median)")
    else:
        logger.info(f"Using user-specified tumor ploidy: {tumor_ploidy}")

    cna_query_segment_df, effective_amplification_threshold = _annotate_amplifications(
        cna_query_segment_df, tumor_ploidy,
        amp_threshold_absolute, amp_threshold_relative, threshold_mode, logger)
    cna_query_segment_df = _annotate_gains(
        cna_query_segment_df, tumor_ploidy,
        gain_threshold_absolute, gain_threshold_relative, threshold_mode, logger)
    cna_query_segment_df, baseline_cn = _annotate_deletions(
        cna_query_segment_df, tumor_ploidy,
        del_threshold_absolute, del_threshold_relative, threshold_mode, sex, logger)
    cna_query_segment_df = _assign_variant_class(cna_query_segment_df)
    cna_query_segment_df = _annotate_loh(cna_query_segment_df, baseline_cn, sex, logger)

    ## Append actionability evidence to input amplifications (column 'biomarker_match')
    #cna_query_segment_df = cna_query_segment_df.merge(
    #    cna_actionable_df, left_on=["aberration_key"], right_on=["aberration_key"], how="left")
    cna_query_segment_df = cna_query_segment_df.merge(
        biomarkers['actionable_df'], left_on=["aberration_key"], right_on=["aberration_key"], how="left")
    cna_query_segment_df.drop(['amp_cond', 'gain_cond', 'hetloss_cond', 'homloss_cond', 'hemloss_cond', 'aberration_key'], axis=1, inplace=True)
    cna_query_segment_df.loc[cna_query_segment_df['biomarker_match'].isnull(),"biomarker_match"] = '.'
    
    ## remove all temporary files
    for fname in temp_files:
        remove_file(fname)
    
    cna_query_segment_df.columns = [col.upper() for col in cna_query_segment_df.columns]
    #cna_query_segment_df.columns = map(str.upper, cna_query_segment_df.columns)

    segments_out = cna_query_segment_df[['CHROMOSOME','SEGMENT_START','SEGMENT_END']].copy()
    segments_out['SEGMENT_NAME'] = (
        cna_query_segment_df['SEGMENT_ID'].astype(str)
        .str.cat(cna_query_segment_df['N_MAJOR'].astype(str), sep='|')
        .str.cat(cna_query_segment_df['N_MINOR'].astype(str), sep='|')
        .str.cat(cna_query_segment_df['CHROMOSOME_ARM'].astype(str), sep='|')
        .str.cat(cna_query_segment_df['CYTOBAND'].astype(str), sep='|')
        .str.cat(cna_query_segment_df['EVENT_TYPE'].astype(str), sep='|')
        .str.cat(cna_query_segment_df['VARIANT_CLASS'].astype(str), sep='|')
        .str.cat(cna_query_segment_df['LOH'].astype(str).replace('.', ''), sep='|')
    )
    segments_out = segments_out.rename(columns={'CHROMOSOME': 'CHROM'})
    segments_out.drop_duplicates(inplace=True)
    segments_out.to_csv(output_segment_fname, sep="\t", header=True, index=False)

    cna_query_segment_df.rename(columns = {
        'CHROMOSOME':'CHROM',
        'SEGMENT_ID':'VAR_ID',
        'N_MAJOR':'CN_MAJOR',
        'N_MINOR':'CN_MINOR'
    }, inplace = True)
    cna_query_segment_df['VAR_ID'] = \
        cna_query_segment_df['VAR_ID'].str.cat(
            cna_query_segment_df['CN_MAJOR'].astype(str), sep=":").str.cat(
                cna_query_segment_df['CN_MINOR'].astype(str), sep=":")
    
    cna_query_segment_df['SAMPLE_ID'] = sample_id
    cna_query_segment_df['TUMOR_PLOIDY'] = float(tumor_ploidy)
    cna_query_segment_df['TUMOR_PLOIDY_SOURCE'] = tumor_ploidy_source
    cna_query_segment_df['TUMOR_PURITY'] = tumor_purity if tumor_purity is not None else '.'
    if expression_data is not None:
        cna_query_segment_df = integrate_variant_expression(cna_query_segment_df, expression_data, logger)
    else:
        logger.info("No expression data provided. Skipping CNA-expression integration")
        cna_query_segment_df['TPM_GENE'] = '.'
        cna_query_segment_df['TPM'] = '.'


    cna_query_segment_df.to_csv(output_segment_gene_fname, sep="\t", header=True, index=False)

    # Generate OncoKB input file
    generate_oncokb_cna_input(cna_query_segment_df, oncokb_input_fname, logger)

    return {
        'status': 0,
        'amp_threshold_effective': float(effective_amplification_threshold),
        'tumor_ploidy': float(tumor_ploidy),
        'tumor_ploidy_source': tumor_ploidy_source,
    }


def annotate_cytoband(
    cna_segments_bt: BedTool, 
    output_dir: str, refdata_assembly_dir: str, 
    logger: Optional[logging.Logger] = None) -> pd.DataFrame:
    
    pybedtools.set_tempdir(output_dir)    
    temp_files = []
    
    # BED file with cytoband annotations
    cytoband_bed_fname = \
        os.path.join(refdata_assembly_dir, 'misc','bed','cytoband', 'cytoband.bed.gz')
    
    cytoband_annotated_segments = pd.DataFrame()
    
    check_file_exists(cytoband_bed_fname, logger = logger)
    cytoband_bed = pybedtools.BedTool(cytoband_bed_fname)
        
    tmp = cna_segments_bt.intersect(cytoband_bed, loj=True)        
    cytoband_intersect = pd.read_csv(
        tmp.fn, sep="\t", header=None, 
        names=['chromosome','segment_start','segment_end','segment_name',
                'cytoband_chrom','cytoband_start','cytoband_end','cytoband_name'])
    temp_files.append(tmp.fn)
    tmp1 = cytoband_intersect.groupby(['chromosome','segment_start','segment_end']).first()
    segments_first = tmp1.index.to_frame()
    tmp2 = cytoband_intersect.groupby(['chromosome','segment_start','segment_end']).last()
    segments_last = tmp2.index.to_frame()
        
    cytoband_first = pd.concat([segments_first, tmp1], axis = 1)
    cytoband_first = cytoband_first.reset_index(drop = True)
    cytoband_first = \
        cytoband_first[['chromosome','segment_start','segment_end','segment_name','cytoband_name']]
    cytoband_first.rename({'cytoband_name': 'cytoband_first'}, axis=1, inplace=True)
        
    cytoband_last = pd.concat([segments_last, tmp2], axis = 1)
    cytoband_last = cytoband_last.reset_index(drop = True)
    cytoband_last = \
        cytoband_last[['chromosome','segment_start','segment_end','segment_name','cytoband_name']]
    cytoband_last.rename({'cytoband_name': 'cytoband_last'}, axis=1, inplace=True)
        
    segments_cytoband = cytoband_first.merge(
        cytoband_last, left_on=["chromosome","segment_start","segment_end","segment_name"], 
        right_on=["chromosome","segment_start","segment_end","segment_name"], how="left")
        
    cytoband_first_annotations = segments_cytoband.cytoband_first.str.split('\\|', expand = True)
    cytoband_first_annotations.columns = ['first_cytoband','first_arm','first_arm_length','first_focal_threshold']
    cytoband_first_annotations = cytoband_first_annotations.astype({'first_focal_threshold':'int'})
    cytoband_last_annotations = segments_cytoband.cytoband_last.str.split('\\|', expand = True)
    cytoband_last_annotations.columns = ['last_cytoband','last_arm','last_arm_length','last_focal_threshold']
        
    cytoband_all = pd.concat([segments_cytoband, cytoband_first_annotations, cytoband_last_annotations], axis = 1)
    
    cytoband_all['segment_start'] = cytoband_all['segment_start'].astype(int)
    cytoband_all['segment_end'] = cytoband_all['segment_end'].astype(int)       
    cytoband_all.loc[:,'segment_length'] = cytoband_all['segment_end'] - cytoband_all['segment_start']
    cytoband_all['event_type'] =  'broad'
    cytoband_all.loc[cytoband_all['segment_length'] < cytoband_all['first_focal_threshold'], "event_type"] = "focal"
        
    cytoband_all['cytoband'] = cytoband_all['first_cytoband'].where(
    cytoband_all['first_cytoband'] == cytoband_all['last_cytoband'], 
    cytoband_all['first_cytoband'] + "-" + cytoband_all['last_cytoband'])
        
    cytoband_all['segment_name'] = cytoband_all['segment_name'].str.cat(cytoband_all['first_arm'], sep = "|").str.cat(
        cytoband_all['cytoband'],sep="|").str.cat(cytoband_all['event_type'], sep="|")    
    cytoband_annotated_segments = cytoband_all[['chromosome','segment_start','segment_end','segment_name']]
    
    ## remove all temporary files
    for fname in temp_files:
        remove_file(fname)
            
    return cytoband_annotated_segments


def annotate_transcripts(cna_segments_bt: BedTool, 
                         output_dir: str,
                         refdata_assembly_dir: str, 
                         transcript_overlap_fraction: float,
                         logger: logging.Logger) -> pd.DataFrame:
    
    """
    Annotate the CNA segments with gene transcripts and return the annotated segments as a DataFrame.
    
    Parameters:
    - cna_segments_bt: A BedTool object representing the CNA segments.
    - output_dir: A string specifying the output directory (for temporary files generation).
    - refdata_assembly_dir: A string specifying the directory of the build-specific PCGR data bundle.
    - transcript_overlap_fraction: A float representing the fraction of overlap required for annotation.
    
    Returns:
    - cna_segments_annotated: A DataFrame containing the annotated CNA segments.
    """
    
    pybedtools.set_tempdir(output_dir)
    temp_files = []
        
    # BED file with protein-coding transcripts
    gene_transcript_bed_fname = \
        os.path.join(refdata_assembly_dir, 'gene','bed','gene_transcript_xref', 'gene_transcript_xref_pc_nopad.bed.gz')
    gene_xref_tsv_fname = \
        os.path.join(refdata_assembly_dir, "gene", "tsv", "gene_transcript_xref", "gene_transcript_xref.tsv.gz")
        
    cna_segments_annotated = pd.DataFrame()

    if os.path.exists(gene_transcript_bed_fname):
        gene_transcript_bed = pybedtools.BedTool(gene_transcript_bed_fname)               
        
        ## Intersect query segments with annotated gene transcripts
        transcript_cna_annotations = gene_transcript_bed.intersect(cna_segments_bt, wao=True, f=transcript_overlap_fraction)
        temp_files.append(transcript_cna_annotations.fn)
        
        if os.path.exists(str(transcript_cna_annotations.fn)) and os.path.getsize(str(transcript_cna_annotations.fn)) > 0:
            colnames = ['chromosome','transcript_start','transcript_end','transcript_annotations',
                        'chromosome2','segment_start','segment_end','segment_name','bp_overlap']
            cna_transcript_annotations = \
                pd.read_csv(str(transcript_cna_annotations.fn), sep="\t", na_values=".", 
                            names=list(colnames), header=None, low_memory=False)
            
            ## ignore transcripts that do not overlap with any copy number segment
            cna_transcript_annotations = \
                cna_transcript_annotations[cna_transcript_annotations['segment_start'] != -1]

            if cna_transcript_annotations.empty is True:                
                return cna_segments_annotated
            
            if cna_transcript_annotations.empty is True:
                #warn_msg = f"Could not find any protein-coding gene annotations for input CNA segments - returning."
                #warn_message(warn_msg, logger)
                return cna_segments_annotated
        
            ## calculate fraction of overlap between segments and transcripts
            cna_transcript_annotations['transcript_overlap_percent'] = \
                cna_transcript_annotations['transcript_end'] - cna_transcript_annotations['transcript_start']
            cna_transcript_annotations['transcript_overlap_percent'] = ((cna_transcript_annotations['bp_overlap'] / \
                cna_transcript_annotations['transcript_overlap_percent']) * 100).astype(float).round(2)
                               
            ## pull out segment identifier, major and minor copy number from the name column, and cytoband information
            major_minor_annotations_df = cna_transcript_annotations.segment_name.str.split('\\|', expand = True)
            major_minor_annotations_df.columns = ['segment_id','n_major','n_minor','chromosome_arm','cytoband','event_type']
            cna_transcript_annotations = pd.concat([cna_transcript_annotations, major_minor_annotations_df], axis = 1)
            
            ## select only relevant columns
            core_segment_annotations = \
                cna_transcript_annotations[['chromosome','segment_start','segment_end','segment_id',
                                            'n_major','n_minor','chromosome_arm','cytoband','event_type',
                                            'transcript_overlap_percent',
                                            'transcript_start','transcript_end']]
            core_segment_annotations = core_segment_annotations.astype({'n_major':'int','n_minor':'int'})
            ## pull out gene cross-reference annotations
            gene_annotations = cna_transcript_annotations.transcript_annotations.str.split('\\|', expand = True)
            gene_annotations.columns = ["ensembl_transcript_id", "ensembl_gene_id","symbol","entrezgene",
                                      "refseq_transcript_id","actionable_gene", "tsg", "tsg_support","tsg_rank",
                                      "oncogene","oncogene_support","oncogene_rank"]
            for elem in ['actionable_gene','oncogene','tsg']:
                gene_annotations = gene_annotations.astype({elem:'string'})
                gene_annotations.loc[gene_annotations[elem] == "", elem] = "FALSE"
                gene_annotations.loc[gene_annotations[elem] == "1", elem] = "TRUE"
            
            for elem in ['symbol','ensembl_gene_id','entrezgene','refseq_transcript_id','tsg_support',
                         'tsg_rank','oncogene_support','oncogene_rank']:
                gene_annotations = gene_annotations.astype({elem:'string'})
                gene_annotations.loc[gene_annotations[elem] == "", elem] = "."
            
            cna_segments_annotated = pd.concat([core_segment_annotations, gene_annotations], axis = 1)
            
            if {'genename'}.issubset(cna_segments_annotated.columns):
                cna_segments_annotated.drop('genename', inplace=True, axis=1)
            
            if os.path.exists(gene_xref_tsv_fname):
                gene_xref_df = pd.read_csv(
                    str(gene_xref_tsv_fname), sep="\t", na_values=".", usecols=["entrezgene","name"], low_memory=False)
                gene_xref_df = gene_xref_df[gene_xref_df['entrezgene'].notnull()].drop_duplicates()
                gene_xref_df["entrezgene"] = gene_xref_df["entrezgene"].astype("int64").astype("string")
                gene_xref_df.rename(columns = {'name':'genename'}, inplace = True)                                        
                cna_segments_annotated = cna_segments_annotated.merge(
                    gene_xref_df, left_on=["entrezgene"], right_on=["entrezgene"], how="left")
                cna_segments_annotated["entrezgene"] = cna_segments_annotated['entrezgene'].str.replace("\\.[0-9]{1,}$", "", regex = True)
                cna_segments_annotated = cna_segments_annotated.fillna('.')

                # Drop records with missing Entrezgene identifier
                n_before = len(cna_segments_annotated)
                cna_segments_annotated = cna_segments_annotated[
                    cna_segments_annotated['entrezgene'] != '.'
                ].copy()
                n_missing_entrezgene = n_before - len(cna_segments_annotated)
                if n_missing_entrezgene > 0:
                    logger.info(
                        f"Removed {n_missing_entrezgene} transcript records with missing Entrezgene identifier"
                    )

                # When a single Entrezgene maps to more than one Ensembl gene ID within
                # the same segment, keep ALL transcripts of the best Ensembl gene ID and
                # discard the others. The best gene ID is the one whose transcripts have the
                # highest maximum transcript_overlap_percent.
                n_before_dedup = len(cna_segments_annotated)

                # For each (segment_id, entrezgene, ensembl_gene_id), get best overlap
                gene_max_overlap = (
                    cna_segments_annotated
                    .groupby(['segment_id', 'entrezgene', 'ensembl_gene_id'],
                             sort=False)['transcript_overlap_percent']
                    .max()
                    .reset_index(name='_max_overlap')
                )
                # For each (segment_id, entrezgene), select the ensembl_gene_id with
                # the highest max overlap (ties broken by first occurrence)
                best_gene_id = (
                    gene_max_overlap
                    .sort_values('_max_overlap', ascending=False)
                    .drop_duplicates(subset=['segment_id', 'entrezgene'], keep='first')
                    [['segment_id', 'entrezgene', 'ensembl_gene_id']]
                )
                # Keep only transcripts belonging to the selected ensembl_gene_id
                cna_segments_annotated = (
                    cna_segments_annotated
                    .merge(best_gene_id,
                           on=['segment_id', 'entrezgene', 'ensembl_gene_id'],
                           how='inner')
                    .reset_index(drop=True)
                )
                n_removed_dups = n_before_dedup - len(cna_segments_annotated)
                if n_removed_dups > 0:
                    logger.info(
                        f"Removed {n_removed_dups} transcript records where Entrezgene mapped "
                        f"to multiple Ensembl gene IDs within a segment"
                    )
                    logger.info(
                        f"Retained all transcripts "
                        f"of the representative gene ID (highest transcript overlap percent)"
                    )

                # Identify transcripts that span multiple segments (for debugging)
                # At this point each (ensembl_transcript_id, segment_id) is already unique,
                # so any transcript_id appearing more than once crosses a segment boundary.
                transcript_segment_counts = (
                    cna_segments_annotated
                    .groupby('ensembl_transcript_id')['segment_id']
                    .nunique()
                )
                multi_segment_transcripts = transcript_segment_counts[transcript_segment_counts > 1].index
                if len(multi_segment_transcripts) > 0:
                    logger.debug(
                        f"{len(multi_segment_transcripts)} transcript(s) span multiple CNA segments"
                    )
                    multi_seg_df = (
                        cna_segments_annotated[
                            cna_segments_annotated['ensembl_transcript_id'].isin(multi_segment_transcripts)
                        ][['symbol', 'ensembl_transcript_id', 'segment_id',
                           'chromosome', 'segment_start', 'segment_end',
                           'transcript_overlap_percent', 'n_major', 'n_minor']]
                        .sort_values(['symbol', 'transcript_overlap_percent'], ascending=[True, False])
                    )
                    for tx_id, grp in multi_seg_df.groupby('ensembl_transcript_id'):
                        symbol = grp['symbol'].iloc[0]
                        best = grp.iloc[0]
                        others = grp.iloc[1:]
                        logger.debug(
                            f"  {symbol} ({tx_id}): majority segment "
                            f"{best['chromosome']}:{best['segment_start']}-{best['segment_end']} "
                            f"[{best['transcript_overlap_percent']:.1f}% overlap, "
                            f"CN={int(best['n_major'])+int(best['n_minor'])}]; "
                            f"also overlaps: "
                            + ", ".join(
                                f"{r['chromosome']}:{r['segment_start']}-{r['segment_end']} "
                                f"[{r['transcript_overlap_percent']:.1f}%]"
                                for _, r in others.iterrows()
                            )
                        )

                # For transcripts spanning multiple segments, keep only the majority segment
                cna_segments_annotated = (
                    cna_segments_annotated
                    .sort_values('transcript_overlap_percent', ascending=False)
                    .drop_duplicates(subset=['ensembl_transcript_id'], keep='first')
                    .reset_index(drop=True)
                )
            else:
                err_msg = f"Could not find {gene_xref_tsv_fname} needed for gene name annotation - exiting"
                error_message(err_msg, logger)                
        else:
            err_msg = f"{transcript_cna_annotations.fn} is empty or not found - exiting"
            error_message(err_msg, logger)
    else:
        err_msg = f"{gene_transcript_bed_fname} is empty or not found - exiting"
        error_message(err_msg, logger)
    
    ## remove all temporary files
    for fname in temp_files:
        remove_file(fname)
    
    return cna_segments_annotated
