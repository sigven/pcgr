#!/usr/bin/env python

import os
import pandas as pd
from typing import Optional
import logging

from pcgr.utils import check_file_exists, remove_file
from pcgr.pcgr_vars import ONCOKB_COLS


def append_oncokb_snv_annotations(
        tsv_gz_fname: str,
        oncokb_maf_hgvsp: Optional[str] = None,
        oncokb_maf_hgvsg: Optional[str] = None,
        logger: Optional[logging.Logger] = None) -> None:
    """
    Merge OncoKB SNV/InDel annotations into the final PCGR TSV.

    Uses VAR_ID (CHROM_POS_REF_ALT, VCF-format) as the join key.
    HGVSp hits are preferred; HGVSg hits fill gaps for variants that
    were only matched by genomic coordinates.

    If OncoKB was not run (both MAF paths are None/missing), the four
    *_OKB columns are added with NA values so downstream R code can
    always rely on their presence.

    Args:
        tsv_gz_fname:    Path to the gzipped variant TSV (read + overwritten)
        oncokb_maf_hgvsp: Path to OncoKB MafAnnotator HGVSp output
        oncokb_maf_hgvsg: Path to OncoKB MafAnnotator HGVSg output
        logger:          Logger instance
    """
    if logger is None:
        logger = logging.getLogger("pcgr-oncokb-snv-append")

    okb_source_cols = list(ONCOKB_COLS.keys())
    okb_renamed_cols = list(ONCOKB_COLS.values())

    if not os.path.exists(tsv_gz_fname):
        logger.warning(f"Variant TSV not found: {tsv_gz_fname} - skipping OncoKB SNV merge")
        return

    tsv_df = pd.read_csv(tsv_gz_fname, sep="\t", low_memory=False)

    # If VAR_ID is missing from TSV (e.g. older run), add NA columns and return
    if 'VAR_ID' not in tsv_df.columns:
        logger.warning("VAR_ID column absent from variant TSV - cannot merge OncoKB SNV annotations")
        for col in okb_renamed_cols:
            tsv_df[col] = pd.NA
        tsv_df.to_csv(tsv_gz_fname, sep="\t", compression="gzip", index=False)
        return

    def _read_oncokb_maf(path: str, label: str) -> Optional[pd.DataFrame]:
        """Read an OncoKB-annotated MAF, return VAR_ID + annotation cols or None."""
        if not path or not os.path.exists(path):
            return None
        try:
            df = pd.read_csv(path, sep="\t", low_memory=False)
        except Exception as e:
            logger.warning(f"Could not read OncoKB MAF ({label}): {e}")
            return None
        if 'VAR_ID' not in df.columns:
            logger.warning(f"VAR_ID absent from OncoKB MAF ({label}) - cannot use for merge")
            return None
        present = [c for c in okb_source_cols if c in df.columns]
        if not present:
            return None
        keep = ['VAR_ID'] + present
        return df[keep].drop_duplicates(subset='VAR_ID')

    hgvsp_df = _read_oncokb_maf(oncokb_maf_hgvsp, "HGVSp")
    hgvsg_df = _read_oncokb_maf(oncokb_maf_hgvsg, "HGVSg")

    # Combine: HGVSp takes priority; HGVSg fills remaining gaps
    if hgvsp_df is not None and hgvsg_df is not None:
        oncokb_df = pd.concat([hgvsp_df, hgvsg_df], ignore_index=True) \
                      .drop_duplicates(subset='VAR_ID', keep='first')
    elif hgvsp_df is not None:
        oncokb_df = hgvsp_df
    elif hgvsg_df is not None:
        oncokb_df = hgvsg_df
    else:
        oncokb_df = None

    # Drop any pre-existing OKB columns to make this function idempotent
    tsv_df = tsv_df.drop(columns=[c for c in okb_renamed_cols if c in tsv_df.columns])

    if oncokb_df is not None and len(oncokb_df) > 0:
        # Ensure all source cols are present (fill absent ones with NA)
        for col in okb_source_cols:
            if col not in oncokb_df.columns:
                oncokb_df[col] = pd.NA
        oncokb_df = oncokb_df.rename(columns=ONCOKB_COLS)
        tsv_df = tsv_df.merge(oncokb_df[['VAR_ID'] + okb_renamed_cols],
                              on='VAR_ID', how='left')
        n_hits = tsv_df[okb_renamed_cols[0]].notna().sum()
        logger.info(f"OncoKB SNV annotations merged: {n_hits} / {len(tsv_df)} variants annotated")
    else:
        logger.info("No OncoKB SNV annotations available - adding NA columns")
        for col in okb_renamed_cols:
            tsv_df[col] = pd.NA

    tsv_df.to_csv(tsv_gz_fname, sep="\t", compression="gzip", index=False)
    logger.info(f"OncoKB SNV columns written to {os.path.basename(tsv_gz_fname)}")


def update_maf(maf_tmp_fname: str,
               maf_fname: str, 
               allelic_support_tags: dict,
               logger = None,
               update_allelic_support = False,
               debug = False):

    """
    Update MAF file from vcf2maf.pl with allelic support data (t_depth, t_ref_count, t_alt_count etc).
    
    Args:
        maf_tmp_fname (str): File name of the temporary MAF file.
        maf_fname (str): File name of the final MAF file.
        allelic_support_tags (dict): Dictionary of allelic support tags (encoded in VCF INFO field and retained in MAF).
        logger: Logger object for logging messages.
        update_allelic_support (bool): Flag indicating whether to update allelic support.
    Returns:
        int: 0 if successful.
    """
    
    # Read MAF file generated with vcf2maf.pl
    check_file_exists(maf_tmp_fname, logger)
    
    header_line = "#version 2.4"
    with open(maf_tmp_fname) as f:
        header_line = f.readline().strip('\n')
    f.close()
        
    raw_maf_data = pd.read_csv(maf_tmp_fname, sep="\t", header=1, dtype='string',na_values=['.'], low_memory=False)
    if update_allelic_support is False:
        # write to file
        os.rename(maf_tmp_fname, maf_fname)
    else:
    
        if 'tumor_dp_tag' in allelic_support_tags:
            if allelic_support_tags['tumor_dp_tag'] != "_NA_":                               
                if {allelic_support_tags['tumor_dp_tag']}.issubset(raw_maf_data.columns):
                    if raw_maf_data[raw_maf_data[allelic_support_tags['tumor_dp_tag']].isna()].empty:

                        raw_maf_data.loc[:,"t_depth"] = raw_maf_data.loc[:,allelic_support_tags['tumor_dp_tag']]
                        
                        if 'tumor_af_tag' in allelic_support_tags: 
                            if allelic_support_tags['tumor_af_tag'] != "_NA_":    
                                if {allelic_support_tags['tumor_af_tag']}.issubset(raw_maf_data.columns):                                   
                                    if raw_maf_data[raw_maf_data[allelic_support_tags['tumor_af_tag']].isna()].empty:
                                        raw_maf_data['t_alt_count'] = None
                                        raw_maf_data.loc[:,"t_alt_count"] = \
                                            raw_maf_data.loc[:,allelic_support_tags['tumor_af_tag']].astype(float) * raw_maf_data.loc[:,"t_depth"].astype(int)
                                    
                                        raw_maf_data.loc[:,"t_alt_count"] = round(raw_maf_data.loc[:,"t_alt_count"].astype(float),0).astype(int)
                                        raw_maf_data['t_ref_count'] = None
                                        raw_maf_data.loc[:,"t_ref_count"] = \
                                            raw_maf_data.loc[:,"t_depth"].astype(int) - raw_maf_data.loc[:,"t_alt_count"]
        
        if 'control_dp_tag' in allelic_support_tags:
            if allelic_support_tags['control_dp_tag'] != "_NA_":                    
                if {allelic_support_tags['control_dp_tag']}.issubset(raw_maf_data.columns):
                    if raw_maf_data[raw_maf_data[allelic_support_tags['control_dp_tag']].isna()].empty:
                        raw_maf_data.loc[:,"n_depth"] = raw_maf_data.loc[:,allelic_support_tags['control_dp_tag']]
                        if 'control_af_tag' in allelic_support_tags: 
                            if allelic_support_tags['control_af_tag'] != "_NA_":    
                                if {allelic_support_tags['control_af_tag']}.issubset(raw_maf_data.columns):
                                    if raw_maf_data[raw_maf_data[allelic_support_tags['control_af_tag']].isna()].empty:
                                        raw_maf_data['n_alt_count'] = None
                                        raw_maf_data.loc[:,"n_alt_count"] = \
                                            raw_maf_data.loc[:,allelic_support_tags['control_af_tag']].astype(float) * raw_maf_data.loc[:,"n_depth"].astype(int)
                                    
                                        raw_maf_data.loc[:,"n_alt_count"] = round(raw_maf_data.loc[:,"n_alt_count"].astype(float),0).astype(int)
                                        raw_maf_data['n_ref_count'] = None
                                        raw_maf_data.loc[:,"n_ref_count"] = \
                                            raw_maf_data.loc[:,"n_depth"].astype(int) - raw_maf_data.loc[:,"n_alt_count"]
        
        #raw_maf_data = raw_maf_data.fillna("")
        with open(maf_fname, 'w') as f:
            f.write(f'{header_line}\n')
        f.close()
        raw_maf_data.to_csv(maf_fname, sep="\t", index=False, mode='a')
        if not debug:
            remove_file(maf_tmp_fname)                     
                
                
        
            
    return 0


def construct_hgvsg(chr_val, start, end, ref, alt):
    """
    Construct HGVSg notation from variant coordinates and alleles.

    Args:
        chr_val: Chromosome (e.g., "1", "X", "MT")
        start: Start position (1-based)
        end: End position (1-based)
        ref: Reference allele
        alt: Alternate allele

    Returns:
        HGVSg notation string (e.g., "1:g.12345A>T") or empty string if invalid
    """
    # Return empty string if any required field is missing or NA
    if pd.isna(chr_val) or pd.isna(start) or pd.isna(end) or pd.isna(ref) or pd.isna(alt):
        return ""

    # Convert to strings to handle properly
    chr_val = str(chr_val)
    start_str = str(start)
    end_str = str(end)
    ref = str(ref)
    alt = str(alt)

    # Handle NaN strings
    if chr_val == 'nan' or start_str == 'nan' or end_str == 'nan' or ref == 'nan' or alt == 'nan':
        return ""

    # Handle empty strings
    if not chr_val or not start_str or not end_str or not ref or not alt:
        return ""

    # SNV: single position, reference != alt, both single base
    if start_str == end_str and len(ref) == 1 and len(alt) == 1:
        return f"{chr_val}:g.{start_str}{ref}>{alt}"

    # Deletion: reference longer than alt
    if len(ref) > len(alt):
        if len(alt) == 0 or alt == "-":
            # Pure deletion
            if start_str == end_str:
                return f"{chr_val}:g.{start_str}del"
            else:
                return f"{chr_val}:g.{start_str}_{end_str}del"
        else:
            # Deletion with replacement (delins)
            if start_str == end_str:
                return f"{chr_val}:g.{start_str}delins{alt}"
            else:
                return f"{chr_val}:g.{start_str}_{end_str}delins{alt}"

    # Insertion: alt longer than reference
    if len(alt) > len(ref):
        if len(ref) == 0 or ref == "-":
            # Pure insertion
            return f"{chr_val}:g.{start_str}_{end_str}ins{alt}"
        else:
            # Complex indel (delins)
            if start_str == end_str:
                return f"{chr_val}:g.{start_str}delins{alt}"
            else:
                return f"{chr_val}:g.{start_str}_{end_str}delins{alt}"

    # MNV or delins: same length but different
    if len(ref) == len(alt) and len(ref) > 1:
        if start_str == end_str:
            return f"{chr_val}:g.{start_str}delins{alt}"
        else:
            return f"{chr_val}:g.{start_str}_{end_str}delins{alt}"

    # Default case
    return ""


def add_var_id_to_vcf(input_vcf: str, output_vcf: str, logger=None) -> None:
    """
    Stamp a VAR_ID INFO field (CHROM_POS_REF_ALT, VCF-format) onto every record
    of the VEP-annotated VCF before vcf2maf conversion.

    Doing this while still in VCF space ensures that indel alleles retain their
    anchor base, making VAR_ID directly comparable with the CHROM_POS_REF_ALT key
    constructed later from the vcf2tsv output in variant.py.

    Args:
        input_vcf:  Path to VEP-annotated VCF (input)
        output_vcf: Path to VCF with VAR_ID INFO field added (output)
        logger:     Logger instance
    """
    import cyvcf2

    if logger is None:
        logger = logging.getLogger("pcgr-add-var-id")

    if not os.path.exists(input_vcf):
        logger.warning(f"VCF not found: {input_vcf} - skipping VAR_ID annotation")
        return

    vcf_in = cyvcf2.VCF(input_vcf)
    vcf_in.add_info_to_header({
        'ID':          'VAR_ID',
        'Number':      '1',
        'Type':        'String',
        'Description': 'Variant identifier CHROM_POS_REF_ALT (VCF-format alleles)'
    })

    writer = cyvcf2.Writer(output_vcf, vcf_in)
    n = 0
    for rec in vcf_in:
        alt = rec.ALT[0] if len(rec.ALT) == 1 else ','.join(rec.ALT)
        rec.INFO['VAR_ID'] = f'{rec.CHROM}_{rec.POS}_{rec.REF}_{alt}'
        writer.write_record(rec)
        n += 1
    writer.close()
    vcf_in.close()
    logger.info(f"Stamped VAR_ID on {n} variants -> {os.path.basename(output_vcf)}")


def generate_oncokb_maf_input(maf_fname: str,
                                output_fname: str,
                                oncokb_maf_query_all: bool = False,
                                logger: Optional[logging.Logger] = None) -> None:
    """
    Generate OncoKB-compatible MAF input file from full MAF file.

    OncoKB requires a minimal set of columns for variant annotation.
    This function extracts the necessary columns and removes duplicates.

    Args:
        maf_fname: Path to the full MAF file generated by vcf2maf
        output_fname: Path to output OncoKB-compatible MAF file
        oncokb_maf_query_all: If True, skip filtering of non-coding variant classes
        logger: Logger instance
    """
    if logger is None:
        logger = logging.getLogger("pcgr-oncokb-maf")

    # Check if MAF file exists
    if not os.path.exists(maf_fname):
        logger.warning(f"MAF file not found: {maf_fname} - skipping OncoKB MAF generation")
        return

    # Read MAF file (skip version header line)
    try:
        with open(maf_fname) as f:
            first_line = f.readline().strip()

        # Check if first line is version header
        skip_rows = 1 if first_line.startswith('#version') else 0

        maf_df = pd.read_csv(maf_fname, sep="\t", skiprows=skip_rows, dtype='string', na_values=['.'], low_memory=False)
    except Exception as e:
        logger.warning(f"Failed to read MAF file: {e} - skipping OncoKB MAF generation")
        return

    if maf_df.empty:
        logger.warning("Empty MAF file - skipping OncoKB MAF generation")
        return

    # Define OncoKB-required columns in order.
    # VAR_ID (CHROM_POS_REF_ALT in VCF format) is carried through so that the
    # OncoKB-annotated MAF can be joined back to the VCF-derived TSV in pcgrr.
    oncokb_columns = [
        'NCBI_Build',
        'Hugo_Symbol',
        'Variant_Classification',
        'Tumor_Sample_Barcode',
        'HGVSp_Short',
        'HGVSp',
        'HGVSg',
        'Chromosome',
        'Start_Position',
        'End_Position',
        'Reference_Allele',
        'Tumor_Seq_Allele1',
        'Tumor_Seq_Allele2',
        'VAR_ID'
    ]

    # Check which required columns are present
    available_columns = [col for col in oncokb_columns if col in maf_df.columns]

    if 'Hugo_Symbol' not in available_columns or 'Tumor_Sample_Barcode' not in available_columns:
        logger.warning("MAF file missing critical columns (Hugo_Symbol or Tumor_Sample_Barcode) - skipping OncoKB MAF generation")
        return

    # Extract only the OncoKB-required columns that are available
    oncokb_df = maf_df[available_columns].copy()

    # Filter out non-coding variants that are unlikely to be actionable
    # IGR = Intergenic Region, Intron, and flanking regions (3' and 5')
    # Skipped when oncokb_maf_query_all=True (intended for TARGETED/WES with non-coding regions of interest)
    if not oncokb_maf_query_all and 'Variant_Classification' in oncokb_df.columns:
        non_coding_classes = ['IGR', 'Intron', "3'Flank", "3'UTR", "5'UTR", 'RNA', 'lincRNA']
        before_filter = len(oncokb_df)
        oncokb_df = oncokb_df[~oncokb_df['Variant_Classification'].isin(non_coding_classes)]
        after_filter = len(oncokb_df)

        if before_filter > after_filter:
            filtered_count = before_filter - after_filter
            logger.info(f"Filtered out {filtered_count} non-coding variants (IGR, Intron, UTR, 3'Flank regions)")
    elif oncokb_maf_query_all:
        logger.info("oncokb_maf_query_all is enabled - skipping non-coding variant filter, all variant classes will be submitted to OncoKB")

    # Construct HGVSg notation if missing or empty
    # HGVSg is useful for OncoKB API calls using genomic coordinates
    if 'HGVSg' in oncokb_df.columns:
        # Check if HGVSg column is empty or has missing values
        needs_hgvsg = oncokb_df['HGVSg'].isna() | (oncokb_df['HGVSg'] == '') | (oncokb_df['HGVSg'] == '.')

        if needs_hgvsg.any():
            # Construct HGVSg for rows that need it
            required_cols = ['Chromosome', 'Start_Position', 'End_Position', 'Reference_Allele', 'Tumor_Seq_Allele2']
            if all(col in oncokb_df.columns for col in required_cols):
                logger.info(f"Constructing HGVSg notation for {needs_hgvsg.sum()} variants with missing HGVSg")
                oncokb_df.loc[needs_hgvsg, 'HGVSg'] = oncokb_df[needs_hgvsg].apply(
                    lambda row: construct_hgvsg(
                        row['Chromosome'],
                        row['Start_Position'],
                        row['End_Position'],
                        row['Reference_Allele'],
                        row['Tumor_Seq_Allele2']
                    ),
                    axis=1
                )
    elif all(col in oncokb_df.columns for col in ['Chromosome', 'Start_Position', 'End_Position', 'Reference_Allele', 'Tumor_Seq_Allele2']):
        # HGVSg column doesn't exist, create it
        logger.info("HGVSg column not found - constructing from genomic coordinates")
        oncokb_df['HGVSg'] = oncokb_df.apply(
            lambda row: construct_hgvsg(
                row['Chromosome'],
                row['Start_Position'],
                row['End_Position'],
                row['Reference_Allele'],
                row['Tumor_Seq_Allele2']
            ),
            axis=1
        )

    # Remove duplicate rows
    oncokb_df = oncokb_df.drop_duplicates()

    if len(oncokb_df) > 0:
        # Write as plain text TSV file
        oncokb_df.to_csv(output_fname, sep="\t", index=False)
        logger.info(f"Generated OncoKB MAF input file: {os.path.basename(output_fname)} ({len(oncokb_df)} variants)")
    else:
        logger.warning("No valid variants for OncoKB MAF input - skipping file generation")