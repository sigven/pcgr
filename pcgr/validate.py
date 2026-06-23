"""
Shared input validation utilities for PCGR and CPSR.

Contains:
- VCF processing helpers (sort/filter, multiallelic detection, INFO-tag collection)
- Molecular-data format checkers (CNA segments, RNA fusion, RNA expression,
  germline classification file, panel-of-normals VCF)
"""
import csv
import re
import gzip
import os

import pandas as pd
from cyvcf2 import VCF

from pcgr.annoutils import read_infotag_file, read_vcfanno_tag_file
from pcgr.utils import (
    error_message, warn_message,
    random_id_generator, remove_file, check_subprocess,
)
from pcgr import pcgr_vars


def build_chrom_filter():
    """Return a bcftools --regions string covering autosomal/sex chromosomes (with and without chr prefix)."""
    chroms = [str(x) for x in [*range(1, 23), 'X', 'Y']]
    return ','.join([*[f'chr{c}' for c in chroms], *chroms])


def detect_multiallelic_sites(vcf_obj):
    """
    Iterate a cyvcf2 VCF object and return a list of variant IDs for multiallelic sites.
    The caller must supply a freshly-opened cyvcf2.VCF object (the iterator is consumed).
    """
    multiallelic = []
    for rec in vcf_obj:
        if len(rec.ALT) > 1:
            pos = rec.start + 1
            alt = ','.join(str(a) for a in rec.ALT)
            multiallelic.append(f'{rec.CHROM}:{pos}_{rec.REF}->{alt}')
    return multiallelic


def bgzip_sort_filter_vcf(input_vcf, output_vcf, output_dir, sample_id, workflow, logger, debug):
    """
    Sort and chromosome-filter a VCF in three steps:
      1. bgzip + tabix the input
      2. bcftools sort into a compressed, indexed temp file
      3. bcftools view restricted to autosomal/sex chromosomes, strip the 'chr' prefix

    For somatic input (*workflow* = ``'pcgr'``) FORMAT metadata lines are also
    stripped and only the first 8 VCF columns are kept (genotype data is not
    needed downstream).  For germline input (*workflow* = ``'cpsr'``) the
    genotype columns are preserved.

    All intermediate files are created in *output_dir* and cleaned up unless
    *debug* is True.  The final uncompressed, filtered VCF is written to
    *output_vcf*; compression and indexing of that file are left to the caller.
    """
    if workflow not in ('pcgr', 'cpsr'):
        raise ValueError(f"workflow must be 'pcgr' or 'cpsr', got {workflow!r}")

    random_str = random_id_generator(10)
    prefix = f'{sample_id}.{workflow}_validate'

    vcf_bgzipped = os.path.join(output_dir, f'{prefix}.sort.{random_str}.vcf.gz')
    vcf_sorted   = os.path.join(output_dir, f'{prefix}.sort.{random_str}.sorted.vcf.gz')
    sort_log     = os.path.join(output_dir, f'{prefix}.sort.{random_str}.log')

    chrom_to_keep = build_chrom_filter()

    cmd_sort = (
        f'bcftools view {input_vcf} | bgzip -cf > {vcf_bgzipped} '
        f'&& tabix -p vcf {vcf_bgzipped} '
        f'&& bcftools sort --temp-dir {output_dir} -Oz {vcf_bgzipped} > {vcf_sorted} '
        f'2> {sort_log} && tabix -p vcf {vcf_sorted}'
    )

    if workflow == 'pcgr':
        # Somatic: strip FORMAT metadata and drop genotype columns (cols 1-8 only)
        # Note: M/MT variants are skipped — requires additional VEP cache handling
        cmd_filter = (
            f"bcftools view --regions {chrom_to_keep} {vcf_sorted} "
            f"| egrep -v '^##FORMAT=' | cut -f1-8 | sed 's/^chr//' > {output_vcf}"
        )
    else:
        # Germline: keep genotype columns intact
        # Note: M/MT variants are skipped — requires additional VEP cache handling
        cmd_filter = (
            f"bcftools view --regions {chrom_to_keep} {vcf_sorted} "
            f"| sed 's/^chr//' > {output_vcf}"
        )

    logger.info('Extracting variants on autosomal/sex chromosomes only (1-22,X,Y) with bcftools')
    check_subprocess(logger, cmd_sort, debug)
    check_subprocess(logger, cmd_filter, debug)

    if not debug:
        for fn in [vcf_bgzipped, vcf_sorted, sort_log,
                   vcf_bgzipped + '.tbi', vcf_sorted + '.tbi']:
            remove_file(fn)


def collect_workflow_info_tags(refdata_assembly_dir, workflow, logger):
    """
    Build and return the merged dict of all VCF INFO tags that will be appended by
    PCGR or CPSR annotation.  *workflow* must be one of ``'pcgr'`` or ``'cpsr'``.
    """
    if workflow not in ('pcgr', 'cpsr'):
        raise ValueError(f"workflow must be 'pcgr' or 'cpsr', got {workflow!r}")

    tags = read_infotag_file(
        os.path.join(refdata_assembly_dir, 'vcf_infotags_other.tsv'), scope=workflow)
    tags.update(read_infotag_file(
        os.path.join(refdata_assembly_dir, 'vcf_infotags_vep.tsv'), scope='vep'))

    # vcfanno variant tracks — CPSR includes two extra germline-specific tracks
    variant_tracks = ['clinvar', 'tcga', 'gwas', 'dbnsfp']
    if workflow == 'cpsr':
        variant_tracks += ['dbmts', 'gnomad_non_cancer']

    bed_tracks = ['simplerepeat', 'winmsk', 'rmsk', 'gerp']

    tag_fnames = {}
    for t in variant_tracks:
        tag_fnames[t] = os.path.join(
            refdata_assembly_dir, 'variant', 'vcf', t, f'{t}.vcfanno.vcf_info_tags.txt')
    for t in bed_tracks:
        tag_fnames[t] = os.path.join(
            refdata_assembly_dir, 'misc', 'bed', t, f'{t}.vcfanno.vcf_info_tags.txt')
    tag_fnames['gene_transcript_xref'] = os.path.join(
        refdata_assembly_dir, 'gene', 'bed', 'gene_transcript_xref',
        'gene_transcript_xref.vcfanno.vcf_info_tags.txt')

    for fname in tag_fnames.values():
        tags.update(read_vcfanno_tag_file(fname, logger))

    return tags


# ---------------------------------------------------------------------------
# Molecular data format validators
# ---------------------------------------------------------------------------

def is_valid_cna(input_cna_segment_fname, logger):
    """
    Check whether the CNA segment file (tab-separated) has the required columns,
    correct data types, and valid coordinate ranges.
    """
    cna_reader = csv.DictReader(open(input_cna_segment_fname, 'r'), delimiter='\t')

    required_columns = ['Chromosome', 'Start', 'End', 'nMajor', 'nMinor']
    fieldnames = cna_reader.fieldnames or []
    missing = set(required_columns) - set(fieldnames)
    if missing:
        err_msg = (
            f"Copy number segment file ({input_cna_segment_fname}) is missing one or more "
            f"required column(s): {', '.join(sorted(missing))}\n"
            f"Column names present in file: {fieldnames}"
        )
        return error_message(err_msg, logger)

    cna_dataframe = pd.read_csv(input_cna_segment_fname, sep='\t')
    if cna_dataframe.empty:
        err_msg = 'Copy number segment file is empty - contains NO segments'
        return error_message(err_msg, logger)

    for elem in ['Start', 'End']:
        if cna_dataframe[elem].dtype.kind not in 'i':
            err_msg = f'Copy number segment file contains non-integer values for column: "{elem}"'
            return error_message(err_msg, logger)

    for elem in ['nMajor', 'nMinor']:
        if cna_dataframe[elem].dtype.kind not in 'if':
            err_msg = f'Copy number segment file contains non-float/integer values for column: "{elem}"'
            return error_message(err_msg, logger)

    for rec in cna_reader:
        if int(rec['End']) < int(rec['Start']):
            err_msg = (
                f"Detected wrongly formatted chromosomal segment - 'Start' is greater than 'End' "
                f"({rec['Chromosome']}:{rec['Start']}-{rec['End']})"
            )
            return error_message(err_msg, logger)
        if int(rec['End']) < 1 or int(rec['Start']) < 1:
            err_msg = (
                f"Detected wrongly formatted chromosomal segment - 'Start' or 'End' is less than "
                f"or equal to zero ({rec['Chromosome']}:{rec['Start']}-{rec['End']})"
            )
            return error_message(err_msg, logger)

    logger.info(f"Copy number segment file ('{os.path.basename(input_cna_segment_fname)}') adheres to the correct format")
    return 0


def is_valid_rna_fusion(rna_fusion_file, logger):
    """
    Check whether the RNA fusion transcript file has the required columns,
    correct data types, and valid fusion gene / breakpoint formats.
    """
    rna_fusion_dataframe = pd.read_csv(rna_fusion_file, sep='\t')

    fieldnames = list(rna_fusion_dataframe.columns)
    required_cols: set[str] = {'FusionGene', 'LeftBreakpoint', 'RightBreakpoint', 'SplitReads'}
    missing = set(required_cols) - set(fieldnames)
    if missing:
        err_msg = (
            f"RNA fusion file ({rna_fusion_file}) is missing required column(s): "
            f"{', '.join(sorted(missing))}\nColumn names present in file: {fieldnames}"
        )
        return error_message(err_msg, logger)
    if rna_fusion_dataframe.empty:
        return error_message('RNA fusion file is empty - contains NO fusions', logger)
    if rna_fusion_dataframe['FusionGene'].dtype.kind not in 'O':
        return error_message(f"'FusionGene' column cannot be of type '{rna_fusion_dataframe['FusionGene'].dtype}'", logger)
    if rna_fusion_dataframe['LeftBreakpoint'].dtype.kind not in 'O':
        return error_message(f"'LeftBreakpoint' column cannot be of type '{rna_fusion_dataframe['LeftBreakpoint'].dtype}'", logger)
    if rna_fusion_dataframe['RightBreakpoint'].dtype.kind not in 'O':
        return error_message(f"'RightBreakpoint' column cannot be of type '{rna_fusion_dataframe['RightBreakpoint'].dtype}'", logger)
    if rna_fusion_dataframe['SplitReads'].dtype.kind not in 'i':
        return error_message(f"'SplitReads' column cannot be of type '{rna_fusion_dataframe['SplitReads'].dtype}'", logger)
    if rna_fusion_dataframe['SplitReads'].min() < 0:
        return error_message("'SplitReads' column cannot contain negative values", logger)
    if 'Score' in rna_fusion_dataframe.columns:
        if rna_fusion_dataframe['Score'].dtype.kind not in 'f':
            return error_message(f"Optional 'Score' column cannot be of type '{rna_fusion_dataframe['Score'].dtype}'", logger)

    gene_pattern = re.compile(r'^[A-Za-z0-9][A-Za-z0-9._-]*$')
    # Accept both bare (2:42264731) and chr-prefixed (chr2:42264731) formats
    bp_pattern   = re.compile(r'^(chr)?([1-9]|1[0-9]|2[0-2]|X|Y|MT):\d+$')
    empty_partner_fusions = []

    for _, rec in rna_fusion_dataframe.iterrows():
        fusion_gene = rec['FusionGene']
        left_bp     = rec['LeftBreakpoint']
        right_bp    = rec['RightBreakpoint']
        split_reads = rec['SplitReads']

        if '--' not in fusion_gene and '::' not in fusion_gene:
            return error_message(
                f"RNA fusion gene '{fusion_gene}' missing mandatory '--' or '::' separator", logger)

        genes = fusion_gene.split('--') if '--' in fusion_gene else fusion_gene.split('::')
        if len(genes) != 2:
            return error_message(
                f"RNA fusion gene '{fusion_gene}' does not have exactly two parts", logger)

        for g in genes:
            if g in ('', 'v'):
                empty_partner_fusions.append(fusion_gene)
            else:
                if not gene_pattern.match(g):
                    return error_message(f"Invalid gene symbol '{g}' in fusion '{fusion_gene}'", logger)
                if g.isdigit():
                    return error_message(f"Gene name '{g}' cannot be all digits", logger)

        if not bp_pattern.match(str(left_bp)):
            return error_message(
                f"Invalid LeftBreakpoint: '{left_bp}' - expected format <chrom>:<position> "
                f"(e.g. '2:42264731' or 'chr2:42264731')", logger)
        if not bp_pattern.match(str(right_bp)):
            return error_message(
                f"Invalid RightBreakpoint: '{right_bp}' - expected format <chrom>:<position> "
                f"(e.g. '2:29223528' or 'chr2:29223528')", logger)
        if int(split_reads) < 0:
            return error_message(f"SplitReads cannot be negative - found '{split_reads}'", logger)

    if empty_partner_fusions:
        n = len(empty_partner_fusions)
        examples = ', '.join(f"'{f}'" for f in empty_partner_fusions[:2])
        suffix = f" (and {n - 2} more)" if n > 2 else ""
        warn_message(f"Found {n} fusion gene(s) with empty/missing partner: {examples}{suffix}", logger)

    logger.info(f"RNA fusion file ('{os.path.basename(rna_fusion_file)}') adheres to the correct format")
    return 0


def is_valid_rna_expression(rna_exp_file, logger):
    """
    Check whether the bulk-RNA expression file has the required columns and valid values.
    """
    rna_exp_reader = csv.DictReader(open(rna_exp_file, 'r'), delimiter='\t')
    fieldnames = rna_exp_reader.fieldnames or []
    if 'TargetID' not in fieldnames or 'TPM' not in fieldnames:
        err_msg = (
            f"Bulk-RNA expression file ({rna_exp_file}) is missing required column(s): "
            f"'TargetID', 'TPM'\nColumn names present in file: {fieldnames}"
        )
        return error_message(err_msg, logger)

    rna_exp_dataframe = pd.read_csv(rna_exp_file, sep='\t')
    if rna_exp_dataframe.empty:
        return error_message('RNA gene expression file is empty - contains NO gene expression estimates', logger)
    if rna_exp_dataframe['TargetID'].dtype.kind not in 'O':
        return error_message(f"'TargetID' column cannot be of type '{rna_exp_dataframe['TargetID'].dtype}'", logger)
    if rna_exp_dataframe['TPM'].dtype.kind not in 'if':
        return error_message(f"'TPM' column cannot be of type '{rna_exp_dataframe['TPM'].dtype}'", logger)

    for rec in rna_exp_reader:
        if not (float(rec['TPM']) >= 0):
            return error_message(f"'TPM' column cannot contain negative values - value was {rec['TPM']}", logger)

    logger.info(f"RNA expression file ('{os.path.basename(rna_exp_file)}') adheres to the correct format")
    return 0


def is_valid_germline(germline_file, build, logger):
    """
    Check whether the CPSR-classified germline variants file has the expected naming
    convention and required columns.
    """
    if not os.path.isfile(germline_file):
        return error_message(f"Germline variants file ({germline_file}) does not exist", logger)

    if not str(germline_file).endswith(f'cpsr.{build}.classification.tsv.gz'):
        return error_message(
            f"Germline variants file ({germline_file}) does not adhere to the correct naming "
            f"format - wrong build or file type", logger)

    with gzip.open(germline_file, 'rt') as f:
        germline_reader = csv.DictReader(f, delimiter='\t')
        fieldnames = germline_reader.fieldnames or []
        for col in pcgr_vars.germline_input_required_cols:
            if col not in fieldnames:
                return error_message(
                    f"Germline variants file ({germline_file}) is missing required column: {col}", logger)

    logger.info(f"Germline variants file ('{os.path.basename(germline_file)}') adheres to the correct format")
    return 0


def is_valid_pon_vcf(panel_normal_vcf_fname, logger):
    """
    Check that the panel-of-normals VCF carries the required 'PANEL_OF_NORMALS' INFO flag.
    """
    vcf_obj = VCF(panel_normal_vcf_fname)
    for e in vcf_obj.header_iter():
        h = e.info()
        if (h.get('HeaderType') == 'INFO' and h.get('Type') == 'Flag'
                and h.get('ID') == 'PANEL_OF_NORMALS'):
            logger.info("Found 'PANEL_OF_NORMALS' INFO flag in panel-of-normals VCF header")
            return 1

    return error_message("INFO flag 'PANEL_OF_NORMALS' is missing from the panel-of-normals VCF header", logger)
