#!/usr/bin/env python
"""
Genome-wide scores from allele-specific copy number segments:

  - Fraction of genome altered (FGA; cBioPortal-style)

and genomic instability ("HRD scar") scores:

  - HRD-LOH  (Abkevich et al., Br J Cancer 2012)
  - LST      (Large-scale state transitions; Popova et al., Cancer Res 2012)
  - TAI      (Telomeric allelic imbalance; Birkbak et al., Cancer Discov 2012)

This is a Python port (translation) of the algorithms as implemented in the
scarHRD R package (https://github.com/sztup/scarHRD, calc.hrd.R/calc.lst.R/
calc.ai_new.R, shrink.seg.ai.R, preprocess.hrd.R), so that the scores can be
computed without installing scarHRD (which hard-depends on the 'sequenza' R
package, only distributed via Bitbucket). scarHRD is Copyright (c) 2020 Zsofia
Sztupinszki, and is distributed under the MIT License - see the notice below
this docstring, and THIRD_PARTY_LICENSES.md in the root of the PCGR repository.

Reference: Sztupinszki et al., "Migrating the SNP array-based homologous
recombination deficiency measures to next generation sequencing data of
breast cancer", npj Breast Cancer 2018 (https://doi.org/10.1038/s41523-018-0066-6)

Fidelity to scarHRD: the three scores reproduce scarHRD (v0.1.1) exactly on
the scarHRD example file (examples/test2.txt: HRD-LOH 25, TAI 35, LST 33,
sum 93) and on 160 randomly sampled TCGA ASCAT3 profiles (BRCA/OV/PRAD/PAAD;
160/160 identical for all three scores when scarHRD's centromere definition
is used). Two details matter for exact agreement:
  - scarHRD defines the centromere of each chromosome as the span of the
    cytoband 'acen' (centromeric) bands (first start to last end). These are
    wider than, and differ from, the centromere-model coordinates in PCGR's
    chromsize file; using the latter changes LST in ~15% of samples (by up to 2)
    and TAI in ~35% (by up to 4). The scar scores therefore use the 'acen'
    bands, read from the cytoband file in the PCGR data bundle (see
    scarhrd_chrom_arms()); the chromsize file is used for chromosome lengths
    only (FGA/WGD).
  - TAI follows calc.ai_new() including its minimum segment size (1 Mb) and
    its per-chromosome "local ploidy" rule for deciding whether a segment is
    balanced (see calc_tai()).

All three scores are computed on autosomes only (1-22), matching scarHRD,
which drops sex chromosomes entirely before scoring.

Segments flagged in the optional input column 'nMajor_nMinor_Imputed' (total
copy number observed, major/minor split imputed) have an unknown allele-
specific state. They are not dropped for HRD-LOH/LST/TAI (dropping them would
make their neighbours adjacent, creating or hiding state transitions), but
treated as unknown: they never merge with neighbouring segments, are never
counted as LOH/imbalanced, and transitions into or out of them are not
counted. Short imputed segments are still removed by the (state-independent)
minimum segment size filters of LST/TAI, as any other segment. Imputed
segments are excluded from the segment coverage check and the genome doubling
call, but included in the fraction of genome altered (total copy number).
"""

## ----------------------------------------------------------------------------
## Third-party notice: the algorithms implemented in this module (HRD-LOH, LST,
## TAI and the segment shrinking/preprocessing steps) are a Python port of the
## scarHRD R package (https://github.com/sztup/scarHRD), which is licensed as
## follows:
##
## MIT License
##
## Copyright (c) 2020 Zsofia Sztupinszki
##
## Permission is hereby granted, free of charge, to any person obtaining a copy
## of this software and associated documentation files (the "Software"), to deal
## in the Software without restriction, including without limitation the rights
## to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
## copies of the Software, and to permit persons to whom the Software is
## furnished to do so, subject to the following conditions:
##
## The above copyright notice and this permission notice shall be included in all
## copies or substantial portions of the Software.
##
## THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
## IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
## FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
## AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
## LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
## OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
## SOFTWARE.
## ----------------------------------------------------------------------------


import logging
import os
from typing import Optional

import pandas as pd

from pcgr.utils import check_file_exists, error_message, get_cna_imputed_flag

AUTOSOMES = [str(i) for i in range(1, 23)]

HRD_LOH_SIZE_LIMIT_MB = 15.0
LST_SEG_SIZE_LIMIT_MB = 10.0
LST_GAP_LIMIT_MB = 3.0
LST_MIN_ARM_SEGMENT_MB = 3.0

MYRIAD_HRD_POSITIVE_THRESHOLD = 42

## Minimum fraction of the autosomal genome that the copy number segments must
## span for the HRD scar scores to be meaningful. Genome-wide profiles (WGS, WES
## with genome-wide segmentation) span ~90% or more; profiles of targeted gene
## panels span only a tiny fraction.
HRD_MIN_AUTOSOMAL_COVERAGE = 0.5

def read_cytoband_centromeres(cytoband_fname: str, logger: Optional[logging.Logger] = None) -> dict:
    """
    Read centromere coordinates from the cytoband file in the PCGR data bundle
    ('<refdata_assembly_dir>/misc/tsv/cytoband/cytoband.tsv.gz'). As in scarHRD,
    the centromere of a chromosome is the span of its cytoband 'acen'
    (centromeric) bands, i.e. the smallest start and largest end among these.

    Returns a dict keyed by chromosome name (without 'chr' prefix), with
    (centromere_start, centromere_end) tuples, e.g. {'1': (121700000, 125100000), ...}
    """
    if logger is None:
        logger = logging.getLogger("pcgr-hrd")
    check_file_exists(cytoband_fname, logger=logger)

    cytoband = pd.read_csv(cytoband_fname, sep="\t", usecols=['chrom', 'start', 'end', 'gieStain'])
    acen = cytoband[cytoband['gieStain'] == 'acen'].copy()
    acen['chrom'] = acen['chrom'].astype(str).str.replace("^(C|c)hr", "", regex=True)
    spans = acen.groupby('chrom').agg(start=('start', 'min'), end=('end', 'max'))
    return {chrom: (int(row['start']), int(row['end'])) for chrom, row in spans.iterrows()}


def scarhrd_chrom_arms(chrom_arms: dict, cytoband_fname: str, logger: Optional[logging.Logger] = None) -> dict:
    """
    Return a copy of 'chrom_arms' (see read_chrom_arms) in which the centromere
    coordinates are replaced by the cytoband 'acen' spans used by scarHRD (see
    read_cytoband_centromeres); chromosome lengths are unchanged.
    """
    centromeres = read_cytoband_centromeres(cytoband_fname, logger=logger)
    arms = {c: dict(v) for c, v in chrom_arms.items()}
    for chrom, (cstart, cend) in centromeres.items():
        if chrom in arms:
            arms[chrom]['centromere_start'] = cstart
            arms[chrom]['centromere_end'] = cend
    return arms


def read_chrom_arms(chromsizes_fname: str, logger: Optional[logging.Logger] = None) -> dict:
    """
    Read PCGR's '<refdata_assembly_dir>/chromsize.<build>.tsv' reference file
    (columns: chrom, length, centromere_left, centromere_right) into a dict
    keyed by chromosome name (without 'chr' prefix), e.g.:
        {'1': {'length': 248956422, 'centromere_start': 122026460, 'centromere_end': 125184587}, ...}
    """
    if logger is None:
        logger = logging.getLogger("pcgr-hrd")
    check_file_exists(chromsizes_fname, logger=logger)

    df = pd.read_csv(
        chromsizes_fname, sep="\t", skiprows=1,
        names=['Chromosome', 'ChromLength', 'CentromereStart', 'CentromereEnd'])
    df['Chromosome'] = df['Chromosome'].astype(str).str.replace("^(C|c)hr", "", regex=True)

    chrom_arms = {}
    for _, row in df.iterrows():
        chrom_arms[row['Chromosome']] = {
            'length': int(row['ChromLength']),
            'centromere_start': int(row['CentromereStart']),
            'centromere_end': int(row['CentromereEnd']),
        }
    return chrom_arms


def _prepare_segments(cna_df: pd.DataFrame, logger: logging.Logger, exclude_imputed: bool = False) -> pd.DataFrame:
    """
    Common preprocessing shared by all three scores:
      - restrict to autosomes 1-22 (sex chromosomes excluded, as in scarHRD)
      - ensure nMajor >= nMinor per segment (scarHRD swaps if not)
      - sort by chromosome/start
      - add boolean column 'Imputed' from the optional input column
        'nMajor_nMinor_Imputed' (segments are dropped if 'exclude_imputed')
    Expects columns: Chromosome, Start, End, nMajor, nMinor (same schema as
    used elsewhere in pcgr/cna.py for the user-supplied CNA segment file).
    """
    required = {'Chromosome', 'Start', 'End', 'nMajor', 'nMinor'}
    if not required.issubset(cna_df.columns):
        error_message(
            f"CNA segment data frame passed to pcgr.hrd is missing required columns: "
            f"{required - set(cna_df.columns)}", logger)

    seg = cna_df[['Chromosome', 'Start', 'End', 'nMajor', 'nMinor']].copy()
    seg['Imputed'] = get_cna_imputed_flag(cna_df)
    if exclude_imputed:
        seg = seg[~seg['Imputed']].copy()
    seg['Chromosome'] = seg['Chromosome'].astype(str).str.replace("^(C|c)hr", "", regex=True)
    seg = seg[seg['Chromosome'].isin(AUTOSOMES)].copy()

    seg['Start'] = seg['Start'].astype(int)
    seg['End'] = seg['End'].astype(int)
    seg['nMajor'] = seg['nMajor'].round(0).astype(int)
    seg['nMinor'] = seg['nMinor'].round(0).astype(int)

    ## scarHRD (preprocess.hrd.R) swaps alleles if the "minor" column exceeds
    ## the "major" column, to guarantee nMajor >= nMinor throughout
    swap_mask = seg['nMinor'] > seg['nMajor']
    seg.loc[swap_mask, ['nMajor', 'nMinor']] = seg.loc[swap_mask, ['nMinor', 'nMajor']].values

    seg = seg.sort_values(['Chromosome', 'Start']).reset_index(drop=True)
    return seg


def _shrink(seg: pd.DataFrame) -> pd.DataFrame:
    """
    Port of scarHRD's shrink.seg.ai(): merge consecutive segments (assumed
    already sorted by Start, single chromosome/arm) that share an identical
    (nMajor, nMinor) allele-specific copy number state into a single segment
    spanning from the first segment's Start to the last segment's End.
    Imputed segments (unknown allele-specific state) never merge with
    non-imputed neighbours; consecutive imputed segments merge with each other.
    """
    if len(seg) <= 1:
        return seg.reset_index(drop=True)

    seg = seg.reset_index(drop=True)
    group_id = [0]
    for i in range(1, len(seg)):
        if seg.loc[i, 'Imputed'] or seg.loc[i - 1, 'Imputed']:
            same_state = bool(seg.loc[i, 'Imputed'] and seg.loc[i - 1, 'Imputed'])
        else:
            same_state = (seg.loc[i, 'nMajor'] == seg.loc[i - 1, 'nMajor'] and
                          seg.loc[i, 'nMinor'] == seg.loc[i - 1, 'nMinor'])
        group_id.append(group_id[-1] if same_state else group_id[-1] + 1)
    seg['_group'] = group_id

    merged = seg.groupby('_group', as_index=False).agg(
        Chromosome=('Chromosome', 'first'),
        Start=('Start', 'min'),
        End=('End', 'max'),
        nMajor=('nMajor', 'first'),
        nMinor=('nMinor', 'first'),
        Imputed=('Imputed', 'first'))
    return merged.drop(columns=[c for c in merged.columns if c == '_group'], errors='ignore')


def calc_autosomal_segment_coverage(
        cna_df: pd.DataFrame,
        chrom_arms: dict,
        logger: Optional[logging.Logger] = None) -> float:
    """
    Fraction of the autosomal genome (chromosomes 1-22) that is spanned by the
    copy number segments with an observed (non-imputed) allele-specific copy
    number (overlapping segments are counted once, segments are clipped to
    chromosome lengths).
    """
    if logger is None:
        logger = logging.getLogger("pcgr-hrd")
    seg = _prepare_segments(cna_df, logger, exclude_imputed=True)
    if len(seg) == 0:
        return 0.0
    seg = _clip_segments_to_chromosomes(seg, chrom_arms)

    genome_length = sum(
        chrom_arms[c]['length'] for c in AUTOSOMES if c in chrom_arms)
    if genome_length <= 0:
        return 0.0

    covered = 0
    for chrom, chrom_seg in seg.groupby('Chromosome'):
        current_start, current_end = None, None
        for start, end in sorted(zip(chrom_seg['Start'], chrom_seg['End'])):
            if current_end is None or start > current_end:
                if current_end is not None:
                    covered += current_end - current_start
                current_start, current_end = start, end
            else:
                current_end = max(current_end, end)
        if current_end is not None:
            covered += current_end - current_start
    return round(float(covered / genome_length), 4)


def check_segment_coverage(
        cna_df: pd.DataFrame,
        chrom_arms: dict,
        min_coverage: float = HRD_MIN_AUTOSOMAL_COVERAGE,
        logger: Optional[logging.Logger] = None) -> float:
    """
    Check that the copy number segments span (at least 'min_coverage' of) the
    autosomal genome, as required for the HRD scar scores. Raises a ValueError
    if not, e.g. for copy number segments from targeted gene panels.

    Returns the autosomal segment coverage (fraction).
    """
    coverage = calc_autosomal_segment_coverage(cna_df, chrom_arms, logger=logger)
    if coverage < min_coverage:
        raise ValueError(
            f"The copy number segments (with observed, non-imputed allele-specific copy number) "
            f"span only {coverage * 100:.1f}% of the autosomal genome "
            f"(minimum required for genomic instability scoring: {min_coverage * 100:.0f}%) - "
            f"the HRD scar scores (HRD-LOH, LST, TAI) require genome-wide copy number segments "
            f"(e.g. from WGS/WES), and are not meaningful for targeted sequencing")
    return coverage


def calc_hrd_loh(
        cna_df: pd.DataFrame,
        size_limit_mb: float = HRD_LOH_SIZE_LIMIT_MB,
        logger: Optional[logging.Logger] = None) -> int:
    """
    HRD-LOH (Abkevich et al. 2012): the number of LOH regions exceeding
    'size_limit_mb' (default 15 Mb) that do NOT cover an entire chromosome.

    A segment is "LOH" when the minor allele copy number is 0 and the
    (capped) major allele copy number is > 0 - i.e. one parental haplotype
    is entirely lost, but not both (that would be a homozygous deletion,
    which is excluded here, consistent with scarHRD).
    """
    if logger is None:
        logger = logging.getLogger("pcgr-hrd")
    seg = _prepare_segments(cna_df, logger)
    if len(seg) == 0:
        return 0

    size_limit_bp = size_limit_mb * 1e6
    n_loh_regions = 0

    for chrom, chrom_seg in seg.groupby('Chromosome'):
        chrom_seg = chrom_seg.sort_values('Start').reset_index(drop=True)

        ## whole-chromosome LOH (every segment on the chromosome is LOH) is
        ## excluded entirely - attributable to chromosome-level aneuploidy,
        ## not a focal HR-deficiency "scar"
        ## (judged on segments with an observed allele-specific state only)
        known_seg = chrom_seg[~chrom_seg['Imputed']]
        whole_chrom_loh = len(known_seg) > 0 and bool((known_seg['nMinor'] == 0).all())
        if whole_chrom_loh:
            continue

        ## cap nMajor at 1: LOH state is treated identically regardless of
        ## how amplified the retained haplotype is, so that e.g. (5,0) and
        ## (1,0) segments merge together as the same "LOH" state
        chrom_seg = chrom_seg.copy()
        chrom_seg['nMajor'] = chrom_seg['nMajor'].clip(upper=1)

        chrom_seg = _shrink(chrom_seg)

        loh = chrom_seg[(chrom_seg['nMinor'] == 0) & (chrom_seg['nMajor'] != 0) & ~chrom_seg['Imputed']]
        loh = loh[(loh['End'] - loh['Start']) > size_limit_bp]
        n_loh_regions += len(loh)

    return n_loh_regions


def _clip_and_split_arms(chrom_seg: pd.DataFrame, centromere_start: int, centromere_end: int):
    """
    Split a chromosome's (already shrunk-eligible) segments into p-arm and
    q-arm segment lists, per scarHRD's calc.lst(): a segment is assigned to
    an arm if it reaches into that arm's side of the centromere (so a
    segment spanning the centromere is assigned to BOTH lists), each list is
    then independently shrunk, and the single boundary segment nearest the
    centromere in each arm is clipped exactly to the centromere coordinate.
    """
    p_arm = chrom_seg[chrom_seg['Start'] <= centromere_start].copy()
    q_arm = chrom_seg[chrom_seg['End'] >= centromere_end].copy()

    q_arm = _shrink(q_arm)
    if len(q_arm) > 0:
        q_arm.loc[q_arm.index[0], 'Start'] = centromere_end

    if len(p_arm) > 0:
        p_arm = _shrink(p_arm)
        p_arm.loc[p_arm.index[-1], 'End'] = centromere_start

    return p_arm, q_arm


def _drop_short_segments_iteratively(arm_seg: pd.DataFrame, min_len_bp: float) -> pd.DataFrame:
    """
    Port of scarHRD's iterative short-segment removal: repeatedly drop the
    first segment shorter than 'min_len_bp' and re-merge (shrink) neighbours
    that become adjacent as a result, until none remain.
    """
    arm_seg = arm_seg.reset_index(drop=True)
    while len(arm_seg) > 0:
        lengths = arm_seg['End'] - arm_seg['Start']
        short_idx = lengths[lengths < min_len_bp].index
        if len(short_idx) == 0:
            break
        arm_seg = arm_seg.drop(index=short_idx[0]).reset_index(drop=True)
        arm_seg = _shrink(arm_seg)
    return arm_seg


def _count_lst_transitions(arm_seg: pd.DataFrame, seg_size_limit_bp: float, gap_limit_bp: float) -> int:
    if len(arm_seg) < 2:
        return 0
    arm_seg = arm_seg.reset_index(drop=True)
    ## transitions into/out of imputed segments (unknown state) are not counted
    is_long = ((arm_seg['End'] - arm_seg['Start']) >= seg_size_limit_bp) & ~arm_seg['Imputed']
    n_lst = 0
    for k in range(1, len(arm_seg)):
        if is_long.iloc[k] and is_long.iloc[k - 1]:
            gap = arm_seg.loc[k, 'Start'] - arm_seg.loc[k - 1, 'End']
            if gap < gap_limit_bp:
                n_lst += 1
    return n_lst


def calc_lst(
        cna_df: pd.DataFrame,
        chrom_arms: dict,
        seg_size_limit_mb: float = LST_SEG_SIZE_LIMIT_MB,
        gap_limit_mb: float = LST_GAP_LIMIT_MB,
        min_arm_segment_mb: float = LST_MIN_ARM_SEGMENT_MB,
        logger: Optional[logging.Logger] = None) -> int:
    """
    LST (Popova et al. 2012): the number of chromosomal breaks between
    adjacent regions of at least 'seg_size_limit_mb' (default 10 Mb), with a
    gap between them no larger than 'gap_limit_mb' (default 3 Mb). Computed
    per chromosome arm (breaks spanning the centromere are never counted),
    after first removing segments shorter than 'min_arm_segment_mb'
    (default 3 Mb) as noise.
    """
    if logger is None:
        logger = logging.getLogger("pcgr-hrd")
    seg = _prepare_segments(cna_df, logger)
    if len(seg) == 0:
        return 0

    seg_size_limit_bp = seg_size_limit_mb * 1e6
    gap_limit_bp = gap_limit_mb * 1e6
    min_arm_segment_bp = min_arm_segment_mb * 1e6

    n_lst = 0
    for chrom, chrom_seg in seg.groupby('Chromosome'):
        if chrom not in chrom_arms:
            logger.warning(f"pcgr-hrd: no centromere info for chromosome '{chrom}' - skipping for LST")
            continue
        chrom_seg = chrom_seg.sort_values('Start').reset_index(drop=True)
        if len(chrom_seg) < 2:
            continue

        centromere_start = chrom_arms[chrom]['centromere_start']
        centromere_end = chrom_arms[chrom]['centromere_end']
        p_arm, q_arm = _clip_and_split_arms(chrom_seg, centromere_start, centromere_end)

        for arm_seg in (p_arm, q_arm):
            arm_seg = _drop_short_segments_iteratively(arm_seg, min_arm_segment_bp)
            n_lst += _count_lst_transitions(arm_seg, seg_size_limit_bp, gap_limit_bp)

    return n_lst


TAI_MIN_SEGMENT_MB = 1.0


def calc_tai(
        cna_df: pd.DataFrame,
        chrom_arms: dict,
        min_segment_mb: float = TAI_MIN_SEGMENT_MB,
        logger: Optional[logging.Logger] = None) -> int:
    """
    Telomeric allelic imbalance (Birkbak et al. 2012): the number of regions with
    allelic imbalance that extend to a chromosome's telomeric end without
    crossing the centromere. Port of scarHRD's calc.ai_new():

      - segments shorter than 'min_segment_mb' (1 Mb) are dropped, and adjacent
        segments with identical allele-specific copy number are then merged
      - a segment is "balanced" depending on the chromosome's local ploidy,
        defined as scarHRD does: the most common (by length) non-zero value of
        the minor allele copy number on the chromosome. For ploidy 1 or even
        (the common case), a segment is balanced when nMajor == nMinor; for odd
        ploidy >= 3, when nMajor + nMinor == ploidy and nMinor != 0
      - the first/last segment of a chromosome counts as telomeric AI if it is
        imbalanced, the chromosome has more than one segment, and the segment
        does not cross the centromere; a single imbalanced segment spanning the
        whole chromosome is "whole-chromosome AI" and is not counted

    Minor details: ties in the ploidy vote and an all-LOH sample (no non-zero
    minor copy number; scarHRD would fail) are resolved arbitrarily / by plain
    allelic imbalance, respectively. For scarHRD-identical results 'chrom_arms'
    should carry scarHRD's centromere coordinates (see scarhrd_chrom_arms()).
    """
    if logger is None:
        logger = logging.getLogger("pcgr-hrd")
    seg = _prepare_segments(cna_df, logger)
    if len(seg) == 0:
        return 0

    parts = []
    for chrom, chrom_seg in seg.groupby('Chromosome', sort=False):
        chrom_seg = chrom_seg.sort_values('Start').reset_index(drop=True)
        chrom_seg = _shrink(chrom_seg)
        chrom_seg = chrom_seg[(chrom_seg['End'] - chrom_seg['Start']) >= min_segment_mb * 1e6]
        if len(chrom_seg) > 0:
            parts.append(_shrink(chrom_seg.reset_index(drop=True)))
    if len(parts) == 0:
        return 0
    seg = pd.concat(parts, ignore_index=True)

    ## distinct non-zero minor allele copy numbers, in order of appearance
    ## (imputed segments do not take part in the local ploidy vote)
    minor_cn_values = [v for v in dict.fromkeys(seg.loc[~seg['Imputed'], 'nMinor'].tolist()) if v != 0]

    n_tai = 0
    for chrom, chrom_seg in seg.groupby('Chromosome', sort=False):
        if chrom not in chrom_arms:
            logger.warning(f"pcgr-hrd: no centromere info for chromosome '{chrom}' - skipping for TAI")
            continue
        chrom_seg = chrom_seg.reset_index(drop=True)
        if len(chrom_seg) == 1:
            continue

        seg_len = chrom_seg['End'] - chrom_seg['Start']
        length_by_minor_cn = {
            v: float(seg_len[(chrom_seg['nMinor'] == v) & ~chrom_seg['Imputed']].sum()) for v in minor_cn_values}
        if length_by_minor_cn:
            ploidy = max(length_by_minor_cn, key=length_by_minor_cn.get)
        else:
            ploidy = None

        if ploidy is None or ploidy == 1 or ploidy % 2 == 0:
            imbalanced = chrom_seg['nMajor'] != chrom_seg['nMinor']
        else:
            balanced = ((chrom_seg['nMajor'] + chrom_seg['nMinor']) == ploidy) & (chrom_seg['nMinor'] != 0)
            imbalanced = ~balanced
        ## telomeric imputed segments (unknown state) are not counted
        imbalanced = imbalanced & ~chrom_seg['Imputed']

        if imbalanced.iloc[0] and chrom_seg['End'].iloc[0] < chrom_arms[chrom]['centromere_start']:
            n_tai += 1
        if imbalanced.iloc[-1] and chrom_seg['Start'].iloc[-1] > chrom_arms[chrom]['centromere_end']:
            n_tai += 1

    return n_tai


FGA_BASELINE_TOTAL_CN = 2
WGD_MIN_FRACTION = 0.5


def _clip_segments_to_chromosomes(seg: pd.DataFrame, chrom_arms: Optional[dict]) -> pd.DataFrame:
    """Drop segments starting beyond, and clip segment ends to, the chromosome length."""
    if chrom_arms is None:
        return seg
    chrom_len = seg['Chromosome'].map({c: v['length'] for c, v in chrom_arms.items()})
    seg = seg[chrom_len.isna() | (seg['Start'] <= chrom_len)].copy()
    chrom_len = chrom_len.loc[seg.index]
    seg['End'] = seg['End'].where(chrom_len.isna() | (seg['End'] <= chrom_len), chrom_len)
    return seg


def calc_genome_doubling(
        cna_df: pd.DataFrame,
        chrom_arms: Optional[dict] = None,
        logger: Optional[logging.Logger] = None) -> Optional[dict]:
    """
    Whole-genome doubling (WGD) call from allele-specific copy number segments,
    following Bielski et al. (Nat Genet 2018, PMID 30013179): a tumor is
    classified as genome-doubled when at least 50% of its autosomal genome has
    a major allele copy number of two or more.

    Unlike a ploidy cut-off, this criterion uses the allele-specific profile
    and so distinguishes a genome-doubled tumor with extensive losses from a
    near-diploid tumor with many gains (both can have an intermediate ploidy).

    Returns a dict with 'wgd_fraction' (fraction of the segmented autosomal
    genome with major allele copy number >= 2) and 'genome_doubled' (bool),
    or None if no usable segments are present. Segments with an imputed
    allele-specific copy number are excluded (numerator and denominator).
    """
    if logger is None:
        logger = logging.getLogger("pcgr-hrd")
    seg = _clip_segments_to_chromosomes(_prepare_segments(cna_df, logger, exclude_imputed=True), chrom_arms)
    if len(seg) == 0:
        return None
    seg_length = (seg['End'] - seg['Start']).clip(lower=0)
    total_length = seg_length.sum()
    if total_length <= 0:
        return None
    wgd_fraction = round(float(seg_length[seg['nMajor'] >= 2].sum() / total_length), 4)
    return {'wgd_fraction': wgd_fraction, 'genome_doubled': bool(wgd_fraction >= WGD_MIN_FRACTION)}


def calc_fraction_genome_altered(
        cna_df: pd.DataFrame,
        chrom_arms: Optional[dict] = None,
        logger: Optional[logging.Logger] = None) -> Optional[float]:
    """
    Fraction of the (autosomal, segmented) genome with an altered total copy
    number, adapted from the cBioPortal definition (cbioportal-core,
    FractionGenomeAlteredCalculator: segment length with |segment mean log2
    ratio| >= 0.2, divided by the total segment length, per sample). For
    integer total copy numbers this is equivalent to "total copy number != 2"
    (adjacent integer states differ by at least log2(3/2) = 0.58), which is
    what is implemented here. Differences to cBioPortal: the input here is
    purity-adjusted allele-specific integer copy number rather than raw
    segment means (low-level gains/losses in impure tumors may fall below 0.2
    in the latter), and cBioPortal applies no chromosome filter at import.

    Notes:
      - the baseline is diploid (not the tumor ploidy), as in cBioPortal -
        near-tetraploid / genome-doubled tumors therefore tend towards 1.0
      - autosomes only (1-22), as for the HRD scar scores, to avoid
        sex-dependent baselines on X/Y
      - the denominator is the total length of the autosomal segments in the
        input (segments are clipped to chromosome lengths if 'chrom_arms' is
        provided), not the full genome size - for targeted data this is the
        segmented territory only
    """
    if logger is None:
        logger = logging.getLogger("pcgr-hrd")
    seg = _prepare_segments(cna_df, logger)
    if len(seg) == 0:
        return None

    seg = _clip_segments_to_chromosomes(seg, chrom_arms)

    seg_length = (seg['End'] - seg['Start']).clip(lower=0)
    total_length = seg_length.sum()
    if total_length <= 0:
        return None

    altered = (seg['nMajor'] + seg['nMinor']) != FGA_BASELINE_TOTAL_CN
    return round(float(seg_length[altered].sum() / total_length), 4)


def compute_genomic_instability_score(
        input_cna_segment_fname: str,
        chromsizes_fname: str,
        tumor_ploidy: Optional[float] = None,
        cytoband_fname: Optional[str] = None,
        logger: Optional[logging.Logger] = None) -> dict:
    """
    Compute the combined genomic instability ("HRD") score - HRD-LOH + LST +
    TAI - from a user-provided allele-specific CNA segment file, following
    the scarHRD algorithm (see module docstring for details on fidelity).

    Args:
        input_cna_segment_fname: TSV with columns Chromosome, Start, End,
            nMajor, nMinor (same format required elsewhere in PCGR for CNA
            segment input).
        chromsizes_fname: path to PCGR's 'chromsize.<build>.tsv' reference file.
        tumor_ploidy: optional - if given, also reports scarHRD's ploidy-
            adjusted variant of the score (LST - 15.5*ploidy + HRD-LOH + TAI).
            This adjustment is present in scarHRD's source but is not the
            widely-cited/clinically-thresholded score (that is the plain,
            unweighted sum - see MYRIAD_HRD_POSITIVE_THRESHOLD).
        cytoband_fname: path to the cytoband file of the PCGR data bundle
            ('misc/tsv/cytoband/cytoband.tsv.gz'). If given, the cytoband 'acen'
            bands are used as centromeres for LST/TAI, which gives results
            identical to scarHRD. If None, the centromere coordinates in the
            chromsize file are used (results may differ slightly from scarHRD).
        logger: optional logger instance.

    Returns:
        dict with keys: hrd_loh, lst, tai, hrd_sum, hrd_sum_ploidy_adjusted
        (None if tumor_ploidy not given), n_segments_used.

    NOTE: this score is only validated for whole-genome/whole-exome,
    paired tumor-normal data, and (per Myriad's assay) primarily for
    ovarian cancer - see caveats discussed with the user before using this
    beyond a research-use, non-diagnostic context. It is not meaningful for
    tumor-only or targeted panel sequencing.
    """
    if logger is None:
        logger = logging.getLogger("pcgr-hrd")

    check_file_exists(input_cna_segment_fname, logger=logger)
    cna_df = pd.read_csv(input_cna_segment_fname, sep="\t", na_values=".")

    chrom_arms = read_chrom_arms(chromsizes_fname, logger=logger)
    if cytoband_fname is not None:
        chrom_arms = scarhrd_chrom_arms(chrom_arms, cytoband_fname, logger=logger)

    check_segment_coverage(cna_df, chrom_arms, logger=logger)

    hrd_loh = calc_hrd_loh(cna_df, logger=logger)
    lst = calc_lst(cna_df, chrom_arms, logger=logger)
    tai = calc_tai(cna_df, chrom_arms, logger=logger)
    hrd_sum = hrd_loh + lst + tai

    hrd_sum_ploidy_adjusted = None
    if tumor_ploidy is not None:
        hrd_sum_ploidy_adjusted = round(lst - 15.5 * tumor_ploidy + hrd_loh + tai, 2)

    seg_used = _prepare_segments(cna_df, logger)

    result = {
        'hrd_loh': hrd_loh,
        'lst': lst,
        'tai': tai,
        'hrd_sum': hrd_sum,
        'hrd_sum_ploidy_adjusted': hrd_sum_ploidy_adjusted,
        'n_segments_used': int((~seg_used['Imputed']).sum()),
        'n_segments_imputed': int(seg_used['Imputed'].sum()),
    }
    if result['n_segments_imputed'] > 0:
        logger.info(
            f"Genomic instability score: n = {result['n_segments_imputed']} autosomal segment(s) with imputed "
            f"allele-specific copy number treated as unknown state")
    logger.info(
        f"Genomic instability score: HRD-LOH={hrd_loh}, LST={lst}, TAI={tai}, "
        f"sum={hrd_sum} (Myriad HRD-positive threshold: >= {MYRIAD_HRD_POSITIVE_THRESHOLD}, ovarian cancer only)")
    return result


if __name__ == "__main__":
    import argparse
    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser(
        description="Compute HRD-LOH/LST/TAI genomic instability scores from an allele-specific CNA segment file")
    parser.add_argument("cna_segment_file", help="TSV with columns Chromosome, Start, End, nMajor, nMinor")
    parser.add_argument("chromsizes_file", help="PCGR chromsize.<build>.tsv reference file")
    parser.add_argument("--ploidy", type=float, default=None, help="Tumor ploidy (optional)")
    parser.add_argument("--cytoband_file", default=None,
                        help="PCGR cytoband file (misc/tsv/cytoband/cytoband.tsv.gz) - use cytoband 'acen' bands "
                             "as centromeres for LST/TAI, for results identical to scarHRD")
    args = parser.parse_args()

    res = compute_genomic_instability_score(
        args.cna_segment_file, args.chromsizes_file, tumor_ploidy=args.ploidy, cytoband_fname=args.cytoband_file)
    for k, v in res.items():
        print(f"{k}\t{v}")
    fga = calc_fraction_genome_altered(
        pd.read_csv(args.cna_segment_file, sep="\t", na_values="."),
        chrom_arms=read_chrom_arms(args.chromsizes_file))
    print(f"fraction_genome_altered\t{fga}")
    wgd = calc_genome_doubling(
        pd.read_csv(args.cna_segment_file, sep="\t", na_values="."),
        chrom_arms=read_chrom_arms(args.chromsizes_file))
    print(f"genome_doubling\t{wgd}")
