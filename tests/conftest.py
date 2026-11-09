import gzip
import logging

import pandas as pd
import pybedtools
import pytest

from pcgr.cna import (
    _annotate_amplifications, _annotate_deletions, _annotate_gains, _annotate_loh,
    _assign_variant_class, annotate_cytoband)

LOGGER = logging.getLogger('test_cna')

# cytoband name format: cytoband|arm|arm_length|focal_threshold
CYTOBANDS = [
    ('1', 0, 125000000, '1p11.1|1p|125000000|25000000'),
    ('1', 125000000, 248956422, '1q21.1|1q|123956422|25000000'),
    ('13', 0, 114364328, '13q14.2|13q|114364328|25000000'),
    ('22', 0, 50818468, '22q11.23|22q|50818468|10000000'),
    ('X', 0, 156040895, 'Xp11.4|Xp|156040895|30000000'),
    ('Y', 0, 57227415, 'Yq11.221|Yq|57227415|10000000'),
]


@pytest.fixture
def refdata_dir(tmp_path):
    bed_dir = tmp_path / 'misc' / 'bed' / 'cytoband'
    bed_dir.mkdir(parents=True)
    with gzip.open(bed_dir / 'cytoband.bed.gz', 'wt') as fh:
        for chrom, start, end, name in CYTOBANDS:
            fh.write(f'{chrom}\t{start}\t{end}\t{name}\n')
    return str(tmp_path)


def segments_bt(segments):
    """segments: (chromosome, start, end, n_major, n_minor), chromosome without 'chr'."""
    df = pd.DataFrame([
        (c, s, e, f'chr{c}:{s}-{e}|{nmaj}|{nmin}') for c, s, e, nmaj, nmin in segments],
        columns=['Chromosome', 'Start', 'End', 'Name'])
    return pybedtools.BedTool.from_dataframe(df)


@pytest.fixture
def classify_segments(tmp_path, refdata_dir):
    """
    Run the segment-level classification chain of annotate_cna_segments() -
    annotate_cytoband, then amplification/gain/deletion/variant class/LOH - on
    segments given as (chromosome, start, end, n_major, n_minor). Returns a frame
    indexed by segment_id with 'variant_class' and 'loh'.
    """
    def _classify(segments, sex, tumor_ploidy=2.0, imputed=(),
                  amp_abs=5, amp_rel=2.5, gain_abs=3, gain_rel=1.5,
                  del_abs=1, del_rel=0.5, mode='combined'):
        cytoband_df = annotate_cytoband(segments_bt(segments), str(tmp_path), refdata_dir, LOGGER)
        ann = cytoband_df['segment_name'].str.split('\\|', expand=True)
        ann.columns = ['segment_id', 'n_major', 'n_minor', 'chromosome_arm', 'cytoband', 'event_type']
        df = pd.concat(
            [cytoband_df[['chromosome', 'segment_start', 'segment_end']].reset_index(drop=True),
             ann.reset_index(drop=True)], axis=1).astype({'n_major': 'int', 'n_minor': 'int'})
        df['cn_major_minor_imputed'] = df['segment_id'].isin(imputed)
        df, _ = _annotate_amplifications(df, tumor_ploidy, amp_abs, amp_rel, mode, LOGGER)
        df = _annotate_gains(df, tumor_ploidy, gain_abs, gain_rel, mode, LOGGER)
        df, baseline_cn = _annotate_deletions(df, tumor_ploidy, del_abs, del_rel, mode, sex, LOGGER)
        df = _assign_variant_class(df)
        df = _annotate_loh(df, baseline_cn, sex, LOGGER)
        return df.set_index('segment_id')[['variant_class', 'loh']]
    return _classify
