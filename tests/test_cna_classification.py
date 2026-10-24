"""
Expected segment labels (VARIANT_CLASS and LOH) from PCGR's CNA classification,
across sex, chromosome type and copy number - run through annotate_cytoband()
and the segment-level classification steps of annotate_cna_segments() (see
classify_segments in conftest.py).

Each sample is tested both with and without sex-chromosome segments: without
them, annotate_cytoband() once read chromosome back as int64 and every autosomal
loss silently became 'neutral'. chr22 is included on purpose (it was once
missing from pcgr_vars.AUTOSOMES).
"""

import pytest

SEXES = ['MALE', 'FEMALE', 'UNKNOWN']

# (n_major, n_minor) -> (variant_class, loh) at tumor ploidy 2.0 and PCGR default
# thresholds in combined mode: amplification CN >= max(5, 2.0 x 2.5) = 5,
# gain CN >= max(3, 2.0 x 1.5) = 3, hetdel 0 < CN <= min(1, 2.0 x 0.5) = 1,
# LOH baseline round(2.0) = 2
AUTOSOME_DIPLOID = {
    (0, 0): ('homdel', '.'),
    (1, 0): ('hetdel', 'deletion'),
    (1, 1): ('neutral', '.'),
    (2, 0): ('neutral', 'copy-neutral'),
    (2, 1): ('gain', '.'),
    (3, 0): ('gain', 'amplification'),
    (3, 1): ('gain', '.'),
    (4, 1): ('amplification', '.'),
    (6, 2): ('amplification', '.'),
}

# chrX at tumor ploidy 2.0, by sex. Gains/amplifications use the same thresholds
# as autosomes for every sex (no sex-specific baseline); losses and LOH differ:
# FEMALE - XX, labelled like autosomes
# MALE - hemizygous: CN 0 is 'hemdel', CN 1 is normal, no LOH calls
# UNKNOWN - losses and LOH on sex chromosomes are skipped
CHRX_DIPLOID = {
    'FEMALE': AUTOSOME_DIPLOID,
    'MALE': {
        (0, 0): ('hemdel', '.'),
        (1, 0): ('neutral', '.'),
        (1, 1): ('neutral', '.'),
        (2, 0): ('neutral', '.'),
        (2, 1): ('gain', '.'),
        (3, 0): ('gain', '.'),
        (3, 1): ('gain', '.'),
        (4, 1): ('amplification', '.'),
        (6, 2): ('amplification', '.'),
    },
    'UNKNOWN': {
        (0, 0): ('neutral', '.'),
        (1, 0): ('neutral', '.'),
        (1, 1): ('neutral', '.'),
        (2, 0): ('neutral', '.'),
        (2, 1): ('gain', '.'),
        (3, 0): ('gain', '.'),
        (3, 1): ('gain', '.'),
        (4, 1): ('amplification', '.'),
        (6, 2): ('amplification', '.'),
    },
}

# autosomes at tumor ploidy 3.5, thresholds as in a typical polyploid run
# (gain_abs 4, amp_abs 5): amplification CN >= max(5, 3.5 x 2.5 = 8.8) -> 9,
# gain CN >= max(4, round(3.5 x 1.5, 1) = 5.2) -> 6, hetdel CN <= min(1, 1.8) = 1,
# LOH baseline round(3.5) = 4
AUTOSOME_POLYPLOID = {
    (0, 0): ('homdel', '.'),
    (1, 0): ('hetdel', 'deletion'),
    (1, 1): ('neutral', '.'),
    (3, 0): ('neutral', 'deletion'),
    (4, 0): ('neutral', 'copy-neutral'),
    (3, 2): ('neutral', '.'),
    (4, 2): ('gain', '.'),
    (6, 2): ('gain', '.'),
    (6, 3): ('amplification', '.'),
    (9, 0): ('amplification', 'amplification'),
}


def _segments(chrom, cases):
    """One 0.5 Mb segment per (n_major, n_minor) case, 2 Mb apart."""
    segments = []
    for i, (n_major, n_minor) in enumerate(cases):
        start = 1000000 + i * 2000000
        segments.append((chrom, start, start + 500000, n_major, n_minor))
    return segments


def _check(labels, chrom, expected):
    for (chrom_, start, end, n_major, n_minor) in _segments(chrom, expected):
        got = tuple(labels.loc[f'chr{chrom_}:{start}-{end}', ['variant_class', 'loh']])
        assert got == expected[(n_major, n_minor)], \
            f'chr{chrom} {n_major}/{n_minor}: got {got}, expected {expected[(n_major, n_minor)]}'


@pytest.mark.parametrize('sex', SEXES)
@pytest.mark.parametrize('with_sex_chromosomes', [False, True],
                         ids=['autosomes_only', 'with_chrX_chrY'])
def test_labels_diploid(classify_segments, sex, with_sex_chromosomes):
    segments = _segments('1', AUTOSOME_DIPLOID) + _segments('22', AUTOSOME_DIPLOID)
    if with_sex_chromosomes:
        segments += _segments('X', CHRX_DIPLOID[sex]) + [('Y', 2000000, 2500000, 0, 0)]
    labels = classify_segments(segments, sex)

    _check(labels, '1', AUTOSOME_DIPLOID)
    _check(labels, '22', AUTOSOME_DIPLOID)
    if with_sex_chromosomes:
        _check(labels, 'X', CHRX_DIPLOID[sex])
        expected_y = {'MALE': 'hemdel', 'FEMALE': 'homdel', 'UNKNOWN': 'neutral'}[sex]
        assert labels.loc['chrY:2000000-2500000', 'variant_class'] == expected_y


@pytest.mark.parametrize('sex', SEXES)
@pytest.mark.parametrize('with_sex_chromosomes', [False, True],
                         ids=['autosomes_only', 'with_chrX'])
def test_labels_polyploid(classify_segments, sex, with_sex_chromosomes):
    segments = _segments('1', AUTOSOME_POLYPLOID) + _segments('22', AUTOSOME_POLYPLOID)
    if with_sex_chromosomes:
        segments.append(('X', 30000000, 30500000, 2, 2))
    labels = classify_segments(segments, sex, tumor_ploidy=3.5, gain_abs=4, amp_abs=5)

    _check(labels, '1', AUTOSOME_POLYPLOID)
    _check(labels, '22', AUTOSOME_POLYPLOID)


@pytest.mark.parametrize('sex', SEXES)
def test_imputed_major_minor_split_is_not_called_as_loh(classify_segments, sex):
    # rescued amplicons carry an imputed (placeholder) major/minor split - they keep
    # their total-CN class but must never be called LOH
    segments = [('1', 1000000, 1500000, 2, 0), ('1', 3000000, 3500000, 6, 0)]
    imputed = ['chr1:1000000-1500000', 'chr1:3000000-3500000']
    labels = classify_segments(segments, sex, imputed=imputed)
    assert tuple(labels.loc['chr1:1000000-1500000']) == ('neutral', '.')
    assert tuple(labels.loc['chr1:3000000-3500000']) == ('amplification', '.')
