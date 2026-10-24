import pytest

from pcgr import pcgr_vars
from pcgr.cna import annotate_cytoband

from conftest import LOGGER, segments_bt


def test_autosomes_include_chr22():
    assert pcgr_vars.AUTOSOMES == [str(x) for x in range(1, 23)]


def test_annotate_cytoband_keeps_str_chromosome_without_sex_chromosomes(tmp_path, refdata_dir):
    # autosome-only input (typical male sample with no chrX/Y segments) must not be
    # read back as int64 - the AUTOSOMES masks compare against str
    segments = [('1', 1000000, 2000000, 0, 0), ('13', 40000000, 41000000, 1, 1)]
    df = annotate_cytoband(segments_bt(segments), str(tmp_path), refdata_dir, LOGGER)
    assert df['chromosome'].map(type).eq(str).all()
    assert df['chromosome'].isin(pcgr_vars.AUTOSOMES).all()


@pytest.mark.parametrize('sex', ['MALE', 'FEMALE'])
@pytest.mark.parametrize('with_x_segment', [False, True])
def test_autosomal_zero_copy_segments_are_homdel(classify_segments, sex, with_x_segment):
    segments = [
        ('1', 1000000, 2000000, 0, 0),
        ('13', 40000000, 41000000, 1, 1),
        ('22', 23000000, 24000000, 0, 0),
    ]
    if with_x_segment:
        segments.append(('X', 1000000, 2000000, 1, 0))
    variant_class = classify_segments(segments, sex)['variant_class']
    assert variant_class['chr1:1000000-2000000'] == 'homdel'
    assert variant_class['chr22:23000000-24000000'] == 'homdel'
    assert variant_class['chr13:40000000-41000000'] == 'neutral'
