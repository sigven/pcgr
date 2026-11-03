import logging

import pytest

from pcgr.hrd_external import (
    HrdInputError, parse_hrd_generic, read_external_hrd, resolve_hrd_status, write_external_hrd)
from pcgr.pcgr_vars import (
    HRD_EXTERNAL_COLUMNS, HRD_STATUS_DEFICIENT, HRD_STATUS_INDETERMINATE, HRD_STATUS_PROFICIENT)

LOGGER = logging.getLogger('test_hrd_external')


def write_hrd(tmp_path, rows):
    fname = tmp_path / 'hrd.tsv'
    fname.write_text('\n'.join('\t'.join(r) for r in rows) + '\n')
    return str(fname)


def test_status_only(tmp_path):
    fname = write_hrd(tmp_path, [
        ['HRD_CALLER', 'HRD_SCORE_TYPE', 'HRD_STATUS'],
        ['MyAssay', 'genomic_scar', 'HRD-positive']])
    res = resolve_hrd_status(parse_hrd_generic(fname))
    assert res['HRD_STATUS'] == HRD_STATUS_DEFICIENT
    assert res['HRD_STATUS_BASIS'] == 'caller'
    assert res['HRD_SCORE'] is None


@pytest.mark.parametrize('value,expected', [
    ('HRD', HRD_STATUS_DEFICIENT), ('hr deficient', HRD_STATUS_DEFICIENT),
    ('Positive', HRD_STATUS_DEFICIENT), ('HRP', HRD_STATUS_PROFICIENT),
    ('HRD-negative', HRD_STATUS_PROFICIENT), ('equivocal', HRD_STATUS_INDETERMINATE)])
def test_status_synonyms(tmp_path, value, expected):
    fname = write_hrd(tmp_path, [
        ['HRD_CALLER', 'HRD_SCORE_TYPE', 'HRD_STATUS'], ['X', 'other', value]])
    assert parse_hrd_generic(fname)['HRD_STATUS'] == expected


def test_score_with_threshold_derives_status(tmp_path):
    fname = write_hrd(tmp_path, [
        ['HRD_CALLER', 'HRD_SCORE_TYPE', 'HRD_SCORE', 'HRD_LOH', 'HRD_LST', 'HRD_TAI'],
        ['HiSNP', 'Genomic scar', '42', '14', '18', '10']])
    res = resolve_hrd_status(parse_hrd_generic(fname), threshold=42)
    assert res['HRD_SCORE_TYPE'] == 'genomic_scar'
    assert res['HRD_STATUS'] == HRD_STATUS_DEFICIENT
    assert res['HRD_STATUS_BASIS'] == 'derived'
    assert res['HRD_SCORE_THRESHOLD_USED'] == 42.0
    assert (res['HRD_LOH'], res['HRD_LST'], res['HRD_TAI']) == (14, 18, 10)
    assert resolve_hrd_status(parse_hrd_generic(fname), threshold=42.5)['HRD_STATUS'] == HRD_STATUS_PROFICIENT


def test_threshold_takes_precedence_over_status(tmp_path):
    fname = write_hrd(tmp_path, [
        ['HRD_CALLER', 'HRD_SCORE_TYPE', 'HRD_STATUS', 'HRD_SCORE'],
        ['CHORD', 'mutational_signature', 'HRP', '0.7']])
    res = resolve_hrd_status(parse_hrd_generic(fname), threshold=0.5, logger=LOGGER)
    assert res['HRD_STATUS'] == HRD_STATUS_DEFICIENT
    assert res['HRD_STATUS_BASIS'] == 'derived'


def test_score_without_threshold_or_status_fails(tmp_path):
    fname = write_hrd(tmp_path, [
        ['HRD_CALLER', 'HRD_SCORE_TYPE', 'HRD_SCORE'], ['X', 'genomic_scar', '50']])
    with pytest.raises(HrdInputError, match='input_hrd_threshold'):
        resolve_hrd_status(parse_hrd_generic(fname))


def test_threshold_without_score_fails(tmp_path):
    fname = write_hrd(tmp_path, [
        ['HRD_CALLER', 'HRD_SCORE_TYPE', 'HRD_STATUS'], ['X', 'genomic_scar', 'HRD']])
    with pytest.raises(HrdInputError, match='HRD_SCORE'):
        resolve_hrd_status(parse_hrd_generic(fname), threshold=42)


@pytest.mark.parametrize('rows,match', [
    ([['HRD_CALLER', 'HRD_STATUS'], ['X', 'HRD']], 'HRD_SCORE_TYPE'),
    ([['HRD_CALLER', 'HRD_SCORE_TYPE'], ['X', 'genomic_scar']], "'HRD_STATUS' and/or"),
    ([['HRD_CALLER', 'HRD_SCORE_TYPE', 'HRD_STATUS'], ['X', 'gis', 'HRD']], 'Unrecognised HRD score type'),
    ([['HRD_CALLER', 'HRD_SCORE_TYPE', 'HRD_STATUS'], ['X', 'other', 'maybe']], 'Unrecognised HRD status'),
    ([['HRD_CALLER', 'HRD_SCORE_TYPE', 'HRD_SCORE'], ['X', 'other', 'high']], 'must be numeric'),
    ([['HRD_CALLER', 'HRD_SCORE_TYPE', 'HRD_STATUS', 'HRD_LOH'], ['X', 'other', 'HRD', '2.5']], 'non-negative integer'),
    ([['HRD_CALLER', 'HRD_SCORE_TYPE', 'HRD_STATUS', 'HRD_SCORE'], ['X', 'other', '.', '.']], 'cannot both be empty'),
    ([['HRD_CALLER', 'HRD_SCORE_TYPE', 'HRD_STATUS'], ['X', 'other', 'HRD'], ['Y', 'other', 'HRP']], 'single data row'),
])
def test_invalid_files(tmp_path, rows, match):
    with pytest.raises(HrdInputError, match=match):
        parse_hrd_generic(write_hrd(tmp_path, rows))


def test_invalid_file_exits(tmp_path):
    fname = write_hrd(tmp_path, [['HRD_CALLER', 'HRD_STATUS'], ['X', 'HRD']])
    with pytest.raises(SystemExit):
        read_external_hrd(fname, logger=LOGGER)


def test_write_external_hrd(tmp_path):
    fname = write_hrd(tmp_path, [
        ['HRD_CALLER', 'HRD_SCORE_TYPE', 'HRD_STATUS'], ['MyAssay', 'composite', 'HRD']])
    out = tmp_path / 'out.hrd.tsv'
    write_external_hrd(read_external_hrd(fname, logger=LOGGER), 'S1', str(out))
    header, row = out.read_text().splitlines()
    assert header.split('\t') == HRD_EXTERNAL_COLUMNS
    values = dict(zip(HRD_EXTERNAL_COLUMNS, row.split('\t')))
    assert values['SAMPLE_ID'] == 'S1'
    assert values['HRD_STATUS'] == HRD_STATUS_DEFICIENT
    assert values['HRD_SCORE'] == '.'
