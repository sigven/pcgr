import pytest

from pcgr import pcgr_vars
from pcgr.external_input import (
    ExternalInputError, is_missing, match_synonym, read_single_row_tsv, to_count, to_float,
    write_external_result)
from pcgr.hrd_external import HrdInputError
from pcgr.msi import MsiInputError, parse_msi_generic


def write_tsv(tmp_path, rows, name='input.tsv'):
    fname = tmp_path / name
    fname.write_text('\n'.join('\t'.join(r) for r in rows) + '\n')
    return str(fname)


def test_error_classes_share_base():
    assert issubclass(MsiInputError, ExternalInputError)
    assert issubclass(HrdInputError, ExternalInputError)


@pytest.mark.parametrize('value', [None, '', ' ', '.', 'NA', 'N/A', 'NaN', 'nan'])
def test_is_missing(value):
    assert is_missing(value)
    assert to_float(value, 'X') is None


@pytest.mark.parametrize('value', ['inf', '-inf', 'Infinity'])
def test_to_float_rejects_non_finite(value):
    with pytest.raises(MsiInputError, match='finite'):
        to_float(value, 'MSI_SCORE', MsiInputError)


def test_to_count():
    assert to_count('12', 'N') == 12
    assert to_count('12.0', 'N') == 12
    with pytest.raises(HrdInputError, match='non-negative integer'):
        to_count('-1', 'N', HrdInputError)


def test_match_synonym():
    synonyms = {'MSI-H': {'MSI-H', 'UNSTABLE'}}
    assert match_synonym('msi_h', synonyms) == 'MSI-H'
    assert match_synonym(' Unstable ', synonyms) == 'MSI-H'
    assert match_synonym('stable', synonyms) is None


@pytest.mark.parametrize('rows,match', [
    ([['A'], ['x']], r'Missing required column\(s\) in TEST file .*: B'),
    ([['A', 'B', 'A'], ['x', 'y', 'z']], 'Duplicated column names'),
    ([['A', 'B']], 'no data row'),
    ([['A', 'B'], ['x', 'y'], ['x', 'y']], 'single data row'),
    ([['A', 'B'], ['x']], 'data row has 1 field'),
])
def test_read_single_row_tsv_errors(tmp_path, rows, match):
    with pytest.raises(ExternalInputError, match=match):
        read_single_row_tsv(write_tsv(tmp_path, rows), 'TEST', ['A', 'B'])


def test_read_single_row_tsv_missing_hint(tmp_path):
    fname = write_tsv(tmp_path, [['A'], ['x']])
    with pytest.raises(ExternalInputError, match='header was: A'):
        read_single_row_tsv(fname, 'TEST', ['B'], missing_hint=lambda h: f" (header was: {','.join(h)})")


def test_read_single_row_tsv(tmp_path):
    row = read_single_row_tsv(write_tsv(tmp_path, [['A', 'B'], [' x ', 'y']]), 'TEST', ['A'])
    assert row == {'A': 'x', 'B': 'y'}


def test_empty_file(tmp_path):
    fname = tmp_path / 'empty.tsv'
    fname.write_text('')
    with pytest.raises(HrdInputError, match='HRD file .* does not exist or is empty'):
        read_single_row_tsv(str(fname), 'HRD', ['A'], HrdInputError)


def test_msi_generic_rejects_infinite_score(tmp_path):
    fname = write_tsv(tmp_path, [['MSI_STATUS', 'MSI_CALLER', 'MSI_SCORE'], ['MSS', 'MANTIS', 'inf']])
    with pytest.raises(MsiInputError, match='finite'):
        parse_msi_generic(fname)


def test_msi_generic_hint_for_wrong_file(tmp_path):
    fname = write_tsv(tmp_path, [['Total_Number_of_Sites', 'Number_of_Unstable_Sites', '%'], ['10', '1', '10.0']])
    with pytest.raises(MsiInputError, match="--input_msi_caller msisensor-pro"):
        parse_msi_generic(fname)


def test_write_external_result(tmp_path):
    out = tmp_path / 'out.tsv'
    write_external_result({'MSI_STATUS': 'MSS', 'MSI_SCORE': None}, 'S1', str(out),
                          pcgr_vars.MSI_EXTERNAL_COLUMNS)
    header, row = out.read_text().strip('\n').split('\n')
    assert header.split('\t') == pcgr_vars.MSI_EXTERNAL_COLUMNS
    values = dict(zip(pcgr_vars.MSI_EXTERNAL_COLUMNS, row.split('\t')))
    assert values['SAMPLE_ID'] == 'S1'
    assert values['MSI_STATUS'] == 'MSS'
    assert values['MSI_SCORE'] == '.'
