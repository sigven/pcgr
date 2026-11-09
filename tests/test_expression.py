"""Tests for identifier handling in parse_expression (pcgr/expression.py)."""
import logging

import pandas as pd
import pytest

from pcgr.expression import parse_expression

LOGGER = logging.getLogger('test_expression')

GENE_ID = 'ENSG00000121410'
TRANSCRIPT_ID = 'ENST00000263100'
REFSEQ_ID = 'NM_130786.4'
DOTTED_ALIAS = 'KB-226F1.2'


@pytest.fixture
def yaml_data(tmp_path):
    """Minimal reference data, with identifiers formatted as in the gene index of the data bundle."""
    xref_dir = tmp_path / 'gene' / 'tsv' / 'gene_transcript_xref'
    xref_dir.mkdir(parents=True)
    gene_properties = {
        'ENTREZGENE': '1', 'ENSEMBL_GENE_ID': GENE_ID, 'BIOTYPE': 'protein_coding',
        'AMBIGUOUS_ID': False, 'SYMBOL': 'A1BG', 'GENENAME': 'alpha-1-B glycoprotein'}
    pd.DataFrame([
        {'ID': GENE_ID, 'ID_TYPE': 'ensembl_gene_id', **gene_properties},
        {'ID': TRANSCRIPT_ID, 'ID_TYPE': 'ensembl_transcript_id', **gene_properties},
        {'ID': REFSEQ_ID, 'ID_TYPE': 'refseq_transcript_id', **gene_properties},
        {'ID': DOTTED_ALIAS, 'ID_TYPE': 'symbol_alias', **gene_properties},
    ]).to_csv(xref_dir / 'gene_index.tsv.gz', sep='\t', index=False)
    pd.DataFrame({
        'ensembl_transcript_id': [TRANSCRIPT_ID], 'refseq_transcript_id': [REFSEQ_ID]
    }).to_csv(xref_dir / 'gene_transcript_xref.tsv.gz', sep='\t', index=False)
    return {'sample_id': 'SAMPLE', 'reference_data': {'path': str(tmp_path)}}


def parse(tmp_path, yaml_data, records):
    fname = tmp_path / 'expression.tsv'
    pd.DataFrame(records, columns=['TargetID', 'TPM']).to_csv(fname, sep='\t', index=False)
    return parse_expression(str(fname), yaml_data, logger=LOGGER)


def test_ensembl_transcript_version_is_removed(tmp_path, yaml_data):
    result = parse(tmp_path, yaml_data, [(f'{TRANSCRIPT_ID}.5', 7.0)])
    assert result['transcript_identifier'] == 'Ensembl'
    assert result['transcript']['ID'].tolist() == [TRANSCRIPT_ID]
    assert result['transcript']['ENSEMBL_TRANSCRIPT_ID'].tolist() == [TRANSCRIPT_ID]
    assert result['gene']['TPM_GENE'].tolist() == [7.0]


def test_ensembl_gene_version_is_removed(tmp_path, yaml_data):
    result = parse(tmp_path, yaml_data, [(f'{GENE_ID}.13', 3.0)])
    assert result['transcript'] is None
    assert result['gene']['ENSEMBL_GENE_ID'].tolist() == [GENE_ID]
    assert result['gene']['TPM_GENE'].tolist() == [3.0]


def test_versions_of_same_transcript_keep_highest_tpm(tmp_path, yaml_data):
    result = parse(tmp_path, yaml_data,
                   [(f'{TRANSCRIPT_ID}.4', 2.0), (f'{TRANSCRIPT_ID}.5', 7.0)])
    assert result['transcript']['TPM'].tolist() == [7.0]


def test_unversioned_ensembl_identifier_is_unchanged(tmp_path, yaml_data):
    result = parse(tmp_path, yaml_data, [(TRANSCRIPT_ID, 7.0)])
    assert result['transcript']['ID'].tolist() == [TRANSCRIPT_ID]


def test_refseq_version_is_kept(tmp_path, yaml_data):
    """RefSeq transcript identifiers are versioned in the gene index."""
    result = parse(tmp_path, yaml_data, [(REFSEQ_ID, 7.0)])
    assert result['transcript_identifier'] == 'RefSeq'
    assert result['transcript']['ID'].tolist() == [REFSEQ_ID]
    assert result['transcript']['ENSEMBL_TRANSCRIPT_ID'].tolist() == [TRANSCRIPT_ID]


def test_symbol_alias_ending_with_number_is_kept(tmp_path, yaml_data):
    result = parse(tmp_path, yaml_data, [(DOTTED_ALIAS, 7.0)])
    assert result['transcript'] is None
    assert result['gene']['SYMBOL'].tolist() == ['A1BG']
    assert result['gene']['TPM_GENE'].tolist() == [7.0]
