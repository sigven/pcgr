import gzip
import logging
import os

import pytest

from pcgr.biomarker import is_exon_specific_fusion_alias, load_all_biomarkers

VARIANT_COLS = [
    'biomarker_source', 'biomarker_source_datestamp', 'variant_id', 'molecular_profile_id',
    'variant_alias', 'variant_name_primary', 'alteration_type', 'alias_type',
    'variant_consequence', 'variant_exon', 'gene', 'symbol', 'entrezgene']

CLINICAL_COLS = [
    'biomarker_source', 'biomarker_source_datestamp', 'biomarker_entity', 'molecular_profile_id',
    'molecular_profile_name', 'molecular_profile_type', 'variant_id', 'evidence_id',
    'variant_origin', 'primary_site', 'cancer_type', 'evidence_type', 'evidence_level',
    'evidence_direction', 'clinical_significance']


def _fusion_variant(variant_id, alias):
    return {
        'biomarker_source': 'civic', 'biomarker_source_datestamp': '20261005',
        'variant_id': variant_id, 'molecular_profile_id': variant_id,
        'variant_alias': alias, 'variant_name_primary': 'EWSR1::FLI1 Fusion',
        'alteration_type': 'FUSION', 'alias_type': 'other_gene',
        'variant_consequence': 'transcript_fusion', 'variant_exon': '.',
        'gene': 'EWSR1::FLI1', 'symbol': '.', 'entrezgene': '2130::2313'}


def _evidence(variant_id, evidence_id, evidence_type, clinical_significance):
    return {
        'biomarker_source': 'civic', 'biomarker_source_datestamp': '20261005',
        'biomarker_entity': 'TRUE', 'molecular_profile_id': variant_id,
        'molecular_profile_name': 'EWSR1::FLI1', 'molecular_profile_type': 'Single',
        'variant_id': variant_id, 'evidence_id': evidence_id, 'variant_origin': 'Somatic',
        'primary_site': 'Soft Tissue', 'cancer_type': 'Ewing Sarcoma',
        'evidence_type': evidence_type, 'evidence_level': 'B: Clinical evidence',
        'evidence_direction': 'Supports', 'clinical_significance': clinical_significance}


def _write_tsv(fname, cols, rows):
    with gzip.open(fname, mode='wt') as f:
        f.write('\t'.join(cols) + '\n')
        for row in rows:
            f.write('\t'.join(row[c] for c in cols) + '\n')


@pytest.fixture
def refdata_dir(tmp_path):
    tsv_dir = tmp_path / 'biomarker' / 'tsv'
    tsv_dir.mkdir(parents=True)

    ## CIViC: exon-specific variant (164) listed before the generic variant (706)
    civic_variants = [
        _fusion_variant('164', 'EWSR1-FLI1 TYPE 1'),
        _fusion_variant('164', 'EWSR1-FLI1 E7-E6'),
        _fusion_variant('164', 'EWSR1-FLI1'),
        _fusion_variant('706', 'EWS-FLI'),
        _fusion_variant('706', 'EWSR1::FLI1'),
        _fusion_variant('706', 'EWSR1-FLI1'),
    ]
    civic_clinical = [
        _evidence('164', 'EID431', 'Prognostic', 'Better Outcome'),
        _evidence('706', 'EID478', 'Diagnostic', 'Positive'),
    ]
    _write_tsv(tsv_dir / 'civic.variant.tsv.gz', VARIANT_COLS, civic_variants)
    _write_tsv(tsv_dir / 'civic.clinical.tsv.gz', CLINICAL_COLS, civic_clinical)
    _write_tsv(tsv_dir / 'cgi.variant.tsv.gz', VARIANT_COLS, [])
    _write_tsv(tsv_dir / 'cgi.clinical.tsv.gz', CLINICAL_COLS, [])
    return str(tmp_path)


@pytest.mark.parametrize('alias,expected', [
    ('EWSR1-FLI1 E7-E6', True),
    ('EGFR-RAD51 e24-e4', True),
    ('EWSR1-FLI1 TYPE 1', True),
    ('EWSR1-FLI1', False),
    ('EWSR1::FLI1', False),
    ('T(11;22)(Q24;Q12)', False),
    ('EWSR1::v', False),
])
def test_is_exon_specific_fusion_alias(alias, expected):
    assert is_exon_specific_fusion_alias(alias) is expected


def test_fusion_biomarkers_aggregate_generic_variants_and_skip_exon_specific(refdata_dir):
    biomarkers = load_all_biomarkers(
        logging.getLogger('test'), refdata_dir,
        biomarker_vartype='FUSION', biomarker_variant_origin='Somatic')

    actionable = biomarkers['actionable_df'].set_index('aberration_key')['biomarker_match']
    items = actionable['EWSR1::FLI1_transcript_fusion'].split(',')

    ## generic variant 706 is matched even though it is not the first row for the key,
    ## exon-specific variant 164 is skipped, and each variant appears only once
    assert items == [
        'civic|706|EID478:Soft_Tissue:Positive:B:Diagnostic:Somatic|by_fusion']
