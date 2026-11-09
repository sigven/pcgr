import gzip
import logging

import pandas as pd
import pytest

from pcgr.cna import generate_oncokb_cna_input
from pcgr.config import get_oncotree_lineage, verify_oncotree_code

LOGGER = logging.getLogger('test_oncokb')


@pytest.fixture
def oncotree_refdata_dir(tmp_path):
    tsv_dir = tmp_path / 'phenotype' / 'tsv'
    tsv_dir.mkdir(parents=True)
    with gzip.open(tsv_dir / 'phenotype_oncotree.tsv.gz', 'wt') as fh:
        fh.write('ot_code\tprimary_site\tot_name\tot_code_path\n')
        fh.write('BONE\tBone\tBone\tBONE\n')
        fh.write('OS\tBone\tOsteosarcoma\tBONE-OS\n')
        fh.write('PNS\tPeripheral Nervous System\tPeripheral Nervous System\tPNS\n')
        fh.write('NST\tPeripheral Nervous System\tNerve Sheath Tumor\tPNS-NST\n')
        fh.write('MPNST\tPeripheral Nervous System\tMalignant Peripheral Nerve Sheath Tumor\tPNS-NST-MPNST\n')
        fh.write('MFH\tSoft Tissue\tUndifferentiated Pleomorphic Sarcoma\tSOFT_TISSUE-MFH\n')
    return str(tmp_path)


def test_oncotree_code_matching_tumor_site_is_kept(oncotree_refdata_dir):
    assert verify_oncotree_code('OS', 'Bone', oncotree_refdata_dir, LOGGER) == 'OS'


@pytest.mark.parametrize('code,tumor_site', [
    ('MPNST', 'Soft Tissue'),
    ('MFH', 'Bone'),
])
def test_oncotree_code_with_other_primary_site_is_kept_with_warning(
        oncotree_refdata_dir, caplog, code, tumor_site):
    with caplog.at_level(logging.WARNING, logger='test_oncokb'):
        assert verify_oncotree_code(code, tumor_site, oncotree_refdata_dir, LOGGER) == code
    assert f"does not match the selected PCGR tumor site '{tumor_site}'" in caplog.text


def test_unrecognized_oncotree_code_is_dropped(oncotree_refdata_dir):
    assert verify_oncotree_code('NOT_A_CODE', 'Bone', oncotree_refdata_dir, LOGGER) is None


def test_oncotree_lineage_lists_code_and_ancestors(oncotree_refdata_dir):
    assert get_oncotree_lineage('MPNST', oncotree_refdata_dir) == [
        'Peripheral Nervous System', 'Nerve Sheath Tumor', 'Malignant Peripheral Nerve Sheath Tumor']
    ## ancestors missing from the file (SOFT_TISSUE) are skipped
    assert get_oncotree_lineage('MFH', oncotree_refdata_dir) == ['Undifferentiated Pleomorphic Sarcoma']
    assert get_oncotree_lineage('NOT_A_CODE', oncotree_refdata_dir) == []
    assert get_oncotree_lineage(None, oncotree_refdata_dir) == []


def test_oncokb_cna_input_has_amplifications_and_homdels_only(tmp_path):
    cna_df = pd.DataFrame({
        'SAMPLE_ID': 's1',
        'SYMBOL': ['CDK4', 'MDM2', 'CDKN2A', 'FAT1', 'JAK2', 'TP53'],
        'VARIANT_CLASS': ['amplification', 'amplification', 'homdel', 'hetdel', 'hemdel', 'gain'],
        'VAR_ID': [f'v{i}' for i in range(6)],
    })
    out = tmp_path / 'oncokb_input.cna.txt'
    generate_oncokb_cna_input(cna_df, str(out), LOGGER)
    oncokb_df = pd.read_csv(out, sep='\t')
    assert oncokb_df[['Hugo_Symbol', 'Copy_Number_Alteration']].values.tolist() == [
        ['CDK4', 'Amplification'], ['MDM2', 'Amplification'], ['CDKN2A', 'Deletion']]
