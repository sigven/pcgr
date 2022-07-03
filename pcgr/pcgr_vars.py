#!/usr/bin/env python

from pcgr._version import __version__

PCGR_VERSION = __version__
DB_VERSION = '20220203'
VEP_VERSION = '105'
GENCODE_VERSION = '39'
NCBI_BUILD_MAF = 'GRCh38'
VEP_ASSEMBLY = 'GRCh38'
MAX_VARIANTS_FOR_REPORT = 500000

tsites = {
    0: 'Any',
    1: 'Adrenal Gland',
    2: 'Ampulla of Vater',
    3: 'Biliary Tract',
    4: 'Bladder/Urinary Tract',
    5: 'Bone',
    6: 'Breast',
    7: 'Cervix',
    8: 'CNS/Brain',
    9: 'Colon/Rectum',
    10: 'Esophagus/Stomach',
    11: 'Eye',
    12: 'Head and Neck',
    13: 'Kidney',
    14: 'Liver',
    15: 'Lung',
    16: 'Lymphoid',
    17: 'Myeloid',
    18: 'Ovary/Fallopian Tube',
    19: 'Pancreas',
    20: 'Peripheral Nervous System',
    21: 'Peritoneum',
    22: 'Pleura',
    23: 'Prostate',
    24: 'Skin',
    25: 'Soft Tissue',
    26: 'Testis',
    27: 'Thymus',
    28: 'Thyroid',
    29: 'Uterus',
    30: 'Vulva/Vagina'
}

tumor_sites = '\n'.join([f'{k} = {tsites[k]}' for k in tsites]) # for displaying in help

GE_panels = {
      0: "CPSR exploratory cancer predisposition panel (n = 433, Genomics England PanelApp / TCGA Germline Study / Cancer Gene Census / Other)",
      1: "Adult solid tumours cancer susceptibility (Genomics England PanelApp)",
      2: "Adult solid tumours for rare disease (Genomics England PanelApp)",
      3: "Bladder cancer pertinent cancer susceptibility (Genomics England PanelApp)",
      4: "Brain cancer pertinent cancer susceptibility (Genomics England PanelApp)",
      5: "Breast cancer pertinent cancer susceptibility (Genomics England PanelApp)",
      6: "Childhood solid tumours cancer susceptibility (Genomics England PanelApp)",
      7: "Colorectal cancer pertinent cancer susceptibility (Genomics England PanelApp)",
      8: "Endometrial cancer pertinent cancer susceptibility (Genomics England PanelApp)",
      9: "Familial Tumours Syndromes of the central & peripheral Nervous system (Genomics England PanelApp)",
      10: "Familial breast cancer (Genomics England PanelApp)",
      11: "Familial melanoma (Genomics England PanelApp)",
      12: "Familial prostate cancer (Genomics England PanelApp)",
      13: "Familial rhabdomyosarcoma (Genomics England PanelApp)",
      14: "GI tract tumours (Genomics England PanelApp)",
      15: "Genodermatoses with malignancies (Genomics England PanelApp)",
      16: "Haematological malignancies cancer susceptibility (Genomics England PanelApp)",
      17: "Haematological malignancies for rare disease (Genomics England PanelApp)",
      18: "Head and neck cancer pertinent cancer susceptibility (Genomics England PanelApp)",
      19: "Inherited MMR deficiency (Lynch Syndrome) - Genomics England PanelApp",
      20: "Inherited non-medullary thyroid cancer (Genomics England PanelApp)",
      21: "Inherited ovarian cancer (without breast cancer) (Genomics England PanelApp)",
      22: "Inherited pancreatic cancer (Genomics England PanelApp)",
      23: "Inherited polyposis (Genomics England PanelApp)",
      24: "Inherited predisposition to acute myeloid leukaemia (AML) - Genomics England PanelApp",
      25: "Inherited predisposition to GIST (Genomics England PanelApp)",
      26: "Inherited renal cancer (Genomics England PanelApp)",
      27: "Inherited phaeochromocytoma and paraganglioma (Genomics England PanelApp)",
      28: "Melanoma pertinent cancer susceptibility (Genomics England PanelApp)",
      29: "Multiple endocrine tumours (Genomics England PanelApp)",
      30: "Multiple monogenic benign skin tumours (Genomics England PanelApp)",
      31: "Neuroendocrine cancer pertinent cancer susceptibility (Genomics England PanelApp)",
      32: "Neurofibromatosis Type 1 (Genomics England PanelApp)",
      33: "Ovarian cancer pertinent cancer susceptibility (Genomics England PanelApp)",
      34: "Parathyroid Cancer (Genomics England PanelApp)",
      35: "Prostate cancer pertinent cancer susceptibility (Genomics England PanelApp)",
      36: "Renal cancer pertinent cancer susceptibility (Genomics England PanelApp)",
      37: "Rhabdoid tumour predisposition (Genomics England PanelApp)",
      38: "Sarcoma cancer susceptibility (Genomics England PanelApp)",
      39: "Sarcoma susceptbility (Genomics England PanelApp)",
      40: "Thyroid cancer pertinent cancer susceptibility (Genomics England PanelApp)",
      41: "Tumour predisposition - childhood onset (Genomics England PanelApp)",
      42: "Upper gastrointestinal cancer pertinent cancer susceptibility (Genomics England PanelApp)"
}

panels = '\n'.join([f'{k} = {GE_panels[k]}' for k in GE_panels]) # for displaying in help
