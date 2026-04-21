#!/usr/bin/env python

from pcgr._version import __version__

## Version - software and bundle
PCGR_VERSION = __version__
DB_VERSION = '20260417'  # database build version (date-based)

## Miscellaneous settings
NCBI_BUILD_MAF = 'GRCh38'
MAX_VARIANTS_FOR_REPORT = 500_000
#MAX_VARIANTS_FOR_REPORT = 80_000
CODING_EXOME_SIZE_MB = 34.0

## Mutational signature settings
RECOMMENDED_N_MUT_SIGNATURE = 200
MINIMUM_N_MUT_SIGNATURE = 100
MAX_SIGNATURE_PREVALENCE = 20

## GENCODE versions
GENCODE_VERSION = {'grch38': 49,'grch37': 19}

## vcfanno settings
VCFANNO_MAX_PROC = 15

## Default value for integer annotations when NA
NA_INTEGER = -99999
NA_FLOAT = -99999.99

## VEP settings/versions
VEP_VERSION = '115'
VEP_ASSEMBLY = {'grch38': 'GRCh38','grch37': 'GRCh37'}
VEP_MIN_FORKS = 1
VEP_MAX_FORKS = 8
VEP_MIN_BUFFER_SIZE = 50
VEP_MAX_BUFFER_SIZE = 30000
VEP_PICK_CRITERIA = ['mane_select','mane_plus_clinical','canonical','biotype','ccds','rank','tsl','appris','length']

## Gene expression comparative analysis resources
EXPRESSION_DB_SOURCES = ['tcga','depmap','treehouse']

## Autosomes, sex chromosomes and mitochondrial chromosome
AUTOSOMES = [str(x) for x in range(1,22)]
SEX_CHROMOSOMES = ['X','Y']

## Sample identifier length (max/min allowed)
SAMPLE_ID_MAX_LENGTH = 40
SAMPLE_ID_MIN_LENGTH = 3

## GnomAD subpopulation allele frequency tags (exomes/genomes)
GNOMAD_MAIN_EXOME_AF_TAGS = ['gnomADe_SAS_AF','gnomADe_NFE_AF','gnomADe_AFR_AF','gnomADe_AMR_AF','gnomADe_EAS_AF']
GNOMAD_MAIN_GENOME_AF_TAGS = ['gnomADg_SAS_AF','gnomADg_NFE_AF','gnomADg_AFR_AF','gnomADg_AMR_AF','gnomADg_EAS_AF']

## Classified germline variant input (from CPSR) - required columns
germline_input_required_cols = [
    'SAMPLE_ID',
    'VAR_ID',
    'GENOMIC_CHANGE',
    'VARIANT_CLASS',
    'GENOTYPE',
    'ALTERATION',
    'DP_CONTROL',
    'CPSR_CLASSIFICATION_SOURCE',
    'GENENAME',
    'ENTREZGENE',
    'ENSEMBL_GENE_ID',
    'ENSEMBL_TRANSCRIPT_ID',
    'HGVSc',
    'HGVSc_RefSeq',
    'HGVSp',
    'CONSEQUENCE',
    'CDS_CHANGE',
    'SYMBOL',
    'CODING_STATUS',
    'PFAM_DOMAIN',
    'PFAM_DOMAIN_NAME',
    'PROTEIN_CHANGE',
    'LOSS_OF_FUNCTION',
    'NULL_VARIANT',
    'DBSNP_RSID',
    'CLINVAR_MSID',
    'CLINVAR_CLASSIFICATION',
    'CLINVAR_VARIANT_ORIGIN',
    'CLINVAR_PHENOTYPE',
    'CLINVAR_CONFLICTED',
    'CLINVAR_REVIEW_STATUS_STARS',
    'CPSR_CLASSIFICATION',
    'CPSR_PATHOGENICITY_SCORE',
    'CPSR_ACMG_CODE',
    'FINAL_CLASSIFICATION'
]

## Primary tumor sites - PCGR
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

## Genomics England panels - cancer predisposition (PanelApp)
GE_panels = {
      0: "CPSR exploratory cancer predisposition panel (PanelApp genes / TCGA's germline study / Cancer Gene Census / Other)",
      1: "Adult solid tumours cancer susceptibility (GEP)",
      2: "Adult solid tumours for rare disease (GEP)",
      3: "Bladder cancer pertinent cancer susceptibility (GEP)",
      4: "Brain cancer pertinent cancer susceptibility (GEP)",
      5: "Breast cancer pertinent cancer susceptibility (GEP)",
      6: "Childhood solid tumours cancer susceptibility (GEP)",
      7: "Colorectal cancer pertinent cancer susceptibility (GEP)",
      8: "Endometrial cancer pertinent cancer susceptibility (GEP)",
      9: "Familial Tumours Syndromes of the central & peripheral Nervous system (GEP)",
      10: "Familial breast cancer (GEP)",
      11: "Familial melanoma (GEP)",
      12: "Familial prostate cancer (GEP)",
      13: "Familial rhabdomyosarcoma (GEP)",
      14: "GI tract tumours (GEP)",
      15: "Genodermatoses with malignancies (GEP)",
      16: "Haematological malignancies cancer susceptibility (GEP)",
      17: "Haematological malignancies for rare disease (GEP)",
      18: "Head and neck cancer pertinent cancer susceptibility (GEP)",
      19: "Inherited MMR deficiency (Lynch Syndrome) - GEP",
      20: "Inherited non-medullary thyroid cancer (GEP)",
      21: "Inherited ovarian cancer (without breast cancer) (GEP)",
      22: "Inherited pancreatic cancer (GEP)",
      23: "Inherited polyposis and early onset colorectal cancer (GEP)",
      24: "Inherited predisposition to acute myeloid leukaemia (AML) (GEP)",
      25: "Inherited susceptibility to acute lymphoblastoid leukaemia (ALL) (GEP)",
      26: "Inherited predisposition to GIST (GEP)",
      27: "Inherited renal cancer (GEP)",
      28: "Inherited phaeochromocytoma and paraganglioma (GEP)",
      29: "Melanoma pertinent cancer susceptibility (GEP)",
      30: "Multiple endocrine tumours (GEP)",
      31: "Multiple monogenic benign skin tumours (GEP)",
      32: "Neuroendocrine cancer pertinent cancer susceptibility (GEP)",
      33: "Neurofibromatosis Type 1 (GEP)",
      34: "Ovarian cancer pertinent cancer susceptibility (GEP)",
      35: "Parathyroid Cancer (GEP)",
      36: "Prostate cancer pertinent cancer susceptibility (GEP)",
      37: "Renal cancer pertinent cancer susceptibility (GEP)",
      38: "Rhabdoid tumour predisposition (GEP)",
      39: "Sarcoma cancer susceptibility (GEP)",
      40: "Sarcoma susceptibility (GEP)",
      41: "Thyroid cancer pertinent cancer susceptibility (GEP)",
      42: "Tumour predisposition - childhood onset (GEP)",
      43: "Upper gastrointestinal cancer pertinent cancer susceptibility (GEP)",
      44: "DNA repair genes pertinent cancer susceptibility (GEP)"
}

panels = '\n'.join([f'{k} = {GE_panels[k]}' for k in GE_panels]) # for displaying in help

## https://www.ensembl.org/info/genome/variation/prediction/predicted_data.html#consequences
VEP_consequence_rank = {
    'transcript_ablation': 1,
    'splice_acceptor_variant': 2,
    'splice_donor_variant': 3,
    'stop_gained': 4,
    'frameshift_variant': 5,
    'stop_lost': 6,
    'start_lost' : 7,
    'transcript_amplification': 8,
    'feature_elongation': 9,
    'feature_truncation': 10,
    'inframe_insertion': 11,
    'inframe_deletion': 12,
    'missense_variant': 13,
    'protein_altering_variant': 14,
    'splice_donor_5th_base_variant': 15,
    'splice_region_variant': 16,
    'splice_donor_region_variant': 17,
    'splice_polypyrimidine_tract_variant': 18,
    'incomplete_terminal_codon_variant': 19,
    'start_retained_variant': 20,
    'stop_retained_variant': 21,
    'synonymous_variant': 22,
    'coding_sequence_variant': 23,
    'mature_miRNA_variant': 24,
    '5_prime_UTR_variant': 25,
    '3_prime_UTR_variant': 26,
    'non_coding_transcript_exon_variant': 27,
    'intron_variant': 28,
    'NMD_transcript_variant': 29,
    'non_coding_transcript_variant': 30,
    'coding_transcript_variant': 31,
    'upstream_gene_variant': 32,
    'downstream_gene_variant': 33,
    'TFBS_ablation': 34,
    'TFBS_amplification': 35,
    'TF_binding_site_variant': 36,
    'regulatory_region_ablation': 37,
    'regulatory_region_amplification': 38,
    'regulatory_region_variant': 39,
    'intergenic_variant': 40,
    'sequence_variant': 41
}

## Regular expressions on the VEP CSQ (consequence) that targets different types of variants
CSQ_MISSENSE_PATTERN = r"^missense_variant"
CSQ_CODING_PATTERN = \
    r"(stop_(lost|gained)|start_lost|frameshift_|missense_|splice_(donor|acceptor)|protein_altering|inframe_)"
CSQ_CODING_SILENT_PATTERN = \
    r"(stop_(lost|gained)|start_lost|frameshift_|missense_|splice_(donor|acceptor)|protein_altering|inframe_|synonymous|(start|stop)_retained)"
CSQ_NULL_PATTERN = r"^(stop_gained|frameshift_)|&stop_gained"
CSQ_SPLICE_REGION_PATTERN = r"(splice_|intron_variant)"
CSQ_SPLICE_DONOR_PATTERN = \
    r"(splice_region_variant|splice_donor_variant|splice_donor_region_variant|splice_donor_5th_base_variant)"
CSQ_SPLICE_ACCEPTOR_PATTERN = \
    r"(splice_polypyrimidine_tract_variant|splice_acceptor_variant)"
CSQ_LOF_PATTERN = r"(stop_gained|frameshift|splice_acceptor_variant|splice_donor_variant|start_lost)"


## MaxEntScan thresholds (relaxed) for splice site disruption (donor/acceptor) - 
## percent drop(gain)
DONOR_REF_MIN_SCORE = 4
ACCEPTOR_REF_MIN_SCORE = 5
DONOR_DISRUPTION_MES_DROP_CUTOFF = -60
ACCEPTOR_DISRUPTION_MES_DROP_CUTOFF = -60

## TCGA tumor cohorts
DISEASE_COHORTS = ['ACC','BLCA','BRCA','CESC',
                'CHOL','COAD','DLBC','ESCA',
                'GBM','HNSC','KICH','KIRC',
                'KIRP','LAML','LGG','LIHC',
                'LUAD','LUSC','MESO','OV',
                'PAAD','PCPG','PRAD','READ',
                'SARC','SKCM','STAD','TGCT',
                'THCA','THYM','UCEC','UCS','UVM']

## Tumor site to TCGA cohort mapping
SITE_TO_TCGA_COHORT = {
    'Adrenal Gland': ['TCGA_ACC','TCGA_PCPG'],
    'Biliary Tract': ['TCGA_CHOL'],
    'Bladder/Urinary Tract': ['TCGA_BLCA'],
    'Breast': ['TCGA_BRCA'],
    'Cervix': ['TCGA_CESC'],
    'CNS/Brain': ['TCGA_GBM','TCGA_LGG'],
    'Colon/Rectum': ['TCGA_COAD','TCGA_READ'],
    'Esophagus/Stomach': ['TCGA_STAD'],
    'Eye': ['TCGA_UVM'],
    'Head and Neck': ['TCGA_HNSC'],
    'Kidney': ['TCGA_KICH','TCGA_KIRC','TCGA_KIRP'],
    'Liver': ['TCGA_LIHC'],
    'Lung': ['TCGA_LUAD','TCGA_LUSC'],
    'Lymphoid': ['TCGA_DLBC'],
    'Myeloid': ['TCGA_LAML'],
    'Ovary/Fallopian Tube': ['TCGA_OV'],
    'Pancreas': ['TCGA_PAAD'],
    'Peripheral Nervous System': ['TCGA_PCPG'],
    'Pleura': ['TCGA_MESO'],
    'Prostate': ['TCGA_PRAD'],
    'Skin': ['TCGA_SKCM'],
    'Soft Tissue': ['TCGA_SARC'],
    'Testis': ['TCGA_TGCT'],
    'Thymus': ['TCGA_THYM'],
    'Thyroid': ['TCGA_THCA'],
    'Uterus': ['TCGA_UCEC','TCGA_UCS']
}

# Tumor site to OncoTree cancer-type code mapping (fallback when no code is user-provided).
# Uses specific cancer-type codes rather than tissue/organ-level codes so that OncoKB
# biomarker queries return meaningful hits. Each code represents the most prevalent and
# best-annotated malignancy for that site.
SITE_TO_ONCOTREE = {
    'Adrenal Gland': 'ADRENAL_GLAND',           # Adrenocortical Carcinoma
    'Ampulla of Vater': 'AMPULLA_OF_VATER',      # Ampullary Carcinoma
    'Biliary Tract': 'BILIARY_TRACT',          
    'Bladder/Urinary Tract': 'BLADDER',  # Bladder Urothelial Carcinoma
    'Bone': 'BONE',                     
    'Breast': 'BREAST',                  
    'Cervix': 'CERVIX',                 
    'CNS/Brain': 'BRAIN',               
    'Colon/Rectum': 'COADREAD',       # Colorectal Adenocarcinoma
    'Esophagus/Stomach': 'STOMACH',      # Stomach Adenocarcinoma
    'Eye': 'EYE',                     
    'Head and Neck': 'HEAD_NECK',         
    'Kidney': 'KIDNEY',                
    'Liver': 'LIVER',                   # Hepatocellular Carcinoma
    'Lung': 'LUNG',                   # Lung Adenocarcinoma
    'Lymphoid': 'LYMPH',              
    'Myeloid': 'MYELOID',                 
    'Ovary/Fallopian Tube': 'OVARY',  
    'Pancreas': 'PANCREAS',               # Pancreatic Ductal Adenocarcinoma
    'Peripheral Nervous System': 'PNS',  # Pheochromocytoma
    'Peritoneum': 'PERITONEUM',             # Primary Serous Peritoneal Carcinoma
    'Pleura': 'PLEURA',                 # Mesothelioma
    'Prostate': 'PROSTATE',               # Prostate Adenocarcinoma
    'Skin': 'SKIN',                    # Melanoma
    'Soft Tissue': 'SOFT_TISSUE',            # Sarcoma NOS
    'Testis': 'TESTIS',                 # Testicular Germ Cell Tumor
    'Thymus': 'THYMUS',                 # Thymoma
    'Thyroid': 'THYROID',                # Papillary Thyroid Carcinoma
    'Uterus': 'UTERUS',                 # Uterine Endometrioid Carcinoma
    'Vulva/Vagina': 'VULVA_VAGINA',           # Vulvar Squamous Cell Carcinoma
}

## DBNSFP algorithm score to PCGR field name mapping
DBNSFP_ALGORITHMS = {
    'aloft': 'DBNSFP_ALOFT',
    'alphamissense': 'DBNSFP_ALPHA_MISSENSE',
    'bayesdel_addaf': 'DBNSFP_BAYESDEL_ADDAF',
    'cadd': 'DBNSFP_CADD',
    'clinpred': 'DBNSFP_CLINPRED',
    'deogen2': 'DBNSFP_DEOGEN2',
    'esm1b': 'DBNSFP_ESM1B',
    'fathmm_xf': 'DBNSFP_FATHMM_XF',
    'gerp_rs': 'DBNSFP_GERP',       
    'list_s2': 'DBNSFP_LIST_S2',
    'metarnn': 'DBNSFP_META_RNN',
    'mutationassessor': 'DBNSFP_MUTATIONASSESSOR',
    'mutationtaster': 'DBNSFP_MUTATIONTASTER',
    'mutformer': 'DBNSFP_MUTFORMER', 
    'mutpred': 'DBNSFP_MUTPRED',
    'm-cap': 'DBNSFP_M_CAP',
    'polyphen2_hvar': 'DBNSFP_POLYPHEN2_HVAR',
    'phactboost': 'DBNSFP_PHACTBOOST',
    'primateai': 'DBNSFP_PRIMATEAI', 
    'provean': 'DBNSFP_PROVEAN',
    'revel': 'DBNSFP_REVEL',
    'sift': 'DBNSFP_SIFT',
    'splice_site_ada': 'DBNSFP_SPLICE_SITE_ADA',
    'splice_site_rf': 'DBNSFP_SPLICE_SITE_RF',
    'vest4': 'DBNSFP_VEST4'
}

ONCOGENICITY = {
    'gnomAD_extremely_rare_AF': 0.00001,
    'gnomAD_common_AF': 0.01,
    'gnomAD_very_common_AF': 0.05,
    'insilico_pred_min_majority': 8,
    'insilico_pred_max_minority': 2,
}

## VICC/ClinGen oncogenicity scoring thresholds
ONCOGENICITY_THRESHOLDS = {
    'likely_benign_upper':    -1,
    'likely_benign_lower':    -6,
    'benign_upper':           -7,
    'likely_oncogenic_lower':  5,
    'likely_oncogenic_upper':  8,
    'oncogenic_lower':         9,
}

## Canonical display order for ONCOGENICITY_CODE — matches the oncogenicity.tsv row order.
## Any code not listed here sorts to the end (future-proofing).
ONCOGENICITY_CODE_ORDER = [
    'ONCG_OVS1',  'ONCG_OVS1_A', 'ONCG_OVS1_B',
    'ONCG_OS1',   'ONCG_OS2_A',  'ONCG_OS2_B',  'ONCG_OS3',
    'ONCG_OM1',   'ONCG_OM2',    'ONCG_OM3',    'ONCG_OM4',
    'ONCG_OP1',   'ONCG_OP3',    'ONCG_OP4',
    'ONCG_SBVS1', 'ONCG_SBS1',   'ONCG_SBS2_A', 'ONCG_SBS2_B',
    'ONCG_SBP1',  'ONCG_SBP2',
]

## OncoKB-derived oncogenicity criteria — injected in the post-hoc refinement pass.
## Scores follow VICC/ClinGen evidence weighting (very strong=8, strong=±4, moderate=±2).
OKB_ONCOGENICITY_CRITERIA = {
    'ONCG_OS2_A':  {'score':  4.0, 'pole': 'P'},  # ONCOGENICITY_OKB == "Oncogenic"
    'ONCG_OS2_B':  {'score':  2.0, 'pole': 'P'},  # ONCOGENICITY_OKB == "Likely Oncogenic"
    'ONCG_OVS1_A': {'score':  8.0, 'pole': 'P'},  # MUTATION_EFFECT_OKB == "Loss-of-function" in TSG
    'ONCG_OVS1_B': {'score':  4.0, 'pole': 'P'},  # MUTATION_EFFECT_OKB == "Likely Loss-of-function" in TSG
    'ONCG_SBS2_A': {'score': -4.0, 'pole': 'B'},  # ONCOGENICITY_OKB == "Neutral"
    'ONCG_SBS2_B': {'score': -2.0, 'pole': 'B'},  # ONCOGENICITY_OKB == "Likely Neutral"
}

# https://github.com/oncokb/oncokb-annotator/tree/master/data
ONCOKB_MAF_REQUIRED_COLS = {
    "VAR_ID",
    "Hugo_Symbol", "NCBI_Build",
    "Chromosome", "Start_Position", "End_Position",
    "Variant_Classification", "HGVSp", "HGVSp_Short",
    "HGVSg","Reference_Allele",
    "Tumor_Seq_Allele1", "Tumor_Seq_Allele2",
    "Tumor_Sample_Barcode"
}

ONCOKB_FUSION_REQUIRED_COLS = {"Tumor_Sample_Barcode", "Fusion"}
ONCOKB_CNA_REQUIRED_COLS = {"Tumor_Sample_Barcode", "Hugo_Symbol", "Copy_Number_Alteration"}

## OncoKB annotation columns: source name → renamed column in output TSV
## Same mapping applies for SNV/InDel, fusion and CNA annotator outputs.
ONCOKB_COLS = {
    'MUTATION_EFFECT':             'MUTATION_EFFECT_OKB',
    'MUTATION_EFFECT_CITATIONS':   'MUTATION_EFFECT_CITATIONS_OKB',
    'ONCOGENIC':                   'ONCOGENICITY_OKB',
    'MUTATION_EFFECT_DESCRIPTION': 'MUTATION_EFFECT_DESCRIPTION_OKB',
}

## Mapping from PCGR VARIANT_CLASS to OncoKB Copy_Number_Alteration values
VARIANT_CLASS_TO_OKB_CNA = {
    'amplification': 'Amplification',
    'homdel':        'Deletion',
    'hetdel':        'Loss',
    'hemdel':        'Loss',
}