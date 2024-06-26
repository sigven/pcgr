#---- color_palette ----#
color_palette <- list()
for (c in c("pathogenicity",
            "oncogenicity",
            "clinical_evidence",
            "cancer_assoc",
            "gene_expression",
            "exp_increase",
            "exp_decrease",
            "cna_variant_class",
            "tier",
            "report_color",
            "bg_dark",
            "warning",
            "success",
            "none")) {
  color_palette[[c]] <- list()
  color_palette[[c]][["levels"]] <- c()
  color_palette[[c]][["values"]] <- c()

  if (c == "cancer_assoc") {

    color_palette[[c]][["breaks"]] <-
      c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
    color_palette[[c]][["values"]] <-
      c("#b8b8ba",
        "#deebf7",
        "#c6dbef",
        "#9ecae1",
        "#6baed6",
        "#4292c6",
        "#2171b5",
        "#08519c",
        "#08306b")
  }

  if (c == "gene_expression") {
    #rep[["config"]][["disease"]] <- list()
    color_palette[[c]][["breaks"]] <-
      c(0.5, 10, 1000)
    color_palette[[c]][["values"]] <-
      c("#b8b8ba",
        "#e8c4e3",
        "#c37dbd",
        "#9b3297")
  }

  if (c == "pathogenicity") {
    color_palette[[c]][["levels"]] <-
      c("Pathogenic", "Likely_Pathogenic",
        "VUS", "Likely_Benign", "Benign")
    color_palette[[c]][["values"]] <-
      c("#9E0142", "#D53E4F", "#2c313c",
        "#78C679", "#077009")
  }

  if (c == "oncogenicity") {
    color_palette[[c]][["levels"]] <-
      c("Oncogenic", "Likely_Oncogenic",
        "VUS", "Likely_Benign", "Benign")
    color_palette[[c]][["values"]] <-
      c("#9E0142", "#D53E4F", "#2c313c",
        "#78C679", "#077009")
  }
  if (c == "clinical_evidence") {
    color_palette[[c]][["levels"]] <-
      c("A: Validated",
        "A: FDA/NCCN/ELN guidelines",
        "B: Clinical evidence",
        "B1: Clinical evidence: late trials",
        "B2: Clinical evidence: early trials",
        "C: Case study",
        "D: Preclinical evidence",
        "E: Indirect evidence")
    color_palette[[c]][["values"]] <-
      c("#009E73", "#009E73", "#56B4E9", "#56B4E9",
        "#56B4E9", "#0072B2", "#E69F00", "#F0E442")
  }
  if (c == "tier") {
    color_palette[[c]][["levels"]] <-
      c("TIER1", "TIER2", "TIER3",
        "TIER4", "TIER5",
        "COL6", "COL7", "COL8",
        "COL9", "COL10", "COL11",
        "COL12", "COL13", "COL14",
        "COL15", "COL16", "COL17",
        "COL18", "COL19", "COL20",
        "COL21", "COL22", "COL23",
        "COL24", "COL25")
    color_palette[[c]][["values"]] <-
      c("#0073C2", "#EFC000", "#CD534C", "#7AA6DC",
        "#8F7700", "#003C67", "#3B3B3B", "#868686",
        "#593196", "#1B9E77", "#276419", "#D95F02",
        "#000000", "#E31A1C", "#E7298A", "#009999",
        "#003366", "#660033", "#4DAC26",
        "#F1B6DA", "#E69F00", "#F0E442",
        "#F4A582", "#0072B2", "#BABABA")
  }
  if (c == "report_color") {
    color_palette[[c]][["levels"]] <- c("tumor_control", "tumor_only")
    color_palette[[c]][["values"]] <- c("#9B3297", "#0073C2")
  }
  if (c == "cna_variant_class") {
    color_palette[[c]][["levels"]] <- c("gain", "homdel")
    color_palette[[c]][["values"]] <- c("#00a65a", "#CD534C")
  }
  if (c == "warning") {
    color_palette[[c]] <- "#ff7518"
  }
  if (c == "bg_dark") {
    color_palette[[c]] <- "#2c313c"
  }
  if (c == "none") {
    color_palette[[c]] <- "#868686"
  }
  if (c == "success") {
    color_palette[[c]] <- "#00a65a"
  }
  if (c == 'exp_increase'){
    color_palette[[c]][["breaks"]] <- c(92, 94, 96, 98)
    color_palette[[c]][["values"]] <-
      c("#caead2","#a4d9b3","#7dc995","#52b777","#00a65a")
  }

  if (c == 'exp_decrease'){
    color_palette[[c]][["breaks"]] <- c(2, 4, 6, 8)
    color_palette[[c]][["values"]] <-
      c("#CD534C","#dd766b","#ea978c","#f6b7ae","#ffd8d2")
  }

}

usethis::use_data(color_palette, overwrite = T)

biomarker_evidence <- list()
#-----evidence types---------#
biomarker_evidence[['types']] <-
  c("predictive","prognostic","diagnostic",
    "oncogenic","predisposing","functional")
usethis::use_data(
  biomarker_evidence, overwrite = T)

#-----input column names/types-----#
data_coltype_defs <- list()
data_coltype_defs[['cna_somatic_raw']] <- readr::cols_only(
  CHROM = readr::col_character(),
  SEGMENT_START = readr::col_double(),
  SEGMENT_END = readr::col_double(),
  VAR_ID = readr::col_character(),
  VARIANT_CLASS = readr::col_character(),
  N_MAJOR = readr::col_integer(),
  N_MINOR = readr::col_integer(),
  CHROMOSOME_ARM = readr::col_character(),
  CYTOBAND = readr::col_character(),
  EVENT_TYPE = readr::col_character(),
  TRANSCRIPT_OVERLAP_PERCENT = readr::col_number(),
  TRANSCRIPT_START = readr::col_double(),
  TRANSCRIPT_END = readr::col_double(),
  ENSEMBL_TRANSCRIPT_ID = readr::col_character(),
  ENSEMBL_GENE_ID = readr::col_character(),
  SYMBOL = readr::col_character(),
  ENTREZGENE = readr::col_character(),
  GENENAME = readr::col_character(),
  REFSEQ_TRANSCRIPT_ID = readr::col_character(),
  ACTIONABLE_GENE	= readr::col_logical(),
  TSG	= readr::col_logical(),
  TSG_SUPPORT	= readr::col_character(),
  TSG_RANK = readr::col_integer(),
  ONCOGENE	= readr::col_logical(),
  ONCOGENE_SUPPORT = readr::col_character(),
  ONCOGENE_RANK	= readr::col_logical(),
  SEGMENT_LENGTH_MB = readr::col_number(),
  BIOMARKER_MATCH = readr::col_character(),
  SAMPLE_ID = readr::col_character())
  #TPM = readr::col_number(),
  #TPM_GENE = readr::col_number())

data_coltype_defs[['snv_indel_somatic_raw']] <- readr::cols_only(
  CHROM = readr::col_character(),
  POS = readr::col_double(),
  REF = readr::col_character(),
  ALT = readr::col_character(),
  DP_TUMOR = readr::col_integer(),
  VAF_TUMOR = readr::col_number(),
  DP_CONTROL = readr::col_integer(),
  VAF_CONTROL = readr::col_number(),
  CALL_CONFIDENCE = readr::col_character(),
  GENOMIC_CHANGE = readr::col_character(),
  VAR_ID = readr::col_character(),
  VARIANT_CLASS = readr::col_character(),
  CONSEQUENCE = readr::col_character(),
  IMPACT = readr::col_character(),
  LOSS_OF_FUNCTION = readr::col_logical(),
  LOF_FILTER = readr::col_character(),
  SPLICE_DONOR_RELEVANT = readr::col_logical(),
  NULL_VARIANT = readr::col_logical(),
  CODING_STATUS = readr::col_character(),
  EXONIC_STATUS = readr::col_character(),
  PROTEIN_CHANGE = readr::col_character(),
  HGVSp_short = readr::col_character(),
  HGVSc = readr::col_character(),
  HGVSp = readr::col_character(),
  CDS_CHANGE = readr::col_character(),
  CDS_RELATIVE_POSITION = readr::col_character(),
  EXON = readr::col_character(),
  EXON_AFFECTED = readr::col_character(),
  MUTATION_HOTSPOT = readr::col_character(),
  MUTATION_HOTSPOT_CANCERTYPE = readr::col_character(),
  MUTATION_HOTSPOT_MATCH = readr::col_character(),
  BIOMARKER_MATCH = readr::col_character(),
  ONCOGENICITY = readr::col_character(),
  ONCOGENICITY_CODE = readr::col_character(),
  ONCOGENICITY_SCORE = readr::col_integer(),
  PFAM_DOMAIN = readr::col_character(),
  PFAM_DOMAIN_NAME = readr::col_character(),
  SYMBOL = readr::col_character(),
  ENTREZGENE = readr::col_character(),
  GENENAME = readr::col_character(),
  ENSEMBL_GENE_ID = readr::col_character(),
  ENSEMBL_TRANSCRIPT_ID = readr::col_character(),
  ENSEMBL_PROTEIN_ID = readr::col_character(),
  REFSEQ_TRANSCRIPT_ID = readr::col_character(),
  REFSEQ_PROTEIN_ID = readr::col_character(),
  UNIPROT_ACC = readr::col_character(),
  UNIPROT_ID = readr::col_character(),
  APPRIS = readr::col_character(),
  CCDS = readr::col_character(),
  CANONICAL = readr::col_character(),
  BIOTYPE = readr::col_character(),
  TRANSCRIPT_MANE_SELECT = readr::col_character(),
  TRANSCRIPT_MANE_PLUS_CLINICAL = readr::col_character(),
  TSG = readr::col_logical(),
  TSG_RANK = readr::col_integer(),
  TSG_SUPPORT = readr::col_character(),
  ONCOGENE = readr::col_logical(),
  ONCOGENE_RANK = readr::col_integer(),
  ONCOGENE_SUPPORT = readr::col_character(),
  INTOGEN_DRIVER = readr::col_character(),
  INTOGEN_ROLE = readr::col_character(),
  CGC_SOMATIC = readr::col_logical(),
  CGC_GERMLINE = readr::col_logical(),
  CGC_TIER = readr::col_character(),
  #TPM = readr::col_number(),
  #TPM_GENE = readr::col_number(),
  #consTPM = readr::col_number(),
  CLINVAR_TRAITS_ALL = readr::col_character(),
  CLINVAR_MSID = readr::col_character(),
  CLINVAR_CLNSIG = readr::col_character(),
  CLINVAR_CLASSIFICATION = readr::col_character(),
  CLINVAR_CONFLICTED = readr::col_logical(),
  CLINVAR_REVIEW_STATUS_STARS = readr::col_integer(),
  CLINVAR_NUM_SUBMITTERS = readr::col_integer(),
  CLINVAR_VARIANT_ORIGIN = readr::col_character(),
  PANEL_OF_NORMALS = readr::col_logical(),
  DBSNP_RSID = readr::col_character(),
  COSMIC_ID = readr::col_character(),
  TCGA_PANCANCER_COUNT = readr::col_integer(),
  TCGA_FREQUENCY = readr::col_character(),
  REGULATORY_ANNOTATION = readr::col_character(),
  TRANSCRIPTION_FACTORS = readr::col_character(),
  MOTIF_SCORE_CHANGE = readr::col_character(),
  MOTIF_NAME = readr::col_character(),
  RMSK_HIT = readr::col_character(),
  SIMPLEREPEATS_HIT = readr::col_logical(),
  WINMASKER_HIT = readr::col_logical(),
  VEP_ALL_CSQ = readr::col_character(),
  gnomADe_AF = readr::col_number(),
  gnomADe_AMR_AF = readr::col_number(),
  gnomADe_AFR_AF = readr::col_number(),
  gnomADe_EAS_AF = readr::col_number(),
  gnomADe_FIN_AF = readr::col_number(),
  gnomADe_ASJ_AF = readr::col_number(),
  gnomADe_OTH_AF = readr::col_number(),
  gnomADe_NFE_AF = readr::col_number(),
  gnomADe_SAS_AF = readr::col_number(),
  EFFECT_PREDICTIONS = readr::col_character(),
  SAMPLE_ID = readr::col_character(),
  VCF_SAMPLE_ID = readr::col_character(),
  GENOME_VERSION = readr::col_character()
)

data_coltype_defs[['snv_indel_germline_raw']] <- readr::cols_only(
  CHROM = readr::col_character(),
  POS = readr::col_double(),
  REF = readr::col_character(),
  ALT = readr::col_character(),
  GENOMIC_CHANGE = readr::col_character(),
  GENOTYPE = readr::col_character(),
  DP_CONTROL = readr::col_integer(),
  VAR_ID = readr::col_character(),
  VARIANT_CLASS = readr::col_character(),
  CONSEQUENCE = readr::col_character(),
  IMPACT = readr::col_character(),
  LOSS_OF_FUNCTION = readr::col_logical(),
  LOF_FILTER = readr::col_character(),
  SPLICE_DONOR_RELEVANT = readr::col_logical(),
  NULL_VARIANT = readr::col_logical(),
  CODING_STATUS = readr::col_character(),
  EXONIC_STATUS = readr::col_character(),
  PROTEIN_CHANGE = readr::col_character(),
  HGVSp_short = readr::col_character(),
  HGVSc = readr::col_character(),
  HGVSp = readr::col_character(),
  CDS_CHANGE = readr::col_character(),
  CDS_RELATIVE_POSITION = readr::col_character(),
  EXON = readr::col_character(),
  EXON_AFFECTED = readr::col_integer(),
  EXON_POSITION = readr::col_integer(),
  LAST_EXON = readr::col_logical(),
  LAST_INTRON = readr::col_logical(),
  INTRON_POSITION = readr::col_integer(),
  MUTATION_HOTSPOT = readr::col_character(),
  MUTATION_HOTSPOT_CANCERTYPE = readr::col_character(),
  MUTATION_HOTSPOT_MATCH = readr::col_character(),
  BIOMARKER_MATCH = readr::col_character(),
  GWAS_HIT = readr::col_character(),
  PFAM_DOMAIN = readr::col_character(),
  PFAM_DOMAIN_NAME = readr::col_character(),
  SYMBOL = readr::col_character(),
  ENTREZGENE = readr::col_character(),
  GENENAME = readr::col_character(),
  ENSEMBL_GENE_ID = readr::col_character(),
  ENSEMBL_TRANSCRIPT_ID = readr::col_character(),
  ENSEMBL_PROTEIN_ID = readr::col_character(),
  REFSEQ_TRANSCRIPT_ID = readr::col_character(),
  REFSEQ_PROTEIN_ID = readr::col_character(),
  UNIPROT_ACC = readr::col_character(),
  UNIPROT_ID = readr::col_character(),
  PRINCIPAL_ISOFORM_FLAG = readr::col_character(),
  CCDS = readr::col_character(),
  CANONICAL = readr::col_character(),
  BIOTYPE = readr::col_character(),
  TRANSCRIPT_MANE_SELECT = readr::col_character(),
  TRANSCRIPT_MANE_PLUS_CLINICAL = readr::col_character(),
  TSG = readr::col_logical(),
  TSG_RANK = readr::col_integer(),
  TSG_SUPPORT = readr::col_character(),
  ONCOGENE = readr::col_logical(),
  ONCOGENE_RANK = readr::col_integer(),
  ONCOGENE_SUPPORT = readr::col_character(),
  INTOGEN_DRIVER = readr::col_character(),
  INTOGEN_ROLE = readr::col_character(),
  CGC_SOMATIC = readr::col_logical(),
  CGC_GERMLINE = readr::col_logical(),
  CGC_TIER = readr::col_character(),
  CPG_SOURCE = readr::col_character(),
  GERP_SCORE = readr::col_number(),
  CLINVAR_TRAITS_ALL = readr::col_character(),
  CLINVAR_MSID = readr::col_character(),
  CLINVAR_CLNSIG = readr::col_character(),
  CLINVAR_UMLS_CUI = readr::col_character(),
  CLINVAR_CLASSIFICATION = readr::col_character(),
  CLINVAR_CONFLICTED = readr::col_logical(),
  CLINVAR_REVIEW_STATUS_STARS = readr::col_integer(),
  CLINVAR_NUM_SUBMITTERS = readr::col_integer(),
  CLINVAR_VARIANT_ORIGIN = readr::col_character(),
  PANEL_OF_NORMALS = readr::col_logical(),
  DBSNP_RSID = readr::col_character(),
  COSMIC_ID = readr::col_character(),
  TCGA_PANCANCER_COUNT = readr::col_integer(),
  TCGA_FREQUENCY = readr::col_character(),
  REGULATORY_ANNOTATION = readr::col_character(),
  TRANSCRIPTION_FACTORS = readr::col_character(),
  MOTIF_SCORE_CHANGE = readr::col_character(),
  MOTIF_NAME = readr::col_character(),
  RMSK_HIT = readr::col_character(),
  VEP_ALL_CSQ = readr::col_character(),
  DBNSFP_SIFT = readr::col_character(),
  DBNSFP_PROVEAN = readr::col_character(),
  DBNSFP_META_RNN = readr::col_character(),
  DBNSFP_FATHMM = readr::col_character(),
  DBNSFP_MUTATIONTASTER = readr::col_character(),
  DBNSFP_DEOGEN2 = readr::col_character(),
  DBNSFP_PRIMATEAI = readr::col_character(),
  DBNSFP_MUTATIONASSESSOR = readr::col_character(),
  DBNSFP_FATHMM_MKL = readr::col_character(),
  DBNSFP_M_CAP = readr::col_character(),
  DBNSFP_LIST_S2 = readr::col_character(),
  DBNSFP_BAYESDEL_ADDAF = readr::col_character(),
  DBNSFP_SPLICE_SITE_ADA = readr::col_character(),
  DBNSFP_SPLICE_SITE_RF = readr::col_character(),
  gnomADe_AF = readr::col_number(),
  gnomADe_AMR_AF = readr::col_number(),
  gnomADe_AFR_AF = readr::col_number(),
  gnomADe_EAS_AF = readr::col_number(),
  gnomADe_FIN_AF = readr::col_number(),
  gnomADe_ASJ_AF = readr::col_number(),
  gnomADe_OTH_AF = readr::col_number(),
  gnomADe_NFE_AF = readr::col_number(),
  gnomADe_SAS_AF = readr::col_number(),
  gnomADe_non_cancer_AC = readr::col_integer(),
  gnomADe_non_cancer_AN = readr::col_integer(),
  gnomADe_non_cancer_NHOMALT = readr::col_integer(),
  gnomADe_non_cancer_AF = readr::col_number(),
  gnomADe_non_cancer_AFR_AC = readr::col_integer(),
  gnomADe_non_cancer_AFR_AN = readr::col_integer(),
  gnomADe_non_cancer_AFR_NHOMALT = readr::col_integer(),
  gnomADe_non_cancer_AFR_AF = readr::col_number(),
  gnomADe_non_cancer_AMR_AC = readr::col_integer(),
  gnomADe_non_cancer_AMR_AN = readr::col_integer(),
  gnomADe_non_cancer_AMR_NHOMALT = readr::col_integer(),
  gnomADe_non_cancer_AMR_AF = readr::col_number(),
  gnomADe_non_cancer_NFE_AC = readr::col_integer(),
  gnomADe_non_cancer_NFE_AN = readr::col_integer(),
  gnomADe_non_cancer_NFE_NHOMALT = readr::col_integer(),
  gnomADe_non_cancer_NFE_AF = readr::col_number(),
  gnomADe_non_cancer_FIN_AC = readr::col_integer(),
  gnomADe_non_cancer_FIN_AN = readr::col_integer(),
  gnomADe_non_cancer_FIN_NHOMALT = readr::col_integer(),
  gnomADe_non_cancer_FIN_AF = readr::col_number(),
  gnomADe_non_cancer_ASJ_AC = readr::col_integer(),
  gnomADe_non_cancer_ASJ_AN = readr::col_integer(),
  gnomADe_non_cancer_ASJ_NHOMALT = readr::col_integer(),
  gnomADe_non_cancer_ASJ_AF = readr::col_number(),
  gnomADe_non_cancer_OTH_AC = readr::col_integer(),
  gnomADe_non_cancer_OTH_AN = readr::col_integer(),
  gnomADe_non_cancer_OTH_NHOMALT = readr::col_integer(),
  gnomADe_non_cancer_OTH_AF = readr::col_number(),
  gnomADe_non_cancer_SAS_AC = readr::col_integer(),
  gnomADe_non_cancer_SAS_AN = readr::col_integer(),
  gnomADe_non_cancer_SAS_NHOMALT = readr::col_integer(),
  gnomADe_non_cancer_SAS_AF = readr::col_number(),
  gnomADe_non_cancer_EAS_AC = readr::col_integer(),
  gnomADe_non_cancer_EAS_AN = readr::col_integer(),
  gnomADe_non_cancer_EAS_NHOMALT = readr::col_integer(),
  gnomADe_non_cancer_EAS_AF = readr::col_number(),
  DBMTS = readr::col_character(),
  EFFECT_PREDICTIONS = readr::col_character(),
  SAMPLE_ID = readr::col_character(),
  VCF_SAMPLE_ID = readr::col_character(),
  GENOME_VERSION = readr::col_character()
)


usethis::use_data(data_coltype_defs, overwrite = T)

tsv_cols <- list()
tsv_cols[['cna']] <-
  c('SAMPLE_ID',
    'VAR_ID',
    'GENOME_VERSION',
    'CN_MAJOR',
    'CN_MINOR',
    'SEGMENT_LENGTH_MB',
    'CYTOBAND',
    'EVENT_TYPE',
    'VARIANT_CLASS',
    'SYMBOL',
    'ENTREZGENE',
    'GENENAME',
    'ENSEMBL_GENE_ID',
    'GENENAME',
    'TUMOR_SUPPRESSOR',
    'TUMOR_SUPPRESSOR_SUPPORT',
    'ONCOGENE',
    'ONCOGENE_SUPPORT',
    'TRANSCRIPT_OVERLAP',
    'ACTIONABILITY_TIER',
    'ACTIONABILITY',
    'ACTIONABILITY_FRAMEWORK',
    'BIOMARKER_MATCH',
    'TARGETED_INHIBITORS_ALL2')


tsv_cols[['snv_indel']] <-
  c('SAMPLE_ID',
    'GENOMIC_CHANGE',
    'GENOME_VERSION',
    'VARIANT_CLASS',
    'SYMBOL',
    'ENTREZGENE',
    'ENSEMBL_GENE_ID',
    'GENENAME',
    'PROTEIN_CHANGE',
    'CONSEQUENCE',
    'PFAM_DOMAIN_NAME',
    'LOSS_OF_FUNCTION',
    'LOF_FILTER',
    'CDS_CHANGE',
    'CODING_STATUS',
    'EXONIC_STATUS',
    'DP_TUMOR',
    'VAF_TUMOR',
    'DP_CONTROL',
    'VAF_CONTROL',
    'MUTATION_HOTSPOT',
    'MUTATION_HOTSPOT_CANCERTYPE',
    'ACTIONABILITY_TIER',
    'ACTIONABILITY',
    'ACTIONABILITY_FRAMEWORK',
    'ONCOGENICITY',
    'ONCOGENICITY_CODE',
    'ONCOGENICITY_SCORE',
    'HGVSc',
    'HGVSp',
    'CANONICAL',
    'CCDS',
    'UNIPROT_ACC',
    'ENSEMBL_TRANSCRIPT_ID',
    'ENSEMBL_PROTEIN_ID',
    'REFSEQ_TRANSCRIPT_ID',
    'REFSEQ_PROTEIN_ID',
    'TRANSCRIPT_MANE_SELECT',
    'TRANSCRIPT_MANE_PLUS_CLINICAL',
    'CGC_TIER',
    'CGC_GERMLINE',
    'CGC_SOMATIC',
    'ONCOGENE',
    'ONCOGENE_SUPPORT',
    'TUMOR_SUPPRESSOR',
    'TUMOR_SUPPRESSOR_SUPPORT',
    'TARGETED_INHIBITORS2',
    'EFFECT_PREDICTIONS',
    'REGULATORY_ANNOTATION',
    'VEP_ALL_CSQ',
    'gnomADe_AF',
    'DBSNP_RSID',
    'COSMIC_ID',
    'TCGA_FREQUENCY',
    'TCGA_PANCANCER_COUNT',
    'CLINVAR_MSID',
    'CLINVAR_CLASSIFICATION',
    'CLINVAR_VARIANT_ORIGIN',
    'CLINVAR_NUM_SUBMITTERS',
    'CLINVAR_REVIEW_STATUS_STARS',
    'CLINVAR_CONFLICTED',
    'BIOMARKER_MATCH',
    'CALL_CONFIDENCE'
    )

usethis::use_data(tsv_cols, overwrite = T)

dt_display <- list()
dt_display[['cna_gene_actionable']] <-
  c("SYMBOL",
    "GENENAME",
    "VARIANT_CLASS",
    "BIOMARKER_EVIDENCE",
    "EVENT_TYPE",
    "CN_TOTAL",
    "CYTOBAND",
    "ENSEMBL_GENE_ID",
    "CANCERGENE_EVIDENCE",
    "MOLECULAR_ALTERATION",
    "ONCOGENE",
    "TUMOR_SUPPRESSOR",
    "TRANSCRIPT_OVERLAP",
    "TISSUE_ASSOC_RANK",
    "SEGMENT",
    "SEGMENT_LENGTH_MB",
    "CN_MAJOR",
    "CN_MINOR",
    "TARGETED_INHIBITORS_ALL",
    "GENOME_VERSION")

dt_display[['cna_eitem']] <-
  c("MOLECULAR_ALTERATION",
    "CN_TOTAL",
    "BM_CANCER_TYPE",
    "BM_EVIDENCE_LEVEL",
    "BM_CONTEXT",
    "BM_MOLECULAR_PROFILE",
    "BM_REFERENCE",
    "BM_EVIDENCE_TYPE",
    "BM_CLINICAL_SIGNIFICANCE",
    "BM_THERAPEUTIC_CONTEXT",
    "BM_SOURCE_DB",
    "BM_RATING",
    "BM_EVIDENCE_ID",
    "BM_EVIDENCE_DESCRIPTION",
    "BM_EVIDENCE_DIRECTION",
    "BM_DISEASE_ONTOLOGY_ID",
    "BM_PRIMARY_SITE",
    "BM_RESOLUTION"
  )

dt_display[['cna_other_oncogenic']] <-
  c("SYMBOL",
    "GENENAME",
    "VARIANT_CLASS",
    "EVENT_TYPE",
    "CN_TOTAL",
    "CYTOBAND",
    "ENSEMBL_GENE_ID",
    "CANCERGENE_EVIDENCE",
    "ONCOGENE",
    "TUMOR_SUPPRESSOR",
    "TRANSCRIPT_OVERLAP",
    "GLOBAL_ASSOC_RANK",
    "TISSUE_ASSOC_RANK",
    "SEGMENT",
    "SEGMENT_LENGTH_MB",
    "CN_MAJOR",
    "CN_MINOR",
    "TARGETED_INHIBITORS_ALL",
    "GENOME_VERSION")

dt_display[['snv_indel_germline_filtered']] <-
  c('GENOMIC_CHANGE',
    'VARIANT_CLASS',
    'EXCLUSION_CRITERIA',
    'SYMBOL',
    'CONSEQUENCE',
    'gnomADe_AF',
    'DBSNP_RSID',
    'COSMIC_ID',
    'TCGA_FREQUENCY',
    'VAF_TUMOR',
    'DP_TUMOR',
    'VAF_CONTROL',
    'DP_CONTROL',
    'gnomADe_AFR_AF',
    'gnomADe_AMR_AF',
    'gnomADe_ASJ_AF',
    'gnomADe_EAS_AF',
    'gnomADe_FIN_AF',
    'gnomADe_NFE_AF',
    'gnomADe_OTH_AF',
    'gnomADe_SAS_AF')

dt_display[['snv_indel_gene_actionable']] <-
  c('MOLECULAR_ALTERATION',
    'GENENAME',
    'VAF_TUMOR',
    'BIOMARKER_EVIDENCE',
    'PROTEIN_DOMAIN',
    'CDS_CHANGE',
    'BM_TOP_RESOLUTION',
    'MUTATION_HOTSPOT',
    'MUTATION_HOTSPOT_CANCERTYPE',
    'TCGA_FREQUENCY',
    'CONSEQUENCE',
    'VARIANT_CLASS',
    'SYMBOL',
    'HGVSc',
    'HGVSp',
    'PREDICTED_EFFECT',
    'LOSS_OF_FUNCTION',
    'LOF_FILTER',
    'ONCOGENICITY',
    'ONCOGENICITY_CODE',
    'ONCOGENICITY_SCORE',
    'VEP_ALL_CSQ',
    'DBSNP_RSID',
    'COSMIC_ID',
    'CLINVAR',
    'ENSEMBL_GENE_ID',
    'ENSEMBL_TRANSCRIPT_ID',
    'REFSEQ_TRANSCRIPT_ID',
    'CANCERGENE_EVIDENCE',
    'TARGETED_INHIBITORS',
    'TARGETED_INHIBITORS_ALL',
    'CALL_CONFIDENCE',
    'DP_TUMOR',
    'DP_CONTROL',
    'VAF_CONTROL',
    'GENOMIC_CHANGE',
    'GENOME_VERSION')

dt_display[['snv_indel_eitem']] <-
  c('MOLECULAR_ALTERATION',
    'BM_CANCER_TYPE',
    'BM_EVIDENCE_LEVEL',
    'BM_CONTEXT',
    'BM_MOLECULAR_PROFILE',
    'BM_REFERENCE',
    'BM_EVIDENCE_DESCRIPTION',
    'BM_CLINICAL_SIGNIFICANCE',
    'BM_EVIDENCE_TYPE',
    'BM_THERAPEUTIC_CONTEXT',
    'BM_RATING',
    'BM_EVIDENCE_DIRECTION',
    'BM_SOURCE_DB',
    'BM_EVIDENCE_ID',
    'BM_DISEASE_ONTOLOGY_ID',
    'BM_PRIMARY_SITE',
    'BM_RESOLUTION',
    'BM_MATCH')


dt_display[['snv_indel_tier3']] <-
  c('SYMBOL',
    'PROTEIN_CHANGE',
    'GENENAME',
    'CONSEQUENCE',
    'ONCOGENICITY',
    'PROTEIN_DOMAIN',
    'MUTATION_HOTSPOT',
    'COSMIC_ID',
    'CDS_CHANGE',
    'HGVSc',
    'HGVSp',
    'MUTATION_HOTSPOT_CANCERTYPE',
    'LOSS_OF_FUNCTION',
    'LOF_FILTER',
    'TCGA_FREQUENCY',
    'PREDICTED_EFFECT',
    'ONCOGENICITY_CODE',
    'ONCOGENICITY_SCORE',
    'VEP_ALL_CSQ',
    'DBSNP_RSID',
    'CLINVAR',
    'ENSEMBL_TRANSCRIPT_ID',
    'REFSEQ_TRANSCRIPT_ID',
    'TARGETED_INHIBITORS',
    'TARGETED_INHIBITORS_ALL',
    'ONCOGENE',
    'TUMOR_SUPPRESSOR',
    'CANCERGENE_EVIDENCE',
    'CANCER_GENE_CENSUS',
    'GLOBAL_ASSOC_RANK',
    'TISSUE_ASSOC_RANK',
    'CALL_CONFIDENCE',
    'DP_TUMOR',
    'VAF_TUMOR',
    'DP_CONTROL',
    'VAF_CONTROL',
    'GENOMIC_CHANGE',
    'GENOME_VERSION')

dt_display[['tier4']] <-
  c('SYMBOL',
    'PROTEIN_CHANGE',
    'GENENAME',
    'CONSEQUENCE',
    'ONCOGENICITY',
    'PROTEIN_DOMAIN',
    'COSMIC_ID',
    'CDS_CHANGE',
    'TCGA_FREQUENCY',
    'HGVSc',
    'HGVSp',
    'PREDICTED_EFFECT',
    'LOSS_OF_FUNCTION',
    'LOF_FILTER',
    'REGULATORY_ANNOTATION',
    'ONCOGENICITY_CODE',
    'ONCOGENICITY_SCORE',
    'VEP_ALL_CSQ',
    'DBSNP_RSID',
    'ENSEMBL_GENE_ID',
    'ENSEMBL_TRANSCRIPT_ID',
    'ENSEMBL_PROTEIN_ID',
    'REFSEQ_TRANSCRIPT_ID',
    'CANCERGENE_EVIDENCE',
    'CANCER_GENE_CENSUS',
    'GLOBAL_ASSOC_RANK',
    'TISSUE_ASSOC_RANK',
    'CLINVAR',
    'TARGETED_INHIBITORS',
    'TARGETED_INHIBITORS_ALL',
    'CALL_CONFIDENCE',
    'DP_TUMOR',
    'VAF_TUMOR',
    'DP_CONTROL',
    'VAF_CONTROL',
    'GENOMIC_CHANGE',
    'GENOME_VERSION')

# dt_display[['tier5']] <-
#   c('SYMBOL',
#     'GENENAME',
#     'CONSEQUENCE',
#     'COSMIC_ID',
#     'TCGA_FREQUENCY',
#     'DBSNP_RSID',
#     'CLINVAR',
#     'ENSEMBL_TRANSCRIPT_ID',
#     'REFSEQ_TRANSCRIPT_ID',
#     'CANCERGENE_EVIDENCE',
#     'CANCER_GENE_CENSUS',
#     'GLOBAL_ASSOC_RANK',
#     'TISSUE_ASSOC_RANK',
#     'REGULATORY_ANNOTATION',
#     'VEP_ALL_CSQ',
#     'CALL_CONFIDENCE',
#     'DP_TUMOR',
#     'VAF_TUMOR',
#     'DP_CONTROL',
#     'VAF_CONTROL',
#     'GENOMIC_CHANGE',
#     'GENOME_VERSION')
#




usethis::use_data(dt_display, overwrite = T)

#---- variant_db_url ----#
variant_db_url <-
  dplyr::bind_rows(
    data.frame(
      name = "DBSNP_RSID",
      group_by_var = "VAR_ID",
      url_prefix = "https://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=",
      link_key_var = "DBSNP_RSID_RAW",
      rename = TRUE,
      link_display_var = "DBSNP_RSID_RAW"),
    data.frame(
      name = "CLINVAR",
      group_by_var = "VAR_ID",
      url_prefix = "http://www.ncbi.nlm.nih.gov/clinvar/variation/",
      link_key_var = "CLINVAR_MSID",
      rename = FALSE,
      link_display_var = "CLINVAR_TRAITS_ALL"),
    data.frame(
      name = "GENENAME",
      group_by_var = "VAR_ID",
      url_prefix = "https://www.ncbi.nlm.nih.gov/gene/",
      link_key_var = "ENTREZGENE",
      rename = TRUE,
      link_display_var = "GENENAME_RAW"),
    data.frame(
      name = "PROTEIN_DOMAIN",
      group_by_var = "VAR_ID",
      url_prefix = "https://www.ebi.ac.uk/interpro/entry/pfam/",
      link_key_var = "PFAM_DOMAIN",
      rename = FALSE,
      link_display_var = "PFAM_DOMAIN_NAME",
      stringsAsFactors = F),
    data.frame(
      name = "COSMIC_ID",
      group_by_var = "VAR_ID",
      url_prefix = "https://cancer.sanger.ac.uk/cosmic/search?q=",
      link_key_var = "COSMIC_ID_RAW",
      rename = TRUE,
      link_display_var = "COSMIC_ID_RAW",
      stringsAsFactors = F),
    data.frame(
      name = "REFSEQ_TRANSCRIPT_ID",
      group_by_var = "VAR_ID",
      url_prefix = "https://www.ncbi.nlm.nih.gov/nuccore/",
      link_key_var = "REFSEQ_TRANSCRIPT_ID_RAW",
      rename = TRUE,
      link_display_var = "REFSEQ_TRANSCRIPT_ID_RAW",
      stringsAsFactors = F),
    data.frame(
      name = "ENSEMBL_TRANSCRIPT_ID",
      group_by_var = "VAR_ID",
      url_prefix = "https://www.ensembl.org/Homo_sapiens/Transcript/Summary?t=",
      link_key_var = "ENSEMBL_TRANSCRIPT_ID_RAW",
      rename = TRUE,
      link_display_var = "ENSEMBL_TRANSCRIPT_ID_RAW",
      stringsAsFactors = F),
    data.frame(
      name = "ENSEMBL_GENE_ID",
      group_by_var = "VAR_ID",
      url_prefix = "https://www.ensembl.org/Homo_sapiens/Gene/Summary?g=",
      link_key_var = "ENSEMBL_GENE_ID_RAW",
      rename = TRUE,
      link_display_var = "ENSEMBL_GENE_ID_RAW",
      stringsAsFactors = F)
  )

usethis::use_data(variant_db_url, overwrite = T)

#---- effect_prediction_algos ----#
effect_prediction_algos <-
  utils::read.table(file = "data-raw/effect_prediction_algorithms.tsv",
             header = T, sep = "\t", quote = "", stringsAsFactors = F)
usethis::use_data(effect_prediction_algos, overwrite = T)

#---- cancer_phenotypes_regex ----#
cancer_phenotypes_regex <-
  paste0("(cancer|carcinoma|tumor|neoplasm|gangliom|lymphom|leukem|meningiom",
         "|blastoma|melanom|chordom|adenoma|sarcom|mesotheli|ependymom|glioma",
         "|neurofibro|keratoacan|nevus|brca|polyposis|myelodysplastic|cowden",
         "|gardner|noonan|fanconi|carney|bullosa|schwanno|li-fraumeni|xeroderma",
         "|leiomyom|muir-|nijmegen|neoplasia|trichoepithelioma|brooke|turcot",
         "|exostos|lynch|drash|wilm|perlman|fibrofolliculomas|hippel|hamartom",
         "|bloom|werner|peutz|tuberous|angiomyolipoma",
         "|lymphoproliferative|stat3|teratoma|thrombocytop|tp63|weaver",
         "|pheochromo|gorlin|telangiectasia|hemangiom|osteochondro|response",
         "|polg-related|ras-associated|dyskeratosis",
         "|waardenburg|beckwidth|birt-hogg|diamond|frasier",
         "|infantile myofibromatosis|proteus|rothmund|russel)")
usethis::use_data(cancer_phenotypes_regex, overwrite = T)


tcga_cohorts <- as.data.frame(
  TCGAbiolinks::getGDCprojects() |>
  dplyr::filter(
    stringr::str_detect(project_id,"TCGA")
  ) |>
  dplyr::select(
    tumor, name
  ) |>
  dplyr::rename(
    tcga_cancer_code = tumor,
    tcga_cancer_name = name
  )
)

usethis::use_data(tcga_cohorts, overwrite = T)

####---- COSMIC signatures (SBS v3.4) ----####
cosmic_signatures <-
  utils::read.table(file = "data-raw/COSMIC_v3.4_SBS_GRCh38.txt",
             header = T, sep = "\t", quote = "", stringsAsFactors = F)

context_order <-
  c("A[C>A]A","A[C>A]C","A[C>A]G","A[C>A]T",
    "C[C>A]A","C[C>A]C","C[C>A]G","C[C>A]T",
    "G[C>A]A","G[C>A]C","G[C>A]G","G[C>A]T",
    "T[C>A]A","T[C>A]C","T[C>A]G","T[C>A]T",
    "A[C>G]A","A[C>G]C","A[C>G]G","A[C>G]T",
    "C[C>G]A","C[C>G]C","C[C>G]G","C[C>G]T",
    "G[C>G]A","G[C>G]C","G[C>G]G","G[C>G]T",
    "T[C>G]A","T[C>G]C","T[C>G]G","T[C>G]T",
    "A[C>T]A","A[C>T]C","A[C>T]G","A[C>T]T",
    "C[C>T]A","C[C>T]C","C[C>T]G","C[C>T]T",
    "G[C>T]A","G[C>T]C","G[C>T]G","G[C>T]T",
    "T[C>T]A","T[C>T]C","T[C>T]G","T[C>T]T",
    "A[T>A]A","A[T>A]C","A[T>A]G","A[T>A]T",
    "C[T>A]A","C[T>A]C","C[T>A]G","C[T>A]T",
    "G[T>A]A","G[T>A]C","G[T>A]G","G[T>A]T",
    "T[T>A]A","T[T>A]C","T[T>A]G","T[T>A]T",
    "A[T>C]A","A[T>C]C","A[T>C]G","A[T>C]T",
    "C[T>C]A","C[T>C]C","C[T>C]G","C[T>C]T",
    "G[T>C]A","G[T>C]C","G[T>C]G","G[T>C]T",
    "T[T>C]A","T[T>C]C","T[T>C]G","T[T>C]T",
    "A[T>G]A","A[T>G]C","A[T>G]G","A[T>G]T",
    "C[T>G]A","C[T>G]C","C[T>G]G","C[T>G]T",
    "G[T>G]A","G[T>G]C","G[T>G]G","G[T>G]T",
    "T[T>G]A","T[T>G]C","T[T>G]G","T[T>G]T")

cosmic_signatures$Type <- factor(
  cosmic_signatures$Type, levels = context_order)
cosmic_signatures <- cosmic_signatures |>
  dplyr::arrange(Type)
cosmic_signatures$Type <- NULL

artefact_signatures <-
  c("SBS27","SBS43", paste0("SBS",rep(46:60)),
    "SBS95")

cosmic_sbs_signatures_no_artefacts <-
  cosmic_signatures

for(sig in artefact_signatures){
  cosmic_sbs_signatures_no_artefacts[,sig] <- NULL
}


cosmic_sbs_signatures_all <- as.matrix(cosmic_signatures)
cosmic_sbs_signatures_no_artefacts <- as.matrix(cosmic_sbs_signatures_no_artefacts)

cosmic_sbs_signatures <- list()
cosmic_sbs_signatures[['all']] <- cosmic_sbs_signatures_all
cosmic_sbs_signatures[['no_artefacts']] <-
  cosmic_sbs_signatures_no_artefacts
usethis::use_data(cosmic_sbs_signatures, overwrite = T)

#usethis::use_data(cosmic_sbs_signatures_all, overwrite = T)
#usethis::use_data(cosmic_sbs_signatures_no_artefacts, overwrite = T)

# immune_celltypes <- as.data.frame(
#   immunedeconv::cell_type_map |>
#     dplyr::filter(method_dataset == "quantiseq") |>
#     dplyr::select(method_cell_type, cell_type) |>
#     dplyr::mutate(cell_type = dplyr::if_else(
#       !is.na(cell_type) &
#         cell_type == "uncharacterized cell",
#       "Uncharacterized cell",
#       as.character(cell_type)
#     )) |>
#     dplyr::mutate(cell_type = factor(
#       cell_type, levels = cell_type)) |>
#     dplyr::distinct()
# )

immune_celltypes2 <- data.frame(
  method_cell_type = c("B.cells","Macrophages.M1",
                       "Macrophages.M2","Monocytes",
                       "Neutrophils","NK.cells",
                       "T.cells.CD4","T.cells.CD8",
                       "Tregs","Dendritic.cells",
                       "Other"),
  cell_type = c("B cell","Macrophage M1",
                "Macrophage M2","Monocyte",
                "Neutrophil","NK cell",
                "T cell CD4+ (non-regulatory)",
                "T cell CD8+",
                "T cell regulatory (Tregs)",
                "Myeloid dendritic cell",
                "Uncharacterized cell")
)


#usethis::use_data(immune_celltypes, overwrite = T)

germline_filter_levels <-
  c("SOMATIC",
    "GERMLINE_DBSNP",
    "GERMLINE_CLINVAR",
    "GERMLINE_GNOMAD",
    "MULTIPLE FILTERS",
    "GERMLINE_HET",
    "GERMLINE_HOM",
    "GERMLINE_PON"
    )

usethis::use_data(germline_filter_levels, overwrite = T)

rm(cancer_phenotypes_regex,
   data_coltype_defs,
   effect_prediction_algos,
   tcga_cohorts,
   immune_celltypes,
   artefact_signatures,
   context_order,
   cosmic_signatures,
   cosmic_sbs_signatures_no_artefacts,
   cosmic_sbs_signatures_all,
   cosmic_sbs_signatures,
   germline_filter_levels,
   sig,
   biomarker_evidence,
   variant_db_url,
   color_palette,
   dt_display,
   tsv_cols,
   c)
