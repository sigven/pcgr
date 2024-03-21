#---- color_palette ----#
color_palette <- list()
for (c in c("pathogenicity",
            "oncogenicity",
            "clinical_evidence",
            "cancer_assoc",
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
    #color_palette[[c]][["breaks"]] <-
    #  c(0.40, 0.55, 0.70, 0.85)
    #color_palette[[c]][["values"]] <-
    #  c("#b8b8ba", "#BDD7E7", "#6BAED6",
    #    "#3182BD", "#08519C")

    #rep[["config"]][["disease"]] <- list()
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
    color_palette[[c]][["values"]] <- c("#2780e3", "#593196")

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

}

usethis::use_data(color_palette, overwrite = T)

#-----evidence types---------#
evidence_types <- c("predictive","prognostic","diagnostic",
                    "oncogenic","predisposing","functional")
usethis::use_data(evidence_types, overwrite = T)

#-----evidence levels---------#
evidence_levels <- c("any","A_B","C_D_E")
usethis::use_data(evidence_levels, overwrite = T)


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
  SPLICE_DONOR_RELEVANT = readr::col_logical(),
  NULL_VARIANT = readr::col_logical(),
  CODING_STATUS = readr::col_character(),
  EXONIC_STATUS = readr::col_character(),
  PROTEIN_CHANGE = readr::col_character(),
  HGVSp_short = readr::col_character(),
  HGVSc = readr::col_character(),
  HGVSp = readr::col_character(),
  CDS_CHANGE = readr::col_character(),
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
  SPLICE_DONOR_RELEVANT = readr::col_logical(),
  NULL_VARIANT = readr::col_logical(),
  CODING_STATUS = readr::col_character(),
  EXONIC_STATUS = readr::col_character(),
  PROTEIN_CHANGE = readr::col_character(),
  HGVSp_short = readr::col_character(),
  HGVSc = readr::col_character(),
  HGVSp = readr::col_character(),
  CDS_CHANGE = readr::col_character(),
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
tsv_cols[['snv_indel']] <-
  c('GENOMIC_CHANGE',
    'VAR_ID',
    'GENOME_VERSION',
    'TIER',
    'TIER_DESCRIPTION',
    'TIER_FRAMEWORK',
    'SAMPLE_ID',
    'VARIANT_CLASS',
    'SYMBOL',
    'GENENAME',
    'PROTEIN_CHANGE',
    'CONSEQUENCE',
    'PFAM_DOMAIN_NAME',
    'LOSS_OF_FUNCTION',
    'CDS_CHANGE',
    'CODING_STATUS',
    'EXONIC_STATUS',
    'MUTATION_HOTSPOT',
    'MUTATION_HOTSPOT_CANCERTYPE',
    'ONCOGENICITY',
    'ONCOGENICITY_CODE',
    'ONCOGENICITY_SCORE',
    'HGVSc',
    'HGVSp',
    'ENTREZGENE',
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
    'TARGETED_INHIBITORS_ALL2',
    'PREDICTED_EFFECT',
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
    'CALL_CONFIDENCE',
    'DP_TUMOR',
    'VAF_TUMOR',
    'DP_CONTROL',
    'VAF_CONTROL')

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
    "BM_EVIDENCE_TYPE",
    "BM_REFERENCE",
    "BM_CLINICAL_SIGNIFICANCE",
    "BM_THERAPEUTIC_CONTEXT",
    "BM_MOLECULAR_PROFILE_NAME",
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
    'ONCOGENICITY',
    'ONCOGENICITY_CODE',
    'ONCOGENICITY_SCORE',
    'VEP_ALL_CSQ',
    'DBSNP_RSID',
    'COSMIC_ID',
    'CLINVAR',
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
    'BM_CLINICAL_SIGNIFICANCE',
    'BM_EVIDENCE_TYPE',
    'BM_THERAPEUTIC_CONTEXT',
    'BM_REFERENCE',
    'BM_MOLECULAR_PROFILE_NAME',
    'BM_RATING',
    'BM_EVIDENCE_DIRECTION',
    'BM_DESCRIPTION',
    'BM_SOURCE_DB',
    'BM_EVIDENCE_ID',
    'BM_DISEASE_ONTOLOGY_ID',
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
    'REGULATORY_ANNOTATION',
    'ONCOGENICITY_CODE',
    'ONCOGENICITY_SCORE',
    'VEP_ALL_CSQ',
    'DBSNP_RSID',
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

dt_display[['tier5']] <-
  c('SYMBOL',
    'GENENAME',
    'CONSEQUENCE',
    'COSMIC_ID',
    'TCGA_FREQUENCY',
    'DBSNP_RSID',
    'CLINVAR',
    'ENSEMBL_TRANSCRIPT_ID',
    'REFSEQ_TRANSCRIPT_ID',
    'CANCERGENE_EVIDENCE',
    'CANCER_GENE_CENSUS',
    'GLOBAL_ASSOC_RANK',
    'TISSUE_ASSOC_RANK',
    'REGULATORY_ANNOTATION',
    'VEP_ALL_CSQ',
    'CALL_CONFIDENCE',
    'DP_TUMOR',
    'VAF_TUMOR',
    'DP_CONTROL',
    'VAF_CONTROL',
    'GENOMIC_CHANGE',
    'GENOME_VERSION')





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

rm(cancer_phenotypes_regex,
   data_coltype_defs,
   effect_prediction_algos,
   tcga_cohorts,
   evidence_types,
   evidence_levels,
   variant_db_url,
   color_palette,
   dt_display,
   tsv_cols,
   c)
