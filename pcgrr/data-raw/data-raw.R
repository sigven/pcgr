####---- Tumor sites -----####

tumor_sites <-
  c('Any',
    'Adrenal Gland',
    'Ampulla of Vater',
    'Biliary Tract',
    'Bladder/Urinary Tract',
    'Bone',
    'Breast',
    'Cervix',
    'CNS/Brain',
    'Colon/Rectum',
    'Esophagus/Stomach',
    'Eye',
    'Head and Neck',
    'Kidney',
    'Liver',
    'Lung',
    'Lymphoid',
    'Myeloid',
    'Ovary/Fallopian Tube',
    'Pancreas',
    'Peripheral Nervous System',
    'Peritoneum',
    'Pleura',
    'Prostate',
    'Skin',
    'Soft Tissue',
    'Testis',
    'Thymus',
    'Thyroid',
    'Uterus',
    'Vulva/Vagina')

usethis::use_data(
  tumor_sites, overwrite = T)

####---- Color_palette ----####
color_palette <- list()
for (c in c("pathogenicity",
            "pathogenicity_onc",
            "oncogenicity",
            "clinical_evidence",
            "cancer_assoc",
            "expression_outlier_high",
            "expression_outlier_low",
            "cna_variant_class",
            "tier_sensitivity",
            "tier_resistance",
            "prognosis",
            "diagnosis",
            "multi",
            "elevel",
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

  if (c == "pathogenicity") {
    color_palette[[c]][["levels"]] <-
      c("Pathogenic",
        "Likely Pathogenic",
        "VUS",
        "Likely Benign",
        "Benign")
    color_palette[[c]][["values"]] <-
      c("#9E0142",
        "#D53E4F",
        "#2c313c",
        "#78C679",
        "#077009")
  }


  if (c == "tier_sensitivity"){
    color_palette[[c]][["levels"]] <-
      c(1, 2, 3)
    color_palette[[c]][["values"]] <-
      c("#1565C0", "#42A5F5", "#868686")
  }

  if (c == "tier_resistance"){
    color_palette[[c]][["levels"]] <-
      c(1, 2, 3)
    color_palette[[c]][["values"]] <-
      c("#C62828", "#E57373", "#868686")

  }

  if (c == "prognosis"){
    color_palette[[c]][["levels"]] <-
      c("Better Outcome", "Poor Outcome")
    color_palette[[c]][["values"]] <-
      c("#2E7D32","#E65100")
  }
  if (c == "diagnosis"){
    color_palette[[c]] <- "#6A1B9A"
  }


  if (c == "pathogenicity_onc") {
    color_palette[[c]][["levels"]] <-
      c("Pathogenic",
        "Likely Pathogenic",
        "VUS",
        "Likely Benign",
        "Benign")
    color_palette[[c]][["values"]] <-
      c("#FF8790",
        "#FF8790",
        "#2c313c",
        "#65a058",
        "#65a058")
  }


  if (c == "oncogenicity") {
    color_palette[[c]][["levels"]] <-
      c("Oncogenic",
        "Likely Oncogenic",
        "VUS",
        "Likely Benign",
        "Benign")
    color_palette[[c]][["values"]] <-
      c("#9E0142",
        "#D53E4F",
        "#2c313c",
        "#78C679",
        "#077009")
  }

  if(c == "elevel"){

    # Evidence level colors (warm gray sequential)
    color_palette[[c]] <- list(
      A = list(bg = "#37303A", fg = "#ffffff"),
      B = list(bg = "#5D5165", fg = "#ffffff"),
      C = list(bg = "#857490", fg = "#ffffff"),
      D = list(bg = "#AE9AB8", fg = "#1a1a1a"),
      E = list(bg = "#D4C5DC", fg = "#1a1a1a"),
      '1' = list(bg = "#37303A", fg = "#ffffff"),
      '2' = list(bg = "#5D5165", fg = "#ffffff"),
      R1 = list(bg = "#37303A", fg = "#ffffff"),
      R2 = list(bg = "#5D5165", fg = "#ffffff"),
      '3A' = list(bg = "#857490", fg = "#ffffff"),
      '3B' = list(bg = "#AE9AB8", fg = "#1a1a1a"),
      '4' = list(bg = "#D4C5DC", fg = "#1a1a1a")
    )
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
      c("#37303A","#37303A",
        "#5D5165","#5D5165",
        "#5D5165","#857490",
        "#AE9AB8","#D4C5DC")
  }
  if (c == "multi") {
    color_palette[[c]][["levels"]] <-
      c("COL1", "COL2", "COL3",
        "COL4", "COL5", "COL6",
        "COL7", "COL8", "COL9",
        "COL10", "COL11",
        "COL12", "COL13", "COL14",
        "COL15", "COL16", "COL17",
        "COL18", "COL19", "COL20",
        "COL21", "COL22", "COL23",
        "COL24", "COL25")
    color_palette[[c]][["values"]] <-
      c("#0073C2", "#EFC000", "#CD534C",
        "#7AA6DC", "#8F7700", "#003C67",
        "#3B3B3B", "#868686", "#593196",
        "#1B9E77", "#276419",
        "#D95F02", "#000000", "#E31A1C",
        "#E7298A", "#009999", "#003366",
        "#660033", "#4DAC26", "#F1B6DA",
        "#E69F00", "#F0E442", "#F4A582",
        "#0072B2", "#BABABA")
  }
  if (c == "report_color") {
    color_palette[[c]][["levels"]] <-
      c("tumor_control", "tumor_only")
    color_palette[[c]][["values"]] <-
      c("#9B3297",
        "#0073C2")
  }
  if (c == "cna_variant_class") {
    color_palette[[c]][["levels"]] <-
      c("amplification",
        "gain",
        "neutral",
        "hetdel",
        "hemdel",
        "homdel")
    color_palette[[c]][["levels_display"]] <-
      c("Amplification",
        "Gain",
        "No aberration/neutral",
        "Shallow deletion",
        "Deep deletion",
        "Deep deletion")
    color_palette[[c]][["values"]] <-
      c("#8B0000",
        "#DD3D2D",
        "#868686",
        "#2171B5",
        "#08306B",
        "#08306B")
  }
  if (c == "warning") {
    color_palette[[c]] <- "#E69F00"
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
  if (c == 'expression_outlier_high'){
    color_palette[[c]][["breaks"]] <- c(96,97,98,99)
    color_palette[[c]][["values"]] <-
      c("#FDB366",  # lightest orange
        "#F67E4B",  # medium orange
        "#DD3D2D",  # red-orange
        "#B81B1B",  # deep red
        "#8B0000")  # darkest red
  }

  # REVISED: Low expression outliers - cool blue gradient
  if (c == 'expression_outlier_low'){
    color_palette[[c]][["breaks"]] <- c(1,2,3,4)
    color_palette[[c]][["values"]] <-
      c("#08306B",  # darkest blue
        "#08519C",  # deep blue
        "#2171B5",  # medium blue
        "#4292C6",  # lighter blue
        "#6BAED6")  # lightest blue
  }

}


## Callout box styles for "no biomarker hits" messages.
## Each entry has border (accent), bg (light tint), and text (dark shade)
## derived from the same hue as the corresponding palette entry.
color_palette[["bm_callout"]] <- list(
  therapeutic_sensitivity = list(
    border = "#1565C0",   # = tier_sensitivity$values[1]
    bg     = "#e3f2fd",
    text   = "#0d47a1"
  ),
  therapeutic_resistance = list(
    border = "#C62828",   # = tier_resistance$values[1]
    bg     = "#ffebee",
    text   = "#b71c1c"
  ),
  prognostic_poor = list(
    border = "#E65100",   # = prognosis$values[2] (Poor Outcome)
    bg     = "#fff3e0",
    text   = "#bf360c"
  ),
  prognostic_better = list(
    border = "#2E7D32",   # = prognosis$values[1] (Better Outcome)
    bg     = "#e8f5e9",
    text   = "#1b5e20"
  ),
  diagnostic_positive = list(
    border = "#6A1B9A",   # = diagnosis
    bg     = "#f3e5f5",
    text   = "#4a148c"
  )
)

usethis::use_data(color_palette, overwrite = T)

bm_categories <- list()
bm_cats <-
  c("therapeutic_sensitivity",
    "prognostic_poor",
    "prognostic_better",
    "diagnostic_positive",
    "diagnostic_negative",
    "therapeutic_resistance")
bm_cats_etype <-
  c("predictive",
    "prognostic",
    "prognostic",
    "diagnostic",
    "diagnostic",
    "predictive")
bm_cats_clnsig <-
  c("Sensitivity/Response",
    "Poor Outcome",
    "Better Outcome",
    "Positive",
    "Negative")

i <- 1
for(c in bm_cats){
  bm_categories[[c]] <- list()
  if(c == "therapeutic_sensitivity" |
     c == "diagnostic_positive" |
     c == "diagnostic_negative" |
     c == "prognostic_poor" |
     c == "prognostic_better"){
    bm_categories[[c]]
    bm_categories[[c]][["etype"]] <- bm_cats_etype[i]
    bm_categories[[c]][["clnsig"]] <- bm_cats_clnsig[i]
  } else {
    bm_categories[[c]] <- list()
    bm_categories[[c]][["etype"]] <- bm_cats_etype[i]
    bm_categories[[c]][["clnsig"]] <-
      c("Resistance/Non-response",
        "Adverse Response",
        "Toxicity")

  }
  i <- i + 1
}

usethis::use_data(
  bm_categories, overwrite = T)

####--- Biomarker evidence defs ----####

bm_evidence <- list()
#-----evidence types---------#
bm_evidence[['types']] <-
  c("predictive",
    "prognostic",
    "diagnostic",
    "oncogenic",
    "predisposing",
    "functional")
#-----significance levels---------#
bm_evidence[['levels']] <-
  c("A: Validated",
    "A: FDA/NCCN/ELN guidelines",
    "B: Clinical evidence",
    "B1: Clinical evidence: late trials",
    "B2: Clinical evidence: early trials",
    "C: Case study",
    "D: Preclinical evidence",
    "E: Indirect evidence")
#-----clinical significance type----#
bm_evidence[['clinical_significance']] <-
  c("Resistance/Non-response",
    "Sensitivity/Response",
    "Adverse Response",
    "Better Outcome",
    "Poor Outcome",
    "Positive",
    "Negative",
    "Toxicity",
    "Uncertain Significance")
## CIViC strong: A, B; OncoKB therapeutic: 1, 2, 3A, 3B, R1, R2;
## OncoKB diagnostic/prognostic: Dx1, Px1 (highest confidence)
bm_evidence[['strong_regex']] <-
  "^(A|B|(R)?(1|2|(3(A|B)))|(Dx|Px)1)"
## CIViC weak: C, D, E; OncoKB therapeutic: 4;
## OncoKB diagnostic: Dx2, Dx3 (also matched via leading D from CIViC pattern);
## OncoKB prognostic: Px2, Px3 (not covered by CIViC letters, added explicitly)
## Note: Dx1/Px1 also match D/Px here but strong_regex is checked first in
## case_when so they are correctly assigned the strong tier.
bm_evidence[['weak_regex']] <-
  "^(C|D|E|4|Px[23])"

usethis::use_data(
  bm_evidence, overwrite = T)



####-----::Data type definitions::-----####
data_coltype_defs <- list()

####----RNA fusions ----####
data_coltype_defs[['rna_fusion_raw']] <- readr::cols_only(
  SAMPLE_ID = readr::col_character(),
  VAR_ID = readr::col_character(),
  VARIANT_CLASS = readr::col_character(),
  FUSION_GENE = readr::col_character(),
  BREAKPOINT_5P = readr::col_character(),
  BREAKPOINT_3P = readr::col_character(),
  SPLIT_READS = readr::col_integer(),
  VARIANT_SUMMARY_OKB = readr::col_character(),
  TUMOR_TYPE_SUMMARY_OKB = readr::col_character(),
  MUTATION_EFFECT_OKB = readr::col_character(),
  MUTATION_EFFECT_DESCRIPTION_OKB = readr::col_character(),
  MUTATION_EFFECT_CITATIONS_OKB = readr::col_character(),
  ONCOGENICITY_OKB = readr::col_character(),
  BIOMARKER_MATCH = readr::col_character()
)


####---- CNA segments----#####
data_coltype_defs[['cna_somatic_segment_raw']] <- readr::cols_only(
  CHROM = readr::col_character(),
  SEGMENT_START = readr::col_double(),
  SEGMENT_END = readr::col_double(),
  SEGMENT_NAME = readr::col_character()
)

####----CNA genes ----####
data_coltype_defs[['cna_somatic_gene_raw']] <- readr::cols_only(
  CHROM = readr::col_character(),
  SEGMENT_START = readr::col_double(),
  SEGMENT_END = readr::col_double(),
  VAR_ID = readr::col_character(),
  VARIANT_CLASS = readr::col_character(),
  CN_MAJOR = readr::col_integer(),
  CN_MINOR = readr::col_integer(),
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
  ONCOGENE_RANK	= readr::col_integer(),
  SEGMENT_LENGTH_MB = readr::col_number(),
  FOLD_CHANGE = readr::col_number(),
  LOH = readr::col_character(),
  TUMOR_PLOIDY = readr::col_number(),
  TUMOR_PLOIDY_SOURCE = readr::col_character(),
  TUMOR_PURITY = readr::col_number(),
  TWOHIT_CANDIDATE_SOMATIC = readr::col_character(),
  TWOHIT_CANDIDATE_GERMLINE = readr::col_character(),
  BIOMARKER_MATCH = readr::col_character(),
  MUTATION_EFFECT_OKB = readr::col_character(),
  MUTATION_EFFECT_DESCRIPTION_OKB = readr::col_character(),
  MUTATION_EFFECT_CITATIONS_OKB = readr::col_character(),
  ONCOGENICITY_OKB = readr::col_character(),
  TUMOR_TYPE_SUMMARY_OKB = readr::col_character(),
  VARIANT_SUMMARY_OKB = readr::col_character(),
  SAMPLE_ID = readr::col_character())
  #TPM = readr::col_number(),
  #TPM_GENE = readr::col_number())

####----SNV/Indel ----####
data_coltype_defs[['snv_indel_somatic_raw']] <- readr::cols_only(
  CHROM = readr::col_character(),
  POS = readr::col_double(),
  REF = readr::col_character(),
  ALT = readr::col_character(),
  DP_TUMOR = readr::col_integer(),
  VAF_TUMOR = readr::col_number(),
  AD_TUMOR = readr::col_integer(),
  DP_CONTROL = readr::col_integer(),
  VAF_CONTROL = readr::col_number(),
  AD_CONTROL = readr::col_integer(),
  CALL_CONFIDENCE = readr::col_character(),
  GENOMIC_CHANGE = readr::col_character(),
  VAR_ID = readr::col_character(),
  VARIANT_CLASS = readr::col_character(),
  CONSEQUENCE = readr::col_character(),
  IMPACT = readr::col_character(),
  LOSS_OF_FUNCTION = readr::col_logical(),
  LOF_FILTER = readr::col_character(),
  MAXENTSCAN = readr::col_character(),
  MAXENTSCAN2 = readr::col_character(),
  MAXENTSCAN_REF = readr::col_number(),
  MAXENTSCAN_ALT = readr::col_number(),
  MAXENTSCAN_DG = readr::col_number(), #delta gain
  MAXENTSCAN_DL = readr::col_number(), #delta loss
  MAXENTSCAN_DIFF = readr::col_number(),
  MAXENTSCAN_PCT_CHANGE = readr::col_number(),
  NMD = readr::col_character(),
  SPLICE_DONOR_RELEVANT = readr::col_logical(),
  SPLICE_EFFECT_MUTSPLICEDB = readr::col_character(),
  NULL_VARIANT = readr::col_logical(),
  CODING_STATUS = readr::col_character(),
  EXONIC_STATUS = readr::col_character(),
  ALTERATION = readr::col_character(),
  PROTEIN_CHANGE = readr::col_character(),
  GRANTHAM_DISTANCE = readr::col_integer(),
  HGVSp_short = readr::col_character(),
  HGVSc = readr::col_character(),
  HGVSc_RefSeq = readr::col_character(),
  HGVSp = readr::col_character(),
  CODON = readr::col_character(),
  CDS_CHANGE = readr::col_character(),
  EXON_POSITION = readr::col_integer(),
  INTRON_POSITION = readr::col_integer(),
  EXON_INTRON_JUNCTION_SPAN = readr::col_logical(),
  CDS_RELATIVE_POSITION = readr::col_number(),
  PROTEIN_RELATIVE_POSITION = readr::col_number(),
  EXON = readr::col_character(),
  EXON_AFFECTED = readr::col_character(),
  MUTATION_HOTSPOT = readr::col_character(),
  MUTATION_HOTSPOT_CANCERTYPE = readr::col_character(),
  MUTATION_HOTSPOT_MATCH = readr::col_character(),
  BIOMARKER_MATCH = readr::col_character(),
  ONCOGENICITY = readr::col_character(),
  ONCOGENICITY_CODE = readr::col_character(),
  ONCOGENICITY_SCORE = readr::col_integer(),
  KNOWN_ONCOGENIC = readr::col_character(),
  KNOWN_ONCOGENIC_SITE = readr::col_character(),
  MUTATION_EFFECT_OKB = readr::col_character(),
  MUTATION_EFFECT_DESCRIPTION_OKB = readr::col_character(),
  MUTATION_EFFECT_CITATIONS_OKB = readr::col_character(),
  ONCOGENICITY_OKB = readr::col_character(),
  TUMOR_TYPE_SUMMARY_OKB = readr::col_character(),
  VARIANT_SUMMARY_OKB = readr::col_character(),
  PFAM_DOMAIN = readr::col_character(),
  PFAM_DOMAIN_NAME = readr::col_character(),
  PFAM_ENTRY_LOCATIONS = readr::col_character(),
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
  MANE_SELECT = readr::col_character(),
  MANE_SELECT2 = readr::col_character(),
  MANE_PLUS_CLINICAL = readr::col_character(),
  MANE_PLUS_CLINICAL2 = readr::col_character(),
  TSG = readr::col_logical(),
  TSG_RANK = readr::col_integer(),
  TSG_SUPPORT = readr::col_character(),
  ONCOGENE = readr::col_logical(),
  ONCOGENE_RANK = readr::col_integer(),
  ONCOGENE_SUPPORT = readr::col_character(),
  INTOGEN_DRIVER = readr::col_character(),
  INTOGEN_ROLE = readr::col_character(),
  #TPM = readr::col_number(),
  #TPM_GENE = readr::col_number(),
  #consTPM = readr::col_number(),
  CLINVAR_TRAITS_ALL = readr::col_character(),
  CLINVAR_MSID = readr::col_character(),
  CLINVAR_CLNSIG = readr::col_character(),
  CLINVAR_CLASSIFICATION = readr::col_character(),
  CLINVAR_CONFLICTED = readr::col_logical(),
  CLINVAR_GOLD_STARS = readr::col_integer(),
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
  gnomADe_REMAINING_AF = readr::col_number(),
  gnomADe_NFE_AF = readr::col_number(),
  gnomADe_SAS_AF = readr::col_number(),
  gnomADg_AF = readr::col_number(),
  gnomADg_AMR_AF = readr::col_number(),
  gnomADg_AFR_AF = readr::col_number(),
  gnomADg_EAS_AF = readr::col_number(),
  gnomADg_FIN_AF = readr::col_number(),
  gnomADg_ASJ_AF = readr::col_number(),
  gnomADg_REMAINING_AF = readr::col_number(),
  gnomADg_NFE_AF = readr::col_number(),
  gnomADg_SAS_AF = readr::col_number(),
  EFFECT_PREDICTIONS = readr::col_character(),
  SAMPLE_ID = readr::col_character(),
  VCF_SAMPLE_ID = readr::col_character(),
  GENOME_VERSION = readr::col_character()
)

####----SNV/Indel germline ----####
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
  AMINO_ACID_START = readr::col_integer(),
  AMINO_ACID_END = readr::col_integer(),
  CONSEQUENCE = readr::col_character(),
  IMPACT = readr::col_character(),
  LOSS_OF_FUNCTION = readr::col_logical(),
  LOF_FILTER = readr::col_character(),
  MAXENTSCAN = readr::col_character(),
  MAXENTSCAN2 = readr::col_character(),
  MAXENTSCAN_REF = readr::col_number(),
  MAXENTSCAN_ALT = readr::col_number(),
  MAXENTSCAN_DG = readr::col_number(), #delta gain
  MAXENTSCAN_DL = readr::col_number(), #delta loss
  MAXENTSCAN_DIFF = readr::col_number(),
  MAXENTSCAN_PCT_CHANGE = readr::col_number(),
  NMD = readr::col_character(),
  SPLICE_DONOR_RELEVANT = readr::col_logical(),
  SPLICE_EFFECT_MUTSPLICEDB = readr::col_character(),
  NULL_VARIANT = readr::col_logical(),
  CODING_STATUS = readr::col_character(),
  EXONIC_STATUS = readr::col_character(),
  PROTEIN_CHANGE = readr::col_character(),
  GRANTHAM_DISTANCE = readr::col_character(),
  ALTERATION = readr::col_character(),
  HGVSp_short = readr::col_character(),
  HGVSc = readr::col_character(),
  HGVSp = readr::col_character(),
  CODON = readr::col_character(),
  HGVSc_RefSeq = readr::col_character(),
  CDS_CHANGE = readr::col_character(),
  CDS_RELATIVE_POSITION = readr::col_number(),
  PROTEIN_RELATIVE_POSITION = readr::col_number(),
  EXON = readr::col_character(),
  EXON_AFFECTED = readr::col_integer(),
  EXON_POSITION = readr::col_integer(),
  LAST_EXON = readr::col_logical(),
  LAST_INTRON = readr::col_logical(),
  INTRON_POSITION = readr::col_integer(),
  EXON_INTRON_JUNCTION_SPAN = readr::col_logical(),
  MUTATION_HOTSPOT = readr::col_character(),
  MUTATION_HOTSPOT_CANCERTYPE = readr::col_character(),
  MUTATION_HOTSPOT_MATCH = readr::col_character(),
  BIOMARKER_MATCH = readr::col_character(),
  GWAS_HIT = readr::col_character(),
  PFAM_DOMAIN = readr::col_character(),
  PFAM_DOMAIN_NAME = readr::col_character(),
  PFAM_ENTRY_LOCATIONS = readr::col_character(),
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
  MANE_SELECT = readr::col_character(),
  MANE_SELECT2 = readr::col_character(),
  MANE_PLUS_CLINICAL = readr::col_character(),
  MANE_PLUS_CLINICAL2 = readr::col_character(),
  TSG = readr::col_logical(),
  TSG_RANK = readr::col_integer(),
  TSG_SUPPORT = readr::col_character(),
  ONCOGENE = readr::col_logical(),
  ONCOGENE_RANK = readr::col_integer(),
  ONCOGENE_SUPPORT = readr::col_character(),
  INTOGEN_DRIVER = readr::col_character(),
  INTOGEN_ROLE = readr::col_character(),
  CPG_SOURCE = readr::col_character(),
  GERP_SCORE = readr::col_number(),
  CLINVAR_TRAITS = readr::col_character(),
  CLINVAR_TRAITS_ALL = readr::col_character(),
  CLINVAR_MSID = readr::col_character(),
  CLINVAR_CLNSIG = readr::col_character(),
  CLINVAR_CONTRIB_CLNS_GERMLINE = readr::col_character(),
  CLINVAR_UMLS_CUI = readr::col_character(),
  CLINVAR_CLASSIFICATION = readr::col_character(),
  CLINVAR_CONFLICTED = readr::col_logical(),
  CLINVAR_GOLD_STARS = readr::col_integer(),
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
  DBNSFP_ALPHA_MISSENSE = readr::col_character(),
  DBNSFP_BAYESDEL_ADDAF = readr::col_character(),
  DBNSFP_CLINPRED = readr::col_character(),
  DBNSFP_DEOGEN2 = readr::col_character(),
  DBNSFP_ESM1B = readr::col_character(),
  DBNSFP_FATHMM_XF = readr::col_character(),
  DBNSFP_LIST_S2 = readr::col_character(),
  DBNSFP_META_RNN = readr::col_character(),
  DBNSFP_MUTATIONASSESSOR = readr::col_character(),
  DBNSFP_MUTATIONTASTER = readr::col_character(),
  DBNSFP_MUTFORMER = readr::col_character(),
  DBNSFP_MUTPRED = readr::col_character(),
  DBNSFP_M_CAP = readr::col_character(),
  DBNSFP_PHACTBOOST = readr::col_character(),
  DBNSFP_POLYPHEN2_HVAR = readr::col_character(),
  DBNSFP_PRIMATEAI = readr::col_character(),
  DBNSFP_PROVEAN = readr::col_character(),
  DBNSFP_SIFT = readr::col_character(),
  DBNSFP_SPLICE_SITE_ADA = readr::col_character(),
  DBNSFP_SPLICE_SITE_RF = readr::col_character(),
  gnomADe_AF = readr::col_number(),
  gnomADe_AMR_AF = readr::col_number(),
  gnomADe_AFR_AF = readr::col_number(),
  gnomADe_EAS_AF = readr::col_number(),
  gnomADe_FIN_AF = readr::col_number(),
  gnomADe_ASJ_AF = readr::col_number(),
  gnomADe_REMAINING_AF = readr::col_number(),
  gnomADe_NFE_AF = readr::col_number(),
  gnomADe_SAS_AF = readr::col_number(),
  gnomADg_AF = readr::col_number(),
  gnomADg_AMR_AF = readr::col_number(),
  gnomADg_AFR_AF = readr::col_number(),
  gnomADg_EAS_AF = readr::col_number(),
  gnomADg_FIN_AF = readr::col_number(),
  gnomADg_ASJ_AF = readr::col_number(),
  gnomADg_REMAINING_AF = readr::col_number(),
  gnomADg_NFE_AF = readr::col_number(),
  gnomADg_SAS_AF = readr::col_number(),
  gALL_GRPMAX = readr::col_character(),
  gNC_AMR = readr::col_character(),
  gNC_AFR = readr::col_character(),
  gNC_EAS = readr::col_character(),
  gNC_FIN = readr::col_character(),
  gNC_NFE = readr::col_character(),
  gNC_SAS = readr::col_character(),
  gNC = readr::col_character(),
  gNC_FAF = readr::col_character(),
  DBMTS = readr::col_character(),
  EFFECT_PREDICTIONS = readr::col_character(),
  SAMPLE_ID = readr::col_character(),
  VCF_SAMPLE_ID = readr::col_character(),
  GENOME_VERSION = readr::col_character()
)

####----SNV/Indel - CPSR ----####
data_coltype_defs[['snv_indel_germline_cpsr']] <- readr::cols_only(
  GENOMIC_CHANGE = readr::col_character(),
  GENOTYPE = readr::col_character(),
  ALTERATION = readr::col_character(),
  DP_CONTROL = readr::col_integer(),
  VAR_ID = readr::col_character(),
  VARIANT_CLASS = readr::col_character(),
  CONSEQUENCE = readr::col_character(),
  LOSS_OF_FUNCTION = readr::col_logical(),
  CODING_STATUS = readr::col_character(),
  PROTEIN_CHANGE = readr::col_character(),
  HGVSc = readr::col_character(),
  HGVSp = readr::col_character(),
  HGVSc_RefSeq = readr::col_character(),
  CDS_CHANGE = readr::col_character(),
  EFFECT_PREDICTIONS = readr::col_character(),
  SPLICE_EFFECT = readr::col_character(),
  PFAM_DOMAIN_NAME = readr::col_character(),
  SYMBOL = readr::col_character(),
  ENTREZGENE = readr::col_character(),
  GENENAME = readr::col_character(),
  ENSEMBL_GENE_ID = readr::col_character(),
  ENSEMBL_TRANSCRIPT_ID = readr::col_character(),
  REFSEQ_TRANSCRIPT_ID = readr::col_character(),
  PFAM_DOMAIN = readr::col_character(),
  PFAM_DOMAIN_NAME = readr::col_character(),
  CLINVAR_MSID = readr::col_character(),
  CLINVAR_PHENOTYPE = readr::col_character(),
  CLINVAR_CLASSIFICATION = readr::col_character(),
  CLINVAR_CONFLICTED = readr::col_logical(),
  CLINVAR_GOLD_STARS = readr::col_integer(),
  CLINVAR_VARIANT_ORIGIN = readr::col_character(),
  DBSNP_RSID = readr::col_character(),
  ASSERTION_AUTHORITY = readr::col_character(),
  ASSERTION_RATIONALE = readr::col_character(),
  CPSR_CLASSIFICATION = readr::col_character(),
  CPSR_PATHOGENICITY_SCORE = readr::col_double(),
  ACMG_CODE = readr::col_character(),
  CLASSIFICATION = readr::col_character(),
  SAMPLE_ID = readr::col_character(),
  GENOME_VERSION = readr::col_character()
)


usethis::use_data(data_coltype_defs, overwrite = T)

####----::Output TSVs:: ----####

###----CNA-----####
tsv_cols <- list()

oncokb_annotations <-
  c('MUTATION_EFFECT_OKB',
    'MUTATION_EFFECT_CITATIONS_OKB',
    'MUTATION_EFFECT_DESCRIPTION_OKB',
    'ONCOGENICITY_OKB',
    'TUMOR_TYPE_SUMMARY_OKB',
    'VARIANT_SUMMARY_OKB',
    'HOTSPOT_OKB',
    'VUS_OKB'
  )

usethis::use_data(oncokb_annotations, overwrite = T)

tsv_cols[['cna']] <-
  c('SAMPLE_ID',
    'VAR_ID',
    'GENOME_VERSION',
    'CN_MAJOR',
    'CN_MINOR',
    'LOH',
    "TWOHIT_CANDIDATE_SOMATIC",
    "TWOHIT_CANDIDATE_GERMLINE",
    'FOLD_CHANGE',
    'TUMOR_PLOIDY',
    'TUMOR_PLOIDY_SOURCE',
    'TUMOR_PURITY',
    'SEGMENT_LENGTH_MB',
    'CYTOBAND',
    'EVENT_TYPE',
    'VARIANT_CLASS',
    'VARIANT_CLASS_DISPLAY',
    'SYMBOL',
    'ENTREZGENE',
    'GENENAME',
    'ENSEMBL_GENE_ID',
    'TUMOR_SUPPRESSOR',
    'TUMOR_SUPPRESSOR_SUPPORT',
    'ONCOGENE',
    'ONCOGENE_SUPPORT',
    oncokb_annotations,
    'TRANSCRIPT_OVERLAP',
    'TRANSCRIPT_OVERLAP_PERCENT',
    'ACTIONABILITY_TIER',
    'ACTIONABILITY',
    'BIOMARKER_MATCH',
    'TARGETED_INHIBITORS_ALL2')

####----RNA fusions-----#####
tsv_cols[['fusion']] <-
  c("VAR_ID",
    "VARIANT_CLASS",
    "ENTREZGENE",
    "FUSION_GENE",
    "FUSION_GENE2",
    "SPLIT_READS",
    "SCORE",
    "FUSION_GENE_5P",
    "FUSION_GENE_3P",
    "BREAKPOINT_5P",
    "BREAKPOINT_3P",
    "GENENAME_5P",
    "ONCOGENE_5P",
    "ENSEMBL_TRANSCRIPT_ID_5P",
    "TARGETED_INHIBITORS_5P",
    "TARGETED_INHIBITORS_ALL_5P",
    "GENENAME_3P",
    "ONCOGENE_3P",
    "ENSEMBL_TRANSCRIPT_ID_3P",
    "TARGETED_INHIBITORS_3P",
    "TARGETED_INHIBITORS_ALL_3P",
    "SAMPLE_ALTERATION",
    "MITDB_NUM_EVIDENCE",
    "MITDB_EVIDENCE",
    "ACTIONABILITY_TIER",
    "ACTIONABILITY",
    oncokb_annotations
  )


####----SNV/Indel-----#####
tsv_cols[['snv_indel']] <-
  c('SAMPLE_ID',
    'GENOMIC_CHANGE',
    'GENOME_VERSION',
    'VARIANT_CLASS',
    'SYMBOL',
    'ENTREZGENE',
    'ENSEMBL_GENE_ID',
    'GENENAME',
    'ALTERATION',
    'CDS_CHANGE',
    'HGVSc',
    'HGVSc_RefSeq',
    'HGVSp',
    'HGVSP',
    'SPLICE_EFFECT',
    'MAXENTSCAN',
    'EFFECT_PREDICTIONS',
    'EXON',
    'CONSEQUENCE',
    'PFAM_DOMAIN_NAME',
    'LOSS_OF_FUNCTION',
    'LOF_FILTER',
    'CODING_STATUS',
    'EXONIC_STATUS',
    'DP_TUMOR',
    'VAF_TUMOR',
    'AD_TUMOR',
    'CALL_CONFIDENCE',
    'DP_CONTROL',
    'VAF_CONTROL',
    'AD_CONTROL',
    'MUTATION_HOTSPOT',
    'MUTATION_HOTSPOT_CANCERTYPE',
    'ACTIONABILITY_TIER',
    'ACTIONABILITY',
    'ONCOGENICITY',
    'ONCOGENICITY_CODE',
    'ONCOGENICITY_SCORE',
    oncokb_annotations,
    'CANONICAL',
    'CCDS',
    'UNIPROT_ACC',
    'ENSEMBL_TRANSCRIPT_ID',
    'ENSEMBL_PROTEIN_ID',
    'REFSEQ_TRANSCRIPT_ID',
    'REFSEQ_PROTEIN_ID',
    'MANE_SELECT',
    'MANE_PLUS_CLINICAL',
    'ONCOGENE',
    'ONCOGENE_SUPPORT',
    'TUMOR_SUPPRESSOR',
    'TUMOR_SUPPRESSOR_SUPPORT',
    'TARGETED_INHIBITORS2',
    'REGULATORY_ANNOTATION',
    'VEP_ALL_CSQ',
    'gnomADe_AF',
    'gnomADg_AF',
    'DBSNP_RSID',
    'COSMIC_ID',
    'TCGA_FREQUENCY',
    'TCGA_PANCANCER_COUNT',
    'CLINVAR_MSID',
    'CLINVAR_CLASSIFICATION',
    'CLINVAR_VARIANT_ORIGIN',
    'CLINVAR_NUM_SUBMITTERS',
    'CLINVAR_GOLD_STARS',
    'CLINVAR_CONFLICTED'
    )

####----SNV/Indel unfiltered-----#####
tsv_cols[['snv_indel_unfiltered']] <-
  c('SAMPLE_ID',
    'GENOMIC_CHANGE',
    'GENOME_VERSION',
    'VARIANT_CLASS',
    'SOMATIC_CLASSIFICATION',
    'SYMBOL',
    'ENTREZGENE',
    'ENSEMBL_GENE_ID',
    'GENENAME',
    'ALTERATION',
    'CDS_CHANGE',
    'HGVSc',
    'HGVSc_RefSeq',
    'HGVSp',
    'HGVSP',
    'SPLICE_EFFECT',
    'MAXENTSCAN',
    'EFFECT_PREDICTIONS',
    'EXON',
    'CONSEQUENCE',
    'PFAM_DOMAIN_NAME',
    'LOSS_OF_FUNCTION',
    'LOF_FILTER',
    'CODING_STATUS',
    'EXONIC_STATUS',
    'DP_TUMOR',
    'VAF_TUMOR',
    'AD_TUMOR',
    'CALL_CONFIDENCE',
    'MUTATION_HOTSPOT',
    'MUTATION_HOTSPOT_CANCERTYPE',
    'ACTIONABILITY_TIER',
    'ACTIONABILITY',
    'ONCOGENICITY',
    'ONCOGENICITY_CODE',
    'ONCOGENICITY_SCORE',
    oncokb_annotations,
    'CANONICAL',
    'CCDS',
    'UNIPROT_ACC',
    'ENSEMBL_TRANSCRIPT_ID',
    'ENSEMBL_PROTEIN_ID',
    'REFSEQ_TRANSCRIPT_ID',
    'REFSEQ_PROTEIN_ID',
    'MANE_SELECT',
    'MANE_PLUS_CLINICAL',
    'ONCOGENE',
    'ONCOGENE_SUPPORT',
    'TUMOR_SUPPRESSOR',
    'TUMOR_SUPPRESSOR_SUPPORT',
    'TARGETED_INHIBITORS2',
    'REGULATORY_ANNOTATION',
    'VEP_ALL_CSQ',
    'gnomADg_AF',
    'gnomADg_AFR_AF',
    'gnomADg_AMR_AF',
    'gnomADg_ASJ_AF',
    'gnomADg_EAS_AF',
    'gnomADg_FIN_AF',
    'gnomADg_NFE_AF',
    'gnomADg_OTH_AF',
    'gnomADg_SAS_AF',
    'gnomADe_AF',
    'gnomADe_AFR_AF',
    'gnomADe_AMR_AF',
    'gnomADe_EAS_AF',
    'gnomADe_FIN_AF',
    'gnomADe_NFE_AF',
    'gnomADe_SAS_AF',
    'DBSNP_RSID',
    'COSMIC_ID',
    'TCGA_FREQUENCY',
    'TCGA_PANCANCER_COUNT',
    'CLINVAR_MSID',
    'CLINVAR_CLASSIFICATION',
    'CLINVAR_VARIANT_ORIGIN',
    'CLINVAR_NUM_SUBMITTERS',
    'CLINVAR_GOLD_STARS',
    'CLINVAR_CONFLICTED'
  )
usethis::use_data(tsv_cols, overwrite = T)
####----::Data tables display columns:: ----####

table_display_cols <- list()

####----CNA -----#####
table_display_cols[['cna_other_oncogenic']] <-
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
    "ACTIONABLE_GENE",
    "TRANSCRIPT_OVERLAP",
    'TRANSCRIPT_OVERLAP_PERCENT',
    "GLOBAL_ASSOC_RANK",
    "TISSUE_ASSOC_RANK",
    "SEGMENT",
    "SEGMENT_LENGTH_MB",
    "LOH",
    "CN_MAJOR",
    "CN_MINOR",
    "TARGETED_INHIBITORS_ALL",
    "GENOME_VERSION")

####----SNV/Indel - germline filtered-----#####
table_display_cols[['snv_indel_germline_filtered']] <-
  c('GENOMIC_CHANGE',
    'VARIANT_CLASS',
    'EXCLUSION_CRITERIA',
    'SYMBOL',
    'CONSEQUENCE',
    'HGVSc',
    'HGVSc_RefSeq',
    'DBSNP_RSID',
    'COSMIC_ID',
    'TCGA_FREQUENCY',
    'CLINVAR_CLASSIFICATION',
    'VAF_TUMOR',
    'DP_TUMOR',
    'gnomADg_AF',
    'gnomADe_AF',
    'gnomADe_AFR_AF',
    'gnomADe_AMR_AF',
    'gnomADe_EAS_AF',
    'gnomADe_FIN_AF',
    'gnomADe_NFE_AF',
    'gnomADe_SAS_AF')

####----SNV/Indel-----#####
table_display_cols[['snv_indel']] <-
  c('SYMBOL',
    'ALTERATION',
    'GENENAME',
    'CONSEQUENCE',
    'ONCOGENICITY',
    'PROTEIN_DOMAIN',
    'MUTATION_HOTSPOT',
    'COSMIC_ID',
    'PROTEIN_CHANGE',
    'CDS_CHANGE',
    'HGVSc',
    'HGVSc_RefSeq',
    'MUTATION_HOTSPOT_CANCERTYPE',
    'LOSS_OF_FUNCTION',
    'LOF_FILTER',
    'TCGA_FREQUENCY',
    'PREDICTED_EFFECT',
    'SPLICE_EFFECT',
    'ONCOGENICITY_CODE',
    'ONCOGENICITY_SCORE',
    'ONCOGENICITY_DOC',
    "MUTATION_EFFECT_OKB",
    "MUTATION_EFFECT_DESCRIPTION_OKB",
    "ONCOGENICITY_OKB",
    'VARIANT_SUMMARY_OKB',
    'TUMOR_TYPE_SUMMARY_OKB',
    'VUS_OKB',
    'DBSNP_RSID',
    'CLINVAR',
    'CLINVAR_CLASSIFICATION',
    'ENSEMBL_TRANSCRIPT_ID',
    'REFSEQ_TRANSCRIPT_ID',
    'MANE_SELECT',
    'TARGETED_INHIBITORS',
    'TARGETED_INHIBITORS_ALL',
    'ONCOGENE',
    'TUMOR_SUPPRESSOR',
    'CANCERGENE_EVIDENCE',
    'GLOBAL_ASSOC_RANK',
    'TISSUE_ASSOC_RANK',
    'CALL_CONFIDENCE',
    'DP_TUMOR',
    'VAF_TUMOR',
    'DP_CONTROL',
    'VAF_CONTROL',
    'gnomADg_AF',
    'gnomADe_AF',
    'GENOMIC_CHANGE',
    'GENOME_VERSION')

usethis::use_data(table_display_cols, overwrite = T)

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
      url_prefix = "https://platform.opentargets.org/target/",
      link_key_var = "ENSEMBL_GENE_ID",
      rename = TRUE,
      link_display_var = "GENENAME_RAW"),
    # data.frame(
    #   name = "GENENAME",
    #   group_by_var = "VAR_ID",
    #   url_prefix = "https://www.ncbi.nlm.nih.gov/gene/",
    #   link_key_var = "ENTREZGENE",
    #   rename = TRUE,
    #   link_display_var = "GENENAME_RAW"),
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
      name = "MANE_SELECT",
      group_by_var = "VAR_ID",
      url_prefix = "https://www.ncbi.nlm.nih.gov/nuccore/",
      link_key_var = "MANE_SELECT_RAW",
      rename = TRUE,
      link_display_var = "MANE_SELECT_RAW",
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
effect_prediction_algos <- as.data.frame(
  readr::read_tsv(file = "data-raw/effect_prediction_algorithms.tsv",
             show_col_types = F))
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
         "|dihydropyrimidine dehydrogenase deficiency",
         "|lymphoproliferative|stat3|teratoma|thrombocytop|tp63|weaver",
         "|pheochromo|gorlin|telangiectasia|hemangiom|osteochondro|response",
         "|polg-related|ras-associated|dyskeratosis",
         "|tyrosine kinase inhibitor|fluorouracil|methotrexate|capecitabine",
         "|imatinib|erlotinib|gefitinib|sunitinib|vemurafenib|crizotinib|cetuximab",
         "|dabrafenib|osimertinib|afatinib|neratinib|lapatinib|ponatinib|regorafenib",
         "|vandetanib|bosutinib|axitinib|cabozantinib|brigatinib|ripretinib|lenvatinib",
         "|midostaurin|duvelisib|copanlisib|alpelisib|larotrectinib|entrectinib",
         "|programmed death ligand|adp-ribose|doxorubicin|cisplatin|carboplatin",
         "|familial adenomatous polyposis|hereditary breast and ovarian cancer",
         "|hereditary nonpolyposis colorectal cancer|multiple endocrine neoplasia",
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

exonic_filter_levels <-
  c("SOMATIC - EXONIC",
    "SOMATIC - NONEXONIC")

usethis::use_data(germline_filter_levels, overwrite = T)
usethis::use_data(exonic_filter_levels, overwrite = T)

rm(cancer_phenotypes_regex,
   effect_prediction_algos,
   tcga_cohorts,
   immune_celltypes2,
   artefact_signatures,
   context_order,
   cosmic_signatures,
   cosmic_sbs_signatures_no_artefacts,
   cosmic_sbs_signatures_all,
   cosmic_sbs_signatures,
   germline_filter_levels,
   sig,
   variant_db_url,
   color_palette,
   data_coltype_defs,
   bm_evidence,
   table_display_cols,
   tsv_cols,
   c)

oncogenicity_criteria <- readr::read_tsv(
  "data-raw/oncogenicity.tsv", show_col_types = F
)
usethis::use_data(oncogenicity_criteria, overwrite = T)

oncokb_base_api_url <-
  "https://www.oncokb.org/api/v1/annotate/"
usethis::use_data(oncokb_base_api_url, overwrite = T)
