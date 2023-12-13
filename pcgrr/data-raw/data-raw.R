#---- color_palette ----#
color_palette <- list()
for (c in c("pathogenicity", "clinical_evidence", "tier",
           "report_color", "warning", "success", "none")) {
  color_palette[[c]] <- list()
  color_palette[[c]][["levels"]] <- c()
  color_palette[[c]][["values"]] <- c()

  if (c == "pathogenicity") {
    color_palette[[c]][["levels"]] <-
      c("Pathogenic", "Likely_Pathogenic",
        "VUS", "Likely_Benign", "Benign")
    color_palette[[c]][["values"]] <-
      c("#9E0142", "#D53E4F", "#000000", "#78C679", "#077009")
  }
  if (c == "clinical_evidence") {
    color_palette[[c]][["levels"]] <-
      c("A: Validated", "A: FDA/NCCN/ELN guidelines",
        "B: Clinical evidence", "B1: Clinical evidence: late trials",
        "B2: Clinical evidence: early trials", "C: Case study",
        "D: Preclinical evidence", "E: Indirect evidence")
    color_palette[[c]][["values"]] <-
      c("#009E73", "#009E73", "#56B4E9", "#56B4E9",
        "#56B4E9", "#0072B2", "#E69F00", "#F0E442")
  }
  if (c == "tier") {
    color_palette[[c]][["levels"]] <-
      c("TIER1", "TIER2", "TIER3", "TIER4", "NONCODING",
        "COL6", "COL7", "COL8", "COL9", "COL10", "COL11",
        "COL12", "COL13", "COL14", "COL15", "COL16",
        "COL17", "COL18", "COL19", "COL20", "COL21",
        "COL22", "COL23", "COL24", "COL25")
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
  if (c == "warning") {
    color_palette[[c]][["levels"]] <- c("warning")
    color_palette[[c]][["values"]] <- c("#ff7518")
  }
  if (c == "none") {
    color_palette[[c]][["levels"]] <- c("none")
    color_palette[[c]][["values"]] <- c("#868686")
  }
  if (c == "success") {
    color_palette[[c]][["levels"]] <- c("success")
    color_palette[[c]][["values"]] <- c("#00a65a")
  }

}

usethis::use_data(color_palette, overwrite = T)

#-----evidence types---------#
evidence_types <- c("predictive","prognostic","diagnostic",
                    "oncogenic","predisposing","functional")
usethis::use_data(evidence_types)

#-----evidence levels---------#
evidence_levels <- c("any","A_B","C_D_E")
usethis::use_data(evidence_levels)


#-----input column names/types-----#
data_coltype_defs <- list()
data_coltype_defs[['cna_somatic_raw']] <- readr::cols_only(
  CHROM = readr::col_character(),
  SEGMENT_START = readr::col_double(),
  SEGMENT_END = readr::col_double(),
  VAR_ID = readr::col_character(),
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
  SEGMENT_LINK = readr::col_character(),
  SAMPLE_ID = readr::col_character())

data_coltype_defs[['snv_indel_somatic_raw']] <- readr::cols_only(
  CHROM = readr::col_character(),
  POS = readr::col_double(),
  REF = readr::col_character(),
  ALT = readr::col_character(),
  DP_TUMOR = readr::col_integer(),
  AF_TUMOR = readr::col_number(),
  DP_CONTROL = readr::col_integer(),
  AF_CONTROL = readr::col_number(),
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
  ONCOGENICITY_CLASSIFICATION = readr::col_character(),
  ONCOGENICITY_CLASSIFICATION_CODE = readr::col_character(),
  ONCOGENICITY_SCORE = readr::col_number(),
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
  PANEL_OF_NORMALS = readr::col_logical(),
  DBSNPRSID = readr::col_character(),
  COSMIC_MUTATION_ID = readr::col_character(),
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
  gnomAD_AF = readr::col_number(),
  gnomAD_AMR_AF = readr::col_number(),
  gnomAD_AFR_AF = readr::col_number(),
  gnomAD_EAS_AF = readr::col_number(),
  gnomAD_FIN_AF = readr::col_number(),
  gnomAD_ASJ_AF = readr::col_number(),
  gnomAD_OTH_AF = readr::col_number(),
  gnomAD_NFE_AF = readr::col_number(),
  gnomAD_SAS_AF = readr::col_number(),
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
  DBSNPRSID = readr::col_character(),
  COSMIC_MUTATION_ID = readr::col_character(),
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
  gnomAD_AF = readr::col_number(),
  gnomAD_AMR_AF = readr::col_number(),
  gnomAD_AFR_AF = readr::col_number(),
  gnomAD_EAS_AF = readr::col_number(),
  gnomAD_FIN_AF = readr::col_number(),
  gnomAD_ASJ_AF = readr::col_number(),
  gnomAD_OTH_AF = readr::col_number(),
  gnomAD_NFE_AF = readr::col_number(),
  gnomAD_SAS_AF = readr::col_number(),
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

#---- variant_db_url ----#
variant_db_url <-
  data.frame(
    name = "DBSNP",
    group_by_var = "VAR_ID",
    url_prefix = "https://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=",
    link_key_var = "DBSNPRSID",
    link_display_var = "DBSNPRSID")
variant_db_url <-
  dplyr::bind_rows(
    variant_db_url,
    data.frame(name = "CLINVAR",
               group_by_var = "VAR_ID",
               url_prefix = "http://www.ncbi.nlm.nih.gov/clinvar/variation/",
               link_key_var = "CLINVAR_MSID",
               link_display_var = "CLINVAR_TRAITS_ALL"),
    data.frame(name = "OFFICIAL_GENENAME",
               group_by_var = "VAR_ID",
               url_prefix = "https://www.ncbi.nlm.nih.gov/gene/",
               link_key_var = "ENTREZGENE",
               link_display_var = "GENENAME"),
    data.frame(name = "PROTEIN_DOMAIN",
               group_by_var = "VAR_ID",
               url_prefix = "https://www.ebi.ac.uk/interpro/entry/pfam/",
               link_key_var = "PFAM_DOMAIN",
               link_display_var = "PFAM_DOMAIN_NAME",
               stringsAsFactors = F),
    data.frame(name = "COSMIC",
               group_by_var = "VAR_ID",
               url_prefix = "https://cancer.sanger.ac.uk/cosmic/search?q=",
               link_key_var = "COSMIC_MUTATION_ID",
               link_display_var = "COSMIC_MUTATION_ID",
               stringsAsFactors = F),
    data.frame(name = "REFSEQ",
               group_by_var = "VAR_ID",
               url_prefix = "https://www.ncbi.nlm.nih.gov/nuccore/",
               link_key_var = "REFSEQ_TRANSCRIPT_ID",
               link_display_var = "REFSEQ_TRANSCRIPT_ID",
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
         "|exostoses|lynch|drash|wilm|perlman|fibrofolliculomas|hippel|hamartom",
         "|bloom|werner|peutz|legius|tuberous|exostosis|angiomyolipoma",
         "|lymphoproliferative|stat3|teratoma|thrombocytop|tp63|wiskott|weaver",
         "|pheochromo|gorlin|telangiectasia|hemangiom|osteochondro|response",
         "|polg-related|ras-associated|dyskeratosis",
         "|waardenburg|beckwidth|birt-hogg|costello|diamond|cardio-facio|frasier",
         "|hirschsprung|hydrocephalus|hyperparathyroidism|immunodeficiency",
         "|infantile myofibromatosis|leopard|proteus|rothmund|russel)")
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
   c)
