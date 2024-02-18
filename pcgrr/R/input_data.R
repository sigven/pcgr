#' Function that reads and validates a fully annotated CNA file from PCGR
#' pre-report pipeline
#'
#' @param fname Path to file name
#' @param ref_data Object with reference data
#' @param settings Object with PCGR report configuration
#'
#' @export
load_somatic_cna <- function(
    fname,
    ref_data = NULL,
    settings = NULL) {

  log4r_info(paste0(
    "Reading annotated molecular dataset (DNA) - somatic copy number aberrations"))

  callset_cna <- pcgrr::load_dna_variants(
    fname = fname,
    cols = pcgrr::data_coltype_defs$cna_somatic_raw,
    ref_data = ref_data,
    vartype = 'cna',
    primary_site = settings[['conf']][['sample_properties']]$site,
    retained_info_tags = "None",
    variant_origin = "Somatic")

  tumor_site <-
    settings[['conf']][['sample_properties']][['site']]

  if (NROW(callset_cna$variant) > 0) {
    callset_cna[['variant']] <- callset_cna[['variant']] |>
      pcgrr::append_cancer_gene_evidence(
        ref_data = ref_data,
        site = tumor_site,
        pos_var = 'SEGMENT_START') |>
      pcgrr::append_drug_var_link(
        ref_data = ref_data
      ) |>
      dplyr::arrange(
        .data$TIER,
        dplyr::desc(.data$TISSUE_ASSOC_RANK),
        dplyr::desc(.data$GLOBAL_ASSOC_RANK))
  }

  return(callset_cna)

}

#' Function that reads and validates an annotated somatic SNV/InDel
#' file from PCGR pre-reporting pipeline
#'
#' @param fname Path to file name
#' @param ref_data Object with reference data
#' @param settings Object with PCGR report configuration
#'
#' @export
load_somatic_snv_indel <- function(
    fname = NA,
    ref_data = NULL,
    settings = NULL) {

  log4r_info(paste0(
    "Reading annotated molecular dataset (DNA) - somatic SNV/InDels"))

  callset <- pcgrr::load_dna_variants(
    fname = fname,
    cols = pcgrr::data_coltype_defs$snv_indel_somatic_raw,
    ref_data = ref_data,
    vartype = 'snv_indel',
    primary_site = settings[['conf']][['sample_properties']]$site,
    retained_info_tags =
      settings[['conf']][['other']]$retained_vcf_info_tags,
    variant_origin = "Somatic")

  tumor_site <-
    settings[['conf']][['sample_properties']][['site']]


  callset[['variant_unfiltered']] <- data.frame()
  callset[['variant']] <- callset[['variant']] |>
    pcgrr::append_dbnsfp_var_link() |>
    pcgrr::append_dbmts_var_link() |>
    pcgrr::append_tcga_var_link() |>
    pcgrr::append_annotation_links() |>
    pcgrr::append_drug_var_link(ref_data = ref_data) |>
    pcgrr::append_tfbs_annotation() |>
    pcgrr::append_cancer_gene_evidence(ref_data = ref_data,
                                       site = tumor_site)

  if (settings$conf$assay_properties$vcf_tumor_only == 1) {
    callset[['variant_unfiltered']] <- callset[['variant']]
    callset[['variant']] <- callset[['variant']] |>
      ## assign evidence tags for germline/somatic state of variants,
      ## partially based on user-defined criteria
      ## (population allele frequency thresholds)
        pcgrr::assign_somatic_germline_evidence2(
          settings = settings) |>

      ## assign somatic classification based on accumulation
      ## of evidence tags and user-defined options
        pcgrr::assign_somatic_classification(
          settings = settings)
  }

  callset[['variant']] <- callset[['variant']] |>
    dplyr::arrange(.data$TIER,
                   dplyr::desc(.data$ONCOGENICITY_SCORE),
                   dplyr::desc(.data$TISSUE_ASSOC_RANK),
                   dplyr::desc(.data$GLOBAL_ASSOC_RANK))

  # callset <-
  #   pcgrr::expand_biomarker_items(
  #     callset = callset,
  #     variant_origin = "somatic")

  return(callset)


}

#' Function that reads and validates CNA or SNV/InDel TSV files
#' file from PCGR/CPSR pre-report (Python) pipeline
#'
#' @param fname Path to raw input file with DNA aberrations (PCGR/CPSR)
#' @param cols column type definitions of raw input file
#' @param ref_data reference data object
#' @param vartype type of DNA aberrations ('snv_indel','cna')
#' @param primary_site primary site of tumor
#' @param retained_info_tags VCF INFO tags to be retained in output (SNVs/InDels)
#' @param variant_origin Germline/Somatic
#'
#' @export
#'
load_dna_variants <- function(
    fname = NA,
    cols = NULL,
    ref_data = NULL,
    vartype = 'snv_indel',
    primary_site = "Any",
    retained_info_tags = "None",
    variant_origin = "Somatic") {

  pcgrr::check_file_exists(fname)
  calls_raw <- suppressWarnings(
    as.data.frame(
      readr::read_tsv(
        file = fname,
        na = c(".","NA"),
        show_col_types = F,
        progress = F
      )
    )
  )

  ## check that all columns are present among columns
  ## read from file
  compulsary_cols <-
    names(cols$cols)

  raw_col_check <-
    tibble::has_name(calls_raw, compulsary_cols)
  if (FALSE %in% raw_col_check) {
    missing_cols <-
      compulsary_cols[!raw_col_check]
    log4r_fatal(
      paste0("Missing required columns in input file ",
             fname, " - ",
             paste(missing_cols, collapse=", ")))
  }

  cols_including_retained <- cols
  retained_cols <- NULL
  if (retained_info_tags != "None") {
    retained_cols <- stringr::str_split(
      retained_info_tags, pattern = ",")[[1]]
    for(c in retained_cols) {
      if (c %in% colnames(calls_raw)) {
        col_retain <- readr::cols_only(
          !!rlang::sym(c) := readr::col_character()
        )
        cols_including_retained$cols <-
          c(cols_including_retained$cols, col_retain$cols)
      }
    }
  }

  calls <- suppressWarnings(
    as.data.frame(
      readr::read_tsv(
        file = fname,
        col_types =
          cols_including_retained,
        na = c(".","NA"),
        progress = F
      )
    )
  )

  retained_cols_renamed <- c()
  if (!is.null(retained_cols)) {
    for(c in retained_cols) {
      if (c %in% colnames(calls)) {
        new_col <- paste0('VCF_INFO_', c)
        retained_cols_renamed <- c(
          retained_cols_renamed, new_col
        )
        calls[,new_col] <- calls[,c]
        calls[,c] <- NULL
      }
    }
  }

  results <- list()
  results[['variant']] <- calls
  results[['biomarker_evidence']] <- list()
  results[['biomarker_evidence']][['items']] <-
    data.frame()

  ## Rename annotations for more clarity
  if ("TSG" %in% colnames(results[['variant']])) {
    results[['variant']] <-
      results[['variant']] |>
      dplyr::rename(
        TUMOR_SUPPRESSOR = "TSG"
      )
  }

  if ("ONCOGENICITY_CLASSIFICATION" %in% colnames(results[['variant']])) {
    results[['variant']] <-
      results[['variant']] |>
      dplyr::rename(
        ONCOGENICITY = "ONCOGENICITY_CLASSIFICATION"
      )
  }

  if ("VEP_ALL_CSQ" %in% colnames(results[['variant']])) {
    results[['variant']] <-
      results[['variant']] |>
      dplyr::mutate(
        VEP_ALL_CSQ = stringr::str_replace_all(
          .data$VEP_ALL_CSQ, ",",", "
        )
      )
  }

  if ("HGVSp_short" %in% colnames(results[['variant']])) {
    results[['variant']] <-
      results[['variant']] |>
      dplyr::rename(
        HGVSP = "HGVSp_short"
      )
  }

  if ("TSG_RANK" %in% colnames(results[['variant']])) {
    results[['variant']] <-
      results[['variant']] |>
      dplyr::rename(
        TUMOR_SUPPRESSOR_RANK = "TSG_RANK"
      )
  }

  if (vartype == 'cna') {

    results[['variant']] <-
      results[['variant']] |>
      dplyr::mutate(REFSEQ_TRANSCRIPT_ID = dplyr::if_else(
        is.na(.data$REFSEQ_TRANSCRIPT_ID),
        "",
        as.character(.data$REFSEQ_TRANSCRIPT_ID)
      )) |>
      dplyr::mutate(
        TRANSCRIPT_OVERLAP = paste(
          .data$ENSEMBL_TRANSCRIPT_ID,
          .data$REFSEQ_TRANSCRIPT_ID,
          .data$TRANSCRIPT_START,
          .data$TRANSCRIPT_END,
          .data$TRANSCRIPT_OVERLAP_PERCENT, sep="|"
        )) |>
      dplyr::select(
        -c("ENSEMBL_TRANSCRIPT_ID",
           "REFSEQ_TRANSCRIPT_ID",
           "TRANSCRIPT_START",
           "TRANSCRIPT_END")) |>
      dplyr::group_by(
        dplyr::across(-c("TRANSCRIPT_OVERLAP",
                         "TRANSCRIPT_OVERLAP_PERCENT"))) |>
      dplyr::summarise(
        TRANSCRIPT_OVERLAP = paste(.data$TRANSCRIPT_OVERLAP, collapse=", "),
        MAX_TRANSCRIPT_OVERLAP_PERCENT =
          max(.data$TRANSCRIPT_OVERLAP_PERCENT, na.rm = T),
        .groups = "drop"
      )

  }

  if ("BIOMARKER_MATCH" %in% colnames(results[['variant']]) &
     "VAR_ID" %in% colnames(results[['variant']])) {

    biomarker_set <-
      results[['variant']] |>
      dplyr::filter(!is.na(.data$BIOMARKER_MATCH))

    if (NROW(biomarker_set) > 0) {

      citations <- as.data.frame(
        ref_data[['biomarker']][['literature']] |>
          dplyr::select(
            c("EVIDENCE_ID",
              "LINK")
          ) |>
          tidyr::separate_rows(
            c("EVIDENCE_ID"),
            sep=";"
          ) |>
          dplyr::group_by(
            EVIDENCE_ID
          ) |>
          dplyr::summarise(
            CITATION = paste(
              unique(.data$LINK), collapse = ", "
            )
          )
      )
      results[['biomarker_evidence']][['items']] <-
        as.data.frame(
          biomarker_set |>
            dplyr::select(
              c("VAR_ID",
                "VARIANT_CLASS",
                "BIOMARKER_MATCH"),
              ) |>
            dplyr::distinct() |>
            tidyr::separate_rows(
              "BIOMARKER_MATCH", sep=",") |>
            tidyr::separate(
              "BIOMARKER_MATCH",
              into = c("BIOMARKER_SOURCE",
                       "VARIANT_ID",
                       "EVIDENCE_ITEMS",
                       "BIOMARKER_MATCHTYPE"),
              sep = "\\|"
            ) |>
            dplyr::mutate(VARIANT_ID = as.character(.data$VARIANT_ID)) |>
            dplyr::left_join(
              dplyr::select(
                ref_data[['biomarker']][['variant']],
                c("VARIANT_ID", "ENTREZGENE","BIOMARKER_SOURCE")),
              by = c("VARIANT_ID","BIOMARKER_SOURCE")) |>
            dplyr::rename(BIOMARKER_MATCH = .data$BIOMARKER_MATCHTYPE) |>
            dplyr::mutate(BIOMARKER_RESOLUTION = dplyr::case_when(
              stringr::str_detect(.data$BIOMARKER_MATCH,"by_cna_segment") ~ "gene",
              stringr::str_detect(.data$BIOMARKER_MATCH,"by_genomic_coord") ~ "genomic",
              !stringr::str_detect(.data$BIOMARKER_MATCH,"by_genomic_coord") &
                stringr::str_detect(.data$BIOMARKER_MATCH,"by_hgvsp_principal") ~ "hgvsp",
              !stringr::str_detect(.data$BIOMARKER_MATCH,"by_genomic_coord") &
                !stringr::str_detect(.data$BIOMARKER_MATCH,"by_hgvsp_principal") &
                stringr::str_detect(.data$BIOMARKER_MATCH,"by_codon_principal") ~ "codon",
              !stringr::str_detect(.data$BIOMARKER_MATCH,"by_genomic_coord") &
                !stringr::str_detect(.data$BIOMARKER_MATCH,"by_hgvsp_principal") &
                !stringr::str_detect(.data$BIOMARKER_MATCH,"by_codon_principal") &
                stringr::str_detect(.data$BIOMARKER_MATCH,"by_hgvsp_nonprincipal")~ "hgvsp_nonprincipal",
              !stringr::str_detect(.data$BIOMARKER_MATCH,"by_genomic_coord") &
                !stringr::str_detect(.data$BIOMARKER_MATCH,"by_hgvsp_principal") &
                !stringr::str_detect(.data$BIOMARKER_MATCH,"by_codon_principal") &
                stringr::str_detect(.data$BIOMARKER_MATCH,"by_exon_") ~ "exon",
              !stringr::str_detect(.data$BIOMARKER_MATCH,"by_genomic_coord") &
                !stringr::str_detect(.data$BIOMARKER_MATCH,"by_hgvsp_") &
                !stringr::str_detect(.data$BIOMARKER_MATCH,"by_codon_") &
                !stringr::str_detect(.data$BIOMARKER_MATCH,"by_exon_") &
                stringr::str_detect(.data$BIOMARKER_MATCH,"by_gene_") ~ "gene",
              TRUE ~ as.character('other')
            )) |>
            tidyr::separate_rows(
              "EVIDENCE_ITEMS", sep="&") |>
            tidyr::separate(
              "EVIDENCE_ITEMS",
              into = c("EVIDENCE_ID",
                       "PRIMARY_SITE",
                       "CLINICAL_SIGNIFICANCE",
                       "EVIDENCE_LEVEL",
                       "EVIDENCE_TYPE",
                       "VARIANT_ORIGIN"),
              sep = ":"
            ) |>
            dplyr::select(
              -c("CLINICAL_SIGNIFICANCE",
                 "EVIDENCE_LEVEL",
                 "EVIDENCE_TYPE")) |>
            dplyr::mutate(PRIMARY_SITE = stringr::str_replace_all(
              .data$PRIMARY_SITE, "_"," "
            )) |>
            dplyr::left_join(
              citations, by = "EVIDENCE_ID",
              relationship = "many-to-many"
            ) |>
            dplyr::left_join(
              dplyr::select(
                ref_data[['biomarker']][['clinical']],
                c("EVIDENCE_ID",
                "CANCER_TYPE",
                "EVIDENCE_TYPE",
                "THERAPEUTIC_CONTEXT",
                "CLINICAL_SIGNIFICANCE",
                "EVIDENCE_LEVEL",
                "EVIDENCE_DESCRIPTION",
                "MOLECULAR_PROFILE_NAME",
                "MOLECULAR_PROFILE_TYPE",
                "RATING",
                "EVIDENCE_DIRECTION")
              ), by = c("EVIDENCE_ID"),
              relationship = "many-to-many"
            ) |>
            dplyr::rename(
              BM_VARIANT_ID = .data$VARIANT_ID,
              BM_EVIDENCE_ID = .data$EVIDENCE_ID,
              BM_SOURCE = .data$BIOMARKER_SOURCE,
              BM_RESOLUTION = .data$BIOMARKER_RESOLUTION,
              BM_MATCH = .data$BIOMARKER_MATCH,
              BM_PRIMARY_SITE = .data$PRIMARY_SITE,
              BM_EVIDENCE_TYPE = .data$EVIDENCE_TYPE,
              BM_CANCER_TYPE = .data$CANCER_TYPE,
              BM_VARIANT_ORIGIN = .data$VARIANT_ORIGIN,
              BM_EVIDENCE_LEVEL = .data$EVIDENCE_LEVEL,
              BM_EVIDENCE_DESCRIPTION = .data$EVIDENCE_DESCRIPTION,
              BM_THERAPEUTIC_CONTEXT = .data$THERAPEUTIC_CONTEXT,
              BM_CLINICAL_SIGNIFICANCE = .data$CLINICAL_SIGNIFICANCE,
              BM_CITATION = .data$CITATION,
              BM_RATING = .data$RATING,
              BM_EVIDENCE_DIRECTION = .data$EVIDENCE_DIRECTION,
              BM_MOLECULAR_PROFILE_NAME = .data$MOLECULAR_PROFILE_NAME,
              BM_MOLECULAR_PROFILE_TYPE = .data$MOLECULAR_PROFILE_TYPE
            ) |>
            dplyr::select(
              c("VAR_ID",
                "VARIANT_CLASS",
                "ENTREZGENE",
                "BM_SOURCE",
                "BM_VARIANT_ID",
                "BM_EVIDENCE_ID",
                "BM_EVIDENCE_TYPE",
                "BM_EVIDENCE_LEVEL",
                "BM_EVIDENCE_DESCRIPTION",
                "BM_EVIDENCE_DIRECTION",
                "BM_CLINICAL_SIGNIFICANCE",
                "BM_VARIANT_ORIGIN",
                "BM_CANCER_TYPE",
                "BM_PRIMARY_SITE",
                "BM_MATCH",
                "BM_RESOLUTION"),
              dplyr::everything()
            ) |>
            dplyr::filter(
                .data$BM_VARIANT_ORIGIN == variant_origin &
                  .data$BM_MOLECULAR_PROFILE_TYPE == "Any") |>
            dplyr::distinct()
        )
    }

    if (variant_origin == "Somatic") {
      results <- pcgrr::assign_acmg_tiers(
        vartype = vartype,
        variants_df = results$variant,
        primary_site = primary_site,
        biomarker_items =
          results$biomarker_evidence$items
      )
    }

  }else{
    log4r_fatal("Input data does not contain 'BIOMARKER_MATCH' column - fatal")
  }

  results[['retained_info_tags']] <- paste(
    retained_cols_renamed, collapse=","
  )

  return(results)

}