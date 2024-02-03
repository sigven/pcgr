#' Function that reads and validates a fully annotated CNA file from PCGR
#' pre-report pipeline
#'
#' @param fname Path to file name
#' @param ref_data Object with reference data
#'
#' @export
load_somatic_cna <- function(fname, ref_data = NULL){

  log4r_info(paste0(
    "Reading annotated molecular dataset (DNA) - somatic copy number aberrations"))

  callset <- load_dna_variants(
    fname = fname,
    cols = pcgrr::data_coltype_defs$cna_somatic_raw,
    ref_data = ref_data,
    retained_info_tags = "None",
    variant_origin = "Somatic")

  return(callset)

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
    settings = NULL){

  log4r_info(paste0(
    "Reading annotated molecular dataset (DNA) - somatic SNV/InDels"))

  callset <- load_dna_variants(
    fname = fname,
    cols = pcgrr::data_coltype_defs$snv_indel_somatic_raw,
    ref_data = ref_data,
    retained_info_tags =
      settings[['conf']][['other']]$retained_vcf_info_tags,
    variant_origin = "Somatic")

  callset[['variant']] <- callset[['variant']] |>
    pcgrr::append_dbnsfp_var_link() |>
    pcgrr::append_dbmts_var_link() |>
    pcgrr::append_tcga_var_link() |>
    pcgrr::append_annotation_links() |>
    pcgrr::append_drug_var_link(ref_data = ref_data) |>
    pcgrr::append_tfbs_annotation() |>
    pcgrr::append_cancer_gene_evidence(ref_data = ref_data)

  callset <-
    pcgrr::expand_biomarker_items(
      callset = callset,
      variant_origin = "somatic")

  return(callset)


}

#' Function that reads and validates CNA or SNV/InDel TSV files
#' file from PCGR/CPSR pre-report pipeline
#'
#' @param fname Path to raw file with DNA aberrations (PCGR/CPSR)
#' @param cols column type definitions of input
#' @param ref_data reference data object
#' @param retained_info_tags VCF INFO tags to be retained in output (SNVs/InDels)
#' @param variant_origin Germline/Somatic
#'
#' @export
#'
load_dna_variants <- function(
    fname = NA,
    cols = NULL,
    ref_data = NULL,
    retained_info_tags = "None",
    variant_origin = "Somatic"){

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
  if(FALSE %in% raw_col_check){
    missing_cols <-
      compulsary_cols[!raw_col_check]
    log4r_fatal(
      paste0("Missing required columns in input file ",
             fname, " - ",
             paste(missing_cols, collapse=", ")))
  }

  cols_including_retained <- cols
  retained_cols <- NULL
  if(retained_info_tags != "None"){
    retained_cols <- stringr::str_split(
      retained_info_tags, pattern = ",")[[1]]
    for(c in retained_cols){
      if(c %in% colnames(calls_raw)){
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
  if(!is.null(retained_cols)){
    for(c in retained_cols){
      if(c %in% colnames(calls)){
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
  results[['biomarker_evidence']][['all']] <- list()
  for (elevel in pcgrr::evidence_levels) {
    results[['biomarker_evidence']][['all']][[elevel]] <- data.frame()
  }

  for (type in pcgrr::evidence_types) {
    results[['biomarker_evidence']][[type]] <- list()
    for (elevel in pcgrr::evidence_levels) {
      results[['biomarker_evidence']][[type]][[elevel]] <- data.frame()
    }
  }

  results[['retained_info_tags']] <- paste(
    retained_cols_renamed, collapse=","
  )

  if("TSG" %in% colnames(results[['variant']])){
    results[['variant']] <-
      results[['variant']] |>
      dplyr::rename(
        TUMOR_SUPPRESSOR = "TSG"
      )
  }
  if("VEP_ALL_CSQ" %in% colnames(results[['variant']])){
    results[['variant']] <-
      results[['variant']] |>
      dplyr::mutate(
        VEP_ALL_CSQ = stringr::str_replace_all(
          .data$VEP_ALL_CSQ, ",",", "
        )
      )
  }
  if("HGVSp_short" %in% colnames(results[['variant']])){
    results[['variant']] <-
      results[['variant']] |>
      dplyr::rename(
        HGVSP = "HGVSp_short"
      )
  }
  if("TSG_RANK" %in% colnames(results[['variant']])){
    results[['variant']] <-
      results[['variant']] |>
      dplyr::rename(
        TUMOR_SUPPRESSOR_RANK = "TSG_RANK"
      )
  }

  if("BIOMARKER_MATCH" %in% colnames(calls) &
     "VAR_ID" %in% colnames(calls)){

    biomarker_set <-
      calls |>
      dplyr::filter(!is.na(.data$BIOMARKER_MATCH))

    citations <- as.data.frame(
      ref_data[['biomarker']][['literature']] |>
      dplyr::select(
        c("EVIDENCE_ID",
        "LINK")
      ) |>
      tidyr::separate_rows(
        .data$EVIDENCE_ID, sep=";"
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

    if(NROW(biomarker_set) > 0){
      results[['biomarker_evidence']][['all']][['any']] <-
        as.data.frame(
          biomarker_set |>
            dplyr::select(
              c("VAR_ID",
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
            dplyr::rename(BIOMARKER_MATCH = BIOMARKER_MATCHTYPE) |>
            dplyr::mutate(BIOMARKER_RESOLUTION = dplyr::case_when(
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
            dplyr::distinct()
        )

      if(NROW(results[['biomarker_evidence']][['all']][['any']]) > 0){

        for (type in pcgrr::evidence_types) {
          results[['biomarker_evidence']][[type]][["any"]] <-
            results[['biomarker_evidence']][['all']][['any']] |>
            dplyr::filter(
              .data$VARIANT_ORIGIN == variant_origin &
                .data$EVIDENCE_TYPE == stringr::str_to_title(type))
          if (NROW(results[['biomarker_evidence']][[type]][["any"]]) > 0) {
            results[['biomarker_evidence']][[type]][["A_B"]] <-
              results[['biomarker_evidence']][[type]][["any"]] |>
              dplyr::filter(
                stringr::str_detect(
                  .data$EVIDENCE_LEVEL, "^(A|B|B1|B2):"))

            if (NROW(results[['biomarker_evidence']][[type]][["A_B"]]) > 0) {
              results[['biomarker_evidence']][[type]][["A_B"]] <-
                results[['biomarker_evidence']][[type]][["A_B"]] |>
                dplyr::arrange(
                  .data$EVIDENCE_LEVEL,
                  dplyr::desc(
                    .data$RATING))
            }

            results[['biomarker_evidence']][[type]][["C_D_E"]] <-
              results[['biomarker_evidence']][[type]][["any"]] |>
              dplyr::filter(
                stringr::str_detect(
                  .data$EVIDENCE_LEVEL, "^(C|D|E):"))

            if (NROW(results[['biomarker_evidence']][[type]][["C_D_E"]]) > 0) {
              results[['biomarker_evidence']][[type]][["C_D_E"]] <-
                results[['biomarker_evidence']][[type]][["C_D_E"]] |>
                dplyr::arrange(
                  .data$EVIDENCE_LEVEL,
                  dplyr::desc(.data$RATING))
            }
          }
        }
      }
    }
  }else{
    log4r_fatal("Input data does not contain 'BIOMARKER_MATCH' column - fatal")
  }

  return(results)

}
