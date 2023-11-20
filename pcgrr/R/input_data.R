#' Function that reads and validates a fully annotated CNA file from PCGR
#' pre-report pipeline
#'
#' @param fname Path to file name
#' @param ref_data Object with reference data
#'
#' @export
load_somatic_cna <- function(fname, ref_data = NULL){

  load_dna_variants(
    fname = fname,
    cols = pcgrr::data_coltype_defs$cna_somatic_raw,
    ref_data = ref_data)

}

#' Function that reads and validates a fully annotated somatic SNV/InDel
#' file from PCGR pre-report pipeline
#'
#' @param fname Path to file name
#' @param ref_data Object with reference data
#'
#' @export
load_somatic_snv_indel <- function(fname, ref_data = NULL){

  callset <- load_dna_variants(
    fname = fname,
    cols = pcgrr::data_coltype_defs$snv_indel_somatic_raw,
    ref_data = ref_data)

  callset[['variant']] <- callset[['variant']] |>
    pcgrr::append_dbnsfp_var_link() |>
    pcgrr::append_dbmts_var_link() |>
    pcgrr::append_tfbs_annotation() |>
    pcgrr::append_tcga_var_link() |>
    pcgrr::append_annotation_links() |>
    pcgrr::append_drug_var_link()


}

#' Function that reads and validates a fully annotated somatic SNV/InDel
#' file from PCGR pre-report pipeline
#'
#' @param fname Path to raw file with DNA aberrations
#' @param cols
#'
load_dna_variants <- function(fname, cols = NULL, ref_data = NULL){

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

  calls <- suppressWarnings(
    as.data.frame(
      readr::read_tsv(
        file = fname,
        col_types =
          cols,
        na = c(".","NA"),
        progress = F
      )
    )
  )

  results <- list()
  results[['variant']] <- calls
  results[['biomarker']] <- data.frame()

  if("BIOMARKER_MATCH" %in% colnames(calls) &
     "VAR_ID" %in% colnames(calls)){

    biomarker_set <-
      calls |>
      dplyr::filter(!is.na(BIOMARKER_MATCH))

    if(NROW(biomarker_set) > 0){
      biomarker_set <- biomarker_set |>
        dplyr::select(VAR_ID, BIOMARKER_MATCH) |>
        tidyr::separate_rows(BIOMARKER_MATCH, sep=",") |>
        tidyr::separate(
          BIOMARKER_MATCH,
          into = c("BIOMARKER_SOURCE",
                   "VARIANT_ID",
                   "EVIDENCE_ITEMS",
                   "BIOMARKER_MATCHTYPE"),
          sep = "\\|"
        ) |>
        tidyr::separate_rows(
          BIOMARKER_EVIDENCE_ITEMS, sep="&") |>
        tidyr::separate(
          BIOMARKER_EVIDENCE_ITEMS,
          into = c("EVIDENCE_ID",
                   "PRIMARY_SITE",
                   "CLINICAL_SIGNIFICANCE",
                   "EVIDENCE_LEVEL",
                   "EVIDENCE_TYPE",
                   "VARIANT_ORIGIN"),
          sep = ":"
        ) |>
        dplyr::mutate(PRIMARY_SITE = stringr::str_replace_all(
          PRIMARY_SITE, "_"," "
        ))

      results[['biomarker']] <- biomarker_set
    }
  }

  return(results)

}
