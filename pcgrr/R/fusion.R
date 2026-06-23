#' Function that checks for recurrence of RNA fusions in Mitelman database
#' and annotates with clinical and cytogenetic features
#'
#' @param query_fusions data frame with RNA fusion calls
#' @param ref_data PCGR reference data bundle (list)
#'
#' @return data frame with RNA fusion calls annotated with Mitelman database features
#' @export
rna_fusion_recurrence_mitdb <- function(
    query_fusions = NULL,
    ref_data = NULL) {

  assertthat::assert_that(!is.null(query_fusions))
  assertthat::assert_that(is.data.frame(query_fusions))
  assertthat::assert_that(NROW(query_fusions) > 0)
  assertthat::assert_that(!is.null(ref_data))
  assertthat::assert_that(is.list(ref_data))
  assertthat::assert_that(!is.null(ref_data$fusion))
  assertthat::assert_that(!is.null(ref_data$fusion$variant))
  assertthat::assert_that(is.data.frame(ref_data$fusion$variant))
  assertthat::assert_that(NROW(ref_data$fusion$variant) > 0)
  assertthat::assert_that(!is.null(ref_data$fusion$recurrent))
  assertthat::assert_that(is.data.frame(ref_data$fusion$recurrent))
  assertthat::assert_that(NROW(ref_data$fusion$recurrent) > 0)

  assertable::assert_colnames(
    query_fusions,
    c("SAMPLE_ID",
      "VAR_ID",
      "FUSION_GENE2"),
    only_colnames = FALSE,
    quiet = TRUE
  )

  assertable::assert_colnames(
    ref_data$fusion$variant,
    c("VARIANT_ALIAS"),
    only_colnames = FALSE,
    quiet = TRUE
  )

  assertable::assert_colnames(
    ref_data$fusion$recurrent,
    c("FUSION_GENE2",
      "MITDB_PRIMARY_SITE",
      "MITDB_KARYOTYPE",
      "MITDB_EVIDENCE",
      "MITDB_EVIDENCE2",
      "MITDB_NUM_EVIDENCE",
      "MITDB_PRIMARY_DIAGNOSIS",
      "MITDB_SOURCE_ID"),
    only_colnames = FALSE,
    quiet = TRUE
  )


  query_fusions <- query_fusions |>
    dplyr::select(
      c("SAMPLE_ID",
        "VAR_ID",
        "FUSION_GENE2")) |>
    dplyr::distinct()

  mitdb_hits <- query_fusions |>
    dplyr::inner_join(
      dplyr::select(
        ref_data$fusion$variant,
        c("VARIANT_ALIAS"),
      ),
      by = c("FUSION_GENE2" = "VARIANT_ALIAS"),
      relationship = "many-to-many"
    ) |>
    dplyr::distinct()

  if (NROW(mitdb_hits) > 0) {
    ## adds clinical and cytogenetic features of
    ##recurrent fusions in Mitelman database
    ## 1. MITDB_PRIMARY_SITE
    ## 2. MITDB_KARYOTYPE
    ## 3. MITDB_EVIDENCE
    ## 4. MITDB_EVIDENCE2
    ## 5. MITDB_NUM_EVIDENCE
    ## 6. MITDB_PRIMARY_DIAGNOSIS
    ## 7. MITDB_SOURCE_ID
    mitdb_hits <- mitdb_hits |>
      dplyr::left_join(
        ref_data$fusion$recurrent,
        by = "FUSION_GENE2",
        relationship = "many-to-many"
      ) |>
      dplyr::arrange(
        .data$SAMPLE_ID,
        .data$VAR_ID,
        .data$MITDB_PRIMARY_SITE,
        .data$MITDB_KARYOTYPE,
        dplyr::desc(.data$MITDB_SOURCE_ID)
      )
  }else{
    mitdb_hits <- query_fusions |>
      dplyr::mutate(
        MITDB_PRIMARY_SITE = NA_character_,
        MITDB_KARYOTYPE = NA_character_,
        MITDB_EVIDENCE = NA_character_,
        MITDB_EVIDENCE2 = NA_character_,
        MITDB_NUM_EVIDENCE = 0,
        MITDB_PRIMARY_DIAGNOSIS = NA_character_,
        MITDB_SOURCE_ID = NA_character_
      )

  }

  return(mitdb_hits)

}

#' Function that generates fusion data for PCGR report
#'
#' @param ref_data PCGR reference data object
#' @param settings PCGR run/configuration settings
#'
#' @export
generate_report_data_fusion <-
  function(ref_data = NULL,
           settings = NULL) {

    pcg_report_fusion <-
      init_fusion_content()
    pcg_report_fusion[["eval"]] <- TRUE

    if (settings$molecular_data$fname_rna_fusion_tsv != "None" &
       file.exists(settings$molecular_data$fname_rna_fusion_tsv)) {

      pcg_report_fusion[["callset"]] <-
        load_rna_fusions(
          settings = settings,
          ref_data = ref_data)
    }

    return(pcg_report_fusion)
  }




#' Find transcripts covering given splice junction breakpoints
#'
#' This function identifies transcripts that cover specified
#' junction breakpoints.
#'
#' @param bp_junctions A data frame with columns:
#'  - BP_CHROM: Chromosome of the breakpoint junction
#'  - BP_POSITION: Chromosome position of the breakpoint junction
#' @param ref_data PCGR reference data bundle (list)
#' @return A data frame with columns:
#'  - CHROM: Chromosome of the breakpoint junction
#'  - BP_POSITION: Position of the breakpoint junction
#'  - ENSEMBL_TRANSCRIPT_ID: Ensembl transcript ID covering the splice junction
#'  - ENSEMBL_GENE_ID: Ensembl gene ID
#'  - GENE_BIOTYPE: Biotype of the transcript
#'  - TRANSCRIPT_START: Start position of the transcript
#'  - TRANSCRIPT_END: End position of the transcript
#'
#' @export
#'
bp_junction_transcript_overlap <- function(
    bp_junctions = NULL,
    ref_data = NULL) {

  # Input checks
  assertthat::assert_that(!is.null(bp_junctions))
  assertthat::assert_that(NROW(bp_junctions) > 0)

  assertthat::assert_that(is.data.frame(bp_junctions))
  assertthat::assert_that(is.list(ref_data))
  assertthat::assert_that(!is.null(ref_data[['gene']]))
  assertthat::assert_that(!is.null(ref_data[['gene']][['transcript_biotype']]))
  assertthat::assert_that(is.data.frame(
    ref_data[['gene']][['transcript_biotype']]))

  assertable::assert_colnames(
    bp_junctions,
    c("BP_CHROM", "BP_POSITION"),
    only_colnames = FALSE,
    quiet = TRUE
  )

  bp_junctions <- bp_junctions |>
    dplyr::mutate(
      BP_START = .data$BP_POSITION,
      BP_END = .data$BP_POSITION
    ) |>
    dplyr::select(
      c("BP_CHROM",
      "BP_START",
      "BP_END")
    ) |>
    dplyr::distinct()

  assertable::assert_colnames(
    ref_data[['gene']][['transcript_biotype']],
    c("CHROM",
      "TRANSCRIPT_START",
      "TRANSCRIPT_END",
      "ENSEMBL_TRANSCRIPT_ID",
      "ENSEMBL_GENE_ID",
      "GENE_BIOTYPE"),
    only_colnames = FALSE,
    quiet = TRUE
  )

  # Convert data.frames to data.tables
  bp_j <- data.table::as.data.table(bp_junctions)

  transcripts <-
    ref_data[["gene"]][["transcript_biotype"]][, c(
      "CHROM",
      "ENSEMBL_GENE_ID",
      "ENSEMBL_TRANSCRIPT_ID",
      "GENCODE_TRANSCRIPT_BIOTYPE",
      "GENE_BIOTYPE",
      "MANE_SELECT2",
      "TRANSCRIPT_START",
      "TRANSCRIPT_END"
    )]

  transcripts <- unique(transcripts)
  tx <- data.table::as.data.table(transcripts)

  # Check column names match expectations
  data.table::setnames(
    bp_j, c("BP_CHROM", "BP_START","BP_END"))
  data.table::setnames(
    tx, c("CHROM",
          "ENSEMBL_GENE_ID",
          "ENSEMBL_TRANSCRIPT_ID",
          "GENCODE_TRANSCRIPT_BIOTYPE",
          "GENE_BIOTYPE",
          "MANE_SELECT2",
          "TRANSCRIPT_START",
          "TRANSCRIPT_END"))

  # Join on chromosome and transcript start and end
  data.table::setkey(
    tx, "CHROM", "TRANSCRIPT_START", "TRANSCRIPT_END")
  data.table::setkey(
    bp_j, "BP_CHROM", "BP_START", "BP_END")

  # Find overlaps between breakpoint junction and
  # transcript coordinates - using a fast binary-search based
  # overlap join
  tx <- data.table::as.data.table(tx)
  stopifnot(data.table::is.data.table(tx))
  bp_j <- data.table::as.data.table(bp_j)
  stopifnot(data.table::is.data.table(bp_j))

  bp_tx_joined <- data.table::foverlaps(
    x = bp_j,
    y = tx,
    by.x =
      c("BP_CHROM", "BP_START", "BP_END"),
    by.y =
      c("CHROM", "TRANSCRIPT_START", "TRANSCRIPT_END"),
    nomatch = 0L
  )

  final_bp_tx <- data.frame()
  if (NROW(bp_tx_joined) > 0) {

    final_bp_tx <-
      as.data.frame(bp_tx_joined) |>
      dplyr::rename(
        BP_POSITION = .data$BP_START
      ) |>
      dplyr::select(
        c("BP_CHROM",
        "BP_POSITION",
        "ENSEMBL_TRANSCRIPT_ID",
        "ENSEMBL_GENE_ID",
        "GENE_BIOTYPE",
        "MANE_SELECT2",
        "TRANSCRIPT_START",
        "TRANSCRIPT_END")
      ) |>
      dplyr::distinct()

  }

  return(final_bp_tx)
}


#' Identify druggable fusion partners
#'
#' This function identifies druggable fusion partners from a given data frame of fusion genes.
#' It uses a reference data bundle to append targeted drug annotations based on the specified primary site.
#'
#' @param df A data frame containing fusion gene information with columns "VAR_ID" and "FUSION_GENE_5P" or "FUSION_GENE_3P".
#' @param partner A character string indicating whether to analyze the 5' or 3
#' ' fusion partner. Default is "5P".
#' @param ref_data A list containing reference data, including drug annotations.
#' @param variant_display A logical value indicating whether to use variant-level drug annotations. Default
#' is FALSE.
#' @param primary_site A character string specifying the primary site for drug annotation. Default is
#' "Lung".
#' @return A data frame with druggable fusion partners and their associated targeted drug annotations.
#' @export
#'
get_druggable_fusion_partner <- function(
    df = NULL,
    partner = "5P",
    ref_data = NULL,
    variant_display = FALSE,
    primary_site = "Lung") {

  assertthat::assert_that(is.character(partner))
  assertthat::assert_that(partner %in% c("5P", "3P"))
  assertthat::assert_that(!is.null(ref_data))
  assertthat::assert_that(is.list(ref_data))
  assertthat::assert_that(!is.null(df))
  assertthat::assert_that(is.data.frame(df))
  gene_col <- paste0("FUSION_GENE_", partner)
  assertable::assert_colnames(
    df, c("VAR_ID", gene_col),
    only_colnames = FALSE, quiet = TRUE
  )

  if ("SYMBOL" %in% colnames(df)) {
    df$SYMBOL <- NULL
  }

  ## append_targeted_drug_annotations require a 'VAR_ID'
  ## and 'SYMBOL' column
  df <- df |>
    dplyr::rename(SYMBOL = dplyr::all_of(gene_col))

  targeted_drugs <- append_targeted_drug_annotations(
    var_df = df,
    primary_site = primary_site,
    ref_data = ref_data)

  targeted_drugs_display <-
    append_drug_var_link(
      var_df = df,
      primary_site = primary_site,
      ref_data = ref_data
    )

  if (variant_display == TRUE) {
    if ("TARGETED_INHIBITORS_ALL" %in% colnames(targeted_drugs_display) &
       "TARGETED_INHIBITORS" %in% colnames(targeted_drugs_display)) {

      targeted_drugs_display[,paste0("TARGETED_INHIBITORS_ALL_", partner)] <-
        targeted_drugs_display$TARGETED_INHIBITORS_ALL
      targeted_drugs_display[,paste0("TARGETED_INHIBITORS_", partner)] <-
        targeted_drugs_display$TARGETED_INHIBITORS
      targeted_drugs_display[,paste0("FUSION_GENE_", partner)] <-
        targeted_drugs_display$SYMBOL

      for (e in c("SYMBOL",
                 "TARGETED_INHIBITORS_ALL",
                 "TARGETED_INHIBITORS")) {
        if (e %in% colnames(targeted_drugs_display)) {
          targeted_drugs_display[,e] <- NULL
        }
      }
    }

  }else{
    if ("TARGETED_INHIBITORS_ALL2" %in% colnames(targeted_drugs) &
       "TARGETED_INHIBITORS2" %in% colnames(targeted_drugs)) {

      targeted_drugs[,paste0("TARGETED_INHIBITORS_ALL2_", partner)] <-
        targeted_drugs$TARGETED_INHIBITORS_ALL2
      targeted_drugs[,paste0("TARGETED_INHIBITORS2_", partner)] <-
        targeted_drugs$TARGETED_INHIBITORS2
      targeted_drugs[,paste0("FUSION_GENE_", partner)] <-
        targeted_drugs$SYMBOL

      for (e in c("SYMBOL",
                 "TARGETED_INHIBITORS_ALL2",
                 "TARGETED_INHIBITORS2")) {
        if (e %in% colnames(targeted_drugs)) {
          targeted_drugs[,e] <- NULL
        }
      }
    }
  }

  if (variant_display == TRUE) {
    targeted_drugs <- targeted_drugs_display
  }

  return(targeted_drugs)

}
