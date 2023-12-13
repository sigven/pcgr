

#' Function that detects and assigns oncogenes subject to copy number
#' amplifications, and tumor suppressor genes subject to homozygous deletions
#' Also detects other drug targets subject to copy number amplifications
#'
#' @param cna_df data frame with copy number-transcript records
#' @param transcript_overlap_pct minimum level of overlap for transcripts (pct)
#' @param log_r_gain logR threshold for copy number amplifications
#' @param log_r_homdel logR threshold for homozygous deletions
#' @param tumor_type type of tumor
#' @param pcgr_data PCGR data object
#'
#'
#' @export
get_oncogene_tsgene_target_sets <- function(
  cna_df,
  transcript_overlap_pct = 100,
  log_r_gain = 0.8,
  log_r_homdel = -0.8,
  tumor_type = "Any",
  pcgr_data = NULL) {

  invisible(assertthat::assert_that(!is.null(pcgr_data)))
  invisible(assertthat::assert_that(
    is.data.frame(cna_df),
    msg = "Argument 'cna_df' must be of type data.frame"))
  assertable::assert_colnames(
    cna_df, c("ONCOGENE", "TUMOR_SUPPRESSOR",
              "MEAN_TRANSCRIPT_CNA_OVERLAP", "LOG_R",
              "SYMBOL", "KEGG_PATHWAY",
              "SEGMENT_ID", "CHROMOSOME",
              "GENE_NAME","EVENT_TYPE"),
    only_colnames = F, quiet = T)

  onco_ts_sets <- list()
  onco_ts_sets[["oncogene_gain"]] <- data.frame()
  onco_ts_sets[["oncogene_gain"]] <-
    dplyr::filter(cna_df, .data$ONCOGENE == T &
                    #.data$TUMOR_SUPPRESSOR == F &
                    .data$MEAN_TRANSCRIPT_CNA_OVERLAP >= transcript_overlap_pct &
                    .data$LOG_R >= log_r_gain)
  onco_ts_sets[["tsgene_loss"]] <- data.frame()
  onco_ts_sets[["tsgene_loss"]] <-
    dplyr::filter(cna_df, .data$TUMOR_SUPPRESSOR == T &
                    .data$MEAN_TRANSCRIPT_CNA_OVERLAP >= transcript_overlap_pct  &
                    .data$LOG_R <= log_r_homdel)

  onco_ts_sets[["other_target"]] <- data.frame()
  onco_ts_sets[["other_target"]] <-
    dplyr::filter(cna_df, .data$TUMOR_SUPPRESSOR == F &
                    .data$ONCOGENE == F &
                    .data$MEAN_TRANSCRIPT_CNA_OVERLAP >= transcript_overlap_pct  &
                    .data$LOG_R >= log_r_gain)

  drug_target_site <-
    pcgrr::targeted_drugs_pr_ttype(
      tumor_type,
      pcgr_data,
      ignore_on_label_early_phase = T)

  if (nrow(onco_ts_sets[["other_target"]]) > 0) {
    onco_ts_sets[["other_target"]] <- onco_ts_sets[["other_target"]] |>
      dplyr::left_join(drug_target_site, by = "SYMBOL") |>
      dplyr::filter(!is.na(.data$DRUGS_ON_LABEL) |
                      !is.na(.data$DRUGS_OFF_LABEL)) |>
      dplyr::select(-c(.data$TUMOR_SUPPRESSOR,
                       .data$ONCOGENE,
                       .data$TUMOR_SUPPRESSOR_EVIDENCE,
                       .data$ONCOGENE_EVIDENCE)) |>
      dplyr::distinct() |>
      dplyr::select(.data$CHROMOSOME,
                    .data$SYMBOL,
                    .data$GENE_NAME,
                    .data$SEGMENT,
                    .data$EVENT_TYPE,
                    .data$DRUGS_ON_LABEL,
                    .data$DRUGS_OFF_LABEL,
                    .data$SEGMENT_LENGTH_MB,
                    .data$LOG_R,
                    .data$SEGMENT_ID,
                    .data$OPENTARGETS_ASSOCIATIONS,
                    .data$OPENTARGETS_RANK,
                    .data$DRUGS_ON_LABEL_INDICATIONS,
                    .data$DRUGS_OFF_LABEL_INDICATIONS,
                    .data$MEAN_TRANSCRIPT_CNA_OVERLAP,
                    .data$KEGG_PATHWAY,
                    .data$TRANSCRIPTS) |>
      dplyr::mutate(MEAN_TRANSCRIPT_CNA_OVERLAP =
                      paste0(.data$MEAN_TRANSCRIPT_CNA_OVERLAP, "%")) |>
      dplyr::mutate(CNA_TYPE = "gain") |>
      dplyr::arrange(
        dplyr::desc(.data$OPENTARGETS_RANK),
        dplyr::desc(nchar(.data$DRUGS_ON_LABEL)),
        dplyr::desc(nchar(.data$DRUGS_OFF_LABEL))
      )

    pcgrr::log4r_info(
      paste0("Detected ", nrow(onco_ts_sets[["other_target"]]),
             " drug targets to amplification (log(2) ratio >= ",
             log_r_gain, "): ",
             paste0(unique(onco_ts_sets[["other_target"]]$SYMBOL),
                    collapse = ", ")))

  }

  for (t in c("oncogene_gain", "tsgene_loss")) {
    if (nrow(onco_ts_sets[[t]]) > 0) {
      onco_ts_sets[[t]] <- onco_ts_sets[[t]] |>
        dplyr::select(-c(.data$TUMOR_SUPPRESSOR, .data$ONCOGENE)) |>
        dplyr::distinct() |>
        dplyr::select(.data$CHROMOSOME,
                      .data$SYMBOL,
                      .data$GENE_NAME,
                      .data$SEGMENT,
                      .data$LOG_R,
                      .data$EVENT_TYPE,
                      .data$SEGMENT_LENGTH_MB,
                      .data$MEAN_TRANSCRIPT_CNA_OVERLAP,
                      .data$KEGG_PATHWAY,
                      .data$OPENTARGETS_ASSOCIATIONS,
                      .data$OPENTARGETS_RANK,
                      .data$TRANSCRIPTS,
                      .data$SEGMENT_ID) |>
        dplyr::mutate(
          MEAN_TRANSCRIPT_CNA_OVERLAP = paste0(
            .data$MEAN_TRANSCRIPT_CNA_OVERLAP, "%")) |>
        dplyr::mutate(CNA_TYPE =
                        dplyr::if_else(t == "oncogene_gain" |
                                         t == "other_target",
                                       "gain", "loss"))
      if (t == "oncogene_gain") {
        pcgrr::log4r_info(
          paste0("Detected ", nrow(onco_ts_sets[[t]]),
                 " proto-oncogene(s) subject to amplification (log(2) ratio >= ",
                 log_r_gain, "): ",
                 paste0(unique(onco_ts_sets[[t]]$SYMBOL), collapse = ", ")))
        onco_ts_sets[[t]] <- onco_ts_sets[[t]] |>
          dplyr::left_join(drug_target_site, by = "SYMBOL") |>
          dplyr::select(.data$CHROMOSOME,
                        .data$SYMBOL,
                        .data$GENE_NAME,
                        .data$SEGMENT,
                        .data$LOG_R,
                        .data$EVENT_TYPE,
                        .data$DRUGS_ON_LABEL,
                        .data$DRUGS_OFF_LABEL, dplyr::everything()) |>
          dplyr::arrange(
            dplyr::desc(.data$OPENTARGETS_RANK),
            dplyr::desc(nchar(.data$DRUGS_ON_LABEL)),
            dplyr::desc(nchar(.data$DRUGS_OFF_LABEL))
            )
      }else{
          pcgrr::log4r_info(
            paste0("Detected ", nrow(onco_ts_sets[[t]]),
                   " tumor suppressor gene(s) subject to ",
                   "homozygous deletions (log(2) ratio <= ",
                   log_r_homdel, "): ",
                   paste0(unique(onco_ts_sets[[t]]$SYMBOL), collapse = ", ")))
      }
    }else{
      if (t == "tsgene_loss") {
        pcgrr::log4r_info(paste0("Detected 0 tumor suppressor genes subject to",
                          " homozygous deletion (log(2) ratio <= ",
                          log_r_homdel))
      }else{
        if (t == "other_target") {
          pcgrr::log4r_info(paste0("Detected 0 other drug targets subject to ",
                            "amplification (log(2) ratio >= ", log_r_gain))
        }else{
          pcgrr::log4r_info(paste0("Detected 0 proto-oncogenes subject to ",
                            "amplification (log(2) ratio >= ", log_r_gain))
        }
      }
    }
  }

  return(onco_ts_sets)

}

#' Functions that finds all transcripts that overlap with CNA segments. A
#' new data frame with one entry per transcript-CNA overlap entry is returned
#'
#' @param cna_df data frame with copy number segments
#' @param pcgr_data PCGR list object with data
#'
#'
#' @export
get_cna_overlapping_transcripts <- function(cna_df, pcgr_data) {

  cna_gr <-
    GenomicRanges::makeGRangesFromDataFrame(
      cna_df, keep.extra.columns = T,
      seqinfo = pcgr_data[["assembly"]][["seqinfo"]],
      seqnames.field = "chromosome",
      start.field = "segment_start", end.field = "segment_end",
      ignore.strand = T, starts.in.df.are.0based = T)

  hits <- GenomicRanges::findOverlaps(
    cna_gr, pcgr_data[["genomic_ranges"]][["gencode_genes"]],
    type = "any", select = "all")
  ranges <-
    pcgr_data[["genomic_ranges"]][["gencode_genes"]][S4Vectors::subjectHits(hits)]
  S4Vectors::mcols(ranges) <-
    c(S4Vectors::mcols(ranges),
      S4Vectors::mcols(cna_gr[S4Vectors::queryHits(hits)]))

  cna_transcript_df <-
    as.data.frame(
      as.data.frame(
        S4Vectors::mcols(ranges)) |>
        dplyr::mutate(
          segment_start =
            as.integer(BiocGenerics::start(
              IRanges::ranges(cna_gr[S4Vectors::queryHits(hits)])))) |>
        dplyr::mutate(
          segment_end =
            as.integer(BiocGenerics::end(
              IRanges::ranges(cna_gr[S4Vectors::queryHits(hits)])))) |>
        dplyr::mutate(
          transcript_start = BiocGenerics::start(ranges)) |>
        dplyr::mutate(
          transcript_end = BiocGenerics::end(ranges)) |>
        dplyr::mutate(
          chrom = as.character(GenomeInfoDb::seqnames(ranges))) |>
        dplyr::mutate(
          ensembl_gene_id = stringr::str_replace(
            .data$ensembl_gene_id, "\\.[0-9]{1,}$", ""
          )
        ) |>
        dplyr::rowwise() |>
        dplyr::mutate(
          transcript_overlap_percent =
            round(as.numeric((min(.data$transcript_end, .data$segment_end) -
                                max(.data$segment_start, .data$transcript_start)) /
                               (.data$transcript_end - .data$transcript_start)) * 100,
                  digits = 2))
    ) |>
    pcgrr::sort_chromosomal_segments(
      chromosome_column = "chrom",
      start_segment = "segment_start",
      end_segment = "segment_end")

  return(cna_transcript_df)
}
