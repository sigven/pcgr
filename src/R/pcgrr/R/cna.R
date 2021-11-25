
#' Function that removes copy number segments that go
#' beyond chrosomal lengths for the given assembly
#'
#' @param cna_segments_df genomic segments with CNA
#' @param genome_assembly genome assembly (grch37/grch38)
#' @param bsg BSgenome object
#'
#' @export
get_valid_chromosome_segments <-
  function(cna_segments_df,
           genome_assembly = "grch37",
           bsg = BSgenome.Hsapiens.UCSC.hg19) {

  invisible(
    assertthat::assert_that(
      !is.null(cna_segments_df) & is.data.frame(cna_segments_df),
      msg = "Argument 'cna_segments_df' must be a valid data.frame() object"))
    assertable::assert_colnames(
      cna_segments_df, c("chromosome", "segment_end"),
      only_colnames = F, quiet = T)
  assertable::assert_coltypes(
    cna_segments_df,
    list(chromosome = character(),
         segment_end = integer()),
    quiet = T)

  chromosome_lengths <-
    data.frame(chromosome = utils::head(names(GenomeInfoDb::seqlengths(bsg)), 24),
               chrom_length = utils::head(GenomeInfoDb::seqlengths(bsg), 24),
               stringsAsFactors = F, row.names = NULL)
  cna_segments_df <- as.data.frame(
    cna_segments_df %>%
      dplyr::left_join(chromosome_lengths, by = c("chromosome")) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(segment_error = .data$segment_end > .data$chrom_length)
    )
  if (nrow(dplyr::filter(cna_segments_df, .data$segment_error == T)) > 0) {
    n_removed <- nrow(dplyr::filter(cna_segments_df, .data$segment_error == T))
    log4r_warn(paste0(
      "Skipping ", n_removed,
      " copy number segments that span beyond chromosomal lengths for ",
      genome_assembly,
      " (make sure chromosomal segments are consistent with assembly)"))
  }
  cna_segments_df <- cna_segments_df %>%
    dplyr::filter(.data$segment_error == F) %>%
    dplyr::select(-c(.data$segment_error, .data$chrom_length))

  return(cna_segments_df)
}

#' Function that gets the chromosome bands of copy number segments
#'
#' @param cna_df genomic ranges object with copy number aberrations
#' @param pcgr_data pcgr data bundle object
#'
#' @export
get_cna_cytoband <- function(cna_df, pcgr_data = NULL) {

  cna_gr <-
    GenomicRanges::makeGRangesFromDataFrame(
      cna_df, keep.extra.columns = T,
      seqinfo = pcgr_data[["assembly"]][["seqinfo"]],
      seqnames.field = "chromosome",
      start.field = "segment_start",
      end.field = "segment_end",
      ignore.strand = T, starts.in.df.are.0based = T)

  cytoband_gr <- pcgr_data[["genomic_ranges"]][["cytoband"]]

  invisible(
    assertthat::assert_that("focalCNAthreshold" %in%
                              names(S4Vectors::mcols(cytoband_gr))))
  cyto_hits <- GenomicRanges::findOverlaps(
    cna_gr, cytoband_gr, type = "any", select = "all")
  ranges <- cytoband_gr[S4Vectors::subjectHits(cyto_hits)]
  S4Vectors::mcols(ranges) <-
    c(S4Vectors::mcols(ranges),
      S4Vectors::mcols(cna_gr[S4Vectors::queryHits(cyto_hits)]))
  cyto_df <- as.data.frame(S4Vectors::mcols(ranges)) %>%
    dplyr::mutate(segment_start =
                    BiocGenerics::start(ranges(cna_gr[S4Vectors::queryHits(cyto_hits)]))) %>%
    dplyr::mutate(segment_end =
                    BiocGenerics::end(ranges(cna_gr[S4Vectors::queryHits(cyto_hits)]))) %>%
    dplyr::mutate(segment_length =
                    BiocGenerics::width(ranges(cna_gr[S4Vectors::queryHits(cyto_hits)])))

  cyto_stats <- as.data.frame(
    cyto_df %>%
      dplyr::group_by(.data$SEGMENT_ID, .data$segment_length) %>%
      dplyr::summarise(
        CYTOBAND = paste(.data$name, collapse = ", "),
        chromosome_arm = paste(unique(.data$arm), collapse = ","),
        focalCNAthresholds = paste(unique(.data$focalCNAthreshold), collapse = ","),
        .groups = "drop") %>%
      dplyr::ungroup() %>%
      dplyr::mutate(
        focalCNAthresholds =
          dplyr::if_else(stringr::str_detect(.data$focalCNAthresholds, ","),
                         as.character(NA),
                         as.character(.data$focalCNAthresholds))) %>%
      dplyr::mutate(
        focalCNAthresholds =
          dplyr::if_else(!is.na(.data$focalCNAthresholds),
                         as.numeric(.data$focalCNAthresholds),
                         as.numeric(NA))) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(broad_cnv_event = .data$segment_length > .data$focalCNAthresholds) %>%
      dplyr::mutate(EVENT_TYPE = "broad") %>%
      dplyr::mutate(
        EVENT_TYPE =
          dplyr::if_else(
            !is.na(.data$broad_cnv_event) & .data$broad_cnv_event == F,
            "focal", .data$EVENT_TYPE)) %>%
      dplyr::mutate(EVENT_TYPE = dplyr::if_else(
        is.na(.data$broad_cnv_event),
        as.character("broad"), .data$EVENT_TYPE)) %>%
      dplyr::mutate(
        CYTOBAND =
          stringr::str_replace(.data$CYTOBAND, ", (\\S+, ) {0,}", " - ")) %>%
      tidyr::separate(.data$SEGMENT_ID,
                      sep = ":",
                      into = c("chrom", "start_stop"),
                      remove = F) %>%
      dplyr::mutate(CYTOBAND = paste0(.data$chrom, ":", .data$CYTOBAND)) %>%
      dplyr::select(.data$SEGMENT_ID, .data$CYTOBAND, .data$EVENT_TYPE)

  )
  cna_df <- cna_df %>%
    dplyr::left_join(cyto_stats, by = "SEGMENT_ID")

  return(cna_df)

}

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
    dplyr::filter(cna_df, .data$ONCOGENE == T & .data$TUMOR_SUPPRESSOR == F &
                    .data$MEAN_TRANSCRIPT_CNA_OVERLAP >= transcript_overlap_pct &
                    .data$LOG_R >= log_r_gain)
  onco_ts_sets[["tsgene_loss"]] <- data.frame()
  onco_ts_sets[["tsgene_loss"]] <-
    dplyr::filter(cna_df, .data$TUMOR_SUPPRESSOR == T &
                    .data$MEAN_TRANSCRIPT_CNA_OVERLAP >= transcript_overlap_pct  &
                    .data$LOG_R <= log_r_homdel)

  onco_ts_sets[["other_target"]] <- data.frame()
  onco_ts_sets[["other_target"]] <-
    dplyr::filter(cna_df, .data$TUMOR_SUPPRESSOR == F & .data$ONCOGENE == F &
                    .data$MEAN_TRANSCRIPT_CNA_OVERLAP >= transcript_overlap_pct  &
                    .data$LOG_R >= log_r_gain)

  drug_target_site <-
    pcgrr::targeted_drugs_pr_ttype(
      tumor_type,
      pcgr_data,
      ignore_on_label_early_phase = T)

  if (nrow(onco_ts_sets[["other_target"]]) > 0) {
    onco_ts_sets[["other_target"]] <- onco_ts_sets[["other_target"]] %>%
      dplyr::left_join(drug_target_site, by = "SYMBOL") %>%
      dplyr::filter(!is.na(.data$DRUGS_ON_LABEL) | !is.na(.data$DRUGS_OFF_LABEL)) %>%
      dplyr::select(-c(.data$TUMOR_SUPPRESSOR, .data$ONCOGENE,
                       .data$TUMOR_SUPPRESSOR_EVIDENCE, .data$ONCOGENE_EVIDENCE)) %>%
      dplyr::distinct() %>%
      dplyr::select(.data$CHROMOSOME, .data$SYMBOL, .data$GENE_NAME, .data$SEGMENT, .data$EVENT_TYPE,
                    .data$DRUGS_ON_LABEL, .data$DRUGS_OFF_LABEL,
                    .data$SEGMENT_LENGTH_MB, .data$LOG_R, .data$SEGMENT_ID,
                    .data$DRUGS_ON_LABEL_INDICATIONS,
                    .data$DRUGS_OFF_LABEL_INDICATIONS, .data$MEAN_TRANSCRIPT_CNA_OVERLAP,
                    .data$KEGG_PATHWAY, .data$TRANSCRIPTS) %>%
      dplyr::mutate(MEAN_TRANSCRIPT_CNA_OVERLAP =
                      paste0(.data$MEAN_TRANSCRIPT_CNA_OVERLAP, "%")) %>%
      dplyr::mutate(CNA_TYPE = "gain")

    log4r_info(
      paste0("Detected ", nrow(onco_ts_sets[["other_target"]]),
             " drug targets to amplification (log(2) ratio >= ",
             log_r_gain, "): ",
             paste0(unique(onco_ts_sets[["other_target"]]$SYMBOL),
                    collapse = ", ")))

  }

  for (t in c("oncogene_gain", "tsgene_loss")) {
    if (nrow(onco_ts_sets[[t]]) > 0) {
      onco_ts_sets[[t]] <- onco_ts_sets[[t]] %>%
        dplyr::select(-c(.data$TUMOR_SUPPRESSOR, .data$ONCOGENE)) %>%
        dplyr::distinct() %>%
        dplyr::select(.data$CHROMOSOME, .data$SYMBOL, .data$GENE_NAME,
                      .data$SEGMENT, .data$LOG_R, .data$EVENT_TYPE,
                      .data$SEGMENT_LENGTH_MB,
                      .data$MEAN_TRANSCRIPT_CNA_OVERLAP,
                      .data$KEGG_PATHWAY,
                      .data$TRANSCRIPTS, .data$SEGMENT_ID) %>%
        dplyr::mutate(
          MEAN_TRANSCRIPT_CNA_OVERLAP = paste0(
            .data$MEAN_TRANSCRIPT_CNA_OVERLAP, "%")) %>%
        dplyr::mutate(CNA_TYPE =
                        dplyr::if_else(t == "oncogene_gain" |
                                         t == "other_target",
                                       "gain", "loss"))
      if (t == "oncogene_gain") {
        log4r_info(
          paste0("Detected ", nrow(onco_ts_sets[[t]]),
                 " proto-oncogene(s) subject to amplification (log(2) ratio >= ",
                 log_r_gain, "): ",
                 paste0(unique(onco_ts_sets[[t]]$SYMBOL), collapse = ", ")))
        onco_ts_sets[[t]] <- onco_ts_sets[[t]] %>%
          dplyr::left_join(drug_target_site, by = "SYMBOL") %>%
          dplyr::select(.data$CHROMOSOME, .data$SYMBOL,
                        .data$GENE_NAME, .data$SEGMENT,
                        .data$LOG_R, .data$EVENT_TYPE,
                        .data$DRUGS_ON_LABEL,
                        .data$DRUGS_OFF_LABEL, dplyr::everything()) %>%
          dplyr::arrange(.data$DRUGS_ON_LABEL, .data$DRUGS_OFF_LABEL)
      }else{
          log4r_info(
            paste0("Detected ", nrow(onco_ts_sets[[t]]),
                   " tumor suppressor gene(s) subject to ",
                   "homozygous deletions (log(2) ratio <= ",
                   log_r_homdel, "): ",
                   paste0(unique(onco_ts_sets[[t]]$SYMBOL), collapse = ", ")))
      }
    }else{
      if (t == "tsgene_loss") {
        log4r_info(paste0("Detected 0 tumor suppressor genes subject to",
                          " homozygous deletion (log(2) ratio <= ",
                          log_r_homdel))
      }else{
        if (t == "other_target") {
          log4r_info(paste0("Detected 0 other drug targets subject to ",
                            "amplification (log(2) ratio >= ", log_r_gain))
        }else{
          log4r_info(paste0("Detected 0 proto-oncogenes subject to ",
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
        S4Vectors::mcols(ranges)) %>%
        dplyr::mutate(
          segment_start =
            as.integer(BiocGenerics::start(
              ranges(cna_gr[S4Vectors::queryHits(hits)])))) %>%
        dplyr::mutate(
          segment_end =
            as.integer(BiocGenerics::end(
              ranges(cna_gr[S4Vectors::queryHits(hits)])))) %>%
        dplyr::mutate(
          transcript_start = BiocGenerics::start(ranges)) %>%
        dplyr::mutate(
          transcript_end = BiocGenerics::end(ranges)) %>%
        dplyr::mutate(
          chrom = as.character(GenomeInfoDb::seqnames(ranges))) %>%
        dplyr::rowwise() %>%
        dplyr::mutate(
          transcript_overlap_percent =
            round(as.numeric((min(.data$transcript_end, .data$segment_end) -
                                max(.data$segment_start, .data$transcript_start)) /
                               (.data$transcript_end - .data$transcript_start)) * 100,
                  digits = 2))
    ) %>%
    pcgrr::sort_chromosomal_segments(
      chromosome_column = "chrom",
      start_segment = "segment_start",
      end_segment = "segment_end")

  return(cna_transcript_df)
}
