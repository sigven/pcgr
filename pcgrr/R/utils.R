#' Function that checks whether a set of column names are present in two
#' different data frames
#'
#' @param df1 data frame 1
#' @param df2 data frame 2
#' @param cnames character vector with column names to be check for presence
#'
#' @return existing_common_columns T/F
#'
#' @export
check_common_colnames <- function(df1 = NULL, df2 = NULL, cnames = NULL) {

  invisible(assertthat::assert_that(
    is.data.frame(df1),
    msg = "Object df1 is not of type data.frame"))
  invisible(assertthat::assert_that(
    is.data.frame(df2), msg = "Object df2 is not of type data.frame"))

  existing_common_columns <- F
  colnames_in_df1 <- unique(cnames %in% colnames(df1))
  colnames_in_df2 <- unique(cnames %in% colnames(df2))

  if (length(colnames_in_df1) == 1 & length(colnames_in_df2) == 1) {
    if (colnames_in_df1 == T & colnames_in_df2 == T) {
      existing_common_columns <- T
    }
  }
  return(existing_common_columns)

}

#' Function that removes column(s) from data frame
#'
#' @param df data.frame with data
#' @param cnames character vector with column names
#'
#' @return df data.frame with columns removed (if present originally)
#'
#' @export
remove_cols_from_df <- function(df, cnames = NULL) {

  invisible(assertthat::assert_that(
    is.data.frame(df),
    msg = "Object df is not of type data.frame"))
  if (!is.null(cnames)) {
    for (c in cnames) {
      col_name <- as.character(c)
      if (col_name %in% colnames(df)) {
        df[, col_name] <- NULL
      }
    }
  }
  return(df)

}

#' Function that plots a histogram of the the variant allelic
#' support (tumor)
#'
#' @param var_df data frame with somatic mutations
#' @param bin_size size of bins for allelic frequency
#'
#' @return p geom_histogram plot from ggplot2
#'
#' @export
af_distribution <- function(var_df, bin_size = 0.05) {

  af_bin_df <- data.frame()
  assertable::assert_colnames(
    var_df, c("VAF_TUMOR","VARIANT_CLASS"),
    only_colnames = F, quiet = T)
  var_df <- var_df |>
    dplyr::select(c("VAF_TUMOR", "VARIANT_CLASS"))

  # if(NROW(var_df) <= 200){
  #   bin_size = 0.10
  # }
  if(NROW(var_df) > 500){
    bin_size = 0.025
  }

  i <- 1
  num_bins <- as.integer(1 / bin_size)
  bin_start <- 0
  while (i <= num_bins) {
    bin_end <- bin_start + bin_size
    bin_name <- as.character(paste0(bin_start, " - ", bin_end))
    for(e in c('SNV','deletion','insertion','indel','substitution')){
      df <- data.frame(bin_name = bin_name,
                       bin_start = bin_start,
                       bin_end = bin_end,
                       bin = as.integer(i),
                       VARIANT_CLASS = e, stringsAsFactors = F)
      af_bin_df <- dplyr::bind_rows(
        af_bin_df, df)
    }
    bin_start <- bin_end
    i <- i + 1
  }

  var_df_trans <- var_df |>
    dplyr::mutate(
      bin = cut(.data$VAF_TUMOR,
                breaks = seq(0,1,bin_size),
                right = F, include.lowest = T, labels = F))

  var_df_trans_bin <- as.data.frame(
    dplyr::group_by(
      var_df_trans, .data$VARIANT_CLASS, .data$bin) |>
      dplyr::summarise(Count = dplyr::n(),
                       .groups = "drop"))
  af_bin_df <- af_bin_df |>
    dplyr::left_join(
      var_df_trans_bin, by = c("bin", "VARIANT_CLASS")) |>
    dplyr::mutate(Count = dplyr::if_else(
      is.na(.data$Count),
      as.numeric(0),
      as.numeric(.data$Count)))

  af_dist_final <- af_bin_df
  for(vclass in unique(af_bin_df$VARIANT_CLASS)){
    sum_count <- sum(
      af_bin_df$Count[af_bin_df$VARIANT_CLASS == vclass])
    if(sum_count == 0){
      af_dist_final <- af_dist_final |>
        dplyr::filter(
          .data$VARIANT_CLASS != vclass)
    }
  }

  return(af_dist_final)

}


#' Function that plots a histogram of the the variant allelic
#' support (tumor) - grouped by tiers
#'
#' @param tier_df data frame with somatic mutations
#' @param bin_size size of bins for allelic frequency
#'
#' @return p geom_histogram plot from ggplot2
#'
#' @export
tier_af_distribution <- function(tier_df, bin_size = 0.05) {

  af_bin_df <- data.frame()
  assertable::assert_colnames(
    tier_df, c("VAF_TUMOR","ACTIONABILITY_TIER"),
    only_colnames = F, quiet = T)
  tier_df <- tier_df |>
    dplyr::select(c("VAF_TUMOR", "ACTIONABILITY_TIER")) |>
    dplyr::mutate(
      TIER = paste0("TIER ", .data$ACTIONABILITY_TIER))

  if(NROW(tier_df) > 500){
    bin_size = 0.025
  }

  i <- 1
  num_bins <- as.integer(1 / bin_size)
  bin_start <- 0
  while (i <= num_bins) {
    bin_end <- bin_start + bin_size
    bin_name <- as.character(paste0(bin_start, " - ", bin_end))
    j <- 1
    while (j <= 5) {
      TIER <- paste0("TIER ", j)
      df <- data.frame(bin_name = bin_name,
                       bin_start = bin_start,
                       bin_end = bin_end,
                       bin = as.integer(i),
                       TIER = TIER, stringsAsFactors = F)
      af_bin_df <- rbind(af_bin_df, df)
      j <- j + 1
    }
    bin_start <- bin_end
    i <- i + 1
  }

  tier_df_trans <- tier_df |>
    dplyr::mutate(
      bin = cut(.data$VAF_TUMOR,
                breaks = seq(0,1,bin_size),
                right = F, include.lowest = T, labels = F))

  tier_df_trans_bin <- as.data.frame(
    dplyr::group_by(
      tier_df_trans, .data$TIER, .data$bin) |>
      dplyr::summarise(Count = dplyr::n(),
                       .groups = "drop"))
  af_bin_df <- af_bin_df |>
    dplyr::left_join(
      tier_df_trans_bin, by = c("bin", "TIER")) |>
    dplyr::mutate(Count = dplyr::if_else(
      is.na(.data$Count),
      as.numeric(0),
      as.numeric(.data$Count)))

  return(af_bin_df)

}

#' Checks for valid chromosome names in data frame of variants
#'
#' @param vcf_data_df data frame
#' @param chromosome_column name of chromosome column
#' @param bsg BSGenome object
#'
#' @return vcf_data_df valid data frame with valid mutations
#'
#' @export
get_valid_chromosomes <- function(vcf_data_df,
                                  chromosome_column = "CHROM",
                                  bsg = NULL) {
  assertthat::assert_that(
    is.data.frame(vcf_data_df),
    msg = paste0("Argument 'vcf_data_df' must be of type data.frame, not ",
                 class(vcf_data_df)))
  assertthat::assert_that(!is.null(bsg),
                          msg = "Please provide a valid BSgenome.Hsapiens object")
  assertable::assert_colnames(
    vcf_data_df, chromosome_column, only_colnames = FALSE, quiet = TRUE)
  vcf_data_df_valid <- vcf_data_df
  vcf_data_df_valid[, chromosome_column] <-
    factor(vcf_data_df_valid[, chromosome_column])
  levels(vcf_data_df_valid[, chromosome_column]) <-
    sub("^([0-9XY])", "chr\\1", levels(vcf_data_df_valid[, chromosome_column]))
  levels(vcf_data_df_valid[, chromosome_column]) <-
    sub("^MT", "chrM", levels(vcf_data_df_valid[, chromosome_column]))
  levels(vcf_data_df_valid[, chromosome_column]) <-
    sub("^(GL[0-9]+).[0-9]", "chrUn_\\L\\1",
        levels(vcf_data_df_valid[, chromosome_column]), perl = TRUE)
  unknown_regs <-
    levels(vcf_data_df_valid[, chromosome_column])
  unknown_regs <- unknown_regs[which(
    !(levels(vcf_data_df_valid[, chromosome_column]) %in%
        GenomeInfoDb::seqnames(bsg)))]
  if (length(unknown_regs) > 0) {
    unknown_regs <- paste(unknown_regs, collapse = ",\ ")
    log4r_warn(paste(
      "Check chr names -- not all match BSgenome.Hsapiens object:\n",
      unknown_regs, sep = " "))
    vcf_data_df_valid <-
      vcf_data_df_valid[vcf_data_df_valid[, chromosome_column]
                        %in% GenomeInfoDb::seqnames(bsg), ]
  }
  vcf_data_df_valid[, chromosome_column] <-
    as.character(vcf_data_df_valid[, chromosome_column])
  return(vcf_data_df_valid)

}

#' Function that excludes genomic aberrations from non-nuclear chromosomes
#'
#' @param vcf_df data frame
#' @param chrom_var variable name of chromosome in data frame
#' @return vcf_df data frame with mutations from nuclear chromosomes only
#'
#' @export
get_ordinary_chromosomes <- function(vcf_df, chrom_var = "CHROM") {
  invisible(assertthat::assert_that(
    is.data.frame(vcf_df),
    msg = "Argument 'vcf_df' must be of type data.frame"))
  assertable::assert_colnames(
    vcf_df, chrom_var, only_colnames = F, quiet = T)
  vcf_df <- vcf_df |>
    dplyr::mutate(
      !!rlang::sym(chrom_var) := as.character(!!rlang::sym(chrom_var)))
  n_before_exclusion <- nrow(vcf_df)
  nuc_chromosomes_df <- data.frame(c(as.character(seq(1:22)), "X", "Y"),
                                   stringsAsFactors = F)
  colnames(nuc_chromosomes_df) <- c(chrom_var)
  vcf_df <- dplyr::semi_join(vcf_df, nuc_chromosomes_df, by = chrom_var)
  n_after_exclusion <- nrow(vcf_df)
  pcgrr::log4r_info(
    paste0("Excluding ",
           n_before_exclusion - n_after_exclusion,
           " variants from non-nuclear chromosomes/scaffolds"))
  return(vcf_df)

}

#' Function that orders genomic aberrations according to order
#' of chromosomes and chromosomal position
#'
#' @param vcf_df data frame
#' @param chrom_var variable name of chromosome in data frame
#' @param pos_var variable name for chromosomal position
#' @return vcf_df data frame with ordered mutations
#'
#' @export
order_variants <- function(
    vcf_df, chrom_var = "CHROM", pos_var = "POS") {
  stopifnot(is.data.frame(vcf_df) &
              chrom_var %in% colnames(vcf_df) &
              pos_var %in% colnames(vcf_df))
  if (nrow(vcf_df) == 0)return(vcf_df)
  vcf_df |>
    dplyr::mutate(!!rlang::sym(chrom_var) :=
                    factor(!!rlang::sym(chrom_var),
                           ordered = T,
                           levels = c(as.character(seq(1:22)), "X", "Y"))) |>
    dplyr::arrange(!!rlang::sym(chrom_var), !!rlang::sym(pos_var)) |>
    dplyr::mutate(!!rlang::sym(chrom_var) :=
                    as.character(!!rlang::sym(chrom_var)))
}


#' Function that sorts chromosomal segments according to chromosome
#' and chromosomal start/end position
#'
#' @param df data frame with chromosome and start + end segment
#' @param chromosome_column name of column for chromosome name is sigven
#' @param start_segment name of column that indicates start of
#' chromosomal segment
#' @param end_segment name of column that indicates end of chromosomal segment
#' @return df_final data frame with sorted chromosomal segments
#'
#' @export
sort_chromosomal_segments <- function(df,
                                      chromosome_column = "CHROM",
                                      start_segment = "POS",
                                      end_segment = "POS") {

  invisible(assertthat::assert_that(
    !is.null(df),
    msg = "Argument 'df' must be a non-NULL object"))
  invisible(assertthat::assert_that(
    is.data.frame(df),
    msg = paste0("Argument 'df' must be of type data.frame, not ", class(df))))
  assertable::assert_colnames(
    df, c(chromosome_column, start_segment, end_segment),
    only_colnames = F, quiet = T)
  if (nrow(df) == 0) {
    return(df)
  }


  df[, start_segment] <- as.integer(df[, start_segment])
  df[, end_segment] <- as.integer(df[, end_segment])
  df_sorted <- df

  chr_prefix <- FALSE
  chromosome_names <- unique(df[, chromosome_column])
  for (m in chromosome_names) {
    if (startsWith(m, "chr")) {
      chr_prefix <- TRUE
    }
  }

  chr_order <- c(as.character(paste0("chr", c(1:22))), "chrX", "chrY")
  if (chr_prefix == FALSE) {
    chr_order <- c(as.character(c(1:22)), "X", "Y")
  }
  df_sorted[, chromosome_column] <-
    factor(df_sorted[, chromosome_column], levels = chr_order)
  df_sorted <- df_sorted[order(df_sorted[, chromosome_column]), ]

  df_final <- NULL
  for (chrom in chr_order) {
    if (nrow(df_sorted[!is.na(df_sorted[, chromosome_column]) &
                       df_sorted[, chromosome_column] == chrom, ]) > 0) {
      chrom_regions <- df_sorted[df_sorted[, chromosome_column] == chrom, ]
      chrom_regions_sorted <-
        chrom_regions[with(chrom_regions,
                           order(chrom_regions[, start_segment],
                                 chrom_regions[, end_segment])), ]
      df_final <- rbind(df_final, chrom_regions_sorted)
    }
  }
  return(df_final)
}


#' Function that performs stringr::str_replace on strings of multiple
#' string columns of a dataframe
#'
#' @param df data frame
#' @param strings name of columns for which string replace is to be performed
#' @param pattern pattern to replace
#' @param replacement string to replace
#' @param replace_all logical - replace all occurrences
#' @return df
#'
#'
#' @export
df_string_replace <- function(df, strings, pattern,
                              replacement, replace_all = F) {
  stopifnot(is.data.frame(df))
  for (column_name in strings) {
    if (column_name %in% colnames(df)) {
      if (replace_all == F) {
        df[, column_name] <-
          stringr::str_replace(df[, column_name],
                               pattern = pattern,
                               replacement = replacement)
      }else{
        df[, column_name] <-
          stringr::str_replace_all(df[, column_name],
                                   pattern = pattern,
                                   replacement = replacement)
      }
    }
  }
  return(df)
}

#' Function that generate stats for a given variant set, considering
#' number of variants/genes affected across tiers, types of variants
#' ()
#'
#' @param callset list object with callset data (CNA or SNVs/InDels)
#' @param name type of variant statistic
#' @param vartype type of variant ('snv_indel', 'cna')
#'
#' @export
variant_stats_report <- function(
    callset = NULL,
    name = "vstats",
    vartype = "snv_indel") {

  call_stats <- list()
  call_stats[[name]] <- list()

  if(vartype == 'snv_indel'){
    for (n in c("n",
                "n_snv",
                "n_sub",
                "n_indel",
                "n_coding",
                "n_noncoding")) {
      call_stats[[name]][[n]] <- 0
    }

    if("VARIANT_CLASS" %in% colnames(callset$variant)){

      call_stats[[name]][["n"]] <-
        callset$variant |>
        nrow()
      if(call_stats[[name]][["n"]] > 0){
        call_stats[[name]][["n_snv"]] <-
          callset$variant |>
          dplyr::filter(.data$VARIANT_CLASS == "SNV") |>
          nrow()
        call_stats[[name]][["n_indel"]] <-
          callset$variant |>
          dplyr::filter(
            .data$VARIANT_CLASS == "insertion" |
              .data$VARIANT_CLASS == "deletion" |
              .data$VARIANT_CLASS == "indel") |>
          nrow()
        call_stats[[name]][["n_sub"]] <-
          callset$variant |>
          dplyr::filter(
            .data$VARIANT_CLASS == "substituion") |>
          nrow()
        call_stats[[name]][["n_coding"]] <-
          callset$variant |>
          dplyr::filter(.data$CODING_STATUS == "coding") |>
          nrow()
        call_stats[[name]][["n_noncoding"]] <-
          callset$variant |>
          dplyr::filter(.data$CODING_STATUS == "noncoding") |>
          nrow()
      }
    }
  }

  if("ACTIONABILITY_TIER" %in% colnames(callset$variant)){

    for (n in c("n_actionable_tier1",
                "n_actionable_tier2",
                "n_actionable_tier3",
                "n_tier4",
                "n_tier5")) {
      call_stats[[name]][[n]] <- 0
    }


    ## only consider biomarker variants with a certain level of
    ## resolution (i.e. not gene level) when calculating tier statistics
    ## for SNVs/InDels

    if("biomarker_evidence" %in% names(callset)){
      biomarkers_for_stats <- data.frame()
      if(NROW(callset$biomarker_evidence$items) > 0){
        biomarkers_for_stats <-
          dplyr::select(
            callset$biomarker_evidence$items,
            c("ACTIONABILITY_TIER", "BM_RESOLUTION",
              "VAR_ID","ENTREZGENE")) |>
          dplyr::filter(!is.na(.data$ACTIONABILITY_TIER)) |>
          dplyr::distinct()

        if(NROW(biomarkers_for_stats) > 0){
          if(vartype == "snv_indel"){
            biomarkers_for_stats <-
              biomarkers_for_stats |>
              dplyr::filter(
                !is.na(.data$ACTIONABILITY_TIER) &
                  .data$BM_RESOLUTION != "gene") |>
              dplyr::select(c("VAR_ID","ACTIONABILITY_TIER","ENTREZGENE")) |>
              dplyr::distinct()
          }
          if(NROW(biomarkers_for_stats) > 0){
            for(i in 1:2){
              call_stats[[name]][[paste0("n_actionable_tier",i)]] <-
                callset$variant |>
                dplyr::filter(
                  !is.na(.data$ACTIONABILITY_TIER) &
                    .data$ACTIONABILITY_TIER == i) |>
                dplyr::inner_join(
                  biomarkers_for_stats,
                  by = c("ACTIONABILITY_TIER","VAR_ID","ENTREZGENE")) |>
                NROW()
            }
          }
        }
      }
      call_stats[[name]][["n_actionable_tier3"]] <-
        callset$variant |>
        dplyr::filter(
          !is.na(.data$ACTIONABILITY_TIER) &
            .data$ACTIONABILITY_TIER == 3) |>
        NROW()

      for(i in 4:5){
        call_stats[[name]][[paste0("n_tier",i)]] <-
          callset$variant |>
          dplyr::filter(
            !is.na(.data$ACTIONABILITY_TIER) &
              .data$ACTIONABILITY_TIER == i) |>
          NROW()
      }
      if (vartype == "cna"){
        call_stats[[name]][["n_tier4"]] <- NULL
        call_stats[[name]][["n_tier5"]] <- NULL
      }

      if(vartype == "snv_indel"){
        if(call_stats[[name]][["n_actionable_tier3"]] > 0){
          call_stats[[name]][["n_actionable_tier3_tsg"]] <-
            callset$variant |>
            dplyr::filter(
              !is.na(.data$ACTIONABILITY_TIER) &
                .data$ACTIONABILITY_TIER == 3) |>
            dplyr::filter(
              .data$TUMOR_SUPPRESSOR == TRUE &
                .data$ONCOGENE == FALSE) |>
            NROW()

          call_stats[[name]][["n_actionable_tier3_oncogene"]] <-
            callset$variant |>
            dplyr::filter(
              !is.na(.data$ACTIONABILITY_TIER) &
                .data$ACTIONABILITY_TIER == 3) |>
            dplyr::filter(
              .data$TUMOR_SUPPRESSOR == FALSE &
                .data$ONCOGENE == TRUE) |>
            NROW()

          call_stats[[name]][["n_actionable_tier3_dualrole"]] <-
            callset$variant |>
            dplyr::filter(
              !is.na(.data$ACTIONABILITY_TIER) &
                .data$ACTIONABILITY_TIER == 3) |>
            dplyr::filter(
              .data$TUMOR_SUPPRESSOR == TRUE &
                .data$ONCOGENE == TRUE) |>
            NROW()

        }else{
          call_stats[[name]][["n_actionable_tier3_tsg"]] <- 0
          call_stats[[name]][["n_actionable_tier3_oncogene"]] <- 0
          call_stats[[name]][["n_actionable_tier3_dualrole"]] <- 0
        }
      }

      for (n in c("n_eitems_diagnostic_tier1",
                  "n_eitems_predictive_tier1",
                  "n_eitems_prognostic_tier1",
                  "n_eitems_diagnostic_tier2",
                  "n_eitems_predictive_tier2",
                  "n_eitems_prognostic_tier2",
                  "n_genes_tier1",
                  "n_genes_tier2",
                  "n_genes_tier3")) {
        call_stats[[name]][[n]] <- 0
      }

      if(NROW(callset$biomarker_evidence$items) > 0){
        if("BM_EVIDENCE_TYPE" %in% colnames(callset$biomarker_evidence$items) &
           "ACTIONABILITY_TIER" %in% colnames(callset$biomarker_evidence$items) &
           "ENTREZGENE" %in% colnames(callset$biomarker_evidence$items) &
           "BM_RESOLUTION" %in% colnames(callset$biomarker_evidence$items)){

          for(tier in c(1,2)){
            for(etype in c("Diagnostic","Predictive","Prognostic")){
              stat <- paste0("n_eitems_",tolower(etype),"_tier",tier)

              if(vartype == "snv_indel"){
                call_stats[[name]][[stat]] <-
                  callset$biomarker_evidence$items |>
                  dplyr::filter(
                    .data$BM_RESOLUTION != "gene" &
                      !is.na(.data$ACTIONABILITY_TIER) &
                      .data$ACTIONABILITY_TIER == tier &
                      .data$BM_EVIDENCE_TYPE == etype) |>
                  NROW()
              }else{
                call_stats[[name]][[stat]] <-
                  callset$biomarker_evidence$items |>
                  dplyr::filter(
                    !is.na(.data$ACTIONABILITY_TIER) &
                      .data$ACTIONABILITY_TIER == tier &
                      .data$BM_EVIDENCE_TYPE == etype) |>
                  NROW()
              }
            }
            stat <- paste0("n_genes_tier",tier)
            if(vartype == "snv_indel"){
              call_stats[[name]][[stat]] <-
                callset$biomarker_evidence$items |>
                dplyr::filter(
                  .data$BM_RESOLUTION != "gene" &
                    .data$ACTIONABILITY_TIER == tier) |>
                dplyr::select(.data$ENTREZGENE) |>
                dplyr::distinct() |>
                NROW()
            }else{
              call_stats[[name]][[stat]] <-
                callset$biomarker_evidence$items |>
                dplyr::filter(
                  .data$ACTIONABILITY_TIER == tier) |>
                dplyr::select(.data$ENTREZGENE) |>
                dplyr::distinct() |>
                NROW()
            }
          }
        }
      }
    }
  }

  if(vartype == 'cna' &
     "ACTIONABILITY_TIER" %in% colnames(callset$variant) &
     "VAR_ID" %in% colnames(callset$variant) &
     "VARIANT_CLASS" %in% colnames(callset$variant)){
    for (n in c("n_tsg_loss",
                "n_oncogene_gain",
                "n_other_drugtarget_gain")) {
      call_stats[[name]][[n]] <- 0
    }
    call_stats[[name]][["n_tsg_loss"]] <-
      callset$variant |>
      dplyr::filter(
        !is.na(.data$ACTIONABILITY_TIER) &
          .data$ACTIONABILITY_TIER == 3 &
          .data$VARIANT_CLASS == "homdel") |>
      nrow()
    call_stats[[name]][["n_oncogene_gain"]] <-
      callset$variant |>
      dplyr::filter(
        !is.na(.data$ACTIONABILITY_TIER) &
          .data$ACTIONABILITY_TIER == 3 &
          .data$VARIANT_CLASS == "gain") |>
      NROW()

    if("TARGETED_INHIBITORS_ALL2" %in% colnames(callset$variant)){
      call_stats[[name]][["n_other_drugtarget_gain"]] <-
        callset$variant |>
        dplyr::filter(
          is.na(.data$ACTIONABILITY_TIER) &
            .data$VARIANT_CLASS == "gain" &
            !is.na(.data$TARGETED_INHIBITORS_ALL2)) |>
        NROW()
    }


    call_stats[[name]][["n_genes_tier3"]] <-
      callset$variant |>
      dplyr::filter(
        is.na(.data$ACTIONABILITY_TIER) &
          .data$ACTIONABILITY_TIER == 3) |>
      NROW()

    call_stats[[name]][["n_segments_gain"]] <-
      callset$variant |>
      dplyr::filter(.data$VARIANT_CLASS == "gain") |>
      dplyr::select(.data$VAR_ID) |>
      dplyr::distinct() |>
      NROW()

    call_stats[[name]][["n_segments_loss"]] <-
      callset$variant |>
      dplyr::filter(.data$VARIANT_CLASS == "homdel") |>
      dplyr::select(.data$VAR_ID) |>
      dplyr::distinct() |>
      NROW()
  }

  if(vartype == 'snv_indel' &
     "FINAL_CLASSIFICATION" %in% colnames(callset$variant)){
    if("BM_EVIDENCE_TYPE" %in% colnames(callset$variant) &
       "GENOMIC_CHANGE" %in% colnames(callset$variant)){

      call_stats[[name]][['n_eitems_predictive']] <- callset$variant |>
        dplyr::filter(.data$BM_EVIDENCE_TYPE == "Predictive") |>
        NROW()
      call_stats[[name]][['n_eitems_prognostic']] <- callset$variant |>
        dplyr::filter(.data$BM_EVIDENCE_TYPE == "Prognostic") |>
        NROW()
      call_stats[[name]][['n_eitems_diagnostic']] <- callset$variant |>
        dplyr::filter(.data$BM_EVIDENCE_TYPE == "Diagnostic") |>
        NROW()
      call_stats[[name]][['n_eitems_predisposing']] <- callset$variant |>
        dplyr::filter(.data$BM_EVIDENCE_TYPE == "Predisposing") |>
        NROW()

      call_stats[[name]][['n_var_eitems']] <- callset$variant |>
        dplyr::filter(!is.na(.data$BM_EVIDENCE_TYPE)) |>
        dplyr::select("GENOMIC_CHANGE") |>
        dplyr::distinct() |>
        NROW()
    }

    call_stats[[name]][['n_p']] <- callset$variant |>
      dplyr::filter(.data$FINAL_CLASSIFICATION == "Pathogenic") |>
      NROW()
    call_stats[[name]][['n_lp']] <- callset$variant |>
      dplyr::filter(.data$FINAL_CLASSIFICATION == "Likely_Pathogenic") |>
      NROW()
    call_stats[[name]][['n_vus']] <- callset$variant |>
      dplyr::filter(.data$FINAL_CLASSIFICATION == "VUS") |>
      NROW()
    call_stats[[name]][['n_lb']] <- callset$variant |>
      dplyr::filter(.data$FINAL_CLASSIFICATION == "Likely_Benign") |>
      NROW()
    call_stats[[name]][['n_b']] <- callset$variant |>
      dplyr::filter(.data$FINAL_CLASSIFICATION == "Benign") |>
      NROW()

  }

  return(call_stats)
}


#' Function that filters variant set on (depth, allelic fraction)
#' for tumor and normal and filters according to settings
#'
#' @param vcf_df data frame with variants
#' @param config list with PCGR configuration settings
#' @param precision number of significant digits for allelic fraction estimation
#'
#' @return vcf_df
#'
#' @export
filter_read_support <- function(vcf_df, config = NULL, precision = 3) {

  pcgrr::log4r_info(paste0(
    paste0("Filtering tumor variants based on allelic depth/fraction (min_dp_tumor=",
           config$allelic_support$tumor_dp_min,
           ", min_af_tumor=",
           config$somatic_snv$allelic_support$tumor_af_min, ")"))
  )

  pcgrr::log4r_info(paste0(
    "Filtering tumor variants based on allelic depth/fraction (min_dp_control=",
    config$somatic_snv$allelic_support$control_dp_min,
    ", max_af_control=",
    config$allelic_support$control_af_max, ")"))

  n_before_dp_af_filtering <- nrow(vcf_df)
  if (!any(is.na(vcf_df$DP_TUMOR))) {
    vcf_df <- dplyr::filter(
      vcf_df,
      .data$DP_TUMOR >= config$somatic_snv$allelic_support$tumor_dp_min)
  }
  if (!any(is.na(vcf_df$VAF_TUMOR))) {
    vcf_df <- dplyr::filter(
      vcf_df,
      .data$VAF_TUMOR >= config$somatic_snv$allelic_support$tumor_af_min)
  }
  if (!any(is.na(vcf_df$AF_CONTROL))) {
    vcf_df <- dplyr::filter(
      vcf_df,
      .data$AF_CONTROL <= config$somatic_snv$allelic_support$control_af_max)
  }
  if (!any(is.na(vcf_df$DP_CONTROL))) {
    vcf_df <- dplyr::filter(
      vcf_df,
      .data$DP_CONTROL >= config$somatic_snv$allelic_support$control_dp_min)
  }
  n_removed <- n_before_dp_af_filtering - nrow(vcf_df)
  percentage <- round(as.numeric((n_removed / n_before_dp_af_filtering) * 100),
                      digits = 2)
  pcgrr::log4r_info(
    paste0("Removed ", n_removed,
           " tumor variants (", percentage,
           "%) based on thresholds for allelic depth/fraction"
    )
  )

  return(vcf_df)
}

#' Function that writes a VCF intended for mutational signature analysis
#'
#' @param calls data frame with calls
#' @param sample_name sample name
#' @param output_directory Output directory for output file
#' @param vcf_fname filename for VCF
#' @param snv_only logical, if TRUE only SNVs are written to VCF
#'
#'
#' @export
write_processed_vcf <- function(calls,
                                sample_name = NULL,
                                output_directory = NULL,
                                vcf_fname = NULL,
                                snv_only = TRUE) {

  pcgrr::log4r_info("Writing VCF file with input calls for signature analysis")

  header_lines <-
    c("##fileformat=VCFv4.2",
      paste0("##INFO=<ID=SAMPLE_ID,Number=1,Type=String,",
             "Description=\"Sample identifier\">"),
      "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")

  vcf_df <- calls |>
    dplyr::select(c("CHROM","POS","REF","ALT"))

  if(snv_only == TRUE){
    vcf_df <- vcf_df |>
      dplyr::filter(nchar(.data$REF) == 1 & nchar(.data$ALT) == 1) |>
      dplyr::filter(stringr::str_detect(.data$REF,"^(A|C|T|G)$") |
                    stringr::str_detect(.data$ALT,"^(A|C|T|G)$"))
  }
  vcf_df <- vcf_df |>
    dplyr::mutate(
      QUAL = ".",
      FILTER = "PASS",
      ID = ".",
      INFO = paste0(
        "SAMPLE_ID=", sample_name)) |>
    dplyr::distinct()

  sample_vcf_content_fname <-
    file.path(output_directory,
              paste0(sample_name, ".",
                     stringi::stri_rand_strings(
                        1, 15, pattern = "[A-Za-z0-9]"),
                     ".vcf_content.tsv"))
  write(header_lines, file = vcf_fname, sep = "\n")

  sample_vcf <- vcf_df[, c("CHROM", "POS", "ID", "REF",
                          "ALT", "QUAL", "FILTER", "INFO")] |>
    dplyr::mutate(REF = as.character(.data$REF)) |>
    dplyr::mutate(ALT = as.character(.data$ALT))

  options(scipen = 999)
  utils::write.table(sample_vcf, file =
                sample_vcf_content_fname,
              sep = "\t", col.names = F,
              quote = F, row.names = F)


  system(paste0("cat ", sample_vcf_content_fname,
                " | egrep -v \"^[XYM]\" | sort -k1,1n -k2,2n -k4,4 -k5,5  >> ",
                vcf_fname))
  system(paste0("cat ", sample_vcf_content_fname,
                " | egrep \"^[XYM]\" | sort -k1,1 -k2,2n -k4,4 -k5,5 >> ",
                vcf_fname))

  system(paste0("bgzip -f ", vcf_fname))
  system(paste0("tabix -p vcf ", vcf_fname, ".gz"))

  system(paste0("rm -f ", sample_vcf_content_fname))
}


#' A function that detects whether the sample name in
#' variant data frame is unique (as present in column name VCF_SAMPLE_ID), throws an error if
#' multiple sample names are present for the CPSR workflow
#'
#' @param df VCF data frame
#' @param sample_name name of sample identifier
#' @param cpsr logical indicating CPSR workflow
#' @return df Vranges object
#'
#' @export
detect_vcf_sample_name <- function(df, sample_name = NULL, cpsr = FALSE) {
  stopifnot(is.data.frame(df) & !is.null(sample_name))
  if ("VCF_SAMPLE_ID" %in% colnames(df)) {
    unique_sample_names <- unique(df$VCF_SAMPLE_ID)
    pcgrr::log4r_info(paste0("Found the following VCF sample names: ",
                             paste(unique_sample_names, collapse = ", ")))

    if (length(unique_sample_names) > 1 & cpsr == T) {
      pcgrr::log4r_info(paste0("Found more than one sample name - VCF with somatic ",
                     "calls? Expecting single sample germline VCF for CPSR"))
      stop()
    }
  }
  df <- df |>
    dplyr::mutate(VCF_SAMPLE_ID = sample_name)
  return(df)
}

#' Write messages to logs at a given priority level
#'
#' See [log4r::levellog()]
#' @param msg Message to log.
#'
#' @export
log4r_info <- function(msg) {
  log4r_logger <- getOption("PCGRR_LOG4R_LOGGER")
  log4r::info(log4r_logger, msg)
}

#' Write messages to logs at a given priority level
#'
#' See [log4r::levellog()]
#' @inheritParams log4r_info
#' @export
log4r_debug <- function(msg) {
  log4r_logger <- getOption("PCGRR_LOG4R_LOGGER")
  log4r::debug(log4r_logger, msg)
}

#' Write messages to logs at a given priority level
#'
#' See [log4r::levellog()]
#' @inheritParams log4r_info
#' @export
log4r_warn <- function(msg) {
  log4r_logger <- getOption("PCGRR_LOG4R_LOGGER")
  log4r::warn(log4r_logger, msg)
}

#' Write messages to logs at a given priority level
#'
#' See [log4r::levellog()]
#' @inheritParams log4r_info
#' @export
log4r_fatal <- function(msg) {
  log4r_logger <- getOption("PCGRR_LOG4R_LOGGER")
  log4r::fatal(log4r_logger, msg)
  stop()
}

#' Get BSgenome Object
#'
#' Gets BSgenome object given a human genome assembly string.
#'
#' @param genome Human genome assembly string: hg38 or hg19.
#'
#' @return BSgenome object.
#'
#' @examples
#' \dontrun{
#' get_genome_obj("hg38")
#' }
#' @export
get_genome_obj <- function(genome) {
  bsgenome <- c(
    grch37 = "BSgenome.Hsapiens.UCSC.hg19",
    grch38 = "BSgenome.Hsapiens.UCSC.hg38",
    hg19 = "BSgenome.Hsapiens.UCSC.hg19",
    hg38 = "BSgenome.Hsapiens.UCSC.hg38"
  )
  pkg <- bsgenome[genome]
  assertthat::assert_that(
    genome %in% names(bsgenome),
    msg = glue::glue(
      "Instead of '{genome}', pick one of: ",
      "{paste(names(bsgenome), collapse = ', ')}"
    )
  )
  if (!pkg_exists(pkg)) {
    stop(glue::glue(
      "{pkg} is not installed on your system.\n",
      "Please install with:\n'BiocManager::install(\"{pkg}\")'\n",
      "(or use 'mamba install -c bioconda bioconductor-bsgenome.hsapiens.ucsc.hgXX' ",
      "if inside a conda environment)."
    ))
  }
  return(eval(parse(text = glue::glue("{pkg}::{pkg}"))))
}

#' Does R Package Exist
#'
#' Checks if the specified R package exists on the local system.
#'
#' @param p The R package to check for.
#' @return TRUE if package exists, FALSE otherwise.
#'
pkg_exists <- function(p) {
  assertthat::assert_that(is.character(p))
  nzchar(system.file(package = p))
}

#' Function that checks the existence of a file
#'
#' @param fname Name of file to check
#'
#' @export
check_file_exists <- function(fname) {

  if (file.exists(fname)) {
    if (file.size(fname) == 0) {
      log4r_fatal(
        paste0("File ", fname, " has zero size - exiting")
      )
    }
  }else{
    log4r_fatal(
      paste0("File ", fname, " does not exist - exiting")
    )
  }
}


#' Create directory
#'
#' @param d Directory to create.
#'
#' @export
mkdir <- function(d) {
  if (!dir.exists(d)) {
    dir.create(d, recursive = TRUE)
  }
  TRUE
}


#' Strip HTML tags
#'
#' Remove HTML tags and comments from text.
#' From https://github.com/yihui/xfun/blob/ccee26/R/string.R#L329.
#'
#' @param x A character vector.
#' @return A character vector with HTML tags and comments stripped off.
#' @export
#' @examples
#' strip_html('<a href="#">Hello <!-- comment -->world!</a>')
strip_html <- function(x) {
    x <- gsub("<!--.*?-->", "", x)
    x <- gsub("<[^>]+>", "", x)
    x
}

#' Export Quarto Environment Variables
#'
#' Export quarto environment variables
#' required when using conda.
#'
#' @param x Path to conda/envs/pcgrr/etc/conda/activate.d/quarto.sh
#' @export
export_quarto_evars <- function(x) {
  vars <- readr::read_delim(
    x, delim = "=", comment = "#", col_names = c("name", "value"), col_types = "cc"
  ) |>
    dplyr::mutate(name = sub("(.* )(QUARTO_.*)", "\\2", .data$name)) |>
    tibble::deframe() |>
    as.list()
  do.call(Sys.setenv, vars)
}
