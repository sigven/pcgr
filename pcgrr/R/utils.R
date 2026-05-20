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

  # if (NROW(var_df) <= 200) {
  #   bin_size = 0.10
  # }
  if (NROW(var_df) > 500) {
    bin_size = 0.025
  }

  i <- 1
  num_bins <- as.integer(1 / bin_size)
  bin_start <- 0
  while (i <= num_bins) {
    bin_end <- bin_start + bin_size
    bin_name <- as.character(paste0(bin_start, " - ", bin_end))
    for (e in c('SNV','deletion','insertion','indel','substitution')) {
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
  for (vclass in unique(af_bin_df$VARIANT_CLASS)) {
    sum_count <- sum(
      af_bin_df$Count[af_bin_df$VARIANT_CLASS == vclass])
    if (sum_count == 0) {
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

  if (NROW(tier_df) > 500) {
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
exclude_non_chrom_variants <- function(vcf_df, chrom_var = "CHROM") {
  invisible(assertthat::assert_that(
    is.data.frame(vcf_df),
    msg = "Argument 'vcf_df' must be of type data.frame"))
  assertable::assert_colnames(
    vcf_df, chrom_var, only_colnames = F, quiet = T)
  vcf_df <- vcf_df |>
    dplyr::mutate(
      !!rlang::sym(chrom_var) := as.character(!!rlang::sym(chrom_var)))
  n_before_exclusion <- nrow(vcf_df)
  nuc_chromosomes_df <- data.frame(
    c(as.character(seq(1:22)), "X", "Y"),
    stringsAsFactors = F)
  colnames(nuc_chromosomes_df) <- c(chrom_var)
  vcf_df <- dplyr::semi_join(
    vcf_df, nuc_chromosomes_df, by = chrom_var)
  n_after_exclusion <- nrow(vcf_df)
  if (n_before_exclusion - n_after_exclusion > 0) {
    log4r_info(
      paste0("Excluding n = ",
             n_before_exclusion - n_after_exclusion,
             " variant(s) from non-nuclear chromosomes/scaffolds"))
  }
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
  if (nrow(vcf_df) == 0) {
    return(vcf_df)
  }
  vcf_df <- vcf_df |>
    dplyr::mutate(
      !!rlang::sym(chrom_var) :=
        factor(!!rlang::sym(chrom_var),
               ordered = T,
               levels = c(as.character(seq(1:22)), "X", "Y"))) |>
    dplyr::arrange(
      !!rlang::sym(chrom_var),
      !!rlang::sym(pos_var)) |>
    dplyr::mutate(
      !!rlang::sym(chrom_var) :=
        as.character(!!rlang::sym(chrom_var)))

  return(vcf_df)
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


#' Function that filters variant set on depth and allelic fraction
#' according to settings provided by user (tumor and control)
#'
#' @param vcf_df data frame with variants
#' @param config list with PCGR configuration settings
#'
#' @return vcf_df
#'
#' @export
filter_read_support <- function(vcf_df, config = NULL) {

  n_before_dp_af_filtering <- nrow(vcf_df)
  actionable_or_hotspot_vars <- data.frame()

  if ("ACTIONABILITY_TIER" %in% colnames(vcf_df) &
     "MUTATION_HOTSPOT" %in% colnames(vcf_df) &
     "VAR_ID" %in% colnames(vcf_df) &
     "SYMBOL" %in% colnames(vcf_df) &
     "CONSEQUENCE" %in% colnames(vcf_df) &
     "HGVSP" %in% colnames(vcf_df)) {
    actionable_or_hotspot_vars <-
      vcf_df |>
      dplyr::filter(
        (!is.na(.data$ACTIONABILITY_TIER) &
           .data$ACTIONABILITY_TIER <= 2) |
          (!is.na(.data$MUTATION_HOTSPOT))) |>
      dplyr::select(c("VAR_ID","SYMBOL","CONSEQUENCE","HGVSP")) |>
      dplyr::distinct()
  }

  warns_tumor <- 0
  tumor_filtering_set <- FALSE
  for (col in c("DP_TUMOR","VAF_TUMOR","AD_TUMOR")) {
    if (col %in% colnames(vcf_df)) {
      if (!is.null(config$somatic_snv$allelic_support$tumor_dp_min) |
         !is.null(config$somatic_snv$allelic_support$tumor_af_min) |
         !is.null(config$somatic_snv$allelic_support$tumor_ad_min)) {
        tumor_filtering_set <- TRUE
        if (any(is.na(vcf_df[, col]))) {
          log4r_warn(
            paste0("Column '", col,
                   "' contains <NA> values - unable to perform filtering on depth/allelic support"))
          warns_tumor <- warns_tumor + 1
        }
      }
    }
  }

  warns_control <- 0
  control_filtering_set <- FALSE
  for (col in c("DP_CONTROL","VAF_CONTROL","AD_CONTROL")) {
    if (col %in% colnames(vcf_df)) {
      if (!is.null(config$somatic_snv$allelic_support$control_dp_min) |
         !is.null(config$somatic_snv$allelic_support$control_af_max) |
         !is.null(config$somatic_snv$allelic_support$control_ad_max)) {
        control_filtering_set <- TRUE
        if (any(is.na(vcf_df[, col]))) {
          log4r_warn(
            paste0("Column '", col,
                   "' contains <NA> values - unable to perform filtering on depth/allelic support"))
          warns_control <- warns_control + 1
        }
      }
    }
  }


  if ("DP_TUMOR" %in% colnames(vcf_df)) {
    if (!any(is.na(vcf_df$DP_TUMOR)) &
        !is.null(config$somatic_snv$allelic_support$tumor_dp_min)) {
      #cat("DP_TUMOR - Number of rows before filtering - ", nrow(vcf_df), "\n")
      #cat("DP_TUMOR - Filtering criteria used: ", config$somatic_snv$allelic_support$tumor_dp_min, "\n")
      vcf_df <- dplyr::filter(
        vcf_df,
        .data$DP_TUMOR >= config$somatic_snv$allelic_support$tumor_dp_min)
      #cat("DP_TUMOR - Number of rows after filtering: ", nrow(vcf_df), "\n")
    }
  }
  if ("VAF_TUMOR" %in% colnames(vcf_df)) {
    if (!any(is.na(vcf_df$VAF_TUMOR)) &
        !is.null(config$somatic_snv$allelic_support$tumor_af_min)) {
      #cat("VAF_TUMOR - Number of rows before filtering - ", nrow(vcf_df), "\n")
      #cat("VAF_TUMOR - Filtering criteria used: ", config$somatic_snv$allelic_support$tumor_af_min, "\n")
      vcf_df <- dplyr::filter(
        vcf_df,
        .data$VAF_TUMOR >= config$somatic_snv$allelic_support$tumor_af_min)
      #cat("VAF_TUMOR - Number of rows after filtering: ", nrow(vcf_df), "\n")
    }
  }
  if ("AD_TUMOR" %in% colnames(vcf_df)) {

    if (!any(is.na(vcf_df$AD_TUMOR)) &
        !is.null(config$somatic_snv$allelic_support$tumor_ad_min)) {
      #cat("AD_TUMOR - Number of rows before filtering - ", nrow(vcf_df), "\n")
      #cat("AD_TUMOR - Filtering criteria used: ", config$somatic_snv$allelic_support$tumor_ad_min, "\n")
      vcf_df <- dplyr::filter(
        vcf_df,
        .data$AD_TUMOR >= config$somatic_snv$allelic_support$tumor_ad_min)
      #cat("AD_TUMOR - Number of rows after filtering: ", nrow(vcf_df), "\n")
    }

  }

  if ("VAF_CONTROL" %in% colnames(vcf_df)) {

    if (!any(is.na(vcf_df$VAF_CONTROL)) &
        !is.null(config$somatic_snv$allelic_support$control_af_max)) {
      #cat("VAF_CONTROL - Number of rows before filtering - ", nrow(vcf_df), "\n")
      #cat("VAF_CONTROL - Filtering criteria used: ", config$somatic_snv$allelic_support$control_af_max, "\n")
      vcf_df <- dplyr::filter(
        vcf_df,
        .data$VAF_CONTROL <= config$somatic_snv$allelic_support$control_af_max)
      #cat("VAF_CONTROL - Number of rows after filtering: ", nrow(vcf_df), "\n")
    }

  }
  if ("DP_CONTROL" %in% colnames(vcf_df)) {
    if (!any(is.na(vcf_df$DP_CONTROL)) &
        !is.null(config$somatic_snv$allelic_support$control_dp_min)) {
      #cat("Number of rows before filtering - ", nrow(vcf_df), "\n")
      #cat("DP_CONTROL - Filtering criteria used: ", config$somatic_snv$allelic_support$control_dp_min, "\n")
      vcf_df <- dplyr::filter(
        vcf_df,
        .data$DP_CONTROL >= config$somatic_snv$allelic_support$control_dp_min)
      #cat("DP_CONTROL - Number of rows after filtering: ", nrow(vcf_df), "\n")
    }
  }

  if ("AD_CONTROL" %in% colnames(vcf_df)) {
    if (!any(is.na(vcf_df$AD_CONTROL)) &
        !is.null(config$somatic_snv$allelic_support$control_ad_max)) {
      #cat("AD_CONTROL - Number of rows before filtering - ", nrow(vcf_df), "\n")
      #cat("AD_CONTROL - Filtering criteria used: ", config$somatic_snv$allelic_support$control_ad_max, "\n")
      vcf_df <- dplyr::filter(
        vcf_df,
        .data$AD_CONTROL <= config$somatic_snv$allelic_support$control_ad_max)
      #cat("AD_CONTROL - Number of rows after filtering: ", nrow(vcf_df), "\n")
    }
  }

  if (warns_tumor == 0 & tumor_filtering_set == TRUE) {

    filtering_criteria_used_tumor <- c()
    if (!is.null(config$somatic_snv$allelic_support$tumor_dp_min)) {
      filtering_criteria_used_tumor <-
        c(filtering_criteria_used_tumor,
          paste0("min_dp_tumor = ",
                 config$somatic_snv$allelic_support$tumor_dp_min))
    }
    if (!is.null(config$somatic_snv$allelic_support$tumor_af_min)) {
      filtering_criteria_used_tumor <-
        c(filtering_criteria_used_tumor,
          paste0("min_af_tumor = ",
                 config$somatic_snv$allelic_support$tumor_af_min))
    }
    if (!is.null(config$somatic_snv$allelic_support$tumor_ad_min)) {
      filtering_criteria_used_tumor <-
        c(filtering_criteria_used_tumor,
          paste0("min_ad_tumor = ",
                 config$somatic_snv$allelic_support$tumor_ad_min))
    }


    log4r_info(paste0(
      "Tumor variant filtering based on allelic depth/fraction: ",
      paste(filtering_criteria_used_tumor, collapse = ", ")
    ))
  }

  if (warns_control == 0 & control_filtering_set == TRUE) {

    filtering_criteria_used_control <- c()
    if (!is.null(config$somatic_snv$allelic_support$control_dp_min)) {
      filtering_criteria_used_control <-
        c(filtering_criteria_used_control,
          paste0("min_dp_control = ",
                 config$somatic_snv$allelic_support$control_dp_min))
    }
    if (!is.null(config$somatic_snv$allelic_support$control_af_max)) {
      filtering_criteria_used_control <-
        c(filtering_criteria_used_control,
          paste0("max_af_control = ",
                 config$somatic_snv$allelic_support$control_af_max))
    }
    if (!is.null(config$somatic_snv$allelic_support$control_ad_max)) {
      filtering_criteria_used_control <-
        c(filtering_criteria_used_control,
          paste0("max_ad_control = ",
                 config$somatic_snv$allelic_support$control_ad_max))
    }

    log4r_info(paste0(
      "Tumor variant filtering based on allelic depth/fraction: ",
      paste(filtering_criteria_used_control, collapse = ", ")
    ))
  }


  n_removed <- n_before_dp_af_filtering - nrow(vcf_df)
  percentage <- round(as.numeric((n_removed / n_before_dp_af_filtering) * 100),
                      digits = 2)
  log4r_info(
    paste0("Excluded n = ", n_removed,
           " tumor variants (", percentage,
           "% of total) based on thresholds for AF/AD/DP support"
    )
  )

  ## Issue warning if no variants are left after filtering
  if (NROW(vcf_df) == 0) {
    log4r_warn(
      paste0("NO tumor variants left after filtering on AF/AD/DP support - ",
             "check input VCF and/or filtering settings"))
  }

  if (NROW(actionable_or_hotspot_vars) > 0) {

    missed_actionable_vars <-
      actionable_or_hotspot_vars |>
      dplyr::anti_join(
        vcf_df, by = c("VAR_ID")) |>
      dplyr::distinct()

    if (NROW(missed_actionable_vars) > 0) {
      log4r_warn(
        paste0("N = ", NROW(missed_actionable_vars),
               " actionable variants/mutation hotspots excluded after filtering on AF/AD/DP support - ",
               "check input VCF and/or filtering settings"))
      log4r_warn(
        paste0("Missed actionable variants/mutation hotspots: ",
               paste0(missed_actionable_vars$SYMBOL,
                      " (",missed_actionable_vars$CONSEQUENCE,
                      ":",missed_actionable_vars$HGVSP,")",
                      collapse = ", "))
      )

    }
  }

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

  log4r_info("Writing VCF file with input calls for signature analysis")

  header_lines <-
    c("##fileformat=VCFv4.2",
      paste0("##INFO=<ID=SAMPLE_ID,Number=1,Type=String,",
             "Description=\"Sample identifier\">"),
      "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")

  vcf_df <- calls |>
    dplyr::select(c("CHROM","POS","REF","ALT"))

  if (snv_only == TRUE) {
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

  unlink(sample_vcf_content_fname, force = TRUE)
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
    log4r_info(paste0("Found the following VCF sample names: ",
                             paste(unique_sample_names, collapse = ", ")))

    if (length(unique_sample_names) > 1 & cpsr == T) {
      log4r_info(paste0("Found more than one sample name - VCF with somatic ",
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
  if (is.null(log4r_logger)) { message(msg); return(invisible(NULL)) }
  log4r::info(log4r_logger, msg)
}

#' Write messages to logs at a given priority level
#'
#' See [log4r::levellog()]
#' @inheritParams log4r_info
#' @export
log4r_debug <- function(msg) {
  log4r_logger <- getOption("PCGRR_LOG4R_LOGGER")
  if (is.null(log4r_logger)) { return(invisible(NULL)) }
  log4r::debug(log4r_logger, msg)
}

#' Write messages to logs at a given priority level
#'
#' See [log4r::levellog()]
#' @inheritParams log4r_info
#' @export
log4r_warn <- function(msg) {
  log4r_logger <- getOption("PCGRR_LOG4R_LOGGER")
  if (is.null(log4r_logger)) { warning(msg, call. = FALSE); return(invisible(NULL)) }
  log4r::warn(log4r_logger, msg)
}

#' Write messages to logs at a given priority level
#'
#' See [log4r::levellog()]
#' @inheritParams log4r_info
#' @export
log4r_fatal <- function(msg) {
  log4r_logger <- getOption("PCGRR_LOG4R_LOGGER")
  if (is.null(log4r_logger)) { stop(msg, call. = FALSE) }
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
      "(or use 'conda install -c bioconda bioconductor-bsgenome.hsapiens.ucsc.hgXX' ",
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

#' Convert Hex Color to RGBA
#'
#' @param hex Hex color code (e.g. "#FF5733").
#' @param alpha Alpha transparency value (0 to 1).
#'
#' @return RGBA color string (e.g. "rgba(255,87,51,0.5)").
#' @export
#'
#'
hex_to_rgba <- function(hex, alpha = 1) {
  rgb <- grDevices::col2rgb(hex)
  sprintf(
    "rgba(%d,%d,%d,%.2f)", rgb[1],
    rgb[2], rgb[3], alpha)
}


#' Plotly Pie Chart - variant statistics
#'
#' Function that generates a pie chart using plotly
#' for a given category in a data frame
#'
#' @param df_variant_stats Data frame with variant statistics
#' @param category Category for pie chart (e.g. CODING_STATUS)
#' @param plot_margin_top Top margin
#' @param plot_margin_bottom Bottom margin
#' @param plot_margin_left Left margin
#' @param plot_margin_right Right margin
#' @param font_family Font family
#' @param font_size Font size
#' @param pie_line_width Line width for pie chart segments
#' @param opacity_filtered_categories Opacity for filtered categories
#' @param hole_size_pie Hole size for pie chart (0 to 1)
#' @return Plotly pie chart object
#'
#' @export
#'
plotly_pie_chart <- function(
    df_variant_stats = NULL,
    category = "CODING_STATUS",
    color_palette = pcgrr::color_palette,
    plot_margin_top = 50,
    plot_margin_bottom = 20,
    plot_margin_left = 20,
    plot_margin_right = 20,
    font_family = "Helvetica",
    font_size = 15,
    pie_line_width = 3,
    opacity_filtered_categories = 0.4,
    hole_size_pie = 0.4) {

  invisible(assertthat::assert_that(
    !is.null(df_variant_stats) &
      is.data.frame(df_variant_stats) &
      "N" %in% colnames(df_variant_stats) &
      "Pct" %in% colnames(df_variant_stats) &
      !is.null(category) &
      is.character(category) &
      category %in% colnames(df_variant_stats)
  ))

  data <- df_variant_stats
  data$LABELS <-
    as.character(data[[category]])  # force to character

  p <- plotly::plot_ly(
    data,
    marker = list(
      colors =
        color_palette$multi$values,
      line = list(
        color = "#FFFFFF",
        width = pie_line_width))) |>
    plotly::add_pie(
      labels = ~LABELS,   # safe because already character
      values = ~N,
      textinfo = "Pct",
      hole = hole_size_pie
    ) |>
    plotly::layout(
      ## legend at the bottom
      legend = list(
        orientation = "h",
        font = t,
        xanchor = "center",
        x = 0.5),
      margin = list(
        l = plot_margin_left,  # left margin
        r = plot_margin_right,  # right margin
        b = plot_margin_bottom,  # bottom margin
        t = plot_margin_top   # top margin
      ))

  return(p)

}


